#!/usr/bin/env python3
import numpy as np
import netCDF4 as nc
import os
import csv
import sys
import argparse
import pyproj
from scipy.interpolate import RegularGridInterpolator, griddata
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
from scipy.ndimage import uniform_filter
import time
from copy import deepcopy
from inspect import currentframe, getframeinfo
import datetime
from scipy import optimize
from typing import List, Dict

PROGRAM_NAME = 'geo_eval'
VERSION = "1.3.1"
bad_float = -32767.0
running_time = datetime.datetime.now()
RUNNING_TIME = running_time.strftime("%Y-%m-%dT%H%M%S")
subpixels_default_max = 30
subpixels_default = 15
window_size_defalut = 5
min_number_of_points = 100
cf = currentframe()
filename_script = getframeinfo(cf).filename

def check_for_bad_geolocation(coord):
    if coord[0] == bad_float or np.isnan(coord[0]):
        return True
    if coord[1] == bad_float or np.isnan(coord[1]):
        return True
    return False

def get_cloud_mask(cld_file: str):
    if not os.path.exists(cld_file):
        print(
            f"File {cld_file} does not exist. See line {get_linenumber()} in {filename_script}. Exiting")
        sys.exit(1)
    with nc.Dataset(cld_file) as f:
        return f["geophysical_data/cloud_flag"][:]


def read_file_list(path_chips: str) -> List[str]:
    with open(path_chips) as file_txt:
        file_list_chip = file_txt.read().splitlines()
    for file in file_list_chip:
        if not os.path.exists(file):
            print(
                f"File {file} does not exist. See line {get_linenumber()} in {filename_script}. Exiting")
            sys.exit(1)
    return file_list_chip


def get_scantime(year: List[int], day: List[int], msec: List[int]) -> np.ndarray:
    scantime = np.zeros(len(year))
    for i, (iyear, iday, imsec) in enumerate(zip(year, day, msec)):
        date_scan = datetime.datetime.strptime(f'{iyear} {iday}', '%Y %j')
        scantime[i] = time.mktime(date_scan.timetuple(
        )) + date_scan.microsecond / 1.0e6 + float(imsec) / 1.0e3
    return scantime


def get_viirs_modis_lat_lon_scantime(f: nc.Dataset):
    lat = np.array(f["navigation_data/latitude"][:], dtype=np.float64)
    lon = np.array(f["navigation_data/longitude"][:], dtype=np.float64)
    year: list[int] = list(
        np.array(f["scan_line_attributes/year"][:], dtype=np.int32))
    day: list[int] = list(
        np.array(f["scan_line_attributes/day"][:], dtype=np.int32))
    msec: list[int] = list(
        np.array(f["scan_line_attributes/msec"][:], dtype=np.int32))
    scantime = get_scantime(year, day, msec)
    return lat, lon, scantime


def get_linenumber():
    cf = currentframe()
    return cf.f_back.f_lineno


def read_file_l1b(path: str):
    with nc.Dataset(path) as f:
        if f.instrument == "MODIS":
            data = np.array(
                f["geophysical_data/Lt_645"][:], dtype=np.float64)
            lat, lon, scantime = get_viirs_modis_lat_lon_scantime(f)
            return lat, lon, data, scantime, None
        if f.instrument == "VIIRS":
            if f.platform == "Suomi-NPP":
                wave = 671
            else:
                wave = 667
            data = np.array(
                f[f"geophysical_data/Lt_{wave}"][:], dtype=np.float64)
            l2_flags_var = f["geophysical_data/l2_flags"]
            flag_masks = l2_flags_var.flag_masks
            flag_meanings = l2_flags_var.flag_meanings.split(sep=' ')
            bit_bowtie = 1
            for bit, name in zip(flag_masks, flag_meanings):
                if name == "BOWTIEDEL":
                    bit_bowtie = bit
            mask = (l2_flags_var[:] & bit_bowtie) > 0
            data[mask] = np.NaN
            lat, lon, scantime = get_viirs_modis_lat_lon_scantime(f)
            return lat, lon, data, scantime, None
        if f.instrument == "OCI" or f.instrument == "OCIS":
            red_wv = 650.0
            red_wavelength = f["sensor_band_parameters/red_wavelength"][:]
            index = np.abs(red_wavelength - red_wv).argmin()
            data = np.array(
                f["observation_data/rhot_red"][:][index, :, :], dtype=np.float64)
            lat = np.array(f["geolocation_data/latitude"][:], dtype=np.float64)
            lon = np.array(f["geolocation_data/longitude"]
                           [:], dtype=np.float64)
            scantime = np.array(f["scan_line_attributes/time"]
                                [:], dtype=np.float64)
            time_coverage_start = f.time_coverage_start
            if time_coverage_start[-1] == 'Z':
                time_coverage_start = time_coverage_start[:-1]
            date_scan = datetime.datetime.strptime(
                time_coverage_start, '%Y-%m-%dT%H:%M:%S.%f')
            scantime += time.mktime(datetime.datetime(date_scan.year, date_scan.month, date_scan.day).timetuple(
            ))
            tilt = f["navigation_data/tilt_angle"][:]
            return lat, lon, data, scantime, tilt
        print(
            f"-Error- : instrument {f.instrument} not found. See line {get_linenumber()} in {filename_script}")
        sys.exit(1)


def get_position(lat: np.ndarray, lon: np.ndarray, feature_lat, feature_lon):
    diff = np.abs(lat - feature_lat) + np.abs(lon - feature_lon)
    (x, y) = np.unravel_index(diff.argmin(), diff.shape)
    point = Point(feature_lat, feature_lon)
    w = 1
    while w < 32:
        i_min = max(0, x - w)
        i_max = min(lat.shape[0] - 1, x + w)
        j_min = max(0, y - w)
        j_max = min(lat.shape[1] - 1, y + w)
        for i in range(i_min, i_max):
            for j in range(j_min, j_max):
                p_0 = (lat[i, j], lon[i, j])
                p_1 = (lat[i, j + 1], lon[i, j + 1])
                p_2 = (lat[i + 1, j + 1], lon[i + 1, j + 1])
                p_3 = (lat[i + 1, j], lon[i + 1, j])
                if check_for_bad_geolocation(p_0) or check_for_bad_geolocation(p_1) or check_for_bad_geolocation(p_2) or check_for_bad_geolocation(p_3):
                    print(f"Bad geolocation in lat, lon array for line {i}, pixel {j}")
                    return None
                polygon = Polygon([p_0, p_1, p_2, p_3, p_0])
                # a bit more accurate determination of the scan and pixel number
                if polygon.contains(point):
                    dist_dict = {}
                    min_dist = point.distance(Point(p_0))
                    i_min = i
                    j_min = j
                    dist_dict[(i, j + 1)] = point.distance(Point(p_1))
                    dist_dict[(i + 1, j + 1)] = point.distance(Point(p_2))
                    dist_dict[(i + 1, j)] = point.distance(Point(p_3))
                    for (i_k, j_k), dist in dist_dict.items():
                        if dist < min_dist:
                            min_dist = dist
                            i_min = i_k
                            j_min = j_k
                    return i_min, j_min, p_0, p_1, p_2, p_3, point
        w = w * 2
    print(
        f"-Warning- : Could not locate the chip with lat={feature_lat} and lon={feature_lon} within the granule. See {get_linenumber()} in {filename_script}")
    return None


def lorentzian(latlon, x0, y0, GammaX, GammaY, Factor, alpha, Amplitude, offset):
    lat, lon = latlon
    lat = lat - x0
    lon = lon - y0
    x = lat * np.cos(alpha) + lon * np.sin(alpha)
    y = - lat * np.sin(alpha) + lon * np.cos(alpha)
    val: np.ndarray = offset + Amplitude / \
        ((x - x0) ** 2 / GammaX ** 2 + (y - y0) ** 2 / GammaY ** 2 + Factor)
    return val.ravel()


def read_chip_file(path: str):
    with nc.Dataset(path) as f:
        ecrx = f["ECR_x_coordinate_array"][:]
        ecry = f["ECR_y_coordinate_array"][:]
        ecrz = f["ECR_z_coordinate_array"][:]
        ecef = {'proj': 'geocent', 'ellps': 'WGS84', 'datum': 'WGS84'}
        lla = {'proj': 'latlong', 'ellps': 'WGS84', 'datum': 'WGS84'}
        transformer = pyproj.Transformer.from_crs(ecef, lla)
        lon, lat, _ = transformer.transform(ecrx, ecry, ecrz)
        return f["Band_1"][:], lat, lon, f.FEATURELATITUDE, f.FEATURELONGITUDE


def parse_command_line():
    parser = argparse.ArgumentParser(
        prog=PROGRAM_NAME,
        description='''\
                        Output:
                        scantime - scantime
                        granule - input L1B granules
                        chip - input chip file
                        feature_latitude - latitude of the center of the chip
                        feature_longitude - longitude of the center of the chip
                        delta_latitude - latitude shift,  latitude coordinate of the peak of NCC
                        delta_longitude -  longitude shift, longitude coordinate of the peak of NCC
                        pixel_number - pixel number between 0 and pixels_per_line-1 (or ccd_pixels-1 for OCI)
                        confidence_level - the maximum value of the normalized cross correlation (NCC) within a sliding lat/lon window
                        tilt - tilt angle
                        fit_delta_latitude - latitude shift from the Lorentzian fit of NCC
                        fit_delta_longitude - longitude shift from the Lorentzian fit of NCC
                        r_squared - R_squared of the the Lorentzian fit
                        ''',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''\
                Additional output info:
                scantime is a scantime of the line where the feature pixel was found.
                delta_latitude and delta_longitude mean that a point in the chip with (lat,lon) coordinates correspond to a point in the granule with (lat + delta_latitude, lon + delta_longitude).
                The maximum possible value of confidence_level is 1 which means 100% match. Higher value is better. Values lower 0.5 indicate poor match.
                fit_delta_latitude, fit_delta_longitude and r_squared can be NaN if the fit does not converge.
         ''')
    parser.add_argument("-v", "--version", action="version", version=VERSION)
    parser.add_argument(
        "--ifile", help="input L1B netcdf file", type=str, required=True)
    parser.add_argument(
        "--ilist",
        help="input chip file list, each filename should be separated by a new line terminator. No two filenames can be put in the same line",
        type=str, required=True)
    parser.add_argument("--ofile", help="output csv", type=str, required=True)
    parser.add_argument("--cloudmask", help="input cloudmask file", type=str)
    parser.add_argument("--odebugfile", help="output debug file", type=str)
    parser.add_argument("--subpixels", help="number of subpixels", type=int)
    parser.add_argument(
        "--window_size", help="the size of sliding window in sensor pixels", type=int)
    parser.add_argument(
        "--confidence", help="minimum confidence level to be included. Default is 0.5", type=float, default=0.5)
    parser.add_argument('--debug', action='store_true')
    parser.add_argument('--no-debug', dest='debug', action='store_false')
    parser.set_defaults(debug=False)
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    args = parser.parse_args()
    return args


def ini_input_data(lat, lon, path_chip):
    tm_band, lat_chip, lon_chip, feature_lat, feature_lon = read_chip_file(
        path_chip)
    return get_position(lat, lon, feature_lat, feature_lon), tm_band, lat_chip, lon_chip


def make_rectangular_grid(lat: np.ndarray, lon: np.ndarray, data: np.ndarray):
    len_x, len_y = lat.shape

    lat_dir = lat[1, 0] - lat[0, 0] > 0
    lon_dir = lon[0, 1] - lon[0, 0] > 0
    # find max and min
    lat_max = np.min(lat.max(axis=0))
    lat_min = np.max(lat.min(axis=0))
    lon_max = np.min((lon.max(axis=1)))
    lon_min = np.max((lon.min(axis=1)))
    lat_st = lat_min
    lat_end = lat_max
    lon_st = lon_min
    lon_end = lon_max
    if not lat_dir:
        lat_st = lat_max
        lat_end = lat_min
    if not lon_dir:
        lon_st = lon_max
        lon_end = lon_min

    lat_lin = np.linspace(lat_st, lat_end, len_x, endpoint=True)
    lon_lin = np.linspace(lon_st, lon_end, len_x, endpoint=True)

    lon_grid, lat_grid = np.meshgrid(lon_lin, lat_lin)
    data_gridded = griddata((lat.flatten(), lon.flatten()), data.flatten(),
                            (lat_grid.flatten(), lon_grid.flatten()), method="nearest")
    return data_gridded.reshape(len_x, len_y), lat_lin, lon_lin


def extract_scene(lat, lon, i, j, w, h, image=None):
    sc_lat = deepcopy(lat[i:i + w, j:j + h])
    sc_lon = deepcopy(lon[i:i + w, j:j + h])
    if image is not None:
        sc_image = deepcopy(image[i:i + w, j:j + h])
        return sc_lat, sc_lon, sc_image
    else:
        return sc_lat, sc_lon


def get_match_layers(rbs, point, sc_lon, sc_lat, sc_image, shift_lon,
                     shift_lat):
    sc_lon_shift = sc_lon - point.y - shift_lon
    sc_lat_shift = sc_lat - point.x - shift_lat
    sc_image_copy = deepcopy(sc_image)
    tm_sc_band_n = rbs(np.column_stack(
        [sc_lat_shift.flatten(), sc_lon_shift.flatten()]))
    tm_sc_band_n = tm_sc_band_n.reshape(sc_image_copy.shape)
    sc_image_copy[np.isnan(tm_sc_band_n)] = np.NaN
    tm_sc_band_n[np.isnan(sc_image_copy)] = np.NaN
    return sc_image_copy, tm_sc_band_n


def get_cross_normalized_correlation(rbs, point, sc_lon, sc_lat, sc_image, shift_lon,
                                     shift_lat):
    sc_image_copy, tm_sc_band_n = get_match_layers(rbs, point, sc_lon, sc_lat, sc_image, shift_lon,
                                                   shift_lat)
    # check if it makes sense to compute it :
    valid_pnts = np.sum(np.isfinite(tm_sc_band_n))
    if valid_pnts < min_number_of_points:
        return -1.0
    # shift mean
    sc_image_copy = sc_image_copy - np.nanmean(sc_image_copy)
    tm_sc_band_n = tm_sc_band_n - np.nanmean(tm_sc_band_n)
    sc_image_copy /= np.nanstd(sc_image_copy)
    tm_sc_band_n /= np.nanstd(tm_sc_band_n)
    ncc = np.nansum(sc_image_copy * tm_sc_band_n) / \
        np.sum(~np.isnan(tm_sc_band_n))
    return ncc


def get_lat_lon_fit_shifts(ncc_original, lat_grid, lon_grid, delta_latitude, delta_longitude):
    ncc_fit = 1 / (1 - ncc_original)
    gamma_lat = np.abs(lat_grid[1] - lat_grid[0]) * 100
    gamma_lon = np.abs(lon_grid[1] - lon_grid[0]) * 100
    latlon = np.meshgrid(lat_grid, lon_grid)
    try:
        p0 = delta_latitude, delta_longitude, gamma_lat, gamma_lon, 1, 0, 100, 0
        popt, pcov = optimize.curve_fit(lorentzian, latlon, ncc_fit.ravel(),
                                        p0=p0, bounds=((-0.1, -0.1, -1, -1, -np.inf, -np.pi / 2, 0, -np.inf),
                                                       (0.1, 0.1, 1, 1, np.inf, np.pi / 2, np.inf, np.inf)))
        x0, y0, GammaX, GammaY, Factor, alpha, Amplitude, offset = popt
        fit = lorentzian(latlon, x0, y0, GammaX, GammaY,
                         Factor, alpha, Amplitude, offset)
        fit = fit.reshape(ncc_fit.shape)
        residual = ncc_fit - fit
        residual *= residual
        ss_res = np.sum(residual)
        ss_data = np.sum((ncc_fit - np.mean(ncc_fit)) ** 2)
        r_squared = 1 - ss_res / ss_data
        return x0, y0, GammaX, GammaY, Factor, alpha, Amplitude, offset, fit, r_squared
    except Exception as e:
        print(e)
        return None


def get_lat_lon_shift(l1b_path: str, chips_list: List[str], output_file_name: str, debug_mode=False, **kwargs):
    lat, lon, image, scantime, tilt_angle = read_file_l1b(l1b_path)
    subpixels_set = subpixels_default  # can't be more than 30, 1km // 30 meters == 30
    # the sliding window size in 1km resolution
    window_size_set = window_size_defalut
    coeff_scale = 1
    min_included_ncc = -1
    if "confidence" in kwargs:
        if not kwargs["confidence"] is None:
            min_included_ncc = kwargs["confidence"]
    w = 30  # // half width of the crop
    h = w

    if "subpixels" in kwargs:
        if not kwargs["subpixels"] is None:
            if kwargs["subpixels"] > subpixels_default_max:
                print("Warning: the subpixel count can't exceed {}".format(
                    subpixels_default_max))
            else:
                subpixels_set = kwargs["subpixels"]
    if "window_size" in kwargs:
        if not kwargs["window_size"] is None:
            if kwargs["window_size"] > w:
                print("Warning: the window_size count can't exceed {}".format(w))
            else:
                window_size_set = kwargs["window_size"]
    # set debug mode if odebug
    if "odebugfile" in kwargs:
        if not kwargs["odebugfile"] is None:
            debug_mode = True

    len_ncc = ((subpixels_set * window_size_set) // 2) * 2 + \
        1  # number of steps, window size
    coeff_scale = subpixels_default_max / subpixels_set
    wide = len_ncc // 2  # half width

    # setting the cloud mask
    if "cloud_mask_file" in kwargs:
        if not kwargs["cloud_mask_file"] is None:
            cld_mask = get_cloud_mask(kwargs["cloud_mask_file"])
            image[cld_mask == 1] = np.NaN
    # data in the output file, csv format
    data_out: List[Dict] = []

    for chip_path in chips_list:
        start_time = time.time()
        # getting chip position within the scene
        position_of_the_chip, tm_band, lat_chip, lon_chip = ini_input_data(lat, lon,
                                                                           chip_path)
        if position_of_the_chip is None:
            print(f"Skipping the chip {chip_path} ...")
            continue
        i_cell, j_cell, p0, p1, p2, p3, point = position_of_the_chip
        if i_cell < 2 * w or (lat.shape[0] - i_cell) < 2 * w:
            print(
                f"The chip {chip_path} is too close to the swath edge : line_number={i_cell}. Skipping ...")
            continue
        if j_cell < 2 * h or (lat.shape[1] - j_cell) < 2 * h:
            print(
                f"The chip {chip_path} is too close to the swath edge : pixel_number={j_cell}. Skipping ...")
            continue
        tb_band_gridded, lat_lin, lon_lin = make_rectangular_grid(
            lat=lat_chip, lon=lon_chip, data=tm_band)
        if debug_mode:
            tb_band_gridded_original = deepcopy(tb_band_gridded)
        # apply sensor view (uniform filter)
        tb_band_gridded = uniform_filter(tb_band_gridded, 33)
        step_size_lat = np.abs(lat_lin[1] - lat_lin[0]) * coeff_scale
        step_size_lon = np.abs(lon_lin[1] - lon_lin[0]) * coeff_scale
        shift_lon_array = np.linspace(
            step_size_lon * (-wide), step_size_lon * wide, len_ncc)
        shift_lat_array = np.linspace(
            step_size_lat * (-wide), step_size_lat * wide, len_ncc)
        rbs = RegularGridInterpolator(
            (lat_lin - point.x, lon_lin - point.y), tb_band_gridded, method="linear", bounds_error=False)
        sc_lat, sc_lon, sc_image = extract_scene(
            lat, lon, i_cell - w, j_cell - h, w * 2 + 1, h * 2 + 1, image)
        ncc_array = np.ones((len_ncc, len_ncc)) * -1.0
        # check if sc_image has enough valid pixels
        valid_pnts = np.sum(np.isfinite(sc_image))
        if valid_pnts < min_number_of_points:
            print(
                f"Number of valid (cloud free) points in the scene {valid_pnts} is too low. Skipping the chip {chip_path}")
            continue
        for ind_lat in range(len_ncc):
            for ind_lon in range(len_ncc):
                shift_lat = shift_lat_array[ind_lat]
                shift_lon = shift_lon_array[ind_lon]
                ncc = get_cross_normalized_correlation(rbs, point, sc_lon, sc_lat, sc_image,
                                                       shift_lon,
                                                       shift_lat)
                ncc_array[ind_lat, ind_lon] = ncc

        print("--- Execution time: {} seconds ---".format(time.time() - start_time))
        print(f"Chip {chip_path} has been processed.")
        (ilat, ilon) = np.unravel_index(ncc_array.argmax(), ncc_array.shape)
        delta_latitude = shift_lat_array[ilat]
        delta_longitude = shift_lon_array[ilon]
        ncc_array[np.isnan(ncc_array)] = -1.0
        max_ncc = np.nanmax(ncc_array)
        if max_ncc < min_included_ncc:
            print(
                f"The confidence level is {max_ncc} which is lower than the threshold = {min_included_ncc}. Skipping ..")
            continue
        if tilt_angle is None:
            tilt = 0.0
        else:
            tilt = tilt_angle[i_cell]
        # computing the Lorentzian fit parameters
        result_out = get_lat_lon_fit_shifts(
            ncc_array, shift_lat_array, shift_lon_array, delta_latitude, delta_longitude)
        if result_out:
            x0_fit, y0_fit, GammaX, GammaY, Factor, alpha, Amplitude, offset, fit, r_squared = result_out
        else:
            print(f"Warning: lorentzian fit for {chip_path} has not converged")
            x0_fit = np.NAN
            y0_fit = np.NAN
            r_squared = np.NAN
            fit = np.ones(ncc_array.shape) * np.NAN
            alpha = np.NAN
            Amplitude = np.NAN
            offset = np.NAN
            Factor = np.NAN
            GammaX = np.NAN
            GammaY = np.NAN
        if debug_mode:
            mode = "a"
            out_ncc_filename = f"ncc_{RUNNING_TIME}_{os.path.basename(l1b_path)}.nc"
            out_ncc_filename = os.path.join(
                os.path.dirname(output_file_name), out_ncc_filename)
            if "odebugfile" in kwargs:
                if not kwargs["odebugfile"] is None:
                    out_ncc_filename = kwargs["odebugfile"]
            if not os.path.exists(out_ncc_filename):
                mode = "w"
            with nc.Dataset(out_ncc_filename, mode=mode) as out_ncc:
                out_ncc: nc.Dataset = out_ncc
                if mode == "w":
                    out_ncc.createDimension("x", len_ncc)
                    out_ncc.createDimension("y", len_ncc)
                    out_ncc.createDimension("lat", sc_image.shape[1])
                    out_ncc.createDimension("lon", sc_image.shape[0])
                    out_ncc.version = VERSION
                    out_ncc.runtime = RUNNING_TIME
                    out_ncc.setncattr("l1b_file", l1b_path)
                    out_ncc.setncattr("chips",",".join(chips_list))
                    for key,val in kwargs.items():
                        if val is not None:
                            out_ncc.setncattr(key,val)
                        else:
                            out_ncc.setncattr(key,"")
                grp: nc.Group = out_ncc.createGroup(
                    os.path.basename(chip_path))
                grp.delta_latitude = delta_latitude
                grp.delta_longitude = delta_longitude
                grp.line = i_cell
                grp.pixel = j_cell
                wkt_polygon_str : str = Polygon([p0, p1, p2, p3]).wkt
                wkt_point_str : str = point.wkt
                grp.setncattr("Point",wkt_point_str)
                grp.setncattr("Bounding_polygon",wkt_polygon_str)
                grp.lat_fit = x0_fit
                grp.lon_fit = y0_fit
                grp.alpha = alpha
                grp.Amplitude = Amplitude
                grp.offset = offset
                grp.Factor = Factor
                grp.GammaX = GammaX
                grp.GammaY = GammaY
                grp.r_squared = r_squared
                grp.tilt = tilt
                grp.feature_latitude = point.x
                grp.feature_longitude = point.y
                grp.createDimension("chip_lat_len", lat_lin.size)
                grp.createDimension("chip_lon_len", lon_lin.size)
                var_grp_local = grp.createVariable(varname="lat_lin", dimensions='chip_lat_len',
                                                   datatype='f4', zlib=True)
                var_grp_local[:] = lat_lin
                var_grp_local = grp.createVariable(varname="lon_lin", dimensions='chip_lon_len',
                                                   datatype='f4', zlib=True)
                var_grp_local[:] = lon_lin
                var_grp_local = grp.createVariable(varname="tm_band", dimensions=('chip_lat_len', 'chip_lon_len'),
                                                   datatype='f4', zlib=True)
                var_grp_local[:] = tb_band_gridded_original
                mat_ncc = grp.createVariable(varname="ncc", dimensions=('x', 'y'),
                                             datatype='f4', zlib=True)
                mat_ncc[:] = ncc_array
                mat_fit = grp.createVariable(varname="fit", dimensions=('x', 'y'),
                                             datatype='f4', zlib=True)
                mat_fit[:] = fit
                lat_array = grp.createVariable(varname="latGrid", dimensions=('x'),
                                               datatype='f4', zlib=True)
                lat_array[:] = shift_lat_array

                lon_array = grp.createVariable(varname="lonGrid", dimensions=('y'),
                                               datatype='f4', zlib=True)
                lon_array[:] = shift_lon_array

                sc_image_crop, chip_downsample = get_match_layers(rbs, point, sc_lon, sc_lat, sc_image,
                                                                  delta_longitude,
                                                                  delta_latitude)
                chip = grp.createVariable(varname="chip", dimensions=('lon', 'lat'),
                                          datatype='f4', zlib=True)
                chip[:] = chip_downsample
                crop = grp.createVariable(varname="crop", dimensions=('lon', 'lat'),
                                          datatype='f4', zlib=True)
                crop[:] = sc_image_crop
                scene = grp.createVariable(varname="scene",
                                           dimensions=('lon', 'lat'),
                                           datatype='f4', zlib=True)
                scene[:] = sc_image
                lat_crop = grp.createVariable(varname="lat",
                                              dimensions=('lon', 'lat'),
                                              datatype='f4', zlib=True)
                lat_crop[:] = sc_lat
                lon_crop = grp.createVariable(varname="lon",
                                              dimensions=('lon', 'lat'),
                                              datatype='f4', zlib=True)
                lon_crop[:] = sc_lon
        granule_name_len = len(os.path.basename(l1b_path))
        chip_name_len = len(os.path.basename(chip_path))
        data_out.append({'scantime'.ljust(14): f"{scantime[i_cell]:14.3f}"})
        data_out[-1]["granule".ljust(granule_name_len)
                     ] = os.path.basename(l1b_path)
        data_out[-1]['chip'.ljust(chip_name_len)] = os.path.basename(chip_path)
        data_out[-1]['feature_latitude'] = f"{point.x:16.7f}"
        data_out[-1]['feature_longitude'] = f"{point.y:17.7f}"
        data_out[-1]['delta_latitude'] = f"{delta_latitude:14.7f}"
        data_out[-1]['delta_longitude'] = f"{delta_longitude:15.7f}"
        data_out[-1]['pixel_number'] = f"{j_cell:12d}"
        data_out[-1]['confidence_level'] = f"{max_ncc:16.7f}"
        data_out[-1]['tilt'.rjust(10)] = f"{tilt:10.5f}"
        if abs(x0_fit) > 1.0 or abs(y0_fit) > 1.0:
            print(f"Fit did not converge for {x0_fit} and {y0_fit}")
            r_squared = np.NAN
        if not np.isnan(r_squared):                
            data_out[-1]['fit_delta_latitude'] = f"{x0_fit:18.7f}"
            data_out[-1]['fit_delta_longitude'] = f"{y0_fit:18.7f}"
            data_out[-1]['r_squared'.rjust(15)] = f"{r_squared:15.7f}"
        else:
            data_out[-1]['fit_delta_latitude'] = "NaN".rjust(18)
            data_out[-1]['fit_delta_longitude'] = "NaN".rjust(18)
            data_out[-1]['r_squared'.rjust(15)] = "NaN"
    if len(data_out) > 0:
        fields = []
        for key in data_out[0]:
            fields.append(key)
        with open(output_file_name, 'w', newline='') as output_csv:
            fields = []
            for key in data_out[0]:
                fields.append(key)
            recorder = csv.DictWriter(output_csv, fieldnames=fields)
            recorder.writeheader()
            recorder.writerows(data_out)
        return 0
    else:
        return 110


if __name__ == "__main__":
    print(PROGRAM_NAME, VERSION)
    args = parse_command_line()
    print("You are running {}, version {}. Current time is {}".format(
        PROGRAM_NAME, VERSION, running_time))
    l1b_path = args.ifile
    chips_list = read_file_list(args.ilist)
    ofile = args.ofile
    cloud_mask_file = args.cloudmask
    window_size = args.window_size
    subpixels = args.subpixels
    ret_code = get_lat_lon_shift(
        l1b_path, chips_list, ofile, debug_mode=args.debug, cloud_mask_file=cloud_mask_file, window_size=window_size,
        subpixels=subpixels, confidence=args.confidence, odebugfile=args.odebugfile)
    sys.exit(ret_code)