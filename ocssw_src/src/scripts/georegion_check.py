#!/usr/bin/env python3
from georegion_utils.georegion_utils import determine_projection, get_bounds, read_wkt_file,check_file_type
import argparse
import datetime
import sys
import os
import geopandas
import netCDF4 as nc
import numpy as np
import pyproj
from shapely import contains_xy, wkt
from shapely.ops import transform
VERSION = "v.1.1.0"
PROGRAM = "georegion_check"

running_time = datetime.datetime.now()
RUNNING_TIME = running_time.strftime("%Y-%m-%dT%H%M%S")

# find netcdf geospatial_bounds attribute
def find_geospatial_bounds_polygon(input_nc_file):
    with nc.Dataset(input_nc_file, "r") as ds:
        if "geospatial_bounds" in ds.ncattrs():
            geospatial_bounds_wkt = ds.getncattr("geospatial_bounds")
            shape_data = geopandas.GeoSeries.from_wkt([geospatial_bounds_wkt])
            poly = []
            for geom in shape_data.geometry:
                poly.append(geom)
            return poly
        else:
            print(f"-E-: No 'geospatial_bounds' attribute found in file '{input_nc_file}'")
            sys.exit(1)


if __name__ == "__main__":
    print(PROGRAM, VERSION)
    parser = argparse.ArgumentParser()
    # input L1/L2 netcdf file, ASCII file or WKT string
    parser.add_argument("--input", "-i", help="Input L1/L2 file, ASCII file or WKT string", required=True)
    # input georegion file
    parser.add_argument("--georegion", "-g", help="Input georegion file", required=True)
    parser.add_argument("--version", "-v", action="version", version="{} {}".format(PROGRAM, VERSION))
    parser.add_argument("--verbose", help="Increase output verbosity", action="store_true")
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    print(
        "You are running {}, version {}. Current time is {}".format(
            PROGRAM, VERSION, running_time
        )
    )
    # prints if the region is inside, outside or intersects the georegion. If not intersects (outside), exits with code 110. If intersects or inside, exits with code 0. If fails, returns code 1
    args = parser.parse_args()

    # first, check if input is a file or a WKT string
    input_arg = args.input
    os.path.isfile(input_arg)
    if os.path.isfile(input_arg):
        # check file type
        file_type = check_file_type(input_arg)
        print("Input file type detected as: {}".format(file_type))
        file_type = check_file_type(input_arg)
        poly = []
        if file_type == "ASCII text":
            shape_data = read_wkt_file(input_arg)
            for geom in shape_data.geometry:
                poly.append(geom)
        elif file_type == "Shapefile":
            shape_data = geopandas.read_file(input_arg)
            for geom in shape_data.geometry:
                poly.append(geom)
        elif file_type == "NetCDF":
            poly.extend(find_geospatial_bounds_polygon(input_arg))
        else:
            print(
                f"-E-: Unsupported file type '{file_type}' for file '{input_arg}'"
            )
            sys.exit(1)
    else:
        # assume WKT string
        shape_data = geopandas.GeoSeries.from_wkt([input_arg])
        poly = []
        for geom in shape_data.geometry:
            poly.append(geom)
    # read georegion netcdf file
    georegion_file = args.georegion
    # check if it exists

    if not os.path.isfile(georegion_file):
        print(f"-E-: Georegion file '{georegion_file}' does not exist")
        sys.exit(1)
    # read lat/lon arrays from georegion file
    with nc.Dataset(georegion_file) as ds:
        if "lat" in ds.variables and "lon" in ds.variables:
            lat = ds.variables["lat"][:]
            lon = ds.variables["lon"][:]
        else:
            print(f"-E-: 'lat' and/or 'lon' variables not found in file '{georegion_file}'")
            sys.exit(1)
    projection_list = []
    lat_lon_bounds = []
    for p in poly:
        proj_type = determine_projection(p.wkt)
        projection_list.append(proj_type)
        lon_max, lon_min, lat_max, lat_min = get_bounds(
            p.wkt, proj_type
        )
        print(
            f"Determined bounds for shape {poly.index(p)+1}/{len(poly)}: lon({lon_min}, {lon_max}), lat({lat_min}, {lat_max}) using projection {proj_type}"
        )
        lat_lon_bounds.append((lon_max, lon_min, lat_max, lat_min))
    ds = nc.Dataset(georegion_file)
    # find georegion variable:
    georegion_var = None
    for var_name in ds.variables:
        if "georegion" in var_name.lower():
            georegion_var = ds.variables[var_name]
            break
    if georegion_var is None:
        print(f"-E-: No variable with 'georegion' found in file '{georegion_file}'")
        sys.exit(1)
    # for each determine indeces in georegion lat/lon arrays
    at_least_one_intersection = False
    for i, bounds in enumerate(lat_lon_bounds):
        lon_max, lon_min, lat_max, lat_min = bounds
        # check if (lat_min, lat_max) and (lon_min, lon_max) intersect with georegion lat/lon arrays
        lat_indices = (lat >= lat_min) & (lat <= lat_max)
        lon_indices = (lon >= lon_min) & (lon <= lon_max)
        if not lat_indices.any() or not lon_indices.any():
            if args.verbose:
                print(
                    f"Shape {i+1}/{len(poly)} is outside of the georegion bounding box, defined in file '{georegion_file}'"
                )
            continue
        else:
            # now we need to compute indexes for lat_min, lat_max, lon_min, lon_max
            lat_min_index = np.where(lat >= lat_min)[0][0]
            lat_max_index = np.where(lat <= lat_max)[0][-1]
            lon_min_index = np.where(lon >= lon_min)[0][0]
            lon_max_index = np.where(lon <= lon_max)[0][-1]
            # extract the subarrays from georegeion variables
            georegion_subarray = georegion_var[
                lat_min_index : lat_max_index + 1, lon_min_index : lon_max_index + 1
            ]
            # create a mesh grid of lat/lon for the subarray
            lat_subarray = lat[lat_min_index : lat_max_index + 1]
            lon_subarray = lon[lon_min_index : lon_max_index + 1]
            lon_grid, lat_grid = np.meshgrid(lon_subarray, lat_subarray)
            projector = pyproj.Transformer.from_crs(
            "EPSG:4326", projection_list[i], always_xy=True
        ).transform
            polygon = wkt.loads(poly[i].wkt)
            polygon_projected = transform(projector, polygon)
            lon_mesh_proj, lat_mesh_proj = projector(lon_grid, lat_grid)
            mask = contains_xy(polygon_projected, lon_mesh_proj, lat_mesh_proj)
            # check if any of the georegion_subarray points within the polygon mask are valid (non-zero)
            if np.any(georegion_subarray[mask] != 0):
                if args.verbose:
                    print(
                        f"Shape {i+1}/{len(poly)} intersects the georegion, defined in file '{georegion_file}'"
                    )
                at_least_one_intersection = True
    ds.close()
    if not at_least_one_intersection:
        print(
            f"-W-: None of the shapes intersects the georegion, defined in file '{georegion_file}'"
        )
        sys.exit(110)
    print(f"At least one of the shapes intersects the georegion, defined in file '{georegion_file}'")
    sys.exit(0)