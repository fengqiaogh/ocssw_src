#!/usr/bin/env python3
import datetime
import geopandas
import rasterio
import rasterio.mask
import argparse
import numpy as np
import netCDF4 as nc
import sys
import os
from typing import List
VERSION = "v.1.0.2"
PROGRAM = "georegion_gen"

running_time = datetime.datetime.now()
RUNNING_TIME = running_time.strftime("%Y-%m-%dT%H%M%S")


def read_wkt_file(path_wkt: str):
    with open(path_wkt) as f:
        wkts = f.read().splitlines()
    return geopandas.GeoSeries.from_wkt(wkts)


def generate_georegion(inp_shapefiles: List[str], out_netcdf: str, stepsize: float, **kwargs):
    poly = []
    if "wkt" in kwargs:
        if kwargs["wkt"]:
            for inp_shapefile in inp_shapefiles:
                shape_data = read_wkt_file(inp_shapefile)
                for geom in shape_data.geometry:
                    poly.append(geom)
        else:
            for inp_shapefile in inp_shapefiles:
                shape_data = geopandas.read_file(inp_shapefile)
                for geom in shape_data.geometry:
                    poly.append(geom)
    poly = geopandas.GeoSeries(poly)
    west_lon = np.min(poly.bounds["minx"][:])
    east_lon = np.max(poly.bounds["maxx"][:])
    north_lat = np.max(poly.bounds["maxy"][:])
    south_lat = np.min(poly.bounds["miny"][:])
    if "global_set" in kwargs:
        if kwargs["global_set"]:
            west_lon = -180.0
            east_lon = 180.0
            north_lat = 90.0
            south_lat = -90.0
    scale_y = stepsize
    scale_x = scale_y
    xsize = int((east_lon - west_lon) / scale_x) + 1
    ysize = int((north_lat - south_lat) / scale_y) + 1
    shift_x_array = np.linspace(west_lon, east_lon, xsize)
    shift_y_array = np.linspace(north_lat, south_lat, ysize)
    trs = rasterio.transform.from_origin(west_lon, north_lat, scale_x, scale_y)
    img = rasterio.features.rasterize(
        poly, transform=trs, out_shape=(ysize, xsize))
    print(
        f"Generated georegion with shape {img.shape} along (lat,lon). Saving to {out_netcdf}")
    with nc.Dataset(out_netcdf, "w") as georegion:
        georegion.version = VERSION
        georegion.program = PROGRAM
        georegion.date_created = RUNNING_TIME
        if "location" in kwargs and not kwargs["location"] is None:
            georegion.location = kwargs["location"]
        georegion.createDimension('lat', ysize)
        georegion.createDimension('lon', xsize)
        georegion_var = georegion.createVariable(varname='georegion', datatype=np.int8, dimensions=('lat', 'lon'), zlib=True,
                                             chunksizes=(432, 864), fill_value=-128)
        lat_var = georegion.createVariable(
            varname='lat', datatype=np.float64, dimensions='lat', zlib=True, chunksizes=(432,), fill_value=-32767)
        lon_var = georegion.createVariable(
            varname='lon', datatype=np.float64, dimensions='lon', zlib=True, chunksizes=(864,), fill_value=-32767)
        georegion_var[:] = img
        lat_var[:] = shift_y_array
        lon_var[:] = shift_x_array

        lon_var.standard_name = "longitude"
        lon_var.long_name = "longitude"
        lon_var.units = "degrees_east"
        lon_var.axis = "X"

        lat_var.standard_name = "latitude"
        lon_var.long_name = "latitude"
        lat_var.units = "degrees_north"
        lat_var.axis = "Y"

        georegion_var.long_name = "geographical_binary_mask"
        georegion_var.description = "Binary mask for masking out a specific earth domain "
        georegion_var.comment = "0 = outside, 1 = inside the region"
        georegion_var.valid_min = 0
        georegion_var.valid_max = 1
        # print global attributes
        georegion.Conventions = "CF-1.8 ACDD-1.3"
        georegion.license = "https://science.nasa.gov/earth-science/earth-science-data/data-information-policy/"
        georegion.naming_authority = "gov.nasa.gsfc.sci.oceandata"
        georegion.keywords_vocabulary = "NASA Global Change Master Directory (GCMD) Science Keywords"
        georegion.standard_name_vocabulary = "CF Standard Name Table v79"
        georegion.creator_name = "NASA/GSFC/OBPG"
        georegion.creator_email = "data@oceancolor.gsfc.nasa.gov"
        georegion.creator_url = "https://oceancolor.gsfc.nasa.gov"
        georegion.project = "Ocean Biology Processing Group"
        georegion.publisher_name = "NASA/GSFC/OB.DAAC"
        georegion.publisher_email = "data@oceancolor.gsfc.nasa.gov"
        georegion.publisher_url = "https://oceancolor.gsfc.nasa.gov"
        georegion.institution = "NASA Goddard Space Flight Center, Ocean Biology Processing Group"
        georegion.product_name = os.path.basename(out_netcdf)
        georegion.history = " ".join(sys.argv)
        georegion.geospatial_lat_min = south_lat
        georegion.geospatial_lat_max = north_lat
        georegion.geospatial_lon_min = west_lon
        georegion.geospatial_lon_max = east_lon
        georegion.geospatial_lat_units = "degrees_north"
        georegion.geospatial_lon_units = "degrees_east"
        georegion.geospatial_lat_resolution = scale_y
        georegion.geospatial_lon_resolution = scale_x
        georegion.grid_mapping_name = "latitude_longitude"
    print("done")


if __name__ == "__main__":
    print(PROGRAM, VERSION)
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--ifiles", help="input shape files", type=str, required=True, nargs='+')
    parser.add_argument(
        "--ofile", help="output georegion file", type=str, required=True)
    parser.add_argument(
        "--stepsize", help="stepsize along latitude", type=float, default=1.0 / 16.0 / 16.0)
    # parser.add_argument(
    #     "--location", help="name of the geographic area", type=str)
    parser.add_argument(
        "--pversion", help="pversion", type=str, default="V11")
    parser.add_argument(
        "--global_set", help="global coverage", type=bool, default=False)
    parser.add_argument(
        "--wkt", help="indicate that input files are ASCII files which contains WKT strings", type=bool, default=False)
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    print("You are running {}, version {}. Current time is {}".format(
        PROGRAM, VERSION, running_time))
    args = parser.parse_args()
    # "data/iho/iho.shp"  # atlantic ocean.  https://www.marineregions.org/gazetteer.php?p=details&id=1912
    shapefile_path = args.ifiles
    stepsize = args.stepsize
    out_nc = args.ofile  # "data/georegion_atlantic.nc"
    generate_georegion(shapefile_path, out_nc, stepsize,
                     global_set=args.global_set, wkt=args.wkt)
