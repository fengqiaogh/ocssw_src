#!/usr/bin/env python3
import datetime
import sys
from shapely import wkt
import pyproj
from shapely.ops import transform
import numpy as np
import geopandas
from shapely import contains_xy
import argparse
import numpy as np
import netCDF4 as nc
import os
from typing import List
from georegion_utils.georegion_utils import determine_projection, get_bounds, read_wkt_file,check_file_type
VERSION = "v.2.1.1"
PROGRAM = "georegion_gen"

running_time = datetime.datetime.now()
RUNNING_TIME = running_time.strftime("%Y-%m-%dT%H%M%S")

def generate_georegion(
    inp_shapefiles: List[str],
    out_netcdf: str,
    stepsize: float,
    bound_accuracy: int,
    **kwargs,
):
    poly = []

    # check file types and apply appropriate reading method
    for inp_shapefile in inp_shapefiles:
        file_type = check_file_type(inp_shapefile)
        if file_type == "ASCII text":
            shape_data = read_wkt_file(inp_shapefile)
            for geom in shape_data.geometry:
                poly.append(geom)
        elif file_type == "Shapefile":
            shape_data = geopandas.read_file(inp_shapefile)
            for geom in shape_data.geometry:
                poly.append(geom)
        else:
            print(
                f"-E-: Unsupported file type '{file_type}' for file '{inp_shapefile}'"
            )
            sys.exit(1)

    projection_list = []
    lat_lon_bounds = []
    global_lat_min = 90.0
    global_lat_max = -90.0
    global_lon_min = 180.0
    global_lon_max = -180.0
    for p in poly:
        proj_type = determine_projection(p.wkt)
        projection_list.append(proj_type)
        lon_max, lon_min, lat_max, lat_min = get_bounds(
            p.wkt, proj_type, resolve=bound_accuracy
        )
        print(
            f"Determined bounds for shape {poly.index(p)+1}/{len(poly)}: lon({lon_min}, {lon_max}), lat({lat_min}, {lat_max}) using projection {proj_type}"
        )
        lat_lon_bounds.append((lon_max, lon_min, lat_max, lat_min))
        global_lat_min = min(global_lat_min, lat_min)
        global_lat_max = max(global_lat_max, lat_max)
        global_lon_min = min(global_lon_min, lon_min)
        global_lon_max = max(global_lon_max, lon_max)
    if "longitude_bounds" in kwargs:
        if kwargs["longitude_bounds"]:
            global_lon_min = kwargs["longitude_bounds"][0]
            global_lon_max = kwargs["longitude_bounds"][1]
    if "latitude_bounds" in kwargs:
        if kwargs["latitude_bounds"]:
            global_lat_min = kwargs["latitude_bounds"][0]
            global_lat_max = kwargs["latitude_bounds"][1]
    print(
        f"Determined global bounds: lon({global_lon_min}, {global_lon_max}), lat({global_lat_min}, {global_lat_max})"
    )
    scale_lat = stepsize
    scale_lon = scale_lat
    # if sizes is provided, override stepsize
    if "sizes" in kwargs and not kwargs["sizes"] is None:
        lat_size = kwargs["sizes"][0]
        lon_size = kwargs["sizes"][1]
        scale_lat = (global_lat_max - global_lat_min) / (lat_size - 1)
        scale_lon = (global_lon_max - global_lon_min) / (lon_size - 1)
        print(
            f"Overriding stepsize with provided sizes: lat size {lat_size}, lon size {lon_size}, resulting stepsize lat {scale_lat}, lon {scale_lon}"
        )
    else:
        lat_size = int((global_lat_max - global_lat_min) / scale_lat) + 1
        lon_size = int((global_lon_max - global_lon_min) / scale_lon) + 1
    print(f"Generating georegion with shape ({lat_size}, {lon_size}) along (lat, lon)")
    lat_grid = np.linspace(global_lat_min, global_lat_max, lat_size)
    lon_grid = np.linspace(global_lon_min, global_lon_max, lon_size)
    georegion_mask = np.zeros((lat_size, lon_size), dtype=bool)
    for p, proj_type, lat_lon_bound in zip(poly, projection_list, lat_lon_bounds):
        # ensure that lon/lon bounds are within global bounds
        lat_lon_bound = (
            min(global_lon_max, lat_lon_bound[0]),
            max(global_lon_min, lat_lon_bound[1]),
            min(global_lat_max, lat_lon_bound[2]),
            max(global_lat_min, lat_lon_bound[3]),
        )
        # determine lat/lon indices
        lat_start_idx = int((lat_lon_bound[3] - global_lat_min) / scale_lat)
        lat_end_idx = int((lat_lon_bound[2] - global_lat_min) / scale_lat) + 1
        lon_start_idx = int((lat_lon_bound[1] - global_lon_min) / scale_lon)
        lon_end_idx = int((lat_lon_bound[0] - global_lon_min) / scale_lon) + 1
        print(
            f"Processing shape {poly.index(p)+1}/{len(poly)} with projection {projection_list[poly.index(p)]}: lat indices ({lat_start_idx}, {lat_end_idx}), lon indices ({lon_start_idx}, {lon_end_idx})"
        )
        sub_lon_grid = lon_grid[lon_start_idx:lon_end_idx]
        sub_lat_grid = lat_grid[lat_start_idx:lat_end_idx]
        sub_lon_mesh, sub_lat_mesh = np.meshgrid(sub_lon_grid, sub_lat_grid)
        projector = pyproj.Transformer.from_crs(
            "EPSG:4326", proj_type, always_xy=True
        ).transform
        polygon = wkt.loads(p.wkt)
        polygon_projected = transform(projector, polygon)
        lon_mesh_proj, lat_mesh_proj = projector(sub_lon_mesh, sub_lat_mesh)
        mask = contains_xy(polygon_projected, lon_mesh_proj, lat_mesh_proj)
        georegion_mask[lat_start_idx:lat_end_idx, lon_start_idx:lon_end_idx] |= mask
    georegion_mask = georegion_mask.astype(np.int8)
    print(
        f"Generated georegion with shape {georegion_mask.shape} along (lat,lon). Saving to {out_netcdf}"
    )
    with nc.Dataset(out_netcdf, "w") as georegion:
        georegion.version = VERSION
        georegion.program = PROGRAM
        georegion.date_created = RUNNING_TIME
        if "location" in kwargs and not kwargs["location"] is None:
            georegion.location = kwargs["location"]
        georegion.createDimension("lat", lat_size)
        georegion.createDimension("lon", lon_size)
        georegion_var = georegion.createVariable(
            varname="georegion",
            datatype=np.int8,
            dimensions=("lat", "lon"),
            zlib=True,
            chunksizes=(432, 864),
            fill_value=-128,
        )
        lat_var = georegion.createVariable(
            varname="lat",
            datatype=np.float64,
            dimensions="lat",
            zlib=True,
            chunksizes=(432,),
            fill_value=-32767,
        )
        lon_var = georegion.createVariable(
            varname="lon",
            datatype=np.float64,
            dimensions="lon",
            zlib=True,
            chunksizes=(864,),
            fill_value=-32767,
        )
        georegion_var[:] = georegion_mask
        lat_var[:] = lat_grid
        lon_var[:] = lon_grid

        lon_var.standard_name = "longitude"
        lon_var.long_name = "longitude"
        lon_var.units = "degrees_east"
        lon_var.axis = "X"

        lat_var.standard_name = "latitude"
        lon_var.long_name = "latitude"
        lat_var.units = "degrees_north"
        lat_var.axis = "Y"

        georegion_var.long_name = "geographical_binary_mask"
        georegion_var.description = (
            "Binary mask for masking out a specific earth domain "
        )
        georegion_var.comment = "0 = outside, 1 = inside the region"
        georegion_var.valid_min = 0
        georegion_var.valid_max = 1
        # print global attributes
        georegion.Conventions = "CF-1.8 ACDD-1.3"
        georegion.license = "https://science.nasa.gov/earth-science/earth-science-data/data-information-policy/"
        georegion.naming_authority = "gov.nasa.gsfc.sci.oceandata"
        georegion.keywords_vocabulary = (
            "NASA Global Change Master Directory (GCMD) Science Keywords"
        )
        georegion.standard_name_vocabulary = "CF Standard Name Table v79"
        georegion.creator_name = "NASA/GSFC/OBPG"
        georegion.creator_email = "data@oceancolor.gsfc.nasa.gov"
        georegion.creator_url = "https://oceancolor.gsfc.nasa.gov"
        georegion.project = "Ocean Biology Processing Group"
        georegion.publisher_name = "NASA/GSFC/OB.DAAC"
        georegion.publisher_email = "data@oceancolor.gsfc.nasa.gov"
        georegion.publisher_url = "https://oceancolor.gsfc.nasa.gov"
        georegion.institution = (
            "NASA Goddard Space Flight Center, Ocean Biology Processing Group"
        )
        georegion.product_name = os.path.basename(out_netcdf)
        georegion.history = " ".join(sys.argv)
        georegion.geospatial_lat_min = global_lat_min
        georegion.geospatial_lat_max = global_lat_max
        georegion.geospatial_lon_min = global_lon_min
        georegion.geospatial_lon_max = global_lon_max
        georegion.geospatial_lat_units = "degrees_north"
        georegion.geospatial_lon_units = "degrees_east"
        georegion.geospatial_lat_resolution = scale_lat
        georegion.geospatial_lon_resolution = scale_lon
        georegion.grid_mapping_name = "latitude_longitude"
    print("done")


if __name__ == "__main__":
    print(PROGRAM, VERSION)
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--ifiles",
        help="input shape files or ASCII files with WKT strings",
        type=str,
        required=True,
        nargs="+",
    )
    parser.add_argument(
        "--ofile", help="output georegion file", type=str, required=True
    )
    parser.add_argument(
        "--stepsize",
        help="stepsize along latitude/longitude",
        type=float,
        default=1 / 239.999,
    )  # gebco resolution
    parser.add_argument("--pversion", help="pversion", type=str, default="V11")
    # optinonal arguments,  add image size along lat and lon, a list of two integers
    parser.add_argument(
        "--sizes", help="image sizes along lat and lon", type=int, nargs=2
    )
    parser.add_argument(
        "--bound_accuracy", help="accuracy of bounds", type=int, default=1000
    )
    parser.add_argument(
        "--latitude_bounds", help="latitude bounds", type=float, nargs=2
    )
    parser.add_argument(
        "--longitude_bounds", help="longitude bounds", type=float, nargs=2
    )
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    print(
        "You are running {}, version {}. Current time is {}".format(
            PROGRAM, VERSION, running_time
        )
    )
    args = parser.parse_args()
    # "data/iho/iho.shp"  # atlantic ocean.  https://www.marineregions.org/gazetteer.php?p=details&id=1912
    shapefile_path = args.ifiles
    stepsize = args.stepsize
    out_nc = args.ofile  # "data/georegion_atlantic.nc"
    generate_georegion(
        shapefile_path,
        out_nc,
        stepsize,
        bound_accuracy=args.bound_accuracy,
        sizes=args.sizes if args.sizes else None,
        latitude_bounds=args.latitude_bounds if args.latitude_bounds else None,
        longitude_bounds=args.longitude_bounds if args.longitude_bounds else None,
    )
