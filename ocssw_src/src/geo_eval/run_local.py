#!/usr/bin/env python3
import os
import csv
import numpy as np
from geo_eval import get_lat_lon_shift, read_file_l1b
from scipy.spatial import ConvexHull
from shapely.geometry.polygon import Polygon
from shapely.geometry import Point
import sys
import argparse


def get_polygon(lat, lon):
    coord = np.column_stack((lat.flatten(), lon.flatten()))
    hull = ConvexHull(coord)
    points = []
    for lat_p, lon_p in zip(coord[hull.vertices, 0], coord[hull.vertices, 1]):
        points.append((lat_p, lon_p))
    points.append((coord[hull.vertices, 0][0], coord[hull.vertices, 1][0]))
    polygon = Polygon(points)
    return polygon


def get_matching_chips(feature_file, polygon):

    feature_names = []
    feature_lats = []
    feature_lons = []

    with open(feature_file, mode='r') as file:
        csvFile = csv.DictReader(file)
        for lines in csvFile:
            feature_lat = float(lines["featureLat"])
            feature_lon = float(lines["featureLon"])
            south = float(lines["south"])
            north = float(lines["north"])
            west = float(lines["west"])
            east = float(lines["east"])
            point = Point(feature_lat, feature_lon)
            if polygon.contains(point):
                p_0 = (north, west)
                p_1 = (north, east)
                p_2 = (south, east)
                p_3 = (south, west)
                polygon_chip = Polygon([p_0, p_1, p_2, p_3, p_0])
                if polygon.contains(polygon_chip):
                    feature_names.append(lines["filename"])
                    feature_lats.append(feature_lat)
                    feature_lons.append(feature_lon)
    return feature_names


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--ifile", help="input L1B netcdf file", type=str, required=True)
    parser.add_argument(
        "--idir", help="directory with chips", type=str, required=True)
    parser.add_argument(
        "--odir", help="output directory", type=str, required=True)
    parser.add_argument('--debug', action='store_true')
    parser.add_argument('--no-debug', dest='debug', action='store_false')
    parser.set_defaults(debug=False)
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    args = parser.parse_args()
    l1b_path = args.ifile
    lat, lon, _, _, _ = read_file_l1b(l1b_path)
    polygon = get_polygon(lat, lon)
    dir_chips = args.idir
    feature_file = os.path.join(dir_chips, "chipindex.csv")
    debug_mode = False
    feature_names = get_matching_chips(feature_file, polygon)
    feature_names = [os.path.join(dir_chips, x) for x in feature_names]
    outfile = os.path.join(args.odir, f"gcpm-{os.path.basename(l1b_path)}.csv")
    if args.debug:
        with open("chip_lists.txt","w") as ilist:
             ilist.write('\n'.join(feature_names))
    ret_code = get_lat_lon_shift(
        l1b_path, feature_names, outfile, debug_mode=args.debug)
    sys.exit(ret_code)
