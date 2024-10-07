#!/usr/bin/env python3

import argparse
import numpy as np
import sys

def centroidnp(points):
    lon_coords = [p[0] for p in points]
    lat_coords = [p[1] for p in points]
    length = len(points)
    centroid_lon = sum(lon_coords)/length
    centroid_lat = sum(lat_coords)/length
    return [centroid_lon, centroid_lat]

def monotonic(x):
    dx = np.diff(x)
    return np.all(dx <= 0) or np.all(dx >= 0)

__version__ = "1.0-20211104"

parser = argparse.ArgumentParser(prog='sort_gring',formatter_class=argparse.RawTextHelpFormatter,
        description='This script sorts GRing points to ensure counter-clockwise order',epilog='''
EXIT Status:
0   : All is well in the world
1   : Dunno, something horrible occurred
110 : No reordering necessary''')
parser.add_argument('--version', action='version', version='%(prog)s ' + __version__)
parser.add_argument('-i','--ifile', type=str, metavar="FILE",default=None, help="name of file with gring data from l1info")
parser.add_argument('--lon', type=str, default=None, help='list of gring longitudes')
parser.add_argument('--lat', type=str, default=None, help="list of gring latitudes")
parser.add_argument('--verbose', '-v', action='store_true', default=False)

args = parser.parse_args()

gringpointlongitude = None
gringpointlatitude = None

if args.ifile:
    with open(args.ifile) as ifile:
        lines = ifile.readlines()
        lines = filter(str.rstrip, lines)
    for line in lines:    
        if "=" in line:
            key, value =  line.split('=')
            values = np.array(value.split(',')).astype(float)
            if 'longitude' in key:
                gringpointlongitude = np.squeeze(np.radians(np.array([values])))
            if 'latitude' in key:
                gringpointlatitude = np.squeeze(np.radians(np.array([values])))

if args.lon:
    value = args.lon
    values = np.array(value.split(',')).astype(float)
    gringpointlongitude = np.squeeze(np.radians(np.array([values])))

if args.lat:
    value = args.lat
    values = np.array(value.split(',')).astype(float)
    gringpointlatitude = np.squeeze(np.radians(np.array([values])))

if (gringpointlongitude is not None) and (gringpointlatitude is not None):
    if gringpointlongitude.shape != gringpointlatitude.shape:
        sys.exit("lat and lon elements of different size")

    # make a list of tuples for the points
    points = list(zip(gringpointlongitude, gringpointlatitude))
    # find the middle of the polygon
    center = centroidnp(points)
    # get the distance from each point to the center
    distance = np.atleast_1d(points) - np.atleast_1d(center)
    # compute an angle from these
    angles = np.arctan2(distance[:,0],distance[:,1])
    # convert the angles into degrees in the range of 0-360
    deg_angles = np.degrees(angles)
    idx = np.where(deg_angles < 0)
    deg_angles[idx] += 360
    # sort the angles and reverse since we want counter-clockwise
    gringsequence = np.argsort(deg_angles, axis=0)[::-1]

    if monotonic(gringsequence):
        if args.verbose:
            print("already ordered appropriately")
        sys.exit(110)
    else:
        sorted_lon = np.array2string(np.take_along_axis(np.degrees(gringpointlongitude), gringsequence, axis=0), precision=5, separator=',')
        sorted_lat = np.array2string(np.take_along_axis(np.degrees(gringpointlatitude), gringsequence, axis=0), precision=5, separator=',')
        if args.verbose:
            print("new world order: {}".format(gringsequence))

        print("gringpointlongitude={}".format(sorted_lon[1:-1]))
        print("gringpointlatitude={}".format(sorted_lat[1:-1]))
else:
    parser.print_help()
