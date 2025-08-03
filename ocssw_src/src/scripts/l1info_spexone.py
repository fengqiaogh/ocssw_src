#! /usr/bin/env python3.8

# l1info for SpexOne L1B. Returns the following information:
#   1. Start Date and Time
#   2. End Date and Time
#   3. Gring information (lon-lats) for upper left, middle and right 
#       points and lower left, middle and right points in clockwise
#       direction

import argparse
import datetime
from netCDF4 import Dataset as NetCDF

__version__ = "1.1 (2024-04-25)"


# CONSTANTS
MIN_LAT_BOUND = -50.0 
MAX_LAT_BOUND = 50.0

MIN_DELTA_LAT = 5.0


'''
    Given a lon and lat dataset, process the gring of the datasets
    Return an array of tuples (lon,lat) for the left and right side
    of the gring
'''
def getGring(lon_dataset, lat_dataset):

    # array containing tuples to be returned
    leftArr = []
    rightArr = []

    # tracks prev lon and lat values on the left and right
    prevLonLeft = None
    prevLatLeft = None
    prevLonRight = None
    prevLatRight = None

    # go through each line and get left and right points
    for i in range(num_bins):

        # only update prev values when points are added
        updatePrevVals = False

        # left, right points
        leftPoint = (lon_dataset[i][num_spatial-1], lat_dataset[i][num_spatial-1])
        rightPoint = (lon_dataset[i][0], lat_dataset[i][0])

        # individual left-right lon lats
        leftLon = leftPoint[0]
        leftLat = leftPoint[1]
        rightLon = rightPoint[0]
        rightLat = rightPoint[1]

        # first point, record the lon lat values as prev and append the left, middle and 
        # right points
        if i == 0:

            # upper middle point needs to be first element in left array bc of gring clock-wise printing.
            # Left elements are reversed when printing and the last element to print is the upper middle.
            upperMiddlePoint = (lon_dataset[i][(num_spatial-1)//2], lat_dataset[i][(num_spatial-1)//2])
            leftArr.append(upperMiddlePoint) 
            leftArr.append(leftPoint)
            rightArr.append(rightPoint)
            updatePrevVals = True
           

        # last point, add left, middle and right points
        elif i == num_bins-1:
            leftArr.append(leftPoint)
            rightArr.append(rightPoint)

            # after the right side is done printing, it prints the lower middle
            # point and then prints the left points in reverse order
            lowerMiddlePoint = (lon_dataset[i][(num_spatial-1)//2], lat_dataset[i][(num_spatial-1)//2])
            rightArr.append(lowerMiddlePoint) 
            updatePrevVals = True
           
        
        # above 50.0 and below -50.0 lat bounds, get more points
        elif (leftLat >= MAX_LAT_BOUND or leftLat <= MIN_LAT_BOUND or rightLat >= MAX_LAT_BOUND or rightLat <= MIN_LAT_BOUND):

            # check both delta lon and lat >= 10 on both sides
            leftDeltaLon = abs(leftLon - prevLonLeft)
            leftDeltaLat = abs(leftLat - prevLatLeft)
            rightDeltaLon = abs(rightLon - prevLonRight)
            rightDeltaLat = abs(rightLat - prevLatRight)

            # if the delta lat is at least MIN_DELTA_LAT, grab the point 
            if ((leftDeltaLat >= MIN_DELTA_LAT and rightDeltaLat >= MIN_DELTA_LAT) or 
                        (leftDeltaLon >= 5.0 and rightDeltaLon >= MIN_DELTA_LAT)):

                leftArr.append(leftPoint)
                rightArr.append(rightPoint)
                updatePrevVals = True
  
        # in between, check only lat
        elif (leftLat < 40.0 and leftLat > -40.0) or (rightLat < 40.0 and rightLat > -40.0):
            deltaLat = abs(leftLat - prevLatLeft)
            if (deltaLat >= 20.0):
                leftArr.append(leftPoint)
                rightArr.append(rightPoint)
                updatePrevVals = True
        
        # update prev values
        if updatePrevVals:
            prevLonLeft = leftLon
            prevLatLeft = leftLat
            prevLonRight = rightLon
            prevLatRight = rightLat
    
    return rightArr, leftArr


# Prints the gring info for both lat and lon.
def printGring(rightArr, leftArr, isLon:bool):
    
    # flag depending on if it's printing lon or lat gring
    flag = 0 if isLon else 1

    # print upper right
    if (isLon):
        print("gringpointlongitude=", end="")
    else:
        print("gringpointlatitude=", end="")


    arrSize = len(rightArr)

    # print all right elements and lower middle
    for i in range(arrSize):
        print(rightArr[i][flag], end=",")

    # print all left elements and upper middle
    for i in reversed(range(arrSize)):
        if i == 0:
            print(leftArr[i][flag], end="\n")
        else:
            print(leftArr[i][flag], end=",")
    
            

def main():
    print("l1info_spexone", __version__)

    # parse command line info
    parser = argparse.ArgumentParser(description="Given a SpexOne L1B file, return its time and gring information")
    parser.add_argument("iFile", type=str, help="SpexOne L1B NetCDF file")
    args = parser.parse_args()

    # open netcdf file
    ncFile = NetCDF(args.iFile, "r+", format="NETCDF4")

    # shape of lon and lat dataset. Global for other function access
    global num_bins       # 1st dim
    global num_spatial    # 2nd dim
    num_bins = len(ncFile.dimensions["bins_along_track"])
    num_spatial = len(ncFile.dimensions["spatial_samples_per_image"])

    # get start
    startTime = ncFile.__dict__["time_coverage_start"]
    endDatetime = datetime.datetime(int(startTime[:4]), int(startTime[5:7]), int(startTime[8:10]))

    # get last bin start time, this will be the end time
    var = ncFile.groups["BIN_ATTRIBUTES"]["image_time"]
    # ignore the valid_max for this variable
    var.set_auto_mask(False)
    lastBinSeconds = float(var[num_bins-1])
    endTime = datetime.timedelta(seconds=lastBinSeconds)
    endDatetime += endTime
    endTimeStr = endDatetime.isoformat()
    if(len(endTimeStr) > 23):
        endTimeStr = endTimeStr[:23]

    # get lon and lat dataset
    lon_dataset = ncFile.groups["GEOLOCATION_DATA"]["longitude"]
    lat_dataset = ncFile.groups["GEOLOCATION_DATA"]["latitude"]

    # get gring items for printing
    rightArr, leftArr = getGring(lon_dataset, lat_dataset)


    ### date and time ###
    print(f"Start_Date={startTime[:10]}")
    print(f"Start_Time={startTime[11:23]}")
    print(f"End_Date={endTimeStr[:10]}")
    print(f"End_Time={endTimeStr[11:]}")
    

    ### gringpoint ###
    printGring(rightArr, leftArr, True)
    printGring(rightArr, leftArr, False)
    

    ### gring sequence ###
    print("gringpointsequence=", end="")

    # 6 counts for the upper and lower indicies not in the center arr
    for i in range(len(leftArr) * 2):
        if i == (len(leftArr) * 2)-1: # last index, dont end with ,
            print(str(i+1))
        else:
            print(str(i+1), end=",")


if __name__ == "__main__":
    main()
