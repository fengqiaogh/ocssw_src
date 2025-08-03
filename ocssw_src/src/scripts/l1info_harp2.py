#! /usr/bin/env python3

# l1info for HARP2 L1B. Returns the following information:
#   1.  Start Date and Time
#   2.  End Date and Time
#   3.  Gring information (lon-lats) that encloses all possible lon,lat pairs for the file.
#       Return in clock-wise order.
#   4. Geospatial bounds. The same a Gring, but in counter-clockwise order.
#
#   Return codes:
#           0 if no issues detected
#           100 if bad file type
#           110 if there's bad gring order
#           120 if there's bad geospatial bounds order
#           130 bad file, missing more than 1 corner of the polygon


import argparse
import numpy as np
from pyproj import Geod
from netCDF4 import Dataset as NetCDF

# year, month, day
__version__ = "2.2 (2024-12-17)"


# Class that holds lon and lat points 
class LonLatPoint:

    # constructor
    def __init__(self, lon, lat):
        self.lon = lon
        self.lat = lat


# Class to store the most east and most west point of a single Harp2 Line
class Coordinates:
    def __init__(self, west:LonLatPoint, east:LonLatPoint):
        self.west = west
        self.east = east


# use pyproj's geod class to calculate the area of the lon and lats
# to determine order
def isClockwise(lons, lats):
    geod = Geod(ellps="WGS84")
    area, perimeter = geod.polygon_area_perimeter(lons, lats)
    return area < 0 # negative area indicates clockwise order




# True if the view and line has at least 2 pixels that are not fill value
# False otherwise.
def Has2GoodPixels(view, line):
    global lonData
    global latData
    global numPixels

    goodPixelCount = 0

    for pix in range(numPixels):
        
        # already seen 2, return and stop checking for more
        if (goodPixelCount == 2):
            return True
        if (~lonData[view][line].mask[pix] and ~latData[view][line].mask[pix]):
            goodPixelCount+=1
    
    return False


# True if the line has is not all fill values. False otherwise.
def LineNotAllMaskValues(view, line):
    global lonData
    global latData
    return np.any(~lonData[view][line].mask) and np.any(~latData[view][line].mask)

# Have a pointer at the start and end of the views line.
# Move the start pointer forward if it contains fill values and move the end
# pointer backwards if it contains fill values.
# return the first good line and last good line when both are good
def GetFirstAndLastGoodLinesForView(view):

    global lonData
    global latData
    global numLines

    firstGoodLine = 0
    lastGoodLine = numLines-1

    # make sure they dont cross when searching
    while (firstGoodLine < lastGoodLine):

        # check if the line is not all fill value
        isGoodFirstLine = LineNotAllMaskValues(view, firstGoodLine)
        isGoodLastLine = LineNotAllMaskValues(view, lastGoodLine)

        # both good, then check that it has at least 2 good pixels
        if (isGoodFirstLine and isGoodLastLine):
            
            goodFirstLine = Has2GoodPixels(view, firstGoodLine)
            goodLastLine = Has2GoodPixels(view, lastGoodLine)

            # both good
            if (goodFirstLine and goodLastLine):
                break
            # only first line has at least 2 pixels
            elif (goodFirstLine and not goodLastLine):
                lastGoodLine-= 1
            # only good last line
            elif (not goodFirstLine and goodLastLine):
                firstGoodLine +=1
            else:
                firstGoodLine+=1
                lastGoodLine-=1

        # good first, bad last, move only the last pointer
        elif (isGoodFirstLine and not isGoodLastLine):
            lastGoodLine-=1
        elif (not isGoodFirstLine and isGoodLastLine):
            firstGoodLine+=1
        else:
            firstGoodLine+=1
            lastGoodLine-=1
    
    return firstGoodLine, lastGoodLine


# Determines if there are at least 2 good lines in the view
def HasGoodFirstAndLastLineForView(view):

    global lonData
    global latData
    global numLines

    firstGoodLine = 0
    lastGoodLine = numLines-1

    # make sure they dont cross when searching
    while (firstGoodLine < lastGoodLine):

        # check if the line is not all fill value
        isGoodFirstLine = LineNotAllMaskValues(view, firstGoodLine)
        isGoodLastLine = LineNotAllMaskValues(view, lastGoodLine)

        # both good, then check that it has at least 2 good pixels
        if (isGoodFirstLine and isGoodLastLine):
            
            goodFirstLine = Has2GoodPixels(view, firstGoodLine)
            goodLastLine = Has2GoodPixels(view, lastGoodLine)

            # both good
            if (goodFirstLine and goodLastLine):
                break
            # only first line has at least 2 pixels
            elif (goodFirstLine and not goodLastLine):
                lastGoodLine-= 1
            # only good last line
            elif (not goodFirstLine and goodLastLine):
                firstGoodLine +=1
            else:
                firstGoodLine+=1
                lastGoodLine-=1

        # good first, bad last, move only the last pointer
        elif (isGoodFirstLine and not isGoodLastLine):
            lastGoodLine-=1
        elif (not isGoodFirstLine and isGoodLastLine):
            firstGoodLine+=1
        else:
            firstGoodLine+=1
            lastGoodLine-=1
    
    # case where is there no good lines or when there is only 1 line
    if (firstGoodLine >= lastGoodLine):
        return False
    
    return True

# False if it contains masked values, true otherwise
def CheckPixelForNoMask(view, line, pixel):
    global lonData
    global latData
    return np.any(~lonData[view][line].mask[pixel]) and np.any(~latData[view][line].mask[pixel])


# Given a view number, get the west and east points for a north or south location
def GetWestEastPointsForView(view, isNorth):

    global lonData
    global latData
    global numLines
    global numPixels

    # get first and last good lines of a view. The first line will contain the most
    # southern pixel and the last line will contain the most northern pixel
    firstLine, lastLine = GetFirstAndLastGoodLinesForView(view)

    # if firstLine and lastLine is false, then this view is not good. 
    # try a different view

    try:
        # isNorth == true, then use the last line. Otherwise, use the first line for south
        line = lastLine if isNorth else firstLine

        # assume the first good pixel is most left and last good pixel is most right
        # move the pointers forward or backwards if they are fill values
        firstGoodPixel = 0
        lastGoodPixel = numPixels-1

        # make sure they dont cross when searching for first and last good pixel
        while (firstGoodPixel < lastGoodPixel):

            # get the status of the pixels
            isGoodFirstPixel = CheckPixelForNoMask(view, line, firstGoodPixel)
            isGoodLastPixel = CheckPixelForNoMask(view, line, lastGoodPixel)

            # both good, dont check anymore and exit while loop
            if (isGoodFirstPixel and isGoodLastPixel):
                break

            # good first pixel, but bad last pixel, move only the last pointer
            elif (isGoodFirstPixel and not isGoodLastPixel):
                lastGoodPixel-=1
            
            # bad first pixel, but good last pixel, move only the first pointer
            elif (not isGoodFirstPixel and isGoodLastPixel):
                firstGoodPixel+=1
            
            # both bad, move both pointers
            else:
                firstGoodPixel+=1
                lastGoodPixel-=1

        
        # if no issues, make the lon, lat pairs
        leftPoint = LonLatPoint(lonData[view][line][firstGoodPixel], latData[view][line][firstGoodPixel])
        rightPoint = LonLatPoint(lonData[view][line][lastGoodPixel], latData[view][line][lastGoodPixel])

        # store this line's lon and lat in a coordinate class
        # Harp's left most pixel of a line is east and right most pixel is west
        # Harp's first line of the file is most south and last line is most north
        coordinateObj = Coordinates(east=rightPoint, west=leftPoint)

        return coordinateObj
    
    # line returns None for first and last line. It will cause an exception. Catch it
    # and return Coordinates that are fill values. 
    except Exception as e:
        print(e)
        return Coordinates(east=LonLatPoint(-999, -999), west=LonLatPoint(-999, -999))


# From a list of views, get the lon and lat data for west and east side of the view
# and return a dictionary where the key is the view number and the value is a Coordinate
# object that contains west and east LonLatPoint objects 
def GetLonLatPointForViewsAsDict(viewsList, isNorth):
    viewCoordinates = {}
    for view in viewsList:
        coordinate = GetWestEastPointsForView(view=view, isNorth=isNorth)
        viewCoordinates[view] = coordinate
    
    return viewCoordinates


# Make a gring in clockwise order
# northCoordinates is a dict with view numbers as keys and a Coordinate object as its
# value. Same for southCoordinates
def ConstructGring(northCoordinates, southCoordinates):
    # gring is in clockwise order, do get the lon and lat data in this order based on the
    # static CORNER_VIEWS variables:

    #   (1) NW Views ---------------------------- NE Views (2)
    #   |                                                   |
    #   |                                                   |
    #   (4) reversed(SW View) --------- reversed(SE Views) (3)

    global NW_CORNER_VIEWS
    global NE_CORNER_VIEWS
    global SW_CORNER_VIEWS
    global SE_CORNER_VIEWS

    gringLons = []
    gringLats = []
    
    # (1) North West Points
    for view in NW_CORNER_VIEWS:
        gringLons.append(northCoordinates[view].west.lon)
        gringLats.append(northCoordinates[view].west.lat)
    
    # (2) North East Points
    for view in NE_CORNER_VIEWS:
        gringLons.append(northCoordinates[view].east.lon)
        gringLats.append(northCoordinates[view].east.lat)
    
    # (3) South East Points
    for view in reversed(SE_CORNER_VIEWS):
        gringLons.append(southCoordinates[view].east.lon)
        gringLats.append(southCoordinates[view].east.lat)
    
    # (4) South West Points
    for view in reversed(SW_CORNER_VIEWS):
        gringLons.append(southCoordinates[view].west.lon)
        gringLats.append(southCoordinates[view].west.lat)

    # check that it is in clockwise order
    
    return gringLons, gringLats



# Make Geospatial Bounds in **Counter Clockwise** Order
# northCoordinates is a dict with view numbers as keys and a Coordinate object as its
# value. Same for southCoordinates
def ConstructGeospatialBounds(northCoordinates, southCoordinates):
    # gring is in clockwise order, do get the lon and lat data in this order based on the
    # static CORNER_VIEWS variables:

    #   (1) reversed(NW Views) -------- reversed(NE Views) (4)
    #   |                                                   |
    #   |                                                   |
    #   (2) SW View ----------------------------- SE Views (3)

    global NW_CORNER_VIEWS
    global NE_CORNER_VIEWS
    global SW_CORNER_VIEWS
    global SE_CORNER_VIEWS

    lons = []
    lats = []
    
    # (1) North West Points
    for view in reversed(NW_CORNER_VIEWS):
        lons.append(northCoordinates[view].west.lon)
        lats.append(northCoordinates[view].west.lat)
    
    # (2) South West Points
    for view in SW_CORNER_VIEWS:
        lons.append(southCoordinates[view].west.lon)
        lats.append(southCoordinates[view].west.lat)
    # (3) South East Points
    for view in SE_CORNER_VIEWS:
        lons.append(southCoordinates[view].east.lon)
        lats.append(southCoordinates[view].east.lat)

    # (3) North East Points
    for view in reversed(NE_CORNER_VIEWS):
        lons.append(northCoordinates[view].east.lon)
        lats.append(northCoordinates[view].east.lat)
    
    return lons, lats
    
  
  # Prints the files Gring information. Gring is a non-closed polygon and it prints
# in clockwise direction.
# Point == tuple(lon, lat)

def PrintGring(lons, lats):

    gringLonStr = "gringpointlongitude="
    gringLatStr = "gringpointlatitude="
    # 1 less point because last point is duplicated. gring doesnt need to be closed
    for i in range(len(lons)):

        # dont print bad points
        if (lons[i] == -999 or lats[i] == -999):
            continue
        gringLonStr += f"{lons[i]},"
        gringLatStr += f"{lats[i]},"
    
    # delete the last "," from the final lon and lat string
    gringLonStr = gringLonStr[:-1]
    gringLatStr = gringLatStr[:-1]

    print(gringLonStr)
    print(gringLatStr)
    print(f"gringpointsequence=1,2,3,4")




# Prints the files Geospatial Bounds information. Geospatial Bounds is a closed polygon, where
# the first point is repeated. 
#
# It prints in counter-clockwise direction from the corners: 
# NOTE: upper right repeated to make closed polygon
#     upper right, upper left, lower left, lower right, upper right
#
# Format is different from gring: 
#     1. space to separate lat,lon value pairs
#     2. comma between pairs
# ie:
# lat lon, lat lon, lat lon, ... lat lon

def PrintGeospatialBounds(lons, lats):
    
    geospatialBoundsStr ="geospatial_bounds="
    
    for i in range(len(lons)):
        if (lons[i] == -999 or lats[i] == -999):
            continue
        geospatialBoundsStr += f"{lats[i]} {lons[i]},"

    # remove last comma from string
    geospatialBoundsStr = geospatialBoundsStr[:-1]

    print(geospatialBoundsStr)


# Given all good views that are ordered from West to East, get 4 points from them 
def Get2GoodViews(allViews):
    leftPoint = 0
    rightPoint = len(allViews)-1

    while(leftPoint < rightPoint):
        goodLeftPoint = HasGoodFirstAndLastLineForView(allViews[leftPoint])
        goodRightPoint = HasGoodFirstAndLastLineForView(allViews[rightPoint])

        # both good
        if (goodLeftPoint and goodRightPoint):
            goodLeftViewNum = allViews[leftPoint]
            goodRightViewNum = allViews[rightPoint]
            return [goodLeftViewNum, goodRightViewNum]
        
        # right side is bad
        elif (goodLeftPoint and not goodRightPoint):
            rightPoint-=1
        # left side is bad  
        elif (not goodLeftPoint and goodRightPoint):
            leftPoint += 1
        # both bad
        else:
            leftPoint+=1
            rightPoint-=1
            
    # end of while loop, none is found, this is a bad file without good views to
    # make a polygon
    global emptyCornerCount
    emptyCornerCount+=1
    return []




def main():
    print("l1info_harp2", __version__)
    global status
    status = 0

    # open the netcdf file and extract global variables 
    helpMessage = """Given a Harp2 L1B file, return time, gring and geospatial bounds. \n
    Status Codes:
    0: Everything is fine,
    100: Wrong file type
    110: Bad Geolocation,
    120: Bad Geospacial bounds,
    130: Bad file, missing more than 1 corner of the polygon or wrong file.
    """
   
    # add the help message when users do not give any command line arguments
    parser = argparse.ArgumentParser(
        description=helpMessage, 
        
        # required to keep helpMessage tab and space formatting. Otherwise, it's 1 long string
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("iFile", type=str, help="HARP2 L1B NetCDF file")
    args = parser.parse_args()

    # check that it is a L1B file. L1B for HARP2 contains geolocation data
    # file format is PACE_OCI.date.fileLevel. Slice by "." and index 2 is the file level
    if (args.iFile.split(".")[2].upper() != "L1B"):
        print("Input file is not L1B. Exiting...")
        return 100
    
    ncFile = NetCDF(args.iFile, "r+", format="NETCDF4")

    # get start and end time 
    start_time = ncFile.__dict__["time_coverage_start"]
    end_time = ncFile.__dict__["time_coverage_end"]

    # load lon and lat data. Make it global for easy access to other functions
    global lonData
    global latData
    lonData = ncFile.groups["geolocation_data"]["longitude"][:]
    latData = ncFile.groups["geolocation_data"]["latitude"][:]
    

    # shape of the dataset. make it global so other functions can access
    global numViews, numLines, numPixels 
    numViews = len(ncFile.dimensions["views"])
    numLines = len(ncFile.dimensions["swath_lines"])
    numPixels = len(ncFile.dimensions["swath_pixels"])


    #######
    ## HARP2 views that the lon, lat pairs will come from which will form the 
    ## gring and geospatial bounds. They are static because the lens and the way 
    ## the spacecraft moves will always be the same, no matter what location.
    ##
    ## The view numbers are the best views that gets the corners of HARP2's
    ## lon and lat data. They are listed in order from west to east and for
    ## each consecutive view, the lon-lat point for that view will be more
    ## east than the previous view.
    ##
    ##      NorthWest                                     NorthEast
    ##      (17)--(15)--(13)--(10)---------------(10)--(15)--(18)--(22)
    ##      |                                                          |
    ##      |                                                          |
    ##      |                                                          |
    ##      |                                                          |
    ##      (62)--(63)--(66)--(68)---------------(79)--(63)--(58)--(56)    
    ##      SouthWest                                     SouthEast
    ##
    ##
    #######
    
    ## All good views to construct each corner of the polygon. These views go from west to east.
    ##
    ## ie:  ALL_NW_VIEWS[0] > ALL_NW_VIEWS[1] > ALL_NW_VIEWS[2] ... ALL_NW_VIEWS[n]
    ##      where ALL_NW_VIEWS[0] is more west than ALL_NW_VIEWS[n]
    global ALL_NE_VIEWS, ALL_NE_VIEWS, ALL_SW_VIEWS, ALL_SE_VIEWS
    ALL_NW_VIEWS = [17, 71, 16, 15, 14, 13, 12, 14, 13, 12, 80, 11, 70, 10]
    ALL_NE_VIEWS = [10, 70, 11, 80, 12, 13, 14, 15, 16, 71, 17, 81, 18, 1, 19, 20, 21, 22]
    ALL_SW_VIEWS = [62, 63, 64, 65, 66, 67, 79, 68]
    ALL_SE_VIEWS = [79, 67, 9, 66, 65, 64, 63, 62, 78, 61, 88, 59, 58, 57, 56]

    # increment when getting 2 good views. if it is empty, increment this.
    # if 2 or more, this is a bad file because it will be missing too may corners 
    global emptyCornerCount
    emptyCornerCount = 0

    # Views to grab to make the NorthWest (NW) corners of the gring. Each consecutive view will 
    # have a lon, lat pair than is more east than the previous view. 
    global NW_CORNER_VIEWS
    NW_CORNER_VIEWS = Get2GoodViews(ALL_NW_VIEWS)
    #NW_CORNER_VIEWS = [17, 15, 13, 10]

    # views for the NE corner. Each consecutive view will also have a lon, lat pair
    # than is more east than the previous view  
    global NE_CORNER_VIEWS
    NE_CORNER_VIEWS = Get2GoodViews(ALL_NE_VIEWS)
    #NE_CORNER_VIEWS = [10, 15, 18, 22]
    
    # SW and SE views. Following the same logic as NW and NE
    global SW_CORNER_VIEWS
    SW_CORNER_VIEWS = Get2GoodViews(ALL_SW_VIEWS)
    #SW_CORNER_VIEWS = [62, 64, 67, 68]

    global SE_CORNER_VIEWS
    SE_CORNER_VIEWS = Get2GoodViews(ALL_SE_VIEWS)
    #SE_CORNER_VIEWS = [79, 63, 58, 56]


    # missing corners, it is a bad file
    if emptyCornerCount > 0:
        status = 130


    # for each of the north views, grab the lon and lat data and save it to a map
    # it will be used to construct the gring and geospatial bounds and so that each
    # view's coordinates are not calculated twice
    northViews = list(set(NW_CORNER_VIEWS + NE_CORNER_VIEWS)) # unique views to north side of map
    northViewCoordinates = GetLonLatPointForViewsAsDict(viewsList=northViews, isNorth=True)
    
    # do the same for the south views
    southViews = list(set(SW_CORNER_VIEWS + SE_CORNER_VIEWS)) #
    southViewCoordinates = GetLonLatPointForViewsAsDict(viewsList=southViews, isNorth=False)


    ###### Make Gring ###### 
    # Gring is in clockwise order.
    gringLons, gringLats = ConstructGring(northViewCoordinates, southViewCoordinates)
    

    ###### Make Geospatial Bounds ###### 
    # Geospatial bounds is in counter-clockwise order
    geospatialLons, geospatialLats = ConstructGeospatialBounds(northViewCoordinates, southViewCoordinates)


    #### print the basic time l1info
    print(f"Start_Date={start_time[:10]}")
    print(f"Start_Time={start_time[11:-1]}")
    print(f"End_Date={end_time[:10]}")
    print(f"End_Time={end_time[11:-1]}")

    #### check if gring is in clockwise order
    if not isClockwise(gringLons, gringLats):
        status = 110
    
    #### print gring
    PrintGring(gringLons, gringLats)

    #### print geospatial bounds
    PrintGeospatialBounds(geospatialLons, geospatialLats)

    return status


if __name__ == "__main__":
    main()