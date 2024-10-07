#! /usr/bin/env python3

import argparse
import sys
from netCDF4 import Dataset as NC
from pyhdf.SD import SD as HDF
from pyhdf.SD import SDC
import numpy as np
from datetime import datetime
from datetime import datetime, timedelta

datetime.isoformat

## Convert an L1A CZCS HDF File into a L1A CZCS NetCDF File
## Reads the contents of the input HDF file and outputs it as a NetCDF
## To generate a L1A file, download 

__version__ = "1.0 2024-05-20"

# Constants
HDF_MODE = SDC.READ
NC_MODE = "w"
NC_FORMAT = "NETCDF4"



# python has no switch, so this dict is used to get the NC equal based on HDF4 types
# NetCDF API sets data types for variabled based on NumPy's type.
# list is not extensive, only includes ones used in CZCS 
GET_NC_TYPES = {
    SDC.UINT8: np.uint8,        # ubyte
    SDC.FLOAT32: np.float32,    # float
    SDC.INT16: np.int16,        # short
    SDC.INT32: np.int32         # int
} 


# Global files to make reading and writing easier across functions
hdfFile = None
ncFile = None
errorAt = "" # name of whatever dataset the program was trying to copy/do



def main():
    print(f"l1aconvert_czcs {__version__}")
    global ncFile
    global hdfFile

    # Create an ArgumentParser object
    parser = argparse.ArgumentParser(description="Given a HDF4 Scientific Dataset (SD), convert it into NetCDF format with the same name.")

    # Define arguments
    parser.add_argument("iFile", type=str, help="HDF4 file")
    parser.add_argument("oFile", nargs='?', type=str, help="Optional Output file name")

    # Parse the command-line arguments
    # if oFile name is not given, then use the filename with .nc extension
    args = parser.parse_args()
    fileName = args.iFile
    oFileName = fileName + ".nc" if args.oFile is None else args.oFile
    
    print(f"\nInput file:\t{fileName}")
    print(f"Output file:\t{oFileName}\n")

    # Opening the file
    try:
        hdfFile = HDF(fileName, HDF_MODE)
        print(f"Opening file:\t{fileName}")
    except:
        print(f"\n-E- Error opening file named: {fileName}.\n Make sure the filetype is hdf4.\n")
        exit()
    
    # on successful open, make a netcdf file to be written into
    #ncFileName = fileName + ".nc" if oFileName == "" else oFileName 
    ncFile = NC(oFileName, NC_MODE, NC_FORMAT)

    copyGlobalAttributes(fileName, oFileName)
    copyDimensions()
    copyDatasets()
    closeFiles()

    print("Finished!")



# Runs at the end of the program, closing both files.
def closeFiles():
    global hdfFile
    global ncFile

    print("Closing HDF4 File...")
    hdfFile.end()
    print("Closing NetCDF File...")
    ncFile.close()



# Convert HDF naming to NetCDF naming
# ie: "Data Center" to "data_center"
def convertToNetcdfName(string:str):

    # hdf has 3 and 6 as names, but in NetCDF, give it a proper name
    if (string == "3"):
        return "vector_elements"
    elif (string == "6"):
        return "calibration_elements"
    return string.lower().replace(" ", "_")



# Copy global attributes from the HDF4 file into the new NetCDF file 
def copyGlobalAttributes(iFile:str, oFile:str):
    global hdfFile
    global ncFile
    global errorAt
    
    # Per CF Conventions, put unit descriptions and values in proper places.
    # Items in this list are already in the dimensions or described in the variables.
    # Some items are not needed (ie. Satrt Year, Day, etc. bc of time_coverage_start)
    IGNORE = [                              # already in:
        "Pixel per Scan Line",                  # dimensions
        "Number of Scan Lines",                 # dimensions
        "Number of Pixel Control Points",       # dimensions
        "Number of Scan Control Points",        # dimensions
        "Start Year",                           # made time_coverate_start and
        "Start Day",                            # time_converate_end
        "Start Millisec", 
        "End Year", 
        "End Day", 
        "End Millisec",

        # not sure about these because variable "slope" has this information
        # for each scan already. 
        "Calibration Slope",    
        "Calibration Intercept",
    ]


    try:   
        print("Copying global attributes...")

        gloablAttr = hdfFile.attributes()

        for name, val in gloablAttr.items():
            errorAt = name
            ## Skips start and end Year, Day and Milisec b/c they will be made into an isodate
            if (name in IGNORE):
                continue
            
            valType = type(val)
            
            # strings
            if (isinstance(val, str)):
                ncFile.setncattr(convertToNetcdfName(name), val)
            
            # numbers
            else:
                val = np.float32(val) if isinstance(val, float) else np.int32(val)
                ncFile.setncattr(convertToNetcdfName(name), val)
        

        # after copying all the global attrs, make the isodate time
        # make isodate for start time
        errorAt = "time_coverage_start"
        year = gloablAttr.get("Start Year")
        day = gloablAttr.get("Start Day")
        msec = gloablAttr.get("Start Millisec")

        start_of_year = datetime(year, 1, 1)
        time = start_of_year + timedelta(days=(day-1), seconds=(msec/1000))
        
        ncFile.setncattr("time_coverage_start", str(time.isoformat()))

        # make isodate for end time
        errorAt = "time_coverage_end"
        year = gloablAttr.get("End Year")
        day = gloablAttr.get("End Day")
        msec = gloablAttr.get("End Millisec")

        start_of_year = datetime(year, 1, 1)
        time = start_of_year + timedelta(days=(day-1), seconds=(msec/1000))
        
        ncFile.setncattr("time_coverage_end", str(time.isoformat()))

        # netcdf has history global attribute, so add one here:
        # to convert a file, it's static except for the input and output files
        ncFile.setncattr("history", f"python3 l1acovert_czcs {iFile} {oFile}")

        # add converter version into the global attributes
        ncFile.setncattr("l1aconvert_czcs_version", __version__)


    except:
        print(f"-E- Error copying global attributes. Was processing <{errorAt}> from HDF4 when error was caught.")
        exit()
    errorAt = "" # reset 



# Open band1 dataset to copy number of scan lines and pixles per scan line
# Open cntl_pt_cols for number of pixel control points
# other dimensions are constant
def copyDimensions():
    global hdfFile
    global ncFile
    global errorAt

    try:
        print("Copying dimensions...")

        # these 2 dimensions were not named in HDF4, so giving it a name in NetCDF
        # datasets that uses these 2 are named 6 and 3.
        ncFile.createDimension("calibration_elements", 6)
        ncFile.createDimension("vector_elements", 3)

        # constant dims
        ncFile.createDimension("num_qual", 5)
        ncFile.createDimension("number_of_bands", 6)

        # Copy over number of scan lines and pixel per scan line data
        errorAt = "Number of Scan Lines & Pixel per Scan Line"
        currSet = hdfFile.select("band1")
        dims = currSet.dimensions()
        ncFile.createDimension("number_of_scan_lines", dims["Number of Scan Lines"])
        ncFile.createDimension("pixels_per_scan_line", dims["Pixels per Scan Line"])
        currSet.endaccess()
        
        # copy over number of pixel control points
        errorAt = "Number of Pixel Control Points"
        currSet = hdfFile.select("cntl_pt_cols")
        dims = currSet.dimensions()
        ncFile.createDimension("number_of_pixel_control_points", dims["Number of Pixel Control Points"])
        currSet.endaccess()
    
    except:
        print(f"-E- Error copying dimensions. Was trying to copy the dimension(s) <{errorAt}> when error was caught.")
        exit()



# Given a NetCDF variable, assign the attributes in attrList.
# attrList is the attribute list from the HDF4 dataset and it is copying
# that over to the NetCDF version
#
# When assigning valid range, it is in string format. The dataType parameter
# helps with the slicing and converting it into the right range
def assignNcVarAttr(ncVar, attrDict:dict, dataType):

    for attr in attrDict.keys():
        if (attr == "long_name"):
            ncVar.long_name = attrDict.get(attr)
        if (attr == "units"):
            ncVar.units = attrDict.get(attr)

        # valid range is valid_min and valid_max in netcdf
        if (attr == "valid range" or attr == "valid_range"):

            # valid range is a string tuple, extract the values '(1,3)
            # slice the string to get rid of the Parentheses () 
            validRange = attrDict.get(attr)[1:-2].split(",")
  
            # floats like ranges
            if (dataType == np.float32):
                ncVar.valid_min = np.float32(validRange[0])
                ncVar.valid_max = np.float32(validRange[1])
            elif (dataType == np.uint8):
                ncVar.valid_min = np.uint8(validRange[0])
                ncVar.valid_max = np.uint8(validRange[1])
            elif (dataType == np.int16):
                ncVar.valid_min = np.int16(validRange[0])
                ncVar.valid_max = np.int16(validRange[1])
            else:
                ncVar.valid_min = np.int32(validRange[0])
                ncVar.valid_max = np.int32(validRange[1])



# Retrieves the HDF4 dataset and get all the required information to build
# the NetCDF version
def copyDatasets():
    global hdfFile
    global ncFile
    global errorAt

    try:
        print("Copying datasets/variables...")

        datasetNames = hdfFile.datasets().keys()
        for name in datasetNames:
            errorAt = name
            # hdf4 getting data
            currSet = hdfFile.select(name)
            hdfDataType = currSet.info()[3]; # index 3 gives the datasets data type
            hdfDims = currSet.dimensions().keys() # get dimension names 
            setAttrs = currSet.attributes() # variable attrs like long_name, units, etc.
            data = currSet.get() # retrieves the data 
            hdfDatasetAttr = currSet.attributes()

            # netcdf writing data
            ncDatatype = GET_NC_TYPES.get(hdfDataType)
            ncDims = tuple(map(lambda dim: convertToNetcdfName(dim), hdfDims))
            newVariable = ncFile.createVariable(name, ncDatatype, ncDims)
            newVariable[:] = data


            # netcdf assiging attributes to the variables. NcVar is an object and can be
            # passed by reference, so the function takes care of it
            assignNcVarAttr(newVariable, hdfDatasetAttr, ncDatatype)
            
    except:
        print(f"-E- Error copying datasets/variables. Error occured with HDF4 dataset named <{errorAt}>")
        exit()



if __name__ == '__main__':
    main()