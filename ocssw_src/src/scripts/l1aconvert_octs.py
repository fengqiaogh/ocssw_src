#! /usr/bin/env python3

# Pyhdf use version 0.10.5. 
# Version 0.11.3 has a bug that does not let you read a char dataset.
# The fix was addressed Aug 2023, but the lastest pip package (as of Nov 2023)
# is still Jun 25, 2023. If a new version is released and it is > Aug 2023,
# use that one bc the fix should be merged.
# Fix commit: https://github.com/fhs/pyhdf/commit/e8f87adb3f308fc6dc6105aa631e2c659da8f4e6

import argparse
from netCDF4 import Dataset as NC
from pyhdf.SD import SD as HDF
from pyhdf.SD import SDC
import numpy as np
from datetime import datetime
from datetime import datetime, timedelta

datetime.isoformat

__version__ = "1.0 2024-05-21"

## Convert an L1A OCTS HDF4 File into a L1A NetCDF File
## Reads the contents of the input HDF4 file and outputs it as a NetCDF
## Naming conventions will be changed. ie: Start Year --> start_year

# Constants used to open, read and write the files
HDF_READ = SDC.READ
NC_WRITE = "w"
NC_FORMAT = "NETCDF4"

# Do not copy these to NetCDF. 
# year, day and milisecs will be made into isodate
SKIP_ATTRIBUTES = (
    "Start Time",        # dont need these because of time_coverage_start and end 
    "Start Year", 
    "Start Day", 
    "Start Millisec", 
    "End Time",
    "End Year", 
    "End Day", 
    "End Millisec",
)

# python has no switch, so this dict is used to get the NetCDF type equivalent based on HDF4 types.
# NetCDF API sets data types for variabled based on NumPy's type.
# list may not extensive, only includes ones used in **OCTS**
# key   = HDF4 types
# value = NetCDF types
GET_NC_TYPES = {
    SDC.UINT8: np.uint8,        # ubyte
    SDC.FLOAT32: np.float32,    # float
    SDC.FLOAT64: np.float64,    # double
    SDC.INT16: np.int16,        # short
    SDC.UINT16: np.uint16,      # ushort 
    SDC.UINT32: np.uint32,      # uint
    SDC.INT32: np.int32,        # int,
    SDC.CHAR: 'S1'              # string stored as 8bit char array
} 


# Global files to make reading and writing easier across functions
hdfFile = None
ncFile = None

# name of whatever dataset the program was trying to copy/do
# if error occurs, we'll know what to look at
errorAt = "" 


def main():
    print(f"l1aconvert_octs {__version__}")

    # get global reference to the file variables
    global ncFile
    global hdfFile
    # Create an ArgumentParser object
    parser = argparse.ArgumentParser(description='''
                                     Given an HDF4 Scientific Dataset (SD), convert it
                                     into a NetCDF (NC) file format with the same name.
                                     If the optional argument [oFile] name is given, 
                                     output as the oFile name.
                                     ''')

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

    # Opening the hdf4 file for reading
    try:
        hdfFile = HDF(fileName, HDF_READ)
        print(f"Opening file:\t{fileName}")
    except:
        print(f"\n-E- Error opening file named: {fileName}.\n Make sure the filetype is hdf4.\n")
        exit()

    # on successful open, make a netcdf file to be written into
    ncFile = NC(oFileName, NC_WRITE, NC_FORMAT)

    copyGlobalAttributes(fileName, oFileName)
    copyDimensions()
    copyDatasets()
    closeFiles()
    print("Finished!")




# Convert HDF naming to NetCDF naming
# ie: "Data Center" to "data_center"
def convertToNetcdfName(string:str):
    return string.lower().replace(" ", "_")




# Copy global attributes from the HDF4 file into the new NetCDF file 
def copyGlobalAttributes(iFile:str, oFile:str):
    global hdfFile
    global ncFile
    global errorAt

    try:   
        print("Copying global attributes...")

        gloablAttr = hdfFile.attributes()

        for name, val in gloablAttr.items():
            errorAt = name # if error occurs, we know which global attr is causing it
            # # Skips start and end Year, Day and Milisec b/c they will be made into an isodate
            if (name in SKIP_ATTRIBUTES):
                continue
            
            # handle different data types when copying
            # if string, then do nothing and write normally 

            # float
            if isinstance(val, float):
                val = np.float32(val)
            elif isinstance(val, int):
                val = np.int32(val)
            # a list of ints, without np.int32 type, netcdf sets it as long long
            elif isinstance(val, list) and isinstance(val[0], int):
                val = np.array(val, dtype=np.int32)


            ## ...Add more but these are the only 2 instances I see in OCTS global attributes
            ncFile.setncattr(convertToNetcdfName(name), val)

        
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

        # add converter version
        ncFile.setncattr("l1aconvert_octs_version", __version__)
        
        # netcdf has history global attribute, so add one here:
        # to convert a file, it's static except for the input and output files
        ncFile.setncattr("history", f"python3 l1acovert_octs {iFile} {oFile}")


    except:
        print(f"-E- Error copying global attributes. Was processing <{errorAt}> from HDF4 when error was caught.")
        exit()
    errorAt = "" # reset 


# Given a NetCDF variable, assign the attributes in attrList.
# attrList is the attribute list from the HDF4 dataset and it is copying
# that over to the NetCDF version
def assignNcVarAttr(ncVar, attrDict:dict, dataType):
    for attr in attrDict.keys():
        if (attr == "long_name" or attr == "long name"):
            ncVar.long_name = attrDict.get(attr)
        if (attr == "units"):
            ncVar.units = attrDict.get(attr)

        # valid range is valid_min and valid_max in netcdf
        if (attr == "valid_range" or attr == "valid range"):

            # get the valid range
            validRange = attrDict.get(attr)

            # valid range can be a string tuple, extract the values '(1,3)
            # slice the string to get rid of the Parentheses ()
            # otherwise, it is already a list
            if isinstance(validRange, str):
                validRange = attrDict.get(attr)[1:-2].split(",")
  
            # Based on the nctype, make sure the value is in numpy type
            npValidMin = None
            npValidMax = None
            if (dataType == np.float32): # float
                npValidMin = np.float32(validRange[0]) 
                npValidMax = np.float32(validRange[1])

            elif (dataType == np.float64): # double
                npValidMin = np.float64(validRange[0]) 
                npValidMax = np.float64(validRange[1])

            elif (dataType == np.uint16): # ushort 
                npValidMin = np.uint16(validRange[0]) 
                npValidMax = np.uint16(validRange[1])

            elif (dataType == np.uint32): # uint
                npValidMin = np.uint32(validRange[0])
                npValidMax = np.uint32(validRange[1])

            elif (dataType == np.int16): # short 
                npValidMin = np.int16(validRange[0]) 
                npValidMax = np.int16(validRange[1])

            else: # int
                npValidMin = np.int32(validRange[0])
                npValidMax = np.int32(validRange[1])
            
            # assign it to valid min and max
            ncVar.valid_min = npValidMin
            ncVar.valid_max = npValidMax


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
            data = currSet.get() # retrieves the data 
            hdfDatasetAttr = currSet.attributes() # variable attrs like long_name, units, etc.

            # netcdf writing data

            # get nc data type based on the hdf type
            ncDatatype = GET_NC_TYPES.get(hdfDataType)
            # get the dimentions to be set as a tuple
            ncDims = tuple(map(lambda dim: convertToNetcdfName(dim), hdfDims))
            
            # create the variable and then assign the data
            newVariable = ncFile.createVariable(name, ncDatatype, ncDims)
            newVariable[:] = data # slice entire variable and assign it 
    
            # netcdf assiging attributes to the variables. NcVar is an object and can be
            # passed by reference, so the function takes care of it
            assignNcVarAttr(newVariable, hdfDatasetAttr, ncDatatype)
            
    except Exception as e:
        print(f"-E- Error copying datasets/variables. Error occured with HDF4 dataset named <{errorAt}>")
        print(f"Reason: {e}")
        exit()


# Runs at the end of the program, closing both files.
def closeFiles():
    global hdfFile
    global ncFile

    print("Closing HDF4 File...")
    hdfFile.end()
    print("Closing NetCDF File...")
    ncFile.close()



# Given a list of variables, extract the dimention data from it 
# pyhdf lib does not have a way to get only the dims, you need to 
# go through the variables and get it
def copyDimensions():
    global hdfFile
    global ncFile
    global errorAt

    errorAt = "Copying Dimensions"
    # create a set so that when parsing all the variables with their dims,
    # dims do not get repeated in the set
    dimSet = set()

    try:
        errorAt = "extract dimension details from HDF file"
        # get key, value pair of all the HDF
        hdfDatasets = hdfFile.datasets()
        for var in hdfDatasets.keys():
            errorAt = f"extract {var} from HDF file"
            varInfo = hdfDatasets[var] # returns something like: (('bands',), (8,), 22, 0
            dimNames = varInfo[0]   # tuple of dim names for var
            dimVals = varInfo[1]    # tuple of dim values for the dims

            # go though each dim and save it into the set
            for i in range(len(dimNames)):
                name = dimNames[i]
                val = dimVals[i]
                newDim = (name, val)
                dimSet.add(newDim)
        
        # after getting all the dims, go through and set them in the NetCDF file
        for dim in dimSet:
            errorAt = f"set {dim} in the NetCDF file"
            name = dim[0]
            val = dim[1]
            ncFile.createDimension(name, val)

    except:
        print(f"-E- Error copying dimensions. Was trying to {errorAt} when error was caught.")
        exit()



if __name__ == '__main__':
    main()