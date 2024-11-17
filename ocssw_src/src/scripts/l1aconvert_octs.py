#! /usr/bin/env python3

# Pyhdf use version 0.10.5. 
# Version 0.11.3 has a bug that does not let you read a char dataset.
# The fix was addressed Aug 2023, but the lastest pip package (as of Nov 2023)
# is still Jun 25, 2023. If a new version is released and it is > Aug 2023,
# use that one bc the fix should be merged.
# Fix commit: https://github.com/fhs/pyhdf/commit/e8f87adb3f308fc6dc6105aa631e2c659da8f4e6

import argparse
from netCDF4 import Dataset as NC
import numpy as np
from pathlib import Path
from pyhdf.SD import SD as HDF
from pyhdf.SD import SDC
from netCDF4 import Dataset as NC
from datetime import datetime, timedelta, timezone

datetime.isoformat

__version__ = "1.0 2024-11-01"

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
    "Scene Center Time",
    "Node Crossing Time",
    "Processing Time",
    "Latitude Units",
    "Longitude Units",
    "Scene Center Latitude",
    "Scene Center Longitude",
    "Scene Center Solar Zenith",
    "Upper Left Latitude",
    "Upper Left Longitude",
    "Upper Right Latitude",
    "Upper Right Longitude",
    "Lower Left Latitude",
    "Lower Left Longitude",
    "Lower Right Latitude",
    "Lower Right Longitude",
    "Northernmost Latitude",
    "Southernmost Latitude",
    "Westernmost Longitude",
    "Easternmost Longitude",
    "Start Center Latitude",
    "Start Center Longitude",
    "End Center Latitude",
    "End Center Longitude",
    "Pixels per Scan Line",
    "Number of Scan Lines",
    "Scene Center Scan Line",
    "Processing Log",
    "Filled Scan Lines",
    "FF Missing Frames",
    "Replacement Flag",
    "Processing Control"
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

# Global variables to make reading/writing easier across functions
hdfFile = None
ncFile = None

# Convert HDF naming to NetCDF naming
# ie: "Data Center" to "data_center"
def convertToNetcdfName(hdfname:str):
    return hdfname.lower().replace(" ", "_")

# name of whatever dataset the program was trying to copy/do
# if error occurs, we'll know what to look at
errorAt = "" 

# copy global attributes from the hdf4 file into new NETCDF file
def copyGlobalAttributes(iFile: str, oFile: str):
    global hdfFile
    global ncFile
    global errorAt

    try:
        print("Copying global attributes")
        globalAttr = hdfFile.attributes()

        for name, val in globalAttr.items():
            errorAt = name

            ## Skips start and end Year, Day and Milisec bc they'll be made into an isodate
            if (name in SKIP_ATTRIBUTES): continue
            
            # handle different data types when copying
            # if string, then do nothing and write normally 

            # strings
            if (isinstance(val, str)):
                ncFile.setncattr(convertToNetcdfName(name), val)
            # numbers
            else:
                val = np.float32(val) if isinstance(val, float) else np.int32(val)
                ncFile.setncattr(convertToNetcdfName(name), val)

            ncFile.setncattr(convertToNetcdfName(name), val)

        # convert node crossing time
        # :node_crossing_time = "19961101 15:37:06.301" ;

        errorAt = "node_crossing_time"
        node_crossing_str = globalAttr.get("Node Crossing Time")
        ncd = node_crossing_str[:8]
        nct = node_crossing_str[9:17]
        fsec = node_crossing_str[18:21]
        node_crossing = datetime.strptime(ncd+nct,"%Y%m%d%H:%M:%S")
        ncFile.setncattr("node_crossing_time", str(node_crossing.strftime('%Y-%m-%dT%H:%M:%S.') + fsec[:3] + 'Z'))

        # convert scene center time
        errorAt = "scene_center_time"
        scene_center_str = globalAttr.get("Scene Center Time")
        scd = node_crossing_str[:8]
        sct = node_crossing_str[9:17]
        fsec = node_crossing_str[18:21]
        scene_center = datetime.strptime(scd+sct,"%Y%m%d%H:%M:%S")
        ncFile.setncattr("scene_center_time", str(scene_center.strftime('%Y-%m-%dT%H:%M:%S.') + fsec[:3] + 'Z'))

        # make isodate for start/end time
        errorAt = "time_coverage_start"
        stime_str = globalAttr.get("Start Time")
        scd = stime_str[:8]
        sct = stime_str[9:17]
        fsec = stime_str[18:21]
        stime = datetime.strptime(scd+sct,"%Y%m%d%H:%M:%S")
        ncFile.setncattr("time_coverage_start", str(stime.strftime('%Y-%m-%dT%H:%M:%S.') + fsec[:3] + 'Z'))

        # make isodate for end time
        errorAt = "time_coverage_end"
        etime_str = globalAttr.get("End Time")
        ecd = etime_str[:8]
        ect = etime_str[9:17]
        fsec = etime_str[18:21]
        etime = datetime.strptime(ecd+ect,"%Y%m%d%H:%M:%S")
        ncFile.setncattr("time_coverage_end", str(etime.strftime('%Y-%m-%dT%H:%M:%S.') + fsec[:3] + 'Z'))

        # add converter version
        ncFile.setncattr("l1aconvert_octs_version", __version__)
        
        # netcdf has history global attribute, so add one here:
        # to convert a file, it's static except for the input and output files
        ncFile.setncattr("history", "{}; python3 l1acovert_seawifs {} {}".format(globalAttr.get("Processing Control")[:-1], iFile, oFile))

        #add date_created
        errorAt = "date_created"
        create_date = datetime.now(tz=timezone.utc)

        ncFile.setncattr("date_created",str(create_date.strftime('%Y-%m-%dT%H:%M:%SZ')))

        # some static metadata
        ncFile.setncattr("cdm_data_type","swath")
        ncFile.setncattr("Conventions","CF-1.8, ACDD-1.3")
        ncFile.setncattr("creator_email","data@oceancolor.gsfc.nasa.gov")
        ncFile.setncattr("creator_name","NASA/GSFC/OBPG")
        ncFile.setncattr("creator_url","https://oceancolor.gsfc.nasa.gov")
        ncFile.setncattr("institution","NASA Goddard Space Flight Center, Ocean Biology Processing Group")
        ncFile.setncattr("instrument","OCTS")
        ncFile.setncattr("keywords_vocabulary","NASA Global Change Master Directory (GCMD) Science Keywords")
        ncFile.setncattr("license","https://www.earthdata.nasa.gov/engage/open-data-services-and-software/data-and-information-policy")
        ncFile.setncattr("naming_authority","gov.nasa.gsfc.oceancolor")
        ncFile.setncattr("platform","ADEOS")
        ncFile.setncattr("processing_level","L1A")
        ncFile.setncattr("processing_version","V2")
        ncFile.setncattr("product_name","{}".format(oFile.name))
        ncFile.setncattr("project","Ocean Biology Processing Group")
        ncFile.setncattr("publisher_email","data@oceancolor.gsfc.nasa.gov")
        ncFile.setncattr("publisher_name","NASA/GSFC/OB.DAAC")
        ncFile.setncattr("publisher_url","https://oceancolor.gsfc.nasa.gov")
        ncFile.setncattr("standard_name_vocabulary","CF Standard Name Table v79")

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
            if ncVar.name != 'l1a_data':
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

def getChunking(ncDims):
    global ncFile
    sizes = {
        'bbt': 5,
        'ndatas': 6,
        'bsamp': 5,
        'sline': 7,
        'nsamp': 400,
        'ndatas12': 72,
        'lines': 32,
        'nrec': 1,
        'ntilt': 3,
        'gains': 4,
        'reft': 1,
        'ntbl': 1,
        'pixels': 128,
        'err': 1,
        'vec': 3,
        'instt': 4,
        'instr': 120,
        'points': 1440,
        'sets': 7,
        'mat': 3,
        'detectors': 10,
        'odatas22': 22,
        'rec': 32,
        'tilts': 1,
        'dets': 1,
        'ninf': 1,
        'rads': 1,
        'blines': 200,
        'coefs': 4,
        'pxls': 7,
        'pref': 128,
        'pairs': 2,
        'flag': 24,
        'sdet': 2,
        'ndatas22': 132,
        'bands': 8,
        'ndatas14': 84,
        'bdatas': 2,
        'odatas': 1,
    }
    chunks = []
    for dim in ncDims:
        dimsize = len(ncFile.dimensions[dim])
        if sizes[dim] <= dimsize:
            chunks.append(sizes[dim])
        else:
            chunks.append(1)
    return chunks

# Retrieves the HDF4 dataset and get all the required information to build
# the NetCDF version
def copyDatasets():
    global hdfFile
    global ncFile
    global errorAt
    skipSDS = (
        'ntilts',
        'entry_year',
        'entry_day',
        'ref_year',
        'ref_day',
        'ref_minute'
    )
    try:
        print("Copying datasets/variables...")

        datasetNames = hdfFile.datasets().keys()
        for name in datasetNames:
            if name in skipSDS:
                continue
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
            # get the dimensions to be set as a tuple
            #ncDims = tuple(map(lambda dim: hdfDims))

            # create the variable and then assign the data
            chunks = getChunking(hdfDims)
            newVariable = ncFile.createVariable(convertToNetcdfName(name), ncDatatype, hdfDims, chunksizes=chunks, zlib=True)
            newVariable[:] = data # slice entire variable and assign it 
    
            # netcdf assigning attributes to the variables. NcVar is an object and can be
            # passed by reference, so the function takes care of it
            assignNcVarAttr(newVariable, hdfDatasetAttr, ncDatatype)
            if name == 'msec':
                newVariable.setncattr("units","milliseconds")
                newVariable.setncattr("comment","milliseconds since start of day, day boundary crossing indicated if value is smaller than previous scan")

    except Exception as e:
        print(f"-E- Error copying datasets/variables. Error occurred with HDF4 dataset named <{errorAt}>")
        print(f"Reason: {e}")
        exit()

def addTimeVariable():
    global ncFile

    try:
        print("adding CF-compliant time variable...")

        msecs = ncFile['msec'][:]
        timevar = np.zeros(len(msecs),np.float64)
        time_coverage_start = datetime.fromisoformat(ncFile.time_coverage_start[:23])
        day_start = datetime.strptime(time_coverage_start.strftime('%Y%m%d'),'%Y%m%d')
        for i,msec in enumerate(msecs):
            seconds = msec / 1000.
            scantime = day_start + timedelta(seconds=seconds)
            if scantime < time_coverage_start:
                scantime = scantime+ timedelta(days=1)
            timevar[i] = scantime.strftime('%s.%f')
        chunksize = [32]
        if len(ncFile.dimensions['rec']) < chunksize[0]:
            chunksize[0] = 1
        timeVariable = ncFile.createVariable('time', np.float64, 'rec', chunksizes=chunksize, zlib=True)
        timeVariable[:] = timevar 
        timeVariable.setncattr("long_name","time")
        timeVariable.setncattr("units","seconds since 1970-1-1")

    except Exception as e:
        print(f"-E- Error occurred with adding time dataset...")
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


# Given a list of variables, extract the dimension data from it 
# pyhdf lib does not have a way to get only the dims, you need to 
# go through the variables and get it
def setDimensions():
    global hdfFile
    global ncFile
    global errorAt

    errorAt = "Setting Dimensions"
    # create a set so that when parsing all the variables with their dims,
    # dims do not get repeated in the set
    dimSet = set()
    dimSet.add(("tilts", 1))
    dimSet.add(("bands", 8))
    dimSet.add(("bbt", 5))
    dimSet.add(("ndatas", 6 ))
    # dimSet.add(("bsamp", 5 ))
    # dimSet.add(("sline", 22))
    # dimSet.add(("nsamp", 400))
    dimSet.add(("ndatas12", 72))
	# lines = 2196
    dimSet.add(("nrec", 1))
    dimSet.add(("ntilt", 3))
    dimSet.add(("gains", 4))
    dimSet.add(("reft", 1))
    dimSet.add(("ntbl", 1))
    dimSet.add(("pixels", 2222))
    dimSet.add(("err", 1))
    dimSet.add(("vec", 3))
    dimSet.add(("instt", 4))
    dimSet.add(("instr", 120))
    dimSet.add(("points", 1440))
    # dimSet.add(("sets", 22))
    dimSet.add(("mat", 3))
    dimSet.add(("detectors", 10))
    dimSet.add(("odatas22", 22))
    # dimSet.add(("rec", 1098))
    dimSet.add(("tilts", 1))
    dimSet.add(("dets", 1))
    dimSet.add(("ninf", 1))
    dimSet.add(("rads", 1))
    # dimSet.add(("blines", 10980))
    dimSet.add(("coefs", 4))
    # dimSet.add(("pxls", 21))
    dimSet.add(("pref", 128))
    dimSet.add(("pairs", 2))
    dimSet.add(("flag", 24))
    dimSet.add(("sdet", 2))
    dimSet.add(("ndatas22", 132))
    dimSet.add(("ndatas14", 84))
    dimSet.add(("bdatas", 2))
    dimSet.add(("odatas", 1))
    try:
        errorAt = "extract dimension details from HDF file"
        # get key, value pair of all the HDF
        hdfDatasets = hdfFile.datasets()
        vars = ["l1a_data",
                "lat",
                "blk_data",
                "start_line",
                "samp_table"]
        
        for var in vars:
            errorAt = f"extract {var} from HDF file"
            varInfo = hdfDatasets[var] # returns something like: (('bands',), (8,), 22, 0
            dimNames = varInfo[0]   # tuple of dim names for var
            dimVals = varInfo[1]    # tuple of dim values for the dims

            # go though each dim and save it into the set
            for i in range(len(dimNames)):
                name = dimNames[i]
                val = dimVals[i]
                newDim = (convertToNetcdfName(name), val)
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

def main():
    print(f"l1aconvert_octs {__version__}")

    # get global reference to the file variables
    global ncFile
    global hdfFile

    parser = argparse.ArgumentParser(description='''
                                     Given an HDF4 Scientific Dataset (SD), convert it
                                     into a NetCDF (NC) file format with the same name.
                                     If the optional argument [oFile] name is given, 
                                     output as the oFile name.
                                     ''')

    parser.add_argument("ifile", type=str, help="input HDF4 File")
    parser.add_argument("ofile", nargs='?', type=str, help="Optional output file name; default = input + '.nc'")
    parser.add_argument("--doi","-d",type=str,default=None,help="DOI string for metadata")

    # Parse the command-line arguments
    # if oFile isn't given, then use the filename with .nc extension
    args = parser.parse_args()
    fileName = args.ifile
    oFileName = fileName + ".nc" if args.ofile is None else args.ofile
    oFileName = Path(oFileName)

    print(f"Input File:\t{fileName}")
    print(f"Output File:\t{oFileName}\n")

    # Opening the hdf4 file for reading:
    try:
        hdfFile = HDF(fileName, HDF_READ)
        print(f"Opening file:\t{fileName}")
    except:
        print(f"\n-E- Error opening file named: {fileName}.\n Make sure the filetype is hdf4.\n")
        exit()

    # if opened successfully, create netcdf file to be written into
    ncFile = NC(oFileName, NC_WRITE, NC_FORMAT)
    setDimensions()
    copyGlobalAttributes(fileName, oFileName)
    if args.doi:
        ncFile.setncattr("doi",args.doi)
    copyDatasets()
    addTimeVariable()
    closeFiles()
    print("Conversion of L1A HDF4 to NETCDF successful")

if __name__ == '__main__':
    main()