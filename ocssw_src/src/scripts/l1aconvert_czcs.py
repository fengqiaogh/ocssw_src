#! /usr/bin/env python3

import argparse
import sys
from netCDF4 import Dataset as NC
from pyhdf.SD import SD as HDF
from pyhdf.SD import SDC
import numpy as np
from datetime import datetime, timedelta, timezone
datetime.isoformat

## Convert an L1A CZCS HDF File into a L1A CZCS NetCDF File
## Reads the contents of the input HDF file and outputs it as a NetCDF
## To generate a L1A file, download 

__version__ = "1.0.1 2025-01-16"

# Constants
HDF_READ = SDC.READ
NC_WRITE = "w"
NC_FORMAT = "NETCDF4"

# Do not copy these to NetCDF. 
# year, day and milisecs will be made into isodate
SKIP_ATTRIBUTES = (
    "Product Name",
    "Data Type",
    "Replacement Flag",
    "Data Center",
    "Station Name",
    "Station Latitude",
    "Station Longitude",
    "Mission",
    "Sensor",
    "Software ID",
    "Processing Time",
    "Creation Time",
    "Input Files",
    "Processing Control",
    "Pixels per Scan Line",
    "Number of Pixel Control Points",
    "Number of Scan Control Points",
    "Scene Center Scan Line",
    "Filled Scan Lines",
    "Start Time",
    "End Time",
    "Start Year",
    "Start Day",
    "Start Millisec",
    "End Year",
    "End Day",
    "End Millisec",
    "Latitude Units",
    "Longitude Units",
    "Scene Center Time",
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
    "Start Center Latitude",
    "Start Center Longitude",
    "End Center Latitude",
    "End Center Longitude",
    "Northernmost Latitude",
    "Southernmost Latitude",
    "Westernmost Longitude",
    "Easternmost Longitude"
    "Pixels per Scan Line",
    "Number of Scan Lines"
)


# python has no switch, so this dict is used to get the NetCDF type equivalent based on HDF4 types.
# NetCDF API sets data types for variables based on NumPy's type.
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
errorAt = "" # name of whatever dataset the program was trying to copy/do
# Convert HDF naming to NetCDF naming
# ie: "Data Center" to "data_center"


def convertToNetcdfName(hdfname:str):
    ncname = hdfname.lower().replace(" ", "_")
    if ncname == "number_of_scan_lines": return "scans"
    elif ncname == "pixels_per_scan_line": return "pixels"
    elif ncname == "number_of_pixel_control_points": return "control_points"
    elif ncname == "num_qual": return "quality_flags"
    elif ncname == "number_of_bands": return "bands"
    elif ncname == "start_node": return "startDirection"
    elif ncname == "end_node": return "endDirection"
    elif ncname == "msec": return "scan_time"
    elif ncname.startswith('lac'): return ncname.replace('lac','LAC')
    else: return ncname

# Convert HDF naming to NetCDF naming
# ie: "Data Center" to "data_center"
def convertToNetcdfDimension(dim:str, sds:str):
    # HDF dimensions are not all named, so we need to give them NetCDF names
    if (dim == "3"): return "vector_elements"
    elif (dim == "6"): return "bands"
    else: 
        ncdim = dim.lower().replace(" ", "_")
        if ncdim == "number_of_scan_lines": return "scans"
        elif ncdim == "pixels_per_scan_line": return "pixels"
        elif ncdim == "number_of_pixel_control_points": return "control_points"
        elif ncdim == "num_qual": return "quality_flags"
        elif ncdim == "number_of_bands": return "bands"
        else: return ncdim

# Copy global attributes from the HDF4 file into the new NetCDF file 
def copyGlobalAttributes(iFile:str, oFile:str):
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


            ## ...Add more but these are the only 2 instances I see in SEAWIFS global attributes
            ncFile.setncattr(convertToNetcdfName(name), val)

        # make isodate for start/end time
        errorAt = "time_coverage_start"
        syear = globalAttr.get("Start Year")
        sday = globalAttr.get("Start Day")
        ssec = "{:.3f}".format(float(globalAttr.get("Start Millisec"))/1000.)
        fssec = ssec.split('.')[1]

        start_of_year = datetime(syear, 1, 1,tzinfo=timezone.utc)
        stime = start_of_year + timedelta(days=(sday-1), seconds=float(ssec))
        
        ncFile.setncattr("time_coverage_start", str(stime.strftime('%Y-%m-%dT%H:%M:%S.') + fssec + 'Z'))

        # make isodate for end time
        errorAt = "time_coverage_end"
        eyear = globalAttr.get("End Year")
        eday = globalAttr.get("End Day") - 1
        esec = "{:.3f}".format(float(globalAttr.get("End Millisec"))/1000.)
        fesec = esec.split('.')[1]

        end_of_year = datetime(eyear, 1, 1, tzinfo=timezone.utc)

        etime = end_of_year + timedelta(days=eday, seconds=float(esec))

        ncFile.setncattr("time_coverage_end",str(etime.strftime('%Y-%m-%dT%H:%M:%S.') + fesec + 'Z'))

        # add converter version
        ncFile.setncattr("l1aconvert_czcs_version", __version__)
        
        # netcdf has history global attribute, so add one here:
        # to convert a file, it's static except for the input and output files
        ncFile.setncattr("history", "{}; python3 l1acovert_czcs {} {}".format(globalAttr.get("Processing Control")[:-1], iFile, oFile))

        #add date_created
        errorAt = "date_created"
        create_date = datetime.now(tz=timezone.utc)

        ncFile.setncattr("date_created",str(create_date.strftime('%Y-%m-%dT%H:%M:%SZ')))
        
        #add date_created
        errorAt = "product_name"
        ncFile.setncattr("product_name","{}".format(oFile))

        errorAt = "static metadata"
        # some static metadata
        ncFile.setncattr("cdm_data_type","swath")
        ncFile.setncattr("Conventions","CF-1.8, ACDD-1.3")
        ncFile.setncattr("creator_email","data@oceancolor.gsfc.nasa.gov")
        ncFile.setncattr("creator_name","NASA/GSFC/OBPG")
        ncFile.setncattr("creator_url","https://oceancolor.gsfc.nasa.gov")
        ncFile.setncattr("institution","NASA Goddard Space Flight Center, Ocean Biology Processing Group")
        ncFile.setncattr("keywords_vocabulary","NASA Global Change Master Directory (GCMD) Science Keywords")
        ncFile.setncattr("license","https://www.earthdata.nasa.gov/engage/open-data-services-and-software/data-and-information-policy")
        ncFile.setncattr("naming_authority","gov.nasa.gsfc.oceancolor")
        ncFile.setncattr("platform","Nimbus-7")
        ncFile.setncattr("instrument","CZCS")
        ncFile.setncattr("processing_level","L1A")
        ncFile.setncattr("processing_version","3")
        ncFile.setncattr("project","Ocean Biology Processing Group")
        ncFile.setncattr("publisher_email","data@oceancolor.gsfc.nasa.gov")
        ncFile.setncattr("publisher_name","NASA/GSFC/OB.DAAC")
        ncFile.setncattr("publisher_url","https://oceancolor.gsfc.nasa.gov")
        ncFile.setncattr("standard_name_vocabulary","CF Standard Name Table v79")

    except:
        print(f"-E- Error copying global attributes. Was processing <{errorAt}> from HDF4 when error was caught.")
        exit(1)
    errorAt = "" # reset 

# Given a NetCDF variable, assign the attributes in attrList.
# attrList is the attribute list from the HDF4 dataset and it is copying
# that over to the NetCDF version
#
# When assigning valid range, it is in string format. The dataType parameter
# helps with the slicing and converting it into the right range
def assignNcVarAttr(ncVar, attrDict:dict, dataType):

    for attr in attrDict.keys():
        if (attr == "long_name" or attr == "long name"):
            ncVar.long_name = attrDict.get(attr)
        if (attr == "units"):
            if ncVar.name[:4] != 'band':
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
        'scans': 32,
        'quality_flags': 1,
        'bands': 6,
        'pixels': 128,
        'control_points': 25,
        'vector_elements': 3,
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
        'cntl_pt_rows',
        'slat',
        'slon',
        'clat',
        'clon',
        'elat',
        'elon',
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
            ncDims = tuple(map(lambda dim: convertToNetcdfDimension(dim, name), hdfDims))

            # create the variable and then assign the data
            chunks = getChunking(ncDims)
            newVariable = ncFile.createVariable(convertToNetcdfName(name), ncDatatype, ncDims, chunksizes=chunks, zlib=True)
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
        exit(1)

def addTimeVariable():
    global ncFile

    try:
        print("adding CF-compliant time variable...")

        msecs = ncFile['scan_time'][:]
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
        if len(ncFile.dimensions['scans']) < chunksize[0]:
            chunksize[0] = 1
        timeVariable = ncFile.createVariable('time', np.float64, 'scans', chunksizes=chunksize, zlib=True)
        timeVariable[:] = timevar 
        timeVariable.setncattr("long_name","time")
        timeVariable.setncattr("units","seconds since 1970-1-1")

    except Exception as e:
        print(f"-E- Error copying datasets/variables. Error occurred with adding time dataset...")
        print(f"Reason: {e}")
        exit(1)

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
    dimSet.add(("bands", 6))
    dimSet.add(("vector_elements", 3))
    dimSet.add(("quality_flags",5))

    try:
        errorAt = "extract dimension details from HDF file"
        # get key, value pair of all the HDF
        hdfDatasets = hdfFile.datasets()
        var = "band1"
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
        
        # copy over number of pixel control points
        errorAt = "Number of Pixel Control Points"
        var = "cntl_pt_cols"
        varInfo = hdfDatasets[var]
        dimNames = varInfo[0]   # tuple of dim names for var
        dimVals = varInfo[1]    # tuple of dim values for the dims
        newDim = (convertToNetcdfName(dimNames[0]), dimVals[0])
        dimSet.add(newDim)  

        # after getting all the dims, go through and set them in the NetCDF file
        
        for dim in dimSet:
            errorAt = f"set {dim} in the NetCDF file"
            name = dim[0]
            val = dim[1]
            ncFile.createDimension(name, val)

        

    except:
        print(f"-E- Error copying dimensions. Was trying to {errorAt} when error was caught.")
        exit(1)

def main():
    print(f"l1aconvert_czcs {__version__}")

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


    print(f"Input File:\t{fileName}")
    print(f"Output File:\t{oFileName}\n")

    # Opening the hdf4 file for reading:
    try:
        hdfFile = HDF(fileName, HDF_READ)
        print(f"Opening file:\t{fileName}")
    except:
        print(f"\n-E- Error opening file named: {fileName}.\n Make sure the filetype is hdf4.\n")
        exit(1)

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