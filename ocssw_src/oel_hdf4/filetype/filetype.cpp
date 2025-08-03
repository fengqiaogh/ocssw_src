#include <iostream>
#include <netcdf.h>
#include <hdf5.h>
#include <libgen.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <boost/algorithm/string.hpp>

#include <netcdf>
#include <genutils.h>
#include <sensorDefs.h>
#include <sensorInfo.h>

#include "aviris.h"
#include "epr_api.h"

#include "filetype.h"

#define EOSMETALEN 32768

using namespace netCDF;
using namespace netCDF::exceptions;
using namespace std; 

static file_format chk_oli(char *filename);
static file_format chk_l5tm(char *filename);
static file_format chk_l7etm(char *filename);
static file_format chk_aviris(char *filename, char *hdrfile, char *imgfile, char *navfile, char *gainfile);
static file_format chk_prism(char *filename);
static file_format chk_seabass(char *filename);

// defined in filetypeXml.cpp
extern "C" file_format chk_safe_xml(char *filename);

static int hdf5AttributeReadString(hid_t group_id, const char* name, char* val) {
    int result = 0;
    herr_t status;
    hid_t attribute_id = H5Aopen(group_id, name, H5P_DEFAULT);
    if(attribute_id >= 0) {
        hid_t atype = H5Aget_type(attribute_id);
        H5T_class_t type_class = H5Tget_class(atype);
        if (type_class == H5T_STRING) {
            hid_t atype_mem = H5Tget_native_type(atype, H5T_DIR_ASCEND);
            status = H5Aread(attribute_id, atype_mem, val);
            if(status >= 0) {
                result = 1;
            }
            H5Tclose(atype_mem);
        }
        H5Tclose(atype);
    }
    H5Aclose(attribute_id);
    return result;
}

static int hdf5AttributeStartsWith(hid_t group_id, const char* name, char* val) {
    char buf[1024];
    if(hdf5AttributeReadString(group_id, name, buf)) {
        if(strncmp(buf, val, strlen(val)) == 0) 
            return 1;
    }
    return 0;
}

static file_format chk_hdf5(char *filename) {
    file_format ret = {FT_INVALID, -1, -1};
    hid_t file_id, group_id, group2_id;
    char buf[1024];

    /* Save old error handler */
    H5E_auto_t old_func;
    void *old_client_data;
    H5Eget_auto(H5E_DEFAULT, &old_func, &old_client_data);

    /* Turn off error handling */
    H5Eset_auto(H5E_DEFAULT, NULL, NULL);

    file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    if(file_id >= 0) {

        // check for VIIRS SDR
        if(hdf5AttributeReadString(file_id, "Platform_Short_Name", buf)) {
            group_id = H5Gopen2(file_id, "Data_Products/VIIRS-M1-SDR", H5P_DEFAULT);
            
            if(group_id >= 0) {
                if(hdf5AttributeStartsWith(group_id, "Instrument_Short_Name", "VIIRS")) {

                    if(strstr(buf, "NPP") || strstr(buf, "NPOESS")) {
                        ret.type = FT_VIIRSL1B;
                        ret.sensor_id = VIIRSN;
                        ret.subsensor_id = VIIRS_NPP;
                        if (want_verbose) {
                            printf("Input file %s is a VIIRS NPP SDR L1B HDF 5 file.\n", filename);
                        }

                    } else if(strstr(buf, "JPSS-1") || strstr(buf, "J01")) {
                        ret.type = FT_VIIRSL1B;
                        ret.sensor_id = VIIRSJ1;
                        ret.subsensor_id = VIIRS_J1;
                        if (want_verbose) {
                            printf("Input file %s is a VIIRS JPSS-1 SDR L1B HDF 5 file.\n", filename);
                        }

                    } else if(strstr(buf, "JPSS-2") || strstr(buf, "J02")) {
                        ret.type = FT_VIIRSL1B;
                        ret.sensor_id = VIIRSJ2;
                        ret.subsensor_id = VIIRS_J2;
                        if (want_verbose) {
                            printf("Input file %s is a VIIRS JPSS-2 SDR L1B HDF 5 file.\n", filename);
                        }
                    }

                }
                H5Gclose(group_id);
            }
        }

        // check for HICO
        if(ret.type == FT_INVALID) {
            group_id = H5Gopen2(file_id, "metadata/FGDC/Identification_Information/Platform_and_Instrument_Identification", H5P_DEFAULT);
            
            if(group_id >= 0) {
                if(hdf5AttributeStartsWith(group_id, "Instrument_Short_Name", "hico")) {
                    group2_id = H5Gopen2(file_id, "metadata/FGDC/Identification_Information/Processing_Level", H5P_DEFAULT);
                    
                    if(group2_id >= 0) {
                        if(hdf5AttributeStartsWith(group2_id, "Processing_Level_Identifier", "Level-1B")) {
                            ret.type = FT_HICOL1B;
                            ret.sensor_id = HICO;
                            if (want_verbose) 
                                printf("Input file %s is a HICO L1B HDF 5 file.\n", filename);
                        }
                        H5Gclose(group2_id);
                    }
                }
                H5Gclose(group_id);
            }
        }


        // check for GOCI
        if(ret.type == FT_INVALID) {
            group_id = H5Gopen2(file_id, "HDFEOS/POINTS/Scene Header", H5P_DEFAULT);
            
            if(group_id >= 0) {
                if(hdf5AttributeStartsWith(group_id, "Scene Title", "GOCI Level-1B Data")) {
                    ret.type = FT_GOCIL1B;
                    ret.sensor_id = GOCI;
                    if (want_verbose) 
                        printf("Input file %s is a GOCI L1B HDF 5 file.\n", filename);
                }
                H5Gclose(group_id);
            }
        }

        // check for SGLI
        if(ret.type == FT_INVALID) {
            group_id = H5Gopen2(file_id, "Global_attributes", H5P_DEFAULT);
            
            if(group_id >= 0) {
                if(hdf5AttributeStartsWith(group_id, "Satellite", "Global Change Observation Mission - Climate (GCOM-C)")) {
                    
                    if(hdf5AttributeStartsWith(group_id, "Sensor", "Second-generation Global Imager (SGLI)")) {
                        if(hdf5AttributeStartsWith(group_id, "Product_level", "Level-1B")) {
        
                            if(hdf5AttributeStartsWith(group_id, "Product_name", "Top of atmosphere radiance (reflectance)")) {
                                ret.type = FT_SGLI;
                                ret.sensor_id = SGLI;
                                if (want_verbose) 
                                    printf("Input file %s is a SGLI L1B HDF 5 file.\n", filename);
                            }
                        }
                    }
                }
                H5Gclose(group_id);
           }
        }

        H5Fclose(file_id);
    }

    /* Restore previous error handler */
    H5Eset_auto(H5E_DEFAULT, old_func, old_client_data);

    return ret;
}

file_type getFormatType(char *filename) {
    file_format ret = getFormat(filename);
    return ret.type;
}

// remove the hidden characters from a string read from a NetCDF or HDF4 file
static void sanitizeString(string &str) {
    // get rid of the trailing 0
    str = str.c_str();

    // now delete begin/end whitespace
    boost::trim(str);
}

/**
 * @brief Check input file format 
 * - HDF5
 * - Netcdf
 * - HDF4 
 * - HDF-EOS L1B
 * - MISR
 * - SeaBASS
 * - MERIS
 * - OLI
 * - AVERIS
 * - PRISM
 * - Landsat 5
 * - Landsat 7     
 * 
 * @param filename 
 * @return FT_INVALID or the file format code 
*/
file_format getFormat(char *filename) {
    file_format ret = { FT_INVALID, -1, -1 };
    int result = check_url(filename);
    // Check if file exists
    if (result == 0) {
        if (access(filename, F_OK) || access(filename, R_OK)) {
            printf("-E- %s: Input file '%s' does not exist or cannot open.\n", __FILE__, filename);
            return ret;
        }

        // Check if HDF5
        if ((ret = chk_hdf5(filename)).type != FT_INVALID) {
            return ret;
        }
    }

    // see if it is a netCDF file
    try {
        NcFile ncFile(filename, NcFile::read);

        // catch netCDF errors
        try {

            string platformStr;
            string instrumentStr;
            string titleStr;
            string processingLevelStr;
            NcGroupAtt title; 

            if (!ncFile.getAtt("platform").isNull()) {
                NcGroupAtt platform = ncFile.getAtt("platform");
                platform.getValues(platformStr);
                sanitizeString(platformStr);
                transform(platformStr.begin(), platformStr.end(), platformStr.begin(), [](unsigned char c){ return tolower(c); });
            }

            if (!ncFile.getAtt("instrument").isNull()) {
                NcGroupAtt instrument = ncFile.getAtt("instrument");
                instrument.getValues(instrumentStr);
                sanitizeString(instrumentStr);
            }

            // Check if Netcdf
            // our netCDF files all have "title" attribute
            if (!ncFile.getAtt("title").isNull()) {
                // File is netcdf 
                title = ncFile.getAtt("title");
                title.getValues(titleStr);
                sanitizeString(titleStr);
                transform(titleStr.begin(), titleStr.end(), titleStr.begin(), [](unsigned char c){ return tolower(c); });

                if (titleStr.find("viirs level-1a") != string::npos) {
                    if (platformStr.find("suomi-npp") != string::npos) {
                        ret.type = FT_VIIRSL1A;
                        ret.sensor_id = VIIRSN;
                        ret.subsensor_id = VIIRS_NPP;
                        if (want_verbose) {
                            printf("Input file %s is VIIRS NPP L1A NetCDF4.\n", filename);
                        }
                    }
                    else if (platformStr.find("jpss-1") != string::npos) {
                        ret.type = FT_VIIRSL1A;
                        ret.sensor_id = VIIRSJ1;
                        ret.subsensor_id = VIIRS_J1;
                        if (want_verbose) {
                            printf("Input file %s is VIIRS JPSS-2 L1A NetCDF4.\n", filename);                    
                        }
                    }
                    else if (platformStr.find("jpss-2") != string::npos) {
                        ret.type = FT_VIIRSL1A;
                        ret.sensor_id = VIIRSJ2;
                        ret.subsensor_id = VIIRS_J2;
                        if (want_verbose) {
                            printf("Input file %s is VIIRS JPSS-2 L1A NetCDF4.\n", filename);
                        }
                    }

                    return ret;
                }
                
                if (titleStr.find("viirs m-band") != string::npos) {
                    if (titleStr.find("viirs m-band reflected solar band") != string::npos) {
                        ret.type = FT_VIIRSL1BNC;
                        if (want_verbose) 
                            printf("Input file %s is VIIRS NPP L1B NetCDF4.\n", filename);
                    }
                    else if (titleStr.find("viirs m-band geolocation data") != string::npos) {
                        ret.type = FT_VIIRSGEONC;
                        if (want_verbose) 
                            printf("Input file %s is VIIRS NPP GEO NetCDF4.\n", filename);
                    }
                    if (platformStr == "suomi-npp") {
                        ret.sensor_id = VIIRSN;
                        ret.subsensor_id = VIIRS_NPP;
                    } else if (platformStr == "jpss-1") {
                        ret.sensor_id = VIIRSJ1;
                        ret.subsensor_id = VIIRS_J1;
                    } else if (platformStr == "jpss-2") {
                        ret.sensor_id = VIIRSJ2;
                        ret.subsensor_id = VIIRS_J2;
                    }
                    return ret;
                }

                if (titleStr.find("meris l1b") != string::npos) {
                    ret.type = FT_MERISCC;
                    ret.sensor_id = MERIS;
                    ret.subsensor_id = -1;
                    if (want_verbose) 
                        printf("Input file %s is MERISS CC file.\n", filename);
                    return ret;
                }

                if (titleStr.find("aviris level-1b") != string::npos) {
                    ret.type = FT_AVIRIS;
                    ret.sensor_id = AVIRIS;
                    ret.subsensor_id = -1;
                    if (want_verbose) {
                        printf("Input file %s is AVIRIS L1B NetCDF4 file.\n", filename);
                    }
                    return ret;
                } 

                if (titleStr.find("olci level 1b product") != string::npos) {
                    string subtype;

                    if (titleStr.find("geo coordinates") != string::npos) {
                        ret.type = FT_OLCIGEO;
                        subtype = "GEO";
                    }
                    else {
                        ret.type = FT_OLCI;
                        subtype = "L1B";
                    }
                
                    NcGroupAtt product = ncFile.getAtt("product_name");
                    string product_name;
                    product.getValues(product_name);
                    sanitizeString(product_name);

                    if (product_name.find("S3A") != string::npos) {
                        ret.sensor_id = OLCIS3A;
                        ret.subsensor_id = OLCI_S3A;
                        if (want_verbose) {
                            printf("Input file %s is OLCI S3A %s file.\n", filename, subtype.c_str());
                        }
                    } else if (product_name.find("S3B") != string::npos) {
                        ret.sensor_id = OLCIS3B;
                        ret.subsensor_id = OLCI_S3B;
                        if (want_verbose) {
                            printf("Input file %s is OLCI S3B %s file.\n", filename, subtype.c_str());
                        }
                    } else {
                        if (want_verbose) {
                            printf("Input file %s is OLCI %s file. Unknown platform.\n", filename, subtype.c_str());
                            ret.type = FT_INVALID;
                        }
                    }
                    return ret;
                }

                if (titleStr.find("hawkeye level-1a") != string::npos) {
                    ret.type = FT_HAWKEYEL1A;
                    ret.sensor_id = HAWKEYE;
                    ret.subsensor_id = -1;
                    if (want_verbose) 
                        printf("Input file %s is HAWKEYE L1A file.\n", filename);
                    return ret;
                }

                if (titleStr.find("oci level-1b") != string::npos) {
                    ret.type = FT_OCIL1B;
                    ret.sensor_id = OCI;
                    ret.subsensor_id = -1;
                    if (want_verbose) 
                        printf("Input file %s is PACE L1B file.\n", filename);
                    return ret;
                }

                if ((titleStr.find("hkt level-1a") != string::npos) || (titleStr.find("pace hkt data") != string::npos)) {
                    ret.type = FT_HKT;
                    ret.sensor_id = OCI;
                    ret.subsensor_id = -1;
                    if (want_verbose) 
                        printf("Input file %s is PACE HKT file.\n", filename);
                    return ret;                
                }
            
                if (titleStr.find("pace oci level-1c") != string::npos) {
                    ret.type = FT_L1C;
                    ret.sensor_id = OCI;
                    ret.subsensor_id = -1;
                    if (want_verbose)                      
                        printf("Input file %s is PACE OCI L1C file.\n", filename);
                    return ret;
                }

                if (titleStr.find("pace harp2 level-1c") != string::npos) {
                    ret.type = FT_L1C;
                    ret.sensor_id = HARP2;
                    ret.subsensor_id = -1;
                    if (want_verbose)                   
                        printf("Input file %s is PACE HARP2 L1C file.\n", filename);              
                    return ret;
                }

                if (titleStr.find("pace spexone level-1c") != string::npos) {
                    ret.type = FT_L1C;
                    ret.sensor_id = SPEXONE;
                    ret.subsensor_id = -1;
                    if (want_verbose)                       
                        printf("Input file %s is PACE SPEXone L1C file.\n", filename);
                    return ret;
                }

                if (titleStr.find("l1c ancillary file") != string::npos) {
                    ret.type = FT_L1CANC;
                    ret.sensor_id = OCI;
                    ret.subsensor_id = -1;
                    if (want_verbose) 
                        printf("Input file %s is L1C ancillary file.\n", filename);
                    return ret;
                }

                if (titleStr.find("spexone level-1b") != string::npos) {
                    ret.type = FT_SPEXONE;
                    ret.sensor_id = SPEXONE;
                    ret.subsensor_id = -1;
                    if (want_verbose) 
                        printf("Input file %s is PACE SPEXone file.\n", filename);
                    return ret;
                }

                if (titleStr.find("harp2 level-1b") != string::npos) {
                    ret.type = FT_HARP2;
                    ret.sensor_id = HARP2;
                    ret.subsensor_id = -1;
                    if (want_verbose) 
                        printf("Input file %s is PACE HARP2 file.\n", filename);
                    return ret;
                }

                // NetCDF for Seawifs L1A
                if (titleStr.find("seawifs level-1a data") != string::npos) {
                    ret.type = FT_SEAWIFSL1ANC;
                    ret.sensor_id = SEAWIFS;

                    if (!ncFile.getAtt("date_type").isNull()) {
                        NcGroupAtt data = ncFile.getAtt("data_type");
                        string dataTypeStr;
                        data.getValues(dataTypeStr);
                        sanitizeString(dataTypeStr);

                        if (dataTypeStr == "GAC") {
                            ret.subsensor_id = SEAWIFS_GAC;
                            if (want_verbose) 
                                printf("Input file %s is SeaWiFS Level-1A GAC.\n", filename);
                        } else {
                            ret.subsensor_id = SEAWIFS_LAC;
                            if (want_verbose) 
                                printf("Input file %s is SeaWiFS Level-1A LAC.\n", filename);
                        }
                    } else {
                        fprintf(stderr, "-E- %s Line %d: Unable to obtain attribute 'data_type' %s\n", __FILE__, __LINE__, filename);
                    } 
                    return ret;
                }

                if (titleStr.find("octs level-1a gac data") != string::npos) {
                    ret.type = FT_OCTSL1ANC;
                    ret.sensor_id = OCTS;
                    ret.subsensor_id = -1;
                    if (want_verbose) 
                        printf("Input file %s is OCTS Level-1A GAC netCDF.\n", filename);
                    return ret;
                }

                if (titleStr.find("czcs level-1a data") != string::npos) {
                    ret.type = FT_CZCSL1ANC;
                    ret.sensor_id = CZCS;
                    ret.subsensor_id = -1;
                    if (want_verbose) 
                        printf("Input file %s is CZCS Level-1A netCDF.\n", filename);
                    return ret;
                }
                
                if (!ncFile.getAtt("processing_level").isNull()) {
                    NcGroupAtt processingLevel = ncFile.getAtt("processing_level");
                    processingLevel.getValues(processingLevelStr);
                    sanitizeString(processingLevelStr);
                } else {
                    fprintf(stderr, "-E- %s Line %d: Unable to obtain attribute 'processing_level' %s\n", __FILE__, __LINE__, filename);
                    return ret;
                }
                    
                if (processingLevelStr == "L1B") {
                    ret.type = FT_L1BNCDF;
                } else if (processingLevelStr == "L2") {
                    ret.type = FT_L2NCDF;
                } else if (processingLevelStr == "L3 Binned") {
                    ret.type = FT_L3BIN;
                } else if (processingLevelStr == "L3 Mapped") {
                    ret.type = FT_L3MAP;
                }
                    
                if (ret.type != FT_INVALID) {
                    ret.sensor_id = instrumentPlatform2SensorId(instrumentStr.c_str(), platformStr.c_str());
                    ret.subsensor_id = sensorId2SubsensorId(ret.sensor_id);
                    if (want_verbose) 
                        printf("Input file %s is a NetCDF4 %s %s file.\n", filename, instrumentStr.c_str(), processingLevelStr.c_str());
                    return ret;
                }

                return ret;
            }

            // Check if file type is HDF4
            // Our HDF4 files have "Title" attribute
            else if (!ncFile.getAtt("Title").isNull()) {
                title = ncFile.getAtt("Title");
                title.getValues(titleStr);
                sanitizeString(titleStr);

                if (titleStr.find("Level-3 Binned Data") != string::npos) {
                    
                    NcGroupAtt sensor = ncFile.getAtt("Sensor Name");
                    if (!sensor.isNull()) {
                        string sensor_name;
                        sensor.getValues(sensor_name);
                        sanitizeString(sensor_name);
                        
                        // kludge for VIIRS EDR L3
                        if (sensor_name == "VIIRS") {
                            sensor_name = "VIIRSN";
                        }

                        if ((ret.sensor_id = sensorName2SensorId(sensor_name.c_str())) != FT_INVALID) {
                            ret.subsensor_id = sensorId2SubsensorId(ret.sensor_id);
                            ret.type = FT_L3BIN;
                            if (want_verbose) {
                                printf("Input file %s is %s.\n", filename, titleStr.c_str());
                            }
                        } else {
                            fprintf(stderr, "-E- %s Line %d: Unknown sensor name in Level-3 file %s\n", __FILE__, __LINE__, filename);
                            return ret;
                        }

                    } else {
                        fprintf(stderr, "-E- %s Line %d: No sensor name attribute in Level-3 file %s\n", __FILE__, __LINE__, filename);
                        return ret;
                    }
                    return ret;

                } else if (titleStr.find("Level-2 Data") != string::npos) {

                    NcGroupAtt sensor = ncFile.getAtt("Sensor Name");
                    if (!sensor.isNull()) {

                        string sensor_name;
                        sensor.getValues(sensor_name);
                        sanitizeString(sensor_name);
                        
                        ret.sensor_id = sensorName2SensorId(sensor_name.c_str());
                        if (ret.sensor_id != -1) {
                            ret.subsensor_id = sensorId2SubsensorId(ret.sensor_id);
                            ret.type = FT_L2HDF;
                            if (want_verbose) 
                                printf("Input file %s is %s.\n", filename, titleStr.c_str());
                        } else {
                            fprintf(stderr, "-E- %s Line %d: Unknown sensor name in Level-2 file %s\n", __FILE__, __LINE__, filename);
                            return ret;
                        }
                    } else {
                        fprintf(stderr, "-E- %s Line %d: No sensor name attribute in Level-2 file %s\n", __FILE__, __LINE__, filename);
                        return ret;
                    }
                    return ret;

                } else if (titleStr == "SeaWiFS Level-1A Data") {

                    ret.type = FT_SEAWIFSL1A;
                    ret.sensor_id = SEAWIFS;

                    NcGroupAtt data = ncFile.getAtt("Data Type");
                    if (!data.isNull()) {

                        string data_type; 
                        data.getValues(data_type);
                        sanitizeString(data_type);

                        if (data_type == "GAC") {
                            ret.subsensor_id = SEAWIFS_GAC;
                            if (want_verbose) 
                                printf("Input file %s is SeaWiFS Level-1A GAC.\n", filename);
                        } else if (data_type == "LAC") {
                            ret.subsensor_id = SEAWIFS_LAC;
                            if (want_verbose) 
                                printf("Input file %s is SeaWiFS Level-1A LAC.\n", filename);
                        } else {
                            ret.subsensor_id = SEAWIFS_LAC;
                            if (want_verbose) 
                                printf("Input file %s is assumed to be SeaWiFS Level-1A LAC.\n", filename);
                        }
                    } else {
                        ret.subsensor_id = SEAWIFS_LAC;
                        if (want_verbose) {
                            printf("Input file %s is assumed to be SeaWiFS Level-1A LAC.\n", filename);
                        }
                    }

                } else if (titleStr == "OCTS Level-1A GAC Data") {
                    ret.type = FT_OCTSL1A;
                    ret.sensor_id = OCTS;
                    ret.subsensor_id = -1;
                    if (want_verbose) 
                        printf("Input file %s is %s.\n", filename, titleStr.c_str());

                } else if (titleStr == "CZCS Level-1A Data") {
                    ret.type = FT_CZCSL1A;
                    ret.sensor_id = CZCS;
                    ret.subsensor_id = -1;
                    if (want_verbose) 
                        printf("Input file %s is %s.\n", filename, titleStr.c_str());

                } else if (titleStr == "OCM1 Level-1B (OBPG)") {
                    ret.type = FT_OCML1B;
                    ret.sensor_id = OCM1;
                    ret.subsensor_id = -1;
                    if (want_verbose) 
                        printf("Input file %s is %s.\n", filename, titleStr.c_str());

                } else if (titleStr.substr(0, 12) == "OCM Level-1B") {
                    ret.type = FT_OCML1BDB;
                    ret.sensor_id = OCM1;
                    ret.subsensor_id = -1;
                    if (want_verbose) 
                        printf("Input file %s is %s.\n", filename, titleStr.c_str());
                    
                } else if (titleStr == "Oceansat OCM2 Level-1B Data") {
                    ret.type = FT_OCM2L1B;
                    ret.sensor_id = OCM2;
                    ret.subsensor_id = -1;
                    if (want_verbose) 
                        printf("Input file %s is %s.\n", filename, titleStr.c_str());

                    /* generic L1B format support */
                } else if (titleStr == "SeaWiFS Level-1B") {
                    ret.type = FT_L1HDF;
                    ret.sensor_id = SEAWIFS;
                    ret.subsensor_id = -1;
                    if (want_verbose) 
                        printf("Input file %s is %s.\n", filename, titleStr.c_str());

                } else if (titleStr == "MERIS Level-1B") {
                    ret.type = FT_L1HDF;
                    ret.sensor_id = MERIS;
                    ret.subsensor_id = -1;
                    if (want_verbose) {
                        printf("Input file %s is %s.\n", filename, titleStr.c_str());
                    }
                } else if (titleStr == "VIIRS Level-1B" || titleStr == "VIIRSN Level-1B") {
                    ret.type = FT_L1HDF;
                    ret.sensor_id = VIIRSN;
                    ret.subsensor_id = VIIRS_NPP;
                    if (want_verbose) 
                        printf("Input file %s is %s.\n", filename, titleStr.c_str());
                    
                } else if (titleStr == "VIIRSJ1 Level-1B") {
                    ret.type = FT_L1HDF;
                    ret.sensor_id = VIIRSJ1;
                    ret.subsensor_id = VIIRS_J1;
                    if (want_verbose) 
                        printf("Input file %s is %s.\n", filename, titleStr.c_str());

                } else if (titleStr == "VIIRSJ2 Level-1B") {
                    ret.type = FT_L1HDF;
                    ret.sensor_id = VIIRSJ2;
                    ret.subsensor_id = VIIRS_J2;
                    if (want_verbose) {
                        printf("Input file %s is %s.\n", filename, titleStr.c_str());
                    }
                } else if (titleStr == "OCM2 Level-1B") {
                    ret.type = FT_L1HDF;
                    ret.sensor_id = OCM2;
                    ret.subsensor_id = -1;
                    if (want_verbose) 
                        printf("Input file %s is %s.\n", filename, titleStr.c_str());
                    
                } else if (titleStr == "OCTS Level-1B") {
                    ret.type = FT_L1HDF;
                    ret.sensor_id = OCTS;
                    ret.subsensor_id = -1;
                    if (want_verbose) 
                        printf("Input file %s is %s.\n", filename, titleStr.c_str());
                    
                } else if (titleStr == "OCTS Level-1B LAC Data") {
                    ret.type = FT_OCTSL1B;
                    ret.sensor_id = OCTS;
                    ret.subsensor_id = -1;
                    if (want_verbose) 
                        printf("Input file %s is %s.\n", filename, titleStr.c_str());
                    
                } else if (titleStr == "HMODIST Level-1B" || titleStr == "MODIST Level-1B") {
                    ret.type = FT_L1HDF;
                    ret.sensor_id = MODIST;
                    ret.subsensor_id = MODIS_TERRA;
                    if (want_verbose) {
                        printf("Input file %s is %s.\n", filename, titleStr.c_str());
                        printf("\n\n *** WARNING: polarization can not be computed from generic format ***\n\n");
                    }
                } else if (titleStr == "HMODISA Level-1B" || titleStr == "MODISA Level-1B") {
                    ret.type = FT_L1HDF;
                    ret.sensor_id = MODISA;
                    ret.subsensor_id = MODIS_AQUA;
                    if (want_verbose) {
                        printf("Input file %s is %s.\n", filename, titleStr.c_str());
                        printf("\n\n *** WARNING: polarization can not be computed from generic format ***\n\n");
                    }

                } else if (titleStr == "CZCS Level-1B") {
                    ret.type = FT_L1HDF;
                    ret.sensor_id = CZCS;
                    ret.subsensor_id = -1;
                    if (want_verbose) 
                        printf("Input file %s is %s.\n", filename, titleStr.c_str());

                } else if (titleStr == "Level-1 cross-calibration pixels") {
                    if (ncFile.getAtt("sensorID").isNull()) {
                        fprintf(stderr, "-E- %s Line %d: Unrecognized sensor name, title %s in input HDF file %s\n", __FILE__, __LINE__, titleStr.c_str(), filename);
                        return ret; 
                    }

                    ret.type = FT_L1XCAL;
                    if (want_verbose) 
                        printf("Input file %s is %s.\n", filename, titleStr.c_str());

                } else if (titleStr == "AVHRR") {
                    ret.type = FT_CLASSAVHRR;
                    ret.sensor_id = AVHRR;
                    ret.subsensor_id = -1;
                    if (want_verbose) 
                        printf("Input file %s is %s.\n", filename, titleStr.c_str());

                } else if (titleStr.find("oceansat-1") != string::npos) {
                    if (!ncFile.getAtt("satellite").isNull()) {
                        ret.type = FT_OCML1BDB;
                        ret.sensor_id = OCM1;
                        ret.subsensor_id = -1;
                        if (want_verbose) 
                            printf("Input file %s is OCM1 DB file\n", filename); 
                    }

                } else {
                    fprintf(stderr, "-E- %s Line %d: Unrecognized title %s in input HDF file %s\n", __FILE__, __LINE__, titleStr.c_str(), filename);
                    return ret;
                }

                return ret; 

            // check for HDF-EOS L1B format
            } else if (!ncFile.getAtt("ArchiveMetadata.0").isNull()) {

                NcGroupAtt archive_metadata = ncFile.getAtt("ArchiveMetadata.0");
                string archive_metadata_str;
                archive_metadata.getValues(archive_metadata_str);
                sanitizeString(archive_metadata_str);

                if (archive_metadata_str.find("MODIS/Terra Calibrated Radiances 5-Min L1B Swath 1km") != string::npos) {
                    ret.type = FT_MODISL1B;
                    ret.sensor_id = MODIST;
                    ret.subsensor_id = MODIS_TERRA;

                    if (want_verbose) 
                        printf("Input file %s is MODIS Terra Level-1B HDF-EOS product.\n", filename); 
                    return ret;

                } else if (archive_metadata_str.find("MODIS/Aqua Calibrated Radiances 5-Min L1B Swath 1km") != string::npos) {
                    ret.type = FT_MODISL1B;
                    ret.sensor_id = MODISA;
                    ret.subsensor_id = MODIS_AQUA;

                    if (want_verbose) 
                        printf("Input file %s is MODIS Aqua Level-1B HDF-EOS product.\n", filename);
                    return ret;

                } else if (archive_metadata_str.find("MODIS/Aqua Calibrated Radiances 5-Min L1B Swath 250m") != string::npos) {
                    printf("Input file %s is MODIS Aqua Level-1B HDF-EOS product at 250m Resolution.\nThis file is no longer accepted as input.\nPlease use the 1KM (LAC) file and set resolution=250\n", filename);
                    return ret;
                } else if (archive_metadata_str.find("MODIS/Aqua Calibrated Radiances 5-Min L1B Swath 500m") != string::npos) {
                    printf("Input file %s is MODIS Aqua Level-1B HDF-EOS product at 500m Resolution.\nThis file is no longer accepted as input.\nPlease use the 1KM (LAC) file and set resolution=500\n", filename);
                    return ret;
                } else if (archive_metadata_str.find("MODIS/Terra Calibrated Radiances 5-Min L1B Swath 250m") != string::npos) {
                    printf("Input file %s is MODIS Terra Level-1B HDF-EOS product at 250m Resolution.\nThis file is no longer accepted as input.\nPlease use the 1KM (LAC) file and set resolution=250\n", filename);
                    return ret;
                } else if (archive_metadata_str.find("MODIS/Terra Calibrated Radiances 5-Min L1B Swath 500m") != string::npos) {
                    printf("Input file %s is MODIS Terra Level-1B HDF-EOS product at 500m Resolution.\nThis file is no longer accepted as input.\nPlease use the 1KM (LAC) file and set resolution=500\n", filename);
                    return ret;
                } else if (archive_metadata_str.find("MODIS/Aqua Geolocation Fields") != string::npos) {
                    ret.type = FT_MODISGEO;
                    ret.sensor_id = MODISA;
                    ret.subsensor_id = MODIS_AQUA;
                    if (want_verbose) 
                        cout << "Input file " << filename << "is MODIS Aqua Geolocation Fields.\n" << endl;
                    
                } else if (archive_metadata_str.find("MODIS/Terra Geolocation Fields") != string::npos) {
                    ret.type = FT_MODISGEO;
                    ret.sensor_id = MODIST;
                    ret.subsensor_id = MODIS_TERRA;
                    if (want_verbose) 
                        cout << "Input file " << filename << " is MODIS Terra Geolocation Fields.\n" << endl;
                    
                } else {
                    cerr << "-E-" << __FILE__ << "Line" << __LINE__ << ": Unrecognized input HDF-EOS file" << filename << endl; 

                    return ret;
                }

            /* MISR */
            } else if ((!ncFile.getAtt("Path_number").isNull()) && (!ncFile.getAtt("SOM_parameters.som_ellipsoid.a").isNull())) {
                ret.type = FT_MISR;
                ret.sensor_id = MISR;
            } else {
                fprintf(stderr, "-E- %s Line %d: Unrecognized title %s in input HDF file %s\n", __FILE__, __LINE__, titleStr.c_str(), filename);
                return ret;
            }

        } catch (const NcException &e) {
            std::cerr << "Unknown NetCDF file type - NetCDF error: " << e.what() << std::endl;
            return ret;
        }
        return ret;

    } catch (const NcException &e) {
        // not a netCDF file
    }

    // check for OLCI, MSI and MERIS SAFE format - in case they specified the xml file

    if ((ret = chk_safe_xml(filename)).type != FT_INVALID) {
        return ret;
    } 

    /* SeaBASS Check */
    if ((ret = chk_seabass(filename)).type != FT_INVALID) {
        if (want_verbose)
            printf("Input file %s is a SeaBASS text file.\n", filename);
        return ret;
    }

    /* MERIS Check */
    {
        EPR_SProductId *product_id;

        epr_init_api(e_log_debug, NULL, NULL);
        product_id = epr_open_product(filename);
        if (product_id != NULL) {
            if (product_id->id_string[8] == '1') { /* it is a level 1 file */
                ret.type = FT_MERISL1B;
                ret.sensor_id = MERIS;
                ret.subsensor_id = -1;
                if (want_verbose) {
                    printf("Input file %s is MERIS L1 file.\n", filename);
                }
            }
            epr_close_product(product_id);
            /*remember to close api (epr_close_api();)*/
            return ret;
        }
    }

    // check for OLI
    if ((ret = chk_oli(filename)).type != FT_INVALID) {
        if (want_verbose) {
            printf("Input file %s is a Landsat 8/9 OLI L1B GEOTIFF file.\n", filename);
        }
        return ret;
    }

    // check for AVIRIS
    char hdrfile[FILENAME_MAX], imgfile[FILENAME_MAX], navfile[FILENAME_MAX], gainfile[FILENAME_MAX];
    if ((ret = chk_aviris(filename, hdrfile, imgfile, navfile, gainfile)).type != FT_INVALID) {
        if (want_verbose) {
            printf("Input file %s is an AVIRIS file.\n", filename);
        }
        return ret;
    }

    // check for PRISM
    if ((ret = chk_prism(filename)).type != FT_INVALID) {
        if (want_verbose) 
            printf("Input file %s is a PRISM file.\n", filename);
        return ret;
    }

    // check for Landsat 5 (L5TM)
    if ((ret = chk_l5tm(filename)).type != FT_INVALID) {
        if (want_verbose) 
            printf("Input file %s is a Landsat 5 TM L1B GEOTIFF file.\n", filename);
        return ret;
    }

    // check for Landsat 7 (L7TM)
    if ((ret = chk_l7etm(filename)).type != FT_INVALID) {
        if (want_verbose) 
            printf("Input file %s is a Landsat 7 TM L1B GEOTIFF file.\n", filename);
        return ret;
    }

    return ret;
}


file_format chk_oli(char *filename) {
    /*  ------------------------------------------------------------------------
        chk_oli

        purpose: check a file to see if it is an OLI Landsat8/9 file

        Returns FT_INVALID if not OLI L1B or the format code

        Parameters: (in calling order)
        Type              Name            I/O     Description
        ----              ----            ---     -----------
        char *            filename            I      file to check
        filehandle *      file             I      input file information

        -----------------------------------------------------------------------*/
    file_format ret = {FT_INVALID, -1, -1};
    const int lineSize = 500;
    int i;
    FILE *fp;
    char line[lineSize + 1];
    char *result;

    /* Open file */
    if ((fp = fopen(filename, "re")) == NULL) {
        return ret;
    }

    // skip blank lines
    do {
        result = fgets(line, lineSize, fp);
        if (result == NULL) {
            fclose(fp);
            return ret;
        }
        trimBlanks(line);
    } while (strlen(line) == 0);

    // first line needs to be "GROUP = L1_METADATA_FILE"
    if (strstr(line, "L1_METADATA_FILE") == NULL && strstr(line, "LANDSAT_METADATA_FILE") == NULL) {
        fclose(fp);
        return ret;
    }

    // within 60 lines look for:
    //   SPACECRAFT_ID = "LANDSAT_8/9"  
    //   SENSOR_ID = "OLI"
    int foundSpacecraft = 0;
    int foundSensor = 0;
    int isLandsat8 = 0;
    int isLandsat9 = 0;
    for (i = 0; i < 60; i++) {
        result = fgets(line, lineSize, fp);
        if (result == NULL) {
            fclose(fp);
            return ret;
        }
        if (strstr(line, "SPACECRAFT_ID")) {
            if (strstr(line, "LANDSAT_8")) {
                foundSpacecraft = 1;
                isLandsat8 = 1;
            }
            else if (strstr(line, "LANDSAT_9")) {
                foundSpacecraft = 1;
                isLandsat9 = 1;
            }
        } else if (strstr(line, "SENSOR_ID")) {
            if (strstr(line, "OLI")) {
                foundSensor = 1;
            }
        }

        if (foundSpacecraft && foundSensor) {
            ret.type = FT_OLIL1B;
            if (isLandsat8) {
                ret.subsensor_id = OLI_L8;
                ret.sensor_id = OLIL8;
            } else if (isLandsat9) {
                ret.subsensor_id = OLI_L9;
                ret.sensor_id = OLIL9;
            }
            
            break;
        }
    }

    fclose(fp);
    return ret;
}

file_format chk_oli_geo(char *filename) {
    /*  ------------------------------------------------------------------------
        chk_oli_geo

        purpose: check a file to see if it is an OLI Landsat8/9 GEO file

        Returns FT_INVALID if not OLI L1B or the format code

        Parameters: (in calling order)
        Type              Name            I/O     Description
        ----              ----            ---     -----------
        char *            filename         I      file to check
        filehandle *      file             I      input file information

        -----------------------------------------------------------------------*/
    file_format ret = {FT_INVALID, -1, -1};
    const int lineSize = 500;
    int i;
    FILE *fp;
    char line[lineSize + 1];
    char *result;

    /* Open file */
    if ((fp = fopen(filename, "re")) == NULL) {
        return ret;
    }

    // skip blank lines
    do {
        result = fgets(line, lineSize, fp);
        if (result == NULL) {
            fclose(fp);
            return ret;
        }
        trimBlanks(line);
    } while (strlen(line) == 0);

    // first line needs to be "GROUP = FILE_HEADER"
    if (strstr(line, "FILE_HEADER") == NULL) {
        fclose(fp);
        return ret;
    }

    // within 60 lines look for:
    //   SATELLITE = "LANDSAT_8/9"
    //   BAND_LIST = ...
    int foundSpacecraft = 0;
    int foundBand = 0;
    for (i = 0; i < 60; i++) {
        result = fgets(line, lineSize, fp);
        if (result == NULL) {
            fclose(fp);
            return ret;
        }
        trimBlanks(line);
        if (strstr(line, "SATELLITE") || strstr(line, "SPACECRAFT_ID")) {
            if (strstr(line, "LANDSAT_8") || strstr(line, "LANDSAT_9")) {
                foundSpacecraft = 1;
            }
        } else if (strstr(line, "BAND_LIST")) {
            foundBand = 1;
        }

        if (foundSpacecraft && foundBand) {
            ret.type = FT_OLIL1B;
            break;
        }
    }

    fclose(fp);
    return ret;
}


file_format chk_prism(char *filename) {
    /**
     * @brief check a file to see if it is a PRISM file
     * 
     * @param filename name of file to check
     * @return FT_INVALID or PRISM format code
     * 
     */ 

    file_format ret = {FT_INVALID, -1, -1};
    const int lineSize = 500;
    FILE *fp;
    char line[lineSize + 1];
    char *result;

    if (strstr(filename, "prm") == NULL) {
        return ret;
    }

    /* Open file */
    if ((fp = fopen(filename, "re")) == NULL) {
        return ret;
    }

    // skip blank lines
    do {
        result = fgets(line, lineSize, fp);
        if (result == NULL) {
            fclose(fp);
            return ret;
        }
        trimBlanks(line);
    } while (strlen(line) == 0);

    // first line needs to be "ENVI"
    if (strstr(line, "ENVI") == NULL) {
        fclose(fp);
        return ret;
    }

    ret.type = FT_PRISM;
    ret.sensor_id = PRISM;

    fclose(fp);
    return ret;
}


file_format chk_aviris(char *filename, char *hdrfile, char *imgfile, char *navfile, char *gainfile) {
/**
 * @brief Check a file to see if it is an AVIRIS hdr file
 * 
 * @param filename name of file to check
 * @param hdrfile
 * @param imgfile
 * @param navfile
 * @param gainfile
 * @return FT_INVALID if not the expected AVERIS format code
 * 
*/
    file_format ret = {FT_INVALID, -1, -1};
    const int lineSize = 500;
    int i;
    FILE *fp;
    char line[lineSize + 1];
    char *result;

    if (!checkAvProcessFile(filename, hdrfile, imgfile, navfile, gainfile, FILENAME_MAX)) {
        strncpy(hdrfile, filename, FILENAME_MAX);
    } else {
        ret.type = FT_AVIRIS;
        ret.sensor_id = AVIRIS;
        return ret;
    }

    /* Open file */
    if ((fp = fopen(hdrfile, "re")) == NULL) {
        return ret;
    }

    // skip blank lines
    do {
        result = fgets(line, lineSize, fp);
        if (result == NULL) {
            fclose(fp);
            return ret;
        }
        trimBlanks(line);
    } while (strlen(line) == 0);

    // first line needs to be "ENVI"
    if (strstr(line, "ENVI") == NULL) {
        fclose(fp);
        return ret;
    }

    // within 20 lines look for:
    //   AVIRIS orthocorrected
    int foundAviris = 0;
    int foundFormat = 0;
    int foundOrtho = 0;
    for (i = 0; i < 20; i++) {
        result = fgets(line, lineSize, fp);
        if (result == NULL) {
            fclose(fp);
            return ret;
        }
        trimBlanks(line);
        if (strstr(line, "AVIRIS")) {
            foundAviris = 1;
        }
        if (strstr(line, "orthocorrected")) {
            foundOrtho = 1;
        }
        if (strstr(line, "interleave")) {
            foundFormat = 1;
        }

        if (foundAviris && foundFormat && foundOrtho) {
            ret.type = FT_AVIRIS;
            ret.sensor_id = AVIRIS;
            break;
        }
    }

    fclose(fp);
    return ret;
}

file_format chk_l5tm(char *filename) {
    /**
     * chk_l5tm
     * Purpose: check a file to see if it is a Landsat 5 TM file
     * 
     * @param filename file to check 
     * @return FT_INVALID if not expected format code
     */

    file_format ret = {FT_INVALID, -1, -1};
    const int lineSize = 500;
    int i;
    FILE *fp;
    char line[lineSize + 1];
    char *result;

    /* Open file */
    if ((fp = fopen(filename, "re")) == NULL) {
        return ret;
    }

    // skip blank lines
    do {
        result = fgets(line, lineSize, fp);
        if (result == NULL) {
            fclose(fp);
            return ret;
        }
        trimBlanks(line);
    } while (strlen(line) == 0);

    // first line needs to be "GROUP = L1_METADATA_FILE"
    if (strstr(line, "L1_METADATA_FILE") == NULL) {
        fclose(fp);
        return ret;
    }

    // within 20 lines look for:
    //   SPACECRAFT_ID = "LANDSAT_5"
    //   SENSOR_ID = "TM"
    int foundSpacecraft = 0;
    int foundSensor = 0;
    for (i = 0; i < 20; i++) {
        result = fgets(line, lineSize, fp);
        if (result == NULL) {
            fclose(fp);
            return ret;
        }
        trimBlanks(line);
        if (strstr(line, "SPACECRAFT_ID")) {
            if (strstr(line, "LANDSAT_5")) {
                foundSpacecraft = 1;
            }
        } else if (strstr(line, "SENSOR_ID")){
            if (strstr(line, "TM")) {
                foundSensor = 1;
            }
        }

        if (foundSpacecraft && foundSensor) {
            ret.type = FT_L5TML1B;
            ret.sensor_id = L5TM;
            break;
        }
    }

    fclose(fp);
    return ret;
}

file_format chk_l5tm_geo(char *filename) {
    /*  ------------------------------------------------------------------------
        chk_l5tm_geo

        purpose: check a file to see if it is a Landsat 5 TM GEO file

        Returns FT_INVALID if not LS% L1B or the format code

        Parameters: (in calling order)
        Type              Name            I/O     Description
        ----              ----            ---     -----------
        char *            filename         I      file to check
        filehandle *      file             I      input file information

        -----------------------------------------------------------------------*/
    file_format ret = {FT_INVALID, -1, -1};
    const int lineSize = 500;
    int i;
    FILE *fp;
    char line[lineSize + 1];
    char *result;

    /* Open file */
    if ((fp = fopen(filename, "re")) == NULL) {
        return ret;
    }

    // skip blank lines
    do {
        result = fgets(line, lineSize, fp);
        if (result == NULL) {
            fclose(fp);
            return ret;
        }
        trimBlanks(line);
    } while (strlen(line) == 0);

    // first line needs to be "GROUP = FILE_HEADER"
    if (strstr(line, "FILE_HEADER") == NULL) {
        fclose(fp);
        return ret;
    }

    // within 20 lines look for:
    //   SATELLITE = "LANDSAT_5" for now unconfirmed. Need to confirm this from actual
    //   BAND_LIST = ...   geoletry file if such things exists for Landsat 5.
    int foundSpacecraft = 0;
    int foundBand = 0;
    for (i = 0; i < 20; i++) {
        result = fgets(line, lineSize, fp);
        if (result == NULL) {
            fclose(fp);
            return ret;
        }
        trimBlanks(line);
        if (strstr(line, "SATELLITE")) {
            if (strstr(line, "LANDSAT_5")) {
                foundSpacecraft = 1;
            }
        } else if (strstr(line, "BAND_LIST")) {
            foundBand = 1;
        }

        if (foundSpacecraft && foundBand) {
            ret.type = FT_L5TML1B;
            break;
        }
    }

    fclose(fp);
    return ret;
}

file_format chk_l7etm(char *filename) {
/**
 * chk_l7etm
 * 
 * Purpose: Check a file to see if it is a Landsat 7 file
 * 
 * @param filename file to check
 * @return FT_INVALID if not expected format code
 * 
*/

    file_format ret = {FT_INVALID, -1, -1};
    const int lineSize = 500;
    int i;
    FILE *fp;
    char line[lineSize + 1];
    char *result;

    /* Open file */
    if ((fp = fopen(filename, "re")) == NULL) {
        return ret;
    }

    // skip blank lines
    do {
        result = fgets(line, lineSize, fp);
        if (result == NULL) {
            fclose(fp);
            return ret;
        }
        trimBlanks(line);
    } while (strlen(line) == 0);

    // first line needs to be "GROUP = L1_METADATA_FILE"
    if (strstr(line, "L1_METADATA_FILE") == NULL) {
        fclose(fp);
        return ret;
    }

    // within 20 lines look for:
    //   SPACECRAFT_ID = "LANDSAT_5"
    //   SENSOR_ID = "TM"
    int foundSpacecraft = 0;
    int foundSensor = 0;
    for (i = 0; i < 20; i++) {
        result = fgets(line, lineSize, fp);
        if (result == NULL) {
            fclose(fp);
            return ret;
        }
        trimBlanks(line);
        if (strstr(line, "SPACECRAFT_ID")) {
            if (strstr(line, "LANDSAT_7")) {
                foundSpacecraft = 1;
            }
        } else if (strstr(line, "SENSOR_ID")){
            if (strstr(line, "ETM")) {
                foundSensor = 1;
            }
        }

        if (foundSpacecraft && foundSensor) {
            ret.type = FT_L7ETML1B;
            ret.sensor_id = L7ETMP;
            break;
        }
    }

    fclose(fp);
    return ret;
}

file_format chk_seabass(char *filename) {
/**
 * chk_seabass
 * Purpose: check a file to see if it is a SeaBASS text file
 *
 * @param filename name of file to check
 * @return -1 if not SeaBass or the format code
 * 
*/

    file_format ret = {FT_INVALID, -1, -1};
    FILE *fp;

    if ((fp = fopen(filename, "r")) == NULL) {
        fprintf(stderr, "-E- : input file %s does not exist or is read protected.\n", filename);
        return ret;
    }


    char buffer[2048];
    fgets(buffer, 2048-1, fp);
    if (strncmp(buffer, "/begin_header", 13) == 0) {
        ret.type = FT_SEABASSRRS;

        // Get delimiter
        // while(1) {
        //     if (fgets(buffer, 2048-1, fp) == NULL) {
        //         fprintf(stderr, "-E- : input SeaBASS file %s does not contain delimiter.\n", filename);
        //         fclose(fp);
        //         exit(1);
        //     }
        //     if (strncmp(buffer, "/delimiter=", 11) == 0) {
        //         buffer[strcspn(buffer, "\n")] = 0;
        //         strcpy(file->delimiter, &buffer[11]);
        //         fclose(fp);
        //         break;
        //     }
        // }

        while (1) {
            fgets(buffer, 2048-1, fp);
            if (strncmp(buffer, "/end_header", 13) == 0) {
                break;
            } else if (strncmp(buffer, "/sensor=", 8) == 0) {
                for (int i=8;i<64;i++){
                    if (buffer[i] == '\n'){
                        buffer[i] = '\0';
                        break;
                    }
                }
                ret.sensor_id = sensorName2SensorId(&buffer[8]);
                ret.subsensor_id = sensorId2SubsensorId(ret.sensor_id);
                break;
            }
        } // while loop
    }

    fclose(fp);

    return ret;
}
