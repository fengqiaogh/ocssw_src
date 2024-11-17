#include <netcdf.h>
#include <hdf.h>
#include <hdf5.h>
#include <libgen.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <dfutils.h>
#include <genutils.h>
#include <sensorDefs.h>
#include <sensorInfo.h>

#include <hdf.h>
#include <mfhdf.h>

#include "aviris.h"
#include "epr_api.h"

#include "filetype.h"

#define EOSMETALEN 32768

static file_format chk_oli(char *filename);
static file_format chk_l5tm(char *filename);
static file_format chk_l7etm(char *filename);
static file_format chk_aviris(char *filename, char *hdrfile, char *imgfile, char *navfile, char *gainfile);
static file_format chk_prism(char *filename);
static file_format chk_seabass(char *filename);

// defined in filetypeXml.cpp
file_format chk_safe_xml(char *filename);


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
        if(strncmp(buf, val, strlen(val)) == 0) {
            return 1;
        }
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
                            if (want_verbose) {
                                printf("Input file %s is a HICO L1B HDF 5 file.\n", filename);
                            }
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
                    if (want_verbose) {
                        printf("Input file %s is a GOCI L1B HDF 5 file.\n", filename);
                    }
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
                                if (want_verbose) {
                                    printf("Input file %s is a SGLI L1B HDF 5 file.\n", filename);
                                }
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

file_format getFormat(char *filename) {
    file_format ret = { FT_INVALID, -1, -1 };

    int32_t sd_id;
    char eosmeta[EOSMETALEN] = "";
    char tempstr[32] = "";
    idDS ds_id;

    /* Does the file exist? */
    if (access(filename, F_OK) || access(filename, R_OK)) {
        printf("-E- %s: Input file '%s' does not exist or cannot open.\n", __FILE__, filename);
        return ret;
    }

    // need to do HDF5 files before trying NetCDF
    if ((ret = chk_hdf5(filename)).type != FT_INVALID) {
        return ret;
    }

    /* Is it netCDF */
    ds_id = openDS(filename);
    if (ds_id.fid != FAIL) {
        if (ds_id.fftype == DS_NCDF) {
            char *titleStr = readAttrStr(ds_id, "title");

            if (titleStr) {
                lowcase(titleStr);
                if (strstr(titleStr, "viirs level-1a")) {
                    char *platformStr = readAttrStr(ds_id, "platform");
                    if (platformStr) {
                        lowcase(platformStr);
                        if (strstr(platformStr, "suomi-npp")) {
                            ret.type = FT_VIIRSL1A;
                            ret.sensor_id = VIIRSN;
                            ret.subsensor_id = VIIRS_NPP;
                            if (want_verbose) {
                                printf("Input file %s is VIIRS NPP L1A NetCDF4.\n", filename);
                            }
                        } else if (strstr(platformStr, "jpss-1")) {
                            ret.type = FT_VIIRSL1A;
                            ret.sensor_id = VIIRSJ1;
                            ret.subsensor_id = VIIRS_J1;
                            if (want_verbose) {
                                printf("Input file %s is VIIRS JPSS-1 L1A NetCDF4.\n", filename);
                            }
                        } else if (strstr(platformStr, "jpss-2")) {
                            ret.type = FT_VIIRSL1A;
                            ret.sensor_id = VIIRSJ2;
                            ret.subsensor_id = VIIRS_J2;
                            if (want_verbose) {
                                printf("Input file %s is VIIRS JPSS-2 L1A NetCDF4.\n", filename);
                            }
                        }
                        free(platformStr);
                    }
                    free(titleStr);
                    endDS(ds_id);
                    return ret;
                }
                if (strstr(titleStr, "viirs m-band")) {
                    if (strstr(titleStr, "viirs m-band reflected solar band")) {
                        ret.type = FT_VIIRSL1BNC;
                        if (want_verbose) {
                            printf("Input file %s is VIIRS NPP L1B NetCDF4.\n", filename);
                        }
                    } else if (strstr(titleStr, "viirs m-band geolocation data")) {
                        ret.type = FT_VIIRSGEONC;
                        if (want_verbose) {
                            printf("Input file %s is VIIRS NPP GEO NetCDF4.\n", filename);
                        }
                    }
                    char *platformStr = readAttrStr(ds_id, "platform");
                    if (platformStr) {
                        lowcase(platformStr);
                        if (strstr(platformStr, "suomi-npp")) {
                            ret.sensor_id = VIIRSN;
                            ret.subsensor_id = VIIRS_NPP;
                        } else if (strstr(platformStr, "jpss-1")) {
                            ret.sensor_id = VIIRSJ1;
                            ret.subsensor_id = VIIRS_J1;
                        } else if (strstr(platformStr, "jpss-2")) {
                            ret.sensor_id = VIIRSJ2;
                            ret.subsensor_id = VIIRS_J2;
                        }
                        free(platformStr);
                    }
                    free(titleStr);
                    endDS(ds_id);
                    return ret;
                }
                if (strstr(titleStr, "meris l1b")) {
                    ret.type = FT_MERISCC;
                    ret.sensor_id = MERIS;
                    ret.subsensor_id = -1;
                    if (want_verbose) {
                        printf("Input file %s is MERIS CC file.\n", filename);
                    }
                    free(titleStr);
                    endDS(ds_id);
                    return ret;
                }
                if (strstr(titleStr, "ocia level-1b")) {
                    ret.type = FT_OCIA;
                    ret.sensor_id = OCIA;
                    ret.subsensor_id = -1;
                    if (want_verbose) {
                        printf("Input file %s is OCIA L1B file.\n", filename);
                    }
                    free(titleStr);
                    endDS(ds_id);
                    return ret;
                }
                if (strstr(titleStr, "aviris level-1b")) {
                    ret.type = FT_AVIRIS;
                    ret.sensor_id = AVIRIS;
                    ret.subsensor_id = -1;
                    if (want_verbose) {
                        printf("Input file %s is AVIRIS L1B NetCDF4 file.\n", filename);
                    }
                    free(titleStr);
                    endDS(ds_id);
                    return ret;
                }
                if (strstr(titleStr, "olci level 1b product")) {
                    char subtype[4];
                    if (strstr(titleStr, "geo coordinates")) {
                        ret.type = FT_OLCIGEO;
                        strcpy(subtype,"GEO");
                    } else {
                        ret.type = FT_OLCI;
                        strcpy(subtype,"L1B");
                    }
                    char *productNameStr = readAttrStr(ds_id, "product_name");
                    if(!strncmp(productNameStr, "S3A", 3)) {
                        ret.sensor_id = OLCIS3A;
                        ret.subsensor_id = OLCI_S3A;
                        if (want_verbose) {
                            printf("Input file %s is OLCI S3A %s file.\n", filename, subtype);
                        }
                        free(productNameStr);
                        free(titleStr);
                        endDS(ds_id);
                        return ret;
                    } else if(!strncmp(productNameStr, "S3B", 3)) {
                        ret.sensor_id = OLCIS3B;
                        ret.subsensor_id = OLCI_S3B;
                        if (want_verbose) {
                            printf("Input file %s is OLCI S3B %s file.\n", filename, subtype);
                        }
                        free(productNameStr);
                        free(titleStr);
                        endDS(ds_id);
                        return ret;
                    }
                    free(productNameStr);
                }
                if (strstr(titleStr, "hawkeye level-1a")) {
                    ret.type = FT_HAWKEYEL1A;
                    ret.sensor_id = HAWKEYE;
                    ret.subsensor_id = -1;
                    if (want_verbose) {
                        printf("Input file %s is HAWKEYE L1A file.\n", filename);
                    }
                    free(titleStr);
                    endDS(ds_id);
                    return ret;
                }
                if (strstr(titleStr, "oci level-1b")) {
                    ret.type = FT_OCIL1B;
                    ret.sensor_id = OCI;
                    ret.subsensor_id = -1;
                    if (want_verbose) {
                        printf("Input file %s is PACE L1B file.\n", filename);
                    }
                    free(titleStr);
                    endDS(ds_id);
                    return ret;
                }
                if (strstr(titleStr, "ocis level-1b")) {
                    ret.type = FT_OCIS;
                    ret.sensor_id = OCIS;
                    ret.subsensor_id = -1;
                    if (want_verbose) {
                        printf("Input file %s is PACE L1B Simulated file.\n", filename);
                    }
                    free(titleStr);
                    endDS(ds_id);
                    return ret;
                }

                if (strstr(titleStr, "hkt level-1a") ||
                    strstr(titleStr, "pace hkt data")) {
                    ret.type = FT_HKT;
                    ret.sensor_id = OCIS;
                    ret.subsensor_id = -1;
                    if (want_verbose) {
                        printf("Input file %s is PACE HKT file.\n", filename);
                    }
                    free(titleStr);
                    endDS(ds_id);
                    return ret;                
                }
                if (strstr(titleStr, "pace oci level-1c")) {
                    ret.type = FT_L1C;
                    ret.sensor_id = OCI;
                    ret.subsensor_id = -1;
                    if (want_verbose) {                       
                        printf("Input file %s is PACE OCI L1C file.\n", filename);
                    }
                    free(titleStr);
                    endDS(ds_id);
                    return ret;
                }
                if (strstr(titleStr, "pace ocis level-1c")) {
                    ret.type = FT_L1C;
                    ret.sensor_id = OCIS;
                    ret.subsensor_id = -1;
                    if (want_verbose) {                       
                        printf("Input file %s is PACE OCIS L1C file.\n", filename);
                    }
                    free(titleStr);
                    endDS(ds_id);
                    return ret;
                }
                if (strstr(titleStr, "pace harp2 level-1c")) {
                    ret.type = FT_L1C;
                    ret.sensor_id = HARP2;
                    ret.subsensor_id = -1;
                    if (want_verbose) {                       
                        printf("Input file %s is PACE HARP2 L1C file.\n", filename);
                    }
                    free(titleStr);
                    endDS(ds_id);
                    return ret;
                }
                if (strstr(titleStr, "pace spexone level-1c")) {
                    ret.type = FT_L1C;
                    ret.sensor_id = SPEXONE;
                    ret.subsensor_id = -1;
                    if (want_verbose) {                       
                        printf("Input file %s is PACE SPEXone L1C file.\n", filename);
                    }
                    free(titleStr);
                    endDS(ds_id);
                    return ret;
                }
                if (strstr(titleStr, "l1c ancillary file")) {
                    ret.type = FT_L1CANC;
                    ret.sensor_id = OCI;
                    ret.subsensor_id = -1;
                    if (want_verbose) {
                        printf("Input file %s is L1C ancillary file.\n", filename);
                    }
                    free(titleStr);
                    endDS(ds_id);
                    return ret;
                }
                if (strstr(titleStr, "spexone level-1b")) {
                    ret.type = FT_SPEXONE;
                    ret.sensor_id = SPEXONE;
                    ret.subsensor_id = -1;
                    if (want_verbose) {
                        printf("Input file %s is PACE SPEXone file.\n", filename);
                    }
                    free(titleStr);
                    endDS(ds_id);
                    return ret;
                }

                if (strstr(titleStr, "harp2 level-1b")) {
                    ret.type = FT_HARP2;
                    ret.sensor_id = HARP2;
                    ret.subsensor_id = -1;
                    if (want_verbose) {
                        printf("Input file %s is PACE HARP2 file.\n", filename);
                    }
                    free(titleStr);
                    endDS(ds_id);
                    return ret;
                }

                // NetCDF for Seawifs L1A
                if (strstr(titleStr, "seawifs level-1a data")) {
                    ret.type = FT_SEAWIFSL1ANC;
                    ret.sensor_id = SEAWIFS;

                    char *dataTypeStr = readAttrStr(ds_id, "data_type");
                    if (dataTypeStr) {
                        if (strcmp(dataTypeStr, "GAC") == 0) {
                            ret.subsensor_id = SEAWIFS_GAC;
                            if (want_verbose) {
                                printf("Input file %s is SeaWiFS Level-1A GAC.\n", filename);
                            }
                        }
                        else {
                            ret.subsensor_id = SEAWIFS_LAC;
                            if (want_verbose) {
                                printf("Input file %s is SeaWiFS Level-1A LAC.\n", filename);
                            }
                        }
                    }
                    free(dataTypeStr);
                    free(titleStr);
                    endDS(ds_id);
                    return ret;
                }
                if (strstr(titleStr, "octs level-1a gac data")) {
                    ret.type = FT_OCTSL1ANC;
                    ret.sensor_id = OCTS;
                    ret.subsensor_id = -1;
                    if (want_verbose) {
                        printf("Input file %s is OCTS Level-1A GAC netCDF.\n", filename);
                    }
                    free(titleStr);
                    endDS(ds_id);
                    return ret;
                }

                char *platformStr = readAttrStr(ds_id, "platform");

                if (platformStr) {
                    char *instrumentStr = readAttrStr(ds_id, "instrument");
                    if (instrumentStr) {
                        char *processingLevelStr = readAttrStr(ds_id, "processing_level");
                        if (processingLevelStr) {
                            if (!strcmp(processingLevelStr, "L1B")) {
                                ret.type = FT_L1BNCDF;
                            } else if (!strcmp(processingLevelStr, "L2")) {
                                ret.type = FT_L2NCDF;
                            } else if (!strcmp(processingLevelStr, "L3 Binned")) {
                                ret.type = FT_L3BIN;
                            } else if (!strcmp(processingLevelStr, "L3 Mapped")) {
                                ret.type = FT_L3MAP;
                            }
                            if (ret.type != FT_INVALID) {
                                ret.sensor_id = instrumentPlatform2SensorId(instrumentStr, platformStr);
                                ret.subsensor_id = sensorId2SubsensorId(ret.sensor_id);
                                if (want_verbose) {
                                    printf("Input file %s is a NetCDF4 %s %s file.\n", filename, instrumentStr, processingLevelStr);
                                }
                                free(processingLevelStr);
                                free(instrumentStr);
                                free(platformStr);
                                free(titleStr);
                                endDS(ds_id);
                                return ret;
                            }
                            free(processingLevelStr);
                        } // processingLevel found
                        free(instrumentStr);
                    } // instrument found
                    free(platformStr);
                } // platform found
                free(titleStr);
            } // title found

            // unknown netCDF file
            return ret;
        } // is a NetCDF file
        endDS(ds_id);
    } // data set opened successfully

    /* Is it HDF? */
    sd_id = SDstart(filename, DFACC_RDONLY);
    if (sd_id != FAIL) {

        /* File is HDF. Is it one of ours? */

        char title[255];
        char sensor[80];
        if (SDreadattr(sd_id, SDfindattr(sd_id, "Title"), (VOIDP) title) == 0) {
            if (strstr(title, "Level-3 Binned Data") != NULL) {
                if (SDreadattr(sd_id, SDfindattr(sd_id, "Sensor Name"), (VOIDP) sensor) == 0) {

                    // kludge for VIIRS EDR L3
                    if (strcmp(sensor, "VIIRS") == 0) {
                        strncpy(sensor, "VIIRSN", strlen("VIIRSN")+1);
                    }

                    if ((ret.sensor_id = sensorName2SensorId(sensor)) != FT_INVALID) {
                        ret.subsensor_id = sensorId2SubsensorId(ret.sensor_id);
                        ret.type = FT_L3BIN;
                        if (want_verbose) {
                            printf("Input file %s is %s.\n", filename, title);
                        }
                    } else {
                        fprintf(stderr, "-E- %s Line %d: Unknown sensor name in Level-3 file %s\n", __FILE__, __LINE__, filename);
                        return ret;
                    }
                } else {
                    fprintf(stderr, "-E- %s Line %d: No sensor name attribute in Level-3 file %s\n", __FILE__, __LINE__, filename);
                    return ret;
                }
            } else if (strstr(title, "Level-2 Data") != NULL) {
                if (SDreadattr(sd_id, SDfindattr(sd_id, "Sensor Name"), (VOIDP) sensor) == 0) {
                    if ((ret.sensor_id = sensorName2SensorId(sensor)) != FT_INVALID) {
                        ret.subsensor_id = sensorId2SubsensorId(ret.sensor_id);
                        ret.type = FT_L2HDF;
                        if (want_verbose) {
                            printf("Input file %s is %s.\n", filename, title);
                        }
                    } else {
                        fprintf(stderr, "-E- %s Line %d: Unknown sensor name in Level-2 file %s\n", __FILE__, __LINE__, filename);
                        return ret;
                    }
                } else {
                    fprintf(stderr, "-E- %s Line %d: No sensor name attribute in Level-2 file %s\n", __FILE__, __LINE__, filename);
                    return ret;
                }
            } else if (strcmp(title, "SeaWiFS Level-1A Data") == 0) {
                ret.type = FT_SEAWIFSL1A;
                ret.sensor_id = SEAWIFS;
                if (SDreadattr(sd_id, SDfindattr(sd_id, "Data Type"), (VOIDP) tempstr) == 0) {
                    if (strcmp(tempstr, "GAC") == 0) {
                        ret.subsensor_id = SEAWIFS_GAC;
                        if (want_verbose) {
                            printf("Input file %s is SeaWiFS Level-1A GAC.\n", filename);
                        }
                    } else if (strcmp(tempstr, "LAC") == 0) {
                        ret.subsensor_id = SEAWIFS_LAC;
                        if (want_verbose) {
                            printf("Input file %s is SeaWiFS Level-1A LAC.\n", filename);
                        }
                    } else {
                        ret.subsensor_id = SEAWIFS_LAC;
                        if (want_verbose) {
                            printf("Input file %s is assumed to be SeaWiFS Level-1A LAC.\n", filename);
                        }
                    }
                } else {
                    ret.subsensor_id = SEAWIFS_LAC;
                    if (want_verbose) {
                        printf("Input file %s is assumed to be SeaWiFS Level-1A LAC.\n", filename);
                    }
                }
            } else if (strcmp(title, "OCTS Level-1A GAC Data") == 0) {
                ret.type = FT_OCTSL1A;
                ret.sensor_id = OCTS;
                ret.subsensor_id = -1;
                if (want_verbose) {
                    printf("Input file %s is %s.\n", filename, title);
                }
            } else if (strcmp(title, "OSMI Level-1A Data") == 0) {
                ret.type = FT_OSMIL1A;
                ret.sensor_id = OSMI;
                ret.subsensor_id = -1;
                if (want_verbose) {
                    printf("Input file %s is %s.\n", filename, title);
                }
            } else if (strcmp(title, "CZCS Level-1A Data") == 0) {
                ret.type = FT_CZCSL1A;
                ret.sensor_id = CZCS;
                ret.subsensor_id = -1;
                if (want_verbose) {
                    printf("Input file %s is %s.\n", filename, title);
                }
            } else if (strcmp(title, "OCM1 Level-1B (OBPG)") == 0) {
                ret.type = FT_OCML1B;
                ret.sensor_id = OCM1;
                ret.subsensor_id = -1;
                if (want_verbose) {
                    printf("Input file %s is %s.\n", filename, title);
                }
            } else if (strncmp(title, "OCM Level-1B", 12) == 0) {
                ret.type = FT_OCML1BDB;
                ret.sensor_id = OCM1;
                ret.subsensor_id = -1;
                if (want_verbose) {
                    printf("Input file %s is %s.\n", filename, title);
                }
            } else if (strcmp(title, "Oceansat OCM2 Level-1B Data") == 0) {
                ret.type = FT_OCM2L1B;
                ret.sensor_id = OCM2;
                ret.subsensor_id = -1;
                if (want_verbose) {
                    printf("Input file %s is %s.\n", filename, title);
                }

                /* generic L1B format support */
            } else if (strcmp(title, "SeaWiFS Level-1B") == 0) {
                ret.type = FT_L1HDF;
                ret.sensor_id = SEAWIFS;
                ret.subsensor_id = -1;
                if (want_verbose) {
                    printf("Input file %s is %s.\n", filename, title);
                }

            } else if (strcmp(title, "MERIS Level-1B") == 0) {
                ret.type = FT_L1HDF;
                ret.sensor_id = MERIS;
                ret.subsensor_id = -1;
                if (want_verbose) {
                    printf("Input file %s is %s.\n", filename, title);
                }
            } else if (strcmp(title, "VIIRS Level-1B") == 0 || strcmp(title, "VIIRSN Level-1B") == 0) {
                ret.type = FT_L1HDF;
                ret.sensor_id = VIIRSN;
                ret.subsensor_id = VIIRS_NPP;
                if (want_verbose) {
                    printf("Input file %s is %s.\n", filename, title);
                }
            } else if (strcmp(title, "VIIRSJ1 Level-1B") == 0) {
                ret.type = FT_L1HDF;
                ret.sensor_id = VIIRSJ1;
                ret.subsensor_id = VIIRS_J1;
                if (want_verbose) {
                    printf("Input file %s is %s.\n", filename, title);
                }
            } else if (strcmp(title, "VIIRSJ2 Level-1B") == 0) {
                ret.type = FT_L1HDF;
                ret.sensor_id = VIIRSJ2;
                ret.subsensor_id = VIIRS_J2;
                if (want_verbose) {
                    printf("Input file %s is %s.\n", filename, title);
                }
            } else if (strcmp(title, "OCM2 Level-1B") == 0) {
                ret.type = FT_L1HDF;
                ret.sensor_id = OCM2;
                ret.subsensor_id = -1;
                if (want_verbose) {
                    printf("Input file %s is %s.\n", filename, title);
                }
            } else if (strcmp(title, "OCTS Level-1B") == 0) {
                ret.type = FT_L1HDF;
                ret.sensor_id = OCTS;
                ret.subsensor_id = -1;
                if (want_verbose) {
                    printf("Input file %s is %s.\n", filename, title);
                }
            } else if (strcmp(title, "MOS Level-1B") == 0) {
                ret.type = FT_L1HDF;
                ret.sensor_id = MOS;
                ret.subsensor_id = -1;
                if (want_verbose) {
                    printf("Input file %s is %s.\n", filename, title);
                }
            } else if (strcmp(title, "OCTS Level-1B LAC Data") == 0) {
                ret.type = FT_OCTSL1B;
                ret.sensor_id = OCTS;
                ret.subsensor_id = -1;
                if (want_verbose) {
                    printf("Input file %s is %s.\n", filename, title);
                }
            } else if (strcmp(title, "OSMI Level-1B") == 0) {
                ret.type = FT_L1HDF;
                ret.sensor_id = OSMI;
                ret.subsensor_id = -1;
                if (want_verbose) {
                    printf("Input file %s is %s.\n", filename, title);
                }
            } else if (strcmp(title, "HMODIST Level-1B") == 0 || strcmp(title, "MODIST Level-1B") == 0) {
                ret.type = FT_L1HDF;
                ret.sensor_id = MODIST;
                ret.subsensor_id = MODIS_TERRA;
                if (want_verbose) {
                    printf("Input file %s is %s.\n", filename, title);
                    printf("\n\n *** WARNING: polarization can not be computed from generic format ***\n\n");
                }
            } else if (strcmp(title, "HMODISA Level-1B") == 0 || strcmp(title, "MODISA Level-1B") == 0) {
                ret.type = FT_L1HDF;
                ret.sensor_id = MODISA;
                ret.subsensor_id = MODIS_AQUA;
                if (want_verbose) {
                    printf("Input file %s is %s.\n", filename, title);
                    printf("\n\n *** WARNING: polarization can not be computed from generic format ***\n\n");
                }
            } else if (strcmp(title, "CZCS Level-1B") == 0) {
                ret.type = FT_L1HDF;
                ret.sensor_id = CZCS;
                ret.subsensor_id = -1;
                if (want_verbose) {
                    printf("Input file %s is %s.\n", filename, title);
                }

            } else if (strstr(title, "Level-1 cross-calibration pixels") != NULL) {
                if (SDreadattr(sd_id, SDfindattr(sd_id, "sensorID"), (VOIDP) & (ret.sensor_id)) != 0) {
                    fprintf(stderr, "-E- %s Line %d: Unrecognized sensor name, title %s in input HDF file %s\n", __FILE__, __LINE__, title, filename);
                    return ret;
                }
                ret.type = FT_L1XCAL;
                if (want_verbose) {
                    printf("Input file %s is %s.\n", filename, title);
                }
            } else if (strstr(title, "AVHRR") != NULL) {
                ret.type = FT_CLASSAVHRR;
                ret.sensor_id = AVHRR;
                ret.subsensor_id = -1;
                if (want_verbose) {
                    printf("Input file %s is %s.\n", filename, title);
                }
            } else {
                fprintf(stderr, "-E- %s Line %d: Unrecognized title %s in input HDF file %s\n", __FILE__, __LINE__, title, filename);
                return ret;
            }

        } else if (SDreadattr(sd_id, SDfindattr(sd_id, "title"), (VOIDP) title) == 0) {

            if (strstr(title, "AVHRR") != NULL) {
                ret.type = FT_CLASSAVHRR;
                ret.sensor_id = AVHRR;
                ret.subsensor_id = -1;
                if (want_verbose) {
                    printf("Input file %s is %s.\n", filename, title);
                }
            }

        } else if (SDreadattr(sd_id, SDfindattr(sd_id, "satellite"), (VOIDP) title) == 0) {

            if (strstr(title, "oceansat-1") != NULL) {
                ret.type = FT_OCML1BDB;
                ret.sensor_id = OCM1;
                ret.subsensor_id = -1;
                if (want_verbose) {
                    printf("Input file %s is %s.\n", filename, "OCM1 DB file");
                }
            }

            /* Is it HDF-EOS L1B format? */

        } else if (SDreadattr(sd_id, SDfindattr(sd_id, "ArchiveMetadata.0"), (VOIDP) eosmeta) == 0) {

            if (strstr(eosmeta, "MODIS/Terra Calibrated Radiances 5-Min L1B Swath 1km") != NULL) {
                ret.type = FT_MODISL1B;
                ret.sensor_id = MODIST;
                ret.subsensor_id = MODIS_TERRA;

                if (want_verbose) {
                    printf("Input file %s is MODIS Terra Level-1B HDF-EOS product.\n", filename);
                }
                return ret;
            } else if (strstr(eosmeta, "MODIS/Aqua Calibrated Radiances 5-Min L1B Swath 1km") != NULL) {
                ret.type = FT_MODISL1B;
                ret.sensor_id = MODISA;
                ret.subsensor_id = MODIS_AQUA;

                if (want_verbose) {
                    printf("Input file %s is MODIS Aqua Level-1B HDF-EOS product.\n", filename);
                }
                return ret;
            } else if (strstr(eosmeta, "MODIS/Aqua Calibrated Radiances 5-Min L1B Swath 250m") != NULL) {
                printf("Input file %s is MODIS Aqua Level-1B HDF-EOS product at 250m Resolution.\nThis file is no longer accepted as input.\nPlease use the 1KM (LAC) file and set resolution=250\n", filename);
                return ret;
            } else if (strstr(eosmeta, "MODIS/Aqua Calibrated Radiances 5-Min L1B Swath 500m") != NULL) {
                printf("Input file %s is MODIS Aqua Level-1B HDF-EOS product at 500m Resolution.\nThis file is no longer accepted as input.\nPlease use the 1KM (LAC) file and set resolution=500\n", filename);
                return ret;
            } else if (strstr(eosmeta, "MODIS/Terra Calibrated Radiances 5-Min L1B Swath 250m") != NULL) {
                printf("Input file %s is MODIS Terra Level-1B HDF-EOS product at 250m Resolution.\nThis file is no longer accepted as input.\nPlease use the 1KM (LAC) file and set resolution=250\n", filename);
                return ret;
            } else if (strstr(eosmeta, "MODIS/Terra Calibrated Radiances 5-Min L1B Swath 500m") != NULL) {
                printf("Input file %s is MODIS Terra Level-1B HDF-EOS product at 500m Resolution.\nThis file is no longer accepted as input.\nPlease use the 1KM (LAC) file and set resolution=500\n", filename);
                return ret;
            } else if (strstr(eosmeta, "MODIS/Aqua Geolocation Fields") != NULL) {
                ret.type = FT_MODISGEO;
                ret.sensor_id = MODISA;
                ret.subsensor_id = MODIS_AQUA;
                if (want_verbose) {
                    printf("Input file %s is MODIS Aqua Geolocation Fields.\n", filename);
                }
            } else if (strstr(eosmeta, "MODIS/Terra Geolocation Fields") != NULL) {
                ret.type = FT_MODISGEO;
                ret.sensor_id = MODIST;
                ret.subsensor_id = MODIS_TERRA;
                if (want_verbose) {
                    printf("Input file %s is MODIS Terra Geolocation Fields.\n", filename);
                }
            } else {
                fprintf(stderr, "-E- %s Line %d: Unrecognized HDF-EOS file %s\n", __FILE__, __LINE__, filename);
                return ret;
            }

            /* MISR */
        } else if ((SDfindattr(sd_id, "Path_number") != -1) && (SDfindattr(sd_id, "SOM_parameters.som_ellipsoid.a") != -1)) {
          ret.type = FT_MISR;
          ret.sensor_id = MISR;


            /* Is it MOS L1B HDF standard product? */

        } else if (GetFileDesc(filename) != NULL) {
            ret.type = FT_MOSL1B;
            ret.sensor_id = MOS;
            ret.subsensor_id = -1;
            if (want_verbose) {
                printf("Input file %s is MOS Level-1B standard product.\n", filename);
            }
        } else {
            fprintf(stderr, "-E- %s Line %d: Unrecognized input HDF file %s\n", __FILE__, __LINE__, filename);
            return ret;
        }

        SDend(sd_id);
    } else {

        // check for OLCI, MSI and MERIS SAFE format - in case they specified the xml file
        if ((ret = chk_safe_xml(filename)).type != FT_INVALID) {
            return ret;
        }

        /* check for SeaBASS? */
        if ((ret = chk_seabass(filename)).type != FT_INVALID) {
            if (want_verbose)
                printf("Input file %s is a SeaBASS text file.\n", filename);
            return ret;
        }

        /* Is it MERIS? */
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
            if (want_verbose) {
                printf("Input file %s is a PRISM file.\n", filename);
            }
            return ret;
        }

        // check for Landsat 5 (L5TM)
        if ((ret = chk_l5tm(filename)).type != FT_INVALID) {
            if (want_verbose) {
                printf("Input file %s is a Landsat 5 TM L1B GEOTIFF file.\n", filename);
            }
            return ret;
        }

        // check for Landsat 7 (L7TM)
        if ((ret = chk_l7etm(filename)).type != FT_INVALID) {
            if (want_verbose) {
                printf("Input file %s is a Landsat 7 TM L1B GEOTIFF file.\n", filename);
            }
            return ret;
        }
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

/**
     chk_prism
     purpose: check a file to see if it is PRISM hdr file
     Returns FT_INVALID if not PRISM or the format code

    @param char *filename      - name of file to check
    @param filehandle *file - input file information
*/

file_format chk_prism(char *filename) {
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

/**
     chk_aviris
     purpose: check a file to see if it is an AVIRIS hdr file
     Returns FT_INVALID if not AVIRIS or the format code

    @param char *filename      - name of file to check
    @param filehandle *file - input file information
*/

file_format chk_aviris(char *filename, char *hdrfile, char *imgfile, char *navfile, char *gainfile) {
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
    /*  ------------------------------------------------------------------------
        chk_oli

        purpose: check a file to see if it is an OLI Landsat8 file

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
    /*  ------------------------------------------------------------------------
        chk_oli

        purpose: check a file to see if it is an OLI Landsat8 file

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

/**
 *   chk_seabass
 *   purpose: check a file to see if it is a SeaBASS text file
 *   Returns -1 if not SeaBass or the format code
 *
 * @param char *filename      - name of file to check
 */

file_format chk_seabass(char *filename) {

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
