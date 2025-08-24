#include "l1afile_oci.hpp"
#include <genutils.h>
#include <algorithm>

using namespace std;
using namespace netCDF;
using namespace netCDF::exceptions;

// Chunking the arrays TODO: move to class as a const var

#define CHUNKBANDS 36
#define CHUNKPIXELS 256
#define CHUNKLINES 128

/**
 * Class is used in l1agen_oci.cpp to generate a Level 1 file.
 */

// FIXME
// // we need to use the timeutils function ccsds_to_yds() instead
// int ccsds_sec_to_yds(uint8_t *cctime, int32_t *iyear, int32_t *iday, double *sec) {
//     uint32_t ui32;
//     uint32_t ccsec;

//     memcpy(&ui32, cctime, 4);
//     ccsec = SWAP_4(ui32);

//     double dsec = (double)ccsec;
//     int leap = leapseconds_since_1993(dsec);
//     ccsec -= (leap + 27);

//     *iday = ccsec / 86400;
//     int32_t jday = *iday + 2436205;  // Jan. 1, 1958 is Julian day 2436205
//     jdate(jday, iyear, iday);

//     // Get milliseconds
//     int32_t msec = cctime[4] * 256 + cctime[5];
//     int32_t isec = ccsec % 86400;
//     *sec = isec + msec / 65536.0;

//     return 0;
// }

/**
 * @brief Convert the secondary packet header of a PACE OCI file from TAI time to unix time
 * @param cctime An array of 8 bytes containing the ccsds header
 * @return unix time of this packet
 */
double ccsds_sec_to_unix(uint8_t *cctime) {

    // Get values from memory, swap endianness
    uint32_t taiTime;  // Little-endian version
    swapc_bytes2((char *)cctime, (char *)&taiTime, sizeof(taiTime), 1);

    // convert to unix time
    double utime = tai58_to_unix(taiTime);

    // Get milliseconds of this packet
    int32_t msec = (cctime[4] << 8) + cctime[5];

    // add in miliseconds
    utime += msec / 65536.0;

    return utime;
}


L1aFile::L1aFile() {}

L1aFile::~L1aFile() {}

void L1aFile::freeFile() {
    delete l1afile;
    l1afile = nullptr;
}

string L1aFile::getFileName() {
    return fileName;
}

// write all default NASA global data into the l1afile. ie: creator_name and email
// also write what commands were used to make the l1a file
int L1aFile::writeGlobalAttributes(std::string history, std::string doi, std::string pversion) {
    set_global_attrs(this->l1afile, history, doi, pversion);

    return 0;
}

int L1aFile::findNavigationIndex(string startTimeStr, double startTimeInUnixSecs, double endTimeInUnixSecs,
                                 double *navigationTimeData, size_t navigationDataSize, int &startIndex,
                                 int &endIndex) {
    double currTimeInUnix;
    int year, month, day;
    sscanf(startTimeStr.c_str(), "%4d-%2d-%2d", &year, &month, &day);

    for (size_t i = 0; i < navigationDataSize; i++) {
        currTimeInUnix = ymds2unix(year, month, day, navigationTimeData[i]);
        if ((currTimeInUnix > startTimeInUnixSecs - 10) && (startIndex == 1e6)) {
            startIndex = i;
            break;
        }
    }

    for (size_t i = navigationDataSize - 1; i >= 0; i--) {
        currTimeInUnix = ymds2unix(year, month, day, navigationTimeData[i]);
        if ((currTimeInUnix < endTimeInUnixSecs + 10) && (endIndex == -1)) {
            endIndex = i;
            break;
        }
    }
    return 0;
}

inline int expandEnvVar(std::string *sValue) {
    if ((*sValue).find_first_of("$") == std::string::npos)
        return 0;
    std::string::size_type posEndIdx = (*sValue).find_first_of("/");
    if (posEndIdx == std::string::npos)
        return 0;
    const std::string envVar = sValue->substr(1, posEndIdx - 1);
    char *envVar_str = getenv(envVar.c_str());
    if (envVar_str == 0x0) {
        printf("Environment variable: %s not defined.\n", sValue->c_str());
        exit(1);
    }
    *sValue = envVar_str + (*sValue).substr(posEndIdx);

    return 0;
}

void L1aFile::applyChunkingAndCompressionToSciData(netCDF::NcVar &sciVar,
                                                   std::vector<netCDF::NcDim> &currDims) {
    // if the current dim's are smaller, then use them instead of the "CHUNK_" ones.
    vector<size_t> chunkedSizes = {CHUNK_LINES, CHUNK_BANDS, CHUNK_PIXELS};  // default chunk size
    size_t currNumBands = currDims[1].getSize();
    size_t currNumPixels = currDims[2].getSize();
    if (currNumBands < CHUNK_BANDS) {
        chunkedSizes[1] = currNumBands;
    }
    if (currNumPixels < CHUNK_PIXELS) {
        chunkedSizes[2] = currNumPixels;
    }
    // skipping lines since lines will always be larger

    bool chunkingDone = false;
    bool compressionDone = false;

    try {
        sciVar.setChunking(sciVar.nc_CHUNKED, chunkedSizes);
        chunkingDone = true;

        // default compression to 5
        sciVar.setCompression(true, true, 5);
        compressionDone = true;
    } catch (NcException &e) {
        if (!chunkingDone) {
            cerr << "Error setting chunking for " << sciVar.getName() << endl;
        }

        // report compression error only
        if (chunkingDone && !compressionDone) {
            cerr << "Error setting chunking for " << sciVar.getName() << endl;
        }
        exit(EXIT_FAILURE);
    }
}

// Reference an NcVar and set the fill value based on the variable type. Down cast the double to be smaller
// if the type is not a double
void L1aFile::setVariableFillValue(NcVar &varRef, int &varType, double fillValue) {
    // set vill value and downcast the double to be something smaller if the type is smaller
    if (varType == NC_BYTE) {
        int8_t ncByte = (int8_t)fillValue;
        varRef.setFill(true, (void *)&ncByte);
    } else if (varType == NC_UBYTE) {
        uint8_t unsignedByte = (uint8_t)fillValue;
        varRef.setFill(true, (void *)&unsignedByte);
    } else if (varType == NC_SHORT) {
        int16_t ncShort = (int16_t)fillValue;
        varRef.setFill(true, (void *)&ncShort);
    } else if (varType == NC_USHORT) {
        uint16_t unsignedShort = (uint16_t)fillValue;
        varRef.setFill(true, (void *)&unsignedShort);
    } else if (varType == NC_INT) {
        int32_t ncInt = (int32_t)fillValue;
        varRef.setFill(true, (void *)&ncInt);
    } else if (varType == NC_UINT) {
        uint32_t unsignedInt = (uint32_t)fillValue;
        varRef.setFill(true, (void *)&unsignedInt);
    } else if (varType == NC_FLOAT) {
        float ncFloat = (float)fillValue;
        varRef.setFill(true, (void *)&ncFloat);
    }
    // otherwise, it is a double and dont need to cast
    else {
        varRef.setFill(true, (void *)&fillValue);
    }
}

// create a variable for a group but it contains an array of flag values and does not have
// min and max attributes
void L1aFile::createFlagVariable(NcGroup &group, string varName, int varType, string longName,
                                 vector<NcDim> &varDims, bool hasFillValue, double fillValue,
                                 vector<double> &flagValues, string flagMeaning) {
    try {
        // create the variable with type and dimensions set
        NcVar ncVariable = group.addVar(varName, varType, varDims);

        // set a fill value if given for a flag variable
        if (hasFillValue) {
            setVariableFillValue(ncVariable, varType, fillValue);
        }

        // add variable attributes to it

        ncVariable.putAtt("long_name", longName);

        // NetCDF API will convert the doubles to be same type as varType. So -999.0 will be -999 if
        // the varType == int
        ncVariable.putAtt("flag_values", varType, flagValues.size(), flagValues.data());

        ncVariable.putAtt("flag_meanings", flagMeaning);

    } catch (NcException &e) {
        cerr << "Error making flag variable: " << varName << " Error: " << e.what() << endl;
        exit(EXIT_FAILURE);
    }
}

// Create a variable for the referenced group
void L1aFile::createVariable(NcGroup &group, string varName, int varType, string longName,
                             vector<NcDim> &varDims, double fillValue, double validMin, double validMax,
                             string units, string reference) {
    try {
        // create the variable with type and dimensions set
        NcVar ncVariable = group.addVar(varName, varType, varDims);

        // if the var name contains "sci_" and has 3 dimensions, then set chunking for it

        // find "sci_" and it does not return a failure code npos. if it is npos, it failed to find the str
        bool containsSci = varName.find("sci_") != string::npos;
        bool has3Dims = varDims.size() == 3;
        if (containsSci && has3Dims) {
            // pass in the variable reference and the current dims to be chunked down
            applyChunkingAndCompressionToSciData(ncVariable, varDims);
        }

        // set a fill value if the smallest possible value for a variable is not the same as the fill
        // value. Otherwise, the validMin is the fill value and wont need any filling
        if (validMin != fillValue) {
            setVariableFillValue(ncVariable, varType, fillValue);
        }

        // add variable attributes to it

        ncVariable.putAtt("long_name", longName);

        // NetCDF API will convert the doubles to be same type as varType. So -999.0 will be -999 if
        // the varType == int
        ncVariable.putAtt("valid_min", varType, validMin);
        ncVariable.putAtt("valid_max", varType, validMax);

        // check if units is given. not all var has units defined
        if (!units.empty()) {
            ncVariable.putAtt("units", units);
        }

        if (!reference.empty()) {
            ncVariable.putAtt("reference", reference);
        }

    } catch (NcException &e) {
        cerr << "Error making variable: " << varName << " Error: " << e.what() << endl;
        exit(EXIT_FAILURE);
    }
}

int L1aFile::initializeL1aFile(char *l1aFileName, uint16_t maxScans, uint16_t numCcdPixels,
                               uint16_t numBlueBands, uint16_t numRedBands, uint16_t numSwirPixels,
                               uint16_t numDarkCcdPixels) {
    // attempt to make the NcFile and catch any issues if something is wrong
    try {
        this->l1afile = new NcFile(l1aFileName, NcFile::replace);
    } catch (NcException &e) {
        cerr << e.what() << "\nFailure initializing OCI L1A file: " << l1aFileName << endl;
        exit(1);
    }

    this->fileName = l1aFileName;

    // SET DIMENSIONS and have var reference so you can set variable parameters
    NcDim numScansDim = this->l1afile->addDim("number_of_scans", NC_UNLIMITED);
    NcDim numMceScansDim = this->l1afile->addDim("number_of_mce_scans", NC_UNLIMITED);
    NcDim numScaScansDim = this->l1afile->addDim("number_of_sca_scans", NC_UNLIMITED);
    NcDim numAttRecordsDim = this->l1afile->addDim("att_records", NC_UNLIMITED);
    NcDim orbRecordsDim = this->l1afile->addDim("orb_records", NC_UNLIMITED);
    NcDim tlmPacketsDim = this->l1afile->addDim("tlm_packets", NC_UNLIMITED);
    NcDim ccdPixelsDim = this->l1afile->addDim("ccd_pixels", numCcdPixels);
    NcDim swirPixelsDim = this->l1afile->addDim("SWIR_pixels", numSwirPixels);
    NcDim blueBandsDim = this->l1afile->addDim("blue_bands", numBlueBands);
    NcDim redBandsDim = this->l1afile->addDim("red_bands", numRedBands);
    NcDim swirBandsDim = this->l1afile->addDim("SWIR_bands", 9);
    NcDim dcPixelsDim = this->l1afile->addDim("DC_pixels", numDarkCcdPixels);
    NcDim numTapsDim = this->l1afile->addDim("number_of_taps", 16);
    NcDim spatialZonesDim = this->l1afile->addDim("spatial_zones", 10);
    NcDim quaternionElementsDim = this->l1afile->addDim("quaternion_elements", 4);
    NcDim vectorElementsDim = this->l1afile->addDim("vector_elements", 3);
    NcDim encoderSamplesDim = this->l1afile->addDim("encoder_samples", 200);
    NcDim encoderChannelsDim = this->l1afile->addDim("encoder_channels", 4);
    NcDim tiltSamplesDim = this->l1afile->addDim("tilt_samples", NC_UNLIMITED);
    NcDim mceBlocksDim = this->l1afile->addDim("MCE_block", 480);
    this->l1afile->addDim("sidecar_tlm", 76);
    NcDim ancilTmlDim = this->l1afile->addDim("ancil_tlm", 6);
    NcDim dauTlmDim = this->l1afile->addDim("DAU_tlm", 620);
    NcDim ddcTlmDim = this->l1afile->addDim("DDC_tlm", 524);
    NcDim tcTlmDim = this->l1afile->addDim("TC_tlm", 1216);
    NcDim icduMceTempTlmDim = this->l1afile->addDim("ICDU_MCE_temp_tlm", 76);
    NcDim adcLatDim = this->l1afile->addDim("ADC_lat", 4);
    NcDim icduThermDim = this->l1afile->addDim("ICDU_therm", 74);
    NcDim daucTempsDim = this->l1afile->addDim("DAUC_temps", 69);
    NcDim icduMceTempsDim = this->l1afile->addDim("ICDU_MCE_temps", 16);
    NcDim lineSkipsDim = this->l1afile->addDim("lin_skips", 33);

    // SET GLOBAL ATTRIBUTES
    // set the rest of the global attributes near the end of L1A processing for additiona user input
    this->l1afile->putAtt("title", "PACE OCI Level-1A Data");
    this->l1afile->putAtt("instrument", "OCI");
    this->l1afile->putAtt("platform", "PACE");
    this->l1afile->putAtt("product_name", l1aFileName);
    this->l1afile->putAtt("processing_version", "V1.0");  // default, will update in writeGlobalMetaData
    this->l1afile->putAtt("processing_level", "L1A");
    this->l1afile->putAtt("cdm_data_type", "swath");
    this->l1afile->putAtt("data_collect_mode", "Earth Collect");
    this->l1afile->putAtt("SWIR_data_mode", "Science");
    this->l1afile->putAtt("CDS_mode", "CDS");
    this->l1afile->putAtt("CDL_version_date", "2024-11-07"); // TODO REMOVE AFTER regression test has new test files

    // CREATE GROUPS AND THEIR VARIABLES
    // order of each variable made is following the original .cdl file order

    NcVar var;  // store var ref for multiple variables when the only thing being set is the long_name
    vector<NcDim> tempVarParameter = {};  // store temp variable parameters that are used only once
    int tempVarType =
        -1;  // temp nc type for when writing fill values to vars that only has fillvalue and long name

    // ###############################################
    // ######### scan_line_attributes group ##########
    // ###############################################

    NcGroup scanLineAttGroup = this->l1afile->addGroup("scan_line_attributes");
    vector<NcDim> dimsParameter = {numScansDim};  // dims used for entire group

    // scan_start_time
    createVariable(scanLineAttGroup,                            // group
                   string("scan_start_time"),                   // name
                   NC_DOUBLE,                                   // type
                   string("Scan start time (seconds of day)"),  // long name
                   dimsParameter,                               // dims
                   BAD_FLT,                                     // fill value
                   0.0,                                         // min
                   172800.0,                                    // max
                   string("seconds since "),                    // units
                   string("")                                   // reference
    );

    // scan_start_CCSDS_sec
    createVariable(scanLineAttGroup,                                      // group
                   string("scan_start_CCSDS_sec"),                        // name
                   NC_UINT,                                               // type
                   string("Scan start CCSDS time (seconds since 1958)"),  // long name
                   dimsParameter,                                         // dims
                   0.0,                                                   // fill value
                   1900000000.0,                                          // min
                   2400000000.0,                                          // max
                   string("seconds"),                                     // units
                   string("")                                             // reference
    );

    // scan_start_CCSDS_usec
    createVariable(scanLineAttGroup,                                // group
                   string("scan_start_CCSDS_usec"),                 // name
                   NC_INT,                                          // type
                   string("Scan start CCSDS time (microseconds)"),  // long name
                   dimsParameter,                                   // dims
                   BAD_FLT,                                         // fill value
                   0.0,                                             // min
                   999999.0,                                        // max
                   string("microseconds"),                          // units
                   string("")                                       // reference
    );

    // spin_ID
    createVariable(scanLineAttGroup,                                // group
                   string("spin_ID"),                               // name
                   NC_INT,                                          // type
                   string("Telescope spin counter from power-up"),  // long name
                   dimsParameter,                                   // dims
                   BAD_FLT,                                         // fill value
                   0.0,                                             // min
                   2147483647.0,                                    // max
                   string(""),                                      // units
                   string("")                                       // reference
    );

    // HAM_side
    createVariable(scanLineAttGroup,                  // group
                   string("HAM_side"),                // name
                   NC_UBYTE,                          // type
                   string("Half-angle mirror side"),  // long name
                   dimsParameter,                     // dims
                   255.0,                             // fill value
                   0.0,                               // min
                   1.0,                               // max
                   string(""),                        // units
                   string("")                         // reference
    );

    // pseq_flag
    vector<double> errorNoErrorFlagVals = {0, 1};
    createFlagVariable(scanLineAttGroup,                                     // group
                       string("pseq_flag"),                                  // name
                       NC_BYTE,                                              // type
                       string("Science packet sequence number error flag"),  // long name
                       dimsParameter,                                        // dims
                       false,                                                // has fillValue?
                       -1.0,                                                 // fillValue
                       errorNoErrorFlagVals,                                 // list of flag values
                       string("no_error error")                              // flag value meaning
    );

    // line_flag - has the same flag values as pseq_flag
    createFlagVariable(scanLineAttGroup,                      // group
                       string("line_flag"),                   // name
                       NC_BYTE,                               // type
                       string("CCD line number error flag"),  // long name
                       dimsParameter,                         // dims
                       false,                                 // has fillValue?
                       -1.0,                                  // fillValue
                       errorNoErrorFlagVals,                  // list of flag values
                       string("no_error error")               // flag value meaning
    );

    // #################################################
    // ######### spatial_spectral_modes group ##########
    // #################################################

    NcGroup spatialSpectralModesGroup = this->l1afile->addGroup("spatial_spectral_modes");

    // next few variables will have spatial_zones as their only parameter

    vector<NcDim> spatialZoneParameter = {spatialZonesDim};

    // spatial_zone_data_type
    vector<double> spatialZoneDataTypeFlags = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13};
    string spatialZoneDataTypeFlagMeanings =
        "no_data earth dark_cal solar_daily solar_monthly response_curve lunar diagnostic static "
        "earth_spectrum no_processing external_snapshop_trigger internal_snapshop_trigger "
        "non-baseline_spectral_aggregation";
    createFlagVariable(spatialSpectralModesGroup,                  // group
                       string("spatial_zone_data_type"),           // name
                       NC_SHORT,                                   // type
                       string("CCD spatial aggregation zone ID"),  // long name
                       spatialZoneParameter,                       // dims
                       true,                                       // has fillValue?
                       BAD_FLT,                                    // fillValue
                       spatialZoneDataTypeFlags,                   // list of flag values
                       spatialZoneDataTypeFlagMeanings             // flag value meaning
    );

    // spatial_aggregation
    vector<double> spatialAggFactorFlags = {0, 1, 2, 4, 8};
    createFlagVariable(spatialSpectralModesGroup,                   // group
                       string("spatial_aggregation"),               // name
                       NC_SHORT,                                    // type
                       string("CCD spatial aggregation per zone"),  // long name
                       spatialZoneParameter,                        // dims
                       true,                                        // has fillValue?
                       BAD_FLT,                                     // fillValue
                       spatialAggFactorFlags,                       // list of flag values
                       string("no_data 1to1 2to1 4to1 8to1")        // flag value meaning
    );

    // spatial_zone_lines
    createVariable(spatialSpectralModesGroup,                 // group
                   string("spatial_zone_lines"),              // name
                   NC_SHORT,                                  // type
                   string("CCD lines per aggregation zone"),  // long name
                   spatialZoneParameter,                      // dims
                   BAD_FLT,                                   // fill value
                   0.0,                                       // min
                   32384.0,                                   // max
                   string(""),                                // units
                   string("")                                 // reference
    );

    // parameter for blue and red spectral mode variables
    vector<NcDim> numOfTapsParameter = {numTapsDim};

    // blue_spectral_mode
    createFlagVariable(spatialSpectralModesGroup,                     // group
                       string("blue_spectral_mode"),                  // name
                       NC_SHORT,                                      // type
                       string("Blue CCD spectral aggregation mode"),  // long name
                       numOfTapsParameter,                            // dims
                       true,                                          // has fillValue?
                       BAD_FLT,                                       // fillValue
                       spatialAggFactorFlags,                  // flag val same as spatial aggrigation factors
                       string("disabled 1to1 2to1 4to1 8to1")  // flag value meaning
    );

    // red_spectral_mode
    createFlagVariable(spatialSpectralModesGroup,                    // group
                       string("red_spectral_mode"),                  // name
                       NC_SHORT,                                     // type
                       string("Red CCD spectral aggregation mode"),  // long name
                       numOfTapsParameter,                           // dims
                       true,                                         // has fillValue?
                       BAD_FLT,                                      // fillValue
                       spatialAggFactorFlags,                  // flag val same as spatial aggrigation factors
                       string("disabled 1to1 2to1 4to1 8to1")  // flag value meaning
    );

    // paramater for line skips
    vector<NcDim> lineSkipsParameter = {lineSkipsDim};

    // aux_param_table
    createVariable(
        spatialSpectralModesGroup,                                                      // group
        string("aux_param_table"),                                                      // name
        NC_SHORT,                                                                       // type
        string("Auxilary parameter table containing linearity mode line skip values"),  // long name
        lineSkipsParameter,                                                             // dims
        BAD_FLT,                                                                        // fill value
        0.0,                                                                            // min
        99.0,                                                                           // max
        string(""),                                                                     // units
        string("")                                                                      // reference
    );

    // #################################################
    // ######### engineering_data group ################
    // #################################################

    NcGroup engineeringDataGroup = this->l1afile->addGroup("engineering_data");

    // DAU_telemetry -- only has long name so dont use the createVariable() call
    tempVarParameter = {tlmPacketsDim, dauTlmDim};
    var = engineeringDataGroup.addVar("DAU_telemetry", NC_UBYTE, tempVarParameter);
    var.putAtt("long_name", "OCI DAU raw telemetry");

    // tlm packets only parameter
    vector<NcDim> tlmPacketsParameter = {tlmPacketsDim};

    // DAU_tlm_time
    createVariable(engineeringDataGroup,                                  // group
                   string("DAU_tlm_time"),                                // name
                   NC_DOUBLE,                                             // type
                   string("DAU telemetry packet time (seconds of day)"),  // long name
                   tlmPacketsParameter,                                   // dims
                   BAD_FLT,                                               // fill value
                   0.0,                                                   // min
                   172800.0,                                              // max
                   string("seconds"),                                     // units
                   string("")                                             // reference
    );

    // DAU_spin_ID
    createVariable(engineeringDataGroup,                                    // group
                   string("DAU_spin_ID"),                                   // name
                   NC_INT,                                                  // type
                   string("Telescope spin ID from DAU telemetry packets"),  // long name
                   tlmPacketsParameter,                                     // dims
                   BAD_FLT,                                                 // fill value
                   0.0,                                                     // min
                   2147483647.0,                                            // max
                   string(""),                                              // units
                   string("")                                               // reference
    );

    // DDC_telemetry -- only has long name so dont use the createVariable() call
    tempVarParameter = {tlmPacketsDim, ddcTlmDim};  // tlm packets and ddc tlm parameter
    var = engineeringDataGroup.addVar("DDC_telemetry", NC_UBYTE, tempVarParameter);
    var.putAtt("long_name", "OCI DDC raw telemetry");

    // DDC_tlm_time
    createVariable(engineeringDataGroup,                                  // group
                   string("DDC_tlm_time"),                                // name
                   NC_DOUBLE,                                             // type
                   string("DDC telemetry packet time (seconds of day)"),  // long name
                   tlmPacketsParameter,                                   // dims
                   BAD_FLT,                                               // fill value
                   0.0,                                                   // min
                   172800.0,                                              // max
                   string("seconds"),                                     // units
                   string("")                                             // reference
    );

    // TC_telemetry -- only has long name so dont use the createVariable() call
    tempVarParameter = {tlmPacketsDim, tcTlmDim};  // tlm packets and tc tlm parameter
    var = engineeringDataGroup.addVar("TC_telemetry", NC_UBYTE, tempVarParameter);
    var.putAtt("long_name", "OCI temperature controller raw telemetry");

    // TC_tlm_time
    createVariable(engineeringDataGroup,                                  // group
                   string("TC_tlm_time"),                                 // name
                   NC_DOUBLE,                                             // type
                   string("TTC telemetry packet time (seconds of day)"),  // long name
                   tlmPacketsParameter,                                   // dims
                   BAD_FLT,                                               // fill value
                   0.0,                                                   // min
                   172800.0,                                              // max
                   string("seconds"),                                     // units
                   string("")                                             // reference
    );

    // DAUC_temp_time
    createVariable(engineeringDataGroup,                                         // group
                   string("DAUC_temp_time"),                                     // name
                   NC_DOUBLE,                                                    // type
                   string("FSW DAUC temperature packet time (seconds of day)"),  // long name
                   tlmPacketsParameter,                                          // dims
                   BAD_FLT,                                                      // fill value
                   0.0,                                                          // min
                   172800.0,                                                     // max
                   string("seconds"),                                            // units
                   string("")                                                    // reference
    );

    // ICDU_MCE_temp_tlm -- only has long name so dont use the createVariable() call
    tempVarParameter = {tlmPacketsDim, icduMceTempTlmDim};
    var = engineeringDataGroup.addVar("ICDU_MCE_temp_tlm", NC_UBYTE, tempVarParameter);
    var.putAtt("long_name", "ICDU-MCE raw temperature telemetry");

    // ICDU_MCE_temp_time
    createVariable(engineeringDataGroup,                                             // group
                   string("ICDU_MCE_temp_time"),                                     // name
                   NC_DOUBLE,                                                        // type
                   string("FSW ICDU-MCE temperature packet time (seconds of day)"),  // long name
                   tlmPacketsParameter,                                              // dims
                   BAD_FLT,                                                          // fill value
                   0.0,                                                              // min
                   172800.0,                                                         // max
                   string("seconds"),                                                // units
                   string("")                                                        // reference
    );

    // ancillary_tlm -- only has long name and valid min so dont use the createVariable() call
    tempVarParameter = {numScansDim, ancilTmlDim};
    tempVarType = NC_SHORT;
    var = engineeringDataGroup.addVar("ancillary_tlm", NC_SHORT, tempVarParameter);
    setVariableFillValue(var, tempVarType, BAD_FLT);
    var.putAtt("long_name", "Ancillary telemetry status flags and error counts");

    // MCE_telemetry -- only has long name and valid min so dont use the createVariable() call
    tempVarParameter = {numMceScansDim, mceBlocksDim};
    var = engineeringDataGroup.addVar("MCE_telemetry", NC_UBYTE, tempVarParameter);
    var.putAtt("long_name", "RTA/HAM mechanism control electronics telemetry");

    // number of mce_scans variable parameter
    vector<NcDim> numMceScansParameter = {numMceScansDim};

    // MCE_spin_ID
    createVariable(engineeringDataGroup,                                  // group
                   string("MCE_spin_ID"),                                 // name
                   NC_INT,                                                // type
                   string("Spin ID from RTA/HAM MCE telemetry packets"),  // long name
                   numMceScansParameter,                                  // dims
                   BAD_FLT,                                               // fill value
                   0.0,                                                   // min
                   2147483647.0,                                          // max
                   string(""),                                            // units
                   string("")                                             // reference
    );

    // MCE_encoder_data
    tempVarParameter = {numMceScansDim, encoderSamplesDim, encoderChannelsDim};
    var = engineeringDataGroup.addVar("MCE_encoder_data", NC_SHORT, tempVarParameter);
    tempVarType = NC_SHORT;
    setVariableFillValue(var, tempVarType, BAD_FLT);
    var.putAtt("long_name", "RTA/HAM encoder samples (channel order is 3, 2, 1, 0)");

    // encoder_spin_ID
    createVariable(engineeringDataGroup,                                // group
                   string("encoder_spin_ID"),                           // name
                   NC_INT,                                              // type
                   string("Spin ID from RTA/HAM MCE encoder packets"),  // long name
                   numMceScansParameter,                                // dims
                   BAD_FLT,                                             // fill value
                   0.0,                                                 // min
                   2147483647.0,                                        // max
                   string(""),                                          // units
                   string("")                                           // reference
    );

    // SCA_telemetry
    tempVarParameter = {numScaScansDim, mceBlocksDim};
    var = engineeringDataGroup.addVar("SCA_telemetry", NC_UBYTE, tempVarParameter);
    var.putAtt("long_name", "SCA mechanism control electronics telemetry");

    vector<NcDim> numScaScansParameter = {numScaScansDim};

    // SCA_tlm_time
    createVariable(engineeringDataGroup,                                  // group
                   string("SCA_tlm_time"),                                // name
                   NC_DOUBLE,                                             // type
                   string("SCA telemetry packet time (seconds of day)"),  // long name
                   numScaScansParameter,                                  // dims
                   BAD_FLT,                                               // fill value
                   0.0,                                                   // min
                   172800.0,                                              // max
                   string("seconds"),                                     // units
                   string("")                                             // reference
    );

    // SCA_spin_ID
    createVariable(engineeringDataGroup,                              // group
                   string("SCA_spin_ID"),                             // name
                   NC_INT,                                            // type
                   string("Spin ID from SCA MCE telemetry packets"),  // long name
                   numScaScansParameter,                              // dims
                   BAD_FLT,                                           // fill value
                   0.0,                                               // min
                   2147483647.0,                                      // max
                   string(""),                                        // units
                   string("")                                         // reference
    );

    // SCA_diffuser_position
    createVariable(engineeringDataGroup,                   // group
                   string("SCA_diffuser_position"),        // name
                   NC_FLOAT,                               // type
                   string("SCA diffuser position angle"),  // long name
                   numScaScansParameter,                   // dims
                   BAD_FLT,                                // fill value
                   0.0,                                    // min
                   360.0,                                  // max
                   string("degrees"),                      // units
                   string("")                              // reference
    );

    // SCA_encoder_data
    tempVarParameter = {numScaScansDim, encoderSamplesDim, encoderChannelsDim};
    var = engineeringDataGroup.addVar("SCA_encoder_data", NC_SHORT, tempVarParameter);
    tempVarType = NC_SHORT;
    setVariableFillValue(var, tempVarType, BAD_FLT);
    var.putAtt("long_name", "SCA encoder samples (channel order is 3, 2, 1, 0)");

    // SCA_encoder_spin_ID
    createVariable(engineeringDataGroup,                        // group
                   string("SCA_encoder_spin_ID"),               // name
                   NC_INT,                                      // type
                   string("Spin ID from SCA encoder packets"),  // long name
                   numScaScansParameter,                        // dims
                   BAD_FLT,                                     // fill value
                   0.0,                                         // min
                   2147483647.0,                                // max
                   string(""),                                  // units
                   string("")                                   // reference
    );

    // var parameter containing only number_of_scans
    vector<NcDim> numScansParameter = {numScansDim};

    // agg_control
    createVariable(engineeringDataGroup,                  // group
                   string("agg_control"),                 // name
                   NC_INT,                                // type
                   string("Aggregation control fields"),  // long name
                   numScansParameter,                     // dims
                   BAD_FLT,                               // fill value
                   0.0,                                   // min
                   2047.0,                                // max
                   string(""),                            // units
                   string("")                             // reference
    );

    // blue_agg_error
    var = engineeringDataGroup.addVar("blue_agg_error", NC_USHORT, numScansParameter);
    var.putAtt("long_name", "UVVIS aggregation error");

    // red_agg_error
    var = engineeringDataGroup.addVar("red_agg_error", NC_USHORT, numScansParameter);
    var.putAtt("long_name", "VISNIR aggregation error");

    // dig_card_error
    createVariable(engineeringDataGroup,                 // group
                   string("dig_card_error"),             // name
                   NC_INT,                               // type
                   string("Digital card error status"),  // long name
                   numScansParameter,                    // dims
                   BAD_FLT,                              // fill value
                   0.0,                                  // min
                   1048575.0,                            // max
                   string(""),                           // units
                   string("")                            // reference
    );

    // CDS_disable
    createFlagVariable(engineeringDataGroup,                               // group
                       string("CDS_disable"),                              // name
                       NC_UBYTE,                                           // type
                       string("Correlated double sampling disable flag"),  // long name
                       tlmPacketsParameter,                                // dims
                       true,                                               // has fillValue?
                       255.0,                                              // fillValue
                       errorNoErrorFlagVals,                               // flag val containly only 0 and 1
                       string("Enable Disable")                            // flag value meaning
    );

    // ADC_latency
    tempVarParameter = {tlmPacketsDim, adcLatDim};
    vector<double> adcLatencyFlagValues = {15, 16};
    createFlagVariable(engineeringDataGroup,                                      // group
                       string("ADC_latency"),                                     // name
                       NC_UBYTE,                                                  // type
                       string("ADC latency (red_hi, red_lo, blue_hi, blue_lo)"),  // long name
                       tempVarParameter,                                          // dims
                       true,                                                      // has fillValue?
                       255.0,                                                     // fillValue
                       adcLatencyFlagValues,                                      // flag val of 15 and 16
                       string("Reset Video")                                      // flag value meaning
    );

    // ICDU_thermisters
    tempVarParameter = {tlmPacketsDim, icduThermDim};
    createVariable(engineeringDataGroup,                               // group
                   string("ICDU_thermisters"),                         // name
                   NC_FLOAT,                                           // type
                   string("ICDU thermistor data from FSW TC packet"),  // long name
                   tempVarParameter,                                   // dims
                   BAD_FLT,                                            // fill value
                   -200.0,                                             // min
                   100.0,                                              // max
                   string(""),                                         // units
                   string("OCI-THRM-SPEC-0108 tables 6-1 and 6-2")     // reference
    );

    // DAUC_temperatures
    tempVarParameter = {tlmPacketsDim, daucTempsDim};
    createVariable(engineeringDataGroup,                                  // group
                   string("DAUC_temperatures"),                           // name
                   NC_FLOAT,                                              // type
                   string("DAUC temperatures from FSW DAUCTEMP packet"),  // long name
                   tempVarParameter,                                      // dims
                   BAD_FLT,                                               // fill value
                   -200.0,                                                // min
                   200.0,                                                 // max
                   string(""),                                            // units
                   string("OCI-SYS-DESC-0195")                            // reference
    );

    // ICDU_MCE_temperatures
    tempVarParameter = {tlmPacketsDim, icduMceTempsDim};
    createVariable(engineeringDataGroup,                                     // group
                   string("ICDU_MCE_temperatures"),                          // name
                   NC_FLOAT,                                                 // type
                   string("ICDU-MCE temperatures from FSW ICDUMCE packet"),  // long name
                   tempVarParameter,                                         // dims
                   BAD_FLT,                                                  // fill value
                   -200.0,                                                   // min
                   200.0,                                                    // max
                   string(""),                                               // units
                   string("OCI-SYS-DESC-0195")                               // reference
    );

    // blue_channel_mask
    var = engineeringDataGroup.addVar("blue_channel_mask", NC_UINT, numOfTapsParameter);
    tempVarType = NC_UINT;
    setVariableFillValue(var, tempVarType, 4294967295.0);
    var.putAtt("long_name", "Channel mask for blue CCD taps");
    var.putAtt("reference", "OCI-ELEC-SPEC-0009");

    // red_channel_mask
    var = engineeringDataGroup.addVar("red_channel_mask", NC_UINT, numOfTapsParameter);
    tempVarType = NC_UINT;
    setVariableFillValue(var, tempVarType, 4294967295.0);
    var.putAtt("long_name", "Channel mask for red CCD taps");
    var.putAtt("reference", "OCI-ELEC-SPEC-0009");

    // #################################################
    // ######### navigation_data group #################
    // #################################################

    NcGroup navigationDataGroup = this->l1afile->addGroup("navigation_data");

    // att_time
    tempVarParameter = {numAttRecordsDim};                           // parameter only used once
    createVariable(navigationDataGroup,                              // group
                   string("att_time"),                               // name
                   NC_DOUBLE,                                        // type
                   string("Attitude sample time (seconds of day)"),  // long name
                   tempVarParameter,                                 // dims
                   BAD_FLT,                                          // fill value
                   0.0,                                              // min
                   172800.0,                                         // max
                   string("seconds"),                                // units
                   string("")                                        // reference
    );

    // att_quat
    tempVarParameter = {numAttRecordsDim, quaternionElementsDim};
    createVariable(navigationDataGroup,                                   // group
                   string("att_quat"),                                    // name
                   NC_FLOAT,                                              // type
                   string("Attitude quaternions (J2000 to spacecraft)"),  // long name
                   tempVarParameter,                                      // dims
                   BAD_FLT,                                               // fill value
                   -1.0,                                                  // min
                   1.0,                                                   // max
                   string(""),                                            // units
                   string("")                                             // reference
    );

    // att_rate
    tempVarParameter = {numAttRecordsDim, vectorElementsDim};
    createVariable(navigationDataGroup,                                   // group
                   string("att_rate"),                                    // name
                   NC_FLOAT,                                              // type
                   string("Attitude angular rates in spacecraft frame"),  // long name
                   tempVarParameter,                                      // dims
                   BAD_FLT,                                               // fill value
                   -0.004,                                                // min
                   0.004,                                                 // max
                   string("radians/second"),                              // units
                   string("")                                             // reference
    );

    // orb_time
    tempVarParameter = {orbRecordsDim};
    createVariable(navigationDataGroup,                           // group
                   string("orb_time"),                            // name
                   NC_DOUBLE,                                     // type
                   string("Orbit vector time (seconds of day)"),  // long name
                   tempVarParameter,                              // dims
                   BAD_FLT,                                       // fill value
                   0,                                             // min
                   172800.0,                                      // max
                   string("seconds"),                             // units
                   string("")                                     // reference
    );

    // orb_pos
    vector<NcDim> orbitPosVelParameters = {orbRecordsDim, vectorElementsDim};
    createVariable(navigationDataGroup,                     // group
                   string("orb_pos"),                       // name
                   NC_FLOAT,                                // type
                   string("Orbit position vectors (ECR)"),  // long name
                   orbitPosVelParameters,                   // dims
                   -9999999.0,                              // fill value
                   -7200000.0,                              // min
                   7200000.0,                               // max
                   string("meters"),                        // units
                   string("")                               // reference
    );

    // orb_vel
    createVariable(navigationDataGroup,                     // group
                   string("orb_vel"),                       // name
                   NC_FLOAT,                                // type
                   string("Orbit velocity vectors (ECR)"),  // long name
                   orbitPosVelParameters,                   // dims
                   BAD_FLT,                                 // fill value
                   -7600.0,                                 // min
                   7600.0,                                  // max
                   string("meters/second"),                 // units
                   string("")                               // reference
    );

    // tilt_time
    vector<NcDim> tiltParameters = {tiltSamplesDim};
    createVariable(navigationDataGroup,                   // group
                   string("tilt_time"),                   // name
                   NC_DOUBLE,                             // type
                   string("Tilt time (seconds of day)"),  // long name
                   tiltParameters,                        // dims
                   BAD_FLT,                               // fill value
                   0,                                     // min
                   172800.0,                              // max
                   string("seconds"),                     // units
                   string("")                             // reference
    );

    // tilt - valid min and max has percision so not using createVariable
    var = navigationDataGroup.addVar("tilt", NC_FLOAT, tiltParameters);
    tempVarType = NC_FLOAT;
    setVariableFillValue(var, tempVarType, BAD_FLT);
    var.putAtt("long_name", "Tilt angle");
    var.putAtt("valid_min", NC_FLOAT, -20.5);
    var.putAtt("valid_max", NC_FLOAT, 20.5);
    var.putAtt("units", "degrees");

    // ##########################################################
    // ######### onboard_calibration_data group #################
    // ##########################################################

    NcGroup onboardCalibrationDataGroup = this->l1afile->addGroup("onboard_calibration_data");

    // DC_blue
    tempVarParameter = {numScansDim, blueBandsDim, dcPixelsDim};
    createVariable(
        onboardCalibrationDataGroup,                                                            // group
        string("DC_blue"),                                                                      // name
        NC_USHORT,                                                                              // type
        string("Dark calibration data for blue focal plane (wavelengths in ascending order)"),  // long name
        tempVarParameter,                                                                       // dims
        BAD_UINT,                                                                               // fill value
        0,                                                                                      // min
        65530.0,                                                                                // max
        string("counts"),                                                                       // units
        string("")                                                                              // reference
    );

    // DC_red
    tempVarParameter = {numScansDim, redBandsDim, dcPixelsDim};
    createVariable(
        onboardCalibrationDataGroup,                                                           // group
        string("DC_red"),                                                                      // name
        NC_USHORT,                                                                             // type
        string("Dark calibration data for red focal plane (wavelengths in ascending order)"),  // long name
        tempVarParameter,                                                                      // dims
        BAD_UINT,                                                                              // fill value
        0,                                                                                     // min
        65530.0,                                                                               // max
        string("counts"),                                                                      // units
        string("")                                                                             // reference
    );

    // DC_SWIR
    tempVarParameter = {numScansDim, swirBandsDim, dcPixelsDim};
    createVariable(onboardCalibrationDataGroup,                     // group
                   string("DC_SWIR"),                               // name
                   NC_UINT,                                         // type
                   string("Dark calibration data for SWIR bands"),  // long name
                   tempVarParameter,                                // dims
                   1048575.0,                                       // fill value
                   0,                                               // min
                   1048570.0,                                       // max
                   string("counts"),                                // units
                   string("")                                       // reference
    );

    // frm_type_DC_SWIR
    tempVarParameter = {numScansDim, dcPixelsDim};
    vector<double> frameTypeFlags = {0, 1, 2, 3, 4};
    createFlagVariable(
        onboardCalibrationDataGroup,                                         // group
        string("frm_type_DC_SWIR"),                                          // name
        NC_BYTE,                                                             // type
        string("SWIR DC frame type for non-science modes"),                  // long name
        tempVarParameter,                                                    // dims
        true,                                                                // has fillValue?
        -1.0,                                                                // fillValue
        frameTypeFlags,                                                      // flag val 0-4
        string("chan_0-7 chan_7-15 chan_16-23 chan_24-31 test_pattern_TDI")  // flag value meaning
    );

    // ##########################################################
    // ######### science_data group #############################
    // ##########################################################

    NcGroup sciDataGroup = this->l1afile->addGroup("science_data");

    // sci_blue
    tempVarParameter = {numScansDim, blueBandsDim, ccdPixelsDim};
    createVariable(sciDataGroup,                                                                  // group
                   string("sci_blue"),                                                            // name
                   NC_USHORT,                                                                     // type
                   string("Science data for blue focal plane (wavelengths in ascending order)"),  // long name
                   tempVarParameter,                                                              // dims
                   BAD_UINT,          // fill value
                   0,                 // min
                   65530.0,           // max
                   string("counts"),  // units
                   string("")         // reference
    );

    // sci_red
    tempVarParameter = {numScansDim, redBandsDim, ccdPixelsDim};
    createVariable(sciDataGroup,                                                                 // group
                   string("sci_red"),                                                            // name
                   NC_USHORT,                                                                    // type
                   string("Science data for red focal plane (wavelengths in ascending order)"),  // long name
                   tempVarParameter,                                                             // dims
                   BAD_UINT,                                                                     // fill value
                   0,                                                                            // min
                   65530.0,                                                                      // max
                   string("counts"),                                                             // units
                   string("")                                                                    // reference
    );

    // sci_SWIR
    tempVarParameter = {numScansDim, swirBandsDim, swirPixelsDim};
    createVariable(sciDataGroup,                           // group
                   string("sci_SWIR"),                     // name
                   NC_UINT,                                // type
                   string("Science data for SWIR bands"),  // long name
                   tempVarParameter,                       // dims
                   1048575.0,                              // fill value
                   0,                                      // min
                   1048570.0,                              // max
                   string("counts"),                       // units
                   string("")                              // reference
    );

    // frm_type_SWIR
    tempVarParameter = {numScansDim, swirPixelsDim};
    createFlagVariable(
        sciDataGroup,                                                        // group
        string("frm_type_SWIR"),                                             // name
        NC_BYTE,                                                             // type
        string("SWIR frame type for non-science modes"),                     // long name
        tempVarParameter,                                                    // dims
        true,                                                                // has fillValue?
        -1.0,                                                                // fillValue
        frameTypeFlags,                                                      // flag val 0-4
        string("chan_0-7 chan_7-15 chan_16-23 chan_24-31 test_pattern_TDI")  // flag value meaning
    );

    return 0;
}

int L1aFile::writeScienceData(DimensionShape* dimShape, uint32_t scanNum, uint16_t numBlueBands, uint16_t numRedBands,
                              uint16_t numSwirBands, uint16_t numCcdPixels, uint16_t numSwirPixels,
                              uint16_t **blueScienceData, uint16_t **redScienceData,
                              uint32_t **swirScienceData, int8_t *swirFrameTypesSciData) {
    // Writes one scan at a time

    NcVar blu_bands;
    NcVar red_bands;
    NcVar swir_bands;
    NcVar swir_frms;
    vector<size_t> start;
    vector<size_t> countBlue;
    vector<size_t> countRed;
    vector<size_t> countSwir;
    vector<size_t> countSingleEntry;  // for 1-D fill value variables

    vector<size_t> startFrms;
    vector<size_t> countFrms;

    NcGroup gid = l1afile->getGroup("science_data");
    blu_bands = gid.getVar("sci_blue");
    red_bands = gid.getVar("sci_red");
    swir_bands = gid.getVar("sci_SWIR");
    swir_frms = gid.getVar("frm_type_SWIR");

    // writes line by line. Get the current line to write for science data
    start = {dimShape->scienceScanNum, 0, 0};
    countBlue = {1, numBlueBands, numCcdPixels};
    countRed = {1, numRedBands, numCcdPixels};
    countSwir = {1, numSwirBands, numSwirPixels};
    countSingleEntry = {1, 1, 1};

    // sci_blue/sci_red/sci_SWIR with fill value if number of blue bands is 0
    (countBlue[1] > 0) ? blu_bands.putVar(start, countBlue, &blueScienceData[0][0])
                     : blu_bands.putVar(start, countSingleEntry, &blueScienceData[0][0]);
    (countRed[1] > 0) ? red_bands.putVar(start, countRed, &redScienceData[0][0])
                     : red_bands.putVar(start, countSingleEntry, &redScienceData[0][0]);
    (countSwir[1] > 0) ? swir_bands.putVar(start, countSwir, &swirScienceData[0][0])
                     : swir_bands.putVar(start, countSingleEntry, &swirScienceData[0][0]);


    // frm_type_SWIR
    startFrms = {dimShape->scienceScanNum, 0};
    countFrms = {1, numSwirPixels};

    if (countFrms[1] > 0)
        swir_frms.putVar(startFrms, countFrms, swirFrameTypesSciData);

    // increment for next science write
    dimShape->incrementScienceScanNumShape();

    return 0;
}

int L1aFile::writeProcessingControl(std::string hktList, std::string l0List, std::string granuleStartTime,
                                    std::string maxgap, std::string nametag, std::string swir_loff_set,
                                    std::string outlist, std::string outfile, std::string doi,
                                    std::string pversion, std::string isSPW, std::string VERSION) {
    try {
        // add the processing control group
        netCDF::NcGroup processCtrlGrp = l1afile->addGroup("processing_control");

        // within process control group, add input parameter group as its child
        netCDF::NcGroup inputParaGrp = processCtrlGrp.addGroup("input_parameter");

        // add current input parameter file names to the input parameter group.
        // this only contains the important command line args for l1agen_oci
        inputParaGrp.putAtt("OCI_packet_file", l0List);
        inputParaGrp.putAtt("maxgap", maxgap);
        inputParaGrp.putAtt("hktlist_iFile", hktList);
        inputParaGrp.putAtt("swir_loff_set", swir_loff_set);
        inputParaGrp.putAtt("start_time", granuleStartTime);
        inputParaGrp.putAtt("outlist", outlist);
        inputParaGrp.putAtt("outfile", outfile);
        inputParaGrp.putAtt("nametag", nametag);
        inputParaGrp.putAtt("doi", doi);
        inputParaGrp.putAtt("pversion", pversion);
        inputParaGrp.putAtt("isSPW", isSPW);

        // add current software details
        processCtrlGrp.putAtt("software_name", "l1agen_oci");
        processCtrlGrp.putAtt("software_version", VERSION);

        // given hkt list path, convert it into a vector and add
        // all the hkt files used separated by commas (,)
        vector<std::string> hktFiles = readFileList(hktList);
        string hktString = "";
        for (size_t i = 0; i < hktFiles.size(); i++) {
            //
            if (i < hktFiles.size() - 1) {
                hktString.append(hktFiles[i]);
                hktString.append(", ");
                continue;
            }
            // dont end with a  comma if last file
            hktString.append(hktFiles[i]);
        }
        processCtrlGrp.putAtt("hkt_list", hktString);

        // do the same for l0 file list
        vector<std::string> l0Files = readFileList(l0List);
        string l0String = "";
        for (size_t i = 0; i < l0Files.size(); i++) {
            if (i < l0Files.size() - 1) {
                l0String.append(l0Files[i]);
                l0String.append(", ");
                continue;
            }
            l0String.append(l0Files[i]);
        }
        processCtrlGrp.putAtt("l0_list", l0String);

    } catch (std::exception &e) {
        std::cerr << "--Error-- writing processing control group. what: " << e.what() << endl;
        return 0;
    }

    return 1;
}

// Given a pointer to the calibration data, write them line by line for red, blue and swir bands.
// If the number of bands is 0, then write in 1 line. the data for said line should be fill-values.
int L1aFile::writeCalibrationData(DimensionShape *dimShape, uint32_t isc, uint16_t numBlueBands, uint16_t numRedBands,
                                  uint16_t numSwirBands, uint16_t numDarkCcdPixels,
                                  uint16_t numDarkSwirPixels, uint16_t *blueDarkCalibrationData,
                                  uint16_t *redDarkCalibrationData, uint32_t *swirDarkCalibrationData,
                                  int8_t *swirDarkCalFrameTypeData) {
    
    // NetCDF group for all calibration data. l1afile is a static variable defined at the top of the file
    NcGroup calibrationGroup = l1afile->getGroup("onboard_calibration_data");

    // start netcdf indices for red, blue and swir dark calibration data
    vector<size_t> darkCalDataStart = {dimShape->numScans, 0, 0};

    // If number of bands is 0, then fill 1 line of fill values
    vector<size_t> darkCalDataFillValCount = {1, 1, 1};

    // --- WRITE BLUE CALIBRATION DATA ---
    NcVar blueDarkCalVar = calibrationGroup.getVar("DC_blue");
    vector<size_t> blueDarkCalCount = {isc, numBlueBands, numDarkCcdPixels};

    numBlueBands > 0
        // there's bands, add the data as is
        ? blueDarkCalVar.putVar(darkCalDataStart, blueDarkCalCount, blueDarkCalibrationData)
        // no bands, write 1 line of fill values for this (1D)
        : blueDarkCalVar.putVar(darkCalDataStart, darkCalDataFillValCount, blueDarkCalibrationData);

    // --- WRITE RED CALIBRATION DATA ---
    NcVar redDarkCalVar = calibrationGroup.getVar("DC_red");
    vector<size_t> redDarkCalCount = {isc, numRedBands, numDarkCcdPixels};

    numRedBands > 0 ? redDarkCalVar.putVar(darkCalDataStart, redDarkCalCount, redDarkCalibrationData)
                    : redDarkCalVar.putVar(darkCalDataStart, darkCalDataFillValCount, redDarkCalibrationData);

    // --- WRITE SWIR CALIBRATION DATA ---
    NcVar swirDarkCalVar = calibrationGroup.getVar("DC_SWIR");
    vector<size_t> swirDarkCalCount = {isc, numSwirBands, numDarkSwirPixels};

    numSwirBands > 0
        ? swirDarkCalVar.putVar(darkCalDataStart, swirDarkCalCount, swirDarkCalibrationData)
        : swirDarkCalVar.putVar(darkCalDataStart, darkCalDataFillValCount, swirDarkCalibrationData);

    // --- WRITE SWIR DARK CALIBRATION FRAME DATA ---
    // Write only if number of dark swir pixels is not 0
    if (numDarkSwirPixels > 0) {
        NcVar swirDarkCalFrameTypeVar = calibrationGroup.getVar("frm_type_DC_SWIR");

        // this one is 2D dataset so only requires 2 indices
        vector<size_t> swirDarkCalFrameTypeStart = {dimShape->numScans, 0};
        vector<size_t> swirDarkCalFrameTypeCount = {isc, numDarkSwirPixels};

        swirDarkCalFrameTypeVar.putVar(swirDarkCalFrameTypeStart, swirDarkCalFrameTypeCount,
                                       swirDarkCalFrameTypeData);
    }

    return 0;
}

// write out meta data for scan line attributes. Extract time from the ancillary data packet and report
// any errors
int L1aFile::writeScanMetaData(DimensionShape *dimShape, uint32_t scanNum, uint8_t *ancillaryData, uint8_t *sciPacketSequenceError,
                               int8_t *ccdLineError, int32_t *spinID, AncillaryPktTimeStamp &timeStruct) {
    
    int32_t year = 0;
    int32_t day = 0;
    double startTime = 0.0;

    // Extract and convert times to seconds of day
    // Extract CCSDS scan start times

    vector<double> startTimeArr(scanNum);
    short int timeOffset = 24;

    vector<uint32_t> ccsdsStartScanSeconds(scanNum);
    vector<int32_t> ccsdsStartScanSubSeconds(scanNum);

    double firstGoodTime = BAD_FLT;
    for (size_t i = 0; i < scanNum; i++) {
        uint32_t apid = (ancillaryData[i * ANCSIZE] % 8) * 256 + ancillaryData[i * ANCSIZE + 1];
        if (apid == ANCILLARY_APID) {
            getAncillaryPacketTime(&ancillaryData[i * ANCSIZE], year, day, startTime);
            if (firstGoodTime == BAD_FLT)
                firstGoodTime = startTime;
            if (startTime < firstGoodTime)
                startTimeArr[i] = startTime + SECONDS_IN_DAY;  // Ensure that startTimeArr[i] stays relative
                                                               // to the same day that the granule started
            else
                startTimeArr[i] = startTime;

            uint32_t ui32;
            memcpy(&ui32, &ancillaryData[i * ANCSIZE + timeOffset + 4], 4);
            swapc_bytes((char *)&ui32, sizeof(uint32_t), 1);
            ccsdsStartScanSeconds[i] = ui32;
            memcpy(&ui32, &ancillaryData[i * ANCSIZE + timeOffset + 8], 4);
            swapc_bytes((char *)&ui32, sizeof(uint32_t), 1);
            ccsdsStartScanSubSeconds[i] = ui32 / 4096;
        } else {
            startTimeArr[i] = BAD_FLT;
            ccsdsStartScanSeconds[i] = 0;
            ccsdsStartScanSubSeconds[i] = BAD_INT;
        }
    }

    // group id 
    NcGroup gid = l1afile->getGroup("scan_line_attributes");
    NcVar var;

    // variables for the NetCDF API on where to start and how much to put 
    vector<size_t> start = {dimShape->numScans};
    vector<size_t> count = {scanNum}; // total scans is how much we are inserting

    // Scan start time (seconds of day)
    var = gid.getVar("scan_start_time");
    
    // get the size of the current scan
    var.putVar(start, count, startTimeArr.data());
    
    // TODO: add an if statement so this is only done once. Crossing dateline might 
    // advance the date
    double tmpTime = yds2unix(timeStruct.year, timeStruct.day, timeStruct.second);
    string timeStr = unix2isodate(tmpTime, 'G');
    timeStr = timeStr.substr(0, 10);
    timeStr.insert(0, "seconds since ");
    var.putAtt("units", timeStr);

    // Scan start time (CCSDS)

    // Seconds since 1958
    var = gid.getVar("scan_start_CCSDS_sec");
    var.putVar(start, count, ccsdsStartScanSeconds.data());

    // Microseconds
    var = gid.getVar("scan_start_CCSDS_usec");
    var.putVar(start, count, ccsdsStartScanSubSeconds.data());

    // Extract and write HAM side
    vector<uint8_t> hamSide(scanNum);


    // fill value HAM_side data. writeTelemetryData() is where the actual data is written
    for (size_t i = 0; i < scanNum; i++) {
        // This part is modified to fix bogus values for the last two records in last chunk of L1A
        int ham = 255;
        hamSide[i] = (uint8_t)ham;
    }
    var = gid.getVar("HAM_side");
    var.putVar(start, count, hamSide.data());

    // Extract and write instrument spin ID
    for (size_t i = 0; i < scanNum; i++) {
        uint32_t ui32;
        memcpy(&ui32, &ancillaryData[i * ANCSIZE + timeOffset], 4);
        swapc_bytes((char *)&ui32, sizeof(uint32_t), 1);
        spinID[i] = ui32;
    }
    var = gid.getVar("spin_ID");
    var.putVar(start, count, spinID);

    // Packet and line number sequence error flags
    var = gid.getVar("pseq_flag");
    var.putVar(start, count, sciPacketSequenceError);
    var = gid.getVar("line_flag");
    var.putVar(start, count, ccdLineError);

    return 0;
}

// Write to global attributes in the NetCDF file
int L1aFile::writeGlobalMetaData(stringstream &startTime, stringstream &endTime, bool isFileAppended,
                                 std::string l1aFileName, std::string startDirectionStr,
                                 std::string endDirectionStr, short dataType, uint16_t swirModeIndex,
                                 uint16_t cdsModeIndex) {
                                    
    string DATA_COLLECT_MODE[] = {"",
                                  "Earth Collect",
                                  "Dark Collect",
                                  "Solar Cal",
                                  "SPCA Cal",
                                  "Response Curve",
                                  "Lunar Cal",
                                  "Diagnostic",
                                  "Static",
                                  "Earth Spectral",
                                  "",
                                  "External Snapshot Trigger",
                                  "Internal Snapshot Trigger",
                                  "Earth Collect with Non-Baseline Spectral Bands"};

    string SWIR_DATA_MODES[] = {"Science", "Diagnostic", "Single-image raw", "Test pattern"};

    string CDS_MODES[] = {"CDS", "Reset", "Video"};


    // if the current file is being appended to an earlier file, update the 'end'
    // stuff to be the most recent files while retaining the "start" stuff.
    if (isFileAppended) {

        l1afile->putAtt("time_coverage_end", endTime.str() + "Z");
        l1afile->putAtt("endDirection", endDirectionStr.c_str());
    }
    
    // if not being appended, write everything
    else {
        l1afile->putAtt("time_coverage_start", startTime.str() + "Z");
        l1afile->putAtt("time_coverage_end", endTime.str() + "Z");
        l1afile->putAtt("product_name", l1aFileName.c_str());
        l1afile->putAtt("startDirection", startDirectionStr.c_str());
        l1afile->putAtt("endDirection", endDirectionStr.c_str());
        l1afile->putAtt("data_collect_mode", DATA_COLLECT_MODE[dataType].c_str());
        l1afile->putAtt("SWIR_data_mode", SWIR_DATA_MODES[swirModeIndex].c_str());
        l1afile->putAtt("CDS_mode", CDS_MODES[cdsModeIndex].c_str());
    }
    
    return 0;
}

// Extract datatypes, spatial, spectral data and aggrigation data to the output NetCDF file
int L1aFile::writeAncillaryData(DimensionShape* dimShape, uint32_t scanNum, uint8_t *ancillaryData) {

    // read uint16 values in and call swap bytes function to change byte order
    uint16_t uint16ByteSwapper = 0;

    // Extract ancillary telemetry status flags and error counts
    short offset = 12;
    vector<int16_t> ancillaryTelemetries(6 * scanNum, BAD_INT);
    for (size_t i = 0; i < scanNum; i++) {
        for (size_t j = 0; j < 6; j++) {
            memcpy(&uint16ByteSwapper, &ancillaryData[i * ANCSIZE + offset + j * 2], 2);
            swapc_bytes((char *)&uint16ByteSwapper, sizeof(int16_t), 1);
            ancillaryTelemetries[i * 6 + j] = uint16ByteSwapper;
        }
    }

    // Extract spatial aggregation / data type table
    short extractedDataTypes[10];
    short extractedSpatialAggVals[10];
    short spatialZoneLineNum[10];

    // 1:1, 2:1, 4:1, 8:1
    short spatialAggRatios[4] = {1, 2, 4, 8};
    offset = 36;

    int32_t year, day;
    double startTime;
    getAncillaryPacketTime(&ancillaryData[0], year, day, startTime);
    uint32_t jayDay = jday(year, 1, day);

    for (size_t i = 0; i < 10; i++) {
        extractedDataTypes[i] = ancillaryData[offset + 3] % 16;
        if ((jayDay < 2459000) && (extractedDataTypes[i] == 11))
            extractedDataTypes[i] = 9;  // data type mod for ETU before June 2020
        extractedSpatialAggVals[i] = ancillaryData[offset + 2] % 4;
        if (extractedDataTypes[i] != 0)
            extractedSpatialAggVals[i] = spatialAggRatios[extractedSpatialAggVals[i]];
        spatialZoneLineNum[i] = ancillaryData[offset] * 256 + ancillaryData[offset + 1];
        offset += 4;
    }
    offset += 4;

    // Extract spectral aggregation and compute numbers of bands

    // Tap enable flags
    uint16_t blueTap = ancillaryData[offset + 2] * 256 + ancillaryData[offset + 3];
    uint16_t redTap = ancillaryData[offset + 0] * 256 + ancillaryData[offset + 1];

    // read in uint32 and swap bytes to change the order
    uint32_t uint32ByteSwapper;
    int32_t int32ByteSwapper;

    memcpy(&uint32ByteSwapper, &ancillaryData[offset + 8], 4);
    swapc_bytes((char *)&uint32ByteSwapper, sizeof(uint32_t), 1);
    uint32_t blueAggregationFactor = uint32ByteSwapper;

    memcpy(&uint32ByteSwapper, &ancillaryData[offset + 4], 4);
    swapc_bytes((char *)&uint32ByteSwapper, sizeof(uint32_t), 1);
    uint32_t redAggregationFactor = uint32ByteSwapper;

    int16_t blueSpectralAggregations[16] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    int16_t redSpectralAggregations[16] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    // there is a total of 16 taps, currTapBit is used as bitwise& to check each tap from the least sig bit to
    // the most significant (right to left). After each tap iteration, multiply by 2 to get the next bit
    // 1  = 0000 0000 0000 0001 : checks the right most bit for blue or red tap and checks if it is enabled
    // 2  = 0000 0000 0000 0010
    // 4  = 0000 0000 0000 0100
    // ...
    // 65536 = 1000 0000 0000 0001
    uint16_t currTapBit = 1;

    // works the same as currTapBit, except you are extracting the aggrigation factor which are 2 bits and bc
    // the aggrigation factors are 32 bits. Each iteration gets multiplied by 4.
    // 3  = 0000 0000 0000 0000 0000 0000 0000 0011
    // 12 = 0000 0000 0000 0000 0000 0000 0000 1100
    // 48 = 0000 0000 0000 0000 0000 0000 0011 0000
    // ...
    // 12884901888 = 1100 0000 0000 0000 0000 0000 0000 0000
    uint32_t currAggBit = 3;

    // works with currAggBit. After grabbing the aggigration factor, divide that number by this to extract the
    // bits. This will be factors of 4, so it will move 2 bits to the right
    // so if we have the aggigration factor at 12, the 2nd iteration, extract bit is 4.12/4 = 3
    // 12 = 0000 0000 0000 0000 0000 0000 0000 1100
    // 3  = 0000 0000 0000 0000 0000 0000 0000 0011 : after shifting it 2 bytes to the right
    uint32_t currAggFactorExtractBit = 1;

    for (size_t i = 0; i < 16; i++) {  // Put tap information in ascending spectral order
        uint16_t isBlueCcdTapsEnabled = (blueTap & currTapBit) / currTapBit;
        if (isBlueCcdTapsEnabled)
            blueSpectralAggregations[15 - i] =
                spatialAggRatios[(blueAggregationFactor & currAggBit) / currAggFactorExtractBit];
        uint16_t isRedCcdTapsEnabled = (redTap & currTapBit) / currTapBit;
        if (isRedCcdTapsEnabled)
            redSpectralAggregations[15 - i] =
                spatialAggRatios[(redAggregationFactor & currAggBit) / currAggFactorExtractBit];

        currTapBit *= 2;
        currAggBit *= 4;
        currAggFactorExtractBit *= 4;
    }

    NcGroup ncGroupId;
    NcVar var;
    vector<size_t> start;
    vector<size_t> count;
    start.push_back(0);

    // Write to file

    ncGroupId = l1afile->getGroup("spatial_spectral_modes");

    // start at index 0 and input 10 items (size of spatial_zones)
    // for the next 3 variables
    start = {0}; 
    count = {10};

    // Data type
    var = ncGroupId.getVar("spatial_zone_data_type");
    var.putVar(start, count, extractedDataTypes);

    // Spatial aggregation
    var = ncGroupId.getVar("spatial_aggregation");
    var.putVar(start, count, extractedSpatialAggVals);

    // Number of lines
    var = ncGroupId.getVar("spatial_zone_lines");
    var.putVar(start, count, spatialZoneLineNum);



    // start at index 0, for 16 items (size of number_of_taps) for the next 2 variables
    // NOTE: did not need to reassign start to 0, but to make it more readable I reassigned it
    //       so you do not need to scroll up to see what the start was
    start = {0};
    count = {16};

    // Blue spectral aggregation
    var = ncGroupId.getVar("blue_spectral_mode");
    var.putVar(start, count, blueSpectralAggregations);

    // Red spectral aggregation
    var = ncGroupId.getVar("red_spectral_mode");
    var.putVar(start, count, redSpectralAggregations);

    // Extract MCE data and write to file
    // Use ancillary packets AFTER science data

    ncGroupId = l1afile->getGroup("engineering_data");

    // write ancillary telemetry status flags and error counts
    start = {dimShape->numScans, 0}; // start at the current num of scans, pos 0
    count = {scanNum, 6}; // number of scans row of data with 6 elements in each row
    var = ncGroupId.getVar("ancillary_tlm");
    var.putVar(start, count, ancillaryTelemetries.data());


    // start and count for aggrigation control variables:
    // agg_control, blue_agg_error, red_agg_error and dig_card_error
    start = {dimShape->numScans};
    count = {scanNum};

    // Write aggregation controls
    vector<int32_t> aggregationControls(scanNum);
    for (size_t i = 0; i < scanNum; i++) {
        memcpy(&int32ByteSwapper, &ancillaryData[i * ANCSIZE + offset - 4], 4);
        swapc_bytes((char *)&int32ByteSwapper, sizeof(int32_t), 1);
        aggregationControls[i] = int32ByteSwapper;
    }
    var = ncGroupId.getVar("agg_control");
    var.putVar(start, count, aggregationControls.data());

    // Aggregation errors
    vector<uint16_t> blueAggErrors(scanNum);
    vector<uint16_t> redAggErrors(scanNum);

    for (size_t i = 0; i < scanNum; i++) {
        memcpy(&uint16ByteSwapper, &ancillaryData[(i + 1) * ANCSIZE + offset + 12], 2);
        swapc_bytes((char *)&uint16ByteSwapper, sizeof(int16_t), 1);
        blueAggErrors[i] = uint16ByteSwapper;

        memcpy(&uint16ByteSwapper, &ancillaryData[(i + 1) * ANCSIZE + offset + 14], 2);
        swapc_bytes((char *)&uint16ByteSwapper, sizeof(int16_t), 1);
        redAggErrors[i] = uint16ByteSwapper;
    }

    // BLUE
    var = ncGroupId.getVar("blue_agg_error");
    var.putVar(start, count, blueAggErrors.data());

    // RED
    var = ncGroupId.getVar("red_agg_error");
    var.putVar(start, count, redAggErrors.data());

    // Digital card error status
    vector<int32_t> digitalCardErrors(scanNum);
    for (size_t i = 0; i < scanNum; i++) {
        memcpy(&int32ByteSwapper, &ancillaryData[(i + 1) * ANCSIZE + offset + 20], 4);
        swapc_bytes((char *)&int32ByteSwapper, sizeof(int32_t), 1);
        digitalCardErrors[i] = int32ByteSwapper;
    }

    var = ncGroupId.getVar("dig_card_error");
    var.putVar(start, count, digitalCardErrors.data());

    return 0;
}

int L1aFile::writeTelemetryData(DimensionShape* dimShape, spatialAggTable *spatialAggList, uint32_t numTelemetryPackets,
                                uint8_t (*telemetryData)[TLMSIZE], int32_t *spinID, uint16_t &cdsMode,
                                uint32_t scanNum, const AncillaryPktTimeStamp &starttime) {
    // numTelemetryPackets:  Number of telemetry packets
    // telemetryData: OCI telemetry data packets
    // cdsMode: CDS mode (0 = enabled, 1 = reset, 2 = video)
    // telemetryZoneStartTimes(2): Telmetry zone start times (msec)
    // telemetryZoneDurations(2): Telmetry zone durations (msec)
    // tdiTime: CCD data line TDI time (clock cycles)

    // Set up output array and extract/convert data from DAU packets
    const uint16_t DAU_APID = 723;
    const uint16_t DDC_TELEMETRY_APID = 701;
    const uint16_t RTA_HAM_MCE_TELEMETRY_APID = 713;  // 711->713, changed by Liang Hong, 10/29/2020
    const uint16_t RTA_HAM_MCE_ENCODER_APID = 712;
    const uint16_t SCA_MCE_TELEMETRY_APID = 717;
    const uint16_t SCA_MCE_ENCODER_APID = 716;
    const uint16_t ICDU_TC_TELEMETRY_THERMISTORS_APID = 656;
    const uint16_t DDC_TABLE_TELEMETRY_APID = 703;
    const uint16_t FSW_DAUC_TEMP_APID = 744;
    const uint16_t FSW_ICCU_MCE_TELEMETRY_TEMP_APID = 745;

    // dimension sizes
    const uint16_t NUM_ICDU_THERMISTERS = 74;
    const uint16_t NUM_DAUC_TEMPS = 69;
    const uint16_t NUM_ICDU_TEMPS = 16;
    const uint16_t NUM_ADC_LATENCY = 4;
    
    // dimensions used for writing to file
    const size_t DDC_TLM = 524;
    const size_t MCE_BLOCKS = 480;
    const size_t ENCODER_SAMPLES = 200;
    const size_t ENCODER_CHANNELS = 4;

    uint16_t dauTelemetryTimePktCount = 0;
    uint16_t ddcTelemetryTimePktCount = 0;
    uint16_t mceTelemetryPktCount = 0;
    uint16_t rtaHamMceEncoderPktCount = 0;
    uint16_t scaTelemetryTimePktCount = 0;
    uint16_t scaMceEncoderPktCount = 0;
    uint16_t tcTelemetryTimePktCount = 0;
    uint16_t daucTempTimePktCount = 0;
    bool hasDdcTableTelemetry = false;
    uint16_t icduMceTempTimePktCount = 0;
    vector<int16_t> telemetryZoneStartTimes(2);
    vector<int16_t> telemetryZoneDurations(2);
    uint16_t tdiTime = 0;

    double l1afile_epoch = yds2unix(starttime.year, starttime.day, 0.0);
    string l1afile_epoch_str = unix2isodate(l1afile_epoch, 'G');
    l1afile_epoch_str = l1afile_epoch_str.substr(0, 10);
    l1afile_epoch_str.insert(0, "seconds since ");

    // dimensions required to make the correct vector size
    const int NUM_DAU_TELEMETRY = 620;
    const int NUM_DDC_TELEMETRY = 524;
    const int NUM_MCE_BLOCKS = 480;  // MCE == Mechanism Control Electronics
    const int NUM_ENCODER_SAMPLES = 200;
    const int NUM_ENCODER_CHANNELS = 4;
    const int NUM_LINE_SKIPS = 33;
    const int NUM_TC_TELEMETRY = 1216;  // TC == Temperature Controller
    const int NUM_ICDU_MCE_TEMP_TELEMETRY = 76;

    // initialize arrays to hold the extracted data
    vector<double> dauTelemetryTimes(numTelemetryPackets, BAD_FLT);
    vector<int32_t> dauSpinIds(numTelemetryPackets, BAD_INT);
    vector<uint8_t> dauTelemetries(numTelemetryPackets * NUM_DAU_TELEMETRY, 0);
    vector<double> ddcTelemetryTimes(numTelemetryPackets, BAD_FLT);
    vector<uint8_t> ddcTelemetries(numTelemetryPackets * NUM_DDC_TELEMETRY, 0);
    vector<int32_t> mceSpinIds(numTelemetryPackets, BAD_INT);
    vector<uint8_t> mceTelemetries(numTelemetryPackets * NUM_MCE_BLOCKS, 0);
    vector<int32_t> encoderSpinIds(numTelemetryPackets, BAD_INT);
    vector<int16_t> mceEncoderData(numTelemetryPackets * NUM_ENCODER_CHANNELS * NUM_ENCODER_SAMPLES, BAD_INT);
    vector<int16_t> auxilaryParameters(NUM_LINE_SKIPS, BAD_INT);
    vector<int32_t> scaSpinIds(numTelemetryPackets, BAD_INT);
    vector<uint8_t> scaTelemetries(numTelemetryPackets * NUM_MCE_BLOCKS, 0);
    vector<double> scaTelemetryTimes(numTelemetryPackets, BAD_FLT);
    vector<float> scaDiffuserPositions(numTelemetryPackets, BAD_FLT);
    vector<int32_t> scaEncoderSpinIds(numTelemetryPackets, BAD_INT);
    vector<int16_t> scaEncoderData(numTelemetryPackets * NUM_ENCODER_CHANNELS * NUM_ENCODER_SAMPLES, BAD_INT);
    vector<double> tcTelemetryTime(numTelemetryPackets, BAD_FLT);
    vector<uint8_t> tcTelemetries(numTelemetryPackets * NUM_TC_TELEMETRY, 0);
    vector<double> daucTempTimes(numTelemetryPackets, BAD_FLT);
    vector<double> icduMceTempTime(numTelemetryPackets, BAD_FLT);
    vector<uint8_t> icduMceTempTelemetries(numTelemetryPackets * NUM_ICDU_MCE_TEMP_TELEMETRY, 0);
    vector<float> icduThermisters(numTelemetryPackets * NUM_ICDU_THERMISTERS, BAD_FLT);
    vector<float> daucTemps(numTelemetryPackets * NUM_DAUC_TEMPS, BAD_FLT);
    vector<float> icduMceTemps(numTelemetryPackets * NUM_ICDU_TEMPS, BAD_FLT);
    vector<uint8_t> adcLatencies(numTelemetryPackets * NUM_ADC_LATENCY, 0);
    vector<uint8_t> cdsDisableFlags(numTelemetryPackets, 0);
    vector<uint8_t> hamSide(numTelemetryPackets);

    uint32_t redMask[16];
    uint32_t blueMask[16];
    const double SCA_DIFFUSER_POSITION_SCALE_FACTOR = 360.0 / pow(2, 32);

    for (size_t i = 0; i < numTelemetryPackets; i++) {
        uint32_t apid = ((uint8_t)telemetryData[i][0] % 8) * 256 + (uint8_t)telemetryData[i][1];
        uint8_t ccsdsTimes[8] = {0, 0, 0, 0, 0, 0, 0, 0};

        uint32_t uint32ByteSwapper;  // temp variable to read into and swap the byte order of uint32
        int16_t int16ByteSwapper;    // temp used to store a int16 to swap byte order
        float float32Reverser;       // temp used to store float34 to reverse it
        ccsdsTimes[6] = 0;
        ccsdsTimes[7] = 0;

        switch (apid) {
            case DAU_APID:
                memcpy(ccsdsTimes, &telemetryData[i][6], 6);
                dauTelemetryTimes[dauTelemetryTimePktCount] = ccsds_sec_to_unix(ccsdsTimes) - l1afile_epoch;
                swapc_bytes2((char*)&telemetryData[i][12], (char *)&uint32ByteSwapper, sizeof(uint32_t), 1);
                dauSpinIds[dauTelemetryTimePktCount] = uint32ByteSwapper;
                memcpy(&dauTelemetries[dauTelemetryTimePktCount * NUM_DAU_TELEMETRY], &telemetryData[i][16],
                       NUM_DAU_TELEMETRY);
                dauTelemetryTimePktCount++;
                break;

            case DDC_TELEMETRY_APID:
                // DDC telemetry
                memcpy(ccsdsTimes, &telemetryData[i][6], 6);
                ddcTelemetryTimes[ddcTelemetryTimePktCount] = ccsds_sec_to_unix(ccsdsTimes) - l1afile_epoch;
                memcpy(&ddcTelemetries[ddcTelemetryTimePktCount * 524], &telemetryData[i][12], 524);
                cdsDisableFlags[ddcTelemetryTimePktCount] = telemetryData[i][29];
                memcpy(&adcLatencies[ddcTelemetryTimePktCount * NUM_ADC_LATENCY], &telemetryData[i][176],
                       NUM_ADC_LATENCY);
                if (ddcTelemetryTimePktCount == 0) {
                    for (size_t j = 0; j < 16; j++) {

                        swapc_bytes2((char*)&telemetryData[i][196 + j * 4], (char *)&uint32ByteSwapper, sizeof(uint32_t), 1);
                        redMask[j] = (uint32_t)uint32ByteSwapper;

                        swapc_bytes2((char*)&telemetryData[i][260 + j * 4], (char *)&uint32ByteSwapper, sizeof(uint32_t), 1);
                        blueMask[j] = (uint32_t)uint32ByteSwapper;
                    }
                }
                ddcTelemetryTimePktCount++;
                break;

            case RTA_HAM_MCE_TELEMETRY_APID:
                // RTA/HAM MCE telemetry
                swapc_bytes2((char*)&telemetryData[i][12], (char *)&uint32ByteSwapper, sizeof(uint32_t), 1);
                mceSpinIds[mceTelemetryPktCount] = uint32ByteSwapper;
                memcpy(&mceTelemetries[mceTelemetryPktCount * NUM_MCE_BLOCKS], &telemetryData[i][16],
                       NUM_MCE_BLOCKS);
                hamSide[mceTelemetryPktCount] = (telemetryData[i][49] & 8) / 8;
                mceTelemetryPktCount++;
                break;

            case RTA_HAM_MCE_ENCODER_APID:
                // RTA/HAM MCE encoder data
                swapc_bytes2((char*)&telemetryData[i][12], (char *)&uint32ByteSwapper, sizeof(uint32_t), 1);
                encoderSpinIds[rtaHamMceEncoderPktCount] = uint32ByteSwapper;

                for (size_t j = 0; j < NUM_ENCODER_CHANNELS * NUM_ENCODER_SAMPLES; j++) {
                    swapc_bytes2((char*)&telemetryData[i][16 + 2 * j], (char *)&int16ByteSwapper, sizeof(int16_t), 1);
                    mceEncoderData[rtaHamMceEncoderPktCount * NUM_ENCODER_CHANNELS * NUM_ENCODER_SAMPLES +
                                   j] = int16ByteSwapper;
                }
                rtaHamMceEncoderPktCount++;
                break;

            case SCA_MCE_TELEMETRY_APID:
                // SCA MCE telemetry
                memcpy(ccsdsTimes, &telemetryData[i][6], 6);
                scaTelemetryTimes[scaTelemetryTimePktCount] = ccsds_sec_to_unix(ccsdsTimes) - l1afile_epoch;
                swapc_bytes2((char*)&telemetryData[i][12], (char *)&uint32ByteSwapper, sizeof(uint32_t), 1);
                scaSpinIds[scaTelemetryTimePktCount] = uint32ByteSwapper;
                memcpy(&scaTelemetries[scaTelemetryTimePktCount * NUM_MCE_BLOCKS], &telemetryData[i][16],
                       NUM_MCE_BLOCKS);
                swapc_bytes2((char*)&telemetryData[i][348], (char *)&uint32ByteSwapper, sizeof(uint32_t), 1);
                float sca_abs_enc_count;
                sca_abs_enc_count =
                    (float)uint32ByteSwapper;  // Convert byte array to unsigned long and little-endian
                scaDiffuserPositions[scaTelemetryTimePktCount] =
                    330.0 -
                    sca_abs_enc_count * SCA_DIFFUSER_POSITION_SCALE_FACTOR;  // Perform conversion from Capon
                scaTelemetryTimePktCount++;
                break;

            case SCA_MCE_ENCODER_APID:
                // SCA MCE encoder data
                swapc_bytes2((char*)&telemetryData[i][12], (char *)&uint32ByteSwapper, sizeof(uint32_t), 1);
                scaEncoderSpinIds[scaMceEncoderPktCount] = uint32ByteSwapper;
                for (size_t j = 0; j < 4 * 200; j++) {
                    swapc_bytes2((char*)&telemetryData[i][16 + 2 * j], (char *)&int16ByteSwapper, sizeof(int16_t), 1);
                    scaEncoderData[scaMceEncoderPktCount * NUM_ENCODER_CHANNELS * NUM_ENCODER_SAMPLES + j] =
                        int16ByteSwapper;
                }
                scaMceEncoderPktCount++;
                break;

            case ICDU_TC_TELEMETRY_THERMISTORS_APID:
                // ICDU TC telemetry and thermistors
                memcpy(ccsdsTimes, &telemetryData[i][6], 6);
                tcTelemetryTime[tcTelemetryTimePktCount] = ccsds_sec_to_unix(ccsdsTimes) - l1afile_epoch;
                memcpy(&tcTelemetries[tcTelemetryTimePktCount * NUM_TC_TELEMETRY], &telemetryData[i][12],
                       NUM_TC_TELEMETRY);
                for (size_t j = 0; j < NUM_ICDU_THERMISTERS; j++) {
                    memcpy(&float32Reverser, &telemetryData[i][184 + 4 * j], 4);
                    swapc_bytes2((char*)&telemetryData[i][184 + 4 * j], (char *)&float32Reverser, sizeof(float), 1);
                    icduThermisters[tcTelemetryTimePktCount * NUM_ICDU_THERMISTERS + j] =
                        float32Reverser;
                }
                tcTelemetryTimePktCount++;
                break;

            case DDC_TABLE_TELEMETRY_APID:
                // DDC table telemetry
                // Only need one set of these
                if (!hasDdcTableTelemetry) {  // dont have it yet, grab it
                    for (size_t j = 0; j < 33; j++) {
                        swapc_bytes2((char*)&telemetryData[i][64 + j * 4], (char *)&uint32ByteSwapper, sizeof(uint32_t), 1);
                        auxilaryParameters[j] = (int16_t)uint32ByteSwapper;
                    }
                    hasDdcTableTelemetry = true;  // never grab it again
                }
                break;

            case FSW_DAUC_TEMP_APID:
                // FSW DAUC temperatures
                memcpy(ccsdsTimes, &telemetryData[i][6], 6);
                daucTempTimes[daucTempTimePktCount] = ccsds_sec_to_unix(ccsdsTimes) - l1afile_epoch;
                for (size_t j = 0; j < NUM_DAUC_TEMPS; j++) {
                    swapc_bytes2((char*)&telemetryData[i][12 + 4 * j], (char *)&float32Reverser, sizeof(float), 1);
                    daucTemps[daucTempTimePktCount * NUM_DAUC_TEMPS + j] = float32Reverser;
                }
                daucTempTimePktCount++;
                break;

            case FSW_ICCU_MCE_TELEMETRY_TEMP_APID:
                // FSW ICCU/MCE telemetry and temperatures
                memcpy(ccsdsTimes, &telemetryData[i][6], 6);
                icduMceTempTime[icduMceTempTimePktCount] = ccsds_sec_to_unix(ccsdsTimes) - l1afile_epoch;
                memcpy(&icduMceTempTelemetries[icduMceTempTimePktCount * NUM_ICDU_MCE_TEMP_TELEMETRY],
                       &telemetryData[i][12], NUM_ICDU_MCE_TEMP_TELEMETRY);
                for (size_t j = 0; j < NUM_ICDU_TEMPS; j++) {
                    swapc_bytes2((char*)&telemetryData[i][24 + 4 * j], (char *)&float32Reverser, sizeof(float), 1);
                    icduMceTemps[icduMceTempTimePktCount * NUM_ICDU_TEMPS + j] =
                        float32Reverser;
                }
                icduMceTempTimePktCount++;
                break;
        }
    }  // tlm loop

    NcGroup engineeringDataGroup = l1afile->getGroup("engineering_data");
    NcVar var;  // rererence var from the engineering data group

    // define where and how much to write into the file
    vector<size_t> start;
    vector<size_t> count;

    if (dauTelemetryTimePktCount == 0)
        dauTelemetryTimePktCount = 1;  // 0.99.20

    start = {dimShape->tlmPackets}; // get start position of tlm packet
    count = {dauTelemetryTimePktCount};

    var = engineeringDataGroup.getVar("DAU_tlm_time");
    var.putVar(start, count, dauTelemetryTimes.data());
    var.putAtt("units", l1afile_epoch_str);

    var = engineeringDataGroup.getVar("DAU_spin_ID");
    var.putVar(start, count, dauSpinIds.data());

    // for 2D, the start position of the second dim is always 0 since it's new data on
    // the next empty space
    start = {dimShape->tlmPackets, 0};
    count = {dauTelemetryTimePktCount, 620};

    var = engineeringDataGroup.getVar("DAU_telemetry");
    var.putVar(start, count, dauTelemetries.data());

    // Get telemetry zone start times and durations
    telemetryZoneStartTimes[0] = dauTelemetries[123];
    telemetryZoneStartTimes[1] = dauTelemetries[125];
    telemetryZoneDurations[0] = dauTelemetries[122];
    telemetryZoneDurations[1] = dauTelemetries[124];
    //}

    if (ddcTelemetryTimePktCount == 0) {
        for (size_t i = 0; i < 16; i++) {
            redMask[i] = 4294967295;
            blueMask[i] = 4294967295;
        }
        ddcTelemetryTimePktCount = 1;  // 0.99.20
    }
    
    // DDC and DAU share the same tlm_packets dim. If one is bigger, the smller one
    // will have a gap with fill values.
    // So if 5 vs. 3, then the data will look like this in the NetCDF
    // Index  : 0 1 2 3 4 5 6 ... n
    // DDC (5): 0 0 0 0 0 
    // DAU (3): 0 0 0 _ _ 
    // For the next append, it will start at index 5

    start = {dimShape->tlmPackets};
    count = {ddcTelemetryTimePktCount};

    // DDC tlm time
    var = engineeringDataGroup.getVar("DDC_tlm_time");
    var.putVar(start, count, ddcTelemetryTimes.data());
    var.putAtt("units", l1afile_epoch_str);

    // CDS disable flag
    var = engineeringDataGroup.getVar("CDS_disable");
    var.putVar(start, count, cdsDisableFlags.data());


    // DDC telemetry
    start = {dimShape->tlmPackets, 0};
    count = {ddcTelemetryTimePktCount, DDC_TLM};
    var = engineeringDataGroup.getVar("DDC_telemetry");
    var.putVar(start, count, ddcTelemetries.data());

    // ADC latency
    start = {dimShape->tlmPackets, 0};
    count = {ddcTelemetryTimePktCount, NUM_ADC_LATENCY};
    var = engineeringDataGroup.getVar("ADC_latency");
    var.putVar(start, count, adcLatencies.data());


    // Red and blue channel masks
    start = {0};
    count = {16};
    var = engineeringDataGroup.getVar("blue_channel_mask");
    var.putVar(start, count, blueMask);
    var = engineeringDataGroup.getVar("red_channel_mask");
    var.putVar(start, count, redMask);

    // Get CDS mode for metadata
    cdsMode = cdsDisableFlags[0] * (adcLatencies[0] - 14);

    //  Get TDI time
    uint16_t uint16ByteSwapper;
    memcpy(&uint16ByteSwapper, &ddcTelemetries[346], 2);
    swapc_bytes((char *)&uint16ByteSwapper, sizeof(int16_t), 1);
    tdiTime = uint16ByteSwapper;

    if (mceTelemetryPktCount == 0) {
        mceTelemetryPktCount = 1;  // 0.99.20
    }


    // RTA/HAM MCE spin ID
    start = {dimShape->numMceScans};
    count = {mceTelemetryPktCount};
    var = engineeringDataGroup.getVar("MCE_spin_ID");
    var.putVar(start, count, mceSpinIds.data());

    // HAM side
    vector<uint8_t> mirrorSide(scanNum, 255);
    for (size_t i = 0; i < scanNum; i++) {
        for (size_t j = 0; j < numTelemetryPackets; j++) {
            if ((int)mceSpinIds[j] == spinID[i]) {
                mirrorSide[i] = hamSide[j];
                break;
            }
            // if mcespinId != spinID AND the curr scan is 0, skip and leave it at 255
            // when i==0 this should almost never happen because MCE Packets are always a spin
            // behind the science data -Fred
            else if (i == 0) {
                continue;
            } else {
                // reason i==0 is a case above is due to i-1. This is for all cases
                // where i=1 or more.
                mirrorSide[i] = 1 - mirrorSide[i - 1];
            }
        }
    }

    start = {dimShape->numScans};
    count = {scanNum};
    NcGroup scgid = l1afile->getGroup("scan_line_attributes");
    var = scgid.getVar("HAM_side");
    var.putVar(start, count, mirrorSide.data());

    
    // MCE_telemetry
    start = {dimShape->numMceScans, 0};
    count = {mceTelemetryPktCount, MCE_BLOCKS};
    var = engineeringDataGroup.getVar("MCE_telemetry");
    var.putVar(start, count, mceTelemetries.data());


    if (rtaHamMceEncoderPktCount == 0)
        rtaHamMceEncoderPktCount = 1;  // 0.99.20
        

    // RTA/HAM encoder spin ID
    start = {dimShape->numMceScans};
    count = {rtaHamMceEncoderPktCount};
    var = engineeringDataGroup.getVar("encoder_spin_ID");
    var.putVar(start, count, encoderSpinIds.data());


    // RTA/HAM encoder data
    start = {dimShape->numMceScans, 0, 0}; // start at current mce scan and put data at index 0s
    count = {rtaHamMceEncoderPktCount, ENCODER_SAMPLES, ENCODER_CHANNELS}; 
    var = engineeringDataGroup.getVar("MCE_encoder_data");
    var.putVar(start, count, mceEncoderData.data());
    
    /// --- UP TO THIS POINT, number_of_mce_scans dimension is finished --- 
    /// --- UPDATE THE DIMENSION SHAPE FOR IT BASED ON WHICH IS LARGER
    /// --- BETWEEN rtaHamMceEncoderPktCount and mceTelemetryPktCoun
    dimShape->incrementNumMceScanShape(rtaHamMceEncoderPktCount, mceTelemetryPktCount);
 
    
    if (scaTelemetryTimePktCount == 0) {
        scaTelemetryTimePktCount = 1;  // 1.01.00
    }

    // SCA MCE spin ID
    var = engineeringDataGroup.getVar("SCA_spin_ID");
    start = {dimShape->numScaScans};
    count = {scaTelemetryTimePktCount};
    var.putVar(start, count, scaSpinIds.data());

    // SCA_telemetry
    start = {dimShape->numScaScans, 0};
    count = {scaTelemetryTimePktCount, MCE_BLOCKS};
    var = engineeringDataGroup.getVar("SCA_telemetry");
    var.putVar(start, count, scaTelemetries.data());

    // SCA telemetry time
    start = {dimShape->numScaScans};
    count = {scaTelemetryTimePktCount};
    var = engineeringDataGroup.getVar("SCA_tlm_time");
    var.putVar(start, count, scaTelemetryTimes.data());
    var.putAtt("units", l1afile_epoch_str);

    // SCA diffuser position
    // start and count == SCA telemetry time
    var = engineeringDataGroup.getVar("SCA_diffuser_position");
    var.putVar(start, count, scaDiffuserPositions.data()); 

    if (scaMceEncoderPktCount == 0) {
        scaMceEncoderPktCount = 1;
    }

    // RTA/HAM encoder spin ID
    start = {dimShape->numScaScans};
    count = {scaMceEncoderPktCount};
    var = engineeringDataGroup.getVar("SCA_encoder_spin_ID");
    var.putVar(start, count, scaEncoderSpinIds.data());

    // RTA/HAM encoder data
    start = {dimShape->numScaScans, 0, 0};
    count = {scaMceEncoderPktCount, ENCODER_SAMPLES, ENCODER_CHANNELS};
    var = engineeringDataGroup.getVar("SCA_encoder_data");
    var.putVar(start, count, scaEncoderData.data());


    /// --- UP TO THIS POINT, number_of_mce_scans dimension is finished --- 
    /// --- UPDATE THE DIMENSION SHAPE FOR IT BASED ON WHICH IS LARGER
    dimShape->incrementNumScaScanShape(scaTelemetryTimePktCount, scaMceEncoderPktCount);

    // additional tlm_packet variable

    if (tcTelemetryTimePktCount == 0)
        tcTelemetryTimePktCount = 1;  // 0.99.20
  

    // TC tlm time
    start = {dimShape->tlmPackets};
    count = {tcTelemetryTimePktCount};
    var = engineeringDataGroup.getVar("TC_tlm_time");
    var.putVar(start, count, tcTelemetryTime.data());
    var.putAtt("units", l1afile_epoch_str);

    // TC telemetry
    const size_t TC_TLM = 1216; // Only TC_telemetry uses this
    start = {dimShape->tlmPackets, 0};
    count = {tcTelemetryTimePktCount, TC_TLM};
    var = engineeringDataGroup.getVar("TC_telemetry");
    var.putVar(start, count, tcTelemetries.data());

    // ICDU thermisters
    start = {dimShape->tlmPackets, 0};
    count = {tcTelemetryTimePktCount, NUM_ICDU_THERMISTERS};
    var = engineeringDataGroup.getVar("ICDU_thermisters");
    var.putVar(start, count, icduThermisters.data());


    // Auxilary parameter table
    engineeringDataGroup = l1afile->getGroup("spatial_spectral_modes");
    var = engineeringDataGroup.getVar("aux_param_table");
    start = {0};
    count = {33}; // lin_skips dim
    var.putVar(start, count, auxilaryParameters.data());


    if (daucTempTimePktCount == 0)
        daucTempTimePktCount = 1;

    engineeringDataGroup = l1afile->getGroup("engineering_data");

    // FSW DAUC temperature time
    start = {dimShape->tlmPackets};
    count = {daucTempTimePktCount};
    var = engineeringDataGroup.getVar("DAUC_temp_time");
    var.putVar(start, count, daucTempTimes.data());
    var.putAtt("units", l1afile_epoch_str);


    // FSW DAUC temperatures
    var = engineeringDataGroup.getVar("DAUC_temperatures");
    start = {dimShape->tlmPackets, 0};
    count = {daucTempTimePktCount, NUM_DAUC_TEMPS};
    var.putVar(start, count, daucTemps.data());


    if (icduMceTempTimePktCount == 0)
        icduMceTempTimePktCount = 1;

    // FSW ICDU/MCE temperature time
    start = {dimShape->tlmPackets};
    count = {icduMceTempTimePktCount};
    var = engineeringDataGroup.getVar("ICDU_MCE_temp_time");
    var.putVar(start, count, icduMceTempTime.data());
    var.putAtt("units", l1afile_epoch_str);


    // FSW ICDU/MCE temperature telemetry
    const size_t ANCIL_TLM = 76;
    start = {dimShape->tlmPackets, 0};
    count = {icduMceTempTimePktCount, ANCIL_TLM};
    var = engineeringDataGroup.getVar("ICDU_MCE_temp_tlm");
    var.putVar(start, count, icduMceTempTelemetries.data());

    // FSW ICDU/MCE temperatures
    var = engineeringDataGroup.getVar("ICDU_MCE_temperatures");
    // start is the same as ICDU_MCE_temp_tlm, just diff count
    count = {icduMceTempTimePktCount, NUM_ICDU_TEMPS};
    var.putVar(start, count, icduMceTemps.data());


    /// --- UP TO THIS POINT, tlm_packets dimension is finished --- 
    /// --- UPDATE THE DIMENSION SHAPE FOR IT BASED ON WHICH IS LARGER of the 5:
    /// --- ddcTelemetryTimePktCount, dauTelemetryTimePktCount, tcTelemetryTimePktCount, 
    /// --- daucTempPktCount and icduMceTempTimePktCount

    // find the largest tlm_packet count and increment the next write position by it
    // all variables listed here uses the tlm_packets dim. 
    size_t tlmPacketIncrementVal = max({
        ddcTelemetryTimePktCount, 
        dauTelemetryTimePktCount, 
        tcTelemetryTimePktCount,
        daucTempTimePktCount,
        icduMceTempTimePktCount
    });
    dimShape->incrementTlmPacketsShape(tlmPacketIncrementVal);


    // check the science data and telemetry collection zones for conflicts and write the results
    // as an attribute to the engineering data group.  There are two telemetry collection zones
    string conflict = "No";

    if ((telemetryZoneDurations[0] == 0) && (telemetryZoneDurations[1] == 0)) {
        cout << "No telemetry zone fields in file" << endl;
        engineeringDataGroup.putAtt("science_telemetry_zone_conflict", conflict);
    } else {
        float clock = 136000;                       // Master clock frequency in msec
        float msecPerLine = (tdiTime + 1) / clock;  // msec per line

        // Find no-data zones
        float *zones = new float[2 * 10];
        float *zonesCopy = new float[2 * 10];
        int zoneCount = 0;
        uint16_t line = 0;
        for (size_t i = 0; i < 10; i++) {
            if (((spatialAggList[i].dataType == 0) || (spatialAggList[i].dataType == 10)) &&
                (spatialAggList[i].lines > 0)) {
                zones[zoneCount] = line * msecPerLine;
                zones[1 * 10 + zoneCount] = (line + spatialAggList[i].lines) * msecPerLine;
                // Check for consecutive no-data zones
                if ((zoneCount > 0) && (zones[zoneCount] == zones[1 * 10 + zoneCount - 1])) {
                    zones[1 * 10 + zoneCount - 1] = zones[1 * 10 + zoneCount];
                } else {
                    zoneCount++;
                }
            }
            line += spatialAggList[i].lines;
        }
        memcpy(zonesCopy, zones, 20);
        // zones new dimension: 2*zoneCount
        for (int i = 0; i < zoneCount; i++)
            zones[zoneCount + i] = zones[10 + i];  // zones = zones(*,0:zoneCount-1)

        //  If first no-data zone is at start of spin, extend last zone by that amount.
        if ((spatialAggList[0].dataType == 0) || (spatialAggList[0].dataType == 10))
            zones[1 * zoneCount + zoneCount - 1] += zones[1 * zoneCount + 0];

        if (zoneCount > 0) {
            // Loop through telemetry zones
            for (size_t i = 0; i < 2; i++) {
                // Find overlapping no-data zone
                int16_t telemetryZoneEndTime = telemetryZoneStartTimes[i] + telemetryZoneDurations[i] + 3;
                int k;
                for (k = zoneCount - 1; k >= 0; k--) {
                    // k = max(where(telemetryZoneEndTime gt zones(0,*)))
                    if (telemetryZoneEndTime > zones[k])
                        break;
                }
                if ((k == -1) || (telemetryZoneStartTimes[i] < zones[k]) ||
                    (telemetryZoneEndTime > zones[1 * zoneCount + k])) {
                    conflict = "Yes";
                    cout << "Telemetry zone " << telemetryZoneStartTimes[i] << ", " << telemetryZoneEndTime
                         << endl;
                    cout << "No-data zone " << zones[k] << ", " << zones[1 * zoneCount + k] << endl;
                }
            }
        } else
            conflict = "Yes";

        engineeringDataGroup.putAtt("science_telemetry_zone_conflict", conflict);

        delete[] zones;
        delete[] zonesCopy;
    }  // if (telemetryZoneDurations[0] == 0) && (telemetryZoneDurations[1] == 0) else

    return 0;
}

// if epoch time for ancillary and hkt is different, sync them by adding their difference to time arr
void L1aFile::synchronizeEpochTime(vector<double> &navigationTimeArr, size_t arrSize,
                                   double &hktFileEpochTime, double &ancillaryEpochTime) {
    if (ancillaryEpochTime != hktFileEpochTime) {
        double epochDiff = hktFileEpochTime - ancillaryEpochTime;
        for (size_t i = 0; i < arrSize; i++) {
            navigationTimeArr[i] += epochDiff;
        }
    }
}

// Read the hktlist and open it to get the basic spacecraft navigation data like attitude time, speed etc.
int L1aFile::writeNavigationData(DimensionShape* dimShape, std::string hktList, AncillaryPktTimeStamp &startTime,
                                 AncillaryPktTimeStamp &endTime) {
    uint16_t maxNumNavRecords = 10000;  // max records possible for navigation data
    size_t attitudeRecCount = 0;
    size_t orbitRecCount = 0;
    size_t tiltRecCount = 0;

    // convert times to unix seconds
    double startTimeAsUnixSecs = yds2unix(startTime.year, startTime.day, startTime.second);
    double endTimeAsUnixSecs = yds2unix(endTime.year, endTime.day, endTime.second);
    double epochTimeAsUnixSecs = yds2unix(startTime.year, startTime.day, 0.0);

    // navigation_data group's time units: seconds since YYYY-MM-DD
    string epochIsodateStr = unix2isodate(epochTimeAsUnixSecs, 'G');
    epochIsodateStr = epochIsodateStr.substr(0, 10);
    epochIsodateStr.insert(0, "seconds since ");

    // read HKT file list, loop through HKT files
    ifstream hktFile(hktList);
    string hktFileName;

    double *attitudeTimes = new double[maxNumNavRecords];
    for (size_t i = 0; i < maxNumNavRecords; i++)
        attitudeTimes[i] = BAD_FLT;
    double *orbitTimes = new double[maxNumNavRecords];
    for (size_t i = 0; i < maxNumNavRecords; i++)
        orbitTimes[i] = BAD_FLT;
    double *tiltTimes = new double[maxNumNavRecords];
    for (size_t i = 0; i < maxNumNavRecords; i++)
        tiltTimes[i] = BAD_FLT;

    float *attitudeRates = new float[3 * maxNumNavRecords];
    for (size_t i = 0; i < 3 * maxNumNavRecords; i++)
        attitudeRates[i] = BAD_FLT;
    float *attitudeQuaternions = new float[4 * maxNumNavRecords];
    for (size_t i = 0; i < 4 * maxNumNavRecords; i++)
        attitudeQuaternions[i] = BAD_FLT;
    float *orbitPositions = new float[3 * maxNumNavRecords];
    for (size_t i = 0; i < 3 * maxNumNavRecords; i++)
        orbitPositions[i] = -9999999;
    float *orbitVelocities = new float[3 * maxNumNavRecords];
    for (size_t i = 0; i < maxNumNavRecords; i++)
        orbitVelocities[i] = -9999999;
    for (size_t i = maxNumNavRecords; i < 3 * maxNumNavRecords; i++)
        orbitVelocities[i] = BAD_FLT;
    float *tilts = new float[maxNumNavRecords];
    for (size_t i = 0; i < maxNumNavRecords; i++)
        tilts[i] = BAD_FLT;
    uint8_t *tiltFlags = new uint8_t[maxNumNavRecords];
    for (size_t i = 0; i < maxNumNavRecords; i++)
        tiltFlags[i] = 0;

    NcVar var;

    while (getline(hktFile, hktFileName)) {  // get the HKT file, line by line. Each line is a new file.
        try {
            NcFile hktFile(hktFileName, NcFile::read);

            string timeCoverageStart = "";
            hktFile.getAtt("time_coverage_start").getValues(timeCoverageStart);

            string timeCoverageEnd = "";
            hktFile.getAtt("time_coverage_end").getValues(timeCoverageEnd);

            // extract year, month, day, hour, min and seconds from the time_coverage_start and end
            int16_t hktYear = 0;
            int16_t hktMonth = 0;
            int16_t hktDay = 0;
            int16_t hktHour = 0;
            int16_t hktMin = 0;
            double hktSec = 0;

            // time coverage to unix, then from unix grab the year, month, day, hour, min and sec
            double timeCoverageStartAsUnix = isodate2unix(timeCoverageStart.c_str());
            unix2ymdhms(timeCoverageStartAsUnix, &hktYear, &hktMonth, &hktDay, &hktHour, &hktMin, &hktSec);

            // conver the time data into unix time
            double hktTimeAsSeconds = (hktHour * 3600) + (hktMin * 60) + hktSec;
            double hktStartTimeAsUnix = ymds2unix(hktYear, hktMonth, hktDay, hktTimeAsSeconds);
            double hktFileEpoch = ymds2unix(hktYear, hktMonth, hktDay, 0.0);

            // if hkt start time is outside the end time of the ancillary data, skip it
            if (hktStartTimeAsUnix > endTimeAsUnixSecs + 10) {
                continue;
            }

            // get hkt end time and skip if end time is before the start time of the ancillary data
            // time coverage to unix, then from unix grab the year, month, day, hour, min and sec
            double timeCoverageEndAsUnix = isodate2unix(timeCoverageEnd.c_str());
            unix2ymdhms(timeCoverageEndAsUnix, &hktYear, &hktMonth, &hktDay, &hktHour, &hktMin, &hktSec);

            // conver the time data into unix time
            hktTimeAsSeconds = (hktHour * 3600) + (hktMin * 60) + hktSec;
            double hktEndTimeAsUnix = ymds2unix(hktYear, hktMonth, hktDay, hktTimeAsSeconds);
            hktFileEpoch = ymds2unix(hktYear, hktMonth, hktDay, 0.0);

            // if hkt's end time is before the start time of the ancillary data, skip also
            if (hktEndTimeAsUnix < startTimeAsUnixSecs - 10) {
                continue;
            }

            // ------ navigation_data dimensions used by attitude, orbit and tilt  ------
            int numQuaternionElements = hktFile.getDim("quaternion_elements").getSize();
            int numVectorElements = hktFile.getDim("vector_elements").getSize();

            // ------ navigation_data/attitude  ------
            int attitudeStartIndex = 1e6;
            int attitudeEndIndex = -1;
            size_t numAttitudeRecords = hktFile.getDim("att_records").getSize();  // att_records dimension

            if (numAttitudeRecords > 0) {
                // read var related to attitude
                vector<double> hktAttTime(numAttitudeRecords);
                vector<float> hktAttQuaternions(numAttitudeRecords * numQuaternionElements);
                vector<float> hktAttRates(numAttitudeRecords * numVectorElements);

                hktFile.getGroup("navigation_data").getVar("att_time").getVar(hktAttTime.data());
                findNavigationIndex(timeCoverageStart, startTimeAsUnixSecs, endTimeAsUnixSecs,
                                    hktAttTime.data(), numAttitudeRecords, attitudeStartIndex,
                                    attitudeEndIndex);

                // grab quaternions and rate of attitude if there is data
                int startEndIndexDiff = attitudeEndIndex - attitudeStartIndex;
                if (startEndIndexDiff > 0) {
                    // grab quaternions and rate from the hkt file
                    hktFile.getGroup("navigation_data").getVar("att_quat").getVar(hktAttQuaternions.data());
                    hktFile.getGroup("navigation_data").getVar("att_rate").getVar(hktAttRates.data());

                    // fix the time if l1a and hkt files have different time epoch
                    synchronizeEpochTime(hktAttTime, numAttitudeRecords, hktFileEpoch, epochTimeAsUnixSecs);

                    // concatenate arrays
                    size_t numIndices = attitudeEndIndex - attitudeStartIndex + 1;
                    memcpy(attitudeTimes + attitudeRecCount, hktAttTime.data() + attitudeStartIndex,
                           numIndices * sizeof(double));
                    memcpy(attitudeQuaternions + attitudeRecCount * numQuaternionElements,
                           hktAttQuaternions.data() + attitudeStartIndex * numQuaternionElements,
                           numQuaternionElements * numIndices * sizeof(float));
                    memcpy(attitudeRates + attitudeRecCount * numVectorElements,
                           hktAttRates.data() + attitudeStartIndex * numVectorElements,
                           numVectorElements * numIndices * sizeof(float));
                    attitudeRecCount += numIndices;
                }
            }

            // ------ navigation_data/orbit  ------
            int orbitStartIndex = 1e6;
            int orbitEndIndex = -1;
            size_t numOrbitRecords = hktFile.getDim("orb_records").getSize();  // att_records dimension

            if (numOrbitRecords > 0) {
                // read var related to attitude
                vector<double> hktOrbTime(numOrbitRecords);
                vector<float> hktOrbPosition(numOrbitRecords * numVectorElements);
                vector<float> hktOrbVelocities(numOrbitRecords * numVectorElements);

                hktFile.getGroup("navigation_data").getVar("orb_time").getVar(hktOrbTime.data());
                findNavigationIndex(timeCoverageStart, startTimeAsUnixSecs, endTimeAsUnixSecs,
                                    hktOrbTime.data(), numOrbitRecords, orbitStartIndex, orbitEndIndex);

                // grab quaternions and rate of attitude if there is data
                int startEndIndexDiff = orbitEndIndex - orbitStartIndex;
                if (startEndIndexDiff > 0) {
                    // grab quaternions and rate from the hkt file
                    hktFile.getGroup("navigation_data").getVar("orb_pos").getVar(hktOrbPosition.data());
                    hktFile.getGroup("navigation_data").getVar("orb_vel").getVar(hktOrbVelocities.data());

                    // fix the time if l1a and hkt files have different time epoch
                    synchronizeEpochTime(hktOrbTime, numOrbitRecords, hktFileEpoch, epochTimeAsUnixSecs);

                    // concatenate arrays
                    size_t numIndices = orbitEndIndex - orbitStartIndex + 1;
                    memcpy(orbitTimes + orbitRecCount, hktOrbTime.data() + orbitStartIndex,
                           numIndices * sizeof(double));
                    memcpy(orbitPositions + orbitRecCount * numVectorElements,
                           hktOrbPosition.data() + orbitStartIndex * numVectorElements,
                           numVectorElements * numIndices * sizeof(float));
                    memcpy(orbitVelocities + orbitRecCount * numVectorElements,
                           hktOrbVelocities.data() + orbitStartIndex * numVectorElements,
                           numVectorElements * numIndices * sizeof(float));
                    orbitRecCount += numIndices;
                }
            }

            // ------ navigation_data/tilt  ------
            int tiltStartIndex = 1e6;
            int tiltEndIndex = -1;
            size_t numTiltRecords = hktFile.getDim("tilt_records").getSize();  // att_records dimension

            if (numTiltRecords > 0) {
                // read var related to attitude
                vector<double> hktTiltTime(numTiltRecords);
                vector<float> hktTilt(numTiltRecords);
                vector<float> hktTiltFlags(numTiltRecords);

                hktFile.getGroup("navigation_data").getVar("tilt_time").getVar(hktTiltTime.data());
                findNavigationIndex(timeCoverageStart, startTimeAsUnixSecs, endTimeAsUnixSecs,
                                    hktTiltTime.data(), numTiltRecords, tiltStartIndex, tiltEndIndex);

                // grab quaternions and rate of attitude if there is data
                int startEndIndexDiff = tiltEndIndex - tiltStartIndex;
                if (startEndIndexDiff > 0) {
                    // grab quaternions and rate from the hkt file
                    hktFile.getGroup("navigation_data").getVar("tilt").getVar(hktTilt.data());

                    // fix the time if l1a and hkt files have different time epoch
                    synchronizeEpochTime(hktTiltTime, numTiltRecords, hktFileEpoch, epochTimeAsUnixSecs);

                    // concatenate arrays
                    size_t numIndices = tiltEndIndex - tiltStartIndex + 1;
                    memcpy(tiltTimes + tiltRecCount, hktTiltTime.data() + tiltStartIndex,
                           numIndices * sizeof(double));
                    memcpy(tilts + tiltRecCount, hktTilt.data() + tiltStartIndex, numIndices * sizeof(float));

                    // check to see if there is tilt flag data in the HKT file
                    NcVar tiltFlagVar = hktFile.getGroup("navigation_data").getVar("tilt_flag");
                    if (!tiltFlagVar.isNull()) {
                        // grab and save the flag data and copy it to the outfile tiltFlag array
                        tiltFlagVar.getVar(hktTiltFlags.data());
                        memcpy(tiltFlags + tiltRecCount, hktTiltFlags.data() + tiltStartIndex, numIndices);

                        // set tilt time and tile angle to fill value where the tilt flag is set
                        for (size_t i = tiltRecCount; i < tiltRecCount + numIndices; i++) {
                            if (tiltFlags[i] > 0) {
                                tiltTimes[i] = BAD_FLT;
                                tilts[i] = BAD_FLT;
                            }
                        }
                    }

                    tiltRecCount += numIndices;
                }
            }

            hktFile.close();

        } catch (NcException &e) {
            // e.what();
            cout << "Error reading " << hktFileName << "!" << endl;
            hktFile.close();
            return NC2_ERR;
        }

    }  // while (getline(file, hktFileName))
    hktFile.close();

    // write to L1A file
    vector<size_t> start;
    vector<size_t> count;

    NcGroup l1a_gid = l1afile->getGroup("navigation_data");

    if (attitudeRecCount == 0) {
        attitudeRecCount = 1;
    }

    // 2D dims consts when writing to navigation data
    const size_t QUATERNION_ELEMENTS = 4;
    const size_t VECTOR_ELEMENTS = 3;

    // att_time
    start = {dimShape->attRecords};
    count = {attitudeRecCount};
    var = l1a_gid.getVar("att_time");
    var.putVar(start, count, attitudeTimes);
    var.putAtt("units", epochIsodateStr);

    // att_quat
    start = {dimShape->attRecords, 0};
    count = {attitudeRecCount, QUATERNION_ELEMENTS};
    var = l1a_gid.getVar("att_quat");
    var.putVar(start, count, attitudeQuaternions);

    // att_rate -- same as quat, just 2nd dim is VECTOR_ELEMENTS, keep start the same.
    count = {attitudeRecCount, VECTOR_ELEMENTS};
    var = l1a_gid.getVar("att_rate");
    var.putVar(start, count, attitudeRates);

    /// --- att_records dim is done writing. Update the shape with attitudeRecCount
    dimShape->incrementAttRecordsShape(attitudeRecCount);


    // orb_time
    start = {dimShape->orbRecords};
    count = {orbitRecCount};
    var = l1a_gid.getVar("orb_time");
    var.putVar(start, count, orbitTimes);
    var.putAtt("units", epochIsodateStr);

    // orb_pos
    start = {dimShape->orbRecords, 0};
    count = {orbitRecCount, VECTOR_ELEMENTS};
    var = l1a_gid.getVar("orb_pos");
    var.putVar(start, count, orbitPositions);

    // orb_vel -- // shares same start and count, only 2nd dim is different
    count = {orbitRecCount, VECTOR_ELEMENTS}; 
    var = l1a_gid.getVar("orb_vel");
    var.putVar(start, count, orbitVelocities);

    /// --- update the shape for orb_records dim
    dimShape->incrementOrbRecordsShape(orbitRecCount);


    // tilt_time
    start = {dimShape->tiltSamples};
    count = {tiltRecCount};
    var = l1a_gid.getVar("tilt_time");
    var.putVar(start, count, tiltTimes);
    var.putAtt("units", epochIsodateStr);

    // tilt -- shares the same start and count as tilt_time
    var = l1a_gid.getVar("tilt");
    var.putVar(start, count, tilts);

    /// --- update the shape for tilt_samples dim
    dimShape->incrementTiltSampleShape(tiltRecCount);

    delete[] attitudeTimes;
    delete[] orbitTimes;
    delete[] tiltTimes;
    delete[] attitudeRates;
    delete[] attitudeQuaternions;
    delete[] orbitPositions;
    delete[] orbitVelocities;
    delete[] tilts;
    delete[] tiltFlags;

    return 0;
}

// closes the l1afile when it is done writing
int L1aFile::close() {
    try {
        this->l1afile->close();
    } catch (NcException &e) {
        cout << e.what() << endl;
        cerr << "Failure closing: " + fileName << endl;
        exit(1);
    }

    return 0;
}