#include <stdio.h>
#include <math.h>
#include <libgen.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <netcdf>

#include <cgal_interp.h>
#include <allocate2d.h>
#include <allocate3d.h>
#include <timeutils.h>
#include <genutils.h>
#include <clo.h>
#include <genutils.h>

#include <gsl/gsl_matrix_float.h>
#include <gsl/gsl_blas.h>

#include "copy_var_utils.hpp"

#define NCACHE 20
const double PI = OEL_PI;
const double RADEG = OEL_RADEG;

using namespace std;
using namespace netCDF;
using namespace netCDF::exceptions;

//    Modification history:
//  Programmer     Organization   Date     Ver   Description of change
//  ----------     ------------   ----     ---   ---------------------
//  Joel Gales     SAIC           01/14/24 1.50  Fix type bug in the
//                                               calculation of the water cloud
//                                               average values for the L1C bins
//
//  Joel Gales     SAIC           02/01/24 1.51  Make CLD processing optional
//  Martin Montes    SSAI           07/22/24 1.52  Fix water/ice cloud mask

/**
 * @brief accum calculates the arithmetic average of ice/water cloud fraction for all L1C grid cells
 * @param totLines number of lines for target file
 * @param nPixels number of pixels for target file
 * @param naLong number of along-track L1C grid cells
 * @param naCross number of across-track L1C grid cells
 * @param brow L1C grid row indexes for resulting lines x pixels
 * @param bcol L1C grid col indexes for resulting lines x pixels
 * @param cloudProd generic input cloud product 2-D array
 * @param cloudPhase generic input cloud phase product 2-D array
 * @param binvalIceCloud binned value for ice cloud mask
 * @param nobsIceCloud number of observations per bin for ice cloud mask
 * @param binvalWaterCloud binned value for water cloud mask
 * @param nobsWaterCloud number of observations per bin for water cloud mask
 * @return the output is the number of observations and values per bin
 */
int accum(uint32_t totLines, uint32_t nPixels, uint32_t naLong, uint32_t naCross, short *brow, short *bcol,
          float **cloudProd, int8_t **cloudPhase, float **binvalIceCloud, short **nobsIceCloud,
          float **binvalWaterCloud, short **nobsWaterCloud) {
    // Clear accumulation arrays
    for (size_t i = 0; i < naLong; i++) {
        for (size_t j = 0; j < naCross; j++) {
            binvalIceCloud[i][j] = 0.0;
            binvalWaterCloud[i][j] = 0.0;
            nobsIceCloud[i][j] = 0;
            nobsWaterCloud[i][j] = 0;
        }
    }

    // Accumulate bin values and counts
    short *brptr = brow;
    short *bcptr = bcol;
    for (size_t i = 0; i < totLines; i++) {
        for (size_t j = 0; j < nPixels; j++) {
            if (*brptr != -1) {
                if (cloudProd[i][j] == -32767) {
                    brptr++;
                    bcptr++;
                    continue;
                }

                if (cloudPhase[i][j] == 1) {
                    binvalIceCloud[*brptr][*bcptr] += cloudProd[i][j];
                    nobsIceCloud[*brptr][*bcptr]++;
                } else {
                    binvalWaterCloud[*brptr][*bcptr] += cloudProd[i][j];
                    nobsWaterCloud[*brptr][*bcptr]++;
                }
            }
            brptr++;
            bcptr++;
        }
    }
    return 0;
}

// accumulate the ice and water cloud fraction
int accumFrac(uint32_t totLines, uint32_t nPixels, uint32_t naLong, uint32_t naCross, short *brow,
              short *bcol, float **cloudProd, int8_t **cloudPhase, float **binvalIceCloud,
              short **nobsIceCloud, float **binvalWaterCloud, short **nobsWaterCloud) {
    // Clear accumulation arrays
    for (size_t i = 0; i < naLong; i++) {
        for (size_t j = 0; j < naCross; j++) {
            binvalIceCloud[i][j] = 0.0;
            binvalWaterCloud[i][j] = 0.0;
            nobsIceCloud[i][j] = 0;
            nobsWaterCloud[i][j] = 0;
        }
    }

    // Accumulate bin values and counts
    short *brptr = brow;
    short *bcptr = bcol;
    for (size_t i = 0; i < totLines; i++) {
        for (size_t j = 0; j < nPixels; j++) {
            if (*brptr != -1) {
                if (cloudProd[i][j] != -32767) {
                    nobsIceCloud[*brptr][*bcptr]++;
                    nobsWaterCloud[*brptr][*bcptr]++;

                    if (cloudProd[i][j] == 1) {
                        if (cloudPhase[i][j] == 1) {
                            binvalIceCloud[*brptr][*bcptr]++;
                        } else {
                            binvalWaterCloud[*brptr][*bcptr]++;
                        }
                    }
                }
            }
            brptr++;
            bcptr++;
        }
    }

    return 0;
}

/**
 * @brief accumWm calculates the arithmetic average of watermask fraction for all L1C grid cells
 * @param nWaterMask total number of water mask cells (rowxcol)
 * @param wLval discrete watermask value
 * @param binVal binned watermask value
 * @param nobs number of observations per bin
 * @return same as accum this function do sum of values and counting but for watermask fraction
 */
int accumWm(uint32_t nWaterMask, short *brow, short *bcol, float *wLval, float **binVal, short **nobs) {
    short *brptr = brow;
    short *bcptr = bcol;
    for (size_t i = 0; i < nWaterMask; i++) {
        if (*brptr != -1) {
            if (wLval[i] == -32767) {
                brptr++;
                bcptr++;
                continue;
            }

            binVal[*brptr][*bcptr] += wLval[i];
            nobs[*brptr][*bcptr]++;
        }
        brptr++;
        bcptr++;
    }

    return 0;
}

/**
 * @brief search row/col indexes of L1C grid for all values ncm in target file
 * @param naLong number of along-track lines for L1C grid
 * @param naCross number of across-track lines for L1C grid
 * @param ncm (lines x pixels) for target file
 * @param gridRes L1C grid spatial resolution in km, nominal: 5.2 km
 * @param lonL1C longitude in degrees for L1C grid (naLong x naCross)
 * @param latL1C latitude in degrees for L1C grid (naLong x naCross)
 * @param lonCm longitude in degrees for target file (lines x pixels)
 * @param latCm  latitude in degrees for target file (lines x pixels)
 * @return output found row/col for all values of target file
 */
int lonlat2rowcol(uint32_t naLong, uint32_t naCross, uint32_t ncm, float gridRes, float *lonL1C,
                  float *latL1C, float *lonCm, float *latCm, short *brow, short *bcol);


/**
 * @brief read different types of cloud products
 * @tparam T
 * @param ncGrp netcdf4 group
 * @param fieldName field name
 * @param array generic target array pointing to template
 * @param nPixels
 * @param nLines
 * @return output an array containing the field name data
 */
template <typename T>
int readCLD(NcGroup ncGrp[], const char *fieldName, T *array, uint32_t nPixels, uint32_t nLines[]) {
    vector<size_t> start, count;
    NcVar var;

    uint32_t nLines0 = 0;

    // Read from trailing granule if specified
    if (!ncGrp[0].isNull()) {
        if (nLines[0] < 250) {
            nLines0 = nLines[0];
        } else {
            nLines0 = 250;
        }

        start.clear();
        start.push_back(nLines[0] - nLines0);
        start.push_back(0);

        count.clear();
        count.push_back(nLines0);
        count.push_back(nPixels);

        var = ncGrp[0].getVar(fieldName);
        var.getVar(start, count, &array[0]);
    }

    var = ncGrp[1].getVar(fieldName);
    var.getVar(&array[nLines0 * nPixels]);

    // Read from following granule if specified
    if (!ncGrp[2].isNull()) {
        start.clear();
        start.push_back(0);
        start.push_back(0);

        count.clear();
        if (nLines[2] < 250)
            count.push_back(nLines[2]);
        else
            count.push_back(250);
        count.push_back(nPixels);

        var = ncGrp[2].getVar(fieldName);
        var.getVar(start, count, &array[(nLines0 + nLines[1]) * nPixels]);
    }

    return 0;
}

/**
 * @brief add full path of env variable to target file after
 * @param sValue
 * @return output target file name with expanded env variable
 */
inline int expandEnvVar(std::string *sValue) {
    if ((*sValue).find_first_of("$") == std::string::npos)
        return 0;
    std::string::size_type posEndIdx = (*sValue).find_first_of("/");
    if (posEndIdx == std::string::npos)
        return 0;
    const std::string envVar = sValue->substr(1, posEndIdx - 1);
    char *envVarStr = getenv(envVar.c_str());
    if (envVarStr == 0x0) {
        printf("Environment variable: %s not defined.\n", sValue->c_str());
        exit(1);
    }
    *sValue = envVarStr + (*sValue).substr(posEndIdx);

    return 0;
}

/**
 * @brief
 * @param argc
 * @param argv
 * @return
 */
int main(int argc, char *argv[]) {
    NcVar var, varin, varout;

    clo_optionList_t *list;
    clo_option_t *option;
    char *strVal;
    char keyword[50];
    char parFilename[FILENAME_MAX];
    char targetFilename[FILENAME_MAX];
    char chlorFilename[FILENAME_MAX];
    char profileName[FILENAME_MAX];
    char metFilename[FILENAME_MAX];
    char aerFilename[FILENAME_MAX];
    char geoscfFilename[FILENAME_MAX];
    char cloudMask1[FILENAME_MAX];
    char cloudMask2[FILENAME_MAX];
    char cloudMask3[FILENAME_MAX];
    char cloudProd1[FILENAME_MAX];
    char cloudProd2[FILENAME_MAX];
    char cloudProd3[FILENAME_MAX];
    char albedoFilename[FILENAME_MAX];
    char camsch4Filename[FILENAME_MAX];
    char camsco2Filename[FILENAME_MAX];
    char camsn2oFilename[FILENAME_MAX];
    char gebcoFilename[FILENAME_MAX];
    char ancoutFilename[FILENAME_MAX];

    list = clo_createList();

    option = clo_addOption(list, "targetfile", CLO_TYPE_IFILE, NULL, "Input target file");
    option = clo_addOption(list, "ancfile", CLO_TYPE_OFILE, NULL, "Output ancillary file");
    option = clo_addOption(list, "merraprofile", CLO_TYPE_IFILE, NULL, "Input MERRA2 profile ancillary file");
    option = clo_addOption(list, "merrametfile", CLO_TYPE_IFILE, NULL, "Input MERRA2 MET ancillary file");
    option = clo_addOption(list, "merraaerfile", CLO_TYPE_IFILE, NULL,
                           "Input MERRA2 AER ancillary file (optional)");
    option = clo_addOption(list, "geoscffile", CLO_TYPE_IFILE, NULL, "Input GEOS CF ancillary file");

    option = clo_addOption(list, "cldmask1", CLO_TYPE_IFILE, NULL, "Input trailing cloud mask L2 file (optional)");
    option = clo_addOption(list, "cldmask2", CLO_TYPE_IFILE, NULL, "Input current cloud mask L2 file (optional)");
    option = clo_addOption(list, "cldmask3", CLO_TYPE_IFILE, NULL, "Input following cloud mask L2 file (optional)");

    option = clo_addOption(list, "cldprod1", CLO_TYPE_IFILE, NULL, "Input trailing cloud product L2 file (optional)");
    option = clo_addOption(list, "cldprod2", CLO_TYPE_IFILE, NULL, "Input current cloud product L2 file (optional)");
    option = clo_addOption(list, "cldprod3", CLO_TYPE_IFILE, NULL, "Input following cloud product L2 file (optional)");

    option = clo_addOption(list, "albedofile", CLO_TYPE_IFILE, NULL, "Input albedo ancillary file");
    option = clo_addOption(list, "chlfile", CLO_TYPE_IFILE, NULL, "Input chlor_a L3 map file (optional)");
    option = clo_addOption(list, "ch4file", CLO_TYPE_IFILE, NULL, "Input CAMS CH4 ancillary file");
    option = clo_addOption(list, "co2file", CLO_TYPE_IFILE, NULL, "Input CAMS CO2 ancillary file");
    option = clo_addOption(list, "n2ofile", CLO_TYPE_IFILE, NULL, "Input CAMS N2O ancillary file");

    option = clo_addOption(list, "gebcofile", CLO_TYPE_IFILE, NULL, "Input GEBCO landmask file (optional)");

    clo_setVersion("1.52");

    if (argc == 1) {
        clo_printUsage(list);
        exit(1);
    }

    clo_readArgs(list, argc, argv);

    // below are all REQUIRED parameters for ancgen to run. So if the user does not provide it, it exit(1).

    strVal = clo_getString(list, "targetfile");
    strcpy(targetFilename, strVal);

    strVal = clo_getString(list, "ancfile");
    strcpy(ancoutFilename, strVal);

    strVal = clo_getString(list, "merraprofile");
    strcpy(profileName, strVal);

    strVal = clo_getString(list, "merrametfile");
    strcpy(metFilename, strVal);

    strVal = clo_getString(list, "geoscffile");
    strcpy(geoscfFilename, strVal);

    strVal = clo_getString(list, "albedofile");
    strcpy(albedoFilename, strVal);

    strVal = clo_getString(list, "ch4file");
    strcpy(camsch4Filename, strVal);

    strVal = clo_getString(list, "co2file");
    strcpy(camsco2Filename, strVal);

    strVal = clo_getString(list, "n2ofile");
    strcpy(camsn2oFilename, strVal);


    // options here are all optional, so search for the keyword and if it doesnt exist, set empty string

    int numOptions, optionId;

    numOptions = clo_getNumOptions(list);
    for (optionId = 0; optionId < numOptions; optionId++) {
        option = clo_getOption(list, optionId);

        strcpy(keyword, option->key);

        if (strcmp(keyword, "chlfile") == 0) {
            if (clo_isOptionSet(option))
                parse_file_name(clo_getOptionString(option), chlorFilename);
            else
                strcpy(chlorFilename, "");

        } else if (strcmp(keyword, "merraaerfile") == 0) {
            if (clo_isOptionSet(option))
                parse_file_name(clo_getOptionString(option), aerFilename);
            else
                strcpy(aerFilename, "");

        } else if (strcmp(keyword, "cldmask1") == 0) {
            if (clo_isOptionSet(option))
                parse_file_name(clo_getOptionString(option), cloudMask1);
            else
                strcpy(cloudMask1, "");

        } else if (strcmp(keyword, "cldmask2") == 0) {
            if (clo_isOptionSet(option))
                parse_file_name(clo_getOptionString(option), cloudMask2);
            else
                strcpy(cloudMask2, "");

        } else if (strcmp(keyword, "cldmask3") == 0) {
            if (clo_isOptionSet(option))
                parse_file_name(clo_getOptionString(option), cloudMask3);
            else
                strcpy(cloudMask3, "");

        } else if (strcmp(keyword, "cldprod1") == 0) {
            if (clo_isOptionSet(option))
                parse_file_name(clo_getOptionString(option), cloudProd1);
            else
                strcpy(cloudProd1, "");

        } else if (strcmp(keyword, "cldprod2") == 0) {
            if (clo_isOptionSet(option))
                parse_file_name(clo_getOptionString(option), cloudProd2);
            else
                strcpy(cloudProd2, "");

        } else if (strcmp(keyword, "cldprod3") == 0) {
            if (clo_isOptionSet(option))
                parse_file_name(clo_getOptionString(option), cloudProd3);
            else
                strcpy(cloudProd3, "");
        } else if (strcmp(keyword, "gebcofile") == 0) {
            if (clo_isOptionSet(option))
                parse_file_name(clo_getOptionString(option), gebcoFilename);
            else
                strcpy(gebcoFilename, "$OCDATAROOT/common/gebco_ocssw_v2020.nc");
        }
        else if (strcmp(keyword, "par") == 0) {
            if (clo_isOptionSet(option)) 
                parse_file_name(clo_getOptionString(option), parFilename);
            else
                strcpy(parFilename, "");
        }
    }
    ////////////////////////////////
    ///////// Open L1C file ////////
    ////////////////////////////////
    // Read lon/lat of l1c
    cout << endl << "Opening L1C file" << endl;
    NcFile *L1Cfile = new NcFile(targetFilename, NcFile::read);
    // get the start time
    NcGroupAtt att = L1Cfile->getAtt("time_coverage_start");
    if (att.isNull()) {
        cout << "Error - Could not find time_coverage_start global attribute.\n";
        exit(1);
    }
    string timeCoverageStart;
    att.getValues(timeCoverageStart);

    // get the end time
    att = L1Cfile->getAtt("time_coverage_end");
    if (att.isNull()) {
        cout << "Error - Could not find time_coverage_end global attribute.\n";
        exit(1);
    }
    string timeCoverageEnd;
    att.getValues(timeCoverageEnd);

    // get instrument
    att = L1Cfile->getAtt("instrument");
    if (att.isNull()) {
        cout << "Error - Could not find instrument global attribute.\n";
        exit(1);
    }
    string instrument;
    att.getValues(instrument);

    NcDim acrossDim = L1Cfile->getDim("bins_across_track");
    uint32_t naCross = acrossDim.getSize();
    NcDim alongDim = L1Cfile->getDim("bins_along_track");
    uint32_t naLong = alongDim.getSize();

    float **outL1c2d = allocate2d_float(naLong, naCross);
    float **latL1C = allocate2d_float(naLong, naCross);
    float **lonL1C = allocate2d_float(naLong, naCross);

    int binsAlongTrack = naLong;
    int binsAcrossTrack = naCross;

    size_t nl1c = binsAlongTrack * binsAcrossTrack;

    vector<NcDim> dims;

    float gridRes = 5.2;
    ////////////////////////////////
    /////// Create ANC file ////////
    ////////////////////////////////
    // Create output ancillary file
    cout << "Creating ancillary file" << endl;
    NcFile *ncOutput;
    ncOutput = new NcFile(ancoutFilename, NcFile::replace);

    // Global metadata
    string dateCreated = string(unix2isodate(now(), 'G'));

    ncOutput->putAtt("date_created", dateCreated);
    ncOutput->putAtt("time_coverage_start", timeCoverageStart);
    ncOutput->putAtt("time_coverage_end", timeCoverageEnd);
    ncOutput->putAtt("title", instrument + " L1C ancillary file");
    ncOutput->putAtt("product_name", basename(ancoutFilename));
    ncOutput->putAtt("creator_name", "NASA/GSFC/OBPG");
    ncOutput->putAtt("creator_url", "http://oceandata.sci.gsfc.nasa.gov");
    ncOutput->putAtt("creator_email", "data@oceancolor.gsfc.nasa.gov");
    ncOutput->putAtt("project", "Ocean Biology Processing Group (NASA/GSFC/OBPG)");
    ncOutput->putAtt("publisher_name", "NASA/GSFC/OBPG");
    ncOutput->putAtt("publisher_url", "http://oceandata.sci.gsfc.nasa.gov");
    ncOutput->putAtt("publisher_email", "data@oceancolor.gsfc.nasa.gov");

    // processing control group
    NcGroup processingControlGrp = ncOutput->addGroup("processing_control");
    processingControlGrp.putAtt("software_name", "ancgen");
    processingControlGrp.putAtt("software_version", clo_getVersion());

    // input parameter group - a subgroup of processing control
    NcGroup inputParametersGrp = processingControlGrp.addGroup("input_parameter");
    inputParametersGrp.putAtt("par", parFilename);
    inputParametersGrp.putAtt("targetfile", targetFilename);
    inputParametersGrp.putAtt("ancfile", ancoutFilename);
    inputParametersGrp.putAtt("merraprofile", profileName);
    inputParametersGrp.putAtt("merrametfile", metFilename);
    inputParametersGrp.putAtt("merraaerfile",aerFilename);
    inputParametersGrp.putAtt("geoscffile", geoscfFilename);
    inputParametersGrp.putAtt("cldmask1", cloudMask1);
    inputParametersGrp.putAtt("cldmask2", cloudMask2);
    inputParametersGrp.putAtt("cldmask3", cloudMask3);
    inputParametersGrp.putAtt("cldprod1", cloudProd1);
    inputParametersGrp.putAtt("cldprod2", cloudProd2);
    inputParametersGrp.putAtt("cldprod3", cloudProd3);
    string tmpStr = albedoFilename;
    replaceOCroots(tmpStr);
    inputParametersGrp.putAtt("albedofile", tmpStr);
    inputParametersGrp.putAtt("chlfile", chlorFilename);
    tmpStr = camsch4Filename;
    replaceOCroots(tmpStr);
    inputParametersGrp.putAtt("ch4file", tmpStr);
    tmpStr = camsco2Filename;
    replaceOCroots(tmpStr);
    inputParametersGrp.putAtt("co2file", tmpStr);
    tmpStr = camsn2oFilename;
    replaceOCroots(tmpStr);
    inputParametersGrp.putAtt("n2ofile", tmpStr);
    tmpStr = gebcoFilename;
    replaceOCroots(tmpStr);
    inputParametersGrp.putAtt("gebcofile", tmpStr);


    // Create netCDF dimensions
    NcDim rows = ncOutput->addDim("bins_along_track", binsAlongTrack);
    NcDim cols = ncOutput->addDim("bins_across_track", binsAcrossTrack);

    // Add target lat/lon fields
    dims.push_back(rows);
    dims.push_back(cols);

    // Get L1C lon/lat
    NcGroup gidGEO = L1Cfile->getGroup("geolocation_data");

    float *fptr;

    var = gidGEO.getVar("latitude");
    var.getVar(&latL1C[0][0]);
    varout = ncOutput->addVar("latitude", ncFloat, dims);

    string overridell[] = {"-32767", "=", "=", "=", "="};
    string overrideTypell[] = {"F", "=", "=", "=", "="};
    copyVarAtts(&var, &varout, overridell, overrideTypell);

    // Convert fill values to -32767
    fptr = &latL1C[0][0];
    for (size_t i = 0; i < nl1c; i++) {
        if (*fptr < -900)
            *fptr = -32767;
        fptr++;
    }
    varout.putVar(&latL1C[0][0]);

    var = gidGEO.getVar("longitude");
    var.getVar(&lonL1C[0][0]);
    varout = ncOutput->addVar("longitude", ncFloat, dims);

    copyVarAtts(&var, &varout, overridell, overrideTypell);

    // Convert fill values to -32767
    fptr = &lonL1C[0][0];
    for (size_t i = 0; i < nl1c; i++) {
        if (*fptr < -900)
            *fptr = -32767;
        fptr++;
    }
    varout.putVar(&lonL1C[0][0]);

    // Get lon/lat min/max for L1C granule
    float lonMax = -180, lonMin = +180;
    fptr = &lonL1C[0][0];
    for (size_t i = 0; i < (size_t)nl1c; i++) {
        if (*fptr < lonMin)
            lonMin = *fptr;
        if (*fptr > lonMax)
            lonMax = *fptr;
        fptr++;
    }

    float latMax = -90, latMin = +90;
    fptr = &latL1C[0][0];
    for (size_t i = 0; i < (size_t)nl1c; i++) {
        if (*fptr < latMin)
            latMin = *fptr;
        if (*fptr > latMax)
            latMax = *fptr;
        fptr++;
    }

    vector<size_t> start, count;
    vector<ptrdiff_t> stride;
    ////////////////////////////////
    //////////// ALBEDO ////////////
    ////////////////////////////////
    size_t albedoNpix = 43200;
    size_t albedoNlines = 21600;
    size_t albedoNfields = 5;
    float albedoPixelSize = 0.00833333;

    cout << endl << "Opening ALBEDO file" << endl;
    NcFile *albedoFile = new NcFile(albedoFilename, NcFile::read);

    uint8_t **albedoRow = new uint8_t *[albedoNlines];

    std::vector<string> albedoFields = {"Albedo_Map_0.659", "Albedo_Map_0.858", "Albedo_Map_1.24",
                                        "Albedo_Map_1.64", "Albedo_Map_2.13"};

    int32_t syear;
    int32_t sday;
    int32_t smsec;
    isodate2ydmsec((char *)timeCoverageStart.c_str(), &syear, &sday, &smsec);
    int32_t period8day = sday / 8;
    start.clear();
    start.push_back(period8day);
    start.push_back(0);
    start.push_back(0);

    count.clear();
    count.push_back(1);
    count.push_back(NCACHE);
    count.push_back(albedoNpix);

    // Loop through albedo fields
    for (size_t k = 0; k < albedoNfields; k++) {
        var = albedoFile->getVar(albedoFields[k].c_str());

        for (size_t i = 0; i < albedoNlines; i++)
            albedoRow[i] = NULL;

        int first = -1;
        for (size_t i = 0; i < naLong; i++) {
            for (size_t j = 0; j < naCross; j++) {
                int latRow = (latL1C[i][j] + 90) / albedoPixelSize;

                if (albedoRow[latRow] == NULL) {
                    // align the row to the cache size boundary
                    size_t cacheBoundary = ((size_t)(latRow / NCACHE)) * NCACHE;
                    if (first == -1)
                        first = cacheBoundary;
                    albedoRow[cacheBoundary] = new uint8_t[albedoNpix * NCACHE];
                    for (size_t k = 1; k < NCACHE; k++)
                        albedoRow[cacheBoundary + k] = albedoRow[cacheBoundary] + k * albedoNpix;

                    start[1] = cacheBoundary;
                    var.getVar(start, count, albedoRow[cacheBoundary]);
                }
                int lonCol = (lonL1C[i][j] + 180) / albedoPixelSize;
                uint8_t albedo = albedoRow[latRow][lonCol];
                if (albedo == 255)
                    outL1c2d[i][j] = -32767;
                else
                    outL1c2d[i][j] = albedo * 0.004;
            }
        }

        cout << "Writing " << albedoFields[k].c_str() << endl;
        varout = ncOutput->addVar(albedoFields[k].c_str(), ncFloat, dims);

        string override[] = {"-32767", "", "=", "-32767", "", "=", "="};
        string overridetype[] = {"F", "", "=", "F", "", "=", "="};
        copyVarAtts(&var, &varout, override, overridetype);

        varout.putVar(&outL1c2d[0][0]);

        for (size_t i = first; i < albedoNlines; i += NCACHE)
            if (albedoRow[i] != NULL)
                delete[] albedoRow[i];
    }
    // albedo field loop
    delete[] albedoRow;
    // SW to NE (landmask)

    // Determine if L1C granule crosses the dateline
    // Adjust L1C lon to 0-360 if dateline crossed

    bool dateLineCross = false;
    if (lonMax - lonMin > 180) {
        cout << endl << "Correcting longitude for dataline crossing" << endl;
        dateLineCross = true;
        lonMax = 0;
        lonMin = 360;
        fptr = &lonL1C[0][0];
        for (size_t i = 0; i < (size_t)nl1c; i++) {
            if (*fptr < 0)
                *fptr += 360;
            if (*fptr < lonMin)
                lonMin = *fptr;
            if (*fptr > lonMax)
                lonMax = *fptr;
            fptr++;
        }
    }

    ////////////////////////////////
    //////////// MERRA2 ////////////
    ////////////////////////////////
    // Read MERRA2 ancillary file
    int nlon = 576;
    int nlat = 361;

    cout << endl << "Opening MERRA2 file" << endl;
    NcFile *merraFile = new NcFile(profileName, NcFile::read);

    float **latMerra = allocate2d_float(361, 576);
    float **lonMerra = allocate2d_float(361, 576);

    for (size_t i = 0; i < 361; i++) {
        for (size_t j = 0; j < 576; j++) {
            latMerra[i][j] = i * 0.50 - 90;
            lonMerra[i][j] = j * 0.625 - 180;

            if (dateLineCross && lonMerra[i][j] < 0)
                lonMerra[i][j] += 360;
        }
    }
    int nmerra = nlon * nlat;

    // Get nearest neighbor lon/lat bin numbers for MERRA
    short **latbinL1C = allocate2d_short(naLong, naCross);
    short **lonbinL1C = allocate2d_short(naLong, naCross);

    fptr = &latL1C[0][0];
    short *slatptr = &latbinL1C[0][0];
    for (size_t i = 0; i < (size_t)nl1c; i++) {
        *slatptr = (short)((*fptr + 90) / 0.50);
        fptr++;
        slatptr++;
    }

    fptr = &lonL1C[0][0];
    short *slonptr = &lonbinL1C[0][0];
    for (size_t i = 0; i < (size_t)nl1c; i++) {
        *slonptr = (short)((*fptr + 180) / 0.625);
        fptr++;
        slonptr++;
    }

    cout << "Computing MERRA to L1C interpolation mapping" << endl;
    int nncTagMerra = cgal_nnc(nmerra, &latMerra[0][0], &lonMerra[0][0], nl1c, &latL1C[0][0], &lonL1C[0][0]);

    size_t nlev = 42;

    NcDim levs = ncOutput->addDim("levels", nlev);

    dims.clear();
    dims.push_back(levs);
    dims.push_back(rows);
    dims.push_back(cols);

    float ***outL1c3d = allocate3d_float(nlev, naLong, naCross);
    float ***merraData = allocate3d_float(42, 361, 576);

    // Geo-potential height profile
    // Temperature vertical profile
    // Relative humidity profile
    // Specific humidity profile

    string merraFields[5] = {"H", "T", "RH", "QV", "O3"};

    // Product loop
    for (size_t k = 0; k < 5; k++) {
        varin = merraFile->getVar(merraFields[k].c_str());
        varin.getVar(&merraData[0][0][0]);

        // Convert MERRA2 fill values to -1e14
        float *fptr = &merraData[0][0][0];
        for (size_t i = 0; i < 42 * 361 * 576; i++) {
            if (*fptr > 1e14)
                *fptr = -1e14;
            fptr++;
        }

        // Level loop
        for (size_t i = 0; i < nlev; i++) {
            cgal_interp2(nmerra, &latMerra[0][0], &lonMerra[0][0], &merraData[i][0][0], nl1c,
                         &outL1c3d[i][0][0], nncTagMerra);

            // Convert interpolation failures to -32767
            fptr = &outL1c3d[0][0][0];
            for (size_t i = 0; i < nlev * naLong * naCross; i++) {
                if (*fptr < 0)
                    *fptr = -32767;
                fptr++;
            }

            // Attempt to correct fill values with nearest-neighbor values
            fptr = &outL1c3d[i][0][0];
            slatptr = &latbinL1C[0][0];
            slonptr = &lonbinL1C[0][0];

            for (size_t j = 0; j < (size_t)nl1c; j++) {
                if (*fptr == -32767)
                    *fptr = merraData[i][*slatptr][*slonptr];

                fptr++;
                slatptr++;
                slonptr++;
            }
        }

        cout << "Writing " << merraFields[k] << endl;

        varout = ncOutput->addVar(merraFields[k].c_str(), ncFloat, dims);

        string override[] = {"-32767", "=", "-32767", "=", "="};
        string overrideType[] = {"F", "=", "F", "=", "="};
        copyVarAtts(&varin, &varout, override, overrideType);

        varout.putVar(&outL1c3d[0][0][0]);
    }

    delete merraFile;

    ////////////////////////////////
    ////////// MERRA2 MET //////////
    ////////////////////////////////
    cout << endl << "Opening MERRA2 MET file" << endl;
    merraFile = new NcFile(metFilename, NcFile::read);

    dims.clear();
    dims.push_back(rows);
    dims.push_back(cols);

    string merraMetFields[10] = {"PS",  "QV10M", "SLP",  "T10M",  "TO3",
                                 "TQV", "U10M",  "V10M", "FRSNO", "FRSEAICE"};

    for (size_t k = 0; k < 10; k++) {
        varin = merraFile->getVar(merraMetFields[k].c_str());
        varin.getVar(&merraData[0][0][0]);

        // Convert MERRA2 fill values to -1e14
        float *fptr = &merraData[0][0][0];
        for (size_t i = 0; i < 361 * 576; i++) {
            if (*fptr > 1e14)
                *fptr = -1e14;
            fptr++;
        }

        cgal_interp2(nmerra, &latMerra[0][0], &lonMerra[0][0], &merraData[0][0][0], nl1c, &outL1c3d[0][0][0],
                     nncTagMerra);

        // Convert interpolation failures to -32767
        fptr = &outL1c3d[0][0][0];
        for (size_t i = 0; i < naLong * naCross; i++) {
            if (*fptr < 0)
                *fptr = -32767;
            fptr++;
        }

        // Attempt to correct fill values with nearest-neighbor values
        fptr = &outL1c3d[0][0][0];
        slatptr = &latbinL1C[0][0];
        slonptr = &lonbinL1C[0][0];

        for (size_t j = 0; j < (size_t)nl1c; j++) {
            if (*fptr == -32767)
                *fptr = merraData[0][*slatptr][*slonptr];

            fptr++;
            slatptr++;
            slonptr++;
        }

        cout << "Writing " << merraMetFields[k] << endl;

        varout = ncOutput->addVar(merraMetFields[k].c_str(), ncFloat, dims);

        string override[] = {"-32767", "=", "-32767", "=", "="};
        string overrideType[] = {"F", "=", "F", "=", "="};
        copyVarAtts(&varin, &varout, override, overrideType);
        varout.putVar(&outL1c3d[0][0][0]);
    }

    ////////////////////////////////
    ////////// MERRA2 AER //////////
    ////////////////////////////////
    if (strcmp(aerFilename, "") != 0) {
        cout << endl << "Opening MERRA2 AER_file" << endl;
        merraFile = new NcFile(aerFilename, NcFile::read);

        dims.clear();
        dims.push_back(rows);
        dims.push_back(cols);

        string merraAerFields[13] = {"BCEXTTAU",  "BCSCATAU",  "DUEXTTAU", "DUSCATAU", "SSEXTTAU",
                                     "SSSCATAU",  "SUEXTTAU",  "SUSCATAU", "OCEXTTAU", "OCSCATAU",
                                     "TOTEXTTAU", "TOTSCATAU", "TOTANGSTR"};

        for (size_t k = 0; k < 13; k++) {
            varin = merraFile->getVar(merraAerFields[k].c_str());
            varin.getVar(&merraData[0][0][0]);

            // Convert MERRA2 fill values to -1e14
            float *fptr = &merraData[0][0][0];
            for (size_t i = 0; i < 361 * 576; i++) {
                if (*fptr > 1e14)
                    *fptr = -1e14;
                fptr++;
            }

            cgal_interp2(nmerra, &latMerra[0][0], &lonMerra[0][0], &merraData[0][0][0], nl1c,
                         &outL1c3d[0][0][0], nncTagMerra);

            // Convert interpolation failures to -32767
            fptr = &outL1c3d[0][0][0];
            for (size_t i = 0; i < naLong * naCross; i++) {
                if (*fptr < 0)
                    *fptr = -32767;
                fptr++;
            }

            // Attempt to correct fill values with nearest-neighbor values
            fptr = &outL1c3d[0][0][0];
            slatptr = &latbinL1C[0][0];
            slonptr = &lonbinL1C[0][0];

            for (size_t j = 0; j < (size_t)nl1c; j++) {
                if (*fptr == -32767)
                    *fptr = merraData[0][*slatptr][*slonptr];

                fptr++;
                slatptr++;
                slonptr++;
            }

            cout << "Writing " << merraAerFields[k] << endl;

            varout = ncOutput->addVar(merraAerFields[k].c_str(), ncFloat, dims);

            string override[] = {"-32767", "=", "-32767", "=", "="};
            string overrideType[] = {"F", "=", "F", "=", "="};
            copyVarAtts(&varin, &varout, override, overrideType);
            varout.putVar(&outL1c3d[0][0][0]);
        }

        free2d_short(lonbinL1C);
        free2d_short(latbinL1C);
        free3d_float(merraData);
    }
    /////////////////////////////////
    //////////// CAMS SAT ///////////
    /////////////////////////////////
    NcFile *camsFile = new NcFile(camsch4Filename, NcFile::read);

    NcDim camsLatDim = camsFile->getDim("latitude_bins");
    uint32_t camsLat = camsLatDim.getSize();

    NcDim camsLonDim = camsFile->getDim("longitude_bins");
    uint32_t camsLon = camsLonDim.getSize();

    float **latCams = allocate2d_float(camsLat, camsLon);
    float **lonCams = allocate2d_float(camsLat, camsLon);

    for (size_t i = 0; i < camsLat; i++) {
        for (size_t j = 0; j < camsLon; j++) {
            latCams[i][j] = i * 2 - 89.0;
            lonCams[i][j] = j * 3 - 178.5;

            if (dateLineCross && lonCams[i][j] < 0)
                lonCams[i][j] += 360;
        }
    }
    int ncams = camsLon * camsLat;

    cout << endl << "Computing CAMS SATELLITE to L1C interpolation mapping" << endl;
    int nncTagCams = cgal_nnc(ncams, &latCams[0][0], &lonCams[0][0], nl1c, &latL1C[0][0], &lonL1C[0][0]);

    dims.clear();
    dims.push_back(levs);
    dims.push_back(rows);
    dims.push_back(cols);

    float ***camsData = allocate3d_float(nlev, camsLat, camsLon);

    int month = atoi(timeCoverageStart.substr(5, 2).c_str());

    string s = to_string(month);
    unsigned int numberOfZeros = 2 - s.length();
    s.insert(0, numberOfZeros, '0');
    s.insert(0, "CH4_");

    varin = camsFile->getVar(s.c_str());
    varin.getVar(&camsData[0][0][0]);
    // Convert CAMS fill values to -32767
    fptr = &camsData[0][0][0];
    for (size_t i = 0; i < nlev * camsLat * camsLon; i++) {
        if (*fptr == -999)
            *fptr = -32767;
        fptr++;
    }

    for (size_t l = 0; l < 42; l++) {
        cgal_interp2(ncams, &latCams[0][0], &lonCams[0][0], &camsData[l][0][0], nl1c, &outL1c3d[l][0][0],
                     nncTagCams);
    }
    // Convert interpolation failures to -32767
    fptr = &outL1c3d[0][0][0];
    for (size_t i = 0; i < nlev * naLong * naCross; i++) {
        if (*fptr < 0)
            *fptr = -32767;
        fptr++;
    }

    s.assign("CH4");
    cout << "Writing " << s << endl;

    // add the variable
    varout = ncOutput->addVar(s.c_str(), ncFloat, dims);

    // add its variable attributes by copying over from the source file
    copyVarAttsCams(&varin, &varout);

    // add the data
    varout.putVar(&outL1c3d[0][0][0]);

    delete camsFile;

    free2d_float(latCams);
    free2d_float(lonCams);
    free3d_float(camsData);

    cgal_release_tag(nncTagCams);
    /////////////////////////////////
    //////////// CAMS INST //////////
    /////////////////////////////////

    camsFile = new NcFile(camsco2Filename, NcFile::read);

    camsLatDim = camsFile->getDim("latitude_bins");
    camsLat = camsLatDim.getSize();

    camsLonDim = camsFile->getDim("longitude_bins");
    camsLon = camsLonDim.getSize();

    latCams = allocate2d_float(camsLat, camsLon);
    lonCams = allocate2d_float(camsLat, camsLon);

    for (size_t i = 0; i < camsLat; i++) {
        for (size_t j = 0; j < camsLon; j++) {
            latCams[i][j] = i * 1.8947368421 - 90.0;
            lonCams[i][j] = j * 3.75 - 180.0;

            if (dateLineCross && lonCams[i][j] < 0)
                lonCams[i][j] += 360;
        }
    }
    ncams = camsLon * camsLat;

    cout << endl << "Computing CAMS INST to L1C interpolation mapping" << endl;
    nncTagCams = cgal_nnc(ncams, &latCams[0][0], &lonCams[0][0], nl1c, &latL1C[0][0], &lonL1C[0][0]);

    camsData = allocate3d_float(nlev, camsLat, camsLon);

    s = to_string(month);
    numberOfZeros = 2 - s.length();
    s.insert(0, numberOfZeros, '0');
    s.insert(0, "CO2_");

    varin = camsFile->getVar(s.c_str());
    varin.getVar(&camsData[0][0][0]);

    // Convert CAMS fill values to -32767
    fptr = &camsData[0][0][0];
    for (size_t i = 0; i < nlev * camsLat * camsLon; i++) {
        if (*fptr == -999)
            *fptr = -32767;
        fptr++;
    }

    for (size_t l = 0; l < 42; l++) {
        cgal_interp2(ncams, &latCams[0][0], &lonCams[0][0], &camsData[l][0][0], nl1c, &outL1c3d[l][0][0],
                     nncTagCams);
    }

    // Convert interpolation failures to -32767
    fptr = &outL1c3d[0][0][0];
    for (size_t i = 0; i < nlev * naLong * naCross; i++) {
        if (*fptr < 0)
            *fptr = -32767;
        fptr++;
    }

    s.assign("CO2");
    cout << "Writing " << s << endl;

    varout = ncOutput->addVar(s.c_str(), ncFloat, dims);
    copyVarAttsCams(&varin, &varout);

    varout.putVar(&outL1c3d[0][0][0]);

    delete camsFile;

    camsFile = new NcFile(camsn2oFilename, NcFile::read);

    s = to_string(month);
    numberOfZeros = 2 - s.length();
    s.insert(0, numberOfZeros, '0');
    s.insert(0, "N2O_");

    varin = camsFile->getVar(s.c_str());
    varin.getVar(&camsData[0][0][0]);

    // Convert CAMS fill values to -32767
    fptr = &camsData[0][0][0];
    for (size_t i = 0; i < nlev * camsLat * camsLon; i++) {
        if (*fptr == -999)
            *fptr = -32767;
        fptr++;
    }

    for (size_t l = 0; l < 42; l++) {
        cgal_interp2(ncams, &latCams[0][0], &lonCams[0][0], &camsData[l][0][0], nl1c, &outL1c3d[l][0][0],
                     nncTagCams);
    }

    // Convert interpolation failures to -32767
    fptr = &outL1c3d[0][0][0];
    for (size_t i = 0; i < nlev * naLong * naCross; i++) {
        if (*fptr < 0)
            *fptr = -32767;
        fptr++;
    }

    s.assign("N2O");
    cout << "Writing " << s << endl;

    varout = ncOutput->addVar(s.c_str(), ncFloat, dims);
    copyVarAttsCams(&varin, &varout);

    varout.putVar(&outL1c3d[0][0][0]);

    delete camsFile;

    free2d_float(latCams);
    free2d_float(lonCams);
    free3d_float(camsData);

    cgal_release_tag(nncTagCams);
    /////////////////////////////////
    //////////// GEOS-CF ////////////
    /////////////////////////////////

    // Read and process GEOS-CF data

    cout << endl << "Opening GEOS-CF file" << endl;
    NcFile *geosFile = new NcFile(geoscfFilename, NcFile::read);

    nlon = 1440;
    nlat = 721;
    int ngeos = nlon * nlat;

    float **latGeos = allocate2d_float(721, 1440);
    float **lonGeos = allocate2d_float(721, 1440);

    // Compute lon/lat geos
    for (size_t i = 0; i < 721; i++) {
        for (size_t j = 0; j < 1440; j++) {
            latGeos[i][j] = i * 0.25 - 90;
            lonGeos[i][j] = j * 0.25 - 180;

            if (dateLineCross && lonGeos[i][j] < 0)
                lonGeos[i][j] += 360;
        }
    }

    // Write level field
    dims.clear();
    dims.push_back(levs);
    varout = ncOutput->addVar("levels", ncFloat, dims);

    float levels[nlev];
    varin = geosFile->getVar("lev");
    varin.getVar(levels);
    copyVarAtts(&varin, &varout);
    varout.putVar(levels);

    float ***geosData = allocate3d_float(42, 721, 1440);

    cout << "Computing GEOS to L1C interpolation mapping" << endl;
    int nncTagGeos = cgal_nnc(ngeos, &latGeos[0][0], &lonGeos[0][0], nl1c, &latL1C[0][0], &lonL1C[0][0]);

    string geosFields[5] = {"NO2", "SO2", "TOTCOL_NO2", "STRATCOL_NO2", "TROPCOL_NO2"};

    for (size_t k = 0; k < 5; k++) {
        varin = geosFile->getVar(geosFields[k].c_str());
        varin.getVar(&geosData[0][0][0]);

        if (geosFields[k].compare("NO2") == 0 || geosFields[k].compare("SO2") == 0)
            nlev = 42;
        else
            nlev = 1;

        // Convert GOES fill values to -1e14
        float *fptr = &geosData[0][0][0];
        for (size_t i = 0; i < nlev * 721 * 1440; i++) {
            if (*fptr > 1e14)
                *fptr = -1e14;
            fptr++;
        }

        for (size_t i = 0; i < nlev; i++) {
            cgal_interp2(ngeos, &latGeos[0][0], &lonGeos[0][0], &geosData[i][0][0], nl1c, &outL1c3d[i][0][0],
                         nncTagGeos);
        }

        // Convert interpolation failures to -32767
        fptr = &outL1c3d[0][0][0];
        for (size_t i = 0; i < nlev * naLong * naCross; i++) {
            if (*fptr < 0)
                *fptr = -32767;
            fptr++;
        }

        dims.clear();
        if (geosFields[k].compare("NO2") == 0 || geosFields[k].compare("SO2") == 0) {
            dims.push_back(levs);
            dims.push_back(rows);
            dims.push_back(cols);
        } else {
            dims.push_back(rows);
            dims.push_back(cols);
        }

        cout << "Writing " << geosFields[k].c_str() << endl;

        varout = ncOutput->addVar(geosFields[k].c_str(), ncFloat, dims);

        string overrideGeos[] = {"-32767", "=", "-32767", "=", "="};
        string overrideTypeGeos[] = {"F", "=", "F", "=", "="};
        copyVarAtts(&varin, &varout, overrideGeos, overrideTypeGeos);

        varout.putVar(&outL1c3d[0][0][0]);
    }

    cgal_release_tag(nncTagGeos);
    free2d_float(latGeos);
    free2d_float(lonGeos);
    free3d_float(geosData);
    free3d_float(outL1c3d);

    ///////////////////////////////
    //////////// CLOUD ////////////
    ///////////////////////////////
    if (strcmp(cloudMask2, "") != 0 && strcmp(cloudProd2, "") != 0) {
        cout << endl << "Opening CLDMSK files" << endl;

        NcFile *cloudMskFile[3] = {NULL, NULL, NULL};

        uint32_t nLines[3] = {0, 0, 0};
        uint32_t nPixels;
        NcDim linesDim;
        NcDim pixelDim;

        uint32_t nLines0 = 0;

        // Open trailing/following L2 CLDMASK datasets
        if (strcmp(cloudMask1, "") != 0) {
            cloudMskFile[0] = new NcFile(cloudMask1, NcFile::read);
            linesDim = cloudMskFile[0]->getDim("number_of_lines");
            nLines[0] = linesDim.getSize();
        }

        if (strcmp(cloudMask3, "") != 0) {
            cloudMskFile[2] = new NcFile(cloudMask3, NcFile::read);
            linesDim = cloudMskFile[2]->getDim("number_of_lines");
            nLines[2] = linesDim.getSize();
        }

        // Open current L2 LDMASK dataset
        cloudMskFile[1] = new NcFile(cloudMask2, NcFile::read);
        linesDim = cloudMskFile[1]->getDim("number_of_lines");
        nLines[1] = linesDim.getSize();

        pixelDim = cloudMskFile[1]->getDim("pixels_per_line");
        nPixels = pixelDim.getSize();

        // Use only 250 (or less) lines from trailing/following datasets
        uint32_t totLines = nLines[1];
        if (nLines[0] != 0) {
            if (nLines[0] < 250)
                totLines += nLines[0];
            else
                totLines += 250;
        }

        if (nLines[2] != 0) {
            if (nLines[2] < 250)
                totLines += nLines[2];
            else
                totLines += 250;
        }

        float **latCm = allocate2d_float(totLines, nPixels);
        float **lonCm = allocate2d_float(totLines, nPixels);
        int8_t **adjMask = allocate2d_schar(totLines, nPixels);

        for (size_t i = 0; i < totLines; i++) {
            for (size_t j = 0; j < nPixels; j++) {
                latCm[i][j] = 0.0;
                lonCm[i][j] = 0.0;
                adjMask[i][j] = 0;
            }
        }
        // Read from trailing granule if specified
        NcGroup ncGroup;

        if (cloudMskFile[0] != NULL) {
            if (nLines[0] < 250) {
                nLines0 = nLines[0];
                start[0] = 0;
            } else {
                nLines0 = 250;
                start[0] = nLines[0] - 250;
            }
            count[0] = nLines0;
            start[1] = 0;
            count[1] = nPixels;

            ncGroup = cloudMskFile[0]->getGroup("navigation_data");
            if (ncGroup.isNull()) {
                var = cloudMskFile[0]->getVar("latitude");
                var.getVar(start, count, &latCm[0][0]);

                var = cloudMskFile[0]->getVar("longitude");
                var.getVar(start, count, &lonCm[0][0]);

                var = cloudMskFile[0]->getVar("ADJ_MASK");
                var.getVar(start, count, &adjMask[0][0]);
            } else {
                var = ncGroup.getVar("latitude");
                var.getVar(start, count, &latCm[0][0]);

                var = ncGroup.getVar("longitude");
                var.getVar(start, count, &lonCm[0][0]);

                ncGroup = cloudMskFile[0]->getGroup("geophysical_data");
                var = ncGroup.getVar("cloud_flag");
                var.getVar(start, count, &adjMask[0][0]);
            }
        }

        ncGroup = cloudMskFile[1]->getGroup("navigation_data");
        if (ncGroup.isNull()) {
            var = cloudMskFile[1]->getVar("latitude");
            var.getVar(&latCm[nLines0][0]);

            var = cloudMskFile[1]->getVar("longitude");
            var.getVar(&lonCm[nLines0][0]);

            var = cloudMskFile[1]->getVar("ADJ_MASK");
            var.getVar(&adjMask[nLines0][0]);
        } else {
            var = ncGroup.getVar("latitude");
            var.getVar(&latCm[nLines0][0]);

            var = ncGroup.getVar("longitude");
            var.getVar(&lonCm[nLines0][0]);

            ncGroup = cloudMskFile[1]->getGroup("geophysical_data");
            var = ncGroup.getVar("cloud_flag");
            var.getVar(&adjMask[nLines0][0]);
        }

        // Read from following granule if specified
        if (cloudMskFile[2] != NULL) {
            start[0] = 0;
            if (nLines[2] < 250) {
                count[0] = nLines[2];
            } else {
                count[0] = 250;
            }
            start[1] = 0;
            count[1] = nPixels;

            ncGroup = cloudMskFile[2]->getGroup("navigation_data");
            if (ncGroup.isNull()) {
                var = cloudMskFile[2]->getVar("latitude");
                var.getVar(start, count, &latCm[nLines0 + nLines[1]][0]);

                var = cloudMskFile[2]->getVar("longitude");
                var.getVar(start, count, &lonCm[nLines0 + nLines[1]][0]);

                var = cloudMskFile[2]->getVar("ADJ_MASK");
                var.getVar(start, count, &adjMask[nLines0 + nLines[1]][0]);
            } else {
                var = ncGroup.getVar("latitude");
                var.getVar(start, count, &latCm[nLines0 + nLines[1]][0]);

                var = ncGroup.getVar("longitude");
                var.getVar(start, count, &lonCm[nLines0 + nLines[1]][0]);

                ncGroup = cloudMskFile[2]->getGroup("geophysical_data");
                var = ncGroup.getVar("cloud_flag");
                var.getVar(start, count, &adjMask[nLines0 + nLines[1]][0]);
            }
        }

        uint32_t ncm = totLines * nPixels;

        // Determine L1C row/col of L2 CLD lon/lat
        // Ported from code developed by F.Patt

        short *brow = new short[ncm];
        short *bcol = new short[ncm];
        lonlat2rowcol(naLong, naCross, ncm, gridRes, &lonL1C[0][0], &latL1C[0][0], &lonCm[0][0], &latCm[0][0],
                      brow, bcol);

        // Open CLD products file
        // Note: CLD files uses same lon/lat as CLDMSK

        NcFile *cloudFile[3];
        NcGroup ncGrp[3];

        if (strcmp(cloudProd1, "") != 0) {
            cloudFile[0] = new NcFile(cloudProd1, NcFile::read);
            ncGrp[0] = cloudFile[0]->getGroup("geophysical_data");
        }

        if (strcmp(cloudProd3, "") != 0) {
            cloudFile[2] = new NcFile(cloudProd3, NcFile::read);
            ncGrp[2] = cloudFile[2]->getGroup("geophysical_data");
        }

        cloudFile[1] = new NcFile(cloudProd2, NcFile::read);
        ncGrp[1] = cloudFile[1]->getGroup("geophysical_data");

        int8_t **cloudPhase = allocate2d_schar(totLines, nPixels);

        // Read cloud phase
        readCLD<int8_t>(ncGrp, "cld_phase_21", &cloudPhase[0][0], nPixels, nLines);

        for (size_t i = 0; i < totLines; i++)
            for (size_t j = 0; j < nPixels; j++)
                if (cloudPhase[i][j] == 3)
                    cloudPhase[i][j] = 1;
                else
                    cloudPhase[i][j] = 0;

        float **binvalIceCloud = allocate2d_float(naLong, naCross);
        float **binvalWaterCloud = allocate2d_float(naLong, naCross);
        short **nobsIceCloud = allocate2d_short(naLong, naCross);
        short **nobsWaterCloud = allocate2d_short(naLong, naCross);

        float **cloudProd = allocate2d_float(totLines, nPixels);

        // Average ADJ_MASK for water/ice cloud pixels

        // Convert BYTE to FLOAT
        for (size_t i = 0; i < totLines; i++)
            for (size_t j = 0; j < nPixels; j++)
                if (adjMask[i][j] == -128)
                    cloudProd[i][j] = -32767;
                else
                    cloudProd[i][j] = (float)adjMask[i][j];

        free2d_schar(adjMask);

        accumFrac(totLines, nPixels, naLong, naCross, brow, bcol, cloudProd, cloudPhase, binvalIceCloud,
                  nobsIceCloud, binvalWaterCloud, nobsWaterCloud);

        // Compute average bin values for ice clouds
        for (size_t i = 0; i < naLong; i++) {
            for (size_t j = 0; j < naCross; j++) {
                if (nobsIceCloud[i][j] != 0)
                    outL1c2d[i][j] = binvalIceCloud[i][j] / nobsIceCloud[i][j];
                else
                    outL1c2d[i][j] = -32767;
            }
        }

        dims.clear();
        dims.push_back(rows);
        dims.push_back(cols);

        cout << "Writing ice cloud fraction" << endl;
        varout = ncOutput->addVar("ice_cloud_fraction", ncFloat, dims);

        string override_icf[] = {"-32767", "", "", "", "", "Ice cloud fraction", "", "1.0", "0.0"};
        string overridetype_icf[] = {"F", "", "", "", "", "=", "", "F", "F"};
        copyVarAtts(&var, &varout, override_icf, overridetype_icf);

        varout.putVar(&outL1c2d[0][0]);

        // Compute average bin values for water clouds
        for (size_t i = 0; i < naLong; i++) {
            for (size_t j = 0; j < naCross; j++) {
                if (nobsWaterCloud[i][j] != 0)
                    outL1c2d[i][j] = binvalWaterCloud[i][j] / nobsWaterCloud[i][j];
                else
                    outL1c2d[i][j] = -32767;
            }
        }

        cout << "Writing water cloud fraction" << endl;
        varout = ncOutput->addVar("water_cloud_fraction", ncFloat, dims);

        string override_wcf[] = {"-32767", "", "", "", "", "Water cloud fraction", "", "1.0", "0.0"};
        string overridetype_wcf[] = {"F", "", "", "", "", "=", "", "F", "F"};
        copyVarAtts(&var, &varout, override_wcf, overridetype_wcf);

        varout.putVar(&outL1c2d[0][0]);
        // Cloud Top Pressure (CTP)
        // Cloud Top Temperature (CTT)
        // Cloud Top Height (CTH)
        // Particle effective radius cer_21 (micron)
        // COT cot_21
        string cloudFields[5] = {"ctp", "ctt", "cth", "cer_21", "cot_21"};

        for (size_t k = 0; k < 5; k++) {
            cout << "Writing " << cloudFields[k].c_str() << endl;

            var = ncGrp[1].getVar(cloudFields[k].c_str());

            for (size_t i = 0; i < totLines; i++)
                for (size_t j = 0; j < nPixels; j++)
                    cloudProd[i][j] = 0.0;

            readCLD<float>(ncGrp, cloudFields[k].c_str(), &cloudProd[0][0], nPixels, nLines);

            accum(totLines, nPixels, naLong, naCross, brow, bcol, cloudProd, cloudPhase, binvalIceCloud,
                  nobsIceCloud, binvalWaterCloud, nobsWaterCloud);

            // Compute average bin values for ice clouds
            for (size_t i = 0; i < naLong; i++) {
                for (size_t j = 0; j < naCross; j++) {
                    if (nobsIceCloud[i][j] != 0)
                        outL1c2d[i][j] = binvalIceCloud[i][j] / nobsIceCloud[i][j];
                    else
                        outL1c2d[i][j] = -32767;
                }
            }

            string outName;

            outName.assign(cloudFields[k]);
            outName.append("_ice_cloud");

            varout = ncOutput->addVar(outName.c_str(), ncFloat, dims);
            copyVarAtts(&var, &varout);
            varout.putVar(&outL1c2d[0][0]);

            // Compute average bin values for water clouds
            for (size_t i = 0; i < naLong; i++) {
                for (size_t j = 0; j < naCross; j++) {
                    if (nobsWaterCloud[i][j] != 0)
                        outL1c2d[i][j] = binvalWaterCloud[i][j] / nobsWaterCloud[i][j];
                    else
                        outL1c2d[i][j] = -32767;
                }
            }

            outName.assign(cloudFields[k]);
            outName.append("_water_cloud");

            varout = ncOutput->addVar(outName.c_str(), ncFloat, dims);
            copyVarAtts(&var, &varout);
            varout.putVar(&outL1c2d[0][0]);
        }

        delete[] brow;
        delete[] bcol;

        free2d_float(latCm);
        free2d_float(lonCm);
        free2d_float(cloudProd);
        free2d_float(binvalIceCloud);
        free2d_float(binvalWaterCloud);
        free2d_short(nobsIceCloud);
        free2d_short(nobsWaterCloud);
    }  // if CLD processing
    ///////////////////////////////
    ////////// WATERMASK //////////
    ///////////////////////////////

    float wmPixelSize = 360.0 / 86400;
    int wmLatRow = (latMin + 90) / wmPixelSize;
    int wmLonCol = (lonMin + 180) / wmPixelSize;

    start.clear();
    start.push_back(wmLatRow);
    start.push_back(wmLonCol);

    size_t nWmLatRow = ((latMax + 90) / wmPixelSize - wmLatRow + 1) / 4;
    size_t nWmLonCol = ((lonMax + 180) / wmPixelSize - wmLonCol + 1) / 4;

    count.clear();
    count.push_back(nWmLatRow);
    count.push_back(nWmLonCol);

    stride.clear();
    stride.push_back(4);
    stride.push_back(4);

    float *lonWm = new float[nWmLatRow * nWmLonCol];
    float *latWm = new float[nWmLatRow * nWmLonCol];

    int k = 0;
    for (size_t i = 0; i < (size_t)nWmLatRow; i++) {
        for (size_t j = 0; j < (size_t)nWmLonCol; j++) {
            latWm[k] = (4 * i + wmLatRow) * wmPixelSize - 90;
            lonWm[k] = (4 * j + wmLonCol) * wmPixelSize - 180;

            if (dateLineCross && lonWm[k] < 0)
                lonWm[k] += 360;
            k++;
        }
    }

    uint32_t nWaterMask = nWmLatRow * nWmLonCol;

    // Reading landmask
    string gebcoFileString = gebcoFilename;
    expandEnvVar(&gebcoFileString);

    NcFile *waterMaskFile = new NcFile(gebcoFileString.c_str(), NcFile::read);

    float *waterMask = new float[nWmLatRow * nWmLonCol];
    uint8_t *waterMaskByte = new uint8_t[nWmLatRow * nWmLonCol];
    var = waterMaskFile->getVar("watermask");
    if (dateLineCross) {
        uint8_t *wm0 = new uint8_t[nWmLatRow * nWmLonCol];
        count[1] = (86400 - wmLonCol) / 4;
        var.getVar(start, count, stride, &wm0[0]);
        for (size_t i = 0; i < (size_t)nWmLatRow; i++)
            memcpy(&waterMaskByte[i * nWmLonCol], &wm0[i * count[1]], count[1] * sizeof(uint8_t));

        start[1] = 0;
        count[1] = nWmLonCol - count[1];

        var.getVar(start, count, stride, &wm0[0]);
        for (size_t i = 0; i < (size_t)nWmLatRow; i++)
            memcpy(&waterMaskByte[i * nWmLonCol + (86400 - wmLonCol) / 4], &wm0[i * count[1]],
                   count[1] * sizeof(uint8_t));
        delete[] wm0;
    } else {
        var.getVar(start, count, stride, &waterMaskByte[0]);
    }

    // Convert byte to float
    k = 0;
    for (size_t i = 0; i < (size_t)nWmLatRow; i++) {
        for (size_t j = 0; j < (size_t)nWmLonCol; j++) {
            waterMask[k] = (float)waterMaskByte[k];
            k++;
        }
    }

    // Determine L1C row/col of L2 CLD lon/lat
    short *brow = new short[nWaterMask];
    short *bcol = new short[nWaterMask];

    float **binvalWm = allocate2d_float(naLong, naCross);
    short **nobsWm = allocate2d_short(naLong, naCross);

    // Clear accumulation arrays
    for (size_t i = 0; i < naLong; i++) {
        for (size_t j = 0; j < naCross; j++) {
            binvalWm[i][j] = 0.0;
            nobsWm[i][j] = 0;
        }
    }

    size_t N = 2147483648 / (nWmLonCol * naLong);
    size_t M = nWmLatRow / N;
    if (nWmLatRow % N != 0)
        M++;
    for (size_t i = 0; i < M; i++) {
        size_t L = N;
        if (nWmLatRow - i * N < N)
            L = nWmLatRow - i * N;

        lonlat2rowcol(naLong, naCross, L * nWmLonCol, gridRes, &lonL1C[0][0], &latL1C[0][0],
                      &lonWm[i * N * nWmLonCol], &latWm[i * N * nWmLonCol], brow, bcol);

        accumWm(L * nWmLonCol, brow, bcol, &waterMask[i * N * nWmLonCol], binvalWm, nobsWm);
    }

    for (size_t i = 0; i < naLong; i++) {
        for (size_t j = 0; j < naCross; j++) {
            if (nobsWm[i][j] != 0)
                outL1c2d[i][j] = binvalWm[i][j] / nobsWm[i][j];
            else
                outL1c2d[i][j] = -32767;
        }
    }

    cout << "Writing waterfraction" << endl;
    varout = ncOutput->addVar("waterfraction", ncFloat, dims);

    string overrideWaterFraction[] = {"-32767", "=", "Water fraction", "waterfraction", "1.0", "0.0"};
    string overrideTypeWaterFraction[] = {"F", "", "", "=", "F", "F"};
    copyVarAtts(&var, &varout, overrideWaterFraction, overrideTypeWaterFraction);
    varout.putVar(&outL1c2d[0][0]);

    delete[] latWm;
    delete[] lonWm;
    delete[] waterMask;

    delete[] brow;
    delete[] bcol;

    free2d_float(binvalWm);
    free2d_short(nobsWm);

    ///////////////////////////////
    //////////// ChlorOR ////////////
    ///////////////////////////////

    k = 0;

    if (strcmp(chlorFilename, "") != 0) {
        float ChlorPixelSize = 360.0 / 8640;
        int ChlorLatRow = (90 - latMax) / ChlorPixelSize;
        int ChlorLonCol = (lonMin + 180) / ChlorPixelSize;

        cout << latMin << " " << latMax << " " << lonMin << " " << lonMax << endl;
        cout << ChlorLatRow << " " << ChlorLonCol << endl;

        start.clear();
        start.push_back(ChlorLatRow);
        start.push_back(ChlorLonCol);

        int nChlorLatRow = (90 - latMin) / ChlorPixelSize - ChlorLatRow + 1;
        int nChlorLonCol = (lonMax + 180) / ChlorPixelSize - ChlorLonCol + 1;

        count.clear();
        count.push_back(nChlorLatRow);
        count.push_back(nChlorLonCol);

        float *lonChlor = new float[nChlorLatRow * nChlorLonCol];
        float *latChlor = new float[nChlorLatRow * nChlorLonCol];

        cout << nChlorLatRow << " " << nChlorLonCol << endl;

        k = 0;
        for (size_t i = 0; i < (size_t)nChlorLatRow; i++) {
            for (size_t j = 0; j < (size_t)nChlorLonCol; j++) {
                latChlor[k] = 90 - (i + ChlorLatRow) * ChlorPixelSize;
                lonChlor[k] = (j + ChlorLonCol) * ChlorPixelSize - 180;

                if (dateLineCross && lonChlor[k] < 0)
                    lonChlor[k] += 360;
                k++;
            }
        }

        int nChlor = nChlorLatRow * nChlorLonCol;

        cout << "Computing CHLOR_A map to L1C interpolation mapping" << endl;
        int nncTagChlor = cgal_nnc(nChlor, &latChlor[0], &lonChlor[0], nl1c, &latL1C[0][0], &lonL1C[0][0]);

        // Reading chlor_a from L3M file
        NcFile *ChlorFile = new NcFile(chlorFilename, NcFile::read);

        float *Chlor = new float[nChlorLatRow * nChlorLonCol];
        var = ChlorFile->getVar("chlor_a");

        if (dateLineCross) {
            float *Chlor0 = new float[nChlorLatRow * nChlorLonCol];
            count[1] = 8640 - ChlorLonCol;
            var.getVar(start, count, &Chlor0[0]);
            for (size_t i = 0; i < (size_t)nChlorLatRow; i++)
                memcpy(&Chlor[i * nChlorLonCol], &Chlor0[i * count[1]], count[1] * sizeof(float));

            start[1] = 0;
            count[1] = nChlorLonCol + ChlorLonCol - 8640;
            cout << count[1] << endl;
            var.getVar(start, count, &Chlor0[0]);
            for (size_t i = 0; i < (size_t)nChlorLatRow; i++)
                memcpy(&Chlor[i * nChlorLonCol + (8640 - ChlorLonCol)], &Chlor0[i * count[1]],
                       count[1] * sizeof(float));
            delete[] Chlor0;
        } else {
            var.getVar(start, count, &Chlor[0]);
        }

        cgal_interp2(nChlor, &latChlor[0], &lonChlor[0], &Chlor[0], nl1c, &outL1c2d[0][0], nncTagChlor);

        fptr = &outL1c2d[0][0];
        for (size_t i = 0; i < naLong * naCross; i++) {
            if (*fptr < 0)
                *fptr = -32767;
            fptr++;
        }
        // Convert interpolation failures to -32767
        cout << "Writing chlor_a" << endl;
        varout = ncOutput->addVar("chlor_a", ncFloat, dims);
        copyVarAtts(&var, &varout);
        varout.putVar(&outL1c2d[0][0]);

        cgal_release_tag(nncTagChlor);
        delete[] latChlor;
        delete[] lonChlor;
        delete[] Chlor;
    }

    free2d_float(outL1c2d);
    free2d_float(latL1C);
    free2d_float(lonL1C);

    return 0;
}

int lonlat2rowcol(uint32_t naLong, uint32_t naCross, uint32_t ncm, float gridRes, float *lonL1C,
                  float *latL1C, float *lonCm, float *latCm, short *brow, short *bcol) {
    // Determine L1C row/col of L2 CLD lon/lat
    // Ported from code developed by F.Patt
    float *fptrLon, *fptrLat;
    // Convert lon/lat to unit vectors (L1C)
    float ***gvec = allocate3d_float(naLong, naCross, 3);
    fptrLon = lonL1C;
    fptrLat = latL1C;
    for (size_t i = 0; i < naLong; i++) {
        for (size_t j = 0; j < naCross; j++) {
            gvec[i][j][0] = cos(*fptrLon / RADEG) * cos(*fptrLat / RADEG);
            gvec[i][j][1] = sin(*fptrLon / RADEG) * cos(*fptrLat / RADEG);
            gvec[i][j][2] = sin(*fptrLat / RADEG);

            fptrLon++;
            fptrLat++;
        }
    }
    // Convert lon/lat to unit vectors (L2 CLD)
    float **bvec = allocate2d_float(ncm, 3);
    fptrLon = lonCm;
    fptrLat = latCm;
    for (size_t i = 0; i < ncm; i++) {
        bvec[i][0] = cos(*fptrLon / RADEG) * cos(*fptrLat / RADEG);
        bvec[i][1] = sin(*fptrLon / RADEG) * cos(*fptrLat / RADEG);
        bvec[i][2] = sin(*fptrLat / RADEG);

        fptrLon++;
        fptrLat++;
    }

    // Compute normal vectors for L1C rows
    // (cross product of first and last vector in each row)
    float **gnvec = allocate2d_float(naLong, 3);
    float *gnvm = new float[naLong];
    float vecm[3];
    for (size_t i = 0; i < naLong; i++) {
        gnvec[i][0] = gvec[i][naCross - 1][1] * gvec[i][0][2] - gvec[i][naCross - 1][2] * gvec[i][0][1];

        gnvec[i][1] = gvec[i][naCross - 1][2] * gvec[i][0][0] - gvec[i][naCross - 1][0] * gvec[i][0][2];

        gnvec[i][2] = gvec[i][naCross - 1][0] * gvec[i][0][1] - gvec[i][naCross - 1][1] * gvec[i][0][0];

        vecm[0] = gnvec[i][0];
        vecm[1] = gnvec[i][1];
        vecm[2] = gnvec[i][2];

        gnvm[i] = sqrt(vecm[0] * vecm[0] + vecm[1] * vecm[1] + vecm[2] * vecm[2]);
    }
    for (size_t i = 0; i < naLong; i++) {
        gnvec[i][0] /= gnvm[i];
        gnvec[i][1] /= gnvm[i];
        gnvec[i][2] /= gnvm[i];
    }
    delete[] gnvm;

    // Compute dot products of L1B vectors with normal vectors
    gsl_matrix_float_view G = gsl_matrix_float_view_array(&gnvec[0][0], naLong, 3);
    gsl_matrix_float_view B = gsl_matrix_float_view_array(&bvec[0][0], ncm, 3);

    float **bdotgn = allocate2d_float(ncm, naLong);
    gsl_matrix_float_view bdotgnMat = gsl_matrix_float_view_array(&bdotgn[0][0], ncm, naLong);

    gsl_blas_sgemm(CblasNoTrans, CblasTrans, 1.0, &B.matrix, &G.matrix, 0.0, &bdotgnMat.matrix);
    // Determine row and column for each pixel within grid
    // Compute normals to center columns
    uint32_t ic = naCross / 2;
    float **cnvec = allocate2d_float(naLong, 3);
    for (size_t i = 0; i < naLong; i++) {
        cnvec[i][0] = gvec[i][ic][1] * gnvec[i][2] - gvec[i][ic][2] * gnvec[i][1];
        cnvec[i][1] = gvec[i][ic][2] * gnvec[i][0] - gvec[i][ic][0] * gnvec[i][2];
        cnvec[i][2] = gvec[i][ic][0] * gnvec[i][1] - gvec[i][ic][1] * gnvec[i][0];
    }
    // Compute grid row nadir resolution
    float *dcm = new float[naLong];
    for (size_t i = 0; i < naLong; i++) {
        vecm[0] = gvec[i][ic + 1][0] - gvec[i][ic][0];
        vecm[1] = gvec[i][ic + 1][1] - gvec[i][ic][1];
        vecm[2] = gvec[i][ic + 1][2] - gvec[i][ic][2];

        dcm[i] = sqrt(vecm[0] * vecm[0] + vecm[1] * vecm[1] + vecm[2] * vecm[2]);
    }

    float db = gridRes / 6311 / 2;

    for (size_t i = 0; i < ncm; i++) {
        if (bdotgn[i][0] > db || bdotgn[i][naLong - 1] < -db) {
            brow[i] = -1;
            bcol[i] = -1;

        } else {
            brow[i] = -1;
            bcol[i] = -1;

            for (size_t j = 0; j < naLong; j++) {
                if (bdotgn[i][j] > -db && bdotgn[i][j] <= db) {
                    brow[i] = j;

                    float bdotcn = 0.0;
                    for (size_t k = 0; k < 3; k++)
                        bdotcn += (bvec[i][k] * cnvec[j][k]);

                    bcol[i] = (short)round(ic + bdotcn / dcm[j]);

                    if (bcol[i] < 0 || bcol[i] >= (short)naCross) {
                        bcol[i] = -1;
                        brow[i] = -1;
                    }
                    break;
                }
            }
        }
    }

    free3d_float(gvec);
    free2d_float(bvec);
    free2d_float(gnvec);
    free2d_float(cnvec);
    free2d_float(bdotgn);

    delete[] dcm;

    return 0;
}
