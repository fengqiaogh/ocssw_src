/*
 * l3mapgen.cpp
 *
 *  Created on: Jan 15, 2014
 *      Author: dshea
 */
#include "l3mapgen.h"

#include "OutFile.h"
#include "OutFilePgm.h"
#include "OutFilePpm.h"
#include "OutFilePpmRgb.h"
#include "OutFilePng.h"
#include "OutFilePngRgb.h"
#include "OutFileTiff.h"
#include "OutFileTiffRgb.h"
#include "OutFileTiffGray.h"
#include "OutFileTiffColor.h"
#include "OutFileHdf4.h"
#include "OutFileNetcdf4.h"
#include <L3FileSMI.h>

#include <iostream>
#include <iomanip>
#include <fstream>

#include <genutils.h>
#include <timeutils.h>
#include <nc4utils.h>
#include <string>

#include <proj.h>
#include <sensorInfo.h>
#include <productInfo.h>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <unordered_map>
#include <unordered_set>
#include <boost/algorithm/string.hpp>
#include "Wave3DParsing.hpp"
#include "CacheSize.h"

using namespace std;
using namespace l3;

#define NUM_SEARCH_POINTS 51

/** Measurement types for different data processing */
enum MeasurementType { AVG, STDEV, VARIANCE, NOBS, NSCENES, OBS_TIME, BIN_NUM };

/** Interpolation types for mapping data */
enum InterpType { INTERP_NEAREST, INTERP_BIN, INTERP_LINEAR, INTERP_AREA };

/** Constants for East and West Checks */
enum EastWest { NOT_EAST_OR_WEST, IS_EAST, IS_WEST };


/**
 * @brief Gets the modifier string based on the measurement type.
 * @param measure The type of measurement.
 * @return The modifier string.
 */
string getModifierFromMeasure(MeasurementType measure) {
    switch (measure) {
        case AVG:
            return "";
        case STDEV:
            return "_stdev";
        case VARIANCE:
            return "_var";
        case NOBS:
            return "_nobs";
        case NSCENES:
            return "_nscenes";
        case OBS_TIME:
            return "_obs_time";
        case BIN_NUM:
            return "_bin_num";
        default:
            return "";
    }
}
// global params for program
double widthInDeg;
double heightInDeg;
double mapEast; /** Map east boundary */
double mapWest; /** Map west boundary */

int imageWidth = 0; 
int imageHeight = 0;
bool doingQuality = false; /**Flag for quality processing */
bool doingRGB = false; /** Flag for RGB processing */
bool doingTransparency = false; /** Flag for transparency processing */
bool trimNSEW = false; /** Flag for trimming North, South, East, West */
bool writeProjectionText = false; /** Flag for writing projection text */

string productName;
vector<string> productNameList;
vector<MeasurementType> productMeasurementList; /** List of measurement types for products */
string qualityProductName; /** Quality product name */

clo_optionList_t* optionList = NULL;

// landmask
bool applyMask = false; /** Flag to apply landmask */
static grid_info_t* landmaskGrid = {0}; // TODO: apply naming convention to typedef in /oel_util/libnetcdfutils/nc_gridutils.h

/**
 * @brief Checks if the given latitude and longitude is land or water.
 * @param lat Latitude
 * @param lon Longitude
 * @return 1 if water, 0 if land
 */
int isLand(float lat, float lon) {
    if (landmaskGrid != NULL) {
        double value;
        int status = get_bylatlon(landmaskGrid, lat, lon, &value);
        if (!status)
            return ((short)value != 1);  // 1 = water
    }
    return (0);  // assume water if lon, lat not found
}

/**
 * @brief Prints start information for output files.
 * @param outFiles Vector of output files.
 */
void printStartInfo(vector<OutFile*> outFiles) {
    if (want_verbose) {
        meta_l3bType* metaData = outFiles[0]->getMetadata();

        clo_printVersion();
        printf("ifile      : %s\n", clo_getString(optionList, "ifile"));
        printf("ofile      : ");
        for (OutFile* outFile : outFiles) {
            if (outFile == outFiles[0])
                printf("%s", outFile->getFileName().c_str());
            else
                printf(",%s", outFile->getFileName().c_str());
        }
        printf("\n");

        printf("oformat    : %s\n", clo_getString(optionList, "oformat"));
        if (clo_isSet(optionList, "ofile2")) {
            printf("ofile2     : %s\n", clo_getString(optionList, "ofile2"));
            printf("oformat2   : %s\n", clo_getString(optionList, "oformat2"));
        }
        printf("projection : %s\n", clo_getRawString(optionList, "projection"));
        printf("resolution : %.3fm\n", outFiles[0]->getResolution());
        if (clo_isSet(optionList, "width"))
            printf("width      : %d\n", clo_getInt(optionList, "width"));
        if (doingRGB)
            printf("product_rgb: %s\n", productName.c_str());
        else
            printf("product    : %s\n", productName.c_str());
        if (doingQuality)
            printf("qual_prod  : %s\n", qualityProductName.c_str());
        printf("north      : %8.3f\n", metaData->north);
        printf("south      : %8.3f\n", metaData->south);
        printf("east       : %8.3f\n", metaData->east);
        printf("west       : %8.3f\n", metaData->west);
        float tmpf = clo_getFloat(optionList, "central_meridian");
        if (tmpf > -900.0) {
            printf("central_meridian : %8.3f\n", tmpf);
        }
        if (clo_isSet(optionList, "lat_ts")) {
            printf("lat_ts     : %8.3f\n", clo_getFloat(optionList, "lat_ts"));
        }
        if (clo_isSet(optionList, "lat_0")) {
            printf("lat_0      : %8.3f\n", clo_getFloat(optionList, "lat_0"));
        }
        if (clo_isSet(optionList, "azimuth")) {
            printf("azimuth    : %8.3f\n", clo_getFloat(optionList, "azimuth"));
        }
        printf("image size : %d x %d\n", imageHeight, imageWidth);

        printf("\n");
    }
}

/**
 * @brief Prints end information for the output file.
 * @param outFile Pointer to the output file.
 */
void printEndInfo(OutFile* outFile) {
    if (want_verbose) {
        printf("\n\n");
        printf("actual data min       : %f\n", outFile->getFileMinVal());
        printf("actual data max       : %f\n", outFile->getFileMaxVal());
        printf("num filled pixels     : %d\n", outFile->getNumFilledPixels());
        printf("percent filled pixels : %.2f%%\n", outFile->getPercentFilledPixels());
        printf("\n");
    }
}

/**
 * @brief Prints the percent completion of the process.
 * @param percentDone The percentage of completion.
 */
void printPercentDone(float percentDone) {
    static float percentPrev = 0.0;
    static const float percentDelta = 0.01;
    if (want_verbose && (percentDone - percentPrev > percentDelta)) {
        percentPrev = percentDone;
        printf("\r%2d%% complete", (int)(percentDone * 100));
        fflush(stdout);
    }
}

/**
 * @brief Gets the central meridian value from the option list or metadata (mapWest and mapEast).
 * @return The central meridian value.
 */
float getCentralMeridian() {
    int i;
    float centralMeridian = clo_getFloat(optionList, "central_meridian");
    if (centralMeridian > -900.0) {
        i = 0;
        while (centralMeridian < -180.0) {
            centralMeridian += 360.0;
            i++;
            if (i > 5) {
                printf("-E- central meridian is way off\n");
                exit(EXIT_FAILURE);
            }
        }
        i = 0;
        while (centralMeridian > 180.0) {
            centralMeridian -= 360.0;
            i++;
            if (i > 5) {
                printf("-E- central meridian is way off\n");
                exit(EXIT_FAILURE);
            }
        }
    } else {
        centralMeridian = (mapWest + mapEast) / 2.0;
    }
    return constrainLon(centralMeridian);
}


/**
 * @brief Converts interpolation type string to InterpType enum.
 * @param str Interpolation type string.
 * @return Corresponding InterpType enum value.
 */
InterpType interpStr2Type(const char* str) {
    string s = str;
    boost::trim(s);
    boost::to_lower(s);

    if (s == "bin")
        return INTERP_BIN;
    if (s == "linear")
        return INTERP_LINEAR;
    if (s == "area")
        return INTERP_AREA;

    return INTERP_NEAREST;
}


/**
 * @brief Converts InterpType enum to string.
 * @param interp Interpolation type enum.
 * @return Corresponding string value.
 */
const char* interpType2Str(InterpType interp) {
    switch (interp) {
        case INTERP_NEAREST:
            return "nearest";
        case INTERP_BIN:
            return "bin";
        case INTERP_LINEAR:
            return "linear";
        case INTERP_AREA:
            return "area";
        default:
            return "unknown";
    }
}

/**
 * @brief Checks if the dateline is crossed given longitude and delta longitude.
 * @param lon Longitude.
 * @param deltaLon Delta longitude.
 * @return True if dateline is crossed, false otherwise.
 */
bool checkDateLineCrossed(double lon, double deltaLon) {
    float minlon = constrainLon(lon) - (deltaLon / 2.0);
    float maxlon = constrainLon(lon) + (deltaLon / 2.0);

    return (minlon > 0 && maxlon < 0);
}

/**
 * @brief Gets a bounding box given latitude, longitude, delta latitude, and delta longitude.
 * @param lat Latitude.
 * @param lon Longitude.
 * @param deltaLat Delta latitude.
 * @param deltaLon Delta longitude.
 * @param eastwest East/West enum value (default: NOT_EAST_OR_WEST).
 * @return Bounding box.
 */
Box_t getBox(float lat, float lon, float deltaLat, float deltaLon, int eastwest = NOT_EAST_OR_WEST) {
    Point_t pMin;
    Point_t pMax;
    lon = constrainLon(lon);
    if (eastwest == IS_EAST) {
        pMin.set<0>(180.0);
        pMin.set<1>(lat - (deltaLat / 2.0));
        pMax.set<0>(lon + (deltaLon / 2.0));
        pMax.set<1>(lat + (deltaLat / 2.0));
    } else if (eastwest == IS_WEST) {
        pMin.set<0>(lon - (deltaLon / 2.0));
        pMin.set<1>(lat - (deltaLat / 2.0));
        pMax.set<0>(180.0);
        pMax.set<1>(lat + (deltaLat / 2.0));
    } else {
        pMin.set<0>(lon - (deltaLon / 2.0));
        pMin.set<1>(lat - (deltaLat / 2.0));
        pMax.set<0>(lon + (deltaLon / 2.0));
        pMax.set<1>(lat + (deltaLat / 2.0));
    }
    return Box_t(pMin, pMax);
}

/**
 * @brief Gets bins inside a bounding box.
 * @param l3File Pointer to L3File.
 * @param lat Latitude.
 * @param lon Longitude.
 * @param deltaLat Delta latitude.
 * @param deltaLon Delta longitude.
 * @param fudge Fudge factor.
 * @param areaWeighted Area-weighted flag.
 * @param eastwest East/West enum value (default: NOT_EAST_OR_WEST).
 * @return Pointer to L3Bin.
 */
L3Bin* getBoxBins(L3File* l3File, float lat, float lon, float deltaLat, float deltaLon, float fudge,
                  bool areaWeighted, int eastwest = NOT_EAST_OR_WEST) {
    Box_t box = getBox(lat, lon, deltaLat, deltaLon, eastwest);

    L3Bin* l3Bin = l3File->getBinsInside(box, areaWeighted);

    if (!l3Bin && fudge > 1) {  // try again with fudge factor
        Box_t box = getBox(lat, lon, deltaLat * fudge, deltaLon * fudge, eastwest);
        l3Bin = l3File->getBinsInside(box, areaWeighted);
    }
    return l3Bin;
}

/**
 * @brief Figure out if we want and can do quality processing
 * @param l3File Pointer to L3File.
 * @param outFiles Vector of output files.
 * @param outFile2 Pointer to the second output file.
 * @return True if quality processing is set up, else false.
 */
bool setupQualityProcessing(L3File* l3File, vector<OutFile*> outFiles, OutFile* outFile2) {
    doingQuality = true;
    clo_option_t* option = clo_findOption(optionList, "use_quality");
    if (clo_isOptionSet(option)) {
        if (clo_getOptionBool(option)) {
            if (l3File->hasQuality()) {
                doingQuality = true;
            } else {
                printf(
                    "-E- Quality processing was requested, "
                    "but the input file does not have quality data.\n");
                exit(EXIT_FAILURE);
            }
            if (clo_isSet(optionList, "quality_product")) {
                qualityProductName = clo_getRawString(optionList, "quality_product");
            } else {
                qualityProductName = "qual_" + productNameList[0];
            }
        } else {
            doingQuality = false;
        }
    } else {
        if (l3File->hasQuality())
            doingQuality = true;
        else
            doingQuality = false;
    }

    l3File->setQualityProcessing(doingQuality);
    for (OutFile* outFile : outFiles) {
        outFile->setQualityProcessing(doingQuality);
        outFile->setQualityName(qualityProductName);
    }
    if (outFile2) {
        outFile2->setQualityProcessing(doingQuality);
        outFile2->setQualityName(qualityProductName);
    }
    return doingQuality;
}

/**
 * @brief Sets up product information for output files.
 * @param productName Product name.
 * @param measure Measurement type.
 * @param outFile Pointer to the output file.
 * @param outFile2 Pointer to the second output file.
 * @param productAttr Product attribute.
 */
void setupProduct(string& productName, MeasurementType measure, OutFile* outFile, OutFile* outFile2, const ProductL3Attributes & productAttr = ProductL3Attributes()) {
    // get the product info
    productInfo_t* productInfo = allocateProductInfo();
    int sensorId = sensorName2SensorId(outFile->getMetadata()->sensor_name);

    if (sensorId == -1) {
        cerr << "-E -Uknown sensor name " << outFile->getMetadata()->sensor_name << endl;
        freeProductInfo(productInfo); // clean up before exit
        exit(EXIT_FAILURE);
    }
    
    bool isWavelengthProduct = false;
    // check if product has wavelength info
    if (productName.find("_") != string::npos) {
        vector<string> parts;
        boost::split(parts, productName, boost::is_any_of("_"));
        // check if any part is a number that could be a wavelength
        for (const auto& part : parts) {
            if (!part.empty() && std::all_of(part.begin(), part.end(), ::isdigit)) {
                isWavelengthProduct = true;
                break;
            }
        }
    }

    if (!findProductInfo(productName.c_str(), sensorId, productInfo)) {
        if (isWavelengthProduct) {
            setupWavelengthProduct(productName, productInfo, sensorId);
        } else {
        // Create metadata if products aren't found in XML product table.
        printf("-E - product %s not found XML product table\n. Creating metadata.", productName.c_str());
        productInfo->description = strdup("Dynamically generated product");
        productInfo->units = strdup("unitless");
        productInfo->dataType = strdup("float");
        productInfo->validMin = 0;
        productInfo->validMax = 100;
        }
    }

    // now we have to fix the productInfo structure
    // AVG, STDEV, VARIANCE, NOBS, NSCENES, OBS_TIME, BIN_NUM
    string tmpStr;
    switch (measure) {
        case AVG:
            break;
        case STDEV:
            tmpStr = productInfo->description;
            productInfo->description = strdup((tmpStr + " (Standard Deviation)").c_str());
            productInfo->palette = strdup(PRODUCT_DEFAULT_palette);
            productInfo->dataType = strdup("float");
            productInfo->suffix = strdup((string(productInfo->suffix) + "_stdev").c_str());
            productInfo->displayScale = strdup("linear");
            productInfo->addOffset = 0.0;
            productInfo->scaleFactor = 1.0;
            break;
        case VARIANCE:
            tmpStr = productInfo->description;
            productInfo->description = strdup((tmpStr + " (Variance)").c_str());
            productInfo->palette = strdup(PRODUCT_DEFAULT_palette);
            productInfo->dataType = strdup("float");
            productInfo->suffix = strdup((string(productInfo->suffix) + "_var").c_str());
            productInfo->displayScale = strdup("linear");
            productInfo->addOffset = 0.0;
            productInfo->scaleFactor = 1.0;
            break;
        case NOBS:
            tmpStr = productInfo->description;
            productInfo->description = strdup((tmpStr + " (number of observations)").c_str());
            productInfo->units = strdup("counts");
            productInfo->palette = strdup(PRODUCT_DEFAULT_palette);
            productInfo->dataType = strdup("short");
            productInfo->suffix = strdup((string(productInfo->suffix) + "_nobs").c_str());
            productInfo->fillValue = PRODUCT_DEFAULT_fillValue;
            productInfo->validMin = 0;
            productInfo->validMax = 32767;
            productInfo->displayScale = strdup("linear");
            productInfo->displayMin = PRODUCT_DEFAULT_displayMin;
            productInfo->displayMax = PRODUCT_DEFAULT_displayMax;
            productInfo->addOffset = 0.0;
            productInfo->scaleFactor = 1.0;
            break;
        case NSCENES:
            tmpStr = productInfo->description;
            productInfo->description = strdup((tmpStr + " (number of scenes)").c_str());
            productInfo->units = strdup("counts");
            productInfo->palette = strdup(PRODUCT_DEFAULT_palette);
            productInfo->dataType = strdup("short");
            productInfo->suffix = strdup((string(productInfo->suffix) + "_nscenes").c_str());
            productInfo->fillValue = PRODUCT_DEFAULT_fillValue;
            productInfo->validMin = 0;
            productInfo->validMax = 32767;
            productInfo->displayScale = strdup("linear");
            productInfo->displayMin = PRODUCT_DEFAULT_displayMin;
            productInfo->displayMax = PRODUCT_DEFAULT_displayMax;
            productInfo->addOffset = 0.0;
            productInfo->scaleFactor = 1.0;
            break;
        case OBS_TIME:
            tmpStr = productInfo->description;
            productInfo->description = strdup((tmpStr + " (observation time, TAI93)").c_str());
            productInfo->units = strdup("counts");
            productInfo->palette = strdup(PRODUCT_DEFAULT_palette);
            productInfo->dataType = strdup("float");
            productInfo->suffix = strdup((string(productInfo->suffix) + "_obs_time").c_str());
            productInfo->fillValue = PRODUCT_DEFAULT_fillValue;
            productInfo->validMin = PRODUCT_DEFAULT_validMin;
            productInfo->validMax = PRODUCT_DEFAULT_validMax;
            productInfo->displayScale = strdup("linear");
            productInfo->displayMin = PRODUCT_DEFAULT_displayMin;
            productInfo->displayMax = PRODUCT_DEFAULT_displayMax;
            productInfo->addOffset = 0.0;
            productInfo->scaleFactor = 1.0;
            break;
        case BIN_NUM:
            tmpStr = productInfo->description;
            productInfo->description = strdup((tmpStr + " (bin ID number)").c_str());
            productInfo->units = strdup("dimensionless");
            productInfo->palette = strdup(PRODUCT_DEFAULT_palette);
            productInfo->dataType = strdup("int");
            productInfo->suffix = strdup((string(productInfo->suffix) + "_bin_num").c_str());
            productInfo->fillValue = PRODUCT_DEFAULT_fillValue;
            productInfo->validMin = PRODUCT_DEFAULT_validMin;
            productInfo->validMax = PRODUCT_DEFAULT_validMax;
            productInfo->displayScale = strdup("linear");
            productInfo->displayMin = PRODUCT_DEFAULT_displayMin;
            productInfo->displayMax = PRODUCT_DEFAULT_displayMax;
            productInfo->addOffset = 0.0;
            productInfo->scaleFactor = 1.0;
            break;
    }  // switch
    
    // load landmask options before palette
    applyMask = clo_getBool(optionList, "mask_land");
    if (applyMask) {
        char* strVal = clo_getOptionString(clo_findOption(optionList, "land"));
        char landmaskFile[FILENAME_MAX];
        parse_file_name(strVal, landmaskFile);
        static const char* landmaskVars[] = {"watermask", "landmask", "z", NULL};
        landmaskGrid = allocate_gridinfo();
        int status = init_gridinfo(landmaskFile, landmaskVars, landmaskGrid);
        if (status != NC_NOERR) {
            free(landmaskGrid);
            landmaskGrid = NULL;
            cerr << "Error reading file " << landmaskFile << ": " << nc_strerror(status) << "\n"
                 << "Land mask will not be applied.\n";
            applyMask = false;
        }
    }

    // load palette
    if (clo_getBool(optionList, "apply_pal")) {
        // set land_rgb before loading palette
        outFile->setLandRGB(clo_getRawString(optionList, "rgb_land"));
        if (outFile2)
            outFile2->setLandRGB(clo_getRawString(optionList, "rgb_land"));

        clo_option_t* option = clo_findOption(optionList, "palfile");
        if (clo_isOptionSet(option)) {
            outFile->setPalette(clo_getOptionString(option), applyMask);
            if (outFile2)
                outFile2->setPalette(clo_getOptionString(option), applyMask);
        } else {
            outFile->setPalette(productInfo->palette, applyMask);
            if (outFile2)
                outFile2->setPalette(productInfo->palette, applyMask);
        }
    }

    // set the default scale factors for RGB
    if (doingRGB) {
        productInfo->displayMin = 0.01;
        productInfo->displayMax = 0.9;
        productInfo->displayScale = strdup("log");
    }

    // override default scale parameters if set on command line
    if (clo_isSet(optionList, "datamin")) {
        productInfo->displayMin = clo_getFloat(optionList, "datamin");
    }
    if (clo_isSet(optionList, "datamax")) {
        productInfo->displayMax = clo_getFloat(optionList, "datamax");
    }
    if (clo_isSet(optionList, "scale_type")) {
        productInfo->displayScale = strdup(clo_getString(optionList, "scale_type"));
    }

    outFile->addProduct(productInfo, applyMask, productAttr);
    if (outFile2)
        outFile2->addProduct(productInfo, applyMask, productAttr);

    freeProductInfo(productInfo);
}

/**
* @brief Parses product name to find "base" name (eg: rrs) and wavelength value and creates product metadata for wavelength-dependent variables that aren't found in product.xml
**/
void setupWavelengthProduct(const string& productName, productInfo_t* productInfo, int sensorId) {
    vector<string> parts;
    boost::split(parts, productName, boost::is_any_of("_"));
    string baseProduct;
    string wavelength;
    
    for (const auto& part : parts) {
        if (std::all_of(part.begin(), part.end(), ::isdigit)) {
            wavelength = part;
            productInfo->prod_ix = stoi(wavelength); // use prod_ix for storing wavelength
        } else {
            if (!baseProduct.empty()) baseProduct += "_";
            baseProduct += part;
        }
    }
    
    productInfo->productName = strdup(productName.c_str());
    productInfo->description = strdup((baseProduct + " at " + wavelength).c_str());
    productInfo->units = strdup("dimensionless"); 
    productInfo->dataType = strdup("float");
    productInfo->displayScale = strdup("linear");
}


/**
 * @brief Writes raw output file.
 * @param l3File Pointer to L3File.
 * @param outFiles Vector of output files.
 * @param outFile2 Pointer to the second output file.
 */
void writeRawFile(L3File* l3File, vector<OutFile*> outFiles, OutFile* outFile2) {
    L3Bin* l3Bin;
    L3Row* l3Row;
    imageHeight = l3File->getNumRows();
    double resolution = EARTH_CIRCUMFERENCE / (imageHeight * 2);
    imageWidth = imageHeight * 2;
    int32_t start;
    int32_t numBins;
    int64_t baseBin;
    int64_t endBin;
    int64_t binNum;
    int32_t row, col;

    int32_t numFilledPixels = 0;
    meta_l3bType* metaData = outFiles[0]->getMetadata();

    string mapDesc = "Bin";
    string projName = "Integerized Sinusoidal";

    string tmpTitle = metaData->sensor_name + (string)" Level-3 " + mapDesc + (string)" Mapped Image";
    if (clo_isSet(optionList, "suite"))
        tmpTitle += " " + (string)clo_getString(optionList, "suite");
    strcpy(metaData->title, tmpTitle.c_str());

    for (OutFile* outFile : outFiles) {
        strcpy(outFile->getMetadata()->title, metaData->title);
        outFile->setFullLatLon(false);
        outFile->setMapProjection(projName);
        outFile->setResolution(resolution);
        outFile->setSize(imageWidth, imageHeight);
    }
    if (outFile2) {
        strcpy(outFile2->getMetadata()->title, metaData->title);
        outFile2->setFullLatLon(false);
        outFile2->setMapProjection(projName);
        outFile2->setResolution(resolution);
        outFile2->setSize(imageWidth, imageHeight);
    }

    // setup all the product structures
    for (size_t i = 0; i < productNameList.size(); i++) {
        OutFile* outFile;
        if (outFiles.size() == 1)
            outFile = outFiles[0];
        else
            outFile = outFiles[i];
        const ProductL3Attributes & productAttr = readProductL3Attributes(l3File, productNameList[i]);
        setupProduct(productNameList[i], productMeasurementList[i], outFile, outFile2, productAttr);
    }

    printStartInfo(outFiles);
    // set cache size for writing
    if (!set_cache_size_chunk_size_write(l3File,outFiles[0],outFile2,want_verbose)) {
        fprintf(stderr, "-E-: %s:%d Could not set cache size for writing.\n",__FILE__,__LINE__);
        exit(EXIT_FAILURE);
    }
    for (OutFile* outFile : outFiles) {
        if (!outFile->open()) {
            printf("-E- Could not open ofile=\"%s\".\n", outFile->getFileName().c_str());
            exit(EXIT_FAILURE);
        }
    }
    if (outFile2) {
        if (!outFile2->open()) {
            printf("-E- Could not open ofile2=\"%s\".\n", outFile2->getFileName().c_str());
            exit(EXIT_FAILURE);
        }
    }

    float centralMeridian = getCentralMeridian();

    for (row = imageHeight - 1; row >= 0; row--) {
        printPercentDone(1.0 - (float)row / (float)imageHeight);

        l3Row = l3File->getRow(row);
        numBins = l3File->getShape()->getNumCols(row);
        baseBin = l3File->getShape()->getBaseBin(row);
        endBin = baseBin + numBins;
        if (centralMeridian < 0)
            binNum = baseBin + numBins * (centralMeridian + 360.0) / 360.0;
        else
            binNum = baseBin + numBins * centralMeridian / 360.0;
        start = (imageWidth - numBins) / 2;

        // clear out beginning empty pixels
        for (col = 0; col < start; col++) {
            for (OutFile* outFile : outFiles)
                outFile->fillPixel(col);
            if (outFile2)
                outFile2->fillPixel(col);
        }

        // set pixel values
        for (int i = 0; i < numBins; i++) {
            l3Bin = l3Row->getBin(binNum);
            if (l3Bin) {
                if (doingRGB) {
                    outFiles[0]->setPixelRGB(col, l3Bin->getMean(0), l3Bin->getMean(1), l3Bin->getMean(2));
                    if (outFile2)
                        outFile2->setPixelRGB(col, l3Bin->getMean(0), l3Bin->getMean(1), l3Bin->getMean(2));
                } else {
                    for (size_t prod = 0; prod < productNameList.size(); prod++) {
                        float val;
                        switch (productMeasurementList[prod]) {
                            case AVG:
                                val = l3Bin->getMean(prod);
                                break;
                            case STDEV:
                                val = l3Bin->getStdev(prod);
                                break;
                            case VARIANCE:
                                val = l3Bin->getVariance(prod);
                                break;
                            case NOBS:
                                val = l3Bin->getNobs();
                                break;
                            case NSCENES:
                                val = l3Bin->getNscenes();
                                break;
                            case OBS_TIME:
                                val = l3Bin->getObsTime();
                                break;
                            case BIN_NUM:
                                val = l3Bin->getBinNum();
                                break;
                            default:
                                val = l3Bin->getMean(prod);
                        }
                        if (outFiles.size() == 1)
                            outFiles[0]->setPixel(col, val, prod);
                        else
                            outFiles[prod]->setPixel(col, val, 0);
                        if (outFile2)
                            outFile2->setPixel(col, val, prod);
                    }
                }
                numFilledPixels++;
            } else {
                for (OutFile* outFile : outFiles)
                    outFile->missingPixel(col);
                if (outFile2)
                    outFile2->missingPixel(col);
            }
            col++;
            binNum++;
            if (binNum >= endBin)
                binNum = baseBin;
        }

        // clear out trailing empty pixels
        for (; col < imageWidth; col++) {
            for (OutFile* outFile : outFiles)
                outFile->fillPixel(col);
            if (outFile2)
                outFile2->fillPixel(col);
        }

        for (OutFile* outFile : outFiles)
            outFile->writeLine();
        if (outFile2)
            outFile2->writeLine();
    }

    for (OutFile* outFile : outFiles) {
        outFile->setNumFilledPixels(numFilledPixels);
        outFile->close();
    }
    if (outFile2) {
        outFile2->setNumFilledPixels(numFilledPixels);
        outFile2->close();
    }
}

/**
 * @brief Writes SMI output file.
 * @param l3File Pointer to L3File.
 * @param outFiles Vector of output files.
 * @param outFile2 Pointer to the second output file.
 */
void writeSmiFile(L3File* l3File, vector<OutFile*> outFiles, OutFile* outFile2) {
    int32_t numFilledPixels = 0;
    meta_l3bType* metaData = outFiles[0]->getMetadata();
    double resolution = outFiles[0]->getResolution();

    string mapDesc = "Standard";
    string projName = "Equidistant Cylindrical";
    string tmpTitle = (string)metaData->sensor_name + " Level-3 " + mapDesc + " Mapped Image";
    strcpy(metaData->title, tmpTitle.c_str());

    // set up image parameters
    if (imageWidth <= 0) {
        imageWidth = rint(widthInDeg / 360.0 * EARTH_CIRCUMFERENCE / resolution);
        if (imageWidth == 0)
            imageWidth = 1;
    } else {
        resolution = widthInDeg / 360.0 * EARTH_CIRCUMFERENCE / imageWidth;
        for (OutFile* outFile : outFiles)
            outFile->setResolution(resolution);
        if (outFile2)
            outFile2->setResolution(resolution);
    }
    imageHeight = rint(heightInDeg / 360.0 * EARTH_CIRCUMFERENCE / resolution);
    if (imageHeight == 0)
        imageHeight = 1;
    double deltaLon = widthInDeg / imageWidth;
    double deltaLat = heightInDeg / imageHeight;

    for (OutFile* outFile : outFiles) {
        strcpy(outFile->getMetadata()->title, metaData->title);
        outFile->setSize(imageWidth, imageHeight);
        outFile->setMapProjection(projName);
    }
    if (outFile2) {
        strcpy(outFile2->getMetadata()->title, metaData->title);
        outFile2->getMetadata()->east = metaData->east;
        outFile2->getMetadata()->west = metaData->west;
        outFile2->setSize(imageWidth, imageHeight);
        outFile2->setMapProjection(projName);
    }

    // set up all the product structures
    for (size_t i = 0; i < productNameList.size(); i++) {
        OutFile* outFile;
        if (outFiles.size() == 1)
            outFile = outFiles[0];
        else
            outFile = outFiles[i];
        const ProductL3Attributes & productAttr = readProductL3Attributes(l3File, productNameList[i]);
        setupProduct(productNameList[i], productMeasurementList[i], outFile, outFile2, productAttr);
    }
    // set up quality processing
    setupQualityProcessing(l3File, outFiles, outFile2);

    printStartInfo(outFiles);
    // set cache size for writing
    if (!set_cache_size_chunk_size_write(l3File,outFiles[0],outFile2,want_verbose)) {
        fprintf(stderr, "-E-: %s:%d Could not set cache size for writing.\n",__FILE__,__LINE__);
        exit(EXIT_FAILURE);
    }
    for (OutFile* outFile : outFiles) {
        if (!outFile->open()) {
            printf("-E- Could not open ofile=\"%s\".\n", outFile->getFileName().c_str());
            exit(EXIT_FAILURE);
        }
    }

    if (outFile2) {
        if (!outFile2->open()) {
            printf("-E- Could not open ofile2=\"%s\".\n", outFile2->getFileName().c_str());
            exit(EXIT_FAILURE);
        }
    }

    InterpType interp = interpStr2Type(clo_getString(optionList, "interp"));
    float fudge = clo_getFloat(optionList, "fudge");

    // loop through output pixels
    double lat = metaData->north - (deltaLat / 2.0);
    L3Bin* l3Bin;

    for (int row = 0; row < imageHeight; row++) {
        printPercentDone((float)row / (float)imageHeight);
        double lon = metaData->west + (deltaLon / 2.0);

        for (int col = 0; col < imageWidth; col++) {
            if (applyMask && isLand(lat, lon)) {
                for (OutFile* outFile : outFiles)
                    outFile->landPixel(col);
                if (outFile2)
                    outFile2->landPixel(col);

            } else {
                switch (interp) {
                    case INTERP_NEAREST:
                        l3Bin = l3File->getClosestBin(lat, lon);
                        break;
                    case INTERP_BIN:
                    case INTERP_AREA: {
                        bool areaWeighted;

                        if (interp == INTERP_AREA)
                            areaWeighted = true;
                        else
                            areaWeighted = false;

                        l3Bin = getBoxBins(l3File, lat, lon, deltaLat, deltaLon, fudge, areaWeighted);

                        break;
                    }
                    default:
                        printf("-E- interp = %s is not implemented.", interpType2Str(interp));
                        exit(EXIT_FAILURE);
                }

                if (l3Bin) {
                    numFilledPixels++;
                    if (doingRGB) {
                        outFiles[0]->setPixelRGB(col, l3Bin->getMean(0), l3Bin->getMean(1),
                                                 l3Bin->getMean(2));
                        if (outFile2)
                            outFile2->setPixelRGB(col, l3Bin->getMean(0), l3Bin->getMean(1),
                                                  l3Bin->getMean(2));
                    } else {
                        for (size_t prod = 0; prod < productNameList.size(); prod++) {
                            float val;
                            switch (productMeasurementList[prod]) {
                                case AVG:
                                    val = l3Bin->getMean(prod);
                                    break;
                                case STDEV:
                                    val = l3Bin->getStdev(prod);
                                    break;
                                case VARIANCE:
                                    val = l3Bin->getVariance(prod);
                                    break;
                                case NOBS:
                                    val = l3Bin->getNobs();
                                    break;
                                case NSCENES:
                                    val = l3Bin->getNscenes();
                                    break;
                                case OBS_TIME:
                                    val = l3Bin->getObsTime();
                                    break;
                                case BIN_NUM:
                                    val = l3Bin->getBinNum();
                                    break;
                                default:
                                    val = l3Bin->getMean(prod);
                            }
                            if (outFiles.size() == 1)
                                outFiles[0]->setPixel(col, val, prod);
                            else
                                outFiles[prod]->setPixel(col, val, 0);
                            if (outFile2)
                                outFile2->setPixel(col, val, prod);
                        }
                    }
                    if (doingQuality) {
                        for (OutFile* outFile : outFiles)
                            outFile->setQuality(col, l3Bin->getQuality());
                        if (outFile2)
                            outFile2->setQuality(col, l3Bin->getQuality());
                    }
                } else {
                    for (OutFile* outFile : outFiles)
                        outFile->missingPixel(col);
                    if (outFile2)
                        outFile2->missingPixel(col);
                }
            }
            lon += deltaLon;
        }  // for col
        for (OutFile* outFile : outFiles)
            outFile->writeLine();
        if (outFile2)
            outFile2->writeLine();
        lat -= deltaLat;
    }  // for row

    for (OutFile* outFile : outFiles) {
        outFile->setNumFilledPixels(numFilledPixels);
        outFile->close();
    }
    if (outFile2) {
        outFile2->setNumFilledPixels(numFilledPixels);
        outFile2->close();
    }
}

// writing netcdf 4 file
/**
 * @brief Writes projected output file using PROJ.4 library (the writitng is done here).
 * @param l3File Pointer to L3File.
 * @param projectionStr Projection string.
 * @param outFiles Vector of output files.
 * @param outFile2 Pointer to the second output file.
 * @param trimNSEW Flag to trim North, South, East, West.
 */
void writeProj4File(L3File* l3File, char* projectionStr, vector<OutFile*> outFiles, OutFile* outFile2,
                    bool trimNSEW) {
    int32_t numFilledPixels = 0;
    meta_l3bType* metaData = outFiles[0]->getMetadata();
    double resolution = outFiles[0]->getResolution();

    L3Bin* l3Bin;

    // how about circumference of earth + 25%
    double limitMin = EARTH_CIRCUMFERENCE * -0.625;
    double limitMax = EARTH_CIRCUMFERENCE * 0.625;
    // calculate the min and max
    double minX = limitMax;
    double minY = limitMax;
    double maxX = limitMin;
    double maxY = limitMin;
    double lat, lon;
    double x, y;
    double deltaLat = (metaData->north - metaData->south) / (NUM_SEARCH_POINTS - 1);
    double deltaLon = (mapEast - mapWest) / (NUM_SEARCH_POINTS - 1);

    // parse central meridian (lon_0)
    float centralMeridian = getCentralMeridian();
    string cmStr = " +lon_0=" + to_string(centralMeridian);
    string tsStr;
    if (clo_isSet(optionList, "lat_ts")) {
        tsStr = " +lat_ts=" + to_string(clo_getFloat(optionList, "lat_ts"));
    }

    float lat0 = (metaData->north + metaData->south) / 2.0;
    string lat0Str;
    if (clo_isSet(optionList, "lat_0")) {
        lat0 = clo_getFloat(optionList, "lat_0");
    }
    lat0Str = " +lat_0=" + to_string(lat0);

    string lat1Str;
    float lat1Val;
    if (clo_isSet(optionList, "lat_1")) {
        lat1Val = clo_getFloat(optionList, "lat_1");
    } else {
        lat1Val = metaData->south;
    }

    // check lat1 for -90 or 90
    if((strcasecmp(projectionStr, "lambert") == 0)) {
        if(lat1Val >= 90)
            lat1Val = 89.999;
        if(lat1Val <= -90)
            lat1Val = -89.999;
    }
    lat1Str = " +lat_1=" + to_string(lat1Val);

    string lat2Str;
    float lat2Val;
    if (clo_isSet(optionList, "lat_2")) {
        lat2Val = clo_getFloat(optionList, "lat_2");
    } else {
        lat2Val = metaData->north;
    }

    // check lat2 for -90 or 90
    if((strcasecmp(projectionStr, "lambert") == 0)) {
        if(lat2Val >= 90)
            lat2Val = 89.999;
        if(lat2Val <= -90)
            lat2Val = -89.999;
    }
    lat2Str = " +lat_2=" + to_string(lat2Val);

    string aziStr;
    if (clo_isSet(optionList, "azimuth")) {
        aziStr = " +alpha=" + to_string(clo_getFloat(optionList, "azimuth"));
    }
    string utmStr;
    if (clo_isSet(optionList, "utm_zone")) {
        string zone = clo_getString(optionList, "utm_zone");
        utmStr = " +zone=";
        utmStr += zone;
        if (utmStr.find('S') != string::npos) {
            utmStr.pop_back();
            utmStr += " +south";
        }
    }

    // define proj.4 parameters according to shortcut
    string mapDesc;
    string projName;
    string projStr;

    if (strcasecmp(projectionStr, "mollweide") == 0) {
        mapDesc = "Mollweide";
        projName = "Mollweide";
        projStr = "+proj=moll +ellps=WGS84 +datum=WGS84";
        projStr += cmStr;

    } else if (strcasecmp(projectionStr, "lambert") == 0) {
        mapDesc = "Lambert";
        projName = "Lambert";
        projStr = "+proj=lcc +ellps=WGS84 +datum=WGS84";
        projStr += cmStr + lat0Str + lat1Str + lat2Str;

    } else if (strcasecmp(projectionStr, "albersconic") == 0) {
        mapDesc = "Albers Equal Area Conic";
        projName = "Albersconic";
        projStr = "+proj=aea +ellps=WGS84 +datum=WGS84";
        projStr += cmStr + lat0Str + lat1Str + lat2Str;

    } else if (strcasecmp(projectionStr, "aeqd") == 0) {
        mapDesc = "Azimuthal Equidistant";
        projName = "AzimuthalEquidistant";
        projStr = "+proj=aeqd +ellps=WGS84 +datum=WGS84";
        projStr += cmStr + lat0Str;

    } else if (strcasecmp(projectionStr, "mercator") == 0) {
        mapDesc = "Mercator";
        projName = "Mercator";
        projStr = "+proj=merc +ellps=WGS84 +datum=WGS84";
        projStr += cmStr;

    } else if (strcasecmp(projectionStr, "tmerc") == 0) {
        mapDesc = "Transverse Mercator";
        projName = "TransverseMercator";
        projStr = "+proj=tmerc +ellps=WGS84 +datum=WGS84";
        projStr += cmStr + lat0Str;

    } else if (strcasecmp(projectionStr, "utm") == 0) {
        mapDesc = "Universal Transverse Mercator";
        projName = "UTM";
        projStr = "+proj=utm";
        projStr += utmStr;

    } else if (strcasecmp(projectionStr, "obliquemerc") == 0) {
        if (!clo_isSet(optionList, "lat_0") || !clo_isSet(optionList, "azimuth")) {
            printf("-E- lat_0 and azimuth need to be defined for obliquemerc projection");
            exit(EXIT_FAILURE);
        }
        mapDesc = "Oblique Mercator";
        projName = "ObliqueMercator";
        projStr = "+proj=omerc +gamma=0 +k_0=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84";
        projStr += lat0Str;
        projStr += " +lonc=" + to_string(centralMeridian);
        projStr += aziStr;

    } else if (strcasecmp(projectionStr, "ease2") == 0) {
        mapDesc = "Ease Grid 2";
        projName = "Ease2";
        projStr = "EPSG:6933";

    } else if (strcasecmp(projectionStr, "stere") == 0) {
        if (!clo_isSet(optionList, "lat_0") || !clo_isSet(optionList, "lat_ts")) {
            printf("-E- lat_0 and lat_ts need to be defined for stere projection");
            exit(EXIT_FAILURE);
        }
        mapDesc = "Stereographic";
        projName = "Stereo";
        projStr =
            "+proj=stere "
            " +ellps=WGS84 +datum=WGS84 +units=m +no_defs";
        projStr += cmStr + tsStr + lat0Str;

    } else if (strcasecmp(projectionStr, "ortho") == 0) {
        mapDesc = "Orthographic";
        projName = "Ortho";
        projStr =
            "+proj=ortho "
            " +ellps=GRS80 +units=m +no_defs";
        projStr += cmStr + lat0Str;

    } else if (strcasecmp(projectionStr, "conus") == 0) {
        mapDesc = "USA Contiguous Albers Equal Area Conic USGS version";
        projName = "Conus";
        projStr =
            "+proj=aea +lat_1=29.5 +lat_2=45.5"
            " +lat_0=23.0 +lon_0=-96 +x_0=0 +y_0=0"
            " +ellps=GRS80 +datum=NAD83 +units=m +no_defs";

    } else if (strcasecmp(projectionStr, "alaska") == 0) {
        mapDesc = "Alaska Albers Equal Area Conic USGS version";
        projName = "Alaska";
        projStr = "EPSG:3338";

    } else if (strcasecmp(projectionStr, "gibs") == 0) {
        if (((metaData->north + metaData->south) / 2.) > 60.) {
            mapDesc = "Stereographic";
            projName = "GIBS Stereo";
            projStr = "EPSG:3413";

        } else if (((metaData->north + metaData->south) / 2.) < -60.) {
            mapDesc = "Stereographic";
            projName = "GIBS Stereo";
            projStr = "EPSG:3031";

        } else {
            mapDesc = "Equidistant Cylindrical";
            projName = "PlateCarree";
            projStr =
                "+proj=eqc"
                " +lat_ts=0 +lat_0=0 +x_0=0 +y_0=0"
                " +ellps=WGS84 +datum=WGS84 +units=m +no_defs";
            projStr += cmStr;
        }

    } else if (strcasecmp(projectionStr, "platecarree") == 0) {
        mapDesc = "Equidistant Cylindrical";
        projName = "PlateCarree";
        projStr =
            "+proj=eqc"
            " +lat_ts=0 +lat_0=0 +x_0=0 +y_0=0"
            " +ellps=WGS84 +datum=WGS84 +units=m +no_defs";
        projStr += cmStr;

    } else {
        mapDesc = "Proj4";
        projName = projectionStr;
        projStr = projectionStr;
    }

    // save parameters in metadata
    string tmpTitle = (string)metaData->sensor_name + " Level-3 " + mapDesc + " Mapped Image";
    strcpy(metaData->title, tmpTitle.c_str());

    for (OutFile* outFile : outFiles) {
        if (outFile->getMetadata() != metaData)
            strcpy(outFile->getMetadata()->title, metaData->title);
        outFile->setMapProjection(projName);
    }
    if (outFile2) {
        strcpy(outFile2->getMetadata()->title, metaData->title);
        outFile2->getMetadata()->east = metaData->east;
        outFile2->getMetadata()->west = metaData->west;
        outFile2->setMapProjection(projName);
    }

    PJ *pj, *pj_new;
    PJ_COORD c, c_out;

    // init the proj4 projections
    pj = proj_create_crs_to_crs(PJ_DEFAULT_CTX, "EPSG:4326", projStr.c_str(), NULL);
    if (pj == NULL) {
        printf("Error - l3mapgen first PROJ projection failed to init\n");
        exit(1);
    }
    pj_new = proj_normalize_for_visualization(PJ_DEFAULT_CTX, pj);
    if (pj_new == NULL) {
        printf("Error - l3mapgen visualization PROJ projection failed to init\n");
        exit(1);
    }
    proj_destroy(pj);
    pj = NULL;

    // calculate start and delta for grid
    // set default z and t
    c.xyzt.z = 0.0;
    c.xyzt.t = 0.0;

    // define "inBox" region for trimNSEW option using lat/lon arrays based on
    // metaData struct
    array<double, 9> lats{{metaData->south, metaData->south, metaData->south,
                                (metaData->north + metaData->south) / 2., metaData->north, metaData->north,
                                metaData->north, (metaData->north + metaData->south) / 2., metaData->south}};
    array<double, 9> lons{{metaData->west, (metaData->east + metaData->west) / 2., metaData->east,
                                metaData->east, metaData->east, (metaData->east + metaData->west) / 2.,
                                metaData->west, metaData->west, metaData->west}};
    vector<Point_t> points;
    for (int32_t i = 0; i < 9; i++) {
        c.xy.x = lons[i];
        c.xy.y = lats[i];
        c_out = proj_trans(pj_new, PJ_FWD, c);

        if (isfinite(c_out.xy.x) && isfinite(c_out.xy.y)) {
            Point_t point(c_out.xy.x, c_out.xy.y);
            points.push_back(point);
        }
    }

    // Create a polygon object and assign the points to it.
    Polygon_t NSEW;
    boost::geometry::assign_points(NSEW, points);

    if (clo_isSet(optionList, "north")) {
        // ...if user defined the region
        lat = metaData->south;
        for (int j = 0; j < NUM_SEARCH_POINTS; j++) {
            lon = metaData->west;
            for (int i = 0; i < NUM_SEARCH_POINTS; i++) {
                c.xy.x = lon;
                c.xy.y = lat;
                c_out = proj_trans(pj_new, PJ_FWD, c);
                x = c_out.xy.x;
                y = c_out.xy.y;
                if (isfinite(x) && isfinite(y)) {
                    if (x < limitMax && x > limitMin) {
                        if (x < minX)
                            minX = x;
                        if (x > maxX)
                            maxX = x;
                    }
                    if (y < limitMax && y > limitMin) {
                        if (y < minY)
                            minY = y;
                        if (y > maxY)
                            maxY = y;
                    }
                }
                lon += deltaLon;
            }
            lat += deltaLat;
        }
    } else {
        // If no user input, grab extents from the file -
        // not taking the easy route and using metadata because of the square
        // Earth problem

        L3Shape* shape = l3File->getShape();
        for (int row = 0; row < l3File->getNumRows(); row++) {
            L3Row* l3row = l3File->getRow(row);
            if(l3row->getNumBins() == 0)
                continue;
            x = std::numeric_limits<double>::quiet_NaN();
            y = std::numeric_limits<double>::quiet_NaN();
            int32_t binindex = 0;

            // if dateline crossed...search every bin
            if (metaData->east < metaData->west) {
                while (binindex < l3row->getNumBins()) {
                    int64_t bin = l3row->getBinByIndex(binindex)->getBinNum();
                    if (bin > -1) {
                        shape->bin2latlon(bin, lat, lon);
                        c.xy.x = lon;
                        c.xy.y = lat;
                        c_out = proj_trans(pj_new, PJ_FWD, c);
                        x = c_out.xy.x;
                        y = c_out.xy.y;
                        binindex++;
                        if (isfinite(x) && isfinite(y)) {
                            if (x < limitMax && x > limitMin) {
                                if (x < minX) {
                                    minX = x;
                                }
                                if (x > maxX) {
                                    maxX = x;
                                }
                            }
                            if (y < limitMax && y > limitMin) {
                                if (y < minY) {
                                    minY = y;
                                }
                                if (y > maxY) {
                                    maxY = y;
                                }
                            }
                        }
                    }
                }
            } else {
                //First column
                while ((!isfinite(x) || !(isfinite(y))) && binindex < l3row->getNumBins()) {
                    int64_t bin = l3row->getBinByIndex(binindex)->getBinNum();
                    if (bin > -1) {
                        shape->bin2latlon(bin, lat, lon);
                        c.xy.x = lon;
                        c.xy.y = lat;
                        c_out = proj_trans(pj_new, PJ_FWD, c);
                        x = c_out.xy.x;
                        y = c_out.xy.y;
                        binindex++;
                    }
                }
                if (isfinite(x) && isfinite(y)) {
                    if (x < limitMax && x > limitMin) {
                        if (x < minX) {
                            minX = x;
                        }
                        if (x > maxX) {
                            maxX = x;
                        }
                    }
                    if (y < limitMax && y > limitMin) {
                        if (y < minY) {
                            minY = y;
                        }
                        if (y > maxY) {
                            maxY = y;
                        }
                    }
                }
                // Last column
                if ((l3row->getNumBins() > 1) && (binindex != l3row->getNumBins())) {
                    binindex = l3row->getNumBins() - 1;
                    x = std::numeric_limits<double>::quiet_NaN();
                    y = std::numeric_limits<double>::quiet_NaN();

                    // search for last good bin
                    while ((!isfinite(x) || !(isfinite(y))) && (binindex > 0)) {
                        int64_t bin = l3row->getBinByIndex(binindex)->getBinNum();
                        if (bin > -1) {
                            shape->bin2latlon(bin, lat, lon);
                            c.xy.x = lon;
                            c.xy.y = lat;
                            c_out = proj_trans(pj_new, PJ_FWD, c);
                            x = c_out.xy.x;
                            y = c_out.xy.y;
                            binindex--;
                        }
                    }
                    if (isfinite(x) && isfinite(y)) {
                        if (x < limitMax && x > limitMin) {
                            if (x < minX) {
                                minX = x;
                            }
                            if (x > maxX) {
                                maxX = x;
                            }
                        }
                        if (y < limitMax && y > limitMin) {
                            if (y < minY) {
                                minY = y;
                            }
                            if (y > maxY) {
                                maxY = y;
                            }
                        }
                    }
                }
           }
        }

        // add a little padding to get last pixel
        minX -= resolution;
        maxX += resolution;
        minY -= resolution;
        maxY += resolution;

    }  // "north" not set

    // set up image parameters
    if (imageWidth <= 0) {
        imageWidth = rint((maxX - minX) / resolution);
        if (imageWidth <= 0)
            imageWidth = 1;
    } else {
        resolution = (maxX - minX) / imageWidth;
        for (OutFile* outFile : outFiles)
            outFile->setResolution(resolution);
        if (outFile2)
            outFile2->setResolution(resolution);
    }
    imageHeight = rint((maxY - minY) / resolution);
    if (imageHeight <= 0)
        imageHeight = 1;
    for (OutFile* outFile : outFiles)
        outFile->setSize(imageWidth, imageHeight);
    if (outFile2)
        outFile2->setSize(imageWidth, imageHeight);

    for (OutFile* outFile : outFiles) {
        outFile->setProj4Info(projStr, minX, maxY);
    }

    if (outFile2) {
        outFile2->setProj4Info(projStr, minX, maxY);
    }

    // set up all the product structures
    for (size_t i = 0; i < productNameList.size(); i++) {
        OutFile* outFile;
        if (outFiles.size() == 1)
            outFile = outFiles[0];
        else
            outFile = outFiles[i];
        const ProductL3Attributes & productAttr = readProductL3Attributes(l3File, productNameList[i]);
        setupProduct(productNameList[i], productMeasurementList[i], outFile, outFile2, productAttr);
    }

    // set up quality processing
    setupQualityProcessing(l3File, outFiles, outFile2);

    printStartInfo(outFiles);
    // set cache size for writing
    if (!set_cache_size_chunk_size_write(l3File,outFiles[0],outFile2,want_verbose)) {
        fprintf(stderr, "-E-: %s:%d Could not set cache size for writing.\n",__FILE__,__LINE__);
        exit(EXIT_FAILURE);
    }
    for (OutFile* outFile : outFiles) {
        if (!outFile->open()) {
            printf("-E- Could not open ofile=\"%s\".\n", outFile->getFileName().c_str());
            exit(EXIT_FAILURE);
        }
    }

    if (outFile2) {
        if (!outFile2->open()) {
            printf("-E- Could not open ofile2=\"%s\".\n", outFile2->getFileName().c_str());
            exit(EXIT_FAILURE);
        }
    }

    InterpType interp = interpStr2Type(clo_getString(optionList, "interp"));
    float fudge = clo_getFloat(optionList, "fudge");

    // loop through output pixels
    vector<double> tmpX(imageWidth);
    vector<double> tmpY(imageWidth);
    vector<bool> inBox(imageWidth, false);

    deltaLon = widthInDeg / imageWidth;
    deltaLat = heightInDeg / imageHeight;
    double startX = minX + resolution / 2;
    double startY = maxY - resolution / 2;
    y = startY;

    // setting values begins
    for (int row = 0; row < imageHeight; row++) {
        printPercentDone((float)row / (float)imageHeight);
        x = startX;
        fill(inBox.begin(), inBox.end(), false);
        for (int col = 0; col < imageWidth; col++) {
            c.xy.x = x;
            c.xy.y = y;
            c_out = proj_trans(pj_new, PJ_INV, c);
            tmpX[col] = c_out.xy.x;
            tmpY[col] = c_out.xy.y;
            if (!isfinite(tmpX[col]) || !isfinite(tmpY[col])) {
                tmpX[col] = -999.0;  // set to fill value
                tmpY[col] = -999.0;  // set to fill value

            } else if (trimNSEW) {
                Point_t pixel(x, y);
                if (boost::geometry::within(pixel, NSEW)) {
                    inBox[col] = true;
                }
            }
            x += resolution;
        }
        // loop through each pixel in this row
        for (int col = 0; col < imageWidth; col++) {
            lon = tmpX[col];
            lat = tmpY[col];
            if (!isfinite(lon) || !isfinite(lat)) {
                for (OutFile* outFile : outFiles)
                    outFile->fillPixel(col);
                if (outFile2)
                    outFile2->fillPixel(col);
            } else if (trimNSEW && inBox[col] == false) {
                for (OutFile* outFile : outFiles)
                    outFile->fillPixel(col);
                if (outFile2)
                    outFile2->fillPixel(col);
            } else if (applyMask && isLand(lat, lon)) {
                for (OutFile* outFile : outFiles)
                    outFile->landPixel(col);
                if (outFile2)
                    outFile2->landPixel(col);

            } else {
                switch (interp) {
                    case INTERP_NEAREST:
                        l3Bin = l3File->getClosestBin(lat, lon);
                        break;
                    case INTERP_BIN:
                    case INTERP_AREA: {
                        bool areaWeighted;

                        if (interp == INTERP_AREA)
                            areaWeighted = true;
                        else
                            areaWeighted = false;

                        if (checkDateLineCrossed(lon, deltaLon)) {
                            // East
                            l3Bin =
                                getBoxBins(l3File, lat, lon, deltaLat, deltaLon, fudge, areaWeighted, IS_EAST);

                            // West
                            L3Bin* tmpBin =
                                getBoxBins(l3File, lat, lon, deltaLat, deltaLon, fudge, areaWeighted, IS_WEST);
                            if (tmpBin) {
                                *l3Bin += *tmpBin;
                            }

                        } else {
                            l3Bin = getBoxBins(l3File, lat, lon, deltaLat, deltaLon, fudge, areaWeighted);
                        }
                        break;
                    }
                    default:
                        printf("-E- interp = %s is not implemented.", interpType2Str(interp));
                        exit(EXIT_FAILURE);
                }

                if (l3Bin) {
                    numFilledPixels++;
                    if (doingRGB) {
                        outFiles[0]->setPixelRGB(col, l3Bin->getMean(0), l3Bin->getMean(1),
                                                 l3Bin->getMean(2));
                        if (outFile2)
                            outFile2->setPixelRGB(col, l3Bin->getMean(0), l3Bin->getMean(1),
                                                  l3Bin->getMean(2));
                    } else {
                        // it loops through product index
                        for (size_t prod = 0; prod < productNameList.size(); prod++) {
                            float val;
                            switch (productMeasurementList[prod]) {
                                case AVG:
                                    val = l3Bin->getMean(prod);
                                    break;
                                case STDEV:
                                    val = l3Bin->getStdev(prod);
                                    break;
                                case VARIANCE:
                                    val = l3Bin->getVariance(prod);
                                    break;
                                case NOBS:
                                    val = l3Bin->getNobs();
                                    break;
                                case NSCENES:
                                    val = l3Bin->getNscenes();
                                    break;
                                case OBS_TIME:
                                    val = l3Bin->getObsTime();
                                    break;
                                case BIN_NUM:
                                    val = l3Bin->getBinNum();
                                    break;
                                default:
                                    val = l3Bin->getMean(prod);
                            }
                            if (outFiles.size() == 1)
                                outFiles[0]->setPixel(col, val, prod);
                            else
                                outFiles[prod]->setPixel(col, val, 0);
                            if (outFile2)
                                outFile2->setPixel(col, val, prod);
                        }
                    }
                    if (doingQuality) {
                        for (OutFile* outFile : outFiles)
                            outFile->setQuality(col, l3Bin->getQuality());
                        if (outFile2)
                            outFile2->setQuality(col, l3Bin->getQuality());
                    }
                } else {
                    for (OutFile* outFile : outFiles)
                        outFile->missingPixel(col);
                    if (outFile2)
                        outFile2->missingPixel(col);
                }
            }
        }  // for col
        for (OutFile* outFile : outFiles) {
            outFile->setLatLon(tmpY.data(), tmpX.data());
            outFile->writeLine();
        }
        if (outFile2) {
            outFile2->setLatLon(tmpY.data(), tmpX.data());
            outFile2->writeLine();
        }
        y -= resolution;
    }  // for row
       // values set ends here

    proj_destroy(pj_new);

    for (OutFile* outFile : outFiles) {
        if (writeProjectionText) {
            string projtxtfilename;
            projtxtfilename = outFile->getFileName();
            projtxtfilename += ".projtxt";
            ofstream projtxtfile(projtxtfilename);
            if (projtxtfile.is_open()) {
                projtxtfile << "# Projection information for " << outFile->getFileName() << "\n";
                projtxtfile << "proj=" << projStr << "\n";
                projtxtfile << "minX=" << setprecision(11) << minX << "\n";
                projtxtfile << "maxX=" << setprecision(11) << maxX << "\n";
                projtxtfile << "minY=" << setprecision(11) << minY << "\n";
                projtxtfile << "maxY=" << setprecision(11) << maxY << "\n";
                projtxtfile << "north=" << setprecision(11) << metaData->north << "\n";
                projtxtfile << "south=" << setprecision(11) << metaData->south << "\n";
                projtxtfile << "east=" << setprecision(11) << metaData->east << "\n";
                projtxtfile << "west=" << setprecision(11) << metaData->west << "\n";

                projtxtfile << "scale_type=" << outFile->getScaleTypeString() << "\n";
                projtxtfile << "datamin=" << setprecision(11) << outFile->getMinValue() << "\n";
                projtxtfile << "datamax=" << setprecision(11) << outFile->getMaxValue() << "\n";
                projtxtfile << "width=" << imageWidth << "\n";
                projtxtfile << "height=" << imageHeight << "\n";

                projtxtfile.close();
            }
        }

        outFile->setNumFilledPixels(numFilledPixels);
        outFile->close();
    }
    if (outFile2) {
        outFile2->setNumFilledPixels(numFilledPixels);
        outFile2->close();
    }
}

/**
 * @brief Creates an output file object based on format.
 * @param oformatStr Output format string.
 * @param useColor Flag to use color.
 * @return Pointer to the created output file object.
 */
OutFile* makeOutputFile(const char* oformatStr, bool useColor) {
    OutFile* outFile = NULL;

    const char* oformatStr2 = getFileFormatName(oformatStr);

    if (oformatStr2 == NULL) {
        printf("-E- Unknown output file format \"%s\"\n", oformatStr);
        exit(EXIT_FAILURE);
    }

    string oformat = oformatStr2;
    if (oformat.compare("PPM") == 0) {
        if (doingRGB) {
            outFile = new OutFilePpmRgb();
        } else {
            if (useColor) {
                outFile = new OutFilePpm();
            } else {
                outFile = new OutFilePgm();
            }
        }
    } else if (oformat.compare("PNG") == 0) {
        if (doingRGB) {
            outFile = new OutFilePngRgb();
        } else {
            outFile = new OutFilePng(useColor);
        }
    } else if (oformat.compare("TIFF") == 0) {
        if (doingRGB) {
            outFile = new OutFileTiffRgb();
        } else {
            if (useColor) {
                outFile = new OutFileTiffColor();
            } else {
                outFile = new OutFileTiffGray();
            }
        }
    } else if (oformat.compare("HDF4") == 0) {
        outFile = new OutFileHdf4();
    } else if (oformat.compare("netCDF4") == 0) {
        outFile = new OutFileNetcdf4();
    } else {
        printf("-E- Output file type %s not implemented\n", oformat.c_str());
        exit(EXIT_FAILURE);
    }
    if (doingTransparency)
        outFile->setTransparency();

    return outFile;
}

/**
 * @brief Creates a list of output files based on the options.
 * @param optionList Pointer to the option list.
 * @return Vector of pointers to the created output files.
 */
vector<OutFile*> makeOutputFileList(clo_optionList_t* optionList) {
    vector<OutFile*> outFiles;
    OutFile* outFile;
    string oformatStr = getFileFormatName(clo_getString(optionList, "oformat"));
    string originalOfile = clo_getString(optionList, "ofile");
    string tag = clo_getString(optionList, "ofile_product_tag");
    // if netcdf4 format specified just one file will be produced
    // need to put a check : if there 3D vars, then only netcdf is acceptable
    if (oformatStr.compare("netCDF4") == 0 || productNameList.size() == 1 || doingRGB) {
        outFile = makeOutputFile(clo_getString(optionList, "oformat"), clo_getBool(optionList, "apply_pal"));
        outFile->setFileName(originalOfile);
        outFiles.push_back(outFile);
    } else {
        size_t pos = originalOfile.find(tag);
        if (pos == string::npos) {
            printf("Error: ofile_product_tag=%s, not found in ofile=%s\n", tag.c_str(),
                   originalOfile.c_str());
            printf(
                "       and you asked for multiple products with image "
                "oformat\n");
            exit(EXIT_FAILURE);
        }

        for (size_t i = 0; i < productNameList.size(); i++) {
            string newName = originalOfile;
            string prodName_clean = productNameList.at(i);
            string modifier = getModifierFromMeasure(productMeasurementList.at(i));
            newName.replace(pos, tag.size(), prodName_clean + modifier);
            outFile =
                makeOutputFile(clo_getString(optionList, "oformat"), clo_getBool(optionList, "apply_pal"));
            outFile->setFileName(newName);
            outFiles.push_back(outFile);
        }  // loop products
    }

    return outFiles;
}

//------------------------------------------------------------------------------
// main
//------------------------------------------------------------------------------

/**
 * @brief Main function for the l3mapgen program.
 * @param argc Argument count.
 * @param argv Argument vector.
 * @return Exit status.
 */
int main(int argc, char* argv[]) {
    vector<OutFile*> outFiles;
    OutFile* outFile2 = NULL;
    char* ifileName;
    string oformat;
    int i;
    char* tmpStr;

    string softwareVersion = to_string(VERSION_MAJOR) + "." + to_string(VERSION_MINOR) + "." + to_string(VERSION_PATCH) + "-" + GITSHA;

    optionList = clo_createList();

    l3mapgenInitOptions(optionList, softwareVersion.c_str());
    if (argc == 1) {
        clo_printUsage(optionList);
        exit(EXIT_FAILURE);
    }
    l3mapgenReadOptions(optionList, argc, argv);

    if (clo_getBool(optionList, "quiet")) {
        want_verbose = 0;
    }

    ifileName = clo_getString(optionList, "ifile");
    // set proper cache size for reading
     if (!set_cache_size_chunk_size_read(ifileName, want_verbose)) {
         fprintf(stderr, "-E-: %s:%d Could not set cache size for reading.\n",__FILE__,__LINE__);
         exit(EXIT_FAILURE);
     }
    // try SMI input file
    int oldVerbose = want_verbose;
    want_verbose = 0;
    L3File* l3File = new L3FileSMI();
    if (!l3File->open(ifileName)) {
        delete l3File; // open failed, delete current object
        
        // try real L3 bin format
        l3File = new L3File();
        if (!l3File->open(ifileName)) {
            printf("-E- Could not open ifile=\"%s\".\n", ifileName);
            exit(EXIT_FAILURE);
        }
    }
    want_verbose = oldVerbose;

    if (clo_getBool(optionList, "use_transparency")) {
        doingTransparency = true;
    }
    if (clo_getBool(optionList, "use_rgb")) {
        doingRGB = true;
        productName = clo_getRawString(optionList, "product_rgb");
    } else {
        doingRGB = false;
        if (clo_isSet(optionList, "product")) {
            productName = clo_getRawString(optionList, "product");
        } else {
            productName.clear();
            for (int i = 0; i < l3File->getNumProducts(); i++) {
                string tmpName = l3File->getProductName(i);
                if (tmpName != "qual_l3") {
                    if (!productName.empty())
                        productName += ",";
                    productName += tmpName;
                }
            }
        }
    }
    trimNSEW = clo_getBool(optionList, "trimNSEW");
    writeProjectionText = clo_getBool(optionList, "write_projtext");

    // read user specified product names and measurements, read wavelengths if specified, and expand product names with the wavelength specifier if applicable
    getProductNames(productName,productNameList,l3File, optionList);
    std::string cleanProdName;
    vector<string> parts;
    // setup measurements, remove measurement type from product name if specified, and create clean product name list for looking up in file
    for (size_t i = 0; i < productNameList.size(); i++) {
        if (i != 0)
            cleanProdName += ",";
        boost::split(parts, productNameList[i], boost::is_any_of(":"));
        if (parts.size() == 1) {
            cleanProdName += parts[0];
            productMeasurementList.push_back(AVG);
        } else if (parts.size() == 2) {
            productNameList[i] = parts[0];  // get rid of the modifier
            cleanProdName += parts[0];
            if (parts[1].compare("avg") == 0)
                productMeasurementList.push_back(AVG);
            else if (parts[1].compare("stdev") == 0)
                productMeasurementList.push_back(STDEV);
            else if (parts[1].compare("var") == 0)
                productMeasurementList.push_back(VARIANCE);
            else if (parts[1].compare("nobs") == 0)
                productMeasurementList.push_back(NOBS);
            else if (parts[1].compare("nscenes") == 0)
                productMeasurementList.push_back(NSCENES);
            else if (parts[1].compare("obs_time") == 0)
                productMeasurementList.push_back(OBS_TIME);
            else if (parts[1].compare("bin_num") == 0)
                productMeasurementList.push_back(BIN_NUM);
            else {
                EXIT_LOG(
                    printf("-E- measurement type \"%s\" "
                           "not understood for product \"%s\".\n",
                           parts[1].c_str(), parts[0].c_str()));
            }
        } else {
            EXIT_LOG(printf("-E- product name not understood \"%s\".\n", productNameList[i].c_str()));
        }
    }
    if (!l3File->setActiveProductList(cleanProdName.c_str())) {
        EXIT_LOG(
            printf("-E- Could not find product=\"%s\" in file=\"%s\".\n", cleanProdName.c_str(), ifileName));
    }

    l3File->setNumCacheRows(clo_getInt(optionList, "num_cache"));

    // copy the L3 meta data since we will modify it a bit for the output file.

    // output file list is here
    // should be modified first
    outFiles = makeOutputFileList(optionList);

    if (clo_isSet(optionList, "ofile2")) {
        outFile2 =
            makeOutputFile(clo_getString(optionList, "oformat2"), clo_getBool(optionList, "apply_pal"));
        outFile2->setFileName(clo_getString(optionList, "ofile2"));
    }

    // resolution = # of meters across 1 pixel in the center of the scene
    double res;
    if (clo_isSet(optionList, "resolution")) {
        outFiles[0]->setResolution(clo_getString(optionList, "resolution"));
        res = outFiles[0]->getResolution();
    } else {
        res = l3File->getMetaData()->resolution;  // in meters
        if (res < 0) {
            outFiles[0]->setResolution("9km");
            res = outFiles[0]->getResolution();
        }
    }

    for (OutFile* outFile : outFiles) {
        outFile->setResolution(res);
    }

    if (outFile2)
        outFile2->setResolution(res);

    // get imageWidth, if specified
    if (clo_isSet(optionList, "width"))
        imageWidth = clo_getInt(optionList, "width");

    // projection
    clo_option_t* projectionOption = clo_findOption(optionList, "projection");
    char* projectionStr = clo_getOptionRawString(projectionOption);
    meta_l3bType metaData = *l3File->getMetaData();
    // check the metadata of the bin file
    if (metaData.north == metaData.south) {
        printf("-E- north and south metadata are equal.\n");
        exit(110);
    }
    if (metaData.east == metaData.west) {
        printf("-E- east and west metadata are equal.\n");
        exit(110);
    }

    // default to whole globe for SMI files
    if ((strcmp(projectionStr, "smi") == 0) || (strcmp(projectionStr, "raw") == 0)) {
        metaData.north = 90.0;
        metaData.south = -90.0;
        metaData.east = 180.0;
        metaData.west = -180.0;
    }
    // read in north, south, east, west from command line
    float north = clo_getFloat(optionList, "north");
    int32_t nsew = 0;
    if (north > -90.0 && north <= 90.0) {
        metaData.north = north;
        nsew++;
    }
    float south = clo_getFloat(optionList, "south");
    if (south < 90.0 && south >= -90) {
        metaData.south = south;
        nsew++;
    }
    float east = clo_getFloat(optionList, "east");
    if (east > -180.0 && east <= 180.0) {
        metaData.east = east;
        nsew++;
    }
    float west = clo_getFloat(optionList, "west");
    if (west < 180.0 && west >= -180.) {
        metaData.west = west;
        nsew++;
    }
    if ((nsew > 0) && (nsew < 4)) {
        printf("-E- If any of north, south, east or west are provided, ALL need to be provided.\n");
        exit(EXIT_FAILURE);
    }
    if (metaData.north <= metaData.south) {
        printf("-E- north must be greater than south.\n");
        exit(EXIT_FAILURE);
    }
    heightInDeg = metaData.north - metaData.south;
    if (heightInDeg > 180.0) {
        printf("-E- height in degrees must be less than or equal to 180.\n");
        exit(EXIT_FAILURE);
    }
    mapEast = metaData.east;
    mapWest = metaData.west;
    if (mapEast < mapWest)
        mapEast += 360;
    widthInDeg = mapEast - mapWest;

    if (widthInDeg > 360.0) {
        printf("-E- width in degrees must be less than or equal to 360.\n");
        exit(EXIT_FAILURE);
    }

    // set other fields in the metadata
    strcpy(metaData.soft_name, "l3mapgen");
    strcpy(metaData.soft_ver, softwareVersion.c_str());
    if ((tmpStr = strrchr(ifileName, '/')) != NULL)
        tmpStr++;
    else
        tmpStr = ifileName;
    strcpy(metaData.infiles, tmpStr);
    metaData.proc_con[0] = 0;
    for (i = 0; i < argc; i++) {
        strcat(metaData.proc_con, argv[i]);
        strcat(metaData.proc_con, " ");
    }
    strcpy(metaData.pversion, clo_getString(optionList, "pversion"));

    // set input parameters
    metaData.input_parms[0] = 0;
    int numOptions = clo_getNumOptions(optionList);
    for (int i = 0; i < numOptions; i++) {
        clo_option_t* option = clo_getOption(optionList, i);
        if (option) {
            if (strcmp(option->key, "help") == 0)
                continue;
            if (strcmp(option->key, "version") == 0)
                continue;
            if (strstr(option->key, "dump_options"))
                continue;
            char* val = option->valStr;
            if (val == NULL)
                val = option->defaultVal;
            if (val == NULL)
                val = "";
            strcat(metaData.input_parms, option->key);
            strcat(metaData.input_parms, "=");
            strcat(metaData.input_parms, val);
            strcat(metaData.input_parms, "|");
        }
    }

    // set processing time
    get_time(metaData.ptime);

    for (OutFile* outFile : outFiles) {
        outFile->setMetaData(&metaData);
        outFile->setDeflate(clo_getInt(optionList, "deflate"));
        outFile->setFullLatLon(clo_getBool(optionList, "full_latlon"));
    }
    if (outFile2) {
        outFile2->setMetaData(&metaData);
        outFile2->setDeflate(clo_getInt(optionList, "deflate"));
        outFile2->setFullLatLon(clo_getBool(optionList, "full_latlon"));
    }
    // writing netcdf4 file? // assume so
    if (strcasecmp(projectionStr, "raw") == 0) {
        writeRawFile(l3File, outFiles, outFile2);
    } else if (strcasecmp(projectionStr, "smi") == 0) {
        writeSmiFile(l3File, outFiles, outFile2);
    } else {
        writeProj4File(l3File, projectionStr, outFiles, outFile2, trimNSEW);
    }

    // check if no pixels are filled
    if (outFiles[0]->getNumFilledPixels() == 0) {
        printf("\nThere are no filled pixels\n");
        printf("Deleting output file.\n");

        string cmd;
        for (OutFile* outFile : outFiles) {
            cmd = "rm -f ";
            cmd += outFile->getFileName();
            system(cmd.c_str());
        }

        if (outFile2) {
            cmd = "rm -f ";
            cmd += outFile2->getFileName();
            system(cmd.c_str());
        }
        exit(110);
    }

    // check the % filled threshold
    float threshold = clo_getFloat(optionList, "threshold");
    if (threshold > 0.0) {
        if (outFiles[0]->getPercentFilledPixels() < threshold) {
            printf(
                "\nPercent filled pixels (%.1f) "
                "is below the threshold (%.1f)\n",
                outFiles[0]->getPercentFilledPixels(), threshold);
            printf("Deleting output file.\n");

            string cmd;
            for (OutFile* outFile : outFiles) {
                cmd = "rm -f ";
                cmd += outFile->getFileName();
                system(cmd.c_str());
            }

            if (outFile2) {
                cmd = "rm -f ";
                cmd += outFile2->getFileName();
                system(cmd.c_str());
            }
            exit(110);
        }
    }

    printEndInfo(outFiles[0]);
    for (OutFile* outFile : outFiles) {
        delete outFile;
    }
    if (outFile2)
        delete outFile2;
    delete l3File;
    return EXIT_SUCCESS;
}
