/*
 * l3mapgen.cpp
 *
 *  Created on: Jan 15, 2014
 *      Author: dshea
 */
#include "l3mapgen.h"

#include "OutFile.h"
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

#include <boost/lexical_cast.hpp>
#include <unordered_map>
#include <unordered_set>
#include <boost/algorithm/string.hpp>
#define EXIT_LOG(...)                                                                     \
    {                                                                                     \
        __VA_ARGS__;                                                                      \
        std::cerr << "Exiting. See line " << __LINE__ << " in " << __FILE__ << std::endl; \
        exit(EXIT_FAILURE);                                                               \
    }

using namespace std;
using namespace l3;

#define NUM_SEARCH_POINTS 51

enum MeasurementType { Avg, Stdev, Variance, Nobs, Nscenes, ObsTime, BinNum };
enum InterpType { Interp_Nearest, Interp_Bin, Interp_Linear, Interp_Area };
enum EastWest { notEastOrWest, IsEast, IsWest };
// namespace to hold 3d expansion/mapping between 3D and 2D
namespace wv3d {
std::vector<std::string> wavelength_3d_list_separated;
std::unordered_map<std::string, std::vector<int32_t>> wv3d_2d_name_to_3d_expansion;
std::unordered_map<std::string, std::string> wv3d_3d_name_to_2d_name;
std::vector<std::string> output_products_with_3d;
size_t wavelength_3d_size;
}  // namespace wv3d

const std::unordered_map<std::string, std::vector<int32_t>>& get_wv3d_2d_name_to_3d_expansion() {
    return wv3d::wv3d_2d_name_to_3d_expansion;
}

const std::unordered_map<std::string, std::string>& get_wv3d_3d_name_to_2d_name() {
    return wv3d::wv3d_3d_name_to_2d_name;
}

size_t get_len_wv3d() {
    return wv3d::wavelength_3d_size;
}

std::string get_modifier_from_measure(MeasurementType measure) {
    switch (measure) {
        case Avg:
            return "";
        case Stdev:
            return "_stdev";
        case Variance:
            return "_var";
        case Nobs:
            return "_nobs";
        case Nscenes:
            return "_nscenes";
        case ObsTime:
            return "_obs_time";
        case BinNum:
            return "_bin_num";
        default:
            return "";
    }
}
// global params for program
double widthInDeg;
double heightInDeg;
double mapEast;
double mapWest;

int imageWidth = 0;
int imageHeight = 0;
bool doingQuality = false;
bool doingRGB = false;
bool doingTransparency = false;
bool trimNSEW = false;
bool write_projtext = false;

string prodName;
vector<string> prodNameList;
vector<MeasurementType> prodMeasurementList;
string qualName;

clo_optionList_t* optionList = NULL;

// landmask
bool applyMask = false;
static grid_info_t* landmaskGrid = {0};

int isLand(float lat, float lon) {
    if (landmaskGrid != NULL) {
        double value;
        int status = get_bylatlon(landmaskGrid, lat, lon, &value);
        if (!status)
            return ((short)value != 1);  // 1 = water
    }
    return (0);  // assume water if lon, lat not found
}

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
            printf("product_rgb: %s\n", prodName.c_str());
        else
            printf("product    : %s\n", prodName.c_str());
        if (doingQuality)
            printf("qual_prod  : %s\n", qualName.c_str());
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

void printPercentDone(float percentDone) {
    static float percentPrev = 0.0;
    static const float percentDelta = 0.01;
    if (want_verbose && (percentDone - percentPrev > percentDelta)) {
        percentPrev = percentDone;
        printf("\r%2d%% complete", (int)(percentDone * 100));
        fflush(stdout);
    }
}

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

InterpType interpStr2Type(const char* str) {
    string s = str;
    boost::trim(s);
    boost::to_lower(s);

    if (s.compare("bin") == 0)
        return Interp_Bin;
    if (s.compare("linear") == 0)
        return Interp_Linear;
    if (s.compare("area") == 0)
        return Interp_Area;

    return Interp_Nearest;
}

const char* interpType2Str(InterpType interp) {
    switch (interp) {
        case Interp_Nearest:
            return "nearest";
        case Interp_Bin:
            return "bin";
        case Interp_Linear:
            return "linear";
        case Interp_Area:
            return "area";
        default:
            return "unknown";
    }
}

bool checkDateLineCrossed(double lon, double deltaLon) {
    bool crossed = false;

    float minlon = constrainLon(lon) - (deltaLon / 2.0);
    float maxlon = constrainLon(lon) + (deltaLon / 2.0);

    if (minlon > 0 && maxlon < 0)
        crossed = true;

    return crossed;
}

Box_t getBox(float lat, float lon, float deltaLat, float deltaLon, int eastwest = notEastOrWest) {
    Point_t pMin;
    Point_t pMax;
    lon = constrainLon(lon);
    if (eastwest == IsEast) {
        pMin.set<0>(180.0);
        pMin.set<1>(lat - (deltaLat / 2.0));
        pMax.set<0>(lon + (deltaLon / 2.0));
        pMax.set<1>(lat + (deltaLat / 2.0));
    } else if (eastwest == IsWest) {
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
    Box_t box(pMin, pMax);
    return box;
}

L3Bin* getBoxBins(L3File* l3File, float lat, float lon, float deltaLat, float deltaLon, float fudge,
                  bool areaWeighted, int eastwest = notEastOrWest) {
    Box_t box = getBox(lat, lon, deltaLat, deltaLon, eastwest);

    L3Bin* l3Bin = l3File->getBinsInside(box, areaWeighted);

    if (!l3Bin && fudge > 1) {  // try again with fudge factor
        Box_t box = getBox(lat, lon, deltaLat * fudge, deltaLon * fudge, eastwest);
        l3Bin = l3File->getBinsInside(box, areaWeighted);
    }
    return l3Bin;
}

/**
 * figure out if we want and can do quality processing
 * @param l3File input bin file
 * @param outFile output file
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
                qualName = clo_getRawString(optionList, "quality_product");
            } else {
                qualName = "qual_" + prodNameList[0];
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
        outFile->setQualityName(qualName);
    }
    if (outFile2) {
        outFile2->setQualityProcessing(doingQuality);
        outFile2->setQualityName(qualName);
    }
    return doingQuality;
}

void setupProduct(string& prodName, MeasurementType measure, OutFile* outFile, OutFile* outFile2) {
    // get the product info
    productInfo_t* p_info;
    p_info = allocateProductInfo();

    int sensorId = sensorName2SensorId(outFile->getMetadata()->sensor_name);
    if (sensorId == -1) {
        printf("-E- Unknown sensor name %s\n", outFile->getMetadata()->sensor_name);
        exit(EXIT_FAILURE);
    }

    if (!findProductInfo(prodName.c_str(), sensorId, p_info)) {
        printf("-E- product %s not found in XML product table\n", prodName.c_str());
        exit(EXIT_FAILURE);
    }

    // now we have to fix the p_info structure
    // Avg, Stdev, Variance, Nobs, Nscenes, ObsTime, BinNum
    string tmpStr;
    switch (measure) {
        case Avg:
            break;
        case Stdev:
            tmpStr = p_info->description;
            free(p_info->description);
            tmpStr += " (Standard Deviation)";
            p_info->description = strdup(tmpStr.c_str());
            free(p_info->palette);
            p_info->palette = strdup(PRODUCT_DEFAULT_palette);
            free(p_info->dataType);
            p_info->dataType = strdup("float");
            tmpStr = p_info->suffix;
            free(p_info->suffix);
            tmpStr += "_stdev";
            p_info->suffix = strdup(tmpStr.c_str());
            free(p_info->displayScale);
            p_info->displayScale = strdup("linear");
            p_info->addOffset = 0.0;
            p_info->scaleFactor = 1.0;
            break;
        case Variance:
            tmpStr = p_info->description;
            free(p_info->description);
            tmpStr += " (Variance)";
            p_info->description = strdup(tmpStr.c_str());
            free(p_info->palette);
            p_info->palette = strdup(PRODUCT_DEFAULT_palette);
            free(p_info->dataType);
            p_info->dataType = strdup("float");
            tmpStr = p_info->suffix;
            free(p_info->suffix);
            tmpStr += "_var";
            p_info->suffix = strdup(tmpStr.c_str());
            free(p_info->displayScale);
            p_info->displayScale = strdup("linear");
            p_info->addOffset = 0.0;
            p_info->scaleFactor = 1.0;
            break;
        case Nobs:
            tmpStr = p_info->description;
            free(p_info->description);
            tmpStr += " (number of observations)";
            p_info->description = strdup(tmpStr.c_str());
            free(p_info->units);
            p_info->units = strdup("counts");
            free(p_info->palette);
            p_info->palette = strdup(PRODUCT_DEFAULT_palette);
            free(p_info->dataType);
            p_info->dataType = strdup("short");
            tmpStr = p_info->suffix;
            free(p_info->suffix);
            tmpStr += "_nobs";
            p_info->suffix = strdup(tmpStr.c_str());
            p_info->fillValue = PRODUCT_DEFAULT_fillValue;
            p_info->validMin = 0;
            p_info->validMax = 32767;
            free(p_info->displayScale);
            p_info->displayScale = strdup("linear");
            p_info->displayMin = PRODUCT_DEFAULT_displayMin;
            p_info->displayMax = PRODUCT_DEFAULT_displayMax;
            p_info->addOffset = 0.0;
            p_info->scaleFactor = 1.0;
            break;
        case Nscenes:
            tmpStr = p_info->description;
            free(p_info->description);
            tmpStr += " (number of scenes)";
            p_info->description = strdup(tmpStr.c_str());
            free(p_info->units);
            p_info->units = strdup("counts");
            free(p_info->palette);
            p_info->palette = strdup(PRODUCT_DEFAULT_palette);
            free(p_info->dataType);
            p_info->dataType = strdup("short");
            tmpStr = p_info->suffix;
            free(p_info->suffix);
            tmpStr += "_nscenes";
            p_info->suffix = strdup(tmpStr.c_str());
            p_info->fillValue = PRODUCT_DEFAULT_fillValue;
            p_info->validMin = 0;
            p_info->validMax = 32767;
            free(p_info->displayScale);
            p_info->displayScale = strdup("linear");
            p_info->displayMin = PRODUCT_DEFAULT_displayMin;
            p_info->displayMax = PRODUCT_DEFAULT_displayMax;
            p_info->addOffset = 0.0;
            p_info->scaleFactor = 1.0;
            break;
        case ObsTime:
            tmpStr = p_info->description;
            free(p_info->description);
            tmpStr += " (observation time, TAI93)";
            p_info->description = strdup(tmpStr.c_str());
            free(p_info->units);
            p_info->units = strdup("counts");
            free(p_info->palette);
            p_info->palette = strdup(PRODUCT_DEFAULT_palette);
            free(p_info->dataType);
            p_info->dataType = strdup("float");
            tmpStr = p_info->suffix;
            free(p_info->suffix);
            tmpStr += "_obs_time";
            p_info->suffix = strdup(tmpStr.c_str());
            p_info->fillValue = PRODUCT_DEFAULT_fillValue;
            p_info->validMin = PRODUCT_DEFAULT_validMin;
            p_info->validMax = PRODUCT_DEFAULT_validMin;
            free(p_info->displayScale);
            p_info->displayScale = strdup("linear");
            p_info->displayMin = PRODUCT_DEFAULT_displayMin;
            p_info->displayMax = PRODUCT_DEFAULT_displayMax;
            p_info->addOffset = 0.0;
            p_info->scaleFactor = 1.0;
            break;
        case BinNum:
            tmpStr = p_info->description;
            free(p_info->description);
            tmpStr += " (bin ID number)";
            p_info->description = strdup(tmpStr.c_str());
            free(p_info->units);
            p_info->units = strdup("dimensionless");
            free(p_info->palette);
            p_info->palette = strdup(PRODUCT_DEFAULT_palette);
            free(p_info->dataType);
            p_info->dataType = strdup("int");
            tmpStr = p_info->suffix;
            free(p_info->suffix);
            tmpStr += "_bin_num";
            p_info->suffix = strdup(tmpStr.c_str());
            p_info->fillValue = PRODUCT_DEFAULT_fillValue;
            p_info->validMin = PRODUCT_DEFAULT_validMin;
            p_info->validMax = PRODUCT_DEFAULT_validMin;
            free(p_info->displayScale);
            p_info->displayScale = strdup("linear");
            p_info->displayMin = PRODUCT_DEFAULT_displayMin;
            p_info->displayMax = PRODUCT_DEFAULT_displayMax;
            p_info->addOffset = 0.0;
            p_info->scaleFactor = 1.0;
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
            outFile->setPalette(p_info->palette, applyMask);
            if (outFile2)
                outFile2->setPalette(p_info->palette, applyMask);
        }
    }

    // set the default scale factors for RGB
    if (doingRGB) {
        p_info->displayMin = 0.01;
        p_info->displayMax = 0.9;
        if (p_info->displayScale)
            free(p_info->displayScale);
        p_info->displayScale = strdup("log");
    }

    // override default scale parameters if set on command line
    if (clo_isSet(optionList, "datamin")) {
        p_info->displayMin = clo_getFloat(optionList, "datamin");
    }
    if (clo_isSet(optionList, "datamax")) {
        p_info->displayMax = clo_getFloat(optionList, "datamax");
    }
    if (clo_isSet(optionList, "scale_type")) {
        if (p_info->displayScale)
            free(p_info->displayScale);
        p_info->displayScale = strdup(clo_getString(optionList, "scale_type"));
    }

    outFile->addProduct(p_info);
    if (outFile2)
        outFile2->addProduct(p_info);

    freeProductInfo(p_info);
}

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

    if (clo_isSet(optionList, "suite"))
        sprintf(metaData->title, "%s Level-3 %s Mapped Image %s", metaData->sensor_name, mapDesc.c_str(), clo_getString(optionList, "suite"));
    else
        sprintf(metaData->title, "%s Level-3 %s Mapped Image", metaData->sensor_name, mapDesc.c_str());
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
    for (size_t i = 0; i < prodNameList.size(); i++) {
        OutFile* outFile;
        if (outFiles.size() == 1)
            outFile = outFiles[0];
        else
            outFile = outFiles[i];
        setupProduct(prodNameList[i], prodMeasurementList[i], outFile, outFile2);
    }

    printStartInfo(outFiles);

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
                    for (size_t prod = 0; prod < prodNameList.size(); prod++) {
                        float val;
                        switch (prodMeasurementList[prod]) {
                            case Avg:
                                val = l3Bin->getMean(prod);
                                break;
                            case Stdev:
                                val = l3Bin->getStdev(prod);
                                break;
                            case Variance:
                                val = l3Bin->getVariance(prod);
                                break;
                            case Nobs:
                                val = l3Bin->getNobs();
                                break;
                            case Nscenes:
                                val = l3Bin->getNscenes();
                                break;
                            case ObsTime:
                                val = l3Bin->getObsTime();
                                break;
                            case BinNum:
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

void writeSmiFile(L3File* l3File, vector<OutFile*> outFiles, OutFile* outFile2) {
    int32_t numFilledPixels = 0;
    meta_l3bType* metaData = outFiles[0]->getMetadata();
    double resolution = outFiles[0]->getResolution();

    string mapDesc = "Standard";
    string projName = "Equidistant Cylindrical";

    sprintf(metaData->title, "%s Level-3 %s Mapped Image", metaData->sensor_name, mapDesc.c_str());

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
    for (size_t i = 0; i < prodNameList.size(); i++) {
        OutFile* outFile;
        if (outFiles.size() == 1)
            outFile = outFiles[0];
        else
            outFile = outFiles[i];
        setupProduct(prodNameList[i], prodMeasurementList[i], outFile, outFile2);
    }
    // set up quality processing
    setupQualityProcessing(l3File, outFiles, outFile2);

    printStartInfo(outFiles);

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
                    case Interp_Nearest:
                        l3Bin = l3File->getClosestBin(lat, lon);
                        break;
                    case Interp_Bin:
                    case Interp_Area: {
                        bool areaWeighted;

                        if (interp == Interp_Area)
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
                        for (size_t prod = 0; prod < prodNameList.size(); prod++) {
                            float val;
                            switch (prodMeasurementList[prod]) {
                                case Avg:
                                    val = l3Bin->getMean(prod);
                                    break;
                                case Stdev:
                                    val = l3Bin->getStdev(prod);
                                    break;
                                case Variance:
                                    val = l3Bin->getVariance(prod);
                                    break;
                                case Nobs:
                                    val = l3Bin->getNobs();
                                    break;
                                case Nscenes:
                                    val = l3Bin->getNscenes();
                                    break;
                                case ObsTime:
                                    val = l3Bin->getObsTime();
                                    break;
                                case BinNum:
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
// the writitng is done here
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
    bool* inBox;
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
    if (clo_isSet(optionList, "lat_1")) {
        lat1Str = " +lat_1=" + to_string(clo_getFloat(optionList, "lat_1"));
    } else {
        lat1Str = " +lat_1=" + to_string(metaData->south);
    }
    string lat2Str;
    if (clo_isSet(optionList, "lat_2")) {
        lat2Str = " +lat_2=" + to_string(clo_getFloat(optionList, "lat_2"));
    } else {
        lat2Str = " +lat_2=" + to_string(metaData->north);
    }
    string aziStr;
    if (clo_isSet(optionList, "azimuth")) {
        aziStr = " +alpha=" + to_string(clo_getFloat(optionList, "azimuth"));
    }
    string utmStr;
    if (clo_isSet(optionList, "utm_zone")) {
        string zone = clo_getString(optionList, "utm_zone");
        utmStr = " +zone=";
        utmStr += zone;
        if (utmStr.find('S') != std::string::npos) {
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
    sprintf(metaData->title, "%s Level-3 %s Mapped Image", metaData->sensor_name, mapDesc.c_str());
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
    //    pj = proj_create_crs_to_crs(PJ_DEFAULT_CTX, "+proj=longlat +ellps=WGS84 +datum=WGS84",
    //    projStr.c_str(),                                NULL);
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
    // c.xyzt.t = HUGE_VAL;
    c.xyzt.t = 0.0;

    // define "inBox" region for trimNSEW option using lat/lon arrays based on
    // metaData struct
    std::array<double, 9> lats{{metaData->south, metaData->south, metaData->south,
                                (metaData->north + metaData->south) / 2., metaData->north, metaData->north,
                                metaData->north, (metaData->north + metaData->south) / 2., metaData->south}};
    std::array<double, 9> lons{{metaData->west, (metaData->east + metaData->west) / 2., metaData->east,
                                metaData->east, metaData->east, (metaData->east + metaData->west) / 2.,
                                metaData->west, metaData->west, metaData->west}};
    std::vector<Point_t> points;
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
        double latmin = 999., latmax = -999., lonmin = 999., lonmax = -999.;
        double lonmin360 = 999., lonmax360 = -999.;
        L3Shape* shape = l3File->getShape();
        for (int row = 0; row < l3File->getNumRows(); row++) {
            // First column
            L3Row* l3row = l3File->getRow(row);
            x = std::numeric_limits<double>::quiet_NaN();
            y = std::numeric_limits<double>::quiet_NaN();
            int32_t binindex = 0;

            // search for first good bin
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
                if (lon < lonmin)
                    lonmin = lon;
                if (lon > lonmax)
                    lonmax = lon;

                double lon360;
                if (lon < 0)
                    lon360 = lon + 360.0;
                else
                    lon360 = lon;
                if (lon360 < lonmin360)
                    lonmin360 = lon360;
                if (lon360 > lonmax360)
                    lonmax360 = lon360;

                if (lat < latmin)
                    latmin = lat;
                if (lat > latmax)
                    latmax = lat;

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
                    if (lon < lonmin)
                        lonmin = lon;
                    if (lon > lonmax)
                        lonmax = lon;

                    double lon360;
                    if (lon < 0)
                        lon360 = lon + 360.0;
                    else
                        lon360 = lon;
                    if (lon360 < lonmin360)
                        lonmin360 = lon360;
                    if (lon360 > lonmax360)
                        lonmax360 = lon360;

                    if (lat < latmin)
                        latmin = lat;
                    if (lat > latmax)
                        latmax = lat;

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

        // add a little padding to get last pixel
        minX -= resolution;
        maxX += resolution;
        minY -= resolution;
        maxY += resolution;

        metaData->north = latmax;
        metaData->south = latmin;

        double delta = lonmax - lonmin;
        double delta360 = lonmax360 - lonmin360;
        if (delta360 < delta) {
            if (lonmin360 > 180.0)
                metaData->west = lonmin360 - 360.0;
            else
                metaData->west = lonmin360;
            if (lonmax360 > 180.0)
                metaData->east = lonmax360 - 360.0;
            else
                metaData->east = lonmax360;
        } else {
            metaData->west = lonmin;
            metaData->east = lonmax;
        }

        mapEast = metaData->east;
        mapWest = metaData->west;
        if (mapEast < mapWest)
            mapEast += 360;
        widthInDeg = mapEast - mapWest;
        heightInDeg = metaData->north - metaData->south;

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
    for (size_t i = 0; i < prodNameList.size(); i++) {
        OutFile* outFile;
        if (outFiles.size() == 1)
            outFile = outFiles[0];
        else
            outFile = outFiles[i];
        setupProduct(prodNameList[i], prodMeasurementList[i], outFile, outFile2);
    }

    // set up quality processing
    setupQualityProcessing(l3File, outFiles, outFile2);

    printStartInfo(outFiles);

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
    double* tmpX = (double*)allocateMemory(imageWidth * sizeof(double), "tmpX");
    double* tmpY = (double*)allocateMemory(imageWidth * sizeof(double), "tmpY");
    inBox = (bool*)allocateMemory(imageWidth * sizeof(bool), "inBox");

    deltaLon = widthInDeg / imageWidth;
    deltaLat = heightInDeg / imageHeight;
    double startX = minX + resolution / 2;
    double startY = maxY - resolution / 2;
    y = startY;
    // setting values begins
    for (int row = 0; row < imageHeight; row++) {
        printPercentDone((float)row / (float)imageHeight);
        x = startX;
        std::fill(inBox, inBox + imageWidth, false);
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
                    case Interp_Nearest:
                        l3Bin = l3File->getClosestBin(lat, lon);
                        break;
                    case Interp_Bin:
                    case Interp_Area: {
                        bool areaWeighted;

                        if (interp == Interp_Area)
                            areaWeighted = true;
                        else
                            areaWeighted = false;

                        if (checkDateLineCrossed(lon, deltaLon)) {
                            // East
                            l3Bin =
                                getBoxBins(l3File, lat, lon, deltaLat, deltaLon, fudge, areaWeighted, IsEast);

                            // West
                            L3Bin* tmpBin =
                                getBoxBins(l3File, lat, lon, deltaLat, deltaLon, fudge, areaWeighted, IsWest);
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
                        for (size_t prod = 0; prod < prodNameList.size(); prod++) {
                            float val;
                            switch (prodMeasurementList[prod]) {
                                case Avg:
                                    val = l3Bin->getMean(prod);
                                    break;
                                case Stdev:
                                    val = l3Bin->getStdev(prod);
                                    break;
                                case Variance:
                                    val = l3Bin->getVariance(prod);
                                    break;
                                case Nobs:
                                    val = l3Bin->getNobs();
                                    break;
                                case Nscenes:
                                    val = l3Bin->getNscenes();
                                    break;
                                case ObsTime:
                                    val = l3Bin->getObsTime();
                                    break;
                                case BinNum:
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
            outFile->setLatLon(tmpY, tmpX);
            outFile->writeLine();
        }
        if (outFile2) {
            outFile2->setLatLon(tmpY, tmpX);
            outFile2->writeLine();
        }
        y -= resolution;
    }  // for row
       // values set ends here
    free(tmpX);
    free(tmpY);
    free(inBox);

    proj_destroy(pj_new);

    for (OutFile* outFile : outFiles) {
        if (write_projtext) {
            string projtxtfilename;
            projtxtfilename = outFile->getFileName();
            projtxtfilename += ".projtxt";
            ofstream projtxtfile(projtxtfilename);
            if (projtxtfile.is_open()) {
                projtxtfile << "# Projection information for " << outFile->getFileName() << "\n";
                projtxtfile << "proj=" << projStr << "\n";
                projtxtfile << "minX=" << std::setprecision(11) << minX << "\n";
                projtxtfile << "maxX=" << std::setprecision(11) << maxX << "\n";
                projtxtfile << "minY=" << std::setprecision(11) << minY << "\n";
                projtxtfile << "maxY=" << std::setprecision(11) << maxY << "\n";
                projtxtfile << "north=" << std::setprecision(11) << metaData->north << "\n";
                projtxtfile << "south=" << std::setprecision(11) << metaData->south << "\n";
                projtxtfile << "east=" << std::setprecision(11) << metaData->east << "\n";
                projtxtfile << "west=" << std::setprecision(11) << metaData->west << "\n";

                projtxtfile << "scale_type=" << outFile->getScaleTypeString() << "\n";
                projtxtfile << "datamin=" << std::setprecision(11) << outFile->getMinValue() << "\n";
                projtxtfile << "datamax=" << std::setprecision(11) << outFile->getMaxValue() << "\n";
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
            outFile = new OutFile_ppm_rgb();
        } else {
            if (useColor) {
                outFile = new OutFile_ppm();
            } else {
                outFile = new OutFile_pgm();
            }
        }
    } else if (oformat.compare("PNG") == 0) {
        if (doingRGB) {
            outFile = new OutFile_png_rgb();
        } else {
            outFile = new OutFile_png(useColor);
        }
    } else if (oformat.compare("TIFF") == 0) {
        if (doingRGB) {
            outFile = new OutFile_tiff_rgb();
        } else {
            if (useColor) {
                outFile = new OutFile_tiff_color();
            } else {
                outFile = new OutFile_tiff_gray();
            }
        }
    } else if (oformat.compare("HDF4") == 0) {
        outFile = new OutFile_hdf4();
    } else if (oformat.compare("netCDF4") == 0) {
        outFile = new OutFile_netcdf4();
    } else {
        printf("-E- Output file type %s not implemented\n", oformat.c_str());
        exit(EXIT_FAILURE);
    }
    if (doingTransparency)
        outFile->setTransparency();

    return outFile;
}

vector<OutFile*> makeOutputFileList(clo_optionList_t* optionList) {
    vector<OutFile*> outFiles;
    OutFile* outFile;
    string oformatStr = getFileFormatName(clo_getString(optionList, "oformat"));
    string originalOfile = clo_getString(optionList, "ofile");
    string tag = clo_getString(optionList, "ofile_product_tag");
    // if netcdf4 format specified just one file will be produced
    // need to put a check : if there 3D vars, then only netcdf is acceptable
    if (oformatStr.compare("netCDF4") == 0 || prodNameList.size() == 1 || doingRGB) {
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

        for (size_t i = 0; i < prodNameList.size(); i++) {
            std::string newName = originalOfile;
            std::string prodName_clean = prodNameList.at(i);
            std::string modifier = get_modifier_from_measure(prodMeasurementList.at(i));
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

int main(int argc, char* argv[]) {
    vector<OutFile*> outFiles;
    OutFile* outFile2 = NULL;
    char* ifileName;
    string oformat;
    int i;
    char* tmpStr;

    char softwareVersion[200];
    sprintf(softwareVersion, "%d.%d.%d-%s", VERSION_MAJOR, VERSION_MINOR, VERSION_PATCH, GITSHA);

    // init the netCDF file cache
    nc_init_chunk_cache();

    optionList = clo_createList();

    l3mapgen_init_options(optionList, softwareVersion);
    if (argc == 1) {
        clo_printUsage(optionList);
        exit(EXIT_FAILURE);
    }
    l3mapgen_read_options(optionList, argc, argv);

    if (clo_getBool(optionList, "quiet")) {
        want_verbose = 0;
    }

    ifileName = clo_getString(optionList, "ifile");

    // try SMI input file
    int oldVerbose = want_verbose;
    want_verbose = 0;
    L3File* l3File = new L3FileSMI();
    if (!l3File->open(ifileName)) {
        // try real L3 bin format
        delete l3File;
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
        prodName = clo_getRawString(optionList, "product_rgb");
    } else {
        doingRGB = false;
        if (clo_isSet(optionList, "product")) {
            prodName = clo_getRawString(optionList, "product");
        } else {
            prodName.clear();
            for (int i = 0; i < l3File->getNumProducts(); i++) {
                string tmpName = l3File->getProductName(i);
                if (tmpName != "qual_l3") {
                    if (!prodName.empty())
                        prodName += ",";
                    prodName += tmpName;
                }
            }
        }
    }
    trimNSEW = clo_getBool(optionList, "trimNSEW");
    write_projtext = clo_getBool(optionList, "write_projtext");

    boost::split(prodNameList, prodName, boost::is_any_of(","));
    string cleanProdName;
    vector<string> parts;
    meta_l3bType metaData = *l3File->getMetaData();  // get metadata
    // all the wavelength for the sensor
    int32_t* wave_array_of_the_sensor;
    // quick look up for wv
    std::unordered_set<int32_t> look_up_for_wv;
    want_verbose = 0;
    const auto total_num_bands =
        rdsensorinfo(metaData.sensorID, 0, "Lambda", (void**)&wave_array_of_the_sensor);
    // to lookup a rank
    productInfo_t* p_info;
    p_info = allocateProductInfo();
    want_verbose = 1;
    for (int i = 0; i < total_num_bands; i++) {
        look_up_for_wv.insert(wave_array_of_the_sensor[i]);
    }
    if (clo_isSet(optionList, "wavelength_3d")) {
        const std::string wavelengths_3d_list = clo_getRawString(optionList, "wavelength_3d");
        boost::split(wv3d::wavelength_3d_list_separated, wavelengths_3d_list, boost::is_any_of(", "),
                     boost::algorithm::token_compress_on);

        std::vector<std::string> colon_expanded_list;
        for (const auto& wv_par : wv3d::wavelength_3d_list_separated) {
            if (boost::contains(wv_par, ":")) {
                std::vector<std::string> pars;
                boost::split(pars, wv_par, boost::is_any_of(":"));
                if (pars.size() != 2) {
                    EXIT_LOG(std::cerr << "--Error-: Wrong range specifier: " << wv_par << std::endl;)
                    }
                try {
                    int wav_st = boost::lexical_cast<int32_t>(pars.at(0));
                    int wav_end = boost::lexical_cast<int32_t>(pars.at(1));
                    if (look_up_for_wv.count(wav_st) == 0) {
                        EXIT_LOG(std::cerr << "--Error--: The start wavelength " << wav_st
                                           << " is not found in sensor wv list.\n Check "
                                              "the range "
                                           << wv_par << std::endl);
                    }
                    if (look_up_for_wv.count(wav_end) == 0) {
                        EXIT_LOG(std::cerr << "--Error--: The end wavelength " << wav_end
                                           << " is not found in sensor wv list.\n Check "
                                              "the range "
                                           << wv_par << std::endl);
                    }
                    for (int32_t i = wav_st; i <= wav_end; i++) {
                        if (look_up_for_wv.count(i) == 0)
                            continue;
                        colon_expanded_list.push_back(boost::lexical_cast<std::string>(i));
                    }
                } catch (const boost::bad_lexical_cast& e) {
                    EXIT_LOG(std::cerr << e.what() << '\n'; std::cerr
                                                            << "--Error--: Provided wavelength are not valid "
                                                               "numbers. "
                                                            << wv_par << std::endl;)
                }
            } else {
                colon_expanded_list.push_back(wv_par);
            }
        }
        wv3d::wavelength_3d_list_separated = std::move(colon_expanded_list);
    }

    // make temp unordered set for quick lookup;
    std::unordered_set<std::string> look_up_table_product_availiable_list;
    const int number_of_products = l3File->getNumProducts();
    for (int i = 0; i < number_of_products; i++) {
        const std::string& name = l3File->getProductName(i);
        look_up_table_product_availiable_list.insert(name);
    }
    std::vector<std::string> temp_prod_name_list;
    std::unordered_set<std::string> already_set_wv;
    // check that there are no duplicates in the wv list
    for (const auto& wv : wv3d::wavelength_3d_list_separated)
        if (already_set_wv.count(wv) == 0) {
            already_set_wv.insert(wv);
        } else {
            EXIT_LOG(std::cerr << "--Error--: A duplicate found in the wavelength_3d list  " << wv
                               << std::endl);
        }

    for (size_t i = 0; i < prodNameList.size(); i++) {
        std::vector<std::string> names;
        boost::split(names, prodNameList.at(i), boost::is_any_of(":"));
        std::string clean_name = names.at(0);
        int result = findProductInfo(clean_name.c_str(), metaData.sensorID, p_info);
        if (result != 1) {
            EXIT_LOG(std::cerr << "--Error--: Could not find the product: " << clean_name << std::endl);
        }
        int prod_rank = p_info->rank;
        const std::string suffix = p_info->suffix;
        std::string prefix = p_info->prefix;

        if (look_up_table_product_availiable_list.count(clean_name) == 0) {
            if (prod_rank != 3) {
                EXIT_LOG(std::cerr << "--Error--: Non-3D Product " << clean_name
                                   << " is not found in the bin  l3 file\n");
            }
            // check if the list empty
            if (wv3d::wavelength_3d_list_separated.empty()) {
                for (const auto& product_in_l3in : look_up_table_product_availiable_list) {
                    int res = findProductInfo(product_in_l3in.c_str(), metaData.sensorID, p_info);

                    if (res != 1) {
                        EXIT_LOG(std::cerr << "--Error--: Could not read the product info " << product_in_l3in
                                           << std::endl);
                    }
                    const std::string local_suffix = p_info->suffix;
                    const std::string local_name = p_info->productName;
                    if (boost::contains(clean_name, local_name)) {
                        if (!local_suffix.empty()) {
                            if (!boost::contains(clean_name, local_suffix))
                                continue;
                        } else {
                            if (clean_name != local_name)
                                continue;
                        }

                        const std::string wave_length = boost::lexical_cast<std::string>(p_info->prod_ix);
                        if (wave_length.empty()) {
                            EXIT_LOG(std::cerr << "--Error--: Not valid 2D slice of " << clean_name
                                               << "in the l3bin file " << std::endl);
                        }
                        wv3d::wavelength_3d_list_separated.push_back(wave_length);
                    }
                }
                std::sort(wv3d::wavelength_3d_list_separated.begin(),
                          wv3d::wavelength_3d_list_separated.end());
            }
            bool prod_3d_expand_found = false;
            for (size_t i = 0; i < wv3d::wavelength_3d_list_separated.size(); i++) {
                std::string wv = wv3d::wavelength_3d_list_separated.at(i);
                std::string prod_3d_name = prefix + "_" + wv + suffix;
                if (look_up_table_product_availiable_list.count(prod_3d_name) == 0) {
                    EXIT_LOG(std::cerr << "--Error--: Neither product " << clean_name
                                       << " or its wavelength 3d slice " << prod_3d_name
                                       << " are found. \nExiting ... " << std::endl);
                } else {
                    prod_3d_expand_found = true;
                    names.at(0) = prod_3d_name;
                    int32_t wavelength;
                    try {
                        wavelength = boost::lexical_cast<int32_t>(wv);
                    } catch (const boost::bad_lexical_cast& e) {
                        EXIT_LOG(std::cerr << e.what() << '\n';
                                 std::cerr << "--Error--: Provided wavelength are not valid "
                                              "numbers. \nExiting...");
                    }

                    wv3d::wv3d_2d_name_to_3d_expansion[clean_name].push_back(wavelength);
                    wv3d::wv3d_3d_name_to_2d_name[prod_3d_name] = clean_name;
                    std::string temp_prod_3d_name;
                    for (const auto& name : names) {
                        if (!temp_prod_3d_name.empty())
                            temp_prod_3d_name += ":";
                        temp_prod_3d_name += name;
                    }
                    temp_prod_name_list.push_back(temp_prod_3d_name);
                }
            }
            if (!prod_3d_expand_found) {
                EXIT_LOG(std::cerr << "--Error--: Product not found : " << clean_name << std::endl);
            }

        } else {
            if (prod_rank != 2) {
                EXIT_LOG(std::cerr << "--Error--: The product in the bin file " << clean_name
                                   << " is not a 2D prodcut" << std::endl);
            }
            temp_prod_name_list.push_back(prodNameList.at(i));
        }
        wv3d::output_products_with_3d.push_back(clean_name);
    }
    prodNameList = std::move(temp_prod_name_list);
    wv3d::wavelength_3d_size = wv3d::wavelength_3d_list_separated.size();
    if (wv3d::wavelength_3d_size >= 1) {
        if (wv3d::wv3d_2d_name_to_3d_expansion.size() > 0) {
            if (clo_isSet(optionList, "ofile2")) {
                const std::string oformatStr2 = getFileFormatName(clo_getString(optionList, "oformat2"));
                if (oformatStr2.compare("netCDF4") != 0) {
                    EXIT_LOG(std::cerr << "The user supplied a 3D product "
                                       << " and the output format for ofile2 is not netCDF4.\n"
                                       << "Exiting ... " << std::endl);
                }
            }
        }
    }
    // setup measurments
    for (size_t i = 0; i < prodNameList.size(); i++) {
        if (i != 0)
            cleanProdName += ",";
        boost::split(parts, prodNameList[i], boost::is_any_of(":"));
        if (parts.size() == 1) {
            cleanProdName += parts[0];
            prodMeasurementList.push_back(Avg);
        } else if (parts.size() == 2) {
            prodNameList[i] = parts[0];  // get rid of the modifier
            cleanProdName += parts[0];
            if (parts[1].compare("avg") == 0)
                prodMeasurementList.push_back(Avg);
            else if (parts[1].compare("stdev") == 0)
                prodMeasurementList.push_back(Stdev);
            else if (parts[1].compare("var") == 0)
                prodMeasurementList.push_back(Variance);
            else if (parts[1].compare("nobs") == 0)
                prodMeasurementList.push_back(Nobs);
            else if (parts[1].compare("nscenes") == 0)
                prodMeasurementList.push_back(Nscenes);
            else if (parts[1].compare("obs_time") == 0)
                prodMeasurementList.push_back(ObsTime);
            else if (parts[1].compare("bin_num") == 0)
                prodMeasurementList.push_back(BinNum);
            else {
                EXIT_LOG(
                    printf("-E- measurement type \"%s\" "
                           "not understood for product \"%s\".\n",
                           parts[1].c_str(), parts[0].c_str()));
            }
        } else {
            EXIT_LOG(printf("-E- product name not understood \"%s\".\n", prodNameList[i].c_str()));
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
    strcpy(metaData.soft_ver, softwareVersion);
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
