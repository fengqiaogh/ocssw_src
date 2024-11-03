#include <xtiffio.h>
#include <geo_normalize.h>

#include "OutFile.h"

#include <stdio.h>
#include <iostream>
#include <math.h>
#include <timeutils.h>
#include <nc4utils.h>
#include <string>
#include <float.h>
#include <regex>

#include <hdf.h>
#include <mfhdf.h>

#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/case_conv.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string.hpp>
#include <stack>

#include <geo_tiffp.h>

using namespace std;

const std::unordered_map<std::string, std::vector<int32_t>>& get_wv3d_2d_name_to_3d_expansion();
const std::unordered_map<std::string, std::string>& get_wv3d_3d_name_to_2d_name();
size_t get_len_wv3d();

std::string remove_wv_from_long_name(const std::string& long_name_with_wv) {
    std::vector<std::string> key_words;
    const std::string delim = " ";
    boost::split(key_words, long_name_with_wv, boost::is_any_of(delim));
    const std::string start = "at";
    const std::string end = "nm";
    std::string out;
    std::string wv_prefix;
    std::stack<std::string> stack;
    for (const auto& word : key_words) {
        if (stack.empty()) {
            if (word != "at") {
                if (!out.empty())
                    out += delim;
                out += word;
            } else {
                stack.push(word);
            }
        } else {
            if (word != "nm") {
                stack.push(word);
            } else {
                while (!stack.empty()) {
                    if (!wv_prefix.empty())
                        wv_prefix += delim;
                    wv_prefix += stack.top();
                    stack.pop();
                }
            }
        }
    }
    boost::trim(out);
    // {
    //     std::cout << "OUT LONG NAME WITHOUT PREFIX " << out << " THE PREFIX "
    //               << wv_prefix << std::endl;
    // }
    return out;
}

//------------------------------------------------------------------------------
// OutFile::ProductStuff
//------------------------------------------------------------------------------

OutFile::ProductStuff::ProductStuff(int32_t width, const productInfo_t* productInfo, double landPixelValue) {
    this->width = width;
    this->productInfo = allocateProductInfo();
    copyProductInfo(this->productInfo, productInfo);
    dataStorage = UByteDS;
    scaleType = Linear;
    scale = 1.0;
    offset = 0.0;
    minOutputVal = 0;
    maxOutputVal = 255;
    minVal = minOutputVal;
    maxVal = maxOutputVal;
    missingValue = fillPix;
    lineData = (double*)allocateMemory(width * sizeof(double), "OutFile::ProductStuff::lineData");
    this->landPixelValue = landPixelValue;
}

OutFile::ProductStuff::ProductStuff(const OutFile::ProductStuff& pStuff) {
    width = pStuff.width;
    productInfo = allocateProductInfo();
    copyProductInfo(productInfo, pStuff.productInfo);
    dataStorage = pStuff.dataStorage;
    scaleType = pStuff.scaleType;
    scale = pStuff.scale;
    offset = pStuff.offset;
    minOutputVal = pStuff.minOutputVal;
    maxOutputVal = pStuff.maxOutputVal;
    minVal = pStuff.minVal;
    maxVal = pStuff.maxVal;
    missingValue = pStuff.missingValue;
    lineData = (double*)allocateMemory(width * sizeof(double), "OutFile::ProductStuff::lineData");
    landPixelValue = pStuff.landPixelValue;
}

OutFile::ProductStuff::~ProductStuff() {
    freeProductInfo(productInfo);
    free(lineData);
}

/**
 * set the scale factors.  note that default minOutputVal=0, maxOutputVal=255
 * @param min min geophysical value
 * @param max max geophysical value
 * @param log do you want log10 scaling
 */
void OutFile::ProductStuff::setScale(double min, double max, ScaleType scaleType) {
    this->scaleType = scaleType;
    minVal = min;
    maxVal = max;
    switch (scaleType) {
        case Linear:
            offset = min;
            scale = (max - offset) / (maxOutputVal - minOutputVal);
            break;
        case Log:
            offset = log10(min);
            scale = (log10(max) - offset) / (maxOutputVal - minOutputVal);
            break;
        case ArcTan:
            offset = min;
            scale = max;
            minVal = calcPhysicalVal(minOutputVal);
            maxVal = calcPhysicalVal(maxOutputVal);
            break;
        default:
            printf("-E- OutFile::setScale - invalid scaleType = %d\n", (int)scaleType);
            exit(EXIT_FAILURE);
    }
}

void OutFile::ProductStuff::setScale(double min, double max, ScaleType scaleType, double minOutput,
                                     double maxOutput) {
    minOutputVal = minOutput;
    maxOutputVal = maxOutput;
    setScale(min, max, scaleType);
}

/**
 * set the scale factors.  note that default minOutputVal=0, maxOutputVal=255
 * @param scale slope
 * @param offset intercept
 * @param scaleType type of scaling to calculate
 */
void OutFile::ProductStuff::setScaleOffset(double scale, double offset, ScaleType scaleType) {
    this->scaleType = scaleType;
    this->scale = scale;
    this->offset = offset;

    // have to set these so calcPhysicalVal does not limit the physical val
    minVal = 0 - FLT_MAX;
    maxVal = FLT_MAX;
    minVal = calcPhysicalVal(minOutputVal);
    maxVal = calcPhysicalVal(maxOutputVal);
}

void OutFile::ProductStuff::setScaleOffset(double scale, double offset, ScaleType scaleType, double minOutput,
                                           double maxOutput) {
    minOutputVal = minOutput;
    maxOutputVal = maxOutput;
    setScaleOffset(scale, offset, scaleType);
}

double OutFile::ProductStuff::calcOutputVal(double val) const {
    double outVal;

    if (val == badPixelValue)
        return missingValue;

    if (val == landPixelValue)
        return landPix;

    // don't scale if output type is floating point
    if (dataStorage == FloatDS || dataStorage == DoubleDS)
        return val;

    switch (scaleType) {
        case Linear:
            outVal = (val - offset) / scale;
            break;
        case Log:
            if (val < 0)
                return minOutputVal;
            else
                outVal = (log10(val) - offset) / scale;
            break;
        case ArcTan:
            outVal = scale * ((atan(0.5 * val - offset) / atan(offset)) + 1);
            break;
        default:
            printf("-E- OutFile::ProductStuff::calcOutputVal - invalid scaleType = %d\n", (int)scaleType);
            exit(EXIT_FAILURE);
    }

    if (outVal < minOutputVal)
        return minOutputVal;
    if (outVal > maxOutputVal)
        return maxOutputVal;
    return outVal;
}

double OutFile::ProductStuff::calcPhysicalVal(double val) const {
    double physicalVal;

    if (val == missingValue)
        return badPixelValue;

    // don't scale if out output type is floating point
    if (dataStorage == FloatDS || dataStorage == DoubleDS)
        return val;

    switch (scaleType) {
        case Linear:
            physicalVal = val * scale + offset;
            break;
        case Log:
            physicalVal = pow(10, val * scale + offset);
            break;
        case ArcTan:
            physicalVal = (tan((val / scale - 1) * atan(offset)) + offset) / 0.5;
            break;
        default:
            printf("-E- OutFile::ProductStuff::calcPhysicalVal - invalid scaleType = %d\n", (int)scaleType);
            exit(EXIT_FAILURE);
    }

    if (physicalVal < minVal)
        return minVal;
    if (physicalVal > maxVal)
        return maxVal;
    return physicalVal;
}

void OutFile::ProductStuff::calcOutputLineVals(void* lineBuffer) const {
    switch (dataStorage) {
        case ByteDS:
            for (int i = 0; i < width; i++)
                ((int8_t*)lineBuffer)[i] = round(calcOutputVal(lineData[i]));
            break;
        case UByteDS:
            for (int i = 0; i < width; i++)
                ((uint8_t*)lineBuffer)[i] = round(calcOutputVal(lineData[i]));
            break;
        case ShortDS:
            for (int i = 0; i < width; i++)
                ((int16_t*)lineBuffer)[i] = round(calcOutputVal(lineData[i]));
            break;
        case UShortDS:
            for (int i = 0; i < width; i++)
                ((uint16_t*)lineBuffer)[i] = round(calcOutputVal(lineData[i]));
            break;
        case IntDS:
            for (int i = 0; i < width; i++)
                ((int32_t*)lineBuffer)[i] = round(calcOutputVal(lineData[i]));
            break;
        case UIntDS:
            for (int i = 0; i < width; i++)
                ((uint32_t*)lineBuffer)[i] = round(calcOutputVal(lineData[i]));
            break;
        case FloatDS:
            for (int i = 0; i < width; i++)
                ((float*)lineBuffer)[i] = lineData[i];
            break;
        case DoubleDS:
            for (int i = 0; i < width; i++)
                ((double*)lineBuffer)[i] = lineData[i];
            break;
        default:
            printf("-E- OutFile::ProductStuff::calcOutputLineVals - unrecognized data type = %d\n",
                   dataStorage);
            exit(EXIT_FAILURE);
    }
}

//------------------------------------------------------------------------------
// OutFile
//------------------------------------------------------------------------------

OutFile::OutFile() {
    landPixelValue = -32766.0;
    width = 0;
    height = 0;
    qualityData = NULL;
    currentLine = 0;
    colorType = Grayscale;
    fileMinVal = DBL_MAX;
    fileMaxVal = 0 - DBL_MAX;
    resolution = 0;
    deflate = 0;
    fullLatLon = LatLon2D;
    latData = NULL;
    lonData = NULL;

    red = (uint8_t*)allocateMemory(256, "red");
    green = (uint8_t*)allocateMemory(256, "green");
    blue = (uint8_t*)allocateMemory(256, "blue");

    // init to grayscale
    for (int i = 0; i < 256; i++) {
        red[i] = i;
        green[i] = i;
        blue[i] = i;
    }
    rgb_land = (uint8_t*)allocateMemory(3, "rgb_land");
    rgb_land[0] = 160;
    rgb_land[1] = 82;
    rgb_land[2] = 45;
    transparent = false;

    metaData = (meta_l3bType*)allocateMemory(sizeof(meta_l3bType), "OutFile::metaData");
    metaData->north = 90.0;
    metaData->south = -90.0;
    metaData->east = 180.0;
    metaData->west = -180.0;

    mapProjection = "Undefined";
}

OutFile::~OutFile() {
    free(red);
    free(green);
    free(blue);
    free(rgb_land);
    for (size_t i = 0; i < productStuff.size(); i++) {
        delete productStuff[i];
    }
    productStuff.clear();
    free(metaData);
    if (qualityData)
        free(qualityData);
}

std::string OutFile::getScaleTypeString(int32_t prod) {
    switch (productStuff[prod]->scaleType) {
        case Log:
            return "Log";
        case Linear:
            return "Linear";
        case ArcTan:
            return "ArcTan";
    }
    return "Linear";
}

void OutFile::setSize(int32_t width, int32_t height) {
    this->width = width;
    this->height = height;
    if (qualityData) {
        free(qualityData);
        qualityData = (uint8_t*)allocateMemory(width, "OutFile::qualityData");
    }
    currentLine = 0;

    for (size_t i = 0; i < productStuff.size(); i++) {
        delete productStuff[i];
    }
    productStuff.clear();
}

int32_t OutFile::getWidth() const {
    return width;
}

int32_t OutFile::getHeight() const {
    return height;
}

void OutFile::setFileName(string fileName) {
    this->fileName = fileName;
}

void OutFile::setPixel(int32_t x, double val, int32_t prod) {
    if (x < 0 || x >= width) {
        fprintf(stderr, "-E- OutFile::setPixel - x=%d is not within range, width=%d.\n", x, width);
        exit(EXIT_FAILURE);
    }
    productStuff[prod]->lineData[x] = val;
    if (val > fileMaxVal)
        fileMaxVal = val;
    if (val < fileMinVal)
        fileMinVal = val;
}

void OutFile::setPixelRGB(int32_t x, float red, float green, float blue) {
    fprintf(stderr, "-E- OutFile::setPixelRGB - RGB not implemented with this file type.\n");
    exit(EXIT_FAILURE);
}

void OutFile::setTransparency() {
    transparent = true;
}

void OutFile::setQuality(int32_t x, uint8_t val) {
    if (x < 0 || x >= width) {
        fprintf(stderr, "-E- OutFile::setQuality - x=%d is not within range, width=%d.\n", x, width);
        exit(EXIT_FAILURE);
    }
    if (!qualityData) {
        fprintf(stderr, "-E- OutFile::setQuality - qualityData id NULL.\n");
        exit(EXIT_FAILURE);
    }
    qualityData[x] = val;
}

void OutFile::landPixel(int32_t x) {
    if (x < 0 || x >= width) {
        fprintf(stderr, "-E- OutFile::landPixel - x=%d is not within range, width=%d.\n", x, width);
        exit(EXIT_FAILURE);
    }
    for (size_t prod = 0; prod < productStuff.size(); prod++)
        productStuff[prod]->lineData[x] = landPixelValue;
    if (qualityData)
        qualityData[x] = qualityUnused;
}

void OutFile::fillPixel(int32_t x) {
    if (x < 0 || x >= width) {
        fprintf(stderr, "-E- OutFile::fillPixel - x=%d is not within range, width=%d.\n", x, width);
        exit(EXIT_FAILURE);
    }
    for (size_t prod = 0; prod < productStuff.size(); prod++)
        productStuff[prod]->lineData[x] = badPixelValue;
    if (qualityData)
        qualityData[x] = qualityUnused;
}

void OutFile::missingPixel(int32_t x) {
    if (x < 0 || x >= width) {
        fprintf(stderr, "-E- OutFile::missingPixel - x=%d is not within range, width=%d.\n", x, width);
        exit(EXIT_FAILURE);
    }
    for (size_t prod = 0; prod < productStuff.size(); prod++)
        productStuff[prod]->lineData[x] = badPixelValue;
    if (qualityData)
        qualityData[x] = qualityUnused;
}

void OutFile::setLatLon(double* lat, double* lon) {
    if (fullLatLon == LatLon2D) {
        // NcVar::putVar takes care of converting double to float
        latData = lat;
        lonData = lon;
    }
}

bool OutFile::setPalette(const char* paletteName, bool applyMask) {
    char* dataRoot;
    string paletteFileName;
    short r[256], g[256], b[256];

    if ((dataRoot = getenv("OCDATAROOT")) == NULL) {
        printf("OCDATAROOT environment variable is not defined.\n");
        return (EXIT_FAILURE);
    }
    paletteFileName = dataRoot;
    paletteFileName += "/common/palette/";
    paletteFileName += paletteName;
    paletteFileName += ".pal";

    if (getlut_file((char*)paletteFileName.c_str(), r, g, b)) {
        fprintf(stderr, "Error reading palette file %s\n", paletteFileName.c_str());
        return false;
    }
    if (applyMask) {
        r[landPix] = rgb_land[0];
        g[landPix] = rgb_land[1];
        b[landPix] = rgb_land[2];
        r[fillPix] = 0;
        g[fillPix] = 0;
        b[fillPix] = 0;
    }
    for (int i = 0; i < 256; i++) {
        red[i] = r[i];
        green[i] = g[i];
        blue[i] = b[i];
    }

    return true;
}

void OutFile::setLandRGB(const char* rgb_land_string) {
    vector<string> rgb;
    boost::split(rgb, rgb_land_string, boost::is_any_of(","));

    rgb_land[0] = std::stoi(rgb[0]);
    rgb_land[1] = std::stoi(rgb[1]);
    rgb_land[2] = std::stoi(rgb[2]);
}

void OutFile::setMetaData(meta_l3bType* metaData) {
    *this->metaData = *metaData;
}

/**
 * Add a product for display type output files
 * @param productInfo info structure to copy
 * @return the index for the new product
 */
int32_t OutFile::addProduct(productInfo_t* productInfo) {
    ProductStuff* stuff = new ProductStuff(width, productInfo, landPixelValue);

    // setup display scaling
    if (!strcmp(productInfo->displayScale, "linear"))
        stuff->setScale(productInfo->displayMin, productInfo->displayMax, Linear);
    else if (!strcmp(productInfo->displayScale, "log"))
        stuff->setScale(productInfo->displayMin, productInfo->displayMax, Log);
    else if (!strcmp(productInfo->displayScale, "arctan"))
        stuff->setScale(productInfo->displayMin, productInfo->displayMax, ArcTan);
    else {
        printf("-E- OutFile::addProduct - invalid displayScale = %s\n", productInfo->displayScale);
        exit(EXIT_FAILURE);
    }

    productStuff.push_back(stuff);
    return productStuff.size() - 1;
}

int32_t OutFile::addProductNonDisplay(productInfo_t* productInfo) {
    ProductStuff* stuff = new ProductStuff(width, productInfo, landPixelValue);

    if (!strcmp(productInfo->dataType, "byte")) {
        stuff->dataStorage = ByteDS;
        stuff->minOutputVal = SCHAR_MIN;
        stuff->maxOutputVal = SCHAR_MAX;
    } else if (!strcmp(productInfo->dataType, "ubyte")) {
        stuff->dataStorage = UByteDS;
        stuff->minOutputVal = 0;
        stuff->maxOutputVal = UCHAR_MAX;
    } else if (!strcmp(productInfo->dataType, "short")) {
        stuff->dataStorage = ShortDS;
        stuff->minOutputVal = SHRT_MIN;
        stuff->maxOutputVal = SHRT_MAX;
    } else if (!strcmp(productInfo->dataType, "ushort")) {
        stuff->dataStorage = UShortDS;
        stuff->minOutputVal = 0;
        stuff->maxOutputVal = USHRT_MAX;
    } else if (!strcmp(productInfo->dataType, "int")) {
        stuff->dataStorage = IntDS;
        stuff->minOutputVal = INT_MIN;
        stuff->maxOutputVal = INT_MAX;
    } else if (!strcmp(productInfo->dataType, "uint")) {
        stuff->dataStorage = UIntDS;
        stuff->minOutputVal = 0;
        stuff->maxOutputVal = UINT_MAX;
    } else if (!strcmp(productInfo->dataType, "float")) {
        stuff->dataStorage = FloatDS;
        stuff->minOutputVal = 0;
        stuff->maxOutputVal = 0;
    } else if (!strcmp(productInfo->dataType, "double")) {
        stuff->dataStorage = DoubleDS;
        stuff->minOutputVal = 0;
        stuff->maxOutputVal = 0;
    } else {
        printf("-E- OutFile::addProductNonDisplay - invalid data type = %s\n", productInfo->dataType);
        exit(EXIT_FAILURE);
    }

    // setup scaling
    stuff->setScaleOffset(productInfo->scaleFactor, productInfo->addOffset, Linear);
    stuff->missingValue = productInfo->fillValue;

    productStuff.push_back(stuff);
    return productStuff.size() - 1;
}

void OutFile::setMapProjection(string projection) {
    mapProjection = projection;
}

void OutFile::setProj4Info(string projStr, double minX, double maxY) {
    proj4String = projStr;

    for (int i = 0; i < 6; i++)
        tiepoints[i] = 0;
    tiepoints[3] = minX;
    tiepoints[4] = maxY;

    pixscale[0] = resolution;
    pixscale[1] = resolution;
    pixscale[2] = 0;
}

void OutFile::setNumFilledPixels(int32_t num) {
    if (metaData) {
        metaData->data_bins = num;
    }
}

int32_t OutFile::getNumFilledPixels() {
    if (metaData)
        return metaData->data_bins;
    else
        return -1;
}

float OutFile::getPercentFilledPixels() {
    if (metaData) {
        float numPix = width * height;
        return metaData->data_bins / numPix * 100.0;
    } else
        return -1;
}

void OutFile::resetFileMinMax() {
    fileMinVal = DBL_MAX;
    fileMaxVal = 0 - DBL_MAX;
}

void OutFile::setResolution(string resolutionStr) {
    resolution = string2resolution(resolutionStr);
}

void OutFile::setQualityProcessing(bool val) {
    if (val) {
        // if width is not set yet allocate some dummy memory to flag we want
        // to do SST quality processing.
        if (width <= 0) {
            if (!qualityData)
                qualityData = (uint8_t*)allocateMemory(2, "OutFile::qualityData");
        } else {
            if (qualityData)
                free(qualityData);
            qualityData = (uint8_t*)allocateMemory(width, "OutFile::qualityData");
        }
    } else {
        if (qualityData) {
            free(qualityData);
            qualityData = NULL;
        }
    }
}

bool OutFile::getQualityProcessing() {
    if (qualityData)
        return true;
    else
        return false;
}

//------------------------------------------------------------------------------
// OutFile_pgm
//------------------------------------------------------------------------------

OutFile_pgm::OutFile_pgm() : OutFile() {
    outfp = NULL;
    fileData = NULL;
}

OutFile_pgm::~OutFile_pgm() {
    if (outfp)
        fclose(outfp);
    if (fileData)
        free(fileData);
}

void OutFile_pgm::setSize(int32_t width, int32_t height) {
    OutFile::setSize(width, height);
    if (fileData)
        free(fileData);
    fileData = (uint8_t*)allocateMemory(width, "OutFile_pgm::fileData");
}

bool OutFile_pgm::open() {
    outfp = fopen(fileName.c_str(), "w");
    if (!outfp)
        return false;

    /* Write pgm header */
    fprintf(outfp, "P5\n");
    fprintf(outfp, "%d %d\n", width, height);
    fprintf(outfp, "255\n");

    return true;
}

void OutFile_pgm::writeLine() {
    for (int i = 0; i < width; i++)
        fileData[i] = (uint8_t)productStuff[0]->calcOutputVal(productStuff[0]->lineData[i]);
    fwrite(fileData, 1, width, outfp);
    currentLine++;
}

bool OutFile_pgm::close() {
    fclose(outfp);
    outfp = NULL;
    return true;
}

//------------------------------------------------------------------------------
// OutFile_ppm
//------------------------------------------------------------------------------

void OutFile_ppm::setSize(int32_t width, int32_t height) {
    OutFile::setSize(width, height);
    if (fileData)
        free(fileData);
    fileData = (uint8_t*)allocateMemory(width * 3, "OutFile_ppm::fileData");
}

bool OutFile_ppm::open() {
    outfp = fopen(fileName.c_str(), "w");
    if (!outfp)
        return false;

    /*
     * Write ppm file header
     */
    fprintf(outfp, "P6\n");
    fprintf(outfp, "%d %d\n", width, height);
    fprintf(outfp, "255\n");

    return true;
}

void OutFile_ppm::writeLine() {
    int j = 0;
    for (int i = 0; i < width; i++) {
        uint8_t val = (uint8_t)round(productStuff[0]->calcOutputVal(productStuff[0]->lineData[i]));
        fileData[j++] = red[val];
        fileData[j++] = green[val];
        fileData[j++] = blue[val];
    }
    fwrite(fileData, 1, width * 3, outfp);
    currentLine++;
}

//------------------------------------------------------------------------------
// OutFile_ppm_rgb
//------------------------------------------------------------------------------

OutFile_ppm_rgb::OutFile_ppm_rgb() : OutFile() {
    outfp = NULL;
    fileData = NULL;
}

OutFile_ppm_rgb::~OutFile_ppm_rgb() {
    if (outfp)
        fclose(outfp);
    if (fileData)
        free(fileData);
}

void OutFile_ppm_rgb::setSize(int32_t width, int32_t height) {
    OutFile::setSize(width, height);
    if (fileData)
        free(fileData);
    fileData = (uint8_t*)allocateMemory(width * 3, "OutFile_ppm_rgb::fileData");
}

bool OutFile_ppm_rgb::open() {
    outfp = fopen(fileName.c_str(), "w");
    if (!outfp)
        return false;

    /*
     * Write ppm file header
     */
    fprintf(outfp, "P6\n");
    fprintf(outfp, "%d %d\n", width, height);
    fprintf(outfp, "255\n");
    return true;
}

bool OutFile_ppm_rgb::close() {
    fclose(outfp);
    outfp = NULL;
    return true;
}

void OutFile_ppm_rgb::setPixel(int32_t x, double val, int32_t prod) {
    fprintf(stderr,
            "-E- OutFile_ppm_rgb::setPixel - only RGB is implemented with this file type.\n");
    exit(EXIT_FAILURE);
}

void OutFile_ppm_rgb::setPixelRGB(int32_t x, float red, float green, float blue) {
    if (x < 0 || x >= width) {
        fprintf(stderr,
                "-E- OutFile_ppm_rgb::setPixel - x=%d is not within range, width=%d.\n",
                x, width);
        exit(EXIT_FAILURE);
    }

    uint8_t* ptr = fileData + x * 3;

    *ptr = round(productStuff[0]->calcOutputVal(red));
    ptr++;
    *ptr = round(productStuff[0]->calcOutputVal(green));
    ptr++;
    *ptr = round(productStuff[0]->calcOutputVal(blue));

    // do this to keep the file min/max reasonable
    if (red > fileMaxVal)
        fileMaxVal = red;
    if (red < fileMinVal)
        fileMinVal = red;
}

void OutFile_ppm_rgb::landPixel(int32_t x) {
    if (x < 0 || x >= width) {
        fprintf(stderr,
                "-E- OutFile_ppm_rgb::landPixel - x=%d is not within range, width=%d.\n",
                x, width);
        exit(EXIT_FAILURE);
    }
    uint8_t* ptr = fileData + x * 3;
    *ptr = landPix;
    ptr++;
    *ptr = landPix;
    ptr++;
    *ptr = landPix;
}

void OutFile_ppm_rgb::fillPixel(int32_t x) {
    if (x < 0 || x >= width) {
        fprintf(stderr,
                "-E- OutFile_ppm_rgb::fillPixel - x=%d is not within range, width=%d.\n",
                x, width);
        exit(EXIT_FAILURE);
    }
    uint8_t* ptr = fileData + x * 3;
    *ptr = fillPix;
    ptr++;
    *ptr = fillPix;
    ptr++;
    *ptr = fillPix;
}

void OutFile_ppm_rgb::missingPixel(int32_t x) {
    if (x < 0 || x >= width) {
        fprintf(stderr,
                "-E- OutFile_ppm_rgb::missingPixel - x=%d is not within range, width=%d.\n",
                x, width);
        exit(EXIT_FAILURE);
    }
    uint8_t* ptr = fileData + x * 3;
    *ptr = fillPix;
    ptr++;
    *ptr = fillPix;
    ptr++;
    *ptr = fillPix;
}

void OutFile_ppm_rgb::writeLine() {
    fwrite(fileData, 1, width * 3, outfp);
    currentLine++;
}

//------------------------------------------------------------------------------
// OutFile_png
//------------------------------------------------------------------------------

OutFile_png::OutFile_png(bool color) : OutFile() {
    isColor = color;
    outfp = NULL;
    fileData = NULL;
    info_ptr = NULL;
    png_ptr = NULL;
    num_text = 10;
}

OutFile_png::~OutFile_png() {
    if (outfp)
        fclose(outfp);
    if (fileData)
        free(fileData);
}

void OutFile_png::setSize(int32_t width, int32_t height) {
    OutFile::setSize(width, height);
    if (fileData)
        free(fileData);
    fileData = (uint8_t*)allocateMemory(width, "OutFile_png::fileData");
}

bool OutFile_png::open() {
    outfp = fopen(fileName.c_str(), "w");
    if (!outfp)
        return false;

    png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if (!png_ptr) {
        fprintf(stderr, "-E- Unable to create PNG write structure.\n");
        exit(EXIT_FAILURE);
    }

    info_ptr = png_create_info_struct(png_ptr);
    if (!info_ptr) {
        png_destroy_write_struct(&png_ptr, (png_infopp)NULL);
        fprintf(stderr, "-E- Unable to create PNG info structure.\n");
        exit(EXIT_FAILURE);
    }
    if (setjmp(png_jmpbuf(png_ptr))) {
        png_destroy_write_struct(&png_ptr, (png_infopp)NULL);
        fprintf(stderr, "-E- Unable to call PNG setjmp().\n");
        exit(EXIT_FAILURE);
    }
    png_init_io(png_ptr, outfp);

    png_text text_ptr[num_text];
    text_ptr[0].key = "projString";
    text_ptr[0].text = const_cast<char*>(proj4String.c_str());
    text_ptr[0].compression = PNG_TEXT_COMPRESSION_NONE;
    text_ptr[1].key = "resolution";
    text_ptr[1].text = const_cast<char*>(std::to_string(resolution).c_str());
    text_ptr[1].compression = PNG_TEXT_COMPRESSION_NONE;
    text_ptr[2].key = "north";
    text_ptr[2].text = const_cast<char*>(std::to_string(metaData->north).c_str());
    text_ptr[2].compression = PNG_TEXT_COMPRESSION_NONE;
    text_ptr[3].key = "south";
    text_ptr[3].text = const_cast<char*>(std::to_string(metaData->south).c_str());
    text_ptr[3].compression = PNG_TEXT_COMPRESSION_NONE;
    text_ptr[4].key = "east";
    text_ptr[4].text = const_cast<char*>(std::to_string(metaData->east).c_str());
    text_ptr[4].compression = PNG_TEXT_COMPRESSION_NONE;
    text_ptr[5].key = "west";
    text_ptr[5].text = const_cast<char*>(std::to_string(metaData->west).c_str());
    text_ptr[5].compression = PNG_TEXT_COMPRESSION_NONE;
    text_ptr[6].key = "minX";
    text_ptr[6].text = const_cast<char*>(std::to_string(tiepoints[3]).c_str());
    text_ptr[6].compression = PNG_TEXT_COMPRESSION_NONE;
    text_ptr[7].key = "maxY";
    text_ptr[7].text = const_cast<char*>(std::to_string(tiepoints[4]).c_str());
    text_ptr[7].compression = PNG_TEXT_COMPRESSION_NONE;
    text_ptr[8].key = "height";
    text_ptr[8].text = const_cast<char*>(std::to_string(height).c_str());
    text_ptr[8].compression = PNG_TEXT_COMPRESSION_NONE;
    text_ptr[9].key = "width";
    text_ptr[9].text = const_cast<char*>(std::to_string(width).c_str());
    text_ptr[9].compression = PNG_TEXT_COMPRESSION_NONE;
    png_set_text(png_ptr, info_ptr, text_ptr, num_text);

    if (isColor) {
        // color palette
        png_set_IHDR(png_ptr, info_ptr, width, height, 8, PNG_COLOR_TYPE_PALETTE, PNG_INTERLACE_NONE,
                     PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);

        uint8_t pal[256 * 3];
        for (int i = 0; i < 256; i++) {
            pal[i * 3] = red[i];
            pal[i * 3 + 1] = green[i];
            pal[i * 3 + 2] = blue[i];
        }
        png_set_PLTE(png_ptr, info_ptr, (png_const_colorp)pal, 256);

        if (transparent) {
            uint8_t transPal[256];
            for (int i = 0; i < 255; i++) {
                transPal[i] = 255;
            }
            transPal[255] = 0;
            png_set_tRNS(png_ptr, info_ptr, transPal, 256, NULL);
        }
    } else {
        // Grayscale
        png_set_IHDR(png_ptr, info_ptr, width, height, 8, PNG_COLOR_TYPE_GRAY, PNG_INTERLACE_NONE,
                     PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);
    }
    png_write_info(png_ptr, info_ptr);

    return true;
}

void OutFile_png::writeLine() {
    for (int i = 0; i < width; i++)
        fileData[i] = (uint8_t)round(productStuff[0]->calcOutputVal(productStuff[0]->lineData[i]));
    png_write_row(png_ptr, (png_bytep)fileData);
    currentLine++;
}

bool OutFile_png::close() {
    png_write_end(png_ptr, info_ptr);
    png_destroy_write_struct(&png_ptr, &info_ptr);

    fclose(outfp);
    outfp = NULL;
    return true;
}

//------------------------------------------------------------------------------
// OutFile_png_rgb
//------------------------------------------------------------------------------

OutFile_png_rgb::OutFile_png_rgb() : OutFile() {
    outfp = NULL;
    fileData = NULL;
    info_ptr = NULL;
    png_ptr = NULL;
    num_text = 10;
}

OutFile_png_rgb::~OutFile_png_rgb() {
    if (outfp)
        fclose(outfp);
    if (fileData)
        free(fileData);
}

void OutFile_png_rgb::setSize(int32_t width, int32_t height) {
    OutFile::setSize(width, height);
    if (fileData)
        free(fileData);
    int samplesPerPixel = 3;
    if (transparent)
        samplesPerPixel = 4;
    fileData = (uint8_t*)allocateMemory(width * samplesPerPixel, "OutFile_png_rgb::fileData");
}

bool OutFile_png_rgb::open() {
    outfp = fopen(fileName.c_str(), "w");
    if (!outfp)
        return false;

    png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if (!png_ptr) {
        fprintf(stderr, "-E- Unable to create PNG write structure.\n");
        exit(EXIT_FAILURE);
    }

    info_ptr = png_create_info_struct(png_ptr);
    if (!info_ptr) {
        png_destroy_write_struct(&png_ptr, (png_infopp)NULL);
        fprintf(stderr, "-E- Unable to create PNG info structure.\n");
        exit(EXIT_FAILURE);
    }
    if (setjmp(png_jmpbuf(png_ptr))) {
        png_destroy_write_struct(&png_ptr, (png_infopp)NULL);
        fprintf(stderr, "-E- Unable to call PNG setjmp().\n");
        exit(EXIT_FAILURE);
    }
    png_init_io(png_ptr, outfp);

    png_text text_ptr[num_text];
    text_ptr[0].key = "projString";
    text_ptr[0].text = const_cast<char*>(proj4String.c_str());
    text_ptr[0].compression = PNG_TEXT_COMPRESSION_NONE;
    text_ptr[1].key = "resolution";
    text_ptr[1].text = const_cast<char*>(std::to_string(resolution).c_str());
    text_ptr[1].compression = PNG_TEXT_COMPRESSION_NONE;
    text_ptr[2].key = "north";
    text_ptr[2].text = const_cast<char*>(std::to_string(metaData->north).c_str());
    text_ptr[2].compression = PNG_TEXT_COMPRESSION_NONE;
    text_ptr[3].key = "south";
    text_ptr[3].text = const_cast<char*>(std::to_string(metaData->south).c_str());
    text_ptr[3].compression = PNG_TEXT_COMPRESSION_NONE;
    text_ptr[4].key = "east";
    text_ptr[4].text = const_cast<char*>(std::to_string(metaData->east).c_str());
    text_ptr[4].compression = PNG_TEXT_COMPRESSION_NONE;
    text_ptr[5].key = "west";
    text_ptr[5].text = const_cast<char*>(std::to_string(metaData->west).c_str());
    text_ptr[5].compression = PNG_TEXT_COMPRESSION_NONE;
    text_ptr[6].key = "minX";
    text_ptr[6].text = const_cast<char*>(std::to_string(tiepoints[3]).c_str());
    text_ptr[6].compression = PNG_TEXT_COMPRESSION_NONE;
    text_ptr[7].key = "maxY";
    text_ptr[7].text = const_cast<char*>(std::to_string(tiepoints[4]).c_str());
    text_ptr[7].compression = PNG_TEXT_COMPRESSION_NONE;
    text_ptr[8].key = "height";
    text_ptr[8].text = const_cast<char*>(std::to_string(height).c_str());
    text_ptr[8].compression = PNG_TEXT_COMPRESSION_NONE;
    text_ptr[9].key = "width";
    text_ptr[9].text = const_cast<char*>(std::to_string(width).c_str());
    text_ptr[9].compression = PNG_TEXT_COMPRESSION_NONE;
    png_set_text(png_ptr, info_ptr, text_ptr, num_text);

    if (transparent)
        png_set_IHDR(png_ptr, info_ptr, width, height, 8, PNG_COLOR_TYPE_RGBA, PNG_INTERLACE_NONE,
                     PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);
    else
        png_set_IHDR(png_ptr, info_ptr, width, height, 8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
                     PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);
    png_write_info(png_ptr, info_ptr);

    return true;
}

bool OutFile_png_rgb::close() {
    if (setjmp(png_jmpbuf(png_ptr))) {
        fprintf(stderr, "-E- OutFile_png_rgb::close - Unable to call PNG setjmp().\n");
        exit(EXIT_FAILURE);
    }
    png_write_end(png_ptr, info_ptr);
    png_destroy_write_struct(&png_ptr, &info_ptr);

    fclose(outfp);
    outfp = NULL;
    return true;
}

void OutFile_png_rgb::setPixel(int32_t x, double val, int32_t prod) {
    fprintf(stderr,
            "-E- OutFile_png_rgb::setPixel - only RGB is implemented with this file type.\n");
    exit(EXIT_FAILURE);
}

void OutFile_png_rgb::setPixelRGB(int32_t x, float red, float green, float blue) {
    if (x < 0 || x >= width) {
        fprintf(stderr,
                "-E- OutFile_png_rgb::setPixel - x=%d is not within range, width=%d.\n",
                x, width);
        exit(EXIT_FAILURE);
    }

    int samplesPerPixel = 3;
    if (transparent)
        samplesPerPixel = 4;
    uint8_t* ptr = fileData + x * samplesPerPixel;

    *ptr = round(productStuff[0]->calcOutputVal(red));
    ptr++;
    *ptr = round(productStuff[0]->calcOutputVal(green));
    ptr++;
    *ptr = round(productStuff[0]->calcOutputVal(blue));
    if (transparent) {
        ptr++;
        if (red == badPixelValue || green == badPixelValue || blue == badPixelValue) {
            *ptr = 0;  // transparent
        } else {
            *ptr = 255;
        }
    }

    // do this to keep the file min/max reasonable
    if (red > fileMaxVal)
        fileMaxVal = red;
    if (red < fileMinVal)
        fileMinVal = red;
}

void OutFile_png_rgb::landPixel(int32_t x) {
    if (x < 0 || x >= width) {
        fprintf(stderr,
                "-E- OutFile_png_rgb::landPixel - x=%d is not within range, width=%d.\n",
                x, width);
        exit(EXIT_FAILURE);
    }
    int samplesPerPixel = 3;
    if (transparent)
        samplesPerPixel = 4;
    uint8_t* ptr = fileData + x * samplesPerPixel;
    *ptr = landPix;
    ptr++;
    *ptr = landPix;
    ptr++;
    *ptr = landPix;
    if (transparent) {
        ptr++;
        *ptr = 255;
    }
}

void OutFile_png_rgb::fillPixel(int32_t x) {
    if (x < 0 || x >= width) {
        fprintf(stderr,
                "-E- OutFile_png_rgb::fillPixel - x=%d is not within range, width=%d.\n",
                x, width);
        exit(EXIT_FAILURE);
    }
    int samplesPerPixel = 3;
    if (transparent)
        samplesPerPixel = 4;
    uint8_t* ptr = fileData + x * samplesPerPixel;
    *ptr = fillPix;
    ptr++;
    *ptr = fillPix;
    ptr++;
    *ptr = fillPix;
    if (transparent) {
        ptr++;
        *ptr = 0;
    }
}

void OutFile_png_rgb::missingPixel(int32_t x) {
    if (x < 0 || x >= width) {
        fprintf(stderr,
                "-E- OutFile_png_rgb::missingPixel - x=%d is not within range, width=%d.\n",
                x, width);
        exit(EXIT_FAILURE);
    }
    int samplesPerPixel = 3;
    if (transparent)
        samplesPerPixel = 4;
    uint8_t* ptr = fileData + x * samplesPerPixel;
    *ptr = fillPix;
    ptr++;
    *ptr = fillPix;
    ptr++;
    *ptr = fillPix;
    if (transparent) {
        ptr++;
        *ptr = 0;
    }
}

void OutFile_png_rgb::writeLine() {
    if (setjmp(png_jmpbuf(png_ptr))) {
        fprintf(stderr, "-E- OutFile_png_rgb::writeLine - Unable to call PNG setjmp().\n");
        exit(EXIT_FAILURE);
    }

    png_write_row(png_ptr, (png_bytep)fileData);
    currentLine++;
}

//------------------------------------------------------------------------------
// OutFile_tiff
//------------------------------------------------------------------------------

OutFile_tiff::~OutFile_tiff() {
    close();
}

bool OutFile_tiff::open() {
    currentLine = 0;

    // open TIFF file
    tiff = XTIFFOpen(fileName.c_str(), "w");
    if (tiff == NULL) {
        cerr << "-E- Could not open outgoing TIFF image" << endl;
        exit(EXIT_FAILURE);
    }

    // extend the TIFF tags to write pversion
    ttag_t TIFFTAG_GDALMetadata = 42112;
    static const TIFFFieldInfo xtiffFieldInfo[] = {
        {TIFFTAG_GDALMetadata, -1, -1, TIFF_ASCII, FIELD_CUSTOM, TRUE, FALSE, "GDALMetadata"}};
    TIFFMergeFieldInfo(tiff, xtiffFieldInfo, 1);
    std::string tagVal = "<GDALMetadata>\n  <Item name=\"OBPG_version\">";
    tagVal += metaData->pversion;
    tagVal += "</Item>\n";
    tagVal += "</GDALMetadata>\n";
    TIFFSetField(tiff, TIFFTAG_GDALMetadata, tagVal.c_str());

    TIFFSetField(tiff, TIFFTAG_IMAGEWIDTH, width);
    TIFFSetField(tiff, TIFFTAG_IMAGELENGTH, height);
    TIFFSetField(tiff, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    TIFFSetField(tiff, TIFFTAG_COMPRESSION, COMPRESSION_LZW);
    TIFFSetField(tiff, TIFFTAG_PREDICTOR, PREDICTOR_NONE);
    TIFFSetField(tiff, TIFFTAG_ROWSPERSTRIP, height);

    // open GeoTIFF interface
    bool hasGeotiffInfo = false;
    gtif = GTIFNew(tiff);
    if (gtif == NULL) {
        cerr << "-E- Could not create geoTIFF structure" << endl;
        exit(EXIT_FAILURE);
    }

    // define GeoTIFF keys for lat/lon projection
    if (mapProjection == "Equidistant Cylindrical" || mapProjection == "PlateCarree") {
        double north, south, east, west;
        if (metaData) {
            north = metaData->north;
            south = metaData->south;
            east = metaData->east;
            west = metaData->west;
        } else {
            north = 90.0;
            south = -90.0;
            east = 180.0;
            west = -180.0;
        }

        // pixel width, height in degrees
        pixscale[0] = (east - west) / width;
        pixscale[1] = (north - south) / height;

        // top left corner pixel lat, lon
        for (int i = 0; i < 6; i++)
            tiepoints[i] = 0;
        tiepoints[3] = west;
        tiepoints[4] = north;

        GTIFKeySet(gtif, GTModelTypeGeoKey, TYPE_SHORT, 1, ModelGeographic);
        GTIFKeySet(gtif, GeographicTypeGeoKey, TYPE_SHORT, 1, GCS_WGS_84);
        hasGeotiffInfo = true;
    }
    // otherwise, parse the proj4 string
    else
        hasGeotiffInfo = (GTIFSetFromProj4(gtif, proj4String.c_str()));

    if (!hasGeotiffInfo) {
        if (proj4String.find("EPSG:") != string::npos)
            hasGeotiffInfo = true;
    }

    if (hasGeotiffInfo) {
        // write GeoTIFF keys
        GTIFKeySet(gtif, GTRasterTypeGeoKey, TYPE_SHORT, 1, RasterPixelIsArea);
        static const std::string& epsg_regex_str{"EPSG:(\\d+)\\b"};
        static const std::regex epsg_regex{epsg_regex_str};
        std::smatch matches;
        if (std::regex_search(proj4String, matches, epsg_regex)) {
            unsigned short epsg = std::stoi(matches[1]);
            GTIFKeySet(gtif, ProjectedCSTypeGeoKey, TYPE_SHORT, 1, epsg);
        }
        GTIFWriteKeys(gtif);

        // write GeoTIFF tags in map units
        TIFFSetField(tiff, GTIFF_PIXELSCALE, 3, pixscale);
        TIFFSetField(tiff, GTIFF_TIEPOINTS, 6, tiepoints);
    }

    // define colormap
    setTiffColor();
    return true;
}

bool OutFile_tiff::close() {
    if (gtif) {
        GTIFFree(gtif);
        XTIFFClose(tiff);
        tiff = NULL;
        gtif = NULL;
    }
    return true;
}

//----- OutFile_tiff_color -----

OutFile_tiff_color::~OutFile_tiff_color() {
    if (fileData)
        free(fileData);
}

void OutFile_tiff_color::setSize(int32_t width, int32_t height) {
    OutFile::setSize(width, height);
    if (fileData)
        free(fileData);

    // if transparent we need to use RGBA mode instead of indexed color map
    if (transparent) {
        fileData = (uint8_t*)allocateMemory(width * 4, "OutFile_tiff_color::fileData");
    } else {
        fileData = (uint8_t*)allocateMemory(width, "OutFile_tiff_color::fileData");
    }
}

void OutFile_tiff_color::writeLine() {
    if (transparent) {
        for (int i = 0; i < width; i++) {
            uint8_t* ptr = fileData + i * 4;

            uint8_t alpha = 255;  // opaque

            uint8_t scaledVal = (uint8_t)round(productStuff[0]->calcOutputVal(productStuff[0]->lineData[i]));
            if (productStuff[0]->lineData[i] == badPixelValue) {
                alpha = 0;  // transparent
            }
            *ptr = red[scaledVal];
            ptr++;
            *ptr = green[scaledVal];
            ptr++;
            *ptr = blue[scaledVal];
            ptr++;
            *ptr = alpha;
        }
    } else {
        for (int i = 0; i < width; i++)
            fileData[i] = (uint8_t)round(productStuff[0]->calcOutputVal(productStuff[0]->lineData[i]));
    }

    if (TIFFWriteScanline(tiff, (void*)fileData, currentLine) == 0) {
        cerr << "-E- Could not write TIFF image line " << currentLine << endl;
        exit(EXIT_FAILURE);
    }

    currentLine++;
}

void OutFile_tiff_color::setTiffColor() {
    if (transparent) {
        TIFFSetField(tiff, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
        TIFFSetField(tiff, TIFFTAG_BITSPERSAMPLE, 8);
        TIFFSetField(tiff, TIFFTAG_SAMPLESPERPIXEL, 4);
        uint16_t out[1];
        out[0] = EXTRASAMPLE_ASSOCALPHA;
        TIFFSetField(tiff, TIFFTAG_EXTRASAMPLES, 1, &out);
        TIFFSetField(tiff, TIFFTAG_COMPRESSION, COMPRESSION_ADOBE_DEFLATE);
    } else {
        TIFFSetField(tiff, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_PALETTE);
        TIFFSetField(tiff, TIFFTAG_BITSPERSAMPLE, 8);
        TIFFSetField(tiff, TIFFTAG_SAMPLESPERPIXEL, 1);

        int ncolors = 256;  // convert bytes to short
        uint16_t* rr = (uint16_t*)malloc(ncolors * sizeof(uint16_t));
        uint16_t* gg = (uint16_t*)malloc(ncolors * sizeof(uint16_t));
        uint16_t* bb = (uint16_t*)malloc(ncolors * sizeof(uint16_t));
        if (rr == NULL || gg == NULL || bb == NULL) {
            cerr << "-E- Could not allocate memory for TIFF color map" << endl;
            exit(EXIT_FAILURE);
        }
        for (int i = 0; i < ncolors; i++) {
            rr[i] = red[i] << 8;
            gg[i] = green[i] << 8;
            bb[i] = blue[i] << 8;
        }
        TIFFSetField(tiff, TIFFTAG_COLORMAP, rr, gg, bb);
        free(rr);
        free(gg);
        free(bb);
    }
}

//----- OutFile_tiff_gray -----

OutFile_tiff_gray::OutFile_tiff_gray() {
    this->landPixelValue = badPixelValue; // Grayscale outputs a floating point value
}

OutFile_tiff_gray::~OutFile_tiff_gray() {
    if (fileData)
        free(fileData);
}

void OutFile_tiff_gray::setSize(int32_t width, int32_t height) {
    OutFile::setSize(width, height);
    if (fileData)
        free(fileData);
    fileData = (float*)allocateMemory(width * sizeof(float), "OutFile_tiff_gray::fileData");
}

void OutFile_tiff_gray::writeLine() {
    for (int i = 0; i < width; i++)
        fileData[i] = productStuff[0]->lineData[i];

    if (TIFFWriteScanline(tiff, (void*)fileData, currentLine) == 0) {
        cerr << "-E- Could not write TIFF image line " << currentLine << endl;
        exit(EXIT_FAILURE);
    }
    currentLine++;
}

void OutFile_tiff_gray::setTiffColor() {
    TIFFSetField(tiff, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
    TIFFSetField(tiff, TIFFTAG_BITSPERSAMPLE, 32);
    TIFFSetField(tiff, TIFFTAG_SAMPLESPERPIXEL, 1);
    TIFFSetField(tiff, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP);
}

//----- OutFile_tiff_rgb -----

OutFile_tiff_rgb::~OutFile_tiff_rgb() {
    if (fileData)
        free(fileData);
}

void OutFile_tiff_rgb::setSize(int32_t width, int32_t height) {
    OutFile::setSize(width, height);
    if (fileData)
        free(fileData);
    if (transparent) {
        fileData = (uint8_t*)allocateMemory(width * 4, "OutFile_tiff_rgb::fileData");
    } else {
        fileData = (uint8_t*)allocateMemory(width * 3, "OutFile_tiff_rgb::fileData");
    }
}

void OutFile_tiff_rgb::writeLine() {
    if (TIFFWriteScanline(tiff, (void*)fileData, currentLine) == 0) {
        cerr << "-E- Could not write TIFF image line " << currentLine << endl;
        exit(EXIT_FAILURE);
    }
    currentLine++;
}

void OutFile_tiff_rgb::setTiffColor() {
    TIFFSetField(tiff, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
    TIFFSetField(tiff, TIFFTAG_BITSPERSAMPLE, 8);
    if (transparent) {
        TIFFSetField(tiff, TIFFTAG_SAMPLESPERPIXEL, 4);
    } else {
        TIFFSetField(tiff, TIFFTAG_SAMPLESPERPIXEL, 3);
    }
    TIFFSetField(tiff, TIFFTAG_COMPRESSION, COMPRESSION_ADOBE_DEFLATE);
}

void OutFile_tiff_rgb::setPixel(int32_t x, double val, int32_t prod) {
    if (x < 0 || x >= width) {
        fprintf(stderr,
                "-E- OutFile_tiff_rgb::setPixel - x=%d is not within range, width=%d.\n",
                x, width);
        exit(EXIT_FAILURE);
    }

    uint8_t samplesPerPixel = 3;
    if (transparent) {
        samplesPerPixel = 4;
    }
    uint8_t* ptr = fileData + x * samplesPerPixel;

    uint8_t alpha = 255;  // opaque

    uint8_t scaledVal = round(productStuff[0]->calcOutputVal(val));
    if (val == badPixelValue) {
        alpha = 0;  // transparent
    }
    *ptr = red[scaledVal];
    ptr++;
    *ptr = green[scaledVal];
    ptr++;
    *ptr = blue[scaledVal];
    if (transparent) {
        ptr++;
        *ptr = alpha;
    }
}

void OutFile_tiff_rgb::setPixelRGB(int32_t x, float red, float green, float blue) {
    if (x < 0 || x >= width) {
        fprintf(stderr,
                "-E- OutFile_tiff_rgb::setPixel - x=%d is not within range, width=%d.\n",
                x, width);
        exit(EXIT_FAILURE);
    }

    uint8_t alpha = 255;  // opaque

    uint8_t samplesPerPixel = 3;
    if (transparent) {
        samplesPerPixel = 4;
    }
    uint8_t* ptr = fileData + x * samplesPerPixel;

    *ptr = round(productStuff[0]->calcOutputVal(red));
    ptr++;
    *ptr = round(productStuff[0]->calcOutputVal(green));
    ptr++;
    *ptr = round(productStuff[0]->calcOutputVal(blue));
    if (transparent) {
        if (red == badPixelValue || green == badPixelValue || blue == badPixelValue) {
            alpha = 0;  // transparent
        }
        ptr++;
        *ptr = alpha;
    }

    if (red > fileMaxVal)  // keep file min/max reasonable
        fileMaxVal = red;
    if (red < fileMinVal)
        fileMinVal = red;
}

void OutFile_tiff_rgb::landPixel(int32_t x) {
    if (x < 0 || x >= width) {
        fprintf(stderr,
                "-E- OutFile_tiff_rgb::landPixel - x=%d is not within range, width=%d.\n",
                x, width);
        exit(EXIT_FAILURE);
    }

    uint8_t samplesPerPixel = 3;
    if (transparent) {
        samplesPerPixel = 4;
    }
    uint8_t* ptr = fileData + x * samplesPerPixel;
    *ptr = landPix;
    ptr++;
    *ptr = landPix;
    ptr++;
    *ptr = landPix;
    if (transparent) {
        ptr++;
        *ptr = 255;
    }
}

void OutFile_tiff_rgb::fillPixel(int32_t x) {
    if (x < 0 || x >= width) {
        fprintf(stderr,
                "-E- OutFile_tiff_rgb::fillPixel - x=%d is not within range, width=%d.\n",
                x, width);
        exit(EXIT_FAILURE);
    }

    uint8_t samplesPerPixel = 3;
    if (transparent) {
        samplesPerPixel = 4;
    }
    uint8_t* ptr = fileData + x * samplesPerPixel;
    *ptr = fillPix;
    ptr++;
    *ptr = fillPix;
    ptr++;
    *ptr = fillPix;
    if (transparent) {
        ptr++;
        *ptr = 0;
    }
}

void OutFile_tiff_rgb::missingPixel(int32_t x) {
    if (x < 0 || x >= width) {
        fprintf(stderr,
                "-E- OutFile_tiff_rgb::missingPixel - x=%d is not within range, width=%d.\n",
                x, width);
        exit(EXIT_FAILURE);
    }

    uint8_t samplesPerPixel = 3;
    if (transparent) {
        samplesPerPixel = 4;
    }
    uint8_t* ptr = fileData + x * samplesPerPixel;
    *ptr = fillPix;
    ptr++;
    *ptr = fillPix;
    ptr++;
    *ptr = fillPix;
    if (transparent) {
        ptr++;
        *ptr = 0;
    }
}

//------------------------------------------------------------------------------
// OutFile_hdf4
//------------------------------------------------------------------------------

OutFile_hdf4::OutFile_hdf4() : OutFile() {
    fileData = NULL;
    sdfid = -1;
    sdsid = -1;
    quality_sdsid = -1;
    hdfDataType = DFNT_FLOAT32;
}

OutFile_hdf4::~OutFile_hdf4() {
    if (fileData)
        free(fileData);
}

void OutFile_hdf4::setSize(int32_t width, int32_t height) {
    OutFile::setSize(width, height);
    if (fileData)
        free(fileData);
    fileData = NULL;
}

bool OutFile_hdf4::open() {
    const char* tmpStr;
    float tmpFloat;
    int32_t tmpInt32;

    currentLine = 0;

    sdfid = SDstart(fileName.c_str(), DFACC_CREATE);

    if (sdfid < 0) {
        printf("-E- Could not create HDF4 file %s\n", fileName.c_str());
        exit(EXIT_FAILURE);
    }

    get_time(metaData->ptime);

    string prodName;
    size_t pos = fileName.find_last_of('/');
    if (pos == string::npos)
        prodName = fileName;
    else
        prodName = fileName.substr(pos + 1);
    DPTB(SDsetattr(sdfid, "Product Name", DFNT_CHAR, prodName.size() + 1, (VOIDP)prodName.c_str()));
    DPTB(SDsetattr(sdfid, "Sensor Name", DFNT_CHAR, strlen(metaData->sensor_name) + 1,
                   (VOIDP)metaData->sensor_name));
    DPTB(SDsetattr(sdfid, "Sensor", DFNT_CHAR, strlen(metaData->sensor) + 1, (VOIDP)metaData->sensor));
    DPTB(SDsetattr(sdfid, "Title", DFNT_CHAR, strlen(metaData->title) + 1, (VOIDP)metaData->title));
    DPTB(SDsetattr(sdfid, "Data Center", DFNT_CHAR, strlen(metaData->data_center) + 1,
                   (VOIDP)metaData->data_center));
    DPTB(
        SDsetattr(sdfid, "Station Name", DFNT_CHAR, strlen(metaData->station) + 1, (VOIDP)metaData->station));
    DPTB(SDsetattr(sdfid, "Station Latitude", DFNT_FLOAT32, 1, (VOIDP)&metaData->station_lat));
    DPTB(SDsetattr(sdfid, "Station Longitude", DFNT_FLOAT32, 1, (VOIDP)&metaData->station_lon));
    DPTB(SDsetattr(sdfid, "Mission", DFNT_CHAR, strlen(metaData->mission) + 1, (VOIDP)metaData->mission));
    DPTB(SDsetattr(sdfid, "Mission Characteristics", DFNT_CHAR, strlen(metaData->mission_char) + 1,
                   (VOIDP)metaData->mission_char));
    DPTB(SDsetattr(sdfid, "Sensor Characteristics", DFNT_CHAR, strlen(metaData->sensor_char) + 1,
                   (VOIDP)metaData->sensor_char));
    DPTB(SDsetattr(sdfid, "Product Type", DFNT_CHAR, strlen(metaData->prod_type) + 1,
                   (VOIDP)metaData->prod_type));
    DPTB(SDsetattr(sdfid, "Processing Version", DFNT_CHAR, strlen(metaData->pversion) + 1,
                   (VOIDP)metaData->pversion));
    DPTB(SDsetattr(sdfid, "Software Name", DFNT_CHAR, strlen(metaData->soft_name) + 1,
                   (VOIDP)metaData->soft_name));
    DPTB(SDsetattr(sdfid, "Software Version", DFNT_CHAR, strlen(metaData->soft_ver) + 1,
                   (VOIDP)metaData->soft_ver));
    DPTB(SDsetattr(sdfid, "Processing Time", DFNT_CHAR, strlen(metaData->ptime) + 1, (VOIDP)metaData->ptime));
    DPTB(SDsetattr(sdfid, "Input Files", DFNT_CHAR, strlen(metaData->infiles) + 1, (VOIDP)metaData->infiles));
    DPTB(SDsetattr(sdfid, "Processing Control", DFNT_CHAR, strlen(metaData->proc_con) + 1,
                   (VOIDP)metaData->proc_con));
    DPTB(SDsetattr(sdfid, "Input Parameters", DFNT_CHAR, strlen(metaData->input_parms) + 1,
                   (VOIDP)metaData->input_parms));
    DPTB(SDsetattr(sdfid, "L2 Flag Names", DFNT_CHAR, strlen(metaData->flag_names) + 1,
                   (VOIDP)metaData->flag_names));

    short syear, sday, eyear, eday;
    double ssec, esec;
    int32_t smsec, emsec;
    unix2yds(metaData->startTime, &syear, &sday, &ssec);
    smsec = (int32_t)(ssec * 1000.0);
    unix2yds(metaData->endTime, &eyear, &eday, &esec);
    emsec = (int32_t)(esec * 1000.0);
    DPTB(SDsetattr(sdfid, "Period Start Year", DFNT_INT16, 1, (VOIDP)&syear));
    DPTB(SDsetattr(sdfid, "Period Start Day", DFNT_INT16, 1, (VOIDP)&sday));
    DPTB(SDsetattr(sdfid, "Period End Year", DFNT_INT16, 1, (VOIDP)&eyear));
    DPTB(SDsetattr(sdfid, "Period End Day", DFNT_INT16, 1, (VOIDP)&eday));
    tmpStr = ydhmsf(metaData->startTime, 'G');
    DPTB(SDsetattr(sdfid, "Start Time", DFNT_CHAR, strlen(tmpStr) + 1, (VOIDP)tmpStr));
    tmpStr = ydhmsf(metaData->endTime, 'G');
    DPTB(SDsetattr(sdfid, "End Time", DFNT_CHAR, strlen(tmpStr) + 1, (VOIDP)tmpStr));
    DPTB(SDsetattr(sdfid, "Start Year", DFNT_INT16, 1, (VOIDP)&syear));
    DPTB(SDsetattr(sdfid, "Start Day", DFNT_INT16, 1, (VOIDP)&sday));
    DPTB(SDsetattr(sdfid, "Start Millisec", DFNT_INT32, 1, (VOIDP)&smsec));
    DPTB(SDsetattr(sdfid, "End Year", DFNT_INT16, 1, (VOIDP)&eyear));
    DPTB(SDsetattr(sdfid, "End Day", DFNT_INT16, 1, (VOIDP)&eday));
    DPTB(SDsetattr(sdfid, "End Millisec", DFNT_INT32, 1, (VOIDP)&emsec));

    DPTB(SDsetattr(sdfid, "Start Orbit", DFNT_INT32, 1, (VOIDP)&metaData->start_orb));
    DPTB(SDsetattr(sdfid, "End Orbit", DFNT_INT32, 1, (VOIDP)&metaData->end_orb));
    DPTB(SDsetattr(sdfid, "Orbit", DFNT_INT32, 1, (VOIDP)&metaData->orbit));
    DPTB(SDsetattr(sdfid, "Map Projection", DFNT_CHAR, mapProjection.size() + 1, mapProjection.c_str()));
    DPTB(SDsetattr(sdfid, "Latitude Units", DFNT_CHAR, 14, (VOIDP) "degrees North"));
    DPTB(SDsetattr(sdfid, "Longitude Units", DFNT_CHAR, 13, (VOIDP) "degrees East"));
    DPTB(SDsetattr(sdfid, "Northernmost Latitude", DFNT_FLOAT32, 1, (VOIDP)&metaData->north));
    DPTB(SDsetattr(sdfid, "Southernmost Latitude", DFNT_FLOAT32, 1, (VOIDP)&metaData->south));
    DPTB(SDsetattr(sdfid, "Westernmost Longitude", DFNT_FLOAT32, 1, (VOIDP)&metaData->west));
    DPTB(SDsetattr(sdfid, "Easternmost Longitude", DFNT_FLOAT32, 1, (VOIDP)&metaData->east));

    float latStep = (metaData->north - metaData->south) / (float)getHeight();
    DPTB(SDsetattr(sdfid, "Latitude Step", DFNT_FLOAT32, 1, (VOIDP)&latStep));
    float lonStep;
    if (metaData->east < metaData->west)
        lonStep = (360 + metaData->east - metaData->west) / (float)getWidth();
    else
        lonStep = (metaData->east - metaData->west) / (float)getWidth();
    DPTB(SDsetattr(sdfid, "Longitude Step", DFNT_FLOAT32, 1, (VOIDP)&lonStep));

    tmpFloat = metaData->south + latStep / 2.0;
    DPTB(SDsetattr(sdfid, "SW Point Latitude", DFNT_FLOAT32, 1, (VOIDP)&tmpFloat));
    tmpFloat = metaData->west + lonStep / 2.0;
    DPTB(SDsetattr(sdfid, "SW Point Longitude", DFNT_FLOAT32, 1, (VOIDP)&tmpFloat));

    tmpInt32 = metaData->data_bins;
    DPTB(SDsetattr(sdfid, "Data Bins", DFNT_INT32, 1, (VOIDP)&tmpInt32));

    tmpInt32 = getHeight();
    DPTB(SDsetattr(sdfid, "Number of Lines", DFNT_INT32, 1, (VOIDP)&tmpInt32));
    tmpInt32 = getWidth();
    DPTB(SDsetattr(sdfid, "Number of Columns", DFNT_INT32, 1, (VOIDP)&tmpInt32));
    DPTB(SDsetattr(sdfid, "Parameter", DFNT_CHAR, strlen(productStuff[0]->productInfo->description) + 1,
                   (VOIDP)productStuff[0]->productInfo->description));
    DPTB(SDsetattr(sdfid, "Measure", DFNT_CHAR, 5, (VOIDP) "Mean"));
    DPTB(SDsetattr(sdfid, "Units", DFNT_CHAR, strlen(productStuff[0]->productInfo->units) + 1,
                   (VOIDP)productStuff[0]->productInfo->units));

    // we only use linear scaling for data storage
    tmpStr = "linear";
    DPTB(SDsetattr(sdfid, "Scaling", DFNT_CHAR, strlen(tmpStr) + 1, (VOIDP)tmpStr));
    tmpStr = "(Slope*l3m_data) + Intercept = Parameter value";
    DPTB(SDsetattr(sdfid, "Scaling Equation", DFNT_CHAR, strlen(tmpStr) + 1, (VOIDP)tmpStr));

    float tmpScale;
    float tmpOffset;
    const char* imageScalingApplied;

    switch (productStuff[0]->dataStorage) {
        case FloatDS:
        case DoubleDS:
            tmpScale = 1.0;
            tmpOffset = 0.0;
            imageScalingApplied = "No";
            break;
        default:
            tmpScale = productStuff[0]->scale;
            tmpOffset = productStuff[0]->offset;
            if (tmpScale == 1.0 && tmpOffset == 0.0)
                imageScalingApplied = "No";
            else
                imageScalingApplied = "Yes";
    }

    DPTB(SDsetattr(sdfid, "Slope", DFNT_FLOAT32, 1, (VOIDP)&tmpScale));
    DPTB(SDsetattr(sdfid, "Intercept", DFNT_FLOAT32, 1, (VOIDP)&tmpOffset));

    tmpFloat = productStuff[0]->productInfo->validMin;
    DPTB(SDsetattr(sdfid, "Data Minimum", DFNT_FLOAT32, 1, (VOIDP)&tmpFloat));
    tmpFloat = productStuff[0]->productInfo->validMax;
    DPTB(SDsetattr(sdfid, "Data Maximum", DFNT_FLOAT32, 1, (VOIDP)&tmpFloat));
    tmpFloat = productStuff[0]->productInfo->displayMin;
    DPTB(SDsetattr(sdfid, "Suggested Image Scaling Minimum", DFNT_FLOAT32, 1, (VOIDP)&tmpFloat));
    tmpFloat = productStuff[0]->productInfo->displayMax;
    DPTB(SDsetattr(sdfid, "Suggested Image Scaling Maximum", DFNT_FLOAT32, 1, (VOIDP)&tmpFloat));

    tmpStr = strdup(productStuff[0]->productInfo->displayScale);
    upcase((char*)tmpStr);
    DPTB(SDsetattr(sdfid, "Suggested Image Scaling Type", DFNT_CHAR, strlen(tmpStr) + 1, (VOIDP)tmpStr));
    free((void*)tmpStr);

    DPTB(SDsetattr(sdfid, "Suggested Image Scaling Applied", DFNT_CHAR, strlen(imageScalingApplied) + 1,
                   (VOIDP)imageScalingApplied));
    DPTB(SDsetattr(sdfid, "_lastModified", DFNT_CHAR, strlen(metaData->ptime) + 1, (VOIDP)metaData->ptime));

    // delete file data
    if (fileData) {
        free(fileData);
        fileData = NULL;
    }

    int32_t dims[2];
    dims[0] = height;
    dims[1] = width;

    if (!strcmp(productStuff[0]->productInfo->dataType, "byte")) {
        hdfDataType = DFNT_INT8;
        sdsid = SDcreate(sdfid, "l3m_data", hdfDataType, 2, dims);
        int8_t tmp = productStuff[0]->missingValue;
        DPTB(SDsetattr(sdsid, "Fill", hdfDataType, 1, (VOIDP)&tmp));
        fileData = allocateMemory(width * sizeof(int8_t), "OutFile_hdf4::open fileData");
    } else if (!strcmp(productStuff[0]->productInfo->dataType, "ubyte")) {
        hdfDataType = DFNT_UINT8;
        sdsid = SDcreate(sdfid, "l3m_data", hdfDataType, 2, dims);
        uint8_t tmp = productStuff[0]->missingValue;
        DPTB(SDsetattr(sdsid, "Fill", hdfDataType, 1, (VOIDP)&tmp));
        fileData = allocateMemory(width * sizeof(uint8_t), "OutFile_hdf4::open fileData");
    } else if (!strcmp(productStuff[0]->productInfo->dataType, "short")) {
        hdfDataType = DFNT_INT16;
        sdsid = SDcreate(sdfid, "l3m_data", hdfDataType, 2, dims);
        int16_t tmp = productStuff[0]->missingValue;
        DPTB(SDsetattr(sdsid, "Fill", hdfDataType, 1, (VOIDP)&tmp));
        fileData = allocateMemory(width * sizeof(int16_t), "OutFile_hdf4::open fileData");
    } else if (!strcmp(productStuff[0]->productInfo->dataType, "ushort")) {
        hdfDataType = DFNT_UINT16;
        sdsid = SDcreate(sdfid, "l3m_data", hdfDataType, 2, dims);
        uint16_t tmp = productStuff[0]->missingValue;
        DPTB(SDsetattr(sdsid, "Fill", hdfDataType, 1, (VOIDP)&tmp));
        fileData = allocateMemory(width * sizeof(uint16_t), "OutFile_hdf4::open fileData");
    } else if (!strcmp(productStuff[0]->productInfo->dataType, "int")) {
        hdfDataType = DFNT_INT32;
        sdsid = SDcreate(sdfid, "l3m_data", hdfDataType, 2, dims);
        int32_t tmp = productStuff[0]->missingValue;
        DPTB(SDsetattr(sdsid, "Fill", hdfDataType, 1, (VOIDP)&tmp));
        fileData = allocateMemory(width * sizeof(int32_t), "OutFile_hdf4::open fileData");
    } else if (!strcmp(productStuff[0]->productInfo->dataType, "uint")) {
        hdfDataType = DFNT_UINT32;
        sdsid = SDcreate(sdfid, "l3m_data", hdfDataType, 2, dims);
        uint32_t tmp = productStuff[0]->missingValue;
        DPTB(SDsetattr(sdsid, "Fill", hdfDataType, 1, (VOIDP)&tmp));
        fileData = allocateMemory(width * sizeof(uint32_t), "OutFile_hdf4::open fileData");
    } else if (!strcmp(productStuff[0]->productInfo->dataType, "float")) {
        hdfDataType = DFNT_FLOAT32;
        sdsid = SDcreate(sdfid, "l3m_data", hdfDataType, 2, dims);
        float tmp = productStuff[0]->missingValue;
        DPTB(SDsetattr(sdsid, "Fill", hdfDataType, 1, (VOIDP)&tmp));
        fileData = allocateMemory(width * sizeof(float), "OutFile_hdf4::open fileData");
    } else if (!strcmp(productStuff[0]->productInfo->dataType, "double")) {
        hdfDataType = DFNT_FLOAT64;
        sdsid = SDcreate(sdfid, "l3m_data", hdfDataType, 2, dims);
        double tmp = productStuff[0]->missingValue;
        DPTB(SDsetattr(sdsid, "Fill", hdfDataType, 1, (VOIDP)&tmp));
        fileData = allocateMemory(width * sizeof(double), "OutFile_hdf4::open fileData");
    } else {
        printf("-E- Data type %s, not supported\n", productStuff[0]->productInfo->dataType);
        exit(EXIT_FAILURE);
    }

    tmpStr = "linear";
    DPTB(SDsetattr(sdsid, "Scaling", DFNT_CHAR, strlen(tmpStr) + 1, (VOIDP)tmpStr));
    tmpStr = "(Slope*l3m_data) + Intercept = Parameter value";
    DPTB(SDsetattr(sdsid, "Scaling Equation", DFNT_CHAR, strlen(tmpStr) + 1, (VOIDP)tmpStr));
    DPTB(SDsetattr(sdsid, "Slope", DFNT_FLOAT32, 1, (VOIDP)&tmpScale));
    DPTB(SDsetattr(sdsid, "Intercept", DFNT_FLOAT32, 1, (VOIDP)&tmpOffset));

    // create the SST quality data set
    if (qualityData) {
        quality_sdsid = SDcreate(sdfid, "l3m_qual", DFNT_UINT8, 2, dims);
        int32_t validRange[2];
        validRange[0] = 0;
        validRange[1] = 2;
        DPTB(SDsetattr(quality_sdsid, "valid_range", DFNT_INT32, 2, (VOIDP)validRange));
    }

    return true;
}

void OutFile_hdf4::writeLine() {
    productStuff[0]->calcOutputLineVals(fileData);

    int32_t start[2];
    int32_t count[2];

    start[0] = currentLine;
    start[1] = 0;
    count[0] = 1;
    count[1] = width;

    if ((SDwritedata(sdsid, start, NULL, count, (VOIDP)fileData)) < 0) {
        printf("\n-E- OutFile_hdf4::writeLine(): SDwritedata unsuccessful\n");
        exit(EXIT_FAILURE);
    }

    if (qualityData) {
        if ((SDwritedata(quality_sdsid, start, NULL, count, (VOIDP)qualityData)) < 0) {
            printf("\n-E- OutFile_hdf4::writeLine(): SDwritedata unsuccessful\n");
            exit(EXIT_FAILURE);
        }
    }

    currentLine++;
}

bool OutFile_hdf4::close() {
    if (fileData) {
        free(fileData);
        fileData = NULL;
    }

    if (metaData) {
        int32_t tmpInt = metaData->data_bins;
        DPTB(SDsetattr(sdfid, "Data Bins", DFNT_INT32, 1, (VOIDP)&tmpInt));
    }
    float tmpFloat = getFileMinVal();
    DPTB(SDsetattr(sdfid, "Data Minimum", DFNT_FLOAT32, 1, (VOIDP)&tmpFloat));
    tmpFloat = getFileMaxVal();
    DPTB(SDsetattr(sdfid, "Data Maximum", DFNT_FLOAT32, 1, (VOIDP)&tmpFloat));

    if (qualityData) {
        SDendaccess(quality_sdsid);
        quality_sdsid = -1;
    }
    SDendaccess(sdsid);
    sdsid = -1;
    SDend(sdfid);
    sdfid = -1;

    /*-----------------------------------------------------------------------*/
    /*  write map_palette */

    // make sure a palette has been loaded
    if (red && green && blue) {
        uint8_t data[768];
        int j = 0;
        for (int i = 0; i < 256; i++) {
            data[j++] = red[i];
            data[j++] = green[i];
            data[j++] = blue[i];
        }

        if ((DFPaddpal(fileName.c_str(), (VOIDP)data)) < 0) {
            printf("-E- OutFile_hdf4::close - Error writing map_palette.\n");
            return false;
        }

        uint16_t pal_ref;
        if ((pal_ref = DFPlastref()) > 0) {
            if ((DFANputlabel(fileName.c_str(), DFTAG_IP8, pal_ref, (char*)"palette")) < 0) {
                printf("-E- OutFile_hdf4::close - Error writing palette label\n");
                return false;
            }
        }
    }
    return true;
}

int32_t OutFile_hdf4::addProduct(productInfo_t* productInfo) {
    return addProductNonDisplay(productInfo);
}

//------------------------------------------------------------------------------
// OutFile_netcdf4
//------------------------------------------------------------------------------

using namespace netCDF;

OutFile_netcdf4::OutFile_netcdf4() : OutFile() {
    ncFile = NULL;
    fileData = NULL;
    landPixelValue = badPixelValue;
}

OutFile_netcdf4::~OutFile_netcdf4() {
    if (fileData) {
        free(fileData);
        fileData = NULL;
    }
}

void OutFile_netcdf4::setSize(int32_t width, int32_t height) {
    OutFile::setSize(width, height);
    if (fileData) {
        free(fileData);
        fileData = NULL;
    }
}

//  modifications must be done here
bool OutFile_netcdf4::open() {
    const char* tmpStr;
    currentLine = 0;

    try {
        // open file
        ncFile = new NcFile(fileName.c_str(), NcFile::replace);

        // write global metadata
        string prodName;
        size_t pos = fileName.find_last_of('/');
        if (pos == string::npos)
            prodName = fileName;
        else
            prodName = fileName.substr(pos + 1);
        ncFile->putAtt("product_name", prodName.c_str());

        ncFile->putAtt("instrument", metaData->sensor);
        std::string source = "satellite observations from ";
        source += metaData->sensor;
        ncFile->putAtt("title", metaData->title);
        ncFile->putAtt("project", PROJECT);

        if (strcmp(metaData->mission, "") != 0) {
            ncFile->putAtt("platform", metaData->mission);
            source += "-";
            source += metaData->mission;
        }
        ncFile->putAtt("source", source);

        ncFile->putAtt("temporal_range", metaData->prod_type);
        ncFile->putAtt("processing_version", metaData->pversion);

        strcpy(metaData->ptime, unix2isodate(time(NULL), 'G'));
        ncFile->putAtt("date_created", metaData->ptime);
        ncFile->putAtt("history", metaData->proc_con);

        ncFile->putAtt("l2_flag_names", metaData->flag_names);

        ncFile->putAtt("time_coverage_start", unix2isodate(metaData->startTime, 'G'));
        ncFile->putAtt("time_coverage_end", unix2isodate(metaData->endTime, 'G'));
        ncFile->putAtt("start_orbit_number", ncInt, metaData->start_orb);
        ncFile->putAtt("end_orbit_number", ncInt, metaData->end_orb);
        ncFile->putAtt("map_projection", mapProjection.c_str());

        ncFile->putAtt("latitude_units", metaData->lat_units);
        ncFile->putAtt("longitude_units", metaData->lon_units);

        ncFile->putAtt("northernmost_latitude", ncFloat, metaData->north);
        ncFile->putAtt("southernmost_latitude", ncFloat, metaData->south);
        ncFile->putAtt("westernmost_longitude", ncFloat, metaData->west);
        ncFile->putAtt("easternmost_longitude", ncFloat, metaData->east);
        ncFile->putAtt("geospatial_lat_max", ncFloat, metaData->north);
        ncFile->putAtt("geospatial_lat_min", ncFloat, metaData->south);
        ncFile->putAtt("geospatial_lon_max", ncFloat, metaData->east);
        ncFile->putAtt("geospatial_lon_min", ncFloat, metaData->west);

        double latStep = (metaData->north - metaData->south) / (float)getHeight();
        double lonStep;
        if (metaData->east < metaData->west)
            lonStep = (360 + metaData->east - metaData->west) / (float)getWidth();
        else
            lonStep = (metaData->east - metaData->west) / (float)getWidth();
        ncFile->putAtt("latitude_step", ncFloat, latStep);
        ncFile->putAtt("longitude_step", ncFloat, lonStep);
        ncFile->putAtt("sw_point_latitude", ncFloat, metaData->south + latStep / 2.0);
        ncFile->putAtt("sw_point_longitude", ncFloat, metaData->west + lonStep / 2.0);

        string resolutionStr = resolution2string(resolution);
        ncFile->putAtt("spatialResolution", resolutionStr);
        ncFile->putAtt("geospatial_lon_resolution", resolutionStr);
        ncFile->putAtt("geospatial_lat_resolution", resolutionStr);
        ncFile->putAtt("geospatial_lat_units", metaData->lat_units);
        ncFile->putAtt("geospatial_lon_units", metaData->lon_units);

        ncFile->putAtt("number_of_lines", ncInt, height);
        ncFile->putAtt("number_of_columns", ncInt, width);
        ncFile->putAtt("measure", "Mean");  // FIXME

        ncFile->putAtt("suggested_image_scaling_minimum", ncFloat, productStuff[0]->productInfo->displayMin);
        ncFile->putAtt("suggested_image_scaling_maximum", ncFloat, productStuff[0]->productInfo->displayMax);
        if (!strcmp(productStuff[0]->productInfo->displayScale, "log"))
            tmpStr = "LOG";
        else if (!strcmp(productStuff[0]->productInfo->displayScale, "arctan"))
            tmpStr = "ATAN";
        else
            tmpStr = "LINEAR";
        ncFile->putAtt("suggested_image_scaling_type", tmpStr);
        ncFile->putAtt("suggested_image_scaling_applied", "No");

        ncFile->putAtt("_lastModified", metaData->ptime);
        ncFile->putAtt("Conventions", "CF-1.6 ACDD-1.3");
        ncFile->putAtt("institution", INSTITUTION);
        ncFile->putAtt("standard_name_vocabulary", STDNAME_VOCABULARY);
        ncFile->putAtt("naming_authority", NAMING_AUTHORITY);

        // create id
        string id = metaData->pversion;
        if (id == "Unspecified") {
            id = "L3/";
        } else {
            id += "/L3/";
        }
        id += metaData->product_name;
        ncFile->putAtt("id", id);

        ncFile->putAtt("license", LICENSE);
        ncFile->putAtt("creator_name", CREATOR_NAME);
        ncFile->putAtt("publisher_name", PUBLISHER_NAME);
        ncFile->putAtt("creator_email", CREATOR_EMAIL);
        ncFile->putAtt("publisher_email", PUBLISHER_EMAIL);
        ncFile->putAtt("creator_url", CREATOR_URL);
        ncFile->putAtt("publisher_url", PUBLISHER_URL);

        ncFile->putAtt("processing_level", "L3 Mapped");
        ncFile->putAtt("cdm_data_type", "grid");  // deprecated?
        if (proj4String.length() > 0)
            ncFile->putAtt("proj4_string", proj4String);

        // Some missions have DOIs
        if (clo_isSet(optionList, "doi")) {
            ncFile->putAtt("identifier_product_doi_authority", "http://dx.doi.org");
            ncFile->putAtt("identifier_product_doi", clo_getString(optionList, "doi"));
        }
        if (clo_isSet(optionList, "suite")) {
            const char* keywordStr = getGCMDKeywords(clo_getString(optionList, "suite"));
            if (keywordStr) {
                ncFile->putAtt("keywords", keywordStr);
                ncFile->putAtt("keywords_vocabulary", KEYWORDS_VOCABULARY);
            }
        }  // suite was set

        // create group: processing_control
        NcGroup grp1 = ncFile->addGroup("processing_control");
        grp1.putAtt("software_name", metaData->soft_name);
        grp1.putAtt("software_version", metaData->soft_ver);
        if ((tmpStr = strrchr(metaData->infiles, '/')) != NULL)
            tmpStr++;
        else
            tmpStr = metaData->infiles;
        grp1.putAtt("input_sources", tmpStr);
        grp1.putAtt("l2_flag_names", metaData->flag_names);

        // create sub-group: input_parameters
        NcGroup grp2 = grp1.addGroup("input_parameters");
        char buf[2048];
        char* end_str;
        char* token = strtok_r(metaData->input_parms, "|", &end_str);
        while (token != NULL) {
            char* end_token;
            strcpy(buf, token);
            char* name = strtok_r(token, "=", &end_token);
            for (uint32_t i = 0; i < strlen(name); i++) {
                if (name[i] == ' ') {
                    name[i] = 0;
                    break;
                }
            }
            char* val = strtok_r(NULL, "|", &end_token);
            if (val == NULL)
                val = "";
            strcpy(buf, val);
            grp2.putAtt(name, buf);
            token = strtok_r(NULL, "|", &end_str);
        }

        // Define dimensions and coordinate variables
        vector<NcDim> dimIds;
        string coordinates;

        if (mapProjection == "Equidistant Cylindrical"  // SMI
            || mapProjection == "PlateCarree" || (proj4String.find("+proj=eqc") != string::npos)) {
            dimIds.push_back(ncFile->addDim("lat", height));
            dimIds.push_back(ncFile->addDim("lon", width));

            if (fullLatLon != LatLonOff)
                fullLatLon = LatLon1D;
        } else {
            dimIds.push_back(ncFile->addDim("y", height));
            dimIds.push_back(ncFile->addDim("x", width));
            if (fullLatLon != LatLonOff) {
                fullLatLon = LatLon2D;
                coordinates = "lat lon";
            }
        }

        // Define variables
        size_t dataSize = 1;
        vector<NcDim> dimIds_2D;
        vector<NcDim> dimIds_3D;
        std::copy(dimIds.begin(), dimIds.end(), std::back_inserter(dimIds_2D));

        if (!get_wv3d_2d_name_to_3d_expansion().empty()) {
            std::copy(dimIds.begin(), dimIds.end(), std::back_inserter(dimIds_3D));
            dimIds_3D.push_back(ncFile->addDim("wavelength", get_len_wv3d()));
            // { std::cout << "Wavelength 3D is " << get_len_wv3d() << std::endl; }

            vector<NcDim> dimIdswv3d = {dimIds_3D[2]};
            NcVar var = ncFile->addVar("wavelength", ncInt, dimIdswv3d);
            var.putAtt("long_name", "wavelengths");
            var.putAtt("units", "nm");
            var.putAtt("_FillValue", ncInt, -32767);
            var.putAtt("valid_min", ncInt, 0);
            var.putAtt("valid_max", ncInt, 20000);
            auto data = (*get_wv3d_2d_name_to_3d_expansion().begin()).second;
            var.putVar(data.data());
        }
        auto dimIds_to_set = &dimIds_2D;
        size_t count_3d_counter = 0;
        for (size_t prodNum = 0; prodNum < productStuff.size(); prodNum++) {
            // determine maximum product size
            NcType varDataType = getDataType(productStuff[prodNum]->dataStorage);
            dataSize = std::max(dataSize, varDataType.getSize());

            // create variable
            productInfo_t* pInfo = productStuff[prodNum]->productInfo;
            std::string temp_wave_name = pInfo->paramDesignator;
            std::string temp_prefix = pInfo->prefix;
            std::string temp_descp = pInfo->description;

            const auto wv3d_3d_name_to_2d_name = get_wv3d_3d_name_to_2d_name();
            std::string full_prod_name = getProductNameFull(pInfo);
            const std::string suffix_prod = pInfo->suffix;
            if (!suffix_prod.empty()) {
                // std::cout << "Prod " << full_prod_name << " has suffix "
                //           << suffix_prod << std::endl;
                std::size_t pos = full_prod_name.find(suffix_prod);
                // std::cout << "Position is here " << pos << "\n";
                full_prod_name = full_prod_name.substr(0, pos);
                // std::cout << "Name without a suffix " << full_prod_name
                //           << std::endl;
            }
            if (wv3d_3d_name_to_2d_name.count(full_prod_name) > 0) {
                const auto prod_name_no_suffix = wv3d_3d_name_to_2d_name.at(full_prod_name);
                const auto prod_name = prod_name_no_suffix + suffix_prod;

                if (product_3d_already_set.count(prod_name) > 0) {
                    product_3d_already_set.at(prod_name)++;
                    slice_2d_in_wv3d[prodNum] = product_3d_already_set.at(prod_name);
                    index_2d_3d[prodNum] = index_2d_3d.at(prod2d_indexes_last_index.at(prod_name));
                    continue;
                } else {
                    product_3d_already_set[prod_name] = 0;
                    slice_2d_in_wv3d[prodNum] = 0;
                }
                prod2d_indexes_last_index[prod_name] = prodNum;
                // std::cout << "Created 3D nc product name " << prod_name
                //           << std::endl;
                // printf("%s = %d = %s\n", pInfo->prefix, pInfo->prod_ix,
                //        pInfo->suffix);

                // std::cout << "Prod name " << prod_name << " "
                //           << pInfo->productName << std::endl;
                dimIds_to_set = &dimIds_3D;
                if (pInfo->paramDesignator)
                    free(pInfo->paramDesignator);
                if (pInfo->prefix)
                    free(pInfo->prefix);
                const auto new_desc = remove_wv_from_long_name(pInfo->description);
                if (pInfo->description)
                    free(pInfo->description);
                pInfo->paramDesignator = strdup("none");
                pInfo->prefix = strdup(prod_name_no_suffix.c_str());
                pInfo->description = strdup(new_desc.c_str());
                // {
                //     std::cout << "Prefix is now " << pInfo->prefix << std::endl;
                //     std::cout << "New Name is " << getProductNameFull(pInfo)
                //               << std::endl;
                // }
            } else {
                dimIds_to_set = &dimIds_2D;
            }
            NcVar var = createProduct(pInfo, varDataType, *dimIds_to_set);
            // {
            //     std::cout << "New var has been created "
            //               << getProductNameFull(pInfo) << std::endl;
            // }
            index_2d_3d[prodNum] = count_3d_counter;
            count_3d_counter++;
            // reset pInfo to default
            {
                // {
                //     std::cout << "temp_wave_name = " << temp_wave_name
                //               << std::endl;
                //     std::cout << "temp_prefix = " << temp_prefix << std::endl;
                //     std::cout << "temp_descp = " << temp_descp << std::endl;
                // }
                if (wv3d_3d_name_to_2d_name.count(full_prod_name) > 0) {
                    free(pInfo->paramDesignator);
                    free(pInfo->prefix);
                    free(pInfo->description);
                    pInfo->paramDesignator = strdup(temp_wave_name.c_str());
                    pInfo->prefix = strdup(temp_prefix.c_str());
                    pInfo->description = strdup(temp_descp.c_str());
                }
                // std::cout << "prodNum = " << prodNum << "; Old Name is "
                //           << getProductNameFull(pInfo) << std::endl;
            }
            initCompression(var);
            prodVars.push_back(var);

            // add display info
            var.putAtt("display_scale", pInfo->displayScale);
            var.putAtt("display_min", ncFloat, pInfo->displayMin);
            var.putAtt("display_max", ncFloat, pInfo->displayMax);

            // specify coordinates as needed
            if (coordinates.length() > 0)
                var.putAtt("coordinates", coordinates);

        }  // for prodNum
        // {
        //     for (const auto& pair : slice_2d_in_wv3d) {
        //         std::cout << "prod number =  " << pair.first
        //                   << " prod index =  " << pair.second << std::endl;
        //     }
        //     for (const auto& var : prodVars) {
        //         std::cout << var.getName() << " " << var.getDimCount()
        //                   << std::endl;
        //         for (const auto& dim : var.getDims()) {
        //             std::cout << dim.getName() << " " << dim.getSize() << "; ";
        //         }
        //         std::cout << "\n";
        //     }
        //     for (const auto& prodindexes : index_2d_3d) {
        //         std::cout << "2d index " << prodindexes.first
        //                   << " 3d index = " << prodindexes.second << std::endl;
        //     }
        // }
        // allocate buffer big enough for one line of the largest data type
        if (fileData)
            free(fileData);
        fileData = allocateMemory(width * dataSize, "OutFile_netcdf4::open fileData");

        // add the quality variable
        if (qualityData) {
            productInfo_t* qInfo = allocateProductInfo();
            if (qualityName.empty()) {
                qualityName = (string) "qual_" + getProductNameFull(productStuff[0]->productInfo);
            }
            if (!findProductInfo(qualityName.c_str(), metaData->sensorID, qInfo)) {
                cerr << "-E- OutFile_netcdf4::open - cannot find product " << qualityName << " in product.xml"
                     << endl;
                exit(EXIT_FAILURE);
            }
            qInfo->fillValue = NC_FILL_UBYTE;
            qualVar = createProduct(qInfo, ncUbyte, dimIds);
            initCompression(qualVar);
            freeProductInfo(qInfo);
        }

        // add latitude and longitude variables

        productInfo_t* latInfo = allocateProductInfo();
        if (!findProductInfo("lat", metaData->sensorID, latInfo)) {
            cerr << "-E- OutFile_netcdf4::open - "
                    "cannot find \"lat\" in product.xml"
                 << endl;
            exit(EXIT_FAILURE);
        }
        productInfo_t* lonInfo = allocateProductInfo();
        if (!findProductInfo("lon", metaData->sensorID, lonInfo)) {
            cerr << "-E- OutFile_netcdf4::open - "
                    "cannot find \"lon\" in product.xml"
                 << endl;
            exit(EXIT_FAILURE);
        }

        if (fullLatLon == LatLon1D) {
            float* latarray = (float*)allocateMemory(getHeight() * sizeof(float), "latarray");
            for (int i = 0; i < getHeight(); i++)
                latarray[i] = metaData->north - latStep * i - latStep / 2.0;
            vector<NcDim> latDim;
            latDim.push_back(dimIds[0]);
            NcVar lat = createProduct(latInfo, ncFloat, latDim);
            lat.putVar(latarray);
            free(latarray);

            float* lonarray = (float*)allocateMemory(getWidth() * sizeof(float), "lonarray");
            for (int i = 0; i < getWidth(); i++) {
                lonarray[i] = metaData->west + lonStep * i + lonStep / 2.0;
                if(lonarray[i] > 180.0)
                    lonarray[i] -= 360.0;
            }
            vector<NcDim> lonDim;
            lonDim.push_back(dimIds[1]);
            NcVar lon = createProduct(lonInfo, ncFloat, lonDim);
            lon.putVar(lonarray);
            free(lonarray);
        } 
        else if (fullLatLon == LatLon2D) {
            latVar = createProduct(latInfo, ncFloat, dimIds);
            lonVar = createProduct(lonInfo, ncFloat, dimIds);
            initCompression(latVar);
            initCompression(lonVar);
        }

        freeProductInfo(latInfo);
        freeProductInfo(lonInfo);

        // store palette
        if (red && green && blue) {
            uint8_t data[768];
            int j = 0;
            for (int i = 0; i < 256; i++) {
                data[j++] = red[i];
                data[j++] = green[i];
                data[j++] = blue[i];
            }
            vector<NcDim> dimIds;
            dimIds.push_back(ncFile->addDim("rgb", 3));
            dimIds.push_back(ncFile->addDim("eightbitcolor", 256));
            NcVar var = ncFile->addVar("palette", ncUbyte, dimIds);
            var.putVar(data);
        }

    } catch (std::exception const& e) {
        cerr << "Exception: " << e.what() << endl;
    }

    return true;
}
//  modifications must be done here
void OutFile_netcdf4::writeLine() {
    vector<size_t> start(2);
    vector<size_t> count(2);
    start[0] = currentLine;
    start[1] = 0;
    count[0] = 1;
    count[1] = width;
    vector<size_t> start_3D = {start[0], start[1], 0};
    vector<size_t> count_3D = {count[0], count[1], 1};
    // calculate and write one line of data for each product
    for (size_t prodNum = 0; prodNum < productStuff.size(); prodNum++) {
        productStuff[prodNum]->calcOutputLineVals(fileData);
        const auto index_in_products = index_2d_3d.at(prodNum);
        // {
        //     std::cout <<"\nWriting " << currentLine << " " <<  prodNum << "
        //     :" << index_in_products << std::endl;
        // }
        if (slice_2d_in_wv3d.count(prodNum) == 0) {         
            prodVars[index_in_products].putVar(start, count, fileData);
        } else {
            const auto slice_to_put = slice_2d_in_wv3d.at(prodNum);
            start_3D[2] = slice_to_put;
            prodVars[index_in_products].putVar(start_3D, count_3D, fileData);
        }
    }

    // write quality data for this line
    if (qualityData) {
        qualVar.putVar(start, count, qualityData);
    }

    // write lat, lon for this line
    if (fullLatLon == LatLon2D) {
        if (latData)
            latVar.putVar(start, count, latData);
        if (lonData)
            lonVar.putVar(start, count, lonData);
    }

    currentLine++;
}

bool OutFile_netcdf4::close() {
    if (metaData) {
        ncFile->putAtt("data_bins", ncInt64, metaData->data_bins);
        if (metaData->data_bins <= 0) {
            ncFile->putAtt("data_minimum", ncFloat, 0.0);
            ncFile->putAtt("data_maximum", ncFloat, 0.0);
        } else {
            ncFile->putAtt("data_minimum", ncFloat, (float)getFileMinVal());
            ncFile->putAtt("data_maximum", ncFloat, (float)getFileMaxVal());
        }
    }
    ncFile->close();
    prodVars.clear();

    if (fileData) {
        free(fileData);
        fileData = NULL;
    }
    if (qualityData) {
        free(qualityData);
        qualityData = NULL;
    }

    return true;
}

int32_t OutFile_netcdf4::addProduct(productInfo_t* productInfo) {
    return addProductNonDisplay(productInfo);
}

void OutFile_netcdf4::initCompression(NcVar var) {
    // if deflate level below 0 then don't deflate
    if (getDeflate() < 1)
        return;

    vector<NcDim> dims = var.getDims();
    vector<size_t> chunksize;
    /*
     * vary chunk size based on dimensions
     * looking to keep the chunks around 32Kb, hopefully no larger than 200Kb
     */

    switch(dims.size()) {
        case 1:
            chunksize.push_back(MIN(dims[0].getSize(), 2048));
            break;
        case 2:
            chunksize.push_back(MIN(dims[0].getSize(), 512));
            chunksize.push_back(MIN(dims[1].getSize(), 1024));
            break;
        case 3:
            chunksize.push_back(MIN(dims[0].getSize(), 16));
            chunksize.push_back(MIN(dims[1].getSize(), 1024));
            chunksize.push_back(MIN(dims[2].getSize(), 8));
            break;
        default:
            printf("-E- OutFile_netcdf4::initCompression - bad rank for variable %s\n", var.getName().c_str());
            exit(EXIT_FAILURE);
    }

    var.setChunking(NcVar::nc_CHUNKED, chunksize);
    var.setCompression(true, true, getDeflate());
}

NcVar OutFile_netcdf4::createProduct(productInfo_t* pInfo, const NcType& ncType,
                                     const std::vector<netCDF::NcDim> ncDim) {
    // create variable
    NcVar var = ncFile->addVar(getProductNameFull(pInfo), ncType, ncDim);

    // add standard metadata
    if (pInfo->description)
        var.putAtt("long_name", pInfo->description);
    if (pInfo->scaleFactor != 1.0 || pInfo->addOffset != 0.0) {
        var.putAtt("scale_factor", ncFloat, pInfo->scaleFactor);
        var.putAtt("add_offset", ncFloat, pInfo->addOffset);
    }
    if (pInfo->units != NULL && strcmp(pInfo->units, "dimensionless") != 0)
        var.putAtt("units", pInfo->units);
    if (pInfo->standardName != NULL && strcmp(pInfo->standardName, "") != 0)
        var.putAtt("standard_name", pInfo->standardName);
    try {
        if (pInfo->fillValue)
            var.putAtt("_FillValue", ncType, pInfo->fillValue);
    } catch (std::exception const& e) {
        cerr << "FillValue exception: " << e.what() << endl;
    }
    if (pInfo->validMin != PRODUCT_DEFAULT_validMin || pInfo->validMax != PRODUCT_DEFAULT_validMax) {
        var.putAtt("valid_min", ncType, (pInfo->validMin - pInfo->addOffset) / pInfo->scaleFactor);
        var.putAtt("valid_max", ncType, (pInfo->validMax - pInfo->addOffset) / pInfo->scaleFactor);
    }
    if (pInfo->reference != NULL && strcmp(pInfo->reference, "") != 0)
        var.putAtt("reference", pInfo->reference);
    if (pInfo->comment != NULL && strcmp(pInfo->comment, "") != 0)
        var.putAtt("comment", pInfo->comment);

    return var;
}

NcType OutFile_netcdf4::getDataType(DataStorage dataStorage) {
    switch (dataStorage) {
        case ByteDS:
            return NC_BYTE;  // ncByte
        case UByteDS:
            return NC_UBYTE;  // ncUbyte
        case ShortDS:
            return NC_SHORT;  // ncShort
        case UShortDS:
            return NC_USHORT;  // ncUshort
        case IntDS:
            return NC_INT;  // ncInt
        case UIntDS:
            return NC_UINT;  // ncUint
        case FloatDS:
            return NC_FLOAT;  // ncFloat
        case DoubleDS:
            return NC_DOUBLE;  // ncDouble
        default:
            cerr << "-E- OutFile_netcdf4::getDataType - illegal data storage "
                    "type = "
                 << dataStorage << endl;
            exit(EXIT_FAILURE);
    }
}
