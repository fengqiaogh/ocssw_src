#include "OutFile.h"
#include "l3mapgen.h"

#include <boost/algorithm/string.hpp>
#include <fstream>
#include <sstream>
#include <regex>

using namespace std;

//------------------------------------------------------------------------------
// OutFile::ProductStuff
//------------------------------------------------------------------------------

OutFile::ProductStuff::ProductStuff(int32_t width, const productInfo_t* productInfo, double landPixelValue) {
    this->width = width;
    this->productInfo = allocateProductInfo();
    copyProductInfo(this->productInfo, productInfo);
    dataStorage = UBYTE_DS;
    scaleType = LINEAR;
    scale = 1.0;
    offset = 0.0;
    minOutputVal = 0;
    maxOutputVal = 254;
    minVal = minOutputVal;
    maxVal = maxOutputVal;
    missingValue = FILL_PIX;
    lineData = (double*)allocateMemory(width * sizeof(double), "OutFile::ProductStuff::lineData");
    this->landPixelValue = landPixelValue;
}

OutFile::ProductStuff::ProductStuff(const OutFile::ProductStuff& productStuff) {
    width = productStuff.width;
    productInfo = allocateProductInfo();
    copyProductInfo(productInfo, productStuff.productInfo);
    dataStorage = productStuff.dataStorage;
    scaleType = productStuff.scaleType;
    scale = productStuff.scale;
    offset = productStuff.offset;
    minOutputVal = productStuff.minOutputVal;
    maxOutputVal = productStuff.maxOutputVal;
    minVal = productStuff.minVal;
    maxVal = productStuff.maxVal;
    missingValue = productStuff.missingValue;
    lineData = (double*)allocateMemory(width * sizeof(double), "OutFile::ProductStuff::lineData");
    landPixelValue = productStuff.landPixelValue;
}

OutFile::ProductStuff::~ProductStuff() {
    freeProductInfo(productInfo);
    free(lineData);
}

void OutFile::ProductStuff::setScale(double min, double max, ScaleType scaleType) {
    this->scaleType = scaleType;
    minVal = min;
    maxVal = max;
    switch (scaleType) {
        case LINEAR:
            offset = min;
            scale = (max - offset) / (maxOutputVal - minOutputVal);
            break;
        case LOG:
            offset = log10(min);
            scale = (log10(max) - offset) / (maxOutputVal - minOutputVal);
            break;
        case ARCTAN:
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
        return LAND_PIX;

    // don't scale if output type is floating point
    if (dataStorage == FLOAT_DS || dataStorage == DOUBLE_DS)
        return val;

    switch (scaleType) {
        case LINEAR:
            outVal = (val - offset) / scale;
            break;
        case LOG:
            if (val < 0)
                return minOutputVal;
            else
                outVal = (log10(val) - offset) / scale;
            break;
        case ARCTAN:
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
    if (dataStorage == FLOAT_DS || dataStorage == DOUBLE_DS)
        return val;

    switch (scaleType) {
        case LINEAR:
            physicalVal = val * scale + offset;
            break;
        case LOG:
            physicalVal = pow(10, val * scale + offset);
            break;
        case ARCTAN:
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
        case BYTE_DS:
            for (int i = 0; i < width; i++)
                ((int8_t*)lineBuffer)[i] = round(calcOutputVal(lineData[i]));
            break;
        case UBYTE_DS:
            for (int i = 0; i < width; i++)
                ((uint8_t*)lineBuffer)[i] = round(calcOutputVal(lineData[i]));
            break;
        case SHORT_DS:
            for (int i = 0; i < width; i++)
                ((int16_t*)lineBuffer)[i] = round(calcOutputVal(lineData[i]));
            break;
        case USHORT_DS:
            for (int i = 0; i < width; i++)
                ((uint16_t*)lineBuffer)[i] = round(calcOutputVal(lineData[i]));
            break;
        case INT_DS:
            for (int i = 0; i < width; i++)
                ((int32_t*)lineBuffer)[i] = round(calcOutputVal(lineData[i]));
            break;
        case UINT_DS:
            for (int i = 0; i < width; i++)
                ((uint32_t*)lineBuffer)[i] = round(calcOutputVal(lineData[i]));
            break;
        case FLOAT_DS:
            for (int i = 0; i < width; i++)
                ((float*)lineBuffer)[i] = lineData[i];
            break;
        case DOUBLE_DS:
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
    width = 0;
    height = 0;
    qualityData = NULL;
    currentLine = 0;
    colorType = GRAYSCALE;
    fileMinVal = DBL_MAX;
    fileMaxVal = 0 - DBL_MAX;
    resolution = 0;
    deflate = 0;
    fullLatLon = LAT_LON_2D;
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
    rgbLand = (uint8_t*)allocateMemory(3, "rgbLand");
    rgbLand[0] = 160;
    rgbLand[1] = 82;
    rgbLand[2] = 45;
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
    free(rgbLand);
    for (size_t i = 0; i < productStuff.size(); i++) {
        delete productStuff[i];
    }
    productStuff.clear();
    free(metaData);
    if (qualityData)
        free(qualityData);
}

string OutFile::getScaleTypeString(int32_t prod) {
    switch (productStuff[prod]->scaleType) {
        case LOG:
            return "LOG";
        case LINEAR:
            return "LINEAR";
        case ARCTAN:
            return "ARCTAN";
    }
    return "LINEAR";
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
    if (fullLatLon == LAT_LON_2D) {
        // NcVar::putVar takes care of converting double to float
        latData = lat;
        lonData = lon;
    }
}

bool OutFile::setPalette(const char* paletteName, bool applyMask) {
    string dataRoot = getenv("OCDATAROOT") ? getenv("OCDATAROOT") : "";
    string paletteFileName;
    regex pattern("/");
    bool isBinary = false;
    short r[256], g[256], b[256];
    struct Color {
        unsigned char r, g, b;
    };
    vector<Color> palette(256);

    if (dataRoot.empty()) {
        printf("OCDATAROOT environment variable is not defined.\n");
        return (EXIT_FAILURE);
    }

    if(!regex_search(paletteName, pattern)) {
        paletteFileName = dataRoot;
        paletteFileName += "/common/palette/";
        paletteFileName += paletteName;

        ifstream file(paletteFileName, ios::binary);

        if(file.is_open()) {
            vector<char> buffer(512); // Read the first 512 bytes
            file.read(buffer.data(), buffer.size());
            std::streamsize bytesRead = file.gcount();

            for (int i = 0; i < bytesRead; ++i) {
                unsigned char c = static_cast<unsigned char>(buffer[i]);
                // Check for characters that are not common printable ASCII or whitespace
                if (!std::isprint(c) && !std::isspace(c) && c != '\t' && c != '\r' && c != '\n') {
                    isBinary = true;
                    break;
                }
            }
        } else {
            paletteFileName += ".pal";
        }
    } else {
        paletteFileName = paletteName;

        ifstream file(paletteFileName, ios::binary);

        if(file.is_open()) {
            vector<char> buffer(512); // Read the first 512 bytes
            file.read(buffer.data(), buffer.size());
            std::streamsize bytesRead = file.gcount();

            for (int i = 0; i < bytesRead; ++i) {
                unsigned char c = static_cast<unsigned char>(buffer[i]);
                // Check for characters that are not common printable ASCII or whitespace
                if (!std::isprint(c) && !std::isspace(c) && c != '\t' && c != '\r' && c != '\n') {
                    isBinary = true;
                    break;
                }
            }
        } else {
            cerr << "palette file " <<  paletteFileName << " not found" << endl;
            return false;
        }
    }

    ifstream file(paletteFileName, ios::binary);
    if(!file.is_open()) {
        cerr << "Error reading palette file " <<  paletteFileName << endl;
        return false;
    }

    if(isBinary) {
        //read binary
        file.read(reinterpret_cast<char*>(palette.data()), 256 * sizeof(Color));
    } else {
        //read txt
        string line;
        int i = 0;
        int rr, gg, bb;
        while(std::getline(file, line)){
            std::istringstream iss(line);
            iss >> rr >> gg >> bb;
            palette[i].r = rr;
            palette[i].g = gg;
            palette[i].b = bb;
            i++;
        }
    }
    file.close();

    for(int i = 0; i < 256; ++i){
        r[i] = (int)palette[i].r;
        g[i] = (int)palette[i].g;
        b[i] = (int)palette[i].b;
    }

    if (applyMask) {
        r[LAND_PIX] = rgbLand[0];
        g[LAND_PIX] = rgbLand[1];
        b[LAND_PIX] = rgbLand[2];
        r[FILL_PIX] = 0;
        g[FILL_PIX] = 0;
        b[FILL_PIX] = 0;
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

    rgbLand[0] = stoi(rgb[0]);
    rgbLand[1] = stoi(rgb[1]);
    rgbLand[2] = stoi(rgb[2]);
}

void OutFile::setMetaData(meta_l3bType* metaData) {
    *this->metaData = *metaData;
}

int32_t OutFile::addProduct(productInfo_t* productInfo, bool applyMask) {
    ProductStuff* stuff = new ProductStuff(width, productInfo, landPixelValue);
    if(applyMask) {
        stuff->maxOutputVal--;
    }

    // setup display scaling
    if (!strcmp(productInfo->displayScale, "linear"))
        stuff->setScale(productInfo->displayMin, productInfo->displayMax, LINEAR);
    else if (!strcmp(productInfo->displayScale, "log"))
        stuff->setScale(productInfo->displayMin, productInfo->displayMax, LOG);
    else if (!strcmp(productInfo->displayScale, "arctan"))
        stuff->setScale(productInfo->displayMin, productInfo->displayMax, ARCTAN);
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
        stuff->dataStorage = BYTE_DS;
        stuff->minOutputVal = SCHAR_MIN;
        stuff->maxOutputVal = SCHAR_MAX;
    } else if (!strcmp(productInfo->dataType, "ubyte")) {
        stuff->dataStorage = UBYTE_DS;
        stuff->minOutputVal = 0;
        stuff->maxOutputVal = UCHAR_MAX;
    } else if (!strcmp(productInfo->dataType, "short")) {
        stuff->dataStorage = SHORT_DS;
        stuff->minOutputVal = SHRT_MIN;
        stuff->maxOutputVal = SHRT_MAX;
    } else if (!strcmp(productInfo->dataType, "ushort")) {
        stuff->dataStorage = USHORT_DS;
        stuff->minOutputVal = 0;
        stuff->maxOutputVal = USHRT_MAX;
    } else if (!strcmp(productInfo->dataType, "int")) {
        stuff->dataStorage = INT_DS;
        stuff->minOutputVal = INT_MIN;
        stuff->maxOutputVal = INT_MAX;
    } else if (!strcmp(productInfo->dataType, "uint")) {
        stuff->dataStorage = UINT_DS;
        stuff->minOutputVal = 0;
        stuff->maxOutputVal = UINT_MAX;
    } else if (!strcmp(productInfo->dataType, "float")) {
        stuff->dataStorage = FLOAT_DS;
        stuff->minOutputVal = 0;
        stuff->maxOutputVal = 0;
    } else if (!strcmp(productInfo->dataType, "double")) {
        stuff->dataStorage = DOUBLE_DS;
        stuff->minOutputVal = 0;
        stuff->maxOutputVal = 0;
    } else {
        printf("-E- OutFile::addProductNonDisplay - invalid data type = %s\n", productInfo->dataType);
        exit(EXIT_FAILURE);
    }

    // setup scaling
    stuff->setScaleOffset(productInfo->scaleFactor, productInfo->addOffset, LINEAR);
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
