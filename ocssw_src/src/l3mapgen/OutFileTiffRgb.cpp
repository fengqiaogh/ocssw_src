#include "OutFileTiffRgb.h"
#include <genutils.h>
#include <stdexcept>

using namespace std;

OutFileTiffRgb::~OutFileTiffRgb() {
    if (fileData)
        free(fileData);
}

void OutFileTiffRgb::setSize(int32_t width, int32_t height) {
    OutFile::setSize(width, height);
    if (fileData)
        free(fileData);
    if (transparent) {
        fileData = (uint8_t*)allocateMemory(width * 4, "OutFileTiffRgb::fileData");
    } else {
        fileData = (uint8_t*)allocateMemory(width * 3, "OutFileTiffRgb::fileData");
    }
}

void OutFileTiffRgb::writeLine() {
    if (TIFFWriteScanline(tiff, (void*)fileData, currentLine) == 0) {
        cerr << "-E- Could not write TIFF image line " << currentLine << endl;
        exit(EXIT_FAILURE);
    }
    currentLine++;
}

void OutFileTiffRgb::setTiffColor() {
    TIFFSetField(tiff, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
    TIFFSetField(tiff, TIFFTAG_BITSPERSAMPLE, 8);
    if (transparent) {
        TIFFSetField(tiff, TIFFTAG_SAMPLESPERPIXEL, 4);
    } else {
        TIFFSetField(tiff, TIFFTAG_SAMPLESPERPIXEL, 3);
    }
    TIFFSetField(tiff, TIFFTAG_COMPRESSION, COMPRESSION_ADOBE_DEFLATE);
}

void OutFileTiffRgb::setPixel(int32_t x, double val, int32_t prod) {
    if (x < 0 || x >= width) {
        fprintf(stderr, "-E- OutFileTiffRgb::setPixel - x=%d is not within range, width=%d.\n", x, width);
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

void OutFileTiffRgb::setPixelRGB(int32_t x, float red, float green, float blue) {
    if (x < 0 || x >= width) {
        fprintf(stderr, "-E- OutFileTiffRgb::setPixel - x=%d is not within range, width=%d.\n", x, width);
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

/**
 * @brief Set a land pixel in the TIFF RGB output file.
 * @param x Pixel x-coordinate.
 */
void OutFileTiffRgb::landPixel(int32_t x) {
    if (x < 0 || x >= width) {
        fprintf(stderr, "-E- OutFileTiffRgb::landPixel - x=%d is not within range, width=%d.\n", x, width);
        exit(EXIT_FAILURE);
    }

    uint8_t samplesPerPixel = 3;
    if (transparent) {
        samplesPerPixel = 4;
    }
    uint8_t* ptr = fileData + x * samplesPerPixel;
    *ptr = LAND_PIX;
    ptr++;
    *ptr = LAND_PIX;
    ptr++;
    *ptr = LAND_PIX;
    if (transparent) {
        ptr++;
        *ptr = 255;
    }
}

/**
 * @brief Set a fill pixel in the TIFF RGB output file.
 * @param x Pixel x-coordinate.
 */
void OutFileTiffRgb::fillPixel(int32_t x) {
    if (x < 0 || x >= width) {
        fprintf(stderr, "-E- OutFileTiffRgb::fillPixel - x=%d is not within range, width=%d.\n", x, width);
        exit(EXIT_FAILURE);
    }

    uint8_t samplesPerPixel = 3;
    if (transparent) {
        samplesPerPixel = 4;
    }
    uint8_t* ptr = fileData + x * samplesPerPixel;
    *ptr = FILL_PIX;
    ptr++;
    *ptr = FILL_PIX;
    ptr++;
    *ptr = FILL_PIX;
    if (transparent) {
        ptr++;
        *ptr = 0;
    }
}

/**
 * @brief Set a pixel as missing in the TIFF RGB output file.
 * @param x Pixel x-coordinate.
 */
void OutFileTiffRgb::missingPixel(int32_t x) {
    if (x < 0 || x >= width) {
        fprintf(stderr, "-E- OutFileTiffRgb::missingPixel - x=%d is not within range, width=%d.\n", x,
                width);
        exit(EXIT_FAILURE);
    }

    uint8_t samplesPerPixel = 3;
    if (transparent) {
        samplesPerPixel = 4;
    }
    uint8_t* ptr = fileData + x * samplesPerPixel;
    *ptr = FILL_PIX;
    ptr++;
    *ptr = FILL_PIX;
    ptr++;
    *ptr = FILL_PIX;
    if (transparent) {
        ptr++;
        *ptr = 0;
    }
}
