#include "OutFilePpmRgb.h"
#include <math.h>
#include <cstdio>
#include <cstdlib>
#include <genutils.h>
#include <cstddef>

OutFilePpmRgb::OutFilePpmRgb() : OutFile() {
    outfp = nullptr;
    fileData = nullptr;
}

OutFilePpmRgb::~OutFilePpmRgb() {
    if (outfp)
        fclose(outfp);
    if (fileData)
        free(fileData);
}


void OutFilePpmRgb::setSize(int32_t width, int32_t height) {
    OutFile::setSize(width, height);
    if (fileData)
        free(fileData);
    fileData = (uint8_t*)allocateMemory(width * 3, "OutFilePpmRgb::fileData");
}

bool OutFilePpmRgb::open() {
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

bool OutFilePpmRgb::close() {
    fclose(outfp);
    outfp = nullptr;
    return true;
}

void OutFilePpmRgb::setPixel(int32_t x, double val, int32_t prod) {
    fprintf(stderr, "-E- OutFilePpmRgb::setPixel - only RGB is implemented with this file type.\n");
    exit(EXIT_FAILURE);
}

void OutFilePpmRgb::setPixelRGB(int32_t x, float red, float green, float blue) {
    if (x < 0 || x >= width) {
        fprintf(stderr, "-E- OutFilePpmRgb::setPixel - x=%d is not within range, width=%d.\n", x, width);
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

void OutFilePpmRgb::landPixel(int32_t x) {
    if (x < 0 || x >= width) {
        fprintf(stderr, "-E- OutFilePpmRgb::landPixel - x=%d is not within range, width=%d.\n", x, width);
        exit(EXIT_FAILURE);
    }
    uint8_t* ptr = fileData + x * 3;
    *ptr = LAND_PIX;
    ptr++;
    *ptr = LAND_PIX;
    ptr++;
    *ptr = LAND_PIX;
}

void OutFilePpmRgb::fillPixel(int32_t x) {
    if (x < 0 || x >= width) {
        fprintf(stderr, "-E- OutFilePpmRgb::fillPixel - x=%d is not within range, width=%d.\n", x, width);
        exit(EXIT_FAILURE);
    }
    uint8_t* ptr = fileData + x * 3;
    *ptr = FILL_PIX;
    ptr++;
    *ptr = FILL_PIX;
    ptr++;
    *ptr = FILL_PIX;
}

void OutFilePpmRgb::missingPixel(int32_t x) {
    if (x < 0 || x >= width) {
        fprintf(stderr, "-E- OutFilePpmRgb::missingPixel - x=%d is not within range, width=%d.\n", x,
                width);
        exit(EXIT_FAILURE);
    }
    uint8_t* ptr = fileData + x * 3;
    *ptr = FILL_PIX;
    ptr++;
    *ptr = FILL_PIX;
    ptr++;
    *ptr = FILL_PIX;
}

void OutFilePpmRgb::writeLine() {
    fwrite(fileData, 1, width * 3, outfp);
    currentLine++;
}
