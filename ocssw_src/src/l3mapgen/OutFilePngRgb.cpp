#include "OutFilePngRgb.h"
#include <genutils.h>
#include <math.h>

using namespace std;

OutFilePngRgb::OutFilePngRgb() : OutFile() {
    outfp = nullptr;
    fileData = nullptr;
    infoPtr = nullptr;
    pngPtr = nullptr;
    numText = 10;
}

OutFilePngRgb::~OutFilePngRgb() {
    if (outfp)
        fclose(outfp);
    if (fileData)
        free(fileData);
}

void OutFilePngRgb::setSize(int32_t width, int32_t height) {
    OutFile::setSize(width, height);
    if (fileData)
        free(fileData);
    int samplesPerPixel = 3;
    if (transparent)
        samplesPerPixel = 4;
    fileData = (uint8_t*)allocateMemory(width * samplesPerPixel, "OutFilePngRgb::fileData");
}

bool OutFilePngRgb::open() {
    outfp = fopen(fileName.c_str(), "w");
    if (!outfp)
        return false;

    pngPtr = png_create_write_struct(PNG_LIBPNG_VER_STRING, nullptr, nullptr, nullptr);
    if (!pngPtr) {
        fprintf(stderr, "-E- Unable to create PNG write structure.\n");
        exit(EXIT_FAILURE);
    }

    infoPtr = png_create_info_struct(pngPtr);
    if (!infoPtr) {
        png_destroy_write_struct(&pngPtr, (png_infopp)nullptr);
        fprintf(stderr, "-E- Unable to create PNG info structure.\n");
        exit(EXIT_FAILURE);
    }
    if (setjmp(png_jmpbuf(pngPtr))) {
        png_destroy_write_struct(&pngPtr, (png_infopp)nullptr);
        fprintf(stderr, "-E- Unable to call PNG setjmp().\n");
        exit(EXIT_FAILURE);
    }
    png_init_io(pngPtr, outfp);

    png_text text_ptr[numText];
    text_ptr[0].key = "projString";
    text_ptr[0].text = const_cast<char*>(proj4String.c_str());
    text_ptr[0].compression = PNG_TEXT_COMPRESSION_NONE;
    text_ptr[1].key = "resolution";
    text_ptr[1].text = const_cast<char*>(to_string(resolution).c_str());
    text_ptr[1].compression = PNG_TEXT_COMPRESSION_NONE;
    text_ptr[2].key = "north";
    text_ptr[2].text = const_cast<char*>(to_string(metaData->north).c_str());
    text_ptr[2].compression = PNG_TEXT_COMPRESSION_NONE;
    text_ptr[3].key = "south";
    text_ptr[3].text = const_cast<char*>(to_string(metaData->south).c_str());
    text_ptr[3].compression = PNG_TEXT_COMPRESSION_NONE;
    text_ptr[4].key = "east";
    text_ptr[4].text = const_cast<char*>(to_string(metaData->east).c_str());
    text_ptr[4].compression = PNG_TEXT_COMPRESSION_NONE;
    text_ptr[5].key = "west";
    text_ptr[5].text = const_cast<char*>(to_string(metaData->west).c_str());
    text_ptr[5].compression = PNG_TEXT_COMPRESSION_NONE;
    text_ptr[6].key = "minX";
    text_ptr[6].text = const_cast<char*>(to_string(tiepoints[3]).c_str());
    text_ptr[6].compression = PNG_TEXT_COMPRESSION_NONE;
    text_ptr[7].key = "maxY";
    text_ptr[7].text = const_cast<char*>(to_string(tiepoints[4]).c_str());
    text_ptr[7].compression = PNG_TEXT_COMPRESSION_NONE;
    text_ptr[8].key = "height";
    text_ptr[8].text = const_cast<char*>(to_string(height).c_str());
    text_ptr[8].compression = PNG_TEXT_COMPRESSION_NONE;
    text_ptr[9].key = "width";
    text_ptr[9].text = const_cast<char*>(to_string(width).c_str());
    text_ptr[9].compression = PNG_TEXT_COMPRESSION_NONE;
    png_set_text(pngPtr, infoPtr, text_ptr, numText);

    if (transparent)
        png_set_IHDR(pngPtr, infoPtr, width, height, 8, PNG_COLOR_TYPE_RGBA, PNG_INTERLACE_NONE,
                     PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);
    else
        png_set_IHDR(pngPtr, infoPtr, width, height, 8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
                     PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);
    png_write_info(pngPtr, infoPtr);

    return true;
}

bool OutFilePngRgb::close() {
    if (setjmp(png_jmpbuf(pngPtr))) {
        fprintf(stderr, "-E- OutFilePngRgb::close - Unable to call PNG setjmp().\n");
        exit(EXIT_FAILURE);
    }
    png_write_end(pngPtr, infoPtr);
    png_destroy_write_struct(&pngPtr, &infoPtr);

    fclose(outfp);
    outfp = nullptr;
    return true;
}

void OutFilePngRgb::setPixel(int32_t x, double val, int32_t prod) {
    fprintf(stderr, "-E- OutFilePngRgb::setPixel - only RGB is implemented with this file type.\n");
    exit(EXIT_FAILURE);
}

void OutFilePngRgb::setPixelRGB(int32_t x, float red, float green, float blue) {
    if (x < 0 || x >= width) {
        fprintf(stderr, "-E- OutFilePngRgb::setPixel - x=%d is not within range, width=%d.\n", x, width);
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

void OutFilePngRgb::landPixel(int32_t x) {
    if (x < 0 || x >= width) {
        fprintf(stderr, "-E- OutFilePngRgb::landPixel - x=%d is not within range, width=%d.\n", x, width);
        exit(EXIT_FAILURE);
    }
    int samplesPerPixel = 3;
    if (transparent)
        samplesPerPixel = 4;
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

void OutFilePngRgb::fillPixel(int32_t x) {
    if (x < 0 || x >= width) {
        fprintf(stderr, "-E- OutFilePngRgb::fillPixel - x=%d is not within range, width=%d.\n", x, width);
        exit(EXIT_FAILURE);
    }
    int samplesPerPixel = 3;
    if (transparent)
        samplesPerPixel = 4;
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

void OutFilePngRgb::missingPixel(int32_t x) {
    if (x < 0 || x >= width) {
        fprintf(stderr, "-E- OutFilePngRgb::missingPixel - x=%d is not within range, width=%d.\n", x,
                width);
        exit(EXIT_FAILURE);
    }
    int samplesPerPixel = 3;
    if (transparent)
        samplesPerPixel = 4;
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

void OutFilePngRgb::writeLine() {
    if (setjmp(png_jmpbuf(pngPtr))) {
        fprintf(stderr, "-E- OutFilePngRgb::writeLine - Unable to call PNG setjmp().\n");
        exit(EXIT_FAILURE);
    }

    png_write_row(pngPtr, (png_bytep)fileData);
    currentLine++;
}
