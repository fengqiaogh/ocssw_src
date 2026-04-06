#include "OutFilePng.h"
#include <genutils.h>
#include <string>
#include <math.h>

using namespace std;

OutFilePng::OutFilePng(bool color) : OutFile() {
    isColor = color;
    outfp = nullptr;
    fileData = nullptr;
    infoPtr = nullptr;
    pngPtr = nullptr;
    numText = 10;
}

OutFilePng::~OutFilePng() {
    if (outfp)
        fclose(outfp);
    if (fileData)
        free(fileData);
}

void OutFilePng::setSize(int32_t width, int32_t height) {
    OutFile::setSize(width, height);
    if (fileData)
        free(fileData);
    fileData = (uint8_t*)allocateMemory(width, "OutFilePng::fileData");
}

bool OutFilePng::open() {
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

    if (isColor) {
        // color palette
        png_set_IHDR(pngPtr, infoPtr, width, height, 8, PNG_COLOR_TYPE_PALETTE, PNG_INTERLACE_NONE,
                     PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);

        uint8_t pal[256 * 3];
        for (int i = 0; i < 256; i++) {
            pal[i * 3] = red[i];
            pal[i * 3 + 1] = green[i];
            pal[i * 3 + 2] = blue[i];
        }
        png_set_PLTE(pngPtr, infoPtr, (png_const_colorp)pal, 256);

        if (transparent) {
            uint8_t transPal[256];
            for (int i = 0; i < 255; i++) {
                transPal[i] = 255;
            }
            transPal[255] = 0;
            png_set_tRNS(pngPtr, infoPtr, transPal, 256, nullptr);
        }
    } else {
        // GRAYSCALE
        png_set_IHDR(pngPtr, infoPtr, width, height, 8, PNG_COLOR_TYPE_GRAY, PNG_INTERLACE_NONE,
                     PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);
    }
    png_write_info(pngPtr, infoPtr);

    return true;
}

void OutFilePng::writeLine() {
    for (int i = 0; i < width; i++)
        fileData[i] = (uint8_t)round(productStuff[0]->calcOutputVal(productStuff[0]->lineData[i]));
    png_write_row(pngPtr, (png_bytep)fileData);
    currentLine++;
}

bool OutFilePng::close() {
    png_write_end(pngPtr, infoPtr);
    png_destroy_write_struct(&pngPtr, &infoPtr);

    fclose(outfp);
    outfp = nullptr;
    return true;
}
