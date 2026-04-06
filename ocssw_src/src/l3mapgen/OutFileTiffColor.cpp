#include "OutFileTiffColor.h"
#include <genutils.h>
#include <cstdlib>

OutFileTiffColor::~OutFileTiffColor() {
    if (fileData)
        free(fileData);
}

void OutFileTiffColor::setSize(int32_t width, int32_t height) {
    OutFile::setSize(width, height);
    if (fileData)
        free(fileData);

    // if transparent we need to use RGBA mode instead of indexed color map
    if (transparent) {
        fileData = (uint8_t*)allocateMemory(width * 4, "OutFileTiffColor::fileData");
    } else {
        fileData = (uint8_t*)allocateMemory(width, "OutFileTiffColor::fileData");
    }
}

void OutFileTiffColor::writeLine() {
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
        std::cerr << "-E- Could not write TIFF image line " << currentLine << std::endl;
        exit(EXIT_FAILURE);
    }

    currentLine++;
}

void OutFileTiffColor::setTiffColor() {
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
            std::cerr << "-E- Could not allocate memory for TIFF color map" << std::endl;
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
