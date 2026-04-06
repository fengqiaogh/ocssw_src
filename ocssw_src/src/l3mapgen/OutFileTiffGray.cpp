#include "OutFileTiffGray.h"
#include <genutils.h>
#include <stdexcept>

using namespace std;

OutFileTiffGray::~OutFileTiffGray() {
    if (fileData)
        free(fileData);
}

void OutFileTiffGray::setSize(int32_t width, int32_t height) {
    OutFile::setSize(width, height);
    if (fileData)
        free(fileData);
    fileData = (float*)allocateMemory(width * sizeof(float), "OutFileTiffGray::fileData");
}

void OutFileTiffGray::writeLine() {
    for (int i = 0; i < width; i++)
        fileData[i] = productStuff[0]->lineData[i];

    if (TIFFWriteScanline(tiff, (void*)fileData, currentLine) == 0) {
        cerr << "-E- Could not write TIFF image line " << currentLine << endl;
        exit(EXIT_FAILURE);
    }
    currentLine++;
}

void OutFileTiffGray::setTiffColor() {
    TIFFSetField(tiff, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
    TIFFSetField(tiff, TIFFTAG_BITSPERSAMPLE, 32);
    TIFFSetField(tiff, TIFFTAG_SAMPLESPERPIXEL, 1);
    TIFFSetField(tiff, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP);
}
