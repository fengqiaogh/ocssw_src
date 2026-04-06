#include "OutFilePgm.h"
#include <cstdio>
#include <cstdlib>
#include <genutils.h>
#include <cstddef>

OutFilePgm::OutFilePgm() : OutFile() {
    outfp = nullptr;
    fileData = nullptr;
}

OutFilePgm::~OutFilePgm() {
    if (outfp)
        fclose(outfp);
    if (fileData)
        free(fileData);
}

void OutFilePgm::setSize(int32_t width, int32_t height) {
    OutFile::setSize(width, height);
    if (fileData)
        free(fileData);
    fileData = (uint8_t*)allocateMemory(width, "OutFilePgm::fileData");
}

bool OutFilePgm::open() {
    outfp = fopen(fileName.c_str(), "w");
    if (!outfp)
        return false;

    /* Write pgm header */
    fprintf(outfp, "P5\n");
    fprintf(outfp, "%d %d\n", width, height);
    fprintf(outfp, "255\n");

    return true;
}

void OutFilePgm::writeLine() {
    for (int i = 0; i < width; i++)
        fileData[i] = (uint8_t)productStuff[0]->calcOutputVal(productStuff[0]->lineData[i]);
    fwrite(fileData, 1, width, outfp);
    currentLine++;
}

bool OutFilePgm::close() {
    fclose(outfp);
    outfp = nullptr;
    return true;
}
