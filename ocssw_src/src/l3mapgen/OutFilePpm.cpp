#include "OutFilePpm.h"
#include <cstdio>
#include <cstdlib>
#include <genutils.h>
#include <cstddef>
#include <math.h>

void OutFilePpm::setSize(int32_t width, int32_t height) {
    OutFile::setSize(width, height);
    if (fileData)
        free(fileData);
    fileData = (uint8_t*)allocateMemory(width * 3, "OutFilePpm::fileData");
}

bool OutFilePpm::open() {
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

void OutFilePpm::writeLine() {
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

