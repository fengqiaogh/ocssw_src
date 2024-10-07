#ifndef  _L1_AVIRIS_H
#define  _L1_AVIRIS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

extern "C" int close_prism(prism4ocia_t *data);
extern "C" prism4ocia_t* open_prism(char *filename, prism4ocia_t **data);
extern "C" int read_prism(prism4ocia_t *data, int32_t recnum);
extern "C" char* getinbasename_av(char *file);
extern "C" int endianess(void);
extern "C" float checkTagLine_f(char *linein, char* tag);
extern "C" char *checkTagLine(char *linein, char* tag);
//extern "C" int readBinScanLine_float(float *Lt, int32_t recnum, int32_t npix, double *gain, int nbands, int numBands, int interleave, int swap, FILE *ptr);

#endif
