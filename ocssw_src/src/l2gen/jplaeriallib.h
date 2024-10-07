/*
 * jplaeriallib.h
 *
 *  Created on: Jun 12, 2015
 *      Author: rhealy
 */

#ifndef SRC_L2GEN_JPLAERIALLIB_H_
#define SRC_L2GEN_JPLAERIALLIB_H_
#include "l1_struc.h"


#define BIP 0
#define BIL 1
#define BSQ 2

#define DEG_TO_RAD  .0174532925199432958

#define BAD_FLT -32767.0
static const int itemSize = 500;

double getValidAngle(double *ang, int32_t npix, int32_t skip);
int readBinScanLine_int2(float *Lt, int32_t recnum, int32_t npix, double *gain, int nbands, int interleave, int swap, FILE *ptr);
int readBinScanLine_float(float *Lt, int32_t recnum, int32_t npix, double *gain, int nbands, int numBands, int interleave, int swap, FILE *ptr);
int readScanLine_jpl(l1str *l1rec, void* data, int32_t recnum);
char* getinbasename(char *file);
char* getinbasename_av(char *file);
char* checkTagLine(char *line, char* tag);
char* checknspTagLine(char *line, char* tag);
void readNextLine_jpl(FILE* fp, char* tag, char* val);
void readWavInfo_jpl(FILE* fp, char* tag, char* val);
int swapc_bytes(char *in, int nbyte, int ntime);
void trimBlanks(char* str);
void getPosVec(float lat, float lon, float alt, double *pos);

#endif /* SRC_L2GEN_JPLAERIALLIB_H_ */
