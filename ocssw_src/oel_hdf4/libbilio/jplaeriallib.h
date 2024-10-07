/*
 * jplaeriallib.h
 *
 *  Created on: Jun 12, 2015
 *      Author: rhealy
 */

#ifndef JPLAERIALLIB_H_
#define JPLAERIALLIB_H_

#include <stdio.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define BIP 0
#define BIL 1
#define BSQ 2

#define DEG_TO_RAD  .0174532925199432958

static const int itemSize = 500;

char* getinbasename_av(char *file);
char* getinbasename(char *file);
void readWavInfo_jpl(FILE* fp, char* tag, char* val);
void readNextLine_jpl(FILE* fp, char* tag, char* val);
int readBinScanLine4ocia_int2(float *Lt, int32_t recnum, int32_t npix, double *gain, int nbands, int numBands, int interleave, int swap, FILE *ptr);
int readBinScanLine4Ocip_float(float *Lt, int32_t recnum, int32_t npix, double *gain, int nbands, int numBands, int interleave, int swap, FILE *ptr);
int readBinScanLine_sub_int2(float *Lt, int32_t recnum, int32_t npix, double *gain, int nbands, int numBands, int interleave, int swap, FILE *ptr);
int readBinScanLine_int2(float *Lt, int32_t recnum, int32_t npix, double *gain, int nbands, int interleave, int swap, FILE *ptr);
int readBinScanLine_float(float *Lt, int32_t recnum, int32_t npix, double *gain, int nbands, int numBands, int interleave, int swap, FILE *ptr);
double getValidAngle(double *ang, int32_t npix, int32_t skip);
char* checkTagLine(char *linein, char* tag);
float checkTagLine_f(char *linein, char* tag);
int checkTagLine_i(char *linein, char* tag);
char* checkTagLine_m(char *linein, char *line, char* tag);
char* checknspTagLine(char *linein, char* tag);
void getPosVec(float lat, float lon, float alt, double *pos);
void getPosVecR(float lat, float lon, float alt, double *pos);

#ifdef __cplusplus
}
#endif


#endif /* JPLAERIALLIB_H_ */
