/*
 * l2extract.h
 *
 *  Created on: May 6, 2014
 *      Author: dshea
 */

#ifndef L2EXTRACT_H_
#define L2EXTRACT_H_

#ifndef BOUNDS_ERROR
#define BOUNDS_ERROR 110
#endif

#include <clo.h>
#include <genutils.h>
#include "version.h"

extern clo_optionList_t* optionList;

int extractNetCDF(const char* inFile, const char* outFile, int startPixel, int endPixel,
                  int startScan, int endScan, const char* productList, const char* waveList);

#ifdef __cplusplus
extern "C" {
#endif

int l2extractInitOptions(clo_optionList_t* list, const char* softwareVersion);
int l2extractReadOptions(clo_optionList_t* list, int argc, char* argv[]);

#ifdef __cplusplus
}
#endif

#endif /* L2EXTRACT_H_ */
