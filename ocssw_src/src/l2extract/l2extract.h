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

int extractNetCDF(const char* infile, const char* outfile, int spix, int epix,
                  int sscan, int escan, const char* prodlist, const char* wavelist);

int l2extract_init_options(clo_optionList_t* list, const char* softwareVersion);
int l2extract_read_options(clo_optionList_t* list, int argc, char* argv[]);

#endif /* L2EXTRACT_H_ */
