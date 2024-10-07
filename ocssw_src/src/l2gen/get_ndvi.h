/**
 * @name get_ndvi.h
 * @brief header file for NDVI/EVI
 * @authors Jakob Lindo (SSAI)
 * @date Feb 2024
 */

#include "l12_proto.h"

/**
 * @short Main entry point for getting NDVI/EVI producs
 * @param l1rec A level 1 record
 * @param prodnum A catalogue index indicating the desired product
 * @param prod An array provided by the caller to hold the product
 */
void get_ndvi_evi(l1str *l1rec, int prodnum, float prod[]);