#ifndef MATH_UTILS_H
#define MATH_UTILS_H
#include "stdlib.h"
#include "stdio.h"
#include "math.h"

/**
 * @brief Matrix transpose
 * 
 * @param inp input
 * @param out output 
 */
void transpose3d(const float inp[3][3], float out[3][3]);
#endif