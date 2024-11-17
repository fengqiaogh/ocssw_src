//
// Created by A.Semenov on 11/25/23.
//
#include "libnav.h"

/**
 * @brief computes cross product v3 = v1 x v2
 * @param v1 first vector
 * @param v2 second vector
 * @param v3 result
 * @return
 */
int crossp_(float *v1, float *v2, float *v3) {
    v3[0] = v1[1] * v2[2] - v1[2] * v2[1];
    v3[1] = v1[2] * v2[0] - v1[0] * v2[2];
    v3[2] = v1[0] * v2[1] - v1[1] * v2[0];
    return 0;
}