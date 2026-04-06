#ifndef COMPUTE_ALPHA_H
#define COMPUTE_ALPHA_H

#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>


/**
 * @brief computes polarization correction frame rotation angle
 * algorithm provided by: F. S. Patt, SAIC, February 2003.
 * @param lon 
 * @param lat 
 * @param senz 
 * @param sena 
 * @param mnorm - mirror normal
 * @param numPix
 * @param alpha - output alpha angle
 */
void compute_alpha(float lon[], float lat[], float senz[], float sena[],
    double mnorm[3], int numPix, float alpha[]);


/**
 * @brief Take attitude quaternion and convert them to equivalent direction cos matrix
 * @param attQuat 1 scan's att_quat (should be length of 4)
 * @param output result to save into 
 * Provided by Fred Patt. From qtom.pro
 */
void convertQuatToCosineMatrix(double attQuat[4], gsl_matrix* output);


/**
 * @brief calculate euler sequence from euler angles in XYZ order.
 * @param eulerAngles 
 * @param output 
 * Provided by Fred Patt. From euler.pro
 */
void computeEulerSequence(float eulerAngles[3], gsl_matrix* output);

#endif