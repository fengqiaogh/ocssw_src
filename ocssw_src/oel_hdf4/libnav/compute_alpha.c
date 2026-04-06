#include "compute_alpha.h"

/* ----------------------------------------------------------------------------- */
/* compute_alpha() - computes polarization correction frame rotation angle       */
/*                                                                               */
/* algorithm provided by: F. S. Patt, SAIC, February 2003.                       */
/* ----------------------------------------------------------------------------- */
void compute_alpha(float lon[],
                   float lat[],
                   float senz[],
                   float sena[],
                   double mnorm[3],
                   int npix,
                   float alpha[]) {
    static double radeg = 57.29577951;

    int ip, i;

    double slon, clon;
    double slat, clat;
    double szen, czen;
    double sazi, cazi;

    double e[3], n[3], v[3], r[3], s[3];
    double sdvcr, vdr, sdr, sdv;

    /* invert mirror normal */
    for (i = 0; i < 3; i++)
        r[i] = -mnorm[i];

    for (ip = 0; ip < npix; ip++) {
        slon = sin(lon[ip] / radeg);
        clon = cos(lon[ip] / radeg);
        slat = sin(lat[ip] / radeg);
        clat = cos(lat[ip] / radeg);
        szen = sin(senz[ip] / radeg);
        czen = cos(senz[ip] / radeg);
        sazi = sin(sena[ip] / radeg);
        cazi = cos(sena[ip] / radeg);

        /* pixel coordinate system (north, east, vertical) in ECR */
        e[0] = -slon;
        e[1] = clon;
        e[2] = 0.0;

        n[0] = -slat * clon;
        n[1] = -slat * slon;
        n[2] = clat;

        v[0] = clat * clon;
        v[1] = clat * slon;
        v[2] = slat;

        /* sensor view vector in ECR */
        for (i = 0; i < 3; i++)
            s[i] = e[i] * szen * sazi + n[i] * szen * cazi + v[i] * czen;

        /* compute rotation angle (alpha) from pixel normal (v) to mirror */
        /* normal (r) about sensor view vector (s)  (Wertz, p. 757)       */
        sdvcr = s[0] * (v[1] * r[2] - v[2] * r[1])
            + s[1] * (v[2] * r[0] - v[0] * r[2])
            + s[2] * (v[0] * r[1] - v[1] * r[0]);

        vdr = v[0] * r[0] + v[1] * r[1] + v[2] * r[2];
        sdr = s[0] * r[0] + s[1] * r[1] + s[2] * r[2];
        sdv = v[0] * s[0] + v[1] * s[1] + v[2] * s[2];

        /* negated to be consistent with Gordon et al. */
        alpha[ip] = -(radeg * atan2(sdvcr, (vdr - sdr * sdv)) - 90.0);
    }
}

/**
 * @brief Take attitude quaternion and convert them to equivalent direction cos matrix
 * @param attQuat 1 scan's att_quat (should be length of 4)
 * @param output result to save into 
 * Provided by Fred Patt. From qtom.pro
 */
void convertQuatToCosineMatrix(double attQuat[4], gsl_matrix* output) {

    double squareAtIndex0 = (attQuat[0]*attQuat[0]);
    double squareAtIndex1 = (attQuat[1]*attQuat[1]);
    double squareAtIndex2 = (attQuat[2]*attQuat[2]);
    double squareAtIndex3 = (attQuat[3]*attQuat[3]);

    gsl_matrix_set(output, 0, 0, squareAtIndex0 - squareAtIndex1 - squareAtIndex2 + squareAtIndex3);
    gsl_matrix_set(output, 1, 0, 2.0 * ((attQuat[0] * attQuat[1]) + (attQuat[2] * attQuat[3])));
    gsl_matrix_set(output, 2, 0, 2.0 * ((attQuat[0] * attQuat[2]) - (attQuat[1] * attQuat[3])));

    gsl_matrix_set(output, 0, 1, 2.0 * ((attQuat[0] * attQuat[1]) - (attQuat[2] * attQuat[3])));
    gsl_matrix_set(output, 1, 1, (-attQuat[0] * attQuat[0]) + squareAtIndex1 - squareAtIndex2 + squareAtIndex3);
    gsl_matrix_set(output, 2, 1, 2.0 * ((attQuat[1] * attQuat[2]) + (attQuat[0] * attQuat[3])));

    gsl_matrix_set(output, 0, 2, 2.0 * ((attQuat[0] * attQuat[2]) + (attQuat[1] * attQuat[3])));
    gsl_matrix_set(output, 1, 2, 2.0 * ((attQuat[1] * attQuat[2]) - (attQuat[0] * attQuat[3])));
    gsl_matrix_set(output, 2, 2, (-attQuat[0] * attQuat[0]) - squareAtIndex1 + squareAtIndex2 + squareAtIndex3);

    return;
}


/**
 * @brief calculate euler sequence from euler angles in XYZ order.
 * @param eulerAngles 
 * @param output 
 * Provided by Fred Patt. From euler.pro
 */
void computeEulerSequence(float eulerAngles[3], gsl_matrix* output) {
    static double radeg = 57.29577951;

    float cosX = cos(eulerAngles[0]/radeg);
    float sinX = sin(eulerAngles[0]/radeg);
    float cosY = cos(eulerAngles[1]/radeg);
    float sinY = sin(eulerAngles[1]/radeg);
    float cosZ = cos(eulerAngles[2]/radeg);
    float sinZ = sin(eulerAngles[2]/radeg);

    // note that these are set in row major order, so compared to IDL version
    // the row and columns are swapped

    gsl_matrix *rotationX = gsl_matrix_calloc(3, 3);
    gsl_matrix *rotationY = gsl_matrix_calloc(3, 3);
    gsl_matrix *rotationZ = gsl_matrix_calloc(3, 3);

    gsl_matrix_set(rotationX, 0, 0, 1.0); 
    gsl_matrix_set(rotationX, 1, 1, cosX);
    gsl_matrix_set(rotationX, 2, 2, cosX);
    gsl_matrix_set(rotationX, 2, 1, sinX);  // ie: idl == rotationX(1, 2) = sinX
    gsl_matrix_set(rotationX, 1, 2, -sinX); // ie: idl == rotationX(2, 1) = -sinX

    gsl_matrix_set(rotationY, 1, 1, 1.0);
    gsl_matrix_set(rotationY, 0, 0, cosY);
    gsl_matrix_set(rotationY, 2, 2, cosY);
    gsl_matrix_set(rotationY, 0, 2, sinY);
    gsl_matrix_set(rotationY, 2, 0, -sinY);

    gsl_matrix_set(rotationZ, 2, 2, 1.0);
    gsl_matrix_set(rotationZ, 1, 1, cosZ);
    gsl_matrix_set(rotationZ, 0, 0, cosZ);
    gsl_matrix_set(rotationZ, 1, 0, sinZ);
    gsl_matrix_set(rotationZ, 0, 1, -sinZ);


    // Compute total rotation as rotationZ * rotationY * rotationX
    // this variable holds the intermediate value of the first matrix multiplication
    gsl_matrix *rotationYX = gsl_matrix_calloc(3, 3);

    // YX = (Y * X) in IDL, so in row major order, do (X * Y)
    // alpha = 1 scales matrix A. (alpha)(A*B) = C 
    // beta > 1 weights the product. (beta)C + (beta)(A * B) = C. Keep at 0.
    // CblasNoTrans == dont flip the matrix
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, rotationX, rotationY, 0.0, rotationYX);
    
    // Z * (Y * X) in IDL, so in row major order, do: (Y * X) * Z
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, rotationYX, rotationZ, 0.0, output);

    gsl_matrix_free(rotationX);
    gsl_matrix_free(rotationY);
    gsl_matrix_free(rotationZ);
    gsl_matrix_free(rotationYX);

    return;
}
