#ifndef IAS_MATH_H
#define IAS_MATH_H

#include "ias_const.h"
#include "ias_structures.h"

void ias_math_compute_3dvec_cross
(
    const IAS_VECTOR *vec1,     /* I: Input vector number one       */
    const IAS_VECTOR *vec2,     /* I: Input vector number two       */
    IAS_VECTOR *vec3            /* O: Output vector (cross product) */
);

double ias_math_compute_3dvec_dot
(
    const IAS_VECTOR *vec1,     /* I: Vector one to be multiplied  */
    const IAS_VECTOR *vec2      /* I: Vector two to be multiplied  */
);

double ias_math_compute_vector_length
(
    const IAS_VECTOR *vec       /* I: Vector to find the length of */
);

int ias_math_compute_unit_vector
(
    const IAS_VECTOR *vec,      /* I: Input vector */
    IAS_VECTOR *unit_vector     /* O: Unit vector of the input vector */
);

double ias_math_eval_legendre
(
    double x,                   /* I: point at which to compute value */
    const double *coefficients, /* I: Array of Coefficients */
    int num_coefficients        /* I: Number of Coefficients in array */
);

double ias_math_interpolate_lagrange
(
    const double *p_YY, /* I: Pointer to the array of input values */
    const double *p_XX, /* I: Pointer to the array of times closest to the
                              requested time */
    int n_pts,          /* I: Number of points for the interpolation */
    double in_time      /* I: Requested time for interpolation */
);

void ias_math_interpolate_lagrange_3dvec
( 
    const IAS_VECTOR *p_YY, /* I: Y-axis or values to be interpolated */
    const double *p_XX,     /* I: X-axis of value to be interpolated */
    int num_points,    /* I: Number of points to use in the interpolation */
    double in_time,    /* I: X-value for which Y value is to be calculated */
    IAS_VECTOR *output /* O: Output result */
);

/* math constants */
double ias_math_get_pi();
double ias_math_get_arcsec_to_radian_conversion();
double ias_math_get_radians_per_degree();
double ias_math_get_degrees_per_radian();

#endif
