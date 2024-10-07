/*****************************************************************************
NAME:    ias_math_interpolate_lagrange_3dvec

PURPOSE:
    Perform lagrange interpolation on a vector.

RETURN VALUE: None

RETURN PARAMETERS:
    output:
        Type = IAS_VECTOR*
*****************************************************************************/
#include "ias_math.h"

void ias_math_interpolate_lagrange_3dvec
( 
    const IAS_VECTOR *p_YY, /* I: Y-axis or values to be interpolated */
    const double *p_XX,     /* I: X-axis of value to be interpolated */
    int num_points,    /* I: Number of points to use in the interpolation */
    double in_time,    /* I: X-value for which Y value is to be calculated */
    IAS_VECTOR *output /* O: Output result */
)
{
    double term_x, term_y, term_z;
    double sum_x, sum_y, sum_z;
    double s;
    int i, j;

    sum_x = sum_y = sum_z = 0.0;
    for (i = 0; i < num_points; i++)
    {
        term_x = p_YY[i].x;
        term_y = p_YY[i].y;
        term_z = p_YY[i].z;

        for (j = 0; j < num_points; j++)
        {
            if (j != i)
            {
                s = ( in_time - p_XX[j]) / ( p_XX[i] - p_XX[j] );

                term_x = term_x * s;
                term_y = term_y * s;
                term_z = term_z * s;
            }
        }

        sum_x = sum_x + term_x;
        sum_y = sum_y + term_y;
        sum_z = sum_z + term_z;
    }

    output->x = sum_x;
    output->y = sum_y;
    output->z = sum_z;
}

