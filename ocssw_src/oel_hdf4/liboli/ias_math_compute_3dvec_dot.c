/*******************************************************************************
NAME:                 ias_math_compute_3dvec_dot        

PURPOSE:        
Multiply two vectors of type IAS_VECTOR (Dot product)

RETURN VALUE:
Type = double
Value of computed dot product

*******************************************************************************/
#include "ias_math.h"

double ias_math_compute_3dvec_dot
(
    const IAS_VECTOR *vec1,            /* I: Vector one to be multiplied  */
    const IAS_VECTOR *vec2             /* I: Vector two to be multiplied  */
)
{
    double dot_product;

    dot_product = vec1->x * vec2->x + vec1->y * vec2->y + vec1->z * vec2->z;

    return dot_product;
}
