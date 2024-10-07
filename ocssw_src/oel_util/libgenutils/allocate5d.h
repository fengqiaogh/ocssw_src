/** @file allocate5d.h
    @brief Utility functions for allocating and freeing five-dimensional arrays of various types.

    This file was created by allocate5d.pl and should not be edited manually.
*/
#ifndef OEL_UTIL_LIBGENUTILS_ALLOCATE5D_H_
#define OEL_UTIL_LIBGENUTILS_ALLOCATE5D_H_

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/** @brief Allocate a five-dimensional array of type float of a given size.

    The order of the parameters are (now) the same as when the data is being accessed.
        E.g., this is valid (ignoring the printf type specifier):

    float *****array = allocate5d_float(10,1,7,5,1);
    printf("%d\n", array[9][0][6][4][0]);
    @param[in] nq slowest incrimenting dimension of array in memory.
    @param[in] nr really slow incrimenting dimension of array in memory.
    @param[in] nz slow incrimenting dimension of array in memory.
    @param[in] ny fast dimension of array.
    @param[in] nx fastest incrimenting dimension of array.

    @return A malloc'd array or NULL if any malloc fails.
*/
float *****allocate5d_float(size_t nq, size_t nr, size_t nz, size_t ny, size_t nx);

/** @brief Free a five-dimensional array created by allocate5d_float.

    @param[in] p Pointer to array created by allocate5d_float.
*/
void free5d_float(float *****p);

#ifdef __cplusplus
}
#endif

#endif /* OEL_UTIL_LIBGENUTILS_ALLOCATE5D_H_ */
