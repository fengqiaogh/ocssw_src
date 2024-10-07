/** @file allocate4d.h
    @brief Utility functions for allocating and freeing four-dimensional arrays of various types.

    This file was created by allocate4d.pl and should not be edited manually.
*/
#ifndef OEL_UTIL_LIBGENUTILS_ALLOCATE4D_H_
#define OEL_UTIL_LIBGENUTILS_ALLOCATE4D_H_

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/** @brief Allocate a four-dimensional array of type int of a given size.

    The order of the parameters are (now) the same as when the data is being accessed.
        E.g., this is valid (ignoring the printf type specifier):

    int ****array = allocate4d_int(1,7, 5, 1);
    printf("%d\n", array[0][6][4][0]);
    @param[in] nr slowest incrimenting dimension of array in memory.
    @param[in] nz slow incrimenting dimension of array in memory.
    @param[in] ny fast dimension of array.
    @param[in] nx fastest incrimenting dimension of array.

    @return A malloc'd array or NULL if any malloc fails.
*/
int ****allocate4d_int(size_t nr,size_t nz, size_t ny, size_t nx);

/** @brief Free a four-dimensional array created by allocate4d_int.

    @param[in] p Pointer to array created by allocate4d_int.
*/
void free4d_int(int ****p);

/** @brief Allocate a four-dimensional array of type float of a given size.

    The order of the parameters are (now) the same as when the data is being accessed.
        E.g., this is valid (ignoring the printf type specifier):

    float ****array = allocate4d_float(1,7, 5, 1);
    printf("%d\n", array[0][6][4][0]);
    @param[in] nr slowest incrimenting dimension of array in memory.
    @param[in] nz slow incrimenting dimension of array in memory.
    @param[in] ny fast dimension of array.
    @param[in] nx fastest incrimenting dimension of array.

    @return A malloc'd array or NULL if any malloc fails.
*/
float ****allocate4d_float(size_t nr,size_t nz, size_t ny, size_t nx);

/** @brief Free a four-dimensional array created by allocate4d_float.

    @param[in] p Pointer to array created by allocate4d_float.
*/
void free4d_float(float ****p);

/** @brief Allocate a four-dimensional array of type double of a given size.

    The order of the parameters are (now) the same as when the data is being accessed.
        E.g., this is valid (ignoring the printf type specifier):

    double ****array = allocate4d_double(1,7, 5, 1);
    printf("%d\n", array[0][6][4][0]);
    @param[in] nr slowest incrimenting dimension of array in memory.
    @param[in] nz slow incrimenting dimension of array in memory.
    @param[in] ny fast dimension of array.
    @param[in] nx fastest incrimenting dimension of array.

    @return A malloc'd array or NULL if any malloc fails.
*/
double ****allocate4d_double(size_t nr,size_t nz, size_t ny, size_t nx);

/** @brief Free a four-dimensional array created by allocate4d_double.

    @param[in] p Pointer to array created by allocate4d_double.
*/
void free4d_double(double ****p);

#ifdef __cplusplus
}
#endif

#endif /* OEL_UTIL_LIBGENUTILS_ALLOCATE4D_H_ */
