/** @file allocate3d.h
	@brief Utility functions for allocating and freeing three-dimensional arrays of various types.

	This file was created by allocate3d.pl and should not be edited manually.
*/
#ifndef OEL_UTIL_LIBGENUTILS_ALLOCATE3D_H_
#define OEL_UTIL_LIBGENUTILS_ALLOCATE3D_H_

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/** @brief Allocate a three-dimensional array of type short of a given size.

	The order of the parameters are (now) the same as when the data is being accessed.
        E.g., this is valid (ignoring the printf type specifier):

	short ***array = allocate3d_short(7, 5, 1);
	printf("%d\n", array[6][4][0]);
	
	@param[in] nz slowest incrimenting dimension of array in memory.
	@param[in] ny middle dimension of array.
	@param[in] nx fastest incrimenting dimension of array.

	@return A malloc'd array or NULL if any malloc fails.
*/
short ***allocate3d_short(size_t nz, size_t ny, size_t nx);
/** @brief Free a three-dimensional array created by allocate3d_short.

	@param[in] p Pointer to array created by allocate3d_short.
*/
void free3d_short(short ***p);

/** @brief Allocate a three-dimensional array of type int of a given size.

	The order of the parameters are (now) the same as when the data is being accessed.
        E.g., this is valid (ignoring the printf type specifier):

	int ***array = allocate3d_int(7, 5, 1);
	printf("%d\n", array[6][4][0]);
	
	@param[in] nz slowest incrimenting dimension of array in memory.
	@param[in] ny middle dimension of array.
	@param[in] nx fastest incrimenting dimension of array.

	@return A malloc'd array or NULL if any malloc fails.
*/
int ***allocate3d_int(size_t nz, size_t ny, size_t nx);
/** @brief Free a three-dimensional array created by allocate3d_int.

	@param[in] p Pointer to array created by allocate3d_int.
*/
void free3d_int(int ***p);

/** @brief Allocate a three-dimensional array of type float of a given size.

	The order of the parameters are (now) the same as when the data is being accessed.
        E.g., this is valid (ignoring the printf type specifier):

	float ***array = allocate3d_float(7, 5, 1);
	printf("%d\n", array[6][4][0]);
	
	@param[in] nz slowest incrimenting dimension of array in memory.
	@param[in] ny middle dimension of array.
	@param[in] nx fastest incrimenting dimension of array.

	@return A malloc'd array or NULL if any malloc fails.
*/
float ***allocate3d_float(size_t nz, size_t ny, size_t nx);
/** @brief Free a three-dimensional array created by allocate3d_float.

	@param[in] p Pointer to array created by allocate3d_float.
*/
void free3d_float(float ***p);

/** @brief Allocate a three-dimensional array of type double of a given size.

	The order of the parameters are (now) the same as when the data is being accessed.
        E.g., this is valid (ignoring the printf type specifier):

	double ***array = allocate3d_double(7, 5, 1);
	printf("%d\n", array[6][4][0]);
	
	@param[in] nz slowest incrimenting dimension of array in memory.
	@param[in] ny middle dimension of array.
	@param[in] nx fastest incrimenting dimension of array.

	@return A malloc'd array or NULL if any malloc fails.
*/
double ***allocate3d_double(size_t nz, size_t ny, size_t nx);
/** @brief Free a three-dimensional array created by allocate3d_double.

	@param[in] p Pointer to array created by allocate3d_double.
*/
void free3d_double(double ***p);

#ifdef __cplusplus
}
#endif

#endif /* OEL_UTIL_LIBGENUTILS_ALLOCATE3D_H_ */
