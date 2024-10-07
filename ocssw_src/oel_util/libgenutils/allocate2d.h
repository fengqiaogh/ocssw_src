/** @file allocate2d.h
	@brief Utility functions for allocating and freeing two-dimensional arrays of various types.

	This file was created by allocate2d.pl and should not be edited manually.
*/
#ifndef OEL_UTIL_LIBGENUTILS_ALLOCATE2D_H_
#define OEL_UTIL_LIBGENUTILS_ALLOCATE2D_H_

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/** @brief Allocate a two-dimensional array of type char of a given size.

	The order of the parameters are (now) the same as when the data is being accessed.
        E.g., this is valid (ignoring the printf type specifier):

	char **array = allocate2d_char(5, 1);
	printf("%d\n", array[4][0]);
	
	@param[in] h Height of array.
	@param[in] w Width of array.

	@return A malloc'd array or NULL if any malloc fails.
*/
char **allocate2d_char(size_t h, size_t w);
/** @brief Free a two-dimensional array created by allocate2d_char.

	@param[in] p Pointer to array created by allocate2d_char.
*/
void free2d_char(char **p);

/** @brief Allocate a two-dimensional array of type unsigned char of a given size.

	The order of the parameters are (now) the same as when the data is being accessed.
        E.g., this is valid (ignoring the printf type specifier):

	unsigned char **array = allocate2d_uchar(5, 1);
	printf("%d\n", array[4][0]);
	
	@param[in] h Height of array.
	@param[in] w Width of array.

	@return A malloc'd array or NULL if any malloc fails.
*/
unsigned char **allocate2d_uchar(size_t h, size_t w);
/** @brief Free a two-dimensional array created by allocate2d_uchar.

	@param[in] p Pointer to array created by allocate2d_uchar.
*/
void free2d_uchar(unsigned char **p);

/** @brief Allocate a two-dimensional array of type signed char of a given size.

	The order of the parameters are (now) the same as when the data is being accessed.
        E.g., this is valid (ignoring the printf type specifier):

	signed char **array = allocate2d_schar(5, 1);
	printf("%d\n", array[4][0]);
	
	@param[in] h Height of array.
	@param[in] w Width of array.

	@return A malloc'd array or NULL if any malloc fails.
*/
signed char **allocate2d_schar(size_t h, size_t w);
/** @brief Free a two-dimensional array created by allocate2d_schar.

	@param[in] p Pointer to array created by allocate2d_schar.
*/
void free2d_schar(signed char **p);

/** @brief Allocate a two-dimensional array of type short of a given size.

	The order of the parameters are (now) the same as when the data is being accessed.
        E.g., this is valid (ignoring the printf type specifier):

	short **array = allocate2d_short(5, 1);
	printf("%d\n", array[4][0]);
	
	@param[in] h Height of array.
	@param[in] w Width of array.

	@return A malloc'd array or NULL if any malloc fails.
*/
short **allocate2d_short(size_t h, size_t w);
/** @brief Free a two-dimensional array created by allocate2d_short.

	@param[in] p Pointer to array created by allocate2d_short.
*/
void free2d_short(short **p);

/** @brief Allocate a two-dimensional array of type int of a given size.

	The order of the parameters are (now) the same as when the data is being accessed.
        E.g., this is valid (ignoring the printf type specifier):

	int **array = allocate2d_int(5, 1);
	printf("%d\n", array[4][0]);
	
	@param[in] h Height of array.
	@param[in] w Width of array.

	@return A malloc'd array or NULL if any malloc fails.
*/
int **allocate2d_int(size_t h, size_t w);
/** @brief Free a two-dimensional array created by allocate2d_int.

	@param[in] p Pointer to array created by allocate2d_int.
*/
void free2d_int(int **p);

/** @brief Allocate a two-dimensional array of type float of a given size.

	The order of the parameters are (now) the same as when the data is being accessed.
        E.g., this is valid (ignoring the printf type specifier):

	float **array = allocate2d_float(5, 1);
	printf("%d\n", array[4][0]);
	
	@param[in] h Height of array.
	@param[in] w Width of array.

	@return A malloc'd array or NULL if any malloc fails.
*/
float **allocate2d_float(size_t h, size_t w);
/** @brief Free a two-dimensional array created by allocate2d_float.

	@param[in] p Pointer to array created by allocate2d_float.
*/
void free2d_float(float **p);

/** @brief Allocate a two-dimensional array of type double of a given size.

	The order of the parameters are (now) the same as when the data is being accessed.
        E.g., this is valid (ignoring the printf type specifier):

	double **array = allocate2d_double(5, 1);
	printf("%d\n", array[4][0]);
	
	@param[in] h Height of array.
	@param[in] w Width of array.

	@return A malloc'd array or NULL if any malloc fails.
*/
double **allocate2d_double(size_t h, size_t w);
/** @brief Free a two-dimensional array created by allocate2d_double.

	@param[in] p Pointer to array created by allocate2d_double.
*/
void free2d_double(double **p);

#ifdef __cplusplus
}
#endif

#endif /* OEL_UTIL_LIBGENUTILS_ALLOCATE2D_H_ */
