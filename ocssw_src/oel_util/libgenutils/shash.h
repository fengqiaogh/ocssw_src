/** @file shash.h
 @brief A simple dictionary library for storing strings.
 */

#ifndef __SHASH_H__
#define __SHASH_H__

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>

#ifndef SHASH_API_VERSION
/** @brief Current API version.  If already defined, use API defined first. */
#define SHASH_API_VERSION 1002000
/** @brief Current API version as a string.  If already defined, use API defined first. */
#define SHASH_API_VERSION_STR "1.2.0"
#endif

/** @brief Implementation-specific, generic type to store the shash object */
typedef struct shash shash;

/** @brief Initialize a shash object.

 @param[in] options Bitwise OR'd option flags.  Currently unused.

 @return Pointer to malloc'd shash.
 */
shash *shash_create(uint32_t options);

/** @brief Destroy a shash object, free'ing the memory used.

 @param[in] h Pointer to shash to be destroyed.

 @return 0 on success.
 */
int shash_destroy(shash *h);

/** @brief Add or overwrite a pointer, associating it with the given key.

 @param[in] h Pointer to shash to be modified.
 @param[in] key Null-terminated string with which to associate the new pointer.
 @param[in] value New pointer to store..

 @return 0 on add, 1 on overwrite, and -1 on error.
 */
int shash_set(shash *h, const char *key, const char *value);

/** @brief Find a pointer associated with the given string.


 @param[in] h Pointer to shash to be searched.
 @param[in] key Key for which to search.

 @return Pointer to value, or NULL if key is not found.

 */
const char *shash_get(shash *h, const char *key);

/** @brief Remove a pointer associated with the given string.

 @param[in] h Pointer to shash to be modified.
 @param[in] key Key for which to remove.

 @return 0 on success.
 */
int shash_remove(shash *h, const char *key);

/** @brief Rewind iterator for traversing all the keys and values.

 @param[in] h Pointer to shash to be reset.

 @return number of items currently in the shash.
 */
int shash_rewind(shash *h);

/** @brief Retrieves the next key-value pair in the shash.  The order in which
 the pointers are returned shall be consistent but may be arbitrary.

 Consistent refers only to running identical code and rewinding and traversing the data.
 Adding or removing pointers may modify the order in which keys are returned.

 On end-of-shash, the shash is rewind()'d and the state of the key and value inputs
 are set to NULL.

 @param[in] h Pointer to shash to be traversed.
 @param[out] key Pointer to fill with the next key, NULL to ignore.
 @param[out] value Pointer to fill with the next value, NULL to ignore.

 @return -1 on error, positive 1 on end-of-shash, and 0 on success.
 */
int shash_next(shash *h, const char **key, const char **value);

/** @brief Retrieves the number of key-value pairs in the hash object.

 @param[in] h Pointer to shash to be analyzed.
 */
int shash_size(shash *h);

/** @brief Returns the source code version and the implemented API version.

 No assumptions are made about the format of the return value.
 */
const char *shash_version();

#ifdef __cplusplus
}
#endif

#endif // __SHASH_H__
