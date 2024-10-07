/** @file phash.h
 @brief A simple dictionary library for storing pointers.

 A phash doesn't actually store any data for the key or the value, merely
 pointers to them.  The in-memory data must persist until the phash is
 destroyed.  If a memory location referred to by a phash is free'd, or
 non-malloc'd pointers simply go out of scope the phash will almost definitely
 try to access invalid memory locations during a call to any function except
 for phash_destroy.  String literals, depending on the compiler, generally
 get put into some kind of read-only memory and persist through the program,
 so they might be safe.
 */

#ifndef __PHASH_H__
#define __PHASH_H__

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>

#ifndef PHASH_API_VERSION
/** @brief Current API version.  If already defined, use API defined first. */
#define PHASH_API_VERSION 1002000
/** @brief Current API version as a string.  If already defined, use API defined first. */
#define PHASH_API_VERSION_STR "1.2.0"
#endif

#define PHASH_NO_COPY_KEYS 1

/** @brief Implementation-specific, generic type to store the phash object */
typedef struct phash phash;

/** @brief Initialize a phash object.

 @param[in] options Bitwise OR'd option flags.  Currently unused.

 @return Pointer to malloc'd phash.
 */
phash *phash_create(uint32_t options);

/** @brief Destroy a phash object, free'ing the memory used.

 @param[in] h Pointer to phash to be destroyed.

 @return 0 on success.
 */
int phash_destroy(phash *h);

/** @brief Add or overwrite a pointer, associating it with the given key.

 @param[in] h Pointer to phash to be modified.
 @param[in] key Null-terminated string with which to associate the new pointer.
 @param[in] value New pointer to store..

 @return 0 on add, 1 on overwrite, and -1 on error.
 */
int phash_set(phash *h, const char *key, void *value);

/** @brief Find a pointer associated with the given string.


 @param[in] h Pointer to phash to be searched.
 @param[in] key Key for which to search.

 @return Pointer to value, or NULL if key is not found.
 */
void *phash_get(phash *h, const char *key);

/** @brief Remove a pointer associated with the given string.

 @param[in] h Pointer to phash to be modified.
 @param[in] key Key for which to remove.

 @return 0 on success.
 */
int phash_remove(phash *h, const char *key);

/** @brief Rewind iterator for traversing all the keys and values.

 @param[in] h Pointer to phash to be reset.

 @return number of items currently in the phash.
 */
int phash_rewind(phash *h);

/** @brief Retrieves the next key-value pair in the phash.  The order in which
 the pointers are returned shall be consistent but may be arbitrary.

 Consistent refers only to running identical code and rewinding and traversing the data.
 Adding or removing pointers may modify the order in which keys are returned.

 On end-of-hash, the hash is rewind()'d and the state of the key and value inputs
 are set to NULL.

 @param[in] h Pointer to phash to be traversed.
 @param[out] key Pointer to fill with the next key, NULL to ignore.
 @param[out] value Pointer to fill with the next value, NULL to ignore.

 @return -1 on error, positive 1 on end-of-hash, and 0 on success.
 */
int phash_next(phash *h, const char **key, void **value);

/** @brief Retrieves the number of key-value pairs in the hash object.

 @param[in] h Pointer to shash to be analyzed.
 */
int phash_size(phash *h);

/** @brief Returns the source code version and the implemented API version.

 No assumptions are made about the format of the return value.
 */
const char *phash_version();

#ifdef __cplusplus
}
#endif

#endif // __PHASH_H__
