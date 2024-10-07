// implemented version numbers
#define SHASH_IMPLEMENTATION 1002000
#define SHASH_IMPLEMENTATION_STR "1.2.0"
#define SHASH_API_VERSION 1002000
#define SHASH_API_VERSION_STR "1.2.0"
#include "shash.h"

#include <limits.h>		/* For CHAR_BIT.  */
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/** @brief Size to initialize shash to upon creation.

 This has to be 2^n due to some bit math cheating ("index & size-1" to do modulus
 division to wrap the target bucket).
 */
#ifndef SHASH_INITIAL_SIZE
#define SHASH_INITIAL_SIZE 64
#endif

/** @brief How full to get before doubling the size of the shash, 0-99.

 This counts values, not actual buckets filled.
 */
#ifndef SHASH_FILLPCT
#define SHASH_FILLPCT 66
#endif

typedef struct shash_bucket shash_bucket;

struct shash_bucket {
    shash_bucket *next;
    char *key;
    char *value;
    uint32_t shash;
};

struct shash {
    shash_bucket **buckets;
    int max_bucket;
    int filled;
    int iter_i;
    shash_bucket *iter;
};

// Murmur3 32-bit hash stolen from qlibc
// https://github.com/wolkykim/qlibc/blob/master/src/utilities/qhash.c
static uint32_t compute_hash(const char *key) {
    if (key == NULL)
        return 0;
    size_t key_length = strlen(key);
    if (key_length == 0) {
        return 0;
    }
    const uint32_t c1 = 0xcc9e2d51;
    const uint32_t c2 = 0x1b873593;
    const unsigned nblocks = (unsigned) (key_length / 4);
    const uint32_t *blocks = (const uint32_t *) (key);
    const uint8_t *tail = (const uint8_t *) (key + (nblocks * 4));
    uint32_t h = 0;
    unsigned i;
    uint32_t k;
    for (i = 0; i < nblocks; i++) {
        k = blocks[i];
        k *= c1;
        k = (k << 15) | (k >> (32 - 15));
        k *= c2;
        h ^= k;
        h = (h << 13) | (h >> (32 - 13));
        h = (h * 5) + 0xe6546b64;
    }
    k = 0;
    switch (key_length & 3) {
    case 3:
        k ^= tail[2] << 16;
        /* no break */
    case 2:
        k ^= tail[1] << 8;
        /* no break */
    case 1:
        k ^= tail[0];
        k *= c1;
        k = (k << 15) | (k >> (32 - 15));
        k *= c2;
        h ^= k;
    };
    h ^= (uint32_t) key_length;
    h ^= h >> 16;
    h *= 0x85ebca6b;
    h ^= h >> 13;
    h *= 0xc2b2ae35;
    h ^= h >> 16;
    return h;
}

static void shash_zero_out(shash_bucket **buckets, int start, int end) {
    int i;
    for (i = start; i <= end; i++) {
        buckets[i] = NULL;
    }
}

shash *shash_create(uint32_t options) {
    (void) options;
    shash *h = malloc(sizeof(shash));
    if (h == NULL) {
        return NULL;
    }
    h->buckets = (shash_bucket**) malloc(SHASH_INITIAL_SIZE * sizeof(shash_bucket*));
    if (h->buckets == NULL) {
        free(h);
        return NULL;
    }
    h->filled = 0;
    h->max_bucket = SHASH_INITIAL_SIZE - 1;
    shash_rewind(h);
    shash_zero_out(h->buckets, 0, h->max_bucket);
    return h;
}

static int shash_destroy_bucket(shash_bucket *bucket, bool follow) {
    if (bucket != NULL) {
        free(bucket->key);
        free(bucket->value);
        int ret = 1;
        if (follow) {
            ret += shash_destroy_bucket(bucket->next, true);
        }
        free(bucket);
        return ret;
    }
    return 0;
}
int shash_destroy(shash *h) {
    int i, total_freed = 0;
    for (i = 0; i <= h->max_bucket; i++) {
        total_freed += shash_destroy_bucket(h->buckets[i], true);
    }
    free(h->buckets);
    free(h);
    return 0;
}

int shash_remove(shash *h, const char *key) {
    if (h == NULL) {
        return -1;
    }
    uint32_t shash = compute_hash(key);
    shash_bucket **bucket_ptr = &(h->buckets[shash & h->max_bucket]);
    shash_bucket *cur_bucket = *bucket_ptr;

    int is_first;
    for (is_first = 1; cur_bucket != NULL; is_first = 0, bucket_ptr = &cur_bucket->next, cur_bucket = *bucket_ptr) {
        if (cur_bucket->shash == shash && strcmp(cur_bucket->key, key) == 0) {
            *bucket_ptr = cur_bucket->next;
            shash_destroy_bucket(cur_bucket, false);
            if (is_first) {
                h->filled--;
            }
            return 0;
        }
    }
    return 1;
}

static int shash_split(shash *h) {
    int oldsize = h->max_bucket + 1;
    int newsize = oldsize << 1;
    shash_bucket **a = (shash_bucket**) realloc((shash_bucket*) h->buckets, newsize * sizeof(shash_bucket*));
    if (a == NULL) {
        return -1;
    }
    shash_zero_out(a, oldsize, --newsize);

    h->max_bucket = newsize;
    h->buckets = a;

    shash_bucket **b, **oentry, *entry;

    int i;
    for (i = 0; i < oldsize; i++, a++) {
        if (*a == NULL)
            continue;
        b = a + oldsize;
        for (oentry = a, entry = *a; entry != NULL; entry = *oentry) {
            if ((entry->shash & newsize) != (unsigned) i) {
                *oentry = entry->next;
                entry->next = *b;
//				if (*b == NULL)
//					h->filled++;
                *b = entry;
                continue;
            } else {
                oentry = &entry->next;
            }
        }
//		if (*a == NULL)
//			h->filled--;
    }
    return 0;
}

const char *shash_get(shash *h, const char *key) {
    shash_bucket *entry;
    shash_bucket **oentry;

    if (h == NULL) {
        return NULL;
    }

    uint32_t shash = compute_hash(key);
    oentry = &(h->buckets[shash & h->max_bucket]);

    for (entry = *oentry; entry != NULL; entry = entry->next) {
        if (entry->shash == shash && strcmp(entry->key, key) == 0) {
            return entry->value;
        }
    }
    return NULL;
}

int shash_set(shash *h, const char *key, const char *value) {
    shash_bucket *entry;
    shash_bucket **oentry;

    if (h == NULL) {
        return -1;
    }

    uint32_t shash = compute_hash(key);
    oentry = &(h->buckets[shash & h->max_bucket]);

    for (entry = *oentry; entry != NULL; entry = entry->next) {
        if (entry->shash == shash && strcmp(entry->key, key) == 0) {
            free(entry->value);
            entry->value = strdup(value);
            if (entry->value == NULL) {
                shash_remove(h, key);
                return -1;
            }
            return 1;
        }
    }
    entry = (shash_bucket*) malloc(sizeof(shash_bucket));
    if (entry == NULL) {
        return -1;
    }

    entry->key = strdup(key);
    entry->value = strdup(value);
    if (entry->key == NULL || entry->value == NULL) {
        if (entry->key == NULL) {
            free(entry->key);
        }
        if (entry->value == NULL) {
            free(entry->value);
        }
        free(entry);
        return -1;
    }
    entry->shash = shash;
    entry->next = *oentry;
    *oentry = entry;

    h->filled++;
    if ((h->filled * 100 / (h->max_bucket + 1)) > SHASH_FILLPCT) {
        if (shash_split(h)) {
            return -1;
        }
    }

    return 0;
}

int shash_rewind(shash *h) {
    h->iter_i = -1;
    h->iter = NULL;
    return h->filled;
}

int shash_next(shash *h, const char **key, const char **value) {
    if (h == NULL) {
        return -1;
    }
    shash_bucket *b = h->iter;
    do {
        if (b != NULL) {
            b = b->next;
        }
        if (b == NULL) {
            h->iter_i++;
            if (h->iter_i > h->max_bucket) {
                h->iter_i = -1;
                break;
            }
            b = h->buckets[h->iter_i];
        }
    } while (b == NULL);
    h->iter = b;
    if (b == NULL) {
        if (key != NULL) {
            *key = NULL;
        }
        if (value != NULL){
            *value = NULL;
        }
        return 1;
    }
    if (key != NULL) {
        *key = b->key;
    }
    if (value != NULL){
        *value = b->value;
    }
    return 0;
}

int shash_size(shash *h){
    return h->filled;
}

const char *shash_version() {
    return "shash " SHASH_IMPLEMENTATION_STR ", API " SHASH_API_VERSION_STR;
}
