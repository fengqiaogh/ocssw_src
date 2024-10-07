// implemented version numbers
#define PHASH_IMPLEMENTATION 1002001
#define PHASH_IMPLEMENTATION_STR "1.2.0"
#define PHASH_API_VERSION 1002000
#define PHASH_API_VERSION_STR "1.2.0"
#include "phash.h"

#include <limits.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/** @brief Size to initialize hash to upon creation.

 This has to be 2^n due to some bit math cheating ("index & size-1" to do modulus
 division to wrap the target bucket).
 */
#ifndef PHASH_INITIAL_SIZE
#define PHASH_INITIAL_SIZE 64
#endif

/** @brief How full to get before doubling the size of the hash, 0-99.

 This counts values, not actual buckets filled.
 */
#ifndef PHASH_FILLPCT
#define PHASH_FILLPCT 66
#endif

typedef struct phash_bucket phash_bucket;

struct phash_bucket {
    phash *owner;
    phash_bucket *next;
    const char *key;
    void *value;
    uint32_t hash;
};

struct phash {
    phash_bucket **buckets;
    int max_bucket;
    int filled;
    int iter_i;
    unsigned options;
    phash_bucket *iter;
};

// Murmur3 32-bit hash stolen from qlibc
// https://github.com/wolkykim/qlibc/blob/master/src/utilities/qhash.c
static uint32_t compute_hash(const char *key) {
    if (key == NULL){
        return 0;
	}
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

static void phash_zero_out(phash_bucket **buckets, int start, int end) {
    int i;
    for (i = start; i <= end; i++) {
        buckets[i] = NULL;
    }
}

phash *phash_create(uint32_t options) {
    phash *h = malloc(sizeof(phash));
    if (h == NULL) {
        return NULL;
    }
    h->options = options;
    h->buckets = (phash_bucket**) malloc(PHASH_INITIAL_SIZE * sizeof(phash_bucket*));
    if (h->buckets == NULL) {
        free(h);
        return NULL;
    }
    h->filled = 0;
    h->max_bucket = PHASH_INITIAL_SIZE - 1;
    phash_rewind(h);
    phash_zero_out(h->buckets, 0, h->max_bucket);
    return h;
}

static int phash_destroy_bucket(phash_bucket *bucket) {
    if (bucket != NULL) {
        int ret = phash_destroy_bucket(bucket->next);
        if (!(bucket->owner->options & PHASH_NO_COPY_KEYS)) {
            free((char*) bucket->key);
        }
        free(bucket);
        return ret + 1;
    }
    return 0;
}
int phash_destroy(phash *h) {
    int i, total_freed = 0;
    for (i = 0; i <= h->max_bucket; i++) {
        total_freed += phash_destroy_bucket(h->buckets[i]);
    }
    free(h->buckets);
    free(h);
    return 0;
}

int phash_remove(phash *h, const char *key) {
    if (h == NULL) {
        return -1;
    }
    uint32_t hash = compute_hash(key);
    phash_bucket **bucket_ptr = &(h->buckets[hash & h->max_bucket]);
    phash_bucket *cur_bucket = *bucket_ptr;

    int is_first;
    for (is_first = 1; cur_bucket != NULL; is_first = 0, bucket_ptr = &cur_bucket->next, cur_bucket = *bucket_ptr) {
        if (cur_bucket->hash == hash && strcmp(cur_bucket->key, key) == 0) {
            *bucket_ptr = cur_bucket->next;
            free(cur_bucket);
            if (is_first) {
                h->filled--;
            }
            return 0;
        }
    }
    return 1;
}

static int phash_split(phash *h) {
    int oldsize = h->max_bucket + 1;
    int newsize = oldsize << 1;
    phash_bucket **a = (phash_bucket**) realloc((phash_bucket*) h->buckets, newsize * sizeof(phash_bucket*));
    if (a == NULL) {
        return -1;
    }
    phash_zero_out(a, oldsize, --newsize);

    h->max_bucket = newsize;
    h->buckets = a;

    phash_bucket **b, **oentry, *entry;

    int i;
    for (i = 0; i < oldsize; i++, a++) {
        if (*a == NULL){
            continue;
		}
        b = a + oldsize;
        for (oentry = a, entry = *a; entry != NULL; entry = *oentry) {
            if ((entry->hash & newsize) != (unsigned) i) {
                *oentry = entry->next;
                entry->next = *b;
//				if (*b == NULL)
//					h->filled++;
                *b = entry;
            } else {
                oentry = &entry->next;
            }
        }
//		if (*a == NULL)
//			h->filled--;
    }
    return 0;
}

void *phash_get(phash *h, const char *key) {
    phash_bucket *entry;
    phash_bucket **oentry;

    if (h == NULL) {
        return NULL;
    }

    uint32_t hash = compute_hash(key);
    oentry = &(h->buckets[hash & h->max_bucket]);

    for (entry = *oentry; entry != NULL; entry = entry->next) {
        if (entry->hash == hash && strcmp(entry->key, key) == 0) {
            return entry->value;
        }
    }
    return NULL;
}

int phash_set(phash *h, const char *key, void *value) {
    register phash_bucket *entry;
    register phash_bucket **oentry;

    if (h == NULL) {
        return -1;
    }

    uint32_t hash = compute_hash(key);
    oentry = &(h->buckets[hash & h->max_bucket]);

    for (entry = *oentry; entry != NULL; entry = entry->next) {
        if (entry->hash == hash && strcmp(entry->key, key) == 0) {
            entry->value = value;
            return 1;
        }
    }
    entry = (phash_bucket*) malloc(sizeof(phash_bucket));
    if (entry == NULL) {
        return 1;
    }

    entry->owner = h;
    if (!(h->options & PHASH_NO_COPY_KEYS)) {
        size_t key_size = strlen(key) + 1;
        entry->key = malloc(key_size * sizeof(char));
        strncpy((char*) entry->key, key, key_size);
    } else {
        entry->key = key;
    }
    entry->value = value;
    entry->hash = hash;
    entry->next = *oentry;
    *oentry = entry;

    h->filled++;
    if ((h->filled * 100 / (h->max_bucket + 1)) > PHASH_FILLPCT) {
        if (phash_split(h)) {
            return -1;
        }
    }

    return 0;
}

int phash_rewind(phash *h) {
    h->iter_i = -1;
    h->iter = NULL;
    return h->filled;
}

int phash_next(phash *h, const char **key, void **value) {
    if (h == NULL) {
        return -1;
    }
    phash_bucket *b = h->iter;
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
        if (value != NULL) {
            *value = NULL;
        }
        return 1;
    }
    if (key != NULL) {
        *key = b->key;
    }
    if (value != NULL) {
        *value = b->value;
    }
    return 0;
}

int phash_size(phash *h){
    return h->filled;
}

const char *phash_version() {
    return "phash " PHASH_IMPLEMENTATION_STR ", API " PHASH_API_VERSION_STR;
}
