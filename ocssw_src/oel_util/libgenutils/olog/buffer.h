
#ifndef __OLOG_BUFFER_H_
#define __OLOG_BUFFER_H_

#include "olog.h"
#include "phash.h"

#include <stdio.h>
#include <stdlib.h>

typedef struct olog_buffer {
	char *buffer_start, *buffer_current;
	size_t buffer_size_max, buffer_size;
} olog_buffer;

olog_backend olog_backend_buffer(size_t initial_size, int8_t min_log_level, int8_t max_log_level);
void olog_backend_buffer_clear(olog_backend *backend);

#endif /* __OLOG_BUFFER_H_ */
