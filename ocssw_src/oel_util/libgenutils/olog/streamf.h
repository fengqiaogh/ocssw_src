
#ifndef __OLOG_STREAMF_H_
#define __OLOG_STREAMF_H_

#include "olog.h"
#include "phash.h"

#include <stdio.h>
#include <stdlib.h>

typedef struct olog_streamf_arguments {
	char *filename, *file_mode;
	FILE *stream;
	char *format, *strftime_format;
	size_t strftime_size;
	bool persistent_hash, persistent_indexes;
	phash *named_pointers;
	void **indexed_pointers;
} olog_streamf_arguments;

// TODO: figure out how not to expose this stuff

typedef int (*olog_print_callback) (olog *olog, olog_backend *backend, uint8_t severity, va_list args);
typedef int (*olog_print_callback_no_vargs) (olog *olog, olog_backend *backend, uint8_t severity);
typedef int (*olog_print_callback_no_args) (olog *olog, olog_backend *backend);

enum streamf_type_t {
	STREAMF_TYPE_STRING = 0,
	STREAMF_TYPE_FUNCTION = 1,
	STREAMF_TYPE_FUNCTION_NO_VARGS = 2,
	STREAMF_TYPE_FUNCTION_NO_ARGS = 3,

	STREAMF_TYPE_LOOKUP_HASH = 16,
	STREAMF_TYPE_LOOKUP_ARRAY = 32
};
typedef enum streamf_type_t streamf_type;
typedef struct print_command {
	streamf_type type;
	void *data;
} print_command;


typedef struct olog_streamf {
	char *filename;
	FILE *stream;
	char *format, *strftime_format;
	size_t strftime_size;
	phash *named_pointers;
	void **indexed_pointers;
	print_command *commands;
	unsigned command_count;
} olog_streamf;

olog_backend olog_backend_streamf(olog_streamf_arguments *arguments, int8_t min_log_level, int8_t max_log_level);

#endif /* __OLOG_STREAMF_H_ */
