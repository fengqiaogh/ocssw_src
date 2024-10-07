#include "olog.h"
#include "olog/buffer.h"

#include <stdio.h>
#include <stdlib.h>

#define BUFFER_INCREASE 2

static int olog_backend_buffer_print(olog *olog, olog_backend *backend, uint8_t severity, va_list args) {
    olog_buffer *data = (olog_buffer*) (backend->data);
    if (!data || !data->buffer_start) {
        return OLOG_EXPLODED;
    }
    va_list args_tmp;
    va_copy(args_tmp, args);

    const char* format = va_arg(args_tmp, const char*);
    size_t max_to_print = data->buffer_size_max - data->buffer_size;
    int ret = vsnprintf(data->buffer_current, max_to_print, format, args_tmp);
    if (ret < 0) {
        return OLOG_EXPLODED;
    }
    if ((size_t) ret >= max_to_print) {
        void *tmp_ptr;
        tmp_ptr = realloc(data->buffer_start, data->buffer_size_max * BUFFER_INCREASE);
        if (tmp_ptr) {
            data->buffer_start = tmp_ptr;
            data->buffer_current = data->buffer_start + data->buffer_size;
            data->buffer_size_max *= BUFFER_INCREASE;
            return olog_backend_buffer_print(olog, backend, severity, args);
        } else {
            free(data->buffer_start);
            data->buffer_start = NULL;
            data->buffer_current = NULL;
            data->buffer_size_max = 0;
            data->buffer_size = 0;
            return 1;
        }
    } else {
        data->buffer_current += ret;
        data->buffer_size += ret;
    }
    return ret < 0 ? OLOG_EXPLODED : 0;
}

static int olog_backend_buffer_destroy(olog *olog, olog_backend *backend) {
    if (backend->data) {
        olog_buffer *data = (olog_buffer*) (backend->data);
        if (data->buffer_start) {
            free(data->buffer_start);
        }
        free(backend->data);
    }
    return 0;
}
static int olog_backend_has_error(olog *olog, olog_backend *backend) {
    if (backend->data == NULL) {
        return 1;
    }
    olog_buffer *data = (olog_buffer*) (backend->data);
    return data->buffer_start == NULL;
}
static int olog_backend_print_error(olog *olog, olog_backend *backend, FILE* stream) {
    if (backend->data == NULL) {
        fprintf(stream, "Couldn't malloc object\n");
        return 1;
    }
    olog_buffer *data = (olog_buffer*) (backend->data);
    if (data->buffer_start == NULL) {
        fprintf(stream, "Couldn't malloc buffer\n");
        return 1;
    }
    return 0;
}

olog_backend olog_backend_buffer(size_t initial_size, int8_t min_log_level, int8_t max_log_level) {
    olog_buffer *data = malloc(sizeof(olog_buffer));
    data->buffer_start = malloc(sizeof(char) * initial_size);
    data->buffer_current = data->buffer_start;
    data->buffer_size_max = initial_size;
    data->buffer_size = 0;
    *data->buffer_start = '\0';
    olog_backend ret = {
            .min_log_level = min_log_level,
            .max_log_level = max_log_level,
            .data = data,
            .print_callback = olog_backend_buffer_print,
            .destroy_callback = olog_backend_buffer_destroy,
            .has_error_callback = olog_backend_has_error,
            .print_error_callback = olog_backend_print_error
    };
    return ret;
}

void olog_backend_buffer_clear(olog_backend *backend) {
    olog_buffer *data = (olog_buffer*) backend->data;
    data->buffer_current = data->buffer_start;
    data->buffer_size = 0;
    *data->buffer_start = '\0';
}
