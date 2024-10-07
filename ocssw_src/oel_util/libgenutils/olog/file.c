#include "olog.h"
#include "olog/file.h"

#include <stdio.h>
#include <stdlib.h>

static int olog_backend_file_print(olog *olog, olog_backend *backend, uint8_t severity, va_list args) {
    const char* format = va_arg(args, const char*);
    int ret = vfprintf((FILE*) backend->data, format, args);
    return ret < 0 ? OLOG_EXPLODED : 0;
}

static int olog_backend_file_destroy(olog *olog, olog_backend *backend) {
    if (backend->data) {
        fclose(backend->data);
    }
    return 0;
}
static int olog_backend_has_error(olog *olog, olog_backend *backend) {
    return backend->data == NULL;
}
static int olog_backend_print_error(olog *olog, olog_backend *backend, FILE* stream) {
    if (backend->data == NULL) {
        fprintf(stream, "Couldn't open file handle\n");
        return 1;
    }
    return 0;
}

olog_backend olog_backend_file_mode(const char *filename, const char *mode, int8_t min_log_level, int8_t max_log_level) {
    FILE *stream = fopen(filename, mode);
    olog_backend ret = {
            .min_log_level = min_log_level,
            .max_log_level = max_log_level,
            .data = stream,
            .print_callback = olog_backend_file_print,
            .destroy_callback = olog_backend_file_destroy,
            .has_error_callback = olog_backend_has_error,
            .print_error_callback = olog_backend_print_error
    };
    return ret;
}
olog_backend olog_backend_file(const char *filename, int8_t min_log_level, int8_t max_log_level) {
    return olog_backend_file_mode(filename, "w", min_log_level, max_log_level);
}
