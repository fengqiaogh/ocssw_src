#include "olog.h"
#include "olog/stream.h"

#include <stdio.h>
#include <stdlib.h>

static int olog_backend_stream_print(olog *olog, olog_backend *backend, uint8_t severity, va_list args) {
    const char* format = va_arg(args, const char*);
    int ret = vfprintf((FILE*) backend->data, format, args);
    return ret < 0 ? OLOG_EXPLODED : 0;
}

olog_backend olog_backend_stream(FILE *stream, int8_t min_log_level, int8_t max_log_level) {
    olog_backend ret = {
            .min_log_level = min_log_level,
            .max_log_level = max_log_level,
            .data = stream,
            .print_callback = olog_backend_stream_print,
            .destroy_callback = NULL
    };
    return ret;
}
