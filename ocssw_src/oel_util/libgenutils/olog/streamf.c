#include "olog.h"
#include "olog/streamf.h"

#include "phash.h"

#include <errno.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/*
 API notes

 %% - literal %
 %{function_name}f - function from hash
 %{function_name}F - function from hash (no vargs)
 %{string_name}s - string from hash
 %[function_index]f - function from array
 %[function_name]F - function from array (no vargs)
 %[string_index]s - string from array

 %s - printf arguments
 %l - severity lowercase
 %u - severity uppercase
 %m - severity mixed case
 %n - severity number

 %t - time
 time string to pass to strftime passed separately

 */

const char *severity_l[] = { "debug", "info", "notice", "warning", "error", "critical", "alert", "emergency" };
const char *severity_u[] = { "DEBUG", "INFO", "NOTICE", "WARNING", "ERROR", "CRITICAL", "ALERT", "EMERGENCY" };
const char *severity_m[] = { "Debug", "Info", "Notice", "Warning", "Error", "Critical", "Alert", "Emergency" };

static int olog_backend_streamf_printf(olog *olog, olog_backend *backend, uint8_t severity, va_list args) {

//	STREAMF_TYPE_STRING = 0,
//	STREAMF_TYPE_FUNCTION = 1,
//	STREAMF_TYPE_FUNCTION_NO_VARGS = 2,
//	STREAMF_TYPE_FUNCTION_NO_ARGS = 3,
//
//	STREAMF_TYPE_LOOKUP_HASH = 16,
//	STREAMF_TYPE_LOOKUP_ARRAY = 32

//	const char* format = va_arg(args, const char*);
//	int ret = vfprintf((FILE*)((olog_streamf*)backend->data)->stream, format, args);

    olog_streamf *backend_data = (olog_streamf*) backend->data;
    int ret = 0;
    for (unsigned i = 0; i < backend_data->command_count && !ret; i++) {
        streamf_type type = backend_data->commands[i].type;
        void *data = backend_data->commands[i].data;

        if (type & STREAMF_TYPE_LOOKUP_HASH) {
            type ^= STREAMF_TYPE_LOOKUP_HASH;
            data = phash_get(backend_data->named_pointers, data);
        } else if (type & STREAMF_TYPE_LOOKUP_ARRAY) {
            type ^= STREAMF_TYPE_LOOKUP_ARRAY;
            data = backend_data->indexed_pointers[(unsigned long long) data];
        }

        switch (type) {
        case STREAMF_TYPE_STRING:
            ret = (fprintf((FILE*) backend_data->stream, "%s", (char*) data) < 0 ? OLOG_EXPLODED : 0);
            break;
        case STREAMF_TYPE_FUNCTION:
            ret = (((olog_print_callback) data)(olog, backend, severity, args) < 0 ? OLOG_EXPLODED : 0);
            break;
        case STREAMF_TYPE_FUNCTION_NO_VARGS:
            ret = (((olog_print_callback_no_vargs) data)(olog, backend, severity) < 0 ? OLOG_EXPLODED : 0);
            break;
        case STREAMF_TYPE_FUNCTION_NO_ARGS:
            ret = (((olog_print_callback_no_args) data)(olog, backend) < 0 ? OLOG_EXPLODED : 0);
            break;
        default:
            ret = OLOG_EXPLODED;
        }
    }

    return ret < 0 ? OLOG_EXPLODED : 0;
}
static int olog_backend_streamf_destroy(olog *olog, olog_backend *backend) {
    olog_streamf *data = (olog_streamf*) backend->data;
    if (data->filename && data->stream) {
        fclose(data->stream);
    }
    if (data->commands) {
        free(data->commands);
    }
    if (data->format) {
        free(data->format);
    }
    if (data->strftime_format) {
        free(data->strftime_format);
    }
    free(data);
    return 0;
}
static int olog_backend_streamf_func_printf(olog *olog, olog_backend *backend, uint8_t severity, va_list args) {
    const char* format = va_arg(args, const char*);
    int ret = vfprintf((FILE*) ((olog_streamf*) backend->data)->stream, format, args);
    return ret < 0 ? OLOG_EXPLODED : 0;
}
static int olog_backend_streamf_func_sev_lower(olog *olog, olog_backend *backend, uint8_t severity) {
    int ret = fprintf((FILE*) ((olog_streamf*) backend->data)->stream, "%s", severity_l[severity]);
    return ret == EOF ? OLOG_EXPLODED : 0;
}
static int olog_backend_streamf_func_sev_upper(olog *olog, olog_backend *backend, uint8_t severity) {
    int ret = fprintf((FILE*) ((olog_streamf*) backend->data)->stream, "%s", severity_u[severity]);
    return ret == EOF ? OLOG_EXPLODED : 0;
}
static int olog_backend_streamf_func_sev_mixed(olog *olog, olog_backend *backend, uint8_t severity) {
    int ret = fprintf((FILE*) ((olog_streamf*) backend->data)->stream, "%s", severity_m[severity]);
    return ret == EOF ? OLOG_EXPLODED : 0;
}
static int olog_backend_streamf_func_sev_numeric(olog *olog, olog_backend *backend, uint8_t severity) {
    int ret = fprintf((FILE*) ((olog_streamf*) backend->data)->stream, "%u", severity);
    return ret == EOF ? OLOG_EXPLODED : 0;
}

static int olog_backend_streamf_func_time(olog *olog, olog_backend *backend) {
    olog_streamf *data = (olog_streamf*) backend->data;

    time_t t_now = time(0);
    struct tm *tm_now = localtime(&t_now);
    char tstring[data->strftime_size];
    strftime(tstring, data->strftime_size, data->strftime_format, tm_now);

    int ret = fprintf(data->stream, "%s", tstring);
    return ret == EOF ? OLOG_EXPLODED : 0;
}

static olog_streamf *olog_backend_streamf_parse_format(olog_streamf_arguments *arguments) {
    olog_streamf *ret = malloc(sizeof(olog_streamf));
    if (ret == NULL) {
        return NULL;
    }

    ret->strftime_format = NULL;
    if (!arguments->strftime_size) {
        ret->strftime_size = 64;
    } else {
        ret->strftime_size = arguments->strftime_size;
    }

    size_t format_length = strlen(arguments->format);
    if (format_length == 0) {
        fprintf(stderr, "Format unspecified, invalid configuration");
        return NULL;
    }
    ++format_length;
    ret->format = malloc(sizeof(char) * format_length);
    if (ret->format == NULL) {
        return NULL;
    }
    strncpy(ret->format, arguments->format, format_length);

    char *pt = ret->format;
    unsigned max_command_count = 1;

    while (*pt != '\0') {
        if (*pt == '%') {
            ++pt;
            if (*pt == '\0') {
                fprintf(stderr, "failed to compile log format: unfinished escape sequence at end of string");
                return NULL;
            } else if (*pt != '%' && *pt != 's') {
                ++max_command_count;
            }
        }
        ++pt;
    }

    print_command *commands = malloc(sizeof(print_command) * max_command_count * 2);
    if (commands == NULL) {
        return NULL;
    }
    unsigned command_i = 0;

    bool added_user_message = false;
    bool last_char_was_command = true;

    char *first_char;
    first_char = pt = ret->format;
    while (*pt != '\0') {
//		printf("%s\n", pt);
        if (*pt == '%') {
            if (!last_char_was_command) {
                commands[command_i].type = STREAMF_TYPE_STRING;
                commands[command_i++].data = first_char;

                last_char_was_command = true;

                if (*(pt + 1) == '%') {
                    *(pt + 1) = '\0';
                    pt += 2;
                    continue;
                } else {
                    *pt = '\0';
                }
            }

            ++pt;

            if (*pt == '\0') {
                fprintf(stderr, "failed to compile log format: unfinished escape sequence at end of string");
                return NULL;
//			} else if (*pt == '%') {
//				/* skip */
            } else if (*pt == '{' || *pt == '[') {
                char start_char = *pt;
                char *quote_end = NULL;
                if (start_char == '{') {
                    quote_end = strchr(++pt, '}');
                } else {
                    quote_end = strchr(++pt, ']');
                }
                if (quote_end == NULL) {
                    fprintf(stderr, "failed to compile log format: unterminated named pointer name starting at: \"%16s\"", pt);
                    return NULL;
                }
                *quote_end = '\0';
                quote_end++;
                switch (*quote_end) {
                case 's':
                    commands[command_i].type = STREAMF_TYPE_STRING;
                    break;
                case 'f':
                    commands[command_i].type = STREAMF_TYPE_FUNCTION;
                    break;
                case 'F':
                    commands[command_i].type = STREAMF_TYPE_FUNCTION_NO_VARGS;
                    break;
                case 'A':
                    commands[command_i].type = STREAMF_TYPE_FUNCTION_NO_ARGS;
                    break;
                default:
                    fprintf(stderr, "failed to compile log format: named pointer name is not followed by either `f` or `s``");
                    return NULL;
                }
                if (start_char == '{') {
                    if (arguments->persistent_hash) {
                        commands[command_i].type |= STREAMF_TYPE_LOOKUP_HASH;
                        commands[command_i++].data = pt;
                    } else {
                        commands[command_i++].data = phash_get(arguments->named_pointers, pt);
                    }
                } else {
                    char *endptr;
                    errno = 0;
                    unsigned long long index = strtoull(pt, &endptr, 10);
                    if (!errno && *endptr != '\0') {
                        errno = 1;
                    }
                    if (errno) {
                        fprintf(stderr, "not a number");
                        return NULL;
                    } else {
                        if (arguments->persistent_indexes) {
                            commands[command_i].type |= STREAMF_TYPE_LOOKUP_ARRAY;
                            commands[command_i++].data = (void*) index;
                        } else {
                            commands[command_i++].data = arguments->indexed_pointers[index];
                        }
                    }
                }
                pt = quote_end + 1;
            } else {
                switch (*pt++) {
                case 's': {
                    added_user_message = true;
                    commands[command_i].type = STREAMF_TYPE_FUNCTION;
                    commands[command_i++].data = olog_backend_streamf_func_printf;
                }
                    break;
                case 'l':
                    commands[command_i].type = STREAMF_TYPE_FUNCTION_NO_VARGS;
                    commands[command_i++].data = olog_backend_streamf_func_sev_lower;
                    break;
                case 'u':
                    commands[command_i].type = STREAMF_TYPE_FUNCTION_NO_VARGS;
                    commands[command_i++].data = olog_backend_streamf_func_sev_upper;
                    break;
                case 'm':
                    commands[command_i].type = STREAMF_TYPE_FUNCTION_NO_VARGS;
                    commands[command_i++].data = olog_backend_streamf_func_sev_mixed;
                    break;
                case 'n':
                    commands[command_i].type = STREAMF_TYPE_FUNCTION_NO_VARGS;
                    commands[command_i++].data = olog_backend_streamf_func_sev_numeric;
                    break;
                case 't':
                    if (ret->strftime_format == NULL) {
                        if (arguments->strftime_format == NULL) {
                            arguments->strftime_format = "%FT%T%z"; //2017-01-28T22:44:31+00:00
                        }
                        size_t strftime_format_length = strlen(arguments->strftime_format);
                        ret->strftime_format = malloc(sizeof(char) * strftime_format_length + 1);
                        if (ret->strftime_format == NULL) {
                            return NULL;
                        }
                        strncpy(ret->strftime_format, arguments->strftime_format, strftime_format_length);
                        ret->strftime_format[strftime_format_length] = '\0';
                    }
                    commands[command_i].type = STREAMF_TYPE_FUNCTION_NO_ARGS;
                    commands[command_i++].data = olog_backend_streamf_func_time;
                    break;
                default:
                    fprintf(stderr, "failed to compile log format: unknown escape sequence: %%%c", pt[-1]);
                    return NULL;
                }
                continue;
            }
        } else if (last_char_was_command) {
            first_char = pt++;
            last_char_was_command = false;
        } else {
            ++pt;
        }
    }
    if (!last_char_was_command) {
        commands[command_i].type = STREAMF_TYPE_STRING;
        commands[command_i++].data = first_char;
    }
    if (!added_user_message) {
        commands[command_i].type = STREAMF_TYPE_FUNCTION;
        commands[command_i++].data = olog_backend_streamf_func_printf;
    }
    ret->commands = commands;
    ret->command_count = command_i;

    if (arguments->named_pointers && arguments->persistent_hash) {
        ret->named_pointers = arguments->named_pointers;
    }
    if (arguments->indexed_pointers && arguments->persistent_indexes) {
        ret->indexed_pointers = arguments->indexed_pointers;
    }

    if (arguments->filename) {
        ret->filename = arguments->filename;
        ret->stream = fopen(arguments->filename, arguments->file_mode);
        if (ret->stream == NULL) {
            fprintf(stderr, "Error opening file %s: %s\n", arguments->filename, strerror(errno));
            return NULL;
        }
    } else if (arguments->stream) {
        ret->filename = NULL;
        ret->stream = arguments->stream;
    } else {
        fprintf(stderr, "No stream or filename given\n");
        return NULL;
    }

    return ret;
}

olog_backend olog_backend_streamf(olog_streamf_arguments *arguments, int8_t min_log_level, int8_t max_log_level) {
    olog_backend ret = {
            .min_log_level = min_log_level,
            .max_log_level = max_log_level,
            .data = (void*) olog_backend_streamf_parse_format(arguments),
            .print_callback = olog_backend_streamf_printf,
            .destroy_callback = olog_backend_streamf_destroy
    };
    return ret;
}

