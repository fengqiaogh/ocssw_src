#include <stdarg.h>
#include "olog.h"
#include "olog/loader.h"

#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct olog {
    olog_backend *backends;
};

olog *olog_create(olog_backend *backends, FILE *stream) {
    if (backends == NULL) {
        if (stream) {
            fprintf(stream, "backends cannot be NULL\n");
        }
        return NULL;
    }
    olog *o = malloc(sizeof(olog));
    if (o == NULL) {
        if (stream) {
            fprintf(stream, "Couldn't malloc olog\n");
        }
        return NULL;
    }

    for (int i = 0;; i++) {
        if (!(backends[i].print_callback || backends[i].destroy_callback)) {
            size_t backends_size = sizeof(olog_backend) * (i + 1);
            o->backends = malloc(backends_size);
            if (o->backends == NULL) {
                free(o);
                if (stream) {
                    fprintf(stream, "Couldn't malloc backends\n");
                }
                return NULL;
            }
            memcpy(o->backends, backends, backends_size);
            break;
        }
    }

    if (stream) {
        if (olog_print_error(o, stream)) {
            olog_destroy(o);
            return NULL;
        }
    } else {
        if (olog_has_error(o)) {
            olog_destroy(o);
            return NULL;
        }
    }

    return o;
}
void olog_destroy(olog *olog) {
    if (olog == NULL) {
        olog = olog_load_default();
    }
    for (int i = 0; olog->backends[i].print_callback || olog->backends[i].destroy_callback; i++) {
        if (olog->backends[i].destroy_callback) {
            olog->backends[i].destroy_callback(olog, &olog->backends[i]);
        }
    }
    free(olog->backends);
    free(olog);
}

bool olog_has_error(olog *olog) {
    if (olog == NULL) {
        olog = olog_load_default();
    }
    for (int i = 0; olog->backends[i].print_callback || olog->backends[i].destroy_callback; i++) {
        if (olog->backends[i].has_error_callback) {
            if (olog->backends[i].has_error_callback(olog, &olog->backends[i])) {
                return true;
            }
        }
    }
    return false;
}
bool olog_print_error(olog *olog, FILE *stream) {
    if (olog == NULL) {
        olog = olog_load_default();
    }
    bool ret = false;
    for (int i = 0; olog->backends[i].print_callback || olog->backends[i].destroy_callback; i++) {
        if (olog->backends[i].print_error_callback) {
            if (olog->backends[i].print_error_callback(olog, &olog->backends[i], stream)) {
                ret = true;
            }
        }
    }
    return ret;
}
int olog_vprint(olog *olog, int8_t level, va_list args) {
    if (olog == NULL) {
        olog = olog_load_default();
    }
    for (int i = 0; olog->backends[i].print_callback || olog->backends[i].destroy_callback; i++) {
        if (olog->backends[i].print_callback && olog->backends[i].min_log_level <= level && olog->backends[i].max_log_level >= level) {
            va_list args_tmp;
            va_copy(args_tmp, args);
            int print_ret = olog->backends[i].print_callback(olog, &olog->backends[i], level, args_tmp);
            if (print_ret != OLOG_CONTINUE) {
                return print_ret;
            }
        }
    }
    return 0;
}

int olog_print(olog *olog, int level, ...) {
    va_list args;
    va_start(args, level);
    int ret = olog_vprint(olog, level, args);
    va_end(args);
    return ret;
}

#define define_olog_print(name, level) \
int olog_ ## name(olog *olog, ...){ \
	va_list args; \
	va_start(args, olog); \
	int ret = olog_vprint(olog, level, args); \
	va_end(args); \
	return ret; \
}

define_olog_print(debug, OLOG_DEBUG)
define_olog_print(info, OLOG_INFO)
define_olog_print(notice, OLOG_NOTICE)
define_olog_print(warn, OLOG_WARNING)
define_olog_print(error, OLOG_ERROR)
define_olog_print(crit, OLOG_CRITICAL)
define_olog_print(alert, OLOG_ALERT)
define_olog_print(emerg, OLOG_EMERGENCY)

