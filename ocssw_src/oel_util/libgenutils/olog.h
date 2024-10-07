/** @file olog.h
	@brief A simple logger, capable of dispatching log events to multiple end
		points.

	olog is meant to be as simple as possible to use, configure, and extend.

*/

#ifndef __OLOG_H_
#define __OLOG_H_

#include <stdarg.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
/*
https://en.wikipedia.org/wiki/Syslog#Severity_level

Value   Severity        Keyword   Description/Examples
0       Emergency       emerg     System is unusable 	This level should not be used by applications.
1       Alert           alert     Should be corrected immediately 	Loss of the primary ISP connection.
                                  Ski Haus Delta has not reported status within status_timeout (120)
2       Critical        crit      Critical conditions 	A failure in the system's primary application.
                                  Ski Haus Delta reports temperature < low_critical (30)
3       Error           err       Error conditions 	An application has exceeded its file storage limit and attempts to write are failing.
                                  Ski Haus Delta reports temperature < low_error (32)
4       Warning         warn      May indicate that an error will occur if action is not taken. 	A non-root file system has only 2GB remaining.
                                  Ski Haus Delta reports temperature < low_warning (36)
5       Notice          notice    Events that are unusual, but not error conditions.
                                  Ski Haus Delta reports temperature < low_notice(50)
6       Informational   info      Normal operational messages that require no action. 	An application has started, paused or ended successfully.
                                  Ski Haus Delta reports temperature 60
7       Debug           debug     Information useful to developers for debugging the application.

The meaning of severity levels other than Emergency and Debug are relative to the application.
For example, if the purpose of the system is to process transactions to update customer account
balance information, an error in the final step should be assigned Alert level. However, an error
occurring in an attempt to display the ZIP code of the customer may be assigned Error or even
Warning level.
*/

#define OLOG_DEBUG     0
#define OLOG_INFO      1
#define OLOG_NOTICE    2
#define OLOG_WARNING   3
#define OLOG_ERROR     4
#define OLOG_CRITICAL  5
#define OLOG_ALERT     6
#define OLOG_EMERGENCY 7

#define OLOG_CONTINUE   0
#define OLOG_NOERROR    0
#define OLOG_EXPLODED   1
#define OLOG_STOP     255

typedef struct olog olog;
typedef struct olog_backend olog_backend;

typedef int (*olog_print_callback)(olog *olog, olog_backend *backend, uint8_t severity, va_list args);
typedef int (*olog_destroy_callback)(olog *olog, olog_backend *backend);
typedef int (*olog_has_error_callback)(olog *olog, olog_backend *backend);
typedef int (*olog_print_error_callback)(olog *olog, olog_backend *backend, FILE* stream);

typedef struct olog_backend {
    /** @brief Minimum log level to receive. */
    int8_t min_log_level;

    /** @brief Maximum log level to receive. */
    int8_t max_log_level;

    /** @brief Called to check if backend is valid. */
    olog_has_error_callback has_error_callback;

    /** @brief Called to print an error, if one exists. */
    olog_print_error_callback print_error_callback;

    /** @brief Called when given data to print. */
    olog_print_callback print_callback;

    /** @brief Called when olog object is destroyed. */
    olog_destroy_callback destroy_callback;

    /** @brief Storage for backend use. */
    void *data;
} olog_backend;

extern olog *olog_global_log;
extern char *olog_global_config;
extern olog_backend olog_backends_end;

olog *olog_create(olog_backend *backends, FILE *stream);
void olog_destroy(olog *olog);

bool olog_has_error(olog *olog);
bool olog_print_error(olog *olog, FILE *stream);

int olog_print(olog *olog, int level, ...);
int olog_vprint(olog *olog, int8_t level, va_list args);

int olog_debug(olog *olog, ...);
int olog_info(olog *olog, ...);
int olog_notice(olog *olog, ...);
int olog_warn(olog *olog, ...);
int olog_err(olog *olog, ...);
int olog_crit(olog *olog, ...);
int olog_emerg(olog *olog, ...);

#define olog_print_verbose(olog, level, ...) \
	do { \
		olog_print(olog, level, __VA_ARGS__); \
		olog_print(olog, level, "Previous message at %s:%d\n", __FILE__, __LINE__); \
	} while (0)

#endif /* __OLOG_H_ */
