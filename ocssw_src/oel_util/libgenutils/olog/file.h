
#ifndef __OLOG_FILE_H_
#define __OLOG_FILE_H_

#include "olog.h"

olog_backend olog_backend_file_mode(const char *filename, const char *mode, int8_t min_log_level, int8_t max_log_level);
olog_backend olog_backend_file(const char *filename, int8_t min_log_level, int8_t max_log_level);

#endif /* __OLOG_FILE_H_ */
