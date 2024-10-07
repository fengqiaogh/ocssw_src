
#ifndef __OLOG_STREAM_H_
#define __OLOG_STREAM_H_

#include "olog.h"

#include <stdio.h>
#include <stdlib.h>

olog_backend olog_backend_stream(FILE *stream, int8_t min_log_level, int8_t max_log_level);

#endif /* __OLOG_STREAM_H_ */
