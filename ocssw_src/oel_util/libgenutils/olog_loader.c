#include "olog.h"
#include "olog/loader.h"
#include "olog/stream.h"
#include "olog/streamf.h"

#include <stdlib.h>

olog *olog_global_log = NULL;
char *olog_global_config = NULL;
olog_backend olog_backends_end = { 0 };

olog *olog_load_default() {
    if (!olog_global_log) {
        olog_streamf_arguments streamf_args = {
                .stream = stderr,
                .format = "%t [%u] %s"
        };
        olog_backend backends[] = {
                olog_backend_streamf(&streamf_args, OLOG_ERROR, OLOG_EMERGENCY),
                olog_backend_stream(stdout, OLOG_INFO, OLOG_WARNING),
                { 0 }
        };
        olog_global_log = olog_create(backends, stderr);
    }
    return olog_global_log;
}

void olog_destroy_default() {
    if (olog_global_log) {
        olog_destroy(olog_global_log);
    }
}
