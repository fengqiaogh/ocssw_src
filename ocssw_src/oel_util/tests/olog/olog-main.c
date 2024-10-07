#include <olog.h>
#include <olog/loader.h>

#include <olog/buffer.h>
#include <olog/file.h>
#include <olog/stream.h>
#include <olog/streamf.h>

#include <phash.h>

#include <stdio.h>
#include <string.h>

static unsigned message_count = 0;
int get_message_count(olog *olog, olog_backend *backend) {
    olog_streamf *backend_data = (olog_streamf*) backend->data;
    int ret = fprintf(backend_data->stream, "%u", ++message_count);
    return ret < 0 ? OLOG_EXPLODED : 0;
}
int get_message_count2(olog *olog, olog_backend *backend) {
    olog_streamf *backend_data = (olog_streamf*) backend->data;
    int ret = fprintf(backend_data->stream, "%3u", ++message_count);
    return ret < 0 ? OLOG_EXPLODED : 0;
}

int main(int argc, char *argv[]) {
    if (0) {
        phash *p1 = phash_create(0);
        phash_set(p1, "host", "host1");

        void *indexed_pointers[1];
        indexed_pointers[0] = get_message_count;

        olog_streamf_arguments streamf_args = {
                .stream = stdout,
                .named_pointers = p1,
                .persistent_hash = true,
                .indexed_pointers = indexed_pointers,
				//.persistent_indexes = true,
                .format = "[%[0]A] %t [%u] %%%n%n | [%l] | [%{host}s] | %s"
        };
        olog_backend streamf = olog_backend_streamf(&streamf_args, 2, 7);
        olog_backend backends[] = {
                streamf,
                { 0 }
        };
        olog *logger = olog_create(backends, stderr);
        if (logger == NULL) {
            return 1;
        }

        olog_err(logger, "Test string 1\n");
        phash_set(p1, "host", "host2");
        indexed_pointers[0] = get_message_count2;
        olog_err(logger, "Test string 2\n");

        phash_destroy(p1);
        olog_destroy(logger);
    }

    if (0) {
        olog_streamf_arguments streamf_args = {
                .stream = stderr,
                .format = "%t [%u] %s"
        };
        olog_backend backends[] = {
                olog_backend_streamf(&streamf_args, 1, 3),
                olog_backend_stream(stdout, 4, 7),
                { 0 }
        };
        olog *logger = olog_create(backends, stderr);
        if (logger == NULL) {
            return 1;
        }

        olog_err(logger, "Test string 1\n");
        olog_err(logger, "Test string 2\n");
        olog_destroy(logger);
    }

    if (0) {
        olog_backend backends[] = {
                olog_backend_stream(stderr, 3, 7),
                olog_backend_stream(stdout, 2, 3),
                //olog_backend_file("/accounts/jlefler/test.out", 0, 7),
                { 0 }
        };
        olog *logger = olog_create(backends, stderr);
        if (logger == NULL) {
            return 1;
        }

        olog_info(logger, "Test string: %s\n", "qwer");
        olog_err(logger, "Test string 3: %s\n", "asdf");
        olog_destroy(logger);
    }

    if (0) {
        olog_backend buffer = olog_backend_buffer(32, 0, 7);
        olog_backend backends[] = {
                buffer,
                { 0 }
        };
        olog *logger = olog_create(backends, stderr);
        if (logger == NULL) {
            return 1;
        }

        for (int i = 0; i < 150; i++) {
//			if (olog_info(logger, "Test string: %d\n", i)){
//				break;
//			}
            if (olog_info(logger, "%d", i % 10)) {
                break;
            }
        }

        printf("%s\n", ((olog_buffer*) buffer.data)->buffer_start);
        olog_destroy(logger);
    }

    if (0) {
        olog_info(NULL, "Test string 1\n");
        olog_crit(NULL, "Test string 2\n");
        olog_print_verbose(NULL, OLOG_INFO, "Test string %d\n", 15);
        olog_crit(NULL, "%s, %d\n", __FILE__, __LINE__);
        olog_destroy(NULL);
//		olog_destroy_default();
    }

    if (1) {
        olog *log = olog_load_default();
        if (olog_global_log == log) {
            printf("Same\n");
        } else {
            printf("not same\n");
        }
    }

    return 0;
}

