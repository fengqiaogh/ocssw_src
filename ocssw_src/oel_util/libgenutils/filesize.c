#include <sys/types.h>
#include <sys/stat.h>

/* -------------------------------------------------------------- */
/* filesize() - returns file length in bytes                      */

/* -------------------------------------------------------------- */
int32_t filesize(const char *filename) {
    struct stat buf;
    int32_t ret;

    ret = stat(filename, &buf);
    if (ret < 0)
        return (ret);
    else
        return (buf.st_size);
}

int32_t filesize_(const char *filename) {
    return (filesize(filename));
}
