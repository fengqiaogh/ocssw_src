#include <setupflags.h>
#include <string.h>
#include <stdlib.h>
//#include <netcdf.h>
//#include <dfutils.h>

void setupflags(char *flagdef, char *flaguse, uint32_t *flagusemask, uint32_t *required, int *status,
                int32_t *BITS) {
    int bitNum;
    char *tmpFlags;
    char *ptr, *ptr2;

    *status = 0;
    *flagusemask = 0;
    *required = 0;

    bitNum = 0;
    tmpFlags = strdup(flagdef);
    ptr = strtok(tmpFlags, ",");
    while (ptr) {
        if (ptr) {
            if ((ptr2 = strstr(flaguse, ptr))) {
                ptr2--;
                if (*ptr2 == '~')
                    *required = *required | BITS[bitNum];
                else
                    *flagusemask = *flagusemask | BITS[bitNum];
            }
        }
        ptr = strtok(NULL, ",");
        bitNum++;
        if (bitNum > 33) {
            *status = -1;
            break;
        }
    }

    free(tmpFlags);
}
