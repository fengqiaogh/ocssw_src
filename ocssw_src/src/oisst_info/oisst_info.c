#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>
#include "genutils.h"


/* ----------------------------------------------------------------------------------- */
/* oisst_info() - echo meta-data from flat binary Reynolds OIV2 SST                    */
/*                                                                                     */
/* The files were written in IEEE binary (big-endian). Each file contains four FORTRAN */
/* records described as follows:                                                       */
/*                                                                                     */
/*    rec 1: date and version number        (8 4-byte integer words)                   */
/*    rec 2: gridded sst values in degC     (360*180 4-byte real words)                */
/*    rec 3: normalized error variance      (360*180 4-byte real words)                */
/*    rec 4: gridded ice concentration      (360*180 1-byte integer words)             */
/*                                                                                     */
/* B. Franz, SAIC, May 2004.                                                           */

/* ----------------------------------------------------------------------------------- */

int main(int argc, char *argv[]) {
    char *sstfile = argv[1];
    FILE *fp = NULL;
    int32_t syear, smon, sday;
    int32_t eyear, emon, eday;
    int32_t ndays, version;

    if ((fp = fopen(sstfile, "r")) == NULL) {
        printf("Error opening SST reference file %s for reading.\n", sstfile);
        exit(1);
    }

    if (fseek(fp, 4, SEEK_SET) < 0) {
        printf("Error reading SST reference file %s.\n", sstfile);
        exit(1);
    }
    if (fread(&syear, sizeof (int32_t), 1, fp) != 1) {
        printf("Error reading SST reference file %s.\n", sstfile);
        exit(1);
    }
    fread(&smon, sizeof (int32_t), 1, fp);
    fread(&sday, sizeof (int32_t), 1, fp);
    fread(&eyear, sizeof (int32_t), 1, fp);
    fread(&emon, sizeof (int32_t), 1, fp);
    fread(&eday, sizeof (int32_t), 1, fp);
    fread(&ndays, sizeof (int32_t), 1, fp);
    fread(&version, sizeof (int32_t), 1, fp);
    fseek(fp, 4, SEEK_CUR);

    if (endianess() == 1) {
        swapc_bytes((char *) &syear, sizeof (int32_t), 1);
        swapc_bytes((char *) &smon, sizeof (int32_t), 1);
        swapc_bytes((char *) &sday, sizeof (int32_t), 1);
        swapc_bytes((char *) &eyear, sizeof (int32_t), 1);
        swapc_bytes((char *) &emon, sizeof (int32_t), 1);
        swapc_bytes((char *) &eday, sizeof (int32_t), 1);
        swapc_bytes((char *) &ndays, sizeof (int32_t), 1);
        swapc_bytes((char *) &version, sizeof (int32_t), 1);
    }

    printf("OI Reynolds SST file %s\n", sstfile);
    printf("start_date=%.4d-%.2d-%.2d 00:00:00\n", syear, smon, sday);
    printf("end_date=%.4d-%.2d-%.2d 23:59:59\n", eyear, emon, eday);

    exit(0);
}

