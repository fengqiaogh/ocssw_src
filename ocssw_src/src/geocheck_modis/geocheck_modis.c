/*
    Modification history:
    Programmer       Organization      Date      Description of change
    --------------   ------------    --------    ---------------------
    Joel Gales       Futuretech      06/21/05    Original Development


gcc -g -c goodgeo.c -I/home/joel/hdf/include
gcc -g goodgeo.o /home/joel/hdf/lib/libmfhdf.a \
/home/joel/hdf/lib/libdf.a /home/joel/hdf/lib/libjpeg.a \
/home/joel/hdf/lib/libz.a -lm -o goodgeo

scp goodgeo swdev@swdev:/home/swdev/bin

scp goodgeo.c swdev@swdev:/home/swdev/src/goodgeo

 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "hdf.h"
#include "mfhdf.h"

#define VERSION "0.1"

void usage(char *progname) {
    printf("%s %s (%s %s)\n",
            progname, VERSION, __DATE__, __TIME__);

    printf("\nUsage: %s geofile [threshold]\n", progname);
    printf("   geofile    = geolocation filename\n");
    printf("   threshold  = pass threshold in percent (default: 20)\n");
    exit(0);
}

int main(int argc, char *argv[]) {

    int32 i;
    int32 k;
    int32 n;
    int32 sds_id;
    int32 dims[8];
    int32 rank;
    int32 nelem;
    int32 dtype;
    int32 nattrs;
    int32 status = 0;

    static char buffer[2 * 2048];

    float *data;

    int32 sd_id_r;
    int32 sds_index;
    int32 start[3] = {0, 0, 0};
    int32 count;
    double threshold;

    if (argc == 1) {
        usage(argv[0]);
        exit(0);
    }

    if (argc == 2) threshold = 20.0;
    else threshold = atof(argv[2]);

    sd_id_r = SDstart(argv[1], DFACC_RDONLY);
    if (sd_id_r == -1) {
        printf("%s not found.\n", argv[1]);
        return -1;
    }

    sds_index = SDnametoindex(sd_id_r, "Longitude");

    if (sds_index != -1) {
        sds_id = SDselect(sd_id_r, sds_index);
        status = SDgetinfo(sds_id, buffer, &rank, dims, &dtype, &nattrs);

        nelem = dims[0];
        for (k = 1; k < rank; k++) nelem *= dims[k];

        data = (float *) calloc(nelem, DFKNTsize(dtype));

        status = SDreaddata(sds_id, start, NULL, dims, (VOIDP) data);

        count = 0;
        for (i = 0; i < nelem; i++) if (data[i] == -999) count++;
        n = (int32) rintf(100 * ((float32) nelem - count) / nelem);
        printf("Percent Longitude good: %d\n", n);
        if (n < threshold) status = -1;
    } else status = -1;


    sds_index = SDnametoindex(sd_id_r, "Latitude");

    if (sds_index != -1) {
        sds_id = SDselect(sd_id_r, sds_index);
        status = SDgetinfo(sds_id, buffer, &rank, dims, &dtype, &nattrs);

        nelem = dims[0];
        for (k = 1; k < rank; k++) nelem *= dims[k];

        data = (float *) calloc(nelem, DFKNTsize(dtype));

        status = SDreaddata(sds_id, start, NULL, dims, (VOIDP) data);

        count = 0;
        for (i = 0; i < nelem; i++) if (data[i] == -999) count++;
        n = (int32) rintf(100 * ((float32) nelem - count) / nelem);
        printf("Percent Latitude  good: %d\n", n);
        if (n < threshold) status = -1;
    } else status = -1;


    SDend(sd_id_r);

    return status;
}

