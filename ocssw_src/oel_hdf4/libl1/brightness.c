#include <sys/types.h>
#include <unistd.h>

/* ============================================================================ */
/* module brightness.c - convert observed IR radiance to temperature            */
/*                                                                              */
/* Written By: B. Franz, NASA/SIMBIOS, August 2003.                             */
/*                                                                              */
/* Generalized for other sensors, B. Franz, 12/2010                             */
/* ============================================================================ */

#include "l1.h"
#include <hdf.h>
#include <mfhdf.h>

#define NBTDIMS     3
#define NBTTABMAX 100
#define NBTDETMAX  20
#define NBTBANDMAX 10

static int16 ntab;
static int16 ndet;
static int16 nbnd;

static float radtab [NBTBANDMAX][NBTDETMAX][NBTTABMAX];
static float temptab [NBTBANDMAX][NBTDETMAX][NBTTABMAX];
static float radinc [NBTBANDMAX][NBTDETMAX];
static int16 uselog [NBTBANDMAX];

/* ----------------------------------------------------------------------------------- */
/* read_bt_table() - reads the radiance to temperature HDF file                        */
/*                                                                                     */
/* B. Franz, SAIC, August 2003.                                                        */

/* ----------------------------------------------------------------------------------- */
int read_bt_table(char *filename, int nbands, int ndets) {
    int32 sd_id;
    int32 sds_id;
    int32 status;
    int32 rank;
    int32 nt;
    int32 nattrs;
    int32 dims [H4_MAX_VAR_DIMS];
    int32 start [NBTDIMS] = {0, 0, 0};
    int32 end [NBTDIMS] = {1, 1, 1};
    char name [H4_MAX_NC_NAME] = "";
    char sdsname[H4_MAX_NC_NAME] = "";

    int ibnd, idet, itab;


    if (strcmp(filename, "") == 0) {
        printf("\nNo brightness temperature conversion provided for this sensor.\n");
        return (1);
    }

    if (want_verbose)
        printf("Loading radiance to brightness temperature from %s\n", filename);

    // Open the file

    sd_id = SDstart(filename, DFACC_RDONLY);
    if (sd_id == FAIL) {
        fprintf(stderr, "-E- %s line %d: SDstart(%s, %d) failed.\n",
                __FILE__, __LINE__, filename, DFACC_RDONLY);
        return (1);
    }

    // Get dimension attributes

    status = SDreadattr(sd_id, SDfindattr(sd_id, "Number of Bands"), &nbnd);
    if (status != 0) {
        status = SDreadattr(sd_id, SDfindattr(sd_id, "number_of_bands"), &nbnd);
        if (status != 0) {
            printf("-E- %s Line %d:  Error reading global attribute %s from %s.\n",
                    __FILE__, __LINE__, "Number of Bands", filename);
            exit(1);
        }
    }
    status = SDreadattr(sd_id, SDfindattr(sd_id, "Number of Detectors"), &ndet);
    if (status != 0) {
        status = SDreadattr(sd_id, SDfindattr(sd_id, "number_of_detectors"), &ndet);
        if (status != 0) {

            printf("-E- %s Line %d:  Error reading global attribute %s from %s.\n",
                    __FILE__, __LINE__, "Number of Detectors", filename);
            exit(1);
        }
    }
    status = SDreadattr(sd_id, SDfindattr(sd_id, "Number of Radiance Levels"), &ntab);
    if (status != 0) {
        status = SDreadattr(sd_id, SDfindattr(sd_id, "number_of_radiance_levels"), &ntab);
        if (status != 0) {

            printf("-E- %s Line %d:  Error reading global attribute %s from %s.\n",
                    __FILE__, __LINE__, "Number of Radiance Levels", filename);
            exit(1);
        }
    }
    if (ntab > NBTTABMAX || nbnd > NBTBANDMAX || ndet > NBTDETMAX) {
        printf("-E- %s Line %d:  BT table dimensions exceed max allowable (%d %d %d).\n",
                __FILE__, __LINE__, ntab, nbnd, ndet);
        exit(1);
    }
    if (nbnd != nbands || ndet > ndets) {
        printf("-E- %s Line %d:  BT table dimensions does not match sensor attributes (%d %d).\n",
                __FILE__, __LINE__, nbnd, ndet);
        exit(1);
    }

    // Read log flag

    strcpy(sdsname, "Uselog");
    sds_id = SDselect(sd_id, SDnametoindex(sd_id, sdsname));
    status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
    status = SDreaddata(sds_id, start, NULL, dims, (VOIDP) uselog);
    if (status != 0) {
        printf("-E- %s:  Error reading SDS %s from %s.\n", __FILE__, sdsname, filename);
        exit(1);
    }
    status = SDendaccess(sds_id);

    // Read radiance levels

    strcpy(sdsname, "Radiances");
    sds_id = SDselect(sd_id, SDnametoindex(sd_id, sdsname));
    status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
    for (ibnd = 0; ibnd < nbnd; ibnd++) for (idet = 0; idet < ndet; idet++) for (itab = 0; itab < ntab; itab++) {
                start[0] = ibnd;
                start[1] = idet;
                start[2] = itab;
                status = SDreaddata(sds_id, start, NULL, end, (VOIDP) & radtab[ibnd][idet][itab]);
                if (status != 0) {
                    printf("-E- %s:  Error reading %s from %s.\n", __FILE__, sdsname, filename);
                    exit(1);
                }
            }
    status = SDendaccess(sds_id);

    // Read temperature levels

    strcpy(sdsname, "Temperatures");
    sds_id = SDselect(sd_id, SDnametoindex(sd_id, sdsname));
    status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
    for (ibnd = 0; ibnd < nbnd; ibnd++) for (idet = 0; idet < ndet; idet++) for (itab = 0; itab < ntab; itab++) {
                start[0] = ibnd;
                start[1] = idet;
                start[2] = itab;
                status = SDreaddata(sds_id, start, NULL, end, (VOIDP) & temptab[ibnd][idet][itab]);
                if (status != 0) {
                    printf("-E- %s:  Error reading %s from %s.\n", __FILE__, sdsname, filename);
                    exit(1);
                }
            }
    status = SDendaccess(sds_id);

    // Close file

    SDend(sd_id);

    return (0);
}

/* ----------------------------------------------------------------------------------- */
/* radiance2bt() - interpolates brightness temps from table                .           */
/*                                                                                     */
/* B. Franz, SAIC, August 2003.                                                        */

/* ----------------------------------------------------------------------------------- */
void radiance2bt(l1str *l1rec, int resolution) {
    static int firstCall = 1;

    int32_t nbands = l1rec->l1file->nbandsir;
    int32_t ndets = l1rec->l1file->ndets;
    int32_t id, ip, ib, ipb;
    int32_t i;
    float a;
    float rad;

    // on first call, load the table
    if (firstCall) {
        firstCall = 0;
        if (l1_input->btfile[0]) {
            if (read_bt_table(l1_input->btfile, nbands, ndets) != 0) {
                printf("Error loading brightness temperature table.\n");
                exit(1);
            }
        } else {
            printf("Brightness temperature file not specified.");
            exit(1);
        }

        for (ib = 0; ib < nbands; ib++)
            for (id = 0; id < ndets; id++)
                radinc[ib][id] = radtab[ib][id][1] - radtab[ib][id][0];
    }

    // modis-specific handling for dets per resolution: needs generalization
    switch (resolution) {
    case 250: id = l1rec->detnum / 4;
        break;
    case 500: id = l1rec->detnum / 2;
        break;
    default: id = l1rec->detnum / 1;
        break;
    }

    for (ip = 0; ip < l1rec->npix; ip++)
        for (ib = 0; ib < nbands; ib++) {

            ipb = ip * NBANDSIR + ib;

            if (!isfinite(l1rec->Ltir[ipb]) || l1rec->Ltir[ipb] <= 0.0) {
                l1rec->Bt[ipb] = BT_LO;
                continue;
            }

            rad = l1rec->Ltir[ipb]*10.0; /* radiance in W/m^2/nm/sr */

            if (uselog[ib])
                rad = log10(rad);

            if (rad <= radtab[ib][id][0])
                l1rec->Bt[ipb] = BT_LO;
            else if (rad >= radtab[ib][id][ntab - 1])
                l1rec->Bt[ipb] = BT_HI;
            else {
                i = (int32_t) ((rad - radtab[ib][id][0]) / radinc[ib][id]);
                a = (rad - radtab[ib][id][i]) / (radtab[ib][id][i + 1] - radtab[ib][id][i]);
                l1rec->Bt[ipb] = temptab[ib][id][i] + a * (temptab[ib][id][i + 1] - temptab[ib][id][i]);
                l1rec->Bt[ipb] -= 273.15;
            }
        }
}
