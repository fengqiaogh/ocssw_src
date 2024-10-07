/*
   create HDF SDS

   Created by Gary Fu, SAIC GSC, 4/98
 */

#include "hdf.h"
#include "mfhdf.h"

int32 create_sds(int32 sdfid, char *sdsname, int32 nt, int32 rank,
        int32 *dims, int32 vid, VOIDP *data) {
    int32 i, sdsid, sdsref;
    int32 start[3] = {0, 0, 0}, edge[3] = {0, 0, 0};

    if ((sdsid = SDcreate(sdfid, sdsname, nt, rank, dims)) < 0)
        return FAIL;

    for (i = 0; i < rank; i++)
        edge[i] = dims[i];

    if ((SDwritedata(sdsid, start, NULL, edge, (VOIDP) data)) < 0)
        return FAIL;

    if ((sdsref = SDidtoref(sdsid)) < 0) {
        fprintf(stderr, "\nSDidtoref failed for %s", sdsname);
        return FAIL;
    }

    if (vid > 0) {
        if ((Vaddtagref(vid, DFTAG_NDG, sdsref)) < 0) {
            fprintf(stderr, "\nVaddtagref failed while trying to link %s", sdsname);
            return FAIL;
        }
    }

    SDendaccess(sdsid);

    return sdsid;
}

