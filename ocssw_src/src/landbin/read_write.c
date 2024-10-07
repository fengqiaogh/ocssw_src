#include <stdio.h>
#include <fcntl.h>
#include <unistd.h>
#include "mfhdf.h"

int32 open_sds(int32 sd_id, char *sds_name) {
    int32 sds_id, sds_index;
    if ((sds_index = SDnametoindex(sd_id, sds_name)) == -1) {
        fprintf(stderr, "Cannot locate SDS %s\n", sds_name);
        return -1;
    }
    if ((sds_id = SDselect(sd_id, sds_index)) == -1) {
        fprintf(stderr, "Cannot select SDS %s\n", sds_name);
        return -1;
    }
    return sds_id;
}

int read_sds_block(int32 sds_id, int startline, int buflines, void *dataP) {
    int32 start[3] = {0, 0, 0}, dimsize[3];
    int32 numtype, nattrs, rank;
    char sds_name[H4_MAX_NC_NAME];

    if (SDgetinfo(sds_id, sds_name, &rank, dimsize, &numtype, &nattrs) == -1) {
        fprintf(stderr, "Cannot get info for SDS %s\n", sds_name);
        SDendaccess(sds_id);
        return -1;
    }
    start[0] = startline;
    if (buflines > 0) {
        dimsize[0] = buflines;
        if (SDreaddata(sds_id, start, NULL, dimsize, dataP) == -1) {
            fprintf(stderr, "Cannot read line %d to %d from SDS %s\n", startline, startline + buflines - 1, sds_name);
            SDendaccess(sds_id);
            return -1;
        }
    }
    SDendaccess(sds_id);
    return 0;
}

int read_sds(int32 sd_id, char *sds_name, void *dataP) {
    int32 sds_id, start[6], rank, dimsize[6];
    int32 numtype, nattrs, i;
    char name[128];

    if ((sds_id = open_sds(sd_id, sds_name)) == -1) return -1;
    if (SDgetinfo(sds_id, name, &rank, dimsize, &numtype, &nattrs) < 0) {
        fprintf(stderr, "Cannot get info for SDS %s\n", sds_name);
        return -1;
    }
    for (i = 0; i < rank; i++) start[i] = 0;
    if (SDreaddata(sds_id, start, NULL, dimsize, dataP) < 0) {
        fprintf(stderr, "Cannot read data from SDS %s\n", sds_name);
        return -1;
    }
    SDendaccess(sds_id);
    return 0;
}

int read_attr(int32 sd_id, char *attrname, void *attr) {
    int32 hdf_ind;

    if ((hdf_ind = SDfindattr(sd_id, attrname)) < 0) {
        fprintf(stderr, "Cannot find global attribute %s\n", attrname);
        return -1;
    }
    if (SDreadattr(sd_id, hdf_ind, attr) < 0) {
        fprintf(stderr, "Cannot read global attribute %s\n", attrname);
        return -1;
    }
    return 0;
}

int write_sds(int32 sd_id, char *SDSname, void *data, int startline, int buflines, int Nl, int Np, int32 numtype, void *fillvalue, float64 scale_factor, char verbose, char gzip) {
    int32 sds_id, sds_idx, rank, dimsizes[2], edges[2], start[2];
    comp_info c_info;

    rank = 2;
    start[0] = startline;
    start[1] = 0;
    edges[0] = buflines;
    edges[1] = Np;
    dimsizes[0] = Nl;
    dimsizes[1] = Np;
    if ((sds_idx = SDnametoindex(sd_id, SDSname)) == -1) {
        if ((sds_id = SDcreate(sd_id, SDSname, numtype, rank, dimsizes)) == -1) {
            fprintf(stderr, "Cannot create SDS %s\n", SDSname);
            return -1;
        }
        if (verbose) printf("  Writing SDS %s ...\n", SDSname);
    } else {
        if ((sds_id = SDselect(sd_id, sds_idx)) == -1) {
            fprintf(stderr, "Cannot select SDS %s\n", SDSname);
            return -1;
        }
        if (verbose) printf("  Overwriting SDS %s ...\n", SDSname);
    }
    if (fillvalue)
        SDsetfillvalue(sds_id, fillvalue);
    if (scale_factor != 0)
        SDsetcal(sds_id, scale_factor, scale_factor / 2., 0, 0, DFNT_FLOAT32);
    if (gzip) {
        c_info.deflate.level = 6;
        SDsetcompress(sds_id, COMP_CODE_DEFLATE, &c_info);
    }
    if (startline >= 0) {
        /*
            step = 1000000 / Np;
        printf("step %d\n", step);
            edges[0] = step;
            sizeoftype = DFKNTsize(numtype);
            for (k=0; k<buflines; k+=step) {
              start[0] = startline + k;
              if (k + step > buflines) edges[0] = buflines - k;
              if (SDwritedata(sds_id, start, NULL, edges, (VOIDP)((char *)data + k * Np * sizeoftype)) == -1) {
                fprintf(stderr,"Cannot write line %d to %d in SDS %s\n",start[0], start[0] + edges[0] - 1, SDSname);
                return -1;
              }
            }
         */
        if (SDwritedata(sds_id, start, NULL, edges, (VOIDP) data) == -1) {
            fprintf(stderr, "Cannot write line %d to %d in SDS %s\n", start[0], start[0] + edges[0] - 1, SDSname);
            return -1;
        }
#ifndef HDF41r3
        if (sds_idx == -1 &&
                startline + buflines < Nl) {
            {
                int k;
                switch (numtype) {
                case DFNT_UINT8: for (k = 0; k < Np; k++) ((uint8 *) data)[k] = *((uint8 *) fillvalue);
                    break;
                case DFNT_INT16: for (k = 0; k < Np; k++) ((int16 *) data)[k] = *((int16 *) fillvalue);
                    break;
                default: fprintf(stderr, "Cannot initialize data type %d\n", numtype);
                    return -1;
                }
                printf("  Initializing line %d to %d in SDS %s ...\n", startline + buflines, Nl - 1, SDSname);
                edges[0] = 1;
                for (k = startline + buflines; k < Nl; k++) {
                    start[0] = k;
                    if (SDwritedata(sds_id, start, NULL, edges, (VOIDP) data) == -1) {
                        fprintf(stderr, "Cannot write line %d in SDS %s\n", k, SDSname);
                        return -1;
                    }
                }
            }
        }
#endif
    }
    SDendaccess(sds_id);
    return 0;

}
