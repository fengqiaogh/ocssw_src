
#include <stdio.h>
#include <math.h>
#include "hdf.h"
#include "mfhdf.h"
#include "export.h"

void read_landmask(int argc, IDL_VPTR argv[]) {

    int32 i;
    int32 j;
    int32 k;
    int32 sds_id;
    int32 dims[8];
    int32 rank;
    int32 dtype;
    int32 nattrs;
    int32 status;

    static char buffer[120];
    UCHAR *data;

    int32 nfiles;
    int32 scan;
    int32 nscan;
    int32 ncol;
    int32 offset;
    int32 sd_id_r;

    int32 ndatasets;
    int32 nglobal_attr;
    int32 count;
    int32 dim_id_r;

    int32 start[2] = {0, 0};
    int32 edges[2] = {0, 0};

    int16 icol, irow;

    int32 gridrows;
    int32 *numbin;
    int32 *basebin;
    int32 totbins;
    float *latbin;

    IDL_VARIABLE *ret;
    UCHAR dat_typ, n_dims;
    IDL_LONG nelem, dim[8];

    start[0] = (IDL_LONG) argv[0]->value.l;
    edges[0] = (IDL_LONG) argv[1]->value.l;

    sd_id_r = SDstart("/home/simbiosd/data/common/LandWater.ISIN.30arcsec.hdf",
            DFACC_RDONLY);

    sds_id = SDselect(sd_id_r, 0);
    status = SDgetinfo(sds_id, buffer, &rank, dims, &dtype, &nattrs);

    edges[0] = 1;
    edges[1] = dims[1];

    n_dims = 2;
    dim[0] = edges[0];
    dim[1] = dims[1];
    data = (UCHAR *) IDL_MakeTempArray(IDL_TYP_BYTE, n_dims, dim,
            IDL_BARR_INI_ZERO, &ret);

    status = SDreaddata(sds_id, start, NULL, edges, (VOIDP) data);

    SDendaccess(sds_id);
    SDend(sd_id_r);

    IDL_VarCopy(ret, argv[2]);

    return;
}

