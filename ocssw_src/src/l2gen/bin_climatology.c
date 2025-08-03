#include "l12_proto.h"

static int32 *basebin;
static int32 *numbin;
static int32 totbins;
static int32 nrows;

#define MAXNVDATA 50

#define MAXCPROD  3
#define CHLORA    0
#define NLW412    1
#define ANGSTROM  2

static int32 fid;
static int32 sdfid;
static int32 vgid;
static int32 vdata_id[MAXNVDATA];

static int32_t *bins;
static float *data [MAXCPROD];
static int32 n_records;

static char pnames[MAXCPROD][20] = {"chlor_a", "nLw_412", "angstrom"};
static int loaded[MAXCPROD] = {0, 0, 0};

/* ---------------------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------------------------- */
int32 open_l3b(char *l3file, int32_t day) {
    static int nmonths = 12;
    static int firstCall = 1;
    static int clim_day [12] = {1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335};
    static char clim_files[12][19] ={"S001031_L3b_MC.hdf", "S032059_L3b_MC.hdf", "S060090_L3b_MC.hdf",
        "S091120_L3b_MC.hdf", "S121151_L3b_MC.hdf", "S152181_L3b_MC.hdf",
        "S182212_L3b_MC.hdf", "S213243_L3b_MC.hdf", "S244273_L3b_MC.hdf",
        "S274304_L3b_MC.hdf", "S305334_L3b_MC.hdf", "S335365_L3b_MC.hdf"};

    char *tmp_str;

    intn i, j;
    int32 vg_ref;
    int32 tag;
    int32 ref;
    int32 vdid;
    char nam_buf[80];
    char cls_buf[80];
    char file[FILENAME_MAX];


    if (l3file != NULL) {

        strcpy(file, l3file);
        if (firstCall) printf("\nOpening L3 data file %s\n", file);
    } else {

        /* Determine filename from day number */
        for (i = 1; i < nmonths; i++)
            if (clim_day[i] > day)
                break;
        i--;

        if ((tmp_str = getenv("OCDATAROOT")) == NULL) {
            printf("OCDATAROOT environment variable is not defined.\n");
            exit(1);
        }
        strcpy(file, tmp_str);
        strcat(file, "/common/");
        strcat(file, clim_files[i]);

        if (firstCall) printf("\nOpening climatology file %s\n", file);
    }

    firstCall = 0;


    /* Open the file and start the SD interface */
    if ((fid = Hopen(file, DFACC_RDONLY, 0)) < 0) {
        printf("-E- %s line %d: Cannot open input HDF file on Hopen\n", __FILE__, __LINE__);
        printf("-E- Level-3 data file - %s\n", file);
        return (0);
    }
    Vstart(fid);
    if ((sdfid = SDstart(file, DFACC_RDONLY)) < 0) {
        printf("-E- %s line %d: Cannot open input HDF file on SDstart- %s\n", __FILE__, __LINE__, file);
        exit(1);
    }

    /* Locate the Binned Data vgroup */

    vgid = -1;
    for (i = 0; i < MAXNVDATA; i++) {
        vdata_id[i] = -1;
    }

    vg_ref = Vfind(fid, "Level-3 Binned Data");
    if (vg_ref < 0) {
        printf("-E- %s line %d: Cannot locate Level-3 Binned Data in %s\n", __FILE__, __LINE__, file);
        exit(1);
    }

    vgid = Vattach(fid, vg_ref, "r");


    /* Locate all the required vdata in the vgroup */

    j = 0;
    for (i = 0; i < Vntagrefs(vgid); i++) {

        Vgettagref(vgid, i, &tag, &ref);
        vdid = VSattach(fid, ref, "r");

        if (vdid == -1) {
            printf("-E- %s line %d: Problem opening Vdata (reference # %d)\n", __FILE__, __LINE__, ref);
            exit(1);
        }

        VSgetname(vdid, nam_buf);
        VSgetclass(vdid, cls_buf);

        if (strcmp(cls_buf, "DataMain") == 0) {
            vdata_id[0] = vdid;
        }

        if (strcmp(cls_buf, "Index") == 0) {
            vdata_id[1] = vdid;
            nrows = VSelts(vdata_id[1]);
        }

        if (strcmp(cls_buf, "DataSubordinate") == 0) {
            vdata_id[2 + j++] = vdid;
        }
    }

    return 1;
}


/* ---------------------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------------------------- */
int32 close_l3b() {
    intn i;

    if (vgid != -1) {

        Vdetach(vgid);

        for (i = 0; i < MAXNVDATA; i++) {
            if (vdata_id[i] != -1) VSdetach(vdata_id[i]);
        }
    }

    SDend(sdfid);
    Vend(fid);
    Hclose(fid);

    return SUCCEED;
}


/* ---------------------------------------------------------------------------------------- */
/* binary search for nearest bin in bin list                                                */

/* ---------------------------------------------------------------------------------------- */
int32 nearest_bin(int32 bin_num) {
    int32 jl = -1, ju = n_records, jm = 0;
    int ascnd;

    ascnd = (bins[n_records - 1] >= bins[0]);
    while (ju - jl > 1) {
        jm = (ju + jl) / 2;
        if (ascnd == (bin_num >= bins[jm]))
            jl = jm;
        else
            ju = jm;
    }
    if (bin_num == bins[jl]) return (jl);
    if (jl + 1 < n_records && bin_num == bins[jl + 1]) return (jl + 1);
    if (jl > 0 && bin_num == bins[jl - 1]) return (jl - 1);
    if (bin_num == bins[0]) return (0);
    if (bin_num == bins[n_records - 1]) return (n_records - 1);
    return (-1);
}

/* ---------------------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------------------------- */
void initbin(void) {
    int row;
    float32 latbin;

    /* ----------------- */
    /* Set up bin arrays */
    /* ----------------- */
    if ((numbin = (int32 *) calloc(nrows, sizeof (int32))) == NULL) {
        printf("-E- %s:%d Error allocating memory to numbin\n", __FILE__, __LINE__);
        exit(FATAL_ERROR);
    }
    if ((basebin = (int32 *) calloc(nrows + 1, sizeof (int32))) == NULL) {
        printf("-E- %s:%d Error allocating memory to basebin\n", __FILE__, __LINE__);
        exit(FATAL_ERROR);
    }
    basebin[0] = 1;
    for (row = 0; row < nrows; row++) {
        latbin = ((row + 0.5)*180.0 / nrows) - 90.0;
        numbin[row] = (int32) (2 * nrows * cos(latbin * OEL_PI / 180.0) + 0.5);
        if (row > 0) {
            basebin[row] = basebin[row - 1] + numbin[row - 1];
        }
    }
    totbins = basebin[nrows - 1] + numbin[nrows - 1] - 1;
}

/* ---------------------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------------------------- */
void closebin(void) {
    free(numbin);
    free(basebin);
}

/* ---------------------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------------------------- */
double constrain_lat(double lat) {
    if (lat > 90) lat = 90;
    if (lat < -90) lat = -90;
    return (lat);
}

/* ---------------------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------------------------- */
double constrain_lon(double lon) {
    while (lon < -180) lon += 360;
    while (lon > 180) lon -= 360;
    return (lon);
}

/* ---------------------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------------------------- */
int16 lat2row(double lat) {
    int16 row;

    row = (int16) ((90 + lat) * nrows / 180.0);
    if (row >= nrows) row = nrows - 1;
    return (row);
}

/* ---------------------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------------------------- */
int32 rowlon2bin(int16 row, double lon) {
    int16 col;
    int32 bin;

    lon = constrain_lon(lon);
    col = (int16) ((lon + 180.0) * numbin[row] / 360.0);
    if (col >= numbin[row]) col = numbin[row] - 1;
    bin = basebin[row] + col;
    return (bin);
}

/* ---------------------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------------------------- */
int32 latlon2bin(double lat, double lon) {
    int16 row, col;
    int32 bin;

    /* Constrain latitudes to [-90,90] and longitudes to [-180,180]. */
    lat = constrain_lat(lat);
    lon = constrain_lon(lon);

    row = lat2row(lat);
    col = (int16) ((lon + 180.0) * numbin[row] / 360.0);
    if (col >= numbin[row]) col = numbin[row] - 1;
    bin = basebin[row] + col;
    return (bin);
}

/* ---------------------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------------------------- */
float bin_climatology(char *l3file, int32_t day, float lon, float lat, char *prodname) {
    static int firstCall = 1;
    static char lastprod[20] = "";
    static int nofile = 0;
    static int32 prodID = -1;
    static int32 interlace0;

    int32 interlace;
    int32 vdata_size;
    int32 prodinx;
    char fields[60];
    char vdata_name[H4_MAX_NC_NAME];
    char buf[200];
    int32 vsize;
    float32 wgt;

    int status;
    int32 nread;
    int32 nprod;
    int32_t i;
    int32 bin_num;
    int32 bin_idx;

    float32 *sum_buf = NULL;
    uint8 *bin_buf = NULL;


    if (strncmp(prodname, lastprod, 8) != 0 || prodID < 0) {
        prodID = -1;
        for (i = 0; i < MAXCPROD; i++)
            if (strncmp(prodname, pnames[i], 8) == 0) {
                prodID = i;
                break;
            }
        if (prodID == -1) {
            printf("-E- %s line %d: Unknown climatology product\n", __FILE__, __LINE__);
            exit(1);
        }
    }

    if (firstCall) {

        /* get bin information */

        if (!open_l3b(l3file, day)) {
            printf("-E- Level-3 %s needed for the correction\n", pnames[prodID]);
            nofile = 1;
            firstCall = 0;
            exit(1);
        }

        status = VSsetfields(vdata_id[0], "bin_num,nobs,weights,flags_set");
        status = VSsetfields(vdata_id[1], "start_num,begin,extent,max");

        /* allocate space for bin numbers */
        VSinquire(vdata_id[0], &n_records, &interlace0, fields, &vdata_size, vdata_name);
        if ((bins = (int32_t *) malloc(n_records * sizeof (int32_t))) == NULL) {
            printf("-E- %s line %d: Error allocating memory to L3 file bins\n", __FILE__, __LINE__);
            exit(1);
        }

        /* allocate temp space for bin_num (4), nobs (4), weights (4), flags_set (2) */
        vsize = 2 * sizeof (int32) + sizeof (float32) + sizeof (int16);
        if ((bin_buf = (uint8 *) calloc(n_records, vsize)) == NULL) {
            printf("-E- %s line %d: Error allocating memory for L3 file bin info\n", __FILE__, __LINE__);
            exit(1);
        }

        /* get bin info numbers for each filled bin */
        nread = VSread(vdata_id[0], (uint8 *) bin_buf, n_records, interlace0);
        if (nread == -1) {
            printf("-E- %s line %d: Unable to read bin numbers.\n", __FILE__, __LINE__);
            exit(FATAL_ERROR);
        }

        /* copy bin numbers */
        for (i = 0; i < n_records; i++) {
            memcpy((void *) &bins[i], (const void *) &bin_buf[i * vsize], 4);
        }

        close_l3b();
        free(bin_buf);

        firstCall = 0;
    }

    if (nofile) return (-999.0);

    if (!loaded[prodID]) {

        /* load product into static array */

        if (l3file != NULL)
            printf("Loading %s coverage.\n", pnames[prodID]);
        else
            printf("Loading %s climatology.\n", pnames[prodID]);


        if (!open_l3b(l3file, day)) {
            printf("-E- %s line %d: Error opening climatology file\n", __FILE__, __LINE__);
            exit(1);
        }

        initbin();

        /* Calculate the number of L3 products */
        /* ----------------------------------- */
        nprod = MAXNVDATA - 2;
        for (i = 0; i < MAXNVDATA; i++) {
            if (vdata_id[i] == -1) nprod--;
        }

        if (nprod < 1) {
            printf("-E- %s: There are too few products in the L3 file\n", __FILE__);
            exit(1);
        }

        status = VSsetfields(vdata_id[0], "bin_num,nobs,weights,flags_set");
        status = VSsetfields(vdata_id[1], "start_num,begin,extent,max");

        /* locate product vdata & set field names */
        for (i = 0; i < nprod; i++) {
            VSgetname(vdata_id[2 + i], fields);
            if (strncmp(fields, pnames[prodID], 8) == 0) {
                prodinx = i;
                strcpy(buf, fields);
                strcat(buf, "_sum,");
                strcat(buf, fields);
                strcat(buf, "_sum_sq");
                status = VSsetfields(vdata_id[2 + i], buf);
                if (status != 0) {
                    printf("-E- %s: Failed to set field for %s.\n", __FILE__, pnames[prodID]);
                    exit(1);
                }
                break;
            }
        }
        if (i == nprod) {
            printf("-E- %s: Failed to locate %s.\n", __FILE__, pnames[prodID]);
            exit(1);
        }

        /* allocate space for product mean */
        VSinquire(vdata_id[2 + prodinx], &n_records, &interlace, fields, &vdata_size, vdata_name);
        if ((data[prodID] = (float32 *) malloc(n_records * sizeof (float32))) == NULL) {
            printf("-E- %s: Error allocating memory to L3 product %s\n", __FILE__, prodname);
            exit(1);
        }

        /* allocate temp space for sum and sum of squares */
        if ((sum_buf = (float32 *) calloc(2 * n_records, sizeof (float32))) == NULL) {
            printf("-E- %s: Error allocating memory to L3 file products\n", __FILE__);
            exit(1);
        }

        /* allocate temp space for bin numbers, nobs, wts, etc. */
        vsize = 2 * sizeof (int32) + sizeof (float32) + sizeof (int16);
        if ((bin_buf = (uint8 *) calloc(n_records, vsize)) == NULL) {
            printf("-E- %s: Error allocating memory for L3 file bin info\n", __FILE__);
            exit(1);
        }

        /* get bin info for each filled bin */
        nread = VSread(vdata_id[0], (uint8 *) bin_buf, n_records, interlace0);
        if (nread == -1) {
            printf("-E- %s line %d: Unable to read bin numbers.\n", __FILE__, __LINE__);
            exit(1);
        }


        /* get product sum and sum_sq for each filled bin */
        /*
        VSseek(vdata_id[2+prodinx], 0);
         */
        nread = VSread(vdata_id[2 + prodinx], (uint8 *) sum_buf, n_records, interlace);
        if (nread == -1) {
            printf("-E- %s line %d: Unable to read sum/sum_sq\n", __FILE__, __LINE__);
            exit(1);
        }

        for (i = 0; i < n_records; i++) {
            memcpy((void *) &wgt, (const void *) &bin_buf[i * vsize + 6], 4);
            data[prodID][i] = (float32) sum_buf[2 * i] / wgt;
            /*
            printf("%d %f %f %f\n",bins[i],sum_buf[2*i],wgt,data[prodID][i]);
             */
        }

        close_l3b();
        free(bin_buf);
        free(sum_buf);
        closebin();

        loaded[prodID] = 1;
    }


    /* now, find the bin of interest for this product */

    bin_num = latlon2bin(lat, lon);
    bin_idx = nearest_bin(bin_num);

    if (bin_idx >= 0) return (data[prodID][bin_idx]);
    else return (-999.0);
}



