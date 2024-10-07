/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   Function: get_cal_swf

    Reads and applies calibration table

    Return:
        calibrated L1B radiances
	
        Sean Bailey, Futuretech Corporation, 2 Jun 2008 

   Modification history:
        Gene Eplee, SAIC, 17 June 2010		Changed fitting function type
                                                descriptor from sds attribute to
                                                sds fitting parameter variable
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

#include <stdlib.h>
#include "hdf4utils.h"
#include <timeutils.h>
#include "cal_l1a.h"
#include "call1a_proto.h"
#include "genutils.h"
#include "getcal_proto.h"

#include <hdf.h>
#include <mfhdf.h>


float ***alloc3d_float(int m, int n, int o) {
    int i, j;
    float ***p;
    p = malloc(m * sizeof (float *));

    for (i = 0; i < m; i++) {
        p[ i ] = malloc(n * sizeof (float **));
        for (j = 0; j < n; j++) {
            p[ i ][ j ] = malloc(o * sizeof (float ***));
        }
    }

    return (p);
}

float ****alloc4d_float(int l, int m, int n, int o) {
    int i, j, k;
    float ****p;
    p = malloc(l * sizeof (float *));

    for (i = 0; i < l; i++) {
        p[ i ] = malloc(m * sizeof (float **));
        for (j = 0; j < m; j++) {
            p[ i ][ j ] = malloc(n * sizeof (float ***));
            for (k = 0; k < n; k++) {
                p[ i ][ j ][ k ] = malloc(o * sizeof (float ****));

            }

        }
    }

    return (p);
}

int32_t ***alloc3d_long(int m, int n, int o) {
    int i, j;
    int32_t ***p;
    p = malloc(m * sizeof (int32_t *));

    for (i = 0; i < m; i++) {
        p[ i ] = malloc(n * sizeof (int32_t **));
        for (j = 0; j < n; j++) {
            p[ i ][ j ] = malloc(o * sizeof (int32_t ***));
        }
    }

    return (p);
}


#define NBAND 8
#define NTEMP 256
// #define NCOEF 5
#define NCOEF 6
#define NDET 4
#define NGAIN 4
#define NTDI 256
#define NKNEE 5
#define K1 0
#define K3 1
#define K4 2
#define MSIDE 3
#define DARK 5
#define GR 4
#define NK 4

static short ftype[NK];

static short gain3corr;
static float ***radcor;
static int32_t *k1_epoch;
// static int32_t *k1_ftype;
static int num_k1_epoch;

static float ***g3corr;
static int32_t *gr_epoch;
// static int32_t *gr_ftype;
static int num_gr_epoch;

static float ***cnts2rad;
static int32_t ***det_off;

static float ***fptempcor;
static int32_t *k3_epoch;
// static int32_t *k3_ftype;
static int num_k3_epoch;

static float ***scanmod;
static int32_t *k4_epoch;
// static int32_t *k4_ftype;
static int num_k4_epoch;

static float ****msidecor;
static int32_t *ms_epoch;
// static int32_t *ms_ftype;
static int num_ms_epoch;

static float ***dark_restore;
static int32_t *dkrest_epoch;
static int num_dkrest_epoch;
static int haveDarkRestore = 0;

static float fp_temp[NBAND][NTEMP];
static short tdi_list[NDET][NTDI];
static float g_f[NBAND][1024];

static int32_t ref_year;
static int32_t ref_day;
static int32_t ref_min;
static short prev_tdi[NBAND] = {-1, -1, -1, -1, -1, -1, -1, -1};
static short prev_gain[NBAND] = {-1, -1, -1, -1, -1, -1, -1, -1};
static int32_t prev_syear;
static int32_t prev_sday;
static int32_t prev_smsec;
static double ref_jsec;

/* -------------------------------------------------------------------------- */
/* read_caltable() - called once to load cal table static arrays              */

/* -------------------------------------------------------------------------- */

void read_caltable(char *cal_path) {

    static int firstCall = 1;
    char name [H4_MAX_NC_NAME] = "";
    char sdsname[H4_MAX_NC_NAME] = "";

    int32 sd_id;
    int32 sds_id;
    int32 rank;
    int32 nt;
    int32 dims[H4_MAX_VAR_DIMS];
    int32 nattrs, attr_index, count, num_type;
    char attr_name[64];
    int32 start4[4] = {0, 0, 0, 0};
    int32 start[3] = {0, 0, 0};
    int32 status;
    float *r;
    int32_t *dr;

    int32_t i, j, k, l;
    int m = 1;

    if (!firstCall) return;

    /* Open the cal cal_path... */
    sd_id = SDstart(cal_path, DFACC_RDONLY);
    if (sd_id == FAIL) {
        fprintf(stderr, "-E- %s line %d: SDstart(%s, %d) failed.\n",
                __FILE__, __LINE__, cal_path, DFACC_RDONLY);
        exit(1);
    }

    /* Functional type index */
    strcpy(sdsname, "ftype");
    sds_id = SDselect(sd_id, SDnametoindex(sd_id, sdsname));
    status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
    status = SDreaddata(sds_id, start, NULL, dims, ftype);
    if (status != 0) {
        printf("-E- %s Line %d:  Error reading SDS %s from %s.\n",
                __FILE__, __LINE__, sdsname, cal_path);
        exit(1);
    } else {
        status = SDendaccess(sds_id);
    }

    /* Radiometric coefficients */
    strcpy(sdsname, "radiometric_coef");
    sds_id = SDselect(sd_id, SDnametoindex(sd_id, sdsname));
    status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
    radcor = alloc3d_float(dims[0], dims[1], dims[2]);
    m = 1;
    for (i = 0; i < rank; i++)
        m *= dims[i];
    r = (float *) malloc(m * sizeof (float *));
    /* Get epoch start times */
    attr_index = SDfindattr(sds_id, "k1_epochs");
    status = SDattrinfo(sds_id, attr_index, attr_name, &num_type, &count);
    num_k1_epoch = count;
    k1_epoch = (int32_t *) malloc(count * sizeof (int32_t *));
    status = SDreadattr(sds_id, attr_index, k1_epoch);
    /* Get functional type index*/
    /*    attr_index = SDfindattr(sds_id,"ftype index");
    status = SDattrinfo(sds_id,attr_index,attr_name,&num_type, &count);
    k1_ftype = (int32_t *) malloc(count * sizeof(int32_t *));
    status = SDreadattr(sds_id,attr_index,k1_ftype); */

    /* Get the data...*/
    status = SDreaddata(sds_id, start, NULL, dims, r);
    if (status != 0) {
        printf("-E- %s Line %d:  Error reading SDS %s from %s.\n",
                __FILE__, __LINE__, sdsname, cal_path);
        exit(1);
    } else {
        status = SDendaccess(sds_id);
        for (i = 0; i < dims[0]; i++) {
            for (j = 0; j < dims[1]; j++) {
                for (k = 0; k < dims[2]; k++) {
                    m = i * dims[1] * dims[2] + j * dims[2] + k;
                    radcor[i][j][k] = r[m];
                }
            }
        }
    }

    free(r);

    /* Get the dark counts */
    strcpy(sdsname, "dark_counts");
    sds_id = SDselect(sd_id, SDnametoindex(sd_id, sdsname));
    if (sds_id >= 0) {
        haveDarkRestore = 1;
        status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
        dark_restore = alloc3d_float(dims[0], dims[1], dims[2]);
        m = 1;
        for (i = 0; i < rank; i++)
            m *= dims[i];
        r = (float *) malloc(m * sizeof (float *));
        /* Get epoch start times */
        attr_index = SDfindattr(sds_id, "dn0_epochs");
        status = SDattrinfo(sds_id, attr_index, attr_name, &num_type, &count);
        num_dkrest_epoch = count;
        dkrest_epoch = (int32_t *) malloc(count * sizeof (int32_t *));
        status = SDreadattr(sds_id, attr_index, dkrest_epoch);
        status = SDreaddata(sds_id, start, NULL, dims, r);
        if (status != 0) {
            printf("-E- %s Line %d:  Error reading SDS %s from %s.\n",
                    __FILE__, __LINE__, sdsname, cal_path);
            exit(1);
        } else {
            status = SDendaccess(sds_id);
            for (i = 0; i < dims[0]; i++) {
                for (j = 0; j < dims[1]; j++) {
                    for (k = 0; k < dims[2]; k++) {
                        m = i * dims[1] * dims[2] + j * dims[2] + k;
                        dark_restore[i][j][k] = r[m];
                    }
                }
            }
        }
        free(r);

    }


    /* Gain 3:1 correction coefficients */
    strcpy(sdsname, "gainratio_coef");
    sds_id = SDselect(sd_id, SDnametoindex(sd_id, sdsname));
    status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
    if (status == 0) {
        if (want_verbose)
            printf("\nGain 3 correction applied.\n");
        gain3corr = 1;
        g3corr = alloc3d_float(dims[0], dims[1], dims[2]);
        m = 1;
        for (i = 0; i < rank; i++)
            m *= dims[i];
        r = (float *) malloc(m * sizeof (float *));
        /* Get epoch start times */
        attr_index = SDfindattr(sds_id, "gr_epochs");
        status = SDattrinfo(sds_id, attr_index, attr_name, &num_type, &count);
        num_gr_epoch = count;
        gr_epoch = (int32_t *) malloc(count * sizeof (int32_t *));
        status = SDreadattr(sds_id, attr_index, gr_epoch);
        /* Get functional type index*/
        /*	attr_index = SDfindattr(sds_id,"ftype index");
        status = SDattrinfo(sds_id,attr_index,attr_name,&num_type, &count);
        gr_ftype = (int32_t *) malloc(count * sizeof(int32_t *));
        status = SDreadattr(sds_id,attr_index,gr_ftype); */
        /* Get the data...*/
        status = SDreaddata(sds_id, start, NULL, dims, r);
        if (status != 0) {
            printf("-E- %s Line %d:  Error reading SDS %s from %s.\n",
                    __FILE__, __LINE__, sdsname, cal_path);
            exit(1);
        } else {
            status = SDendaccess(sds_id);
            for (i = 0; i < dims[0]; i++) {
                for (j = 0; j < dims[1]; j++) {
                    for (k = 0; k < dims[2]; k++) {
                        m = i * dims[1] * dims[2] + j * dims[2] + k;
                        g3corr[i][j][k] = r[m];
                    }
                }
            }
        }

        free(r);
    } else {
        gain3corr = 0;
    }

    /* Counts to Radiance coefficients */
    strcpy(sdsname, "counts_to_radiance");
    sds_id = SDselect(sd_id, SDnametoindex(sd_id, sdsname));
    status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
    cnts2rad = alloc3d_float(dims[2], dims[1], dims[0]);
    m = 1;
    for (i = 0; i < rank; i++)
        m *= dims[i];
    r = (float *) malloc(m * sizeof (float *));
    /* Get the data...*/
    status = SDreaddata(sds_id, start, NULL, dims, r);
    if (status != 0) {
        printf("-E- %s Line %d:  Error reading SDS %s from %s.\n",
                __FILE__, __LINE__, sdsname, cal_path);
        exit(1);
    } else {
        status = SDendaccess(sds_id);
        for (i = 0; i < dims[0]; i++) {
            for (j = 0; j < dims[1]; j++) {
                for (k = 0; k < dims[2]; k++) {
                    m = i * dims[1] * dims[2] + j * dims[2] + k;
                    cnts2rad[k][j][i] = r[m];
                }
            }
        }
    }

    free(r);


    /* Detector offsets */
    strcpy(sdsname, "detector_offsets");
    sds_id = SDselect(sd_id, SDnametoindex(sd_id, sdsname));
    status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
    det_off = alloc3d_long(dims[2], dims[1], dims[0]);
    m = 1;
    for (i = 0; i < rank; i++)
        m *= dims[i];
    dr = (int32_t *) malloc(m * sizeof (int32_t *));
    /* Get the data...*/
    status = SDreaddata(sds_id, start, NULL, dims, dr);
    if (status != 0) {
        printf("-E- %s Line %d:  Error reading SDS %s from %s.\n",
                __FILE__, __LINE__, sdsname, cal_path);
        exit(1);
    } else {
        status = SDendaccess(sds_id);
        for (i = 0; i < dims[0]; i++) {
            for (j = 0; j < dims[1]; j++) {
                for (k = 0; k < dims[2]; k++) {
                    m = i * dims[1] * dims[2] + j * dims[2] + k;
                    det_off[k][j][i] = dr[m];
                }
            }
        }
    }

    free(dr);

    /* Temperature coefficients */
    strcpy(sdsname, "temperature_coef");
    sds_id = SDselect(sd_id, SDnametoindex(sd_id, sdsname));
    status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
    fptempcor = alloc3d_float(dims[0], dims[1], dims[2]);
    m = 1;
    for (i = 0; i < rank; i++)
        m *= dims[i];
    r = (float *) malloc(m * sizeof (float *));
    /* Get epoch start times */
    attr_index = SDfindattr(sds_id, "k3_epochs");
    status = SDattrinfo(sds_id, attr_index, attr_name, &num_type, &count);
    num_k3_epoch = count;
    k3_epoch = (int32_t *) malloc(count * sizeof (int32_t *));
    status = SDreadattr(sds_id, attr_index, k3_epoch);
    /* Get functional type index*/
    /*    attr_index = SDfindattr(sds_id,"ftype index");
    status = SDattrinfo(sds_id,attr_index,attr_name,&num_type, &count);
    k3_ftype = (int32_t *) malloc(count * sizeof(int32_t *));
    status = SDreadattr(sds_id,attr_index,k3_ftype); */
    /* Get the data...*/
    status = SDreaddata(sds_id, start, NULL, dims, r);
    if (status != 0) {
        printf("-E- %s Line %d:  Error reading SDS %s from %s.\n",
                __FILE__, __LINE__, sdsname, cal_path);
        exit(1);
    } else {
        status = SDendaccess(sds_id);
        for (i = 0; i < dims[0]; i++) {
            for (j = 0; j < dims[1]; j++) {
                for (k = 0; k < dims[2]; k++) {
                    m = i * dims[1] * dims[2] + j * dims[2] + k;
                    fptempcor[i][j][k] = r[m];
                }
            }
        }
    }

    free(r);

    /* Scan Modulation correction*/
    strcpy(sdsname, "scan_mod_coef");
    sds_id = SDselect(sd_id, SDnametoindex(sd_id, sdsname));
    status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
    scanmod = alloc3d_float(dims[0], dims[1], dims[2]);
    m = 1;
    for (i = 0; i < rank; i++)
        m *= dims[i];
    r = (float *) malloc(m * sizeof (float *));
    /* Get epoch start times */
    attr_index = SDfindattr(sds_id, "k4_epochs");
    status = SDattrinfo(sds_id, attr_index, attr_name, &num_type, &count);
    num_k4_epoch = count;
    k4_epoch = (int32_t *) malloc(count * sizeof (int32_t *));
    status = SDreadattr(sds_id, attr_index, k4_epoch);
    /* Get functional type index*/
    /*    attr_index = SDfindattr(sds_id,"ftype index");
    status = SDattrinfo(sds_id,attr_index,attr_name,&num_type, &count);
    k4_ftype = (int32_t *) malloc(count * sizeof(int32_t *));
    status = SDreadattr(sds_id,attr_index,k4_ftype); */
    /* Get the data...*/
    status = SDreaddata(sds_id, start, NULL, dims, r);
    if (status != 0) {
        printf("-E- %s Line %d:  Error reading SDS %s from %s.\n",
                __FILE__, __LINE__, sdsname, cal_path);
        exit(1);
    } else {
        status = SDendaccess(sds_id);
        for (i = 0; i < dims[0]; i++) {
            for (j = 0; j < dims[1]; j++) {
                for (k = 0; k < dims[2]; k++) {
                    m = i * dims[1] * dims[2] + j * dims[2] + k;
                    scanmod[i][j][k] = r[m];
                }
            }
        }
    }

    free(r);

    /* Mirror Side correction*/
    strcpy(sdsname, "mirror_side_coef");
    sds_id = SDselect(sd_id, SDnametoindex(sd_id, sdsname));
    status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
    msidecor = alloc4d_float(dims[0], dims[1], dims[2], dims[3]);
    m = 1;
    for (i = 0; i < rank; i++)
        m *= dims[i];
    r = (float *) malloc(m * sizeof (float *));
    /* Get epoch start times */
    attr_index = SDfindattr(sds_id, "ms_epochs");
    status = SDattrinfo(sds_id, attr_index, attr_name, &num_type, &count);
    num_ms_epoch = count;
    ms_epoch = (int32_t *) malloc(count * sizeof (int32_t *));
    status = SDreadattr(sds_id, attr_index, ms_epoch);
    /* Get functional type index*/
    /*    attr_index = SDfindattr(sds_id,"ftype index");
    status = SDattrinfo(sds_id,attr_index,attr_name,&num_type, &count);
    ms_ftype = (int32_t *) malloc(count * sizeof(int32_t *));
    status = SDreadattr(sds_id,attr_index,ms_ftype); */
    /* Get the data...*/
    status = SDreaddata(sds_id, start4, NULL, dims, r);
    if (status != 0) {
        printf("-E- %s Line %d:  Error reading SDS %s from %s.\n",
                __FILE__, __LINE__, sdsname, cal_path);
        exit(1);
    } else {
        status = SDendaccess(sds_id);
        for (i = 0; i < dims[0]; i++) {
            for (j = 0; j < dims[1]; j++) {
                for (k = 0; k < dims[2]; k++) {
                    for (l = 0; l < dims[3]; l++) {
                        m = i * dims[1] * dims[2] * dims[3] + j * dims[2] * dims[3] + k * dims[3] + l;
                        msidecor[i][j][k][l] = r[m];
                    }
                }
            }
        }
    }

    free(r);

    /* Focal Plane Temperature Array*/
    strcpy(sdsname, "fp_temp");
    sds_id = SDselect(sd_id, SDnametoindex(sd_id, sdsname));
    status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
    /* Get the data...*/
    status = SDreaddata(sds_id, start, NULL, dims, fp_temp);
    if (status != 0) {
        printf("-E- %s Line %d:  Error reading SDS %s from %s.\n",
                __FILE__, __LINE__, sdsname, cal_path);
        exit(1);
    } else {
        status = SDendaccess(sds_id);
    }

    /* TDI List Array*/
    strcpy(sdsname, "TDI_list");
    sds_id = SDselect(sd_id, SDnametoindex(sd_id, sdsname));
    status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
    /* Get the data...*/
    status = SDreaddata(sds_id, start, NULL, dims, tdi_list);
    if (status != 0) {
        printf("-E- %s Line %d:  Error reading SDS %s from %s.\n",
                __FILE__, __LINE__, sdsname, cal_path);
        exit(1);
    } else {
        status = SDendaccess(sds_id);
    }

    /* Get Reference Date info */
    attr_index = SDfindattr(sd_id, "Reference Year");
    status = SDreadattr(sd_id, attr_index, &ref_year);
    if (status != 0) {
        printf("-E- %s Line %d:  Error reading Reference Year from %s.\n",
                __FILE__, __LINE__, cal_path);
        exit(1);
    }
    attr_index = SDfindattr(sd_id, "Reference Day");
    status = SDreadattr(sd_id, attr_index, &ref_day);
    if (status != 0) {
        printf("-E- %s Line %d:  Error reading Reference Day from %s.\n",
                __FILE__, __LINE__, cal_path);
        exit(1);
    }
    attr_index = SDfindattr(sd_id, "Reference Minute");
    status = SDreadattr(sd_id, attr_index, &ref_min);
    if (status != 0) {
        printf("-E- %s Line %d:  Error reading Reference Minute from %s.\n",
                __FILE__, __LINE__, cal_path);
        exit(1);
    }

    ref_jsec = yds2unix((int) ref_year, (int) ref_day, ((double) (ref_min)) / 1440.0);

    /* terminate access to the SD interface and close the cal_path */
    status = SDend(sd_id);

    firstCall = 0;
}

/* --------------------------------------------------------------------------*/
/* apply_time_coef - given the functional form index, coefficients           */
/*                      and "x" parameter, returns "y"                       */

/* --------------------------------------------------------------------------*/
float apply_time_coef(short functype, float coef[NCOEF], double param) {
    float correction_factor;

    switch (functype) {
        // default functype (0) - double exponential
    case 0:
        correction_factor = 1.0 / (coef[0] -
                coef[1] * (1.0 - exp(-coef[2] * param)) -
                coef[3] * (1.0 - exp(-coef[4] * param)));
        break;
        // functype (1) - quadratic
    case 1:
        correction_factor = 1.0 / (coef[0] + coef[1] * param + coef[2] * pow(param, 2));
        break;
        // functype (2) - exponential + linear
    case 2:
        correction_factor = 1.0 / (coef[0] -
                coef[1] * (1.0 - exp(-coef[2] * param)) -
                coef[3] * param);
        break;
    }

    return correction_factor;
}

/* --------------------------------------------------------------------------*/
/* get_epoch_idx - returns the index for the epoch                           */

/* --------------------------------------------------------------------------*/
int get_epoch_idx(int param, double jsec) {
    int epidx = 0;
    int32_t *epochs;
    int cnt, i;

    switch (param) {
    case K1:
        epochs = k1_epoch;
        cnt = num_k1_epoch;
        break;
    case GR:
        epochs = gr_epoch;
        cnt = num_gr_epoch;
        break;
    case K3:
        epochs = k3_epoch;
        cnt = num_k3_epoch;
        break;
    case K4:
        epochs = k4_epoch;
        cnt = num_k4_epoch;
        break;
    case MSIDE:
        epochs = ms_epoch;
        cnt = num_ms_epoch;
        break;
    case DARK:
        epochs = dkrest_epoch;
        cnt = num_dkrest_epoch;
        break;
    }
    for (i = 0; i < cnt; i++) {
        if (jsec > (double) epochs[i]) {
            epidx = i;
        }
    }
    return epidx;
}

/* --------------------------------------------------------------------------*/
/* calc_knees_two - slightly modified version of the original...             */

/* --------------------------------------------------------------------------*/
void calc_knees_two(int16 *tdi, float32 counts[8][4][5], float32 rads[8][4][5]) {

    int16 dets[4];
    int32 i, j, k;
    int32 scnts[4]; /* saturation counts      */
    float32 srads[4]; /* saturation radiance    */
    float32 loc_slopes[4];
    float32 slopes[NBAND][4][4];
    int32 cnts[NBAND][4][4];
    int32 oindex[NDET];


    for (i = 0; i < NBAND; i++)
        for (j = 0; j < 4; j++)
            for (k = 0; k < 4; k++) {
                slopes[i][j][k] = cnts2rad[i][j][k];
                cnts[i][j][k] = det_off[i][j][k];
            }

    for (i = 0; i < NBAND; i++) {
        for (j = 0; j < 4; j++)
            dets[j] = tdi_list[j][tdi[i]] - 1;
        for (j = 0; j < NGAIN; j++) {
            for (k = 0; k < NDET; k++) {
                scnts[k] = 1023 - cnts[i][j][dets[k]];
                srads[k] = scnts[k] * slopes[i][j][dets[k]];
                loc_slopes[k] = slopes[i][j][dets[k]];
            }

            sort_srads(srads, oindex);

            rads[i][j][0] = 0;
            for (k = 1; k < 5; k++)
                rads[i][j][k] = srads[oindex[k - 1]];

            counts[i][j][0] = 0;
            counts[i][j][1] = (scnts[oindex[0]] +
                    srads[oindex[0]] / loc_slopes[oindex[1]] +
                    srads[oindex[0]] / loc_slopes[oindex[2]] +
                    srads[oindex[0]] / loc_slopes[oindex[3]]) / 4.0;

            counts[i][j][2] = (scnts[oindex[0]] + scnts[oindex[1]] +
                    srads[oindex[1]] / loc_slopes[oindex[2]] +
                    srads[oindex[1]] / loc_slopes[oindex[3]]) / 4.0;


            counts[i][j][3] = (scnts[oindex[0]] + scnts[oindex[1]] +
                    scnts[oindex[2]] +
                    srads[oindex[2]] / loc_slopes[oindex[3]]) / 4.0;

            counts[i][j][4] = (scnts[oindex[0]] + scnts[oindex[1]] +
                    scnts[oindex[2]] + scnts[oindex[3]]) / 4.0;

        }
    }
}
/* -------------------------------------------------------------------------- */
/* l1b_rad - returns calibrated L1B radiances                                 */

/* -------------------------------------------------------------------------- */
int32_t l1b_rad(int syear, int sday, int32_t smsec, int32_t msec,
        char *dtype, int32_t nsta, int32_t ninc, int32_t npix,
        float *dark_mean, short *gain, short *tdi,
        short *scan_temp, float *inst_temp, int mside,
        short *l1a_data, float *l1b_data, cal_mod_struc *cal_mod) {
    int i, pixel, knee, n, count;
    int band;
    int l1_data_lo, l1_data_hi;
    int dark;
    float dark_ratio;
    double delta_t;
    int gn[NBAND];
    float k1, gcorr, k3, k4, mirror;
    int16 k1_ftype, gr_ftype, k4_ftype, ms_ftype;
    int epoch_idx[NCOEF];
    double jsec;
    float coef[NCOEF];
    float cal_offset = 0.;
    int gain_flag, tdi_flag;
    float slope;
    short count1, count2;
    short scanpix;
    int called_calc_knees = 0;

    jsec = yds2unix(syear, sday, ((double) (msec)) / 1000.0);
    delta_t = (jsec - (double) ref_jsec) / 86400.;

    // TDI, gain and scan_temp funny business...Gene E. got some 'spalin' to do ;)
    for (band = 0; band < NBAND; band++) {
        if (tdi[band] < 0)
            tdi[band] = 0;
        if (tdi[band] > 255)
            tdi[band] = 255;
        if (scan_temp[band] < 0)
            scan_temp[band] = 0;
        if (scan_temp[band] > 255)
            scan_temp[band] = 255;
        // Fix gain 2 / gain 3 telemetry bit flip
        switch (gain[band]) {
        case 0:
            gn[band] = 0;
            break;
        case 1:
            gn[band] = 2;
            break;
        case 2:
            gn[band] = 1;
            break;
        case 3:
            gn[band] = 3;
            break;
        }

    }

    for (band = 0; band < NBAND && tdi[band] == prev_tdi[band]; band++);
    if (band < NBAND)
        tdi_flag = 1;
    else
        tdi_flag = 0;

    // Calc knees ...

    if (syear != prev_syear || sday != prev_sday
            || smsec != prev_smsec || tdi_flag) {
        prev_syear = syear;
        prev_sday = sday;
        prev_smsec = smsec;

        for (band = 0; band < NBAND; band++)
            prev_tdi[band] = tdi[band];

        calc_knees_two(tdi, cal_counts, cal_rads);
        called_calc_knees = 1;
    }


    for (band = 0; band < NBAND && gn[band] == prev_gain[band]; band++);

    if (band < NBAND)
        gain_flag = 1;
    else
        gain_flag = 0;
    // ... and then generate the 'radiance' LUT...
    if (called_calc_knees || gain_flag) {
        called_calc_knees = 0;
        for (band = 0; band < NBAND; band++)
            prev_gain[band] = gn[band];
        for (band = 0; band < NBAND; band++) {
            for (knee = 1; knee <= 4; knee++) {
                n = 1;
                while (((int16) cal_counts[band][gn[band]][knee] ==
                        (int16) cal_counts[band][gn[band]][knee - n]) && n <= knee)
                    n++;
                count1 = (int16) cal_counts[band][gn[band]][knee - n] + 1;
                count2 = (int16) cal_counts[band][gn[band]][knee];
                if (knee == 1)
                    count1 = 0;
                if (knee == 4)
                    count2 = 1023;
                slope = (cal_rads[band][gn[band]][knee] -
                        cal_rads[band][gn[band]][knee - n]) /
                        (cal_counts[band][gn[band]][knee] -
                        cal_counts[band][gn[band]][knee - n]);
                for (count = count1; count <= count2; count++)
                    g_f[band][count] =
                        slope * (count - cal_counts[band][gn[band]][knee - n]) +
                    cal_rads[band][gn[band]][knee - n];
            }
        }
    }

    // ...and now for the real heavy lifting...
    epoch_idx[0] = get_epoch_idx(K1, jsec);
    epoch_idx[1] = get_epoch_idx(K3, jsec);
    epoch_idx[2] = get_epoch_idx(MSIDE, jsec);
    epoch_idx[3] = get_epoch_idx(K4, jsec);
    epoch_idx[4] = get_epoch_idx(GR, jsec);

    for (band = 0; band < NBAND; band++) {

        // Radiometric decay correction
        for (i = 0; i < NCOEF; i++) {
            coef[i] = radcor[epoch_idx[0]][i][band];
        }
        k1_ftype = (int16) coef[NCOEF - 1];
        k1 = apply_time_coef(k1_ftype, coef, delta_t);

        gcorr = 1.;
        if (gain3corr) {
            for (i = 0; i < NCOEF; i++) {
                coef[i] = g3corr[epoch_idx[4]][i][band];
            }
            // gain3corr is NOT applied in the inverse...so invert it :)
            gr_ftype = (int16) coef[NCOEF - 1];
            gcorr = 1. / apply_time_coef(gr_ftype, coef, delta_t);
        }
        /* Check for override on system gain and offset */
        /* and pass back the gain, offset used  */
        if (cal_mod->flag == 1) {
            k1 = cal_mod->gain[band] * k1;
        } else if (cal_mod->flag == 2) {
            cal_offset = cal_mod->offset[band];
        }
        else if (cal_mod->flag == 3) {
            k1 = cal_mod->gain[band] * k1;
            cal_offset = cal_mod->offset[band];
        } else {
            cal_mod->gain[band] = k1;
            cal_mod->offset[band] = cal_offset;
        }

        // Temperature correction
        for (i = 0; i < NCOEF; i++)
            coef[i] = fptempcor[epoch_idx[1]][i][band];
        // k3_ftype = (int16)fptempcor[epoch_idx[1]][NCOEF][band];
        // k3 = apply_time_coef(k3_ftype,coef,delta_t);
        k3 = (1.0 + coef[0]*(fp_temp[band][scan_temp[band]] - coef[1]));

        // Mirror side correction
        for (i = 0; i < NCOEF; i++)
            coef[i] = msidecor[epoch_idx[2]][i][mside][band];
        ms_ftype = (int16) coef[NCOEF - 1];
        mirror = apply_time_coef(ms_ftype, coef, delta_t);
        if (mside == 1) mirror = 1.0 / mirror;

        // Scan modulation (a.k.a. RVS) correction
        for (i = 0; i < NCOEF; i++)
            coef[i] = scanmod[epoch_idx[3]][i][band];
        k4_ftype = (int16) coef[NCOEF - 1];
        for (pixel = 0; pixel < npix; pixel++) {

            scanpix = nsta + ninc*pixel;
            //  apply scan modulation correction to earth view data only...
            if ((strcmp(dtype, "SOL") != 0) && (strcmp(dtype, "TDI") != 0) &&
                    (strcmp(dtype, "IGC") != 0))
                k4 = apply_time_coef(k4_ftype, coef, scanpix - 643);
            else k4 = 1.0;

            // If the cal table has the dark_count array, use it.
            // If not, use the mean of the scene
            if (haveDarkRestore) {
                int dkidx = get_epoch_idx(DARK, jsec);

                dark = ceil(dark_restore[dkidx][gn[band]][band]);
                dark_ratio = (float) dark - dark_restore[dkidx][gn[band]][band];

            } else {
                dark = ceil(dark_mean[band]);
                dark_ratio = (float) dark - dark_mean[band];
            }

            l1_data_lo = l1a_data[band * npix + pixel] - dark;

            if (l1_data_lo < 0)
                l1_data_lo = 0;
            if (l1_data_lo > 1022)
                l1_data_lo = 1022;

            l1_data_hi = l1_data_lo + 1;

            l1b_data[band * npix + pixel] = (k1 * gcorr * k3 * k4 * mirror * (g_f[band][l1_data_lo]*(1 - dark_ratio) + g_f[band][l1_data_hi] * dark_ratio) + cal_offset);

        }

    }

    return SUCCEED;

}

