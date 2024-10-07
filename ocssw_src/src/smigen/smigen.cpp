//   This program reads a l3_bin land product in a sinusoidal projection and
//   writes either a SMI file in either sinusoidal or platte carre projection or
//   a true color pgm file.

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <libgen.h>
#include <stdlib.h>
#include <stdint.h>
#include <sys/types.h>
#include <sys/utsname.h>

#include "netcdf.h"  // Needs to be first to define netcdf stuff JMG
#include "smiinc.h"
#include "smigen_input.h"
#include "palette.h"
#include "hdf.h"
#include "hdf5.h"
//#include "netcdf.h"
#include "hdf_bin.h"

#include "smi_maplists.h"

#include "genutils.h"
#include <timeutils.h>

#define PI2 1.570796326794897
#define PI  3.141592653589793

#define CMD_ARGS    "p:r:w:g:"   /* Valid commandline options */
#define BYTE    unsigned char

#define MAX_NUM_INPUTROW 1000

#define VERSION L3M_SOFTVER_VAL

using namespace std;

int32 nrows;

int32 L3M_PARAMS = 0;
char **parmname_list;
char **parmname_short;
char **unit_list;
char **scaling_list;
float32 *maximum_list;
float32 *minimum_list;
char **palette_list;
char **precision_list;

extern "C" int32 read_attrs(int32, meta_l3bType *);

extern "C" int getlut_file(char * lut_file, short *rlut, short *glut, short *blut);

//extern "C" int32_t put_smi(char *l3m_path,
int32_t put_smi(char *l3m_path,
        char *l3m_name,
        uint8 *l3m_data,
        int32 *dim_sizes,
        float32 *lat_range,
        float32 *lon_range,
        char *measure,
        char *scale_type,
        float32 *si_used,
        float32 *aminmax,
        char *atype,
        char *aopt,
        char *infiles,
        float32 *l3m_dminmax,
        meta_l3bType *meta_l3b,
        unsigned char *map_palette,
        char *softid,
        char *proc_con,
        instr input,
        char *precision,
        /*              float32 *si8_used,*/
        uint8 *qual_byt,
        uint8 isHDF5,
        VOIDP fill);

extern "C" int smigen_input(int argc, char **argv, instr *input);
extern "C" void set_param_string(instr *input_str);

/*
  Revision 5.20 10/31/16
  Add support for CF compliant metadata
  J. Gales

  Revision 5.10 08/03/15
  Add support for L2C SMAP binfile input
  J. Gales

  Revision 5.04 06/06/15
  Fill "Product Type" attribute for Aquarius
  J. Gales

  Revision 5.03 08/08/14
  Check whether output product exists in input binfile
  If not then exit with error
  J. Gales

  Revision 5.02 03/05/14
  Add support for additional CF metadata
  Modify isHDF5/isCDF4 code to handle new netcdf4 library
  Check for Aquarius mission attribute if possible
  J. Gales

  Revision 5.01 02/28/14
  Add support for l3m data compression
  J. Gales

  Revision 5.00 02/14/14
  Add support for NETCDF4 output
  J. Gales

  Revision 4.44 10/04/13
  Revert default stype to 0
  Force stype to be specified for non-default products
  J. Gales

  Revision 4.43 09/27/13
  Default stype to 1
  Add user-defined resolutions: ukm, udeg
  J. Gales

  Revision 4.42 06/27/13
  Check whether specified product exists in NETCDF4 file
  Call read_attrs()
  J. Gales

  Revision 4.41 06/24/13
  Fix NETCDF4 metadata bug
  Call read_attrs()
  J. Gales

  Revision 4.40 06/14/13
  Add support for NETCDF4 i/o
  Add support for HEALPIX
  J. Gales

  Revision 4.30 01/30/13
  Add CF-compliant metadata (HDF5)
  J. Gales

  Revision 4.29 09/21/12
  Test for HDF5 with error handling off
  J. Gales

  Revision 4.28 08/08/12
  Delare "k" variable as int64_t to handle 1km resolution.
  J. Gales

  Revision 4.27 02/12/12
  Print bad resolution in "Improper resolution type" statement.
  J. Gales

  Revision 4.26 09/21/11
  Add support for arctangent scaling (SSS)
  J. Gales

  Revision 4.25 07/28/11
  Incorporate quality fix for gap-fill from Susan Walsh (Miami)
  J. Gales

  Revision 4.20 07/18/11
  Use hdf_bin.h rather than hdf5_bin.h
  J. Gales

  Revision 4.18 09/03/10
  Exit if bad vdata reads
  J. Gales

  Revision 4.17 06/21/10
  Reduce scl_fac for F precision to 100000 to mininize chance of integer
  overflow for large values.
  J. Gales

  Revision 4.16 06/05/10
  Fix metadata for HDF5 files
  J. Gales

  Revision 4.14 12/17/09
  Trap NaN in sum_buf[k*2]
  J. Gales

  Revision 4.13 10/02/09
  Use endianess() rather than u_name.sysname to check for endian
  J. Gales

  Revision 4.12 09/01/09
  Include the scaling attributes for float
  Change resolution default to 9km (not SMI)
  Add fill_value global attribute
  B. Franz

  Revision 4.11 08/28/09
  Fix bad ifile behavior
  J. Gales

  Revision 4.10 08/11/09
  Add processing version
  Remover rflag
  J. Gales

  Revision 4.00 03/12/09
  Add support for HDF5 input/output
  J. Gales

  Revision 3.60 07/16/07
  Fix 16-bit output for land products
  Skip "white" coastline for 16-bit land products 
  J. Gales

  Revision 3.40 04/12/07
  LONWEST/LONEAST not allowed with SEAMLON
  Handle cases where min_out_col > max_out_col (regions split by dateline)
  J. Gales

  Revision 3.38 04/11/07
  Fix free() bug for seam_lon != -180
  J. Gales

  Revision 3.37 02/07/07
  Set default value of precision parameter to "B"
  J. Gales

  Revision 3.36 11/07/06
  Set palette to default if not specified for non-prodtable product
  Set proddesc to "Not Specified" if not specified for non-prodtable product
  Set precison to "B" if not specified for non-prodtable product
  J. Gales

  Revision 3.35 11/07/06
  Write quality vdata if FLOAT precision
  J. Gales

  Revision 3.33 03/28/06
  Check whether input_row > MAX_NUM_INPUTROW
  J. Gales

  Revision 3.32 03/23/06
  Fix unsigned int print bug for total # of grid points
  J. Gales

  Revision 3.31 03/17/06
  Request palette name if product not found in product list
  J. Gales

 */

void usage(const char *progname) {
    printf("%s %s (%s %s)\n",
            progname, VERSION, __DATE__, __TIME__);

    printf("\nUsage: %s ifile ofile prod-name\n", progname);
    printf("   par        = parameter filename\n");
    printf("   ifile      = input bin filename\n");
    printf("   ofile      = output map filename\n");
    printf("   oformat    = output format: 1 (HDF4 [default]), 2 (netCDF4), 3 (HDF5)\n");
    printf("   deflate    = apply internal compression for netCDF output\n");
    printf("   prod       = product name\n");
    printf("   precision  = output map precision: 'B' (default), 'I', 'F'\n");
    printf("   palfile    = palette filename\n");
    printf("   pversion   = processing version (default is Unspecified)\n");
    printf("   datamin    = minumum value for data scaling (default is prod-specific)\n");
    printf("   datamax    = maximum value for data scaling (default is prod-specific)\n");
    printf("   stype      = scaling type,1=LINEAR,2=LOG (default is prod-specific)\n");
    printf("   meas       = measurement to map, 1=mean, 2=var, 3=stdev, 4=pixels, 5=scenes (default=1)\n");
    printf("   loneast    = Easternmost longitude (default=+180)\n");
    printf("   lonwest    = Westernmost longitude (default=-180)\n");
    printf("   latnorth   = Northernmost latitude (default=+90)\n");
    printf("   latsouth   = Southernmost latitude (default=-90)\n");
    printf("   projection = SIN | RECT (default=RECT)\n");
    printf("   resolution = 36km | 18km | 9km | 4km | 2km | 1km | hkm | qkm\n");
    printf("                1deg (one deg) | hdeg (0.5 deg) | qdeg (0.25 deg)\n");
    printf("                10deg (0.1 deg) |udeg-#.# (#.# deg)| ukm-#.# (#.# km) (default=9km)\n");
    printf("   seam_lon   = Longitude of Left Edge of Map (default=-180)\n");
    printf("   proddesc   = Product Description\n");
    printf("   units      = Product Units\n");
    exit(0);
}

int main(int argc, char **argv) {
    int32 i, j, l, jm1, tk;
    int64_t k;
    int32 out_row;
    int32 prev_input_row;
    int32 input_row;
    int32 input_row_array[MAX_NUM_INPUTROW];
    int32 *in2out;

    int32 fid;
    int32 sdfid;
    int32 sds_id;
    int32 vgid;
    int32 vdata_id[4];
    int32 status;
    int32 dims_in[8];
    int32 rank;
    int32 nattrs;
    int32 dims_out[2];
    int32 dims_out_sub[2];
    int32 type;
    int32 max_rowbuf;
    int32 start[2], edges[2];
    int32 *out_sum;
    int32 *out_num;
    int32 *out_sum_rot;
    int32 *out_num_rot;
    int32 input_col;
    int32 pixmin;
    int32 pixmax;
    int32 binindex_buf[4];
    int32 nread;
    int32 offset;
    int32 bin_num;
    int32 min_out_col, max_out_col;
    int32 min_out_row, max_out_row;
    int32 sub_scan[2] = {0, 0};
    int32 sub_pixl[2];
    int32 prod_num = -1;
    int32 filled_data_bins = 0;
    int32 rotate;

    uint32 tgp;

    float32 ratio;
    float32 cos_theta_output;
    float32 cos_theta_input;
    float32 cos_fac;
    float32 flt_val;
    float32 wgt;
    float32 *sum_buf = NULL;
    float32 scl_fac;
    float32 latf, lonf;

    float32 slope, intercept;

    float64 pix_sz=0;

    int32 *input_databuf;
    int32 input_value;
    int16 nobs;
    int16 nscenes;
    int16 *input_databuf_16;
    int16 cur_ptr;
    int16 last_ptr;
    int16 lon, lat;
    uint16 *ptr_arr;
    uint16 maxval;
    uint8 *input_qualbuf;
    uint8 *qual_buf = NULL;
    uint8 input_qual;
    uint8 *best_qual;
    uint8 *best_qual_rot;

    uint8 *par_byt;
    uint8 *qual_byt = NULL;
    uint16 scale_val_16b;
    char buf[128];
    char prodtablename[128];
    uint8 *binlist_buf = NULL;
    uint8 fill_val = 255;
    uint16 fill_val_int16 = 65535;
    float32 fill_val_float = -32767.0;
    uint8 i8;
    uint8 *ptr_i8;
    VOIDP fill=nullptr;

    uint8 land_input = 0;
    uint8 hdf4_input = 0;
    uint8 hdf5_input = 0;
    uint8 ncdf4_input = 0;
    uint8 healpix_input = 0;
    uint8 default_palfile = 0;
    meta_l3bType meta_l3b;
    meta_l3bType *ptr_meta_l3b;
    instr input;

    float32 si_used[2], aminmax[2], dminmax[2] = {+1e10, -1e10}, f32;

    char atype[SMI_MAX_STR_SHORT];
    char aopt[SMI_MAX_STR_SHORT];
    char scale_type[SMI_MAX_STR_SHORT];
    char ptime[SMI_MAX_STR_SHORT];

    char proc_con[SMI_MAX_STR];
    char softid[SMI_MAX_STR];

    char *tmp_str;
    char *cptr;

    int16 * mixed_buf[360 * 180];
    int16 mixed[128][8];

    short *r, *g, *b;

    static unsigned int pow2[16] = {1, 2, 4, 8, 16, 32, 64, 128, 256,
        512, 1024, 2048, 4096, 8192, 16384, 32768};
    time_t tnow;
    struct tm *tmnow;

    struct utsname u_name;

    FILE *fp;
    FILE *mask_fp;

    hid_t binlist_tid;
    hid_t binindex_tid;
    hid_t bindata_tid;
    hid_t bin_dataset_idx;
    hsize_t one = 1;

    int get_minmax_rowcol(float32 [], float32 [], char *, int32,
            int32 *, int32 *, int32 *, int32 *);

    int32 open_input_hdf(char *, char *, int32 *, int32 *,
            int32 *, int32 [], meta_l3bType *);

    int32 open_input_hdf5(char *hdf5_file, char *pname,
            Hdf::hdf5_bin *input_binfile,
            hid_t * bin_dataset_idx);

    int32 open_input_ncdf4(int ncid, const char *pname,
            size_t *nrows, int *grpid,
            int *binindex_id, int *binlist_id,
            int *bindata_id, int *binqual_id,
            meta_l3bType * meta_l3b);

    /*
      row_inp = Nrow_inp * (0.5 - (theta/180))
      row_out = Nrow_out * (0.5 - (theta/180))

      Solving for theta/180 we get:

      row_out = (Nrow_out/Nrow_inp) * row_inp


      col_inp = Npix_inp * (0.5 + (phi/360) * cos(theta_inp))  "SIN"

      col_out = Npix_out * (0.5 + (phi/360) * cos(theta_out))  "SIN"
      col_out = Npix_out * (0.5 + (phi/360))                   "RECT"

      Solving for phi/360 we get:

      col_out = Npix_out * (0.5 * (1 - fac) + (col_inp/Npix_inp) * fac)

      where fac = cos(theta_out) / cos(theta_inp)              "SIN"
      where fac = 1 / cos(theta_inp)                           "RECT"

     */

    setlinebuf(stdout);

    get_time(ptime);

    if (smigen_input(argc, argv, &input) != 0) {
        usage(argv[0]);
        exit(1);
    }

    printf("%s %s (%s %s)\n", argv[0], VERSION, __DATE__, __TIME__);

    strcpy(proc_con, argv[0]);
    for (i = 1; i < argc; i++) {
        strcat(proc_con, " ");
        strcat(proc_con, argv[i]);
    }

    if (input.loneast < -180 || input.loneast > +180) {
        printf("LONEAST must be between -180 and +180.\n");
        exit(1);
    }

    if (input.lonwest < -180 || input.lonwest > +180) {
        printf("LONWEST must be between -180 and +180.\n");
        exit(1);
    }


    if ((tmp_str = getenv("OCDATAROOT")) == NULL) {
        printf("OCDATAROOT environment variable is not defined.\n");
        return (1);
    }
    strcpy(buf, tmp_str);
    strcat(buf, "/common/smigen_defaults.par");
    fp = fopen(buf, "r");

    if (fp == 0x0) {
        printf("SMIGEN defaults file: %s not found.\n", buf);
        printf("Using defaults for product table and palette directory.\n");

        if (strcmp(input.palfile, "DEFAULT") == 0) {
            strcpy(input.palfile, tmp_str);
            strcat(input.palfile, "/common/palette");
            default_palfile = 1;
        }

        strcpy(prodtablename, tmp_str);
        strcat(prodtablename, "/common/smigen_product_table.dat");

    } else {
        fgets(buf, 128, fp);
        if (strcmp(input.palfile, "DEFAULT") == 0) {
            buf[strlen(buf) - 1] = 0;
            cptr = strchr(buf, '=') + 1;
            strcpy(input.palfile, tmp_str);
            strcat(input.palfile, "/");
            strcat(input.palfile, cptr);
            default_palfile = 1;
        }
        fgets(buf, 128, fp);
        buf[strlen(buf) - 1] = 0;
        cptr = strchr(buf, '=') + 1;
        strcpy(prodtablename, tmp_str);
        strcat(prodtablename, "/");
        strcat(prodtablename, cptr);

        // kludge for adding deflate to the smigen_defaults.par file
        // since this code doesn't use clo...
        fgets(buf, 128, fp);
        buf[strlen(buf) - 1] = 0;
        cptr = strchr(buf, '=') + 1;
        input.deflate = atoi(cptr);
        fclose(fp);
    }


    /* Read product table */
    fp = fopen(prodtablename, "r");
    if (fp == 0x0) {
        printf("SMIGEN product table \"%s\" not found.\n", buf);
        exit(1);
    }
    while (fgets(buf, 128, fp) != NULL) {
        if ((buf[0] >= 0x41) && (buf[0] <= 0x5a)) L3M_PARAMS++;
    }
    fseek(fp, 0, SEEK_SET);

    parmname_list = (char**) calloc(L3M_PARAMS, sizeof (char *));
    parmname_short = (char**) calloc(L3M_PARAMS, sizeof (char *));
    unit_list = (char**) calloc(L3M_PARAMS, sizeof (char *));
    scaling_list = (char**) calloc(L3M_PARAMS, sizeof (char *));
    palette_list = (char**) calloc(L3M_PARAMS, sizeof (char *));
    maximum_list = (float32 *) calloc(L3M_PARAMS, sizeof (float32));
    minimum_list = (float32 *) calloc(L3M_PARAMS, sizeof (float32));
    precision_list = (char **) calloc(L3M_PARAMS, sizeof (char *));

    i = 0;
    while (fgets(buf, 128, fp) != NULL) {
        if ((buf[0] >= 0x41) && (buf[0] <= 0x5a)) {

            cptr = strtok(buf, ":");
            parmname_list[i] = (char*) malloc(strlen(cptr) + 1);
            strcpy(parmname_list[i], cptr);

            cptr = strtok(NULL, ":");
            parmname_short[i] = (char*) malloc(strlen(cptr) + 1);
            strcpy(parmname_short[i], cptr);

            cptr = strtok(NULL, ":");
            unit_list[i] = (char*) malloc(strlen(cptr) + 1);
            strcpy(unit_list[i], cptr);

            cptr = strtok(NULL, ":");
            scaling_list[i] = (char*) malloc(strlen(cptr) + 1);
            strcpy(scaling_list[i], cptr);

            cptr = strtok(NULL, ":");
            minimum_list[i] = (float32) atof(cptr);

            cptr = strtok(NULL, ":");
            maximum_list[i] = (float32) atof(cptr);

            cptr = strtok(NULL, ":");
            precision_list[i] = (char*) malloc(strlen(cptr) + 1);
            strcpy(precision_list[i], cptr);

            cptr = strtok(NULL, "\n");
            palette_list[i] = (char*) malloc(strlen(cptr) + 1);
            strcpy(palette_list[i], cptr);

            i++;
        }
    }
    fclose(fp);


    if (strncmp(input.prod, "refl", 4) == 0) {

        prod_num = -2;
        input.datamin = 0.0;
        if (input.datamax == 0.0) input.datamax = 0.394863;

        intercept = 0.0;
        slope = 0.0;

    } else if (strncmp(input.prod, "rhos", 4) == 0) {

        prod_num = -3;
        input.datamin = 0.0;
        if (input.datamax == 0.0) input.datamax = 0.394863;

        intercept = 0.0;
        slope = 0.0;

    } else {

        for (i = 0; i < L3M_PARAMS; i++) {
            if (strcmp(parmname_short[i], input.prod) == 0) {
                cout << i << " " << parmname_short[i] << endl;
                prod_num = i;

                strcpy(scale_type, scaling_list[i]);
                if (input.stype == 1) strcpy(scale_type, "linear");
                if (input.stype == 2) strcpy(scale_type, "logarithmic");
                if (input.stype == 3) strcpy(scale_type, "arctan");
                if (input.stype == 0) {
                    if (strcmp("linear", scale_type) == 0) {
                        input.stype = 1;
                    } else if (strcmp("logarithmic", scale_type) == 0) {
                        input.stype = 2;
                    } else if (strcmp("arctan", scale_type) == 0) {
                        input.stype = 3;
                    }
                }

                if (input.datamin == 0.0) input.datamin = minimum_list[i];
                if (input.datamax == 0.0) input.datamax = maximum_list[i];

                if (input.precision[0] == 0)
                    strcpy(input.precision, precision_list[i]);

                /* Set Maximum value */
                /* ----------------- */
                if (strcmp(input.precision, "I") == 0) {
                    maxval = 65534;
                } else {
                    maxval = 250;
                }

                if (strcmp(scale_type, "linear") == 0) {

                    strcpy(scale_type, "LINEAR");

                    intercept = input.datamin;
                    slope = (input.datamax - intercept) / maxval;
                }

                if (strcmp(scale_type, "logarithmic") == 0) {

                    strcpy(scale_type, "LOG");

                    intercept = log10(input.datamin);
                    slope = (log10(input.datamax) - intercept) / maxval;
                }

                if (strcmp(scale_type, "arctan") == 0) {
                    strcpy(scale_type, "ATAN");
                }

                if (input.proddesc[0] == 0)
                    strcpy(input.proddesc, parmname_list[i]);

                if (input.units[0] == 0)
                    strcpy(input.units, unit_list[i]);

                /* Read palette file */
                if (default_palfile) {
                    strcat(input.palfile, "/");
                    strcat(input.palfile, palette_list[i]);
                    strcat(input.palfile, ".pal");
                }

                if (!(r = (short *) calloc(256, sizeof (short)))) {
                    fprintf(stderr, "smigen: Error allocating space for red.\n");
                    return -1;
                };
                if (!(g = (short *) calloc(256, sizeof (short)))) {
                    fprintf(stderr, "smigen: Error allocating space for green.\n");
                    return -1;
                };
                if (!(b = (short *) calloc(256, sizeof (short)))) {
                    fprintf(stderr, "smigen: Error allocating space for blue.\n");
                    return -1;
                };

                if (getlut_file(input.palfile, r, g, b)) {
                    fprintf(stderr, "Error reading palette file %s\n", input.palfile);
                    free(r);
                    free(g);
                    free(b);
                    return -1;
                }
                for (i = 0; i < 256; i++) {
                    input.palette[i * 3] = r[i];
                    input.palette[i * 3 + 1] = g[i];
                    input.palette[i * 3 + 2] = b[i];
                }
                free(r);
                free(g);
                free(b);

                break;
            }
        }
    }



    /* Non-Product Table product */
    /* ------------------------- */
    if (prod_num == -1) {
        printf("Product: \"%s\" not found in default product list.\n\n", input.prod);

        if (input.datamin == 0 && input.datamax == 0 && input.stype == 0) {
            printf("Make sure \"datamin\", \"datamax\" and \"stype\"\n");
            printf("are specified on the command line.\n");
            exit(1);
        }
        if (input.proddesc[0] == 0) strcpy(input.proddesc, "Not Specified");
        if (input.precision[0] == 0) strcpy(input.precision, "B");

        if (input.stype == 1) strcpy(scale_type, "linear");
        if (input.stype == 2) strcpy(scale_type, "logarithmic");
        if (input.stype == 3) strcpy(scale_type, "arctan");

        /* Set Maximum value */
        /* ----------------- */
        if (strcmp(input.precision, "I") == 0) {
            maxval = 65534;
        } else {
            maxval = 250;
        }

        if (strcmp(scale_type, "linear") == 0) {

            strcpy(scale_type, "LINEAR");

            intercept = input.datamin;
            slope = (input.datamax - intercept) / maxval;
        }

        if (strcmp(scale_type, "logarithmic") == 0) {

            strcpy(scale_type, "LOG");

            intercept = log10(input.datamin);
            slope = (log10(input.datamax) - intercept) / maxval;
        }

        if (strcmp(scale_type, "arctan") == 0) {
            strcpy(scale_type, "ATAN");
        }

        /* Read palette file */

        if (!(r = (short *) calloc(256, sizeof (short)))) {
            fprintf(stderr, "smigen: Error allocating space for red.\n");
            return -1;
        };
        if (!(g = (short *) calloc(256, sizeof (short)))) {
            fprintf(stderr, "smigen: Error allocating space for green.\n");
            return -1;
        };
        if (!(b = (short *) calloc(256, sizeof (short)))) {
            fprintf(stderr, "smigen: Error allocating space for blue.\n");
            return -1;
        };

        if (strstr(input.palfile, ".pal") == NULL) {
            strcat(input.palfile, "/default.pal");
        }

        if (getlut_file(input.palfile, r, g, b)) {
            fprintf(stderr, "Error reading palette file %s\n", input.palfile);
            free(r);
            free(g);
            free(b);
            return -1;
        }
        for (i = 0; i < 256; i++) {
            input.palette[i * 3] = r[i];
            input.palette[i * 3 + 1] = g[i];
            input.palette[i * 3 + 2] = b[i];
        }
        free(r);
        free(g);
        free(b);
    }

    aminmax[0] = input.datamin;
    aminmax[1] = input.datamax;
    strcpy(atype, scale_type);

    si_used[0] = slope;
    si_used[1] = intercept;


    /* Set slope to 1 and intercept to 0 if float output */
    /* ------------------------------------------------- */
    if (strcmp(input.precision, "F") == 0) {
        strcpy(aopt, "No");
        si_used[0] = 1;
        si_used[1] = 0;
        strcpy(scale_type, "LINEAR");
    } else {
        strcpy(aopt, "Yes");
    }

    if (strcmp(input.resolution, "SMI") == 0) {
        dims_out[0] = 2048;
        dims_out[1] = 4096;
    } else if (strcmp(input.resolution, "SMI4") == 0) {
        dims_out[0] = 4096;
        dims_out[1] = 8192;
    } else if (strcmp(input.resolution, "LAND") == 0) {
        dims_out[0] = 4320;
        dims_out[1] = 8640;
    } else if (strcmp(input.resolution, "9km") == 0) {
        dims_out[0] = 2160;
        dims_out[1] = 4320;
    } else if (strcmp(input.resolution, "4km") == 0) {
        dims_out[0] = 4320;
        dims_out[1] = 8640;
    } else if (strcmp(input.resolution, "2km") == 0) {
        dims_out[0] = 4320 * 2;
        dims_out[1] = 8640 * 2;
    } else if (strcmp(input.resolution, "1km") == 0) {
        dims_out[0] = 4320 * 4;
        dims_out[1] = 8640 * 4;
    } else if (strcmp(input.resolution, "hkm") == 0) {
        dims_out[0] = 4320 * 8;
        dims_out[1] = 8640 * 8;
    } else if (strcmp(input.resolution, "qkm") == 0) {
        dims_out[0] = 4320 * 16;
        dims_out[1] = 8640 * 16;
    } else if (strcmp(input.resolution, "18km") == 0) {
        dims_out[0] = 1080;
        dims_out[1] = 2160;
    } else if (strcmp(input.resolution, "36km") == 0) {
        dims_out[0] = 540;
        dims_out[1] = 1080;
    } else if (strcmp(input.resolution, "90km") == 0) {
        dims_out[0] = 216;
        dims_out[1] = 432;
    } else if (strcmp(input.resolution, "thirddeg") == 0) {
        dims_out[0] = 180 * 3;
        dims_out[1] = 360 * 3;
    } else if (strcmp(input.resolution, "1deg") == 0) {
        dims_out[0] = 180;
        dims_out[1] = 360;
    } else if (strcmp(input.resolution, "hdeg") == 0) {
        dims_out[0] = 360;
        dims_out[1] = 720;
    } else if (strcmp(input.resolution, "qdeg") == 0) {
        dims_out[0] = 720;
        dims_out[1] = 1440;
    } else if (strcmp(input.resolution, "10deg") == 0) {
        dims_out[0] = 1800;
        dims_out[1] = 3600;
    } else if (strncmp(input.resolution, "udeg-", 5) == 0) {
        dims_out[0] = (int32) (180 / atof(&input.resolution[5]));
        dims_out[1] = 2 * dims_out[0];
    } else if (strncmp(input.resolution, "ukm-", 4) == 0) {
        dims_out[0] = (int32) (4320 * 4 / atof(&input.resolution[4]) + 0.5);
        dims_out[1] = 2 * dims_out[0];
    } else {
        printf("Improper resolution type: %s\n", input.resolution);
        exit(1);
    }
    printf("dims_out: %d x %d\n", dims_out[0], dims_out[1]);
    static Hdf::hdf5_bin input_binfile;

    // Check if HDF5 input file
    H5E_auto_t old_func;
    void *old_client_data;
    H5Eget_auto(H5E_DEFAULT, &old_func, &old_client_data);
    H5Eset_auto(H5E_DEFAULT, NULL, NULL); // Turn off error handling
    htri_t ishdf5 = H5Fis_hdf5(input.ifile);
    // Restore previous error handler
    H5Eset_auto(H5E_DEFAULT, old_func, old_client_data);

    int ncid, grpid, binindex_id, binlist_id, bindata_id, binqual_id;
    int nside;

    if (ishdf5 > 0) {
        status = nc_open(input.ifile, 0, &ncid);
        if (status != NC_NOERR) {
            status = open_input_hdf5(input.ifile, input.prod, &input_binfile,
                    &bin_dataset_idx);
            if (status == FAIL) exit(1);
            hdf5_input = 1;
        } else {
            char nam_buf[256];
            nam_buf[0] = 0;
            status = nc_get_att(ncid, NC_GLOBAL, "Mission", nam_buf);
            if (status != NC_NOERR)
                status = nc_get_att(ncid, NC_GLOBAL, "mission", nam_buf);
            if ((strcmp(nam_buf, "SAC-D Aquarius") == 0) ||
                    (strcmp(nam_buf, "SMAP") == 0)) {
                nc_close(ncid);
                status = open_input_hdf5(input.ifile, input.prod, &input_binfile,
                        &bin_dataset_idx);
                if (status == FAIL) exit(1);
                hdf5_input = 1;
            } else {
                status = open_input_ncdf4(ncid, input.prod, (size_t *) & nrows,
                        &grpid, &binindex_id, &binlist_id,
                        &bindata_id, &binqual_id, &meta_l3b);
                ncdf4_input = 1;
                if ((nrows % 2) == 1) {
                    healpix_input = 1;
                    printf("HEALPIX input\n");
                    nside = (nrows + 1) / 4;
                } else {
                    healpix_input = 0;
                }
            }
        }
    } else {
        status = open_input_hdf(input.ifile, input.prod, &fid, &sdfid,
                &vgid, vdata_id, &meta_l3b);

        if (status == FAIL) exit(1);

        sds_id = SDselect(sdfid, SDnametoindex(sdfid, input.prod));
        if (sds_id != -1) land_input = 1;
        else hdf4_input = 1;

        if (hdf4_input == 1 &&
                (strcmp(input.prod, "pixels") != 0) &&
                (strcmp(input.prod, "scenes") != 0)) {
            if (VSfind(fid, input.prod) == 0) {
                printf("Product \"%s\" not found.\n", input.prod);
                exit(-1);
            }
        }
    }

    printf("Input file      : %s\n", input.ifile);
    printf("Projection      : %s\n", input.projection);
    printf("Resolution      : %s\n", input.resolution);
    printf("Gap Fill        : %d\n", input.gap_fill);
    if (input.minobs != 0) printf("Min Obs      : %d\n", input.minobs);
    if (prod_num >= 0) printf("Parameter       : %s\n", parmname_list[prod_num]);
    if (prod_num == -2) printf("Parameter       : %s\n", input.prod);
    if (prod_num == -3) printf("Parameter       : %s\n", input.prod);
    printf("Measure         : %s\n", measure_list[input.meas - 1]);
    printf("Scale Type      : %s\n", scale_type);
    printf("Data Min (abs)  : %8.4f\n", input.datamin);
    printf("Data Max (abs)  : %8.4f\n", input.datamax);
    printf("Precision       : %s\n", input.precision);
    printf("Scale Slope     : %8.4f\n", si_used[0]);
    printf("Scale Intercept : %8.4f\n", si_used[1]);
    printf("Palette File    : %s\n", input.palfile);
    printf("Eastmost Long.  : %8.3f\n", input.loneast);
    printf("Westmost Long.  : %8.3f\n", input.lonwest);
    printf("Northmost Lat.  : %8.3f\n", input.latnorth);
    printf("Southmost Lat.  : %8.3f\n", input.latsouth);
    printf("Seam Longitude  : %8.3f\n\n", input.seam_lon);


    /* LAND INPUT */
    /* ========== */
    if (land_input) {

        /* SDS input */
        /* --------- */
        status = SDgetinfo(sds_id, buf, &rank, dims_in, &type, &nattrs);

        i = SDfindattr(sdfid, "Pixel size (meters)");
        if (i != -1) SDreadattr(sdfid, i, (VOIDP) & pix_sz);
        /* printf("pixsz: %d %f\n", i,pix_sz);*/

        /* 6371007.181 Earth Radius (m) */
        dims_in[0] =
                (int32) ((180 * 6371007.181) / (57.29577951308232087684 * pix_sz) + 0.5);
        dims_in[1] = 2 * dims_in[0];
        /*printf("input dimensions (full map): %d %d\n", dims_in[0], dims_in[1]);*/

        sub_scan[0] = 0;
        sub_pixl[0] = 0;
        sub_scan[1] = dims_in[0] - 1;
        sub_pixl[1] = dims_in[1] - 1;
        i = SDfindattr(sdfid, "Subset Scan Range");
        if (i != -1) SDreadattr(sdfid, i, (VOIDP) sub_scan);
        i = SDfindattr(sdfid, "Subset Pixel Range");
        if (i != -1) SDreadattr(sdfid, i, (VOIDP) sub_pixl);
        /* printf("sub_scan (r/c): %d %d  sub_pixl (r/c): %d %d\n", 
               sub_scan[0], sub_scan[1], sub_pixl[0], sub_pixl[1]);*/

    } else if (hdf4_input == 1) {

        /* BIN INPUT */
        /* ========= */

        /* Vdata (bin file) input */
        /* ---------------------- */
        dims_in[0] = nrows;
        dims_in[1] = 2 * nrows;
        printf("input dimensions (full map): %d %d\n", dims_in[0], dims_in[1]);

        binlist_buf = (uint8 *) calloc(dims_in[1], 12);
        sum_buf = (float32 *) calloc(dims_in[1], 2 * sizeof (float32));
        qual_buf = (uint8 *) calloc(dims_in[1], sizeof (uint8));

        strcpy(buf, input.prod);
        strcat(buf, "_sum,");
        strcat(buf, input.prod);
        strcat(buf, "_sum_sq");
        status = VSsetfields(vdata_id[0], "bin_num,nobs,nscenes,weights");
        status = VSsetfields(vdata_id[1], "start_num,begin,extent,max");
        status = VSsetfields(vdata_id[2], buf);
        if (vdata_id[3] != -1)
            status = VSsetfields(vdata_id[3], "qual_l3");
    } else if (hdf5_input == 1) {
        // HDF5 BinList
        nrows = input_binfile.nrows;
        dims_in[0] = nrows;
        dims_in[1] = 2 * nrows;
        printf("input dimensions (full map): %d %d\n", dims_in[0], dims_in[1]);

        binlist_buf = (uint8 *) calloc(dims_in[1], 12);
        sum_buf = (float32 *) calloc(dims_in[1], 2 * sizeof (float32));

        binlist_tid = H5Tcreate(H5T_COMPOUND, 12);
        H5Tinsert(binlist_tid, "bin_num", 0, H5T_STD_U32LE);
        H5Tinsert(binlist_tid, "nobs", 4, H5T_NATIVE_USHORT);
        H5Tinsert(binlist_tid, "nscenes", 6, H5T_NATIVE_USHORT);
        H5Tinsert(binlist_tid, "weights", 8, H5T_NATIVE_FLOAT);

        // HDF5 Binindex
        binindex_tid = H5Tcreate(H5T_COMPOUND, 16);
        H5Tinsert(binindex_tid, "start_num", 0, H5T_STD_I32LE);
        H5Tinsert(binindex_tid, "begin", 4, H5T_STD_I32LE);
        H5Tinsert(binindex_tid, "extent", 8, H5T_STD_I32LE);
        H5Tinsert(binindex_tid, "max", 12, H5T_STD_I32LE);

        // HDF5 BinProduct
        bindata_tid = H5Tcreate(H5T_COMPOUND, 2 * sizeof (H5T_NATIVE_FLOAT));

        strcpy(buf, input.prod);
        strcat(buf, "_sum");
        H5Tinsert(bindata_tid, buf, 0, H5T_NATIVE_FLOAT);

        strcat(buf, "_sq");
        H5Tinsert(bindata_tid, buf, sizeof (H5T_NATIVE_FLOAT), H5T_NATIVE_FLOAT);

        vdata_id[3] = -1;
    } else {
        // NETCDF4
        dims_in[0] = nrows;
        dims_in[1] = 2 * nrows;

        binlist_buf = (uint8 *) calloc(dims_in[1], 16);
        sum_buf = (float32 *) calloc(dims_in[1], 2 * sizeof (float32));

        vdata_id[3] = -1;
        if (binqual_id != -1) {
            qual_buf = (uint8 *) calloc(dims_in[1], sizeof (uint8));
            vdata_id[3] = 0;
        }
    }

    int binlist_sz;
    if (ncdf4_input)
        binlist_sz = 16;
    else
        binlist_sz = 12;

    /* Adjust for geo-subset with bin file input */
    /* ----------------------------------------- */
    if (hdf4_input || hdf5_input) {
        f32 = input.latnorth;
        input.latnorth = -input.latsouth;
        input.latsouth = -f32;
    }

    get_minmax_rowcol(&input.latnorth, &input.lonwest, input.projection,
            dims_out[0],
            &min_out_row, &max_out_row, &min_out_col, &max_out_col);

    printf("min/max out row/col: %d %d %d %d\n",
            min_out_row, max_out_row, min_out_col, max_out_col);

    /*
    if ((max_out_row - min_out_row) <= 0 || (max_out_col - min_out_col) <= 0) {
      printf("Improper boundary limits: East: %f  West: %f North: %f South: %f\n", 
             input.loneast, input.lonwest, -input.latsouth, -input.latnorth);
      exit(1);
    }
     */


    if (hdf4_input || hdf5_input) {
        f32 = input.latnorth;
        input.latnorth = -input.latsouth;
        input.latsouth = -f32;
    }

    if (land_input && max_out_col < min_out_col) {
        printf("max_out_col must be greater than min_out_col for land input.\n");
        exit(1);
    }

    dims_out_sub[0] = max_out_row - min_out_row + 1;
    if (max_out_col > min_out_col) {
        dims_out_sub[1] = max_out_col - min_out_col + 1;
    } else {
        dims_out_sub[1] = dims_out[1] - (min_out_col - max_out_col + 1);
    }

    if (strcmp(input.precision, "F") == 0)
        par_byt = (uint8 *) calloc(dims_out_sub[0] * dims_out_sub[1],
            4 * sizeof (uint8));
    else if (strcmp(input.precision, "I") == 0)
        par_byt = (uint8 *) calloc(dims_out_sub[0] * dims_out_sub[1],
            2 * sizeof (uint8));
    else
        par_byt = (uint8 *) calloc(dims_out_sub[0] * dims_out_sub[1],
            sizeof (uint8));

    if (vdata_id[3] != -1) {
        qual_byt = (uint8 *) calloc(dims_out_sub[0] * dims_out_sub[1],
                sizeof (uint8));
        memset(qual_byt, 255, dims_out_sub[0] * dims_out_sub[1]);
    }


    /* Setup Arrays */
    /* ------------ */
    out_sum = (int32 *) calloc(dims_out[1], sizeof (int32));
    out_num = (int32 *) calloc(dims_out[1], sizeof (int32));
    best_qual = (uint8 *) calloc(dims_out[1], sizeof (uint8));

    if (input.seam_lon != -180) {
        out_sum_rot = (int32 *) calloc(dims_out[1], sizeof (int32));
        out_num_rot = (int32 *) calloc(dims_out[1], sizeof (int32));
        best_qual_rot = (uint8 *) calloc(dims_out[1], sizeof (uint8));
    }

    ratio = (float32) dims_out[0] / dims_in[0];
    max_rowbuf = (int32) (1 / ratio + 1);
    if (healpix_input) max_rowbuf *= 2;

    /* input_databuf (sinusoidal rojection) */
    input_databuf = (int32 *) calloc(max_rowbuf * dims_in[1], sizeof (int32));
    input_databuf_16 = (int16 *) calloc(max_rowbuf * dims_in[1], sizeof (int16));
    input_qualbuf = (uint8 *) calloc(max_rowbuf * dims_in[1], sizeof (uint8));

    in2out = (int32 *) calloc(dims_in[0], sizeof (int32));
    for (i = 0; i < dims_in[0]; i++) in2out[i] = (int32) ((i + 0.5) * ratio);

    start[1] = 0;
    edges[1] = sub_pixl[1] - sub_pixl[0] + 1;

    prev_input_row = 0;
    input_row = 0;
    int32_t nbinsread = 0;

    if (land_input) scl_fac = 10000;
    if (hdf4_input || hdf5_input || ncdf4_input) scl_fac = 1000000;
    if (strcmp(input.precision, "F") == 0 && hdf4_input) scl_fac = 100000;
    if (strcmp(input.precision, "F") == 0 && hdf5_input) scl_fac = 100000;

    /* For each scan ... (Main Loop) (NORTH to SOUTH) */
    /* ---------------------------------------------- */
    for (out_row = 0; out_row < dims_out[0]; out_row++) {

        if ((out_row % 500) == 0) {
            time(&tnow);
            tmnow = localtime(&tnow);
            printf("out_row:%6d %s", out_row, asctime(tmnow));
        }

        if (healpix_input) {
            float cs = cos((out_row + 1)*(180.0 / dims_out[0])*(PI / 180));
            if (cs > (2. / 3))
                input_row = nside * sqrt(3 * (1 - cs));
            else if (cs >= -(2. / 3))
                input_row = nside * (2 - 1.5 * cs);
            else
                input_row = nside * (4 - sqrt(3 * (1 + cs)));
            for (i = prev_input_row; i < input_row; i++)
                input_row_array[i - prev_input_row] = i;
        } else {
            i = 0;
            while (in2out[input_row] == out_row) {
                input_row_array[i++] = input_row;
                input_row++;
                if (i == MAX_NUM_INPUTROW) {
                    printf("i > MAX_NUM_INPUTROW\n");
                    exit(-1);
                }
            }
        }

        start[0] = prev_input_row - sub_scan[0];
        edges[0] = input_row - prev_input_row;
        /*    printf("prev: %d  sub_scan: %d\n", prev_input_row, sub_scan[0]);*/

        /* Fill Input Data Buffer */
        /* ---------------------- */
        if (land_input) {

            /* (Land) Map File Input */
            /* --------------------- */
            if (start[0] >= 0 && start[0] <= sub_scan[1] - sub_scan[0]) {
                status = SDreaddata(sds_id, start, NULL, edges,
                        (VOIDP) & input_databuf_16[sub_pixl[0]]);

                for (i = 0; i < edges[0] * edges[1]; i++)
                    input_databuf[sub_pixl[0] + i] = input_databuf_16[sub_pixl[0] + i];
            }
        } else {

            /* Bin File Input */
            /* -------------- */

            /* vdata 0: BinList  (bin_num, nobs, nscenes, weights)
               vdata 1: BinIndex (start_num, begin, extent, max)
               vdata 2: Data     (sum, sum_sq)

               where start_num is the starting bin # for each row (1-based)
                     begin is the starting bin # where data actually exist 
                     (0 if row empty)
                     extent is the number of bins in row for which data exists
                     max is the total number of bins in row

             */

            /* Initialize input_databuf */
            /* ------------------------ */
            for (i = 0; i < max_rowbuf * dims_in[1]; i++) input_databuf[i] = -32767;

            /* For each line in databuf ... */
            /* ---------------------------- */
            for (i = 0; i < edges[0]; i++) {

                /* Get bin row info */
                /* ---------------- */
                if (hdf4_input) {
                    VSseek(vdata_id[1], start[0] + i);
                    nread = VSread(vdata_id[1], (uint8 *) binindex_buf,
                            1, FULL_INTERLACE);

                    if (nread == -1) {
                        printf("Problem with binindex read: %d\n", out_row);
                        exit(1);
                    }
                } else if (hdf5_input) {
                    // HDF5 input
                    hid_t dataset = input_binfile.get_index_table();
                    hid_t filespace = H5Dget_space(dataset);
                    hsize_t strt = start[0] + i;
                    status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, &strt, NULL,
                            &one, NULL);
                    hid_t dataspace = H5Screate_simple(1, &one, NULL);
                    status = H5Dread(dataset, binindex_tid, dataspace,
                            filespace, H5P_DEFAULT, binindex_buf);

                    H5Sclose(dataspace);
                    H5Sclose(filespace);
                } else {
                    size_t indexp = start[0] + i;
                    status = nc_get_var1(grpid, binindex_id, &indexp, binindex_buf);
                }

                /* If row has data ... */
                /* ------------------- */
                if (binindex_buf[1] != 0) {

                    /* error check */
                    if (binindex_buf[2] > dims_in[1]) {
                        printf("Bin Index (%d) greater than column dimension (%d) for out row: %d\n",
                                binindex_buf[2], dims_in[1], out_row);
                        printf("%d %d\n", out_row, start[0] + i);
                        exit(1);
                    }

                    /* Get bin #, # of scenes and weights for each filled bin in row */
                    /* ------------------------------------------------------------- */
                    if (hdf4_input) {
                        nread = VSread(vdata_id[0], binlist_buf, binindex_buf[2],
                                FULL_INTERLACE);
                        if (nread == -1) {
                            printf("Problem with binlist read: %d\n", out_row);
                            exit(1);
                        }

                        /* Get data values (sum, sum_sq) for each filled bin in row */
                        /* -------------------------------------------------------- */
                        if (strcmp(input.prod, "pixels") != 0 &&
                                strcmp(input.prod, "scenes") != 0) {
                            nread = VSread(vdata_id[2], (uint8 *) sum_buf, binindex_buf[2],
                                    FULL_INTERLACE);

                            if (nread == -1) {
                                printf("Problem with sum/sum_sqr read: %d\n", out_row);
                                exit(1);
                            }

                            if (vdata_id[3] != -1) {
                                nread = VSread(vdata_id[3], (uint8 *) qual_buf,
                                        binindex_buf[2], FULL_INTERLACE);

                                if (nread == -1) {
                                    printf("Problem with qual_l3 read: %d\n", out_row);
                                    exit(1);
                                }
                            }
                        }
                    } else if (hdf5_input) {
                        // HDF5 input

                        hsize_t strt = nbinsread;
                        hsize_t count = binindex_buf[2];

                        // Get bin #, # of scenes and weights for each filled bin in row
                        hid_t dataset = input_binfile.get_list_table();
                        hid_t filespace = H5Dget_space(dataset);

                        status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, &strt,
                                NULL, &count, NULL);
                        hid_t dataspace = H5Screate_simple(1, &count, NULL);
                        status = H5Dread(dataset, binlist_tid, dataspace,
                                filespace, H5P_DEFAULT, binlist_buf);

                        H5Sclose(dataspace);
                        H5Sclose(filespace);

                        // Get data values (sum, sum_sq) for each filled bin in row
                        if (strcmp(input.prod, "pixels") != 0 &&
                                strcmp(input.prod, "scenes") != 0) {

                            hid_t dataset = input_binfile.get_data_table(bin_dataset_idx);

                            hid_t filespace = H5Dget_space(dataset);

                            status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, &strt,
                                    NULL, &count, NULL);
                            hid_t dataspace = H5Screate_simple(1, &count, NULL);
                            status = H5Dread(dataset, bindata_tid, dataspace,
                                    filespace, H5P_DEFAULT, sum_buf);

                            H5Sclose(dataspace);
                            H5Sclose(filespace);
                        }

                        nbinsread += binindex_buf[2];
                    } else {
                        // CDF4 input
                        size_t strt = nbinsread;
                        size_t count = binindex_buf[2];

                        status = nc_get_vara(grpid, binlist_id, &strt, &count,
                                binlist_buf);

                        if (binqual_id != -1) {
                            status = nc_get_vara(grpid, binqual_id, &strt, &count,
                                    qual_buf);
                        }

                        if (strcmp(input.prod, "pixels") != 0 &&
                                strcmp(input.prod, "scenes") != 0) {
                            status = nc_get_vara(grpid, bindata_id, &strt, &count, sum_buf);
                        }

                        nbinsread += binindex_buf[2];
                    }

                    float heal2sin;
                    if (healpix_input) {
                        cos_theta_input =
                                cos((((input_row_array[i] + 0.5) / dims_in[0]) - 0.5) * PI);
                        heal2sin = (dims_in[1] * cos_theta_input) / binindex_buf[3];
                    } else {
                        offset = (dims_in[1] - binindex_buf[3]) / 2;
                    }


                    /* Populate input databuf with bin data */
                    /* ------------------------------------ */
                    for (k = 0; k < binindex_buf[2]; k++) {
                        memcpy(&bin_num, &binlist_buf[k * binlist_sz], 4);
                        memcpy(&nobs, &binlist_buf[k * binlist_sz + 4], 2);
                        memcpy(&nscenes, &binlist_buf[k * binlist_sz + 6], 2);
                        memcpy(&wgt, &binlist_buf[k * binlist_sz + 8], 4);

                        // Skip if nan
                        if (input.meas == 1 || input.meas == 2) {
                            if (isnan(sum_buf[k * 2])) {
                                continue;
                            }
                        }

                        /* Check if less than minobs */
                        if (input.minobs != 0)
                            if ((uint32_t) nobs < input.minobs) continue;

                        if (healpix_input) {
                            int32_t rot_bin_num = ((bin_num - binindex_buf[0]) +
                                    binindex_buf[3] / 2) % binindex_buf[3];
                            // j = i*dims_in[1] + (bin_num - binindex_buf[0]) * heal2sin
                            //+ (dims_in[1] - binindex_buf[3]*heal2sin) / 2;

                            j = i * dims_in[1] + rot_bin_num * heal2sin
                                    + (dims_in[1] - binindex_buf[3] * heal2sin) / 2;

                        } else {
                            j = i * dims_in[1] + (bin_num - binindex_buf[0]) + offset;
                        }
                        //	    printf("%d %d %f\n", out_row, j,  (bin_num - binindex_buf[0]) * heal2sin);
                        if (j >= (max_rowbuf * dims_in[1])) {
                            printf("Databuf element (%d) greater than buffer size (%d) for out_row: %d\n",
                                    j, max_rowbuf * dims_in[1], out_row);
                            printf("i: %d  k: %d  bin_num: %d  offset: %d\n",
                                    i, (int) k, bin_num, offset);
                            printf("binindex_buf[0]: %d  binindex_buf[1]: %d\n",
                                    binindex_buf[0], binindex_buf[1]);
                            printf("binindex_buf[2]: %d  binindex_buf[3]: %d\n",
                                    binindex_buf[2], binindex_buf[3]);
                            exit(1);
                        }

                        if (j < 0) {
                            printf("Databuf element (%d) is negative for out_row: %d\n",
                                    j, out_row);
                            printf("i: %d  k: %d  bin_num: %d  offset: %d  bin_row: %d\n",
                                    i, (int) k, bin_num, offset, start[0] + i);
                            printf("bin start number for row: %d  data starting bin #: %d\n",
                                    binindex_buf[0], binindex_buf[1]);
                            printf("# of bins in row with data: %d  # of bins in row: %d\n",
                                    binindex_buf[2], binindex_buf[3]);
                            exit(1);
                        }

                        if (vdata_id[3] != -1) input_qualbuf[j] = qual_buf[k];

                        if (input.meas == 1) {

                            flt_val = sum_buf[k * 2] / wgt;
                            input_databuf[j] = (int32) (scl_fac * flt_val + 0.5);

                        } else if ((input.meas == 2 || input.meas == 3) &&
                                (wgt * wgt > nscenes)) {

                            flt_val = sum_buf[k * 2] / wgt;
                            flt_val = (sum_buf[k * 2 + 1] / wgt) - (flt_val * flt_val);

                            flt_val = flt_val * wgt * wgt / (wgt * wgt - nscenes);

                            if (input.meas == 3) flt_val = sqrt(flt_val);

                            input_databuf[j] = (int32) (scl_fac * flt_val + 0.5);

                        } else if (input.meas == 4) {

                            input_databuf[j] = (int32) ((1.0 * nobs) / nscenes + 0.5);

                        } else if (input.meas == 5) {

                            input_databuf[j] = (int32) (nscenes);
                        }

                    } /* for (k=0; */



                    /* Shift bin rows with odd number of pixels by "half" pixel */
                    /* -------------------------------------------------------- */
                    if (binindex_buf[3] % 2 == 1 &&
                            input.meas <= 3 &&
                            getenv("SMIGENNOHALF") == NULL) {

                        /* SMIGENNOHALF only set for debugging */

                        for (k = dims_in[1] - offset - 2; k >= offset; k--) {

                            j = i * dims_in[1] + k;
                            if (k == offset) jm1 = j + binindex_buf[3] - 1;
                            else jm1 = j - 1;

                            if (input_databuf[jm1] != -32767) {
                                if (input_databuf[j] != -32767) {
                                    /*
                                    input_databuf[j] += input_databuf[jm1];
                                    input_databuf[j] /= 2;
                                     */
                                }
                            }

                        } /* k loop */

                        /* Last pixel equal to first */
                        input_databuf[(i + 1) * dims_in[1] - offset - 1] = input_databuf[j];

                    } /* if (binindex_buf[3] ... */

                    /* Fill in gaps in input buffer */
                    /* ---------------------------- */
                    /* 1 pixel gaps */
                    if (input.gap_fill >= 1 && input.meas <= 3) {
                        for (k = 1; k < max_rowbuf * dims_in[1] - 1; k++) {
                            if (input_databuf[k - 1] != -32767 &&
                                    input_databuf[k] == -32767 &&
                                    input_databuf[k + 1] != -32767) {
                                if ((vdata_id[3] != -1)) {
                                    if (input_qualbuf[k - 1] == input_qualbuf[k + 1]) {
                                        /* same quality, avg data and copy qual to filled pixel */
                                        input_databuf[k] = (int32)
                                                (0.5 * (input_databuf[k - 1] + input_databuf[k + 1]) + 0.5);
                                        input_qualbuf[k] = input_qualbuf[k - 1];
                                    } else {
                                        /* use whichever value has the better (lower) quality */
                                        tk = (input_databuf[k - 1] < input_databuf[k + 1]) ? k - 1 : k + 1;
                                        input_qualbuf[k] = input_qualbuf[tk];
                                        input_databuf[k] = input_databuf[tk];
                                    }
                                } else {
                                    input_databuf[k] = (int32)
                                            (0.5 * (input_databuf[k - 1] + input_databuf[k + 1]) + 0.5);
                                }
                            }
                        }
                    }

                    /* 2 pixel gaps */
                    if (input.gap_fill >= 2 && input.meas <= 3) {
                        for (k = 1; k < max_rowbuf * dims_in[1] - 2; k++) {
                            if (input_databuf[k - 1] != -32767 &&
                                    input_databuf[k] == -32767 &&
                                    input_databuf[k + 1] == -32767 &&
                                    input_databuf[k + 2] != -32767) {
                                if ((vdata_id[3] != -1)) {
                                    if (input_qualbuf[k - 1] == input_qualbuf[k + 2]) {
                                        /* same quality, average data and copy quality to filled pixel */
                                        input_databuf[k] = (int32)
                                                ((input_databuf[k + 2] + 2 * input_databuf[k - 1]) / 3 + 0.5);
                                        input_databuf[k + 1] = (int32)
                                                ((2 * input_databuf[k + 2] + input_databuf[k - 1]) / 3 + 0.5);
                                        input_qualbuf[k] = input_qualbuf[k - 1];
                                        input_qualbuf[k + 1] = input_qualbuf[k - 1];
                                    } else {
                                        /* use whichever value has the better (lower) quality */
                                        tk = (input_databuf[k - 1] < input_databuf[k + 2]) ? k - 1 : k + 2;
                                        input_qualbuf[k] = input_qualbuf[tk];
                                        input_databuf[k] = input_databuf[tk];
                                        input_qualbuf[k + 1] = input_qualbuf[tk];
                                        input_databuf[k + 1] = input_databuf[tk];
                                    }
                                } else {
                                    input_databuf[k] = (int32)
                                            ((input_databuf[k + 2] + 2 * input_databuf[k - 1]) / 3 + 0.5);
                                    input_databuf[k + 1] = (int32)
                                            ((2 * input_databuf[k + 2] + input_databuf[k - 1]) / 3 + 0.5);
                                }
                            }
                        }
                    }

                } /* if binindex_buf ... */
            } /* for (i=0; */

        } /* if (land_input) */

        prev_input_row = input_row;

        if (out_row < min_out_row) continue;
        if (out_row > max_out_row) break;

        for (i = 0; i < dims_out[1]; i++) out_sum[i] = -32767;
        for (i = 0; i < dims_out[1]; i++) out_num[i] = 0;
        for (i = 0; i < dims_out[1]; i++) best_qual[i] = 255;

        cos_theta_output = cos((((out_row + 0.5) / dims_out[0]) - 0.5) * PI);

        /* Loop over rows contributing to output row */
        /* ----------------------------------------- */
        for (j = 0; j < edges[0]; j++) {

            cos_theta_input =
                    cos((((input_row_array[j] + 0.5) / dims_in[0]) - 0.5) * PI);

            if (strcmp(input.projection, "SIN") == 0) {
                cos_fac = cos_theta_output / cos_theta_input;
                pixmin = (int32) (0.5 * dims_out[1] * (1 - cos_theta_output));
                pixmax = (int32) (0.5 * dims_out[1] * (1 + cos_theta_output));
            }

            if (strcmp(input.projection, "RECT") == 0) {
                cos_fac = 1 / cos_theta_input;
                pixmin = 0;
                pixmax = dims_out[1];
            }

            if (pixmin < 0) pixmin = 0;
            if (pixmax > dims_out[1]) pixmax = dims_out[1];

            /* For each output pixel accumulate from input pixels */
            /* -------------------------------------------------- */
            for (i = pixmin; i < pixmax; i++) {

                input_col =
                        (int32) ((((float32) (i + 0.5) / dims_out[1]) -
                        0.5 * (1 - cos_fac)) * (dims_in[1] / cos_fac));

                input_value = input_databuf[j * dims_in[1] + input_col];

                if (input_value != -32767) {
                    //	  printf("%d %d %d %d\n", out_row, j, i, input_col);
                    if ((vdata_id[3] != -1)) {
                        input_qual = input_qualbuf[j * dims_in[1] + input_col];
                        if (input_qual < best_qual[i]) {
                            out_num[i] = 0;
                            out_sum[i] = -32767;
                            best_qual[i] = input_qual;
                        } else if (input_qual > best_qual[i]) continue;
                    }

                    out_num[i]++;

                    if (out_sum[i] != -32767)
                        out_sum[i] += input_value;
                    else
                        out_sum[i] = input_value;
                } /* if ... */
            } /* i loop */
        } /* j loop */


        /* Rotate if SEAM_LON != -180 */
        /* -------------------------- */
        if (input.seam_lon != -180 && land_input == 0) {
            rotate = (int32) ((input.seam_lon + 180) * (pixmax - pixmin) / 360 + 0.5);

            for (i = pixmin; i < pixmax; i++) {
                out_sum_rot[i] = out_sum[(i + rotate) % (pixmax - pixmin)];
                out_num_rot[i] = out_num[(i + rotate) % (pixmax - pixmin)];
                best_qual_rot[i] = best_qual[(i + rotate) % (pixmax - pixmin)];
            }

            for (i = pixmin; i < pixmax; i++) {
                out_sum[i] = out_sum_rot[i];
                out_num[i] = out_num_rot[i];
                best_qual[i] = best_qual_rot[i];
            }
        }


        /* Initialize output byte arrays with fill values */
        /* --------------------------------------------- */
        if (strncmp(input.prod, "refl", 4) == 0 ||
                strncmp(input.prod, "rhos", 4) == 0) {
            fill_val = 0;
        } else {
            fill_val = 255;
        }

        for (i = 0; i < dims_out[1]; i++) {

            if (max_out_col > min_out_col) {
                if (i < min_out_col || i > max_out_col) continue;
            } else {
                if (i < min_out_col && i > max_out_col) continue;
            }

            if (max_out_col > min_out_col) {
                if (land_input || healpix_input)
                    k = (out_row - min_out_row) * dims_out_sub[1] + (i - min_out_col);
                else
                    k = (max_out_row - out_row) * dims_out_sub[1] + (i - min_out_col);
            } else {
                if (land_input || healpix_input)
                    k = (out_row - min_out_row) * dims_out_sub[1] +
                        (i - min_out_col) * (i >= min_out_col) +
                    (dims_out[1] - min_out_col + i - 1) * (i < max_out_col);
                else
                    k = (max_out_row - out_row) * dims_out_sub[1] +
                        (i - min_out_col) * (i >= min_out_col) +
                    (dims_out[1] - min_out_col + i - 1) * (i < max_out_col);
            }

            if (k < 0 || k >= (dims_out_sub[0] * dims_out_sub[1])) {
                printf("k: %d i: %d\n", (int) k, i);
                exit(1);
            }

            if (strcmp(input.precision, "F") == 0) {
                memcpy(&par_byt[4 * k], &fill_val_float, 4);
            } else if (strcmp(input.precision, "I") == 0) {
                memcpy(&par_byt[2 * k], &fill_val, 1);
                memcpy(&par_byt[2 * k + 1], &fill_val, 1);
            } else if (strcmp(input.precision, "B") == 0) {
                par_byt[k] = fill_val;
            }
        }

        /* Scale each output pixel */
        /* ----------------------- */
        for (i = 0; i < dims_out[1]; i++) {

            if (max_out_col > min_out_col) {
                if (i < min_out_col || i > max_out_col) continue;
            } else {
                if (i < min_out_col && i > max_out_col) continue;
            }

            if (max_out_col > min_out_col) {
                if (land_input || healpix_input)
                    k = (out_row - min_out_row) * dims_out_sub[1] + (i - min_out_col);
                else
                    k = (max_out_row - out_row) * dims_out_sub[1] + (i - min_out_col);
            } else {
                if (land_input || healpix_input)
                    k = (out_row - min_out_row) * dims_out_sub[1] +
                        (i - min_out_col) * (i >= min_out_col) +
                    (dims_out[1] - min_out_col + i - 1) * (i < max_out_col);
                else
                    k = (max_out_row - out_row) * dims_out_sub[1] +
                        (i - min_out_col) * (i >= min_out_col) +
                    (dims_out[1] - min_out_col + i - 1) * (i < max_out_col);
            }

            if (out_sum[i] != -32767) {

                filled_data_bins++;

                flt_val = (out_sum[i] / out_num[i]) / scl_fac;

                if (strcmp(input.precision, "F") == 0) {
                    if (flt_val < dminmax[0]) dminmax[0] = flt_val;
                    if (flt_val > dminmax[1]) dminmax[1] = flt_val;
                    memcpy(&par_byt[4 * k], &flt_val, 4);
                    if (vdata_id[3] != -1) qual_byt[k] = best_qual[i];
                    continue;
                }

                if (flt_val < dminmax[0]) dminmax[0] = flt_val;
                if (flt_val > dminmax[1]) dminmax[1] = flt_val;

                if (flt_val < input.datamin) flt_val = input.datamin;
                if (flt_val > input.datamax) flt_val = input.datamax;

                if (prod_num == -2) {

                    if (out_sum[i] > 0.0) {
                        flt_val = 396.3 * log10(8.546 * flt_val + 1);
                    }

                } else if (prod_num == -3) {

                    if (flt_val > 0.01)
                        flt_val = rint((log10(flt_val) + 2.0) / (2.0 / 255));
                    else
                        flt_val = 0;

                    if (flt_val < 0) flt_val = 0;
                    if (flt_val > 255) flt_val = 255;

                    /*
                    if (out_sum[i] > 0.0) {
                      flt_val = 400.0*log10(50.0*flt_val+1);
                     */

                } else if (strcmp(scale_type, "LINEAR") == 0) {
                    flt_val = (flt_val - intercept) / slope;
                } else if (strcmp(scale_type, "LOG") == 0) {
                    flt_val = (log10(flt_val) - intercept) / slope;
                }

                if (strcmp(input.precision, "I") == 0) {
                    if (strcmp(input.prod, "pixels") == 0 ||
                            strcmp(input.prod, "scenes") == 0) {
                        scale_val_16b = (uint16) out_sum[i];
                        if (scale_val_16b < dminmax[0]) dminmax[0] = scale_val_16b;
                        if (scale_val_16b > dminmax[1]) dminmax[1] = scale_val_16b;
                    } else {
                        scale_val_16b = (uint16) (flt_val + 0.5);
                    }
                    memcpy(&par_byt[2 * k], &scale_val_16b, 2);
                } else {
                    par_byt[k] = (uint8) (flt_val + 0.5);
                }

                /* Fill quality output buf */
                /* ----------------------- */
                if (vdata_id[3] != -1) qual_byt[k] = best_qual[i];

            }
        } /* i loop */

    } /* out_row loop */

    if (land_input) SDendaccess(sds_id);


    /* LAND MASK */
    /* --------- */
    if (land_input && strcmp(input.projection, "RECT") == 0) {

        tmp_str = getenv("OCDATAROOT");
        if (tmp_str == 0x0) {
            printf("Environment variable OCDATAROOT not defined.\n");
            exit(1);
        }
        strcpy(buf, tmp_str);
        strcat(buf, "/common/landmask.dat");

        printf("Opening: %s for masking\n", buf);
        mask_fp = fopen(buf, "rb");

        if (mask_fp == 0x0) {
            printf("Land mask file: %s not found.\n", buf);
            exit(1);
        }

        ptr_arr = (uint16 *) calloc(1024 * 64, sizeof (uint16));
        for (i = 0; i < 8 * 128; i++) mixed_buf[i] = NULL;

        fseek(mask_fp, 1024 * 2, SEEK_SET);
        fread(ptr_arr, 2, 1024 * 64, mask_fp);

        uname(&u_name);
        //    if (strcmp(u_name.sysname, "Linux") == 0) {
        if (endianess() == 1) {
            ptr_i8 = (uint8 *) ptr_arr;
            for (i = 0; i < 1024 * 64; i++) {
                memcpy(&i8, ptr_i8, 1);
                memcpy(ptr_i8, ptr_i8 + 1, 1);
                memcpy(ptr_i8 + 1, &i8, 1);
                ptr_i8 += 2;
            }
        }

        /* (NORTH to SOUTH) */
        /* ---------------- */
        for (out_row = 0; out_row < dims_out[0]; out_row++) {

            if ((out_row % 500) == 0) {
                time(&tnow);
                tmnow = localtime(&tnow);
                /*	printf("out_row:%6d %s", out_row, asctime(tmnow));*/
            }

            if (out_row < min_out_row) continue;
            if (out_row > max_out_row) break;

            if (strcmp(input.projection, "RECT") == 0) {
                pixmin = 0;
                pixmax = dims_out[1];
            }

            if (pixmin < 0) pixmin = 0;
            if (pixmax > dims_out[1]) pixmax = dims_out[1];

            lat = 180 - (int32) (((out_row + 0.5) / dims_out[0]) * 180) - 1;
            latf = lat - (180 - (((out_row + 0.5) / dims_out[0]) * 180) - 1);

            /* For each output pixel accumulate from input pixels */
            /* -------------------------------------------------- */
            for (i = pixmin; i < pixmax; i++) {

                if (max_out_col > min_out_col) {
                    if (i < min_out_col || i > max_out_col) continue;
                } else {
                    if (i < min_out_col && i > max_out_col) continue;
                }

                lon = (int32) (((i + 0.5) / dims_out[1]) * 360);
                lonf = (((i + 0.5) / dims_out[1]) * 360) - lon;

                cur_ptr = ptr_arr[360 * lat + lon];

                if (cur_ptr == 0) {
                    last_ptr = 0;

                    if (max_out_col > min_out_col) {
                        k = (out_row - min_out_row) * dims_out_sub[1] + (i - min_out_col);
                    } else {
                        k = (out_row - min_out_row) * dims_out_sub[1] +
                                (i - min_out_col) * (i >= min_out_col) +
                                (dims_out[1] - min_out_col + i - 1) * (i < max_out_col);
                    }

                    if (strcmp(input.precision, "F") == 0) {
                        memcpy(&par_byt[4 * k], &fill_val_float, 4);
                    } else if (strcmp(input.precision, "I") == 0) {
                        memcpy(&par_byt[2 * k], &fill_val, 1);
                        memcpy(&par_byt[2 * k + 1], &fill_val, 1);
                    } else if (strcmp(input.precision, "B") == 0) {
                        par_byt[k] = fill_val;
                    }
                    continue;
                }

                if (cur_ptr == 1) {
                    last_ptr = 1;
                    continue;
                }

                if (cur_ptr > 65 && cur_ptr != last_ptr) {
                    last_ptr = cur_ptr;
                    if (mixed_buf[cur_ptr - 66] == NULL) {

                        fseek(mask_fp, cur_ptr * 1024 * 2, SEEK_SET);
                        fread(&mixed[0], 2, 1024, mask_fp);

                        //	    if (strcmp(u_name.sysname, "Linux") == 0) {
                        if (endianess() == 1) {
                            ptr_i8 = (uint8 *) & mixed[0];
                            for (k = 0; k < 1024; k++) {
                                memcpy(&i8, ptr_i8, 1);
                                memcpy(ptr_i8, ptr_i8 + 1, 1);
                                memcpy(ptr_i8 + 1, &i8, 1);
                                ptr_i8 += 2;
                            }
                        }

                        mixed_buf[cur_ptr - 66] = (int16 *) calloc(8 * 128, 2);
                        memcpy(mixed_buf[cur_ptr - 66], &mixed[0], 8 * 128 * 2);

                    } else {
                        memcpy(&mixed[0], mixed_buf[cur_ptr - 66], 8 * 128 * 2);
                    }
                }

                j = (int32) (128 * (1 - latf));
                l = (int32) (128 * lonf);
                if ((mixed[j][l / 16] & pow2[15 - (l % 16)]) == 0) {

                    if (max_out_col > min_out_col) {
                        k = (out_row - min_out_row) * dims_out_sub[1] + (i - min_out_col);
                    } else {
                        k = (out_row - min_out_row) * dims_out_sub[1] +
                                (i - min_out_col) * (i >= min_out_col) +
                                (dims_out[1] - min_out_col + i - 1) * (i < max_out_col);
                    }

                    if (strcmp(input.precision, "F") == 0) {
                        memcpy(&par_byt[4 * k], &fill_val_float, 4);
                    } else if (strcmp(input.precision, "I") == 0) {
                        memcpy(&par_byt[2 * k], &fill_val, 1);
                        memcpy(&par_byt[2 * k + 1], &fill_val, 1);
                    } else if (strcmp(input.precision, "B") == 0) {
                        par_byt[k] = fill_val;
                    }
                }

            } /* pixel loop */

            /* Remove "white" coastline pixels */
            /* ------------------------------- */
            if (strcmp(input.precision, "B") == 0) {

                for (i = pixmin; i < pixmax; i++) {
                    if (max_out_col > min_out_col) {
                        if (i < min_out_col || i > max_out_col) continue;
                    } else {
                        if (i < min_out_col && i > max_out_col) continue;
                    }

                    if (max_out_col > min_out_col) {
                        k = (out_row - min_out_row) * dims_out_sub[1] + (i - min_out_col);
                    } else {
                        k = (out_row - min_out_row) * dims_out_sub[1] +
                                (i - min_out_col) * (i >= min_out_col) +
                                (dims_out[1] - min_out_col + i - 1) * (i < max_out_col);
                    }


                    /* Left Edge */
                    /* --------- */
                    if (i >= min_out_col + 2 && par_byt[k] == 0 && par_byt[k - 1] == 255) {
                        for (j = k + 1; j < dims_out[1]; j++) {
                            if (par_byt[j] != 0) {
                                i8 = par_byt[j];
                                break;
                            }
                        }

                        for (l = j - 1; l >= k; l--) {
                            par_byt[l] = i8;
                        }
                    } /* left edge */

                    /* Right Edge */
                    /* ---------- */
                    if (i <= max_out_col - 2 && par_byt[k] == 0 && par_byt[k + 1] == 255) {
                        for (j = k - 1; j >= 0; j--) {
                            if (par_byt[j] != 0) {
                                i8 = par_byt[j];
                                break;
                            }
                        }

                        for (l = j + 1; l <= k; l++) {
                            par_byt[l] = i8;
                        }
                    } /* right edge */

                } /* pixel loop (coastline) */
            } /* remove white coastline */

        } /* out_row loop */

        free(ptr_arr);
        for (i = 0; i < 360 * 180; i++) {
            if (mixed_buf[i] != NULL) free(mixed_buf[i]);
        }
        fclose(mask_fp);

    } /* End LANDMASK section */


    printf("\n");
    printf("Actual Data Min                 : %8.4f\n", dminmax[0]);
    printf("Actual Data Max                 : %8.4f\n", dminmax[1]);

    if (strcmp(input.precision, "F") == 0) {
        printf("\nData is output in float format (unscaled)\n");
    } else {
        printf("Scaled Data Min                 : %8.4f\n", aminmax[0]);
        printf("Scaled Data Max                 : %8.4f\n\n", aminmax[1]);
    }

    if (hdf4_input)
        printf("\nNumber of input bins            : %ld\n", (long)meta_l3b.data_bins);

    printf("Number of filled grid points    : %d\n", filled_data_bins);

    tgp = dims_out[0] * dims_out[1];
    printf("Total number of grid points     : %u\n", tgp);

    printf("Percentage of filled grid points: %8.4f%%\n",
            (100 * (float32) filled_data_bins) / tgp);
    printf("Output File                     : %s\n", input.ofile);



    /* Write byte data to file */
    /* ----------------------- */
    if (strncmp(input.prod, "refl", 4) == 0 ||
            strncmp(input.prod, "rhos", 4) == 0) {
        fp = fopen(input.ofile, "wb");
        fprintf(fp, "%s\n", "P5");
        fprintf(fp, "%d\n", dims_out_sub[1]);
        fprintf(fp, "%d\n", dims_out_sub[0]);
        fprintf(fp, "%d\n", 255);
        fwrite(&par_byt[0], 1, dims_out_sub[0] * dims_out_sub[1], fp);
        fclose(fp);
    } else {

        if (strcmp(input.precision, "I") == 0) {
            for (i = 0; i < L3M_PARAMS; i++) {
                if (strcmp(parmname_short[i], input.prod) == 0) {

                    prod_num = i;
                    /*
                    si8_used[0] = maximum_list[i];
                    si8_used[1] = minimum_list[i];
                     */
                    break;
                }
            }
        }

        uint8 isHDF5 = 0;
        if (getFileFormatName(input.oformat) == NULL) {
            if (hdf4_input) strcpy(input.oformat, "HDF4");
            if (hdf5_input) strcpy(input.oformat, "HDF5");
            if (ncdf4_input) strcpy(input.oformat, "netCDF4");
        }
        if (strcmp(input.oformat, "HDF5") == 0) isHDF5 = 1;
        if (strcmp(input.oformat, "netCDF4") == 0) isHDF5 = 2;
        //    if (H5Fis_hdf5( input.ifile)) isHDF5 = 1;
        //    if ( ncdf4_input == 1) isHDF5 = 2;
        set_param_string(&input);

        if (strcmp(input.precision, "F") == 0) {
            fill = (VOIDP) & fill_val_float;
        } else if (strcmp(input.precision, "I") == 0) {
            fill = (VOIDP) & fill_val_int16;
        } else if (strcmp(input.precision, "B") == 0) {
            fill = (VOIDP) & fill_val;
        }

        ptr_meta_l3b = &meta_l3b;
        if (isHDF5 == 1) ptr_meta_l3b = &input_binfile.meta_l3b;
        ptr_meta_l3b->data_bins = filled_data_bins;
        put_smi(input.ofile, input.prod, par_byt,
                dims_out_sub, &input.latnorth, &input.lonwest,
                (char*) measure_list[input.meas - 1], scale_type,
                si_used, aminmax, atype, aopt, input.ifile, dminmax,
                ptr_meta_l3b, (unsigned char *) input.palette, softid, proc_con,
                input, input.precision, qual_byt, isHDF5, fill);
    }

    free(in2out);
    free(input_databuf_16);
    free(input_databuf);
    free(input_qualbuf);
    free(out_sum);
    free(out_num);
    free(par_byt);
    free(qual_byt);
    free(best_qual);

    if (input.seam_lon != -180) {
        free(out_sum_rot);
        free(out_num_rot);
        free(best_qual_rot);
    }

    if (binlist_buf != NULL) free(binlist_buf);
    if (sum_buf != NULL) free(sum_buf);
    if (qual_buf != NULL) free(qual_buf);


    for (i = 0; i < L3M_PARAMS; i++) {
        free(parmname_list[i]);
        free(parmname_short[i]);
        free(unit_list[i]);
        free(scaling_list[i]);
        free(palette_list[i]);
        free(precision_list[i]);
    }
    free(parmname_list);
    free(parmname_short);
    free(unit_list);
    free(scaling_list);
    free(palette_list);
    free(maximum_list);
    free(minimum_list);
    free(precision_list);

    return SUCCEED;
}

int32 open_input_hdf(char *hdf_file, char *pname, int32 *fid, int32 *sdfid,
        int32 *vgid, int32 vdata_id[], meta_l3bType *meta_l3b) {
    intn i;
    int32 vg_ref;
    int32 tag;
    int32 ref;
    int32 vdid;
    char nam_buf[80];
    char cls_buf[80];
    char *tmp_str;

    if ((*fid = Hopen(hdf_file, DFACC_RDONLY, 0)) < 0) {
        fprintf(stderr, "Error: Cannot open input HDF file on Hopen - %s\n",
                hdf_file);
        return FAIL;
    }

    Vstart(*fid);

    if ((*sdfid = SDstart(hdf_file, DFACC_RDONLY)) < 0) {
        fprintf(stderr, "Error: Cannot open input HDF file on SDstart- %s\n",
                hdf_file);
        return FAIL;
    }

    read_l3b_meta_hdf4(*sdfid, meta_l3b);

    *vgid = -1;
    vdata_id[0] = -1;
    vdata_id[1] = -1;
    vdata_id[2] = -1;
    vdata_id[3] = -1;

    vg_ref = Vfind(*fid, "Level-3 Binned Data");
    // HDF4 binfile 
    if (vg_ref > 0) {
        *vgid = Vattach(*fid, vg_ref, "r");

        tmp_str = (char *) malloc(strlen(hdf_file) + 1);
        strcpy(tmp_str, hdf_file);

        HXsetdir(dirname(tmp_str));

        for (i = 0; i < Vntagrefs(*vgid); i++) {
            Vgettagref(*vgid, i, &tag, &ref);
            vdid = VSattach(*fid, ref, "r");
            VSgetname(vdid, nam_buf);
            VSgetclass(vdid, cls_buf);

            if (strcmp(cls_buf, "DataMain") == 0) {
                vdata_id[0] = vdid;
            }

            if (strcmp(cls_buf, "Index") == 0) {
                vdata_id[1] = vdid;
                nrows = VSelts(vdata_id[1]);
            }

            if (strcmp(cls_buf, "DataSubordinate") == 0 &&
                    strcmp(nam_buf, pname) == 0) {
                vdata_id[2] = vdid;
            }

            if (strcmp(cls_buf, "DataQuality") == 0) {
                vdata_id[3] = vdid;
            }
        }

        free(tmp_str);

        if (vdata_id[0] == -1) {
            printf("\"DataMain\" Vdata Not Found.\n");
            exit(-1);
        }

        if (vdata_id[1] == -1) {
            printf("\"Index\" Vdata Not Found.\n");
            exit(-1);
        }

    } else {

        /* land metadata */

        strcpy(meta_l3b->sensor_name, "SeaWiFS");
        strcpy(meta_l3b->sensor,
                "Sea-viewing Wide Field-of-view Sensor (SeaWiFS)");

        sprintf(meta_l3b->title, "%s%s", meta_l3b->sensor_name,
                " Level-3 Standard Mapped Image");

        strcpy(meta_l3b->mission, "SeaStar SeaWiFS");

        strcpy(meta_l3b->mission_char, "Nominal orbit: inclination = 98.2 (Sun-synchronous); node = 12 noon local (descending); eccentricity = <0.002; altitude = 705 km; ground speed = 6.75 km/sec");

        strcpy(meta_l3b->sensor, "SeaWiFS");

        strcpy(meta_l3b->sensor_char, "Number of bands = 8; number of active bands = 8; wavelengths per band (nm) = 412, 443, 490, 510, 555, 670, 765, 865; bits per pixel = 10; instantaneous field-of-view = 1.5835 mrad; pixels per scan = 1285; scan rate = 6/sec; sample rate = 7710/sec");

    }

    return SUCCEED;
}

int32 open_input_hdf5(char *hdf5_file, char *pname,
        Hdf::hdf5_bin *input_binfile,
        hid_t *bin_dataset_idx) {
    input_binfile->open(hdf5_file);

    if (strcmp(pname, "pixels") == 0)
        return 0;

    char ds_name[200];
    int i = 0;

    while (1) {
        hid_t dataset = input_binfile->get_data_table(i);
        H5Iget_name(dataset, ds_name, 200);
        if (strcmp(pname, &ds_name[21]) == 0) break;
        i++;

        if (i == input_binfile->nprod()) {
            cout << "Input product: \"" << pname
                    << "\" not found in input binfile." << endl;
            exit(1);
        }

    }
    *bin_dataset_idx = i;
    return 0;
}

int32 open_input_ncdf4(int ncid, const char *pname,
        size_t *nrows, int *grpid,
        int *binindex_id, int *binlist_id,
        int *bindata_id, int *binqual_id,
        meta_l3bType *meta_l3b) {
    int status;
    int dimid;

    status = nc_inq_ncid(ncid, "level-3_binned_data", grpid);
    status = nc_inq_dimid(*grpid, "binIndexDim", &dimid);
    status = nc_inq_dimlen(*grpid, dimid, nrows);

    status = nc_inq_varid(*grpid, "BinIndex", binindex_id);
    status = nc_inq_varid(*grpid, "BinList", binlist_id);

    status = nc_inq_varid(*grpid, pname, bindata_id);
    if (status != NC_NOERR) {
        cout << "Product \"" << pname << "\" not found in input file" << endl;
        exit(1);
    }

    status = nc_inq_varid(*grpid, "qual_l3", binqual_id);
    if (status != NC_NOERR) *binqual_id = -1;

    read_l3b_meta_hdf4(-ncid, meta_l3b);

    int proc_ctl_grp_id;
    status = nc_inq_grp_ncid(ncid, "processing_control", &proc_ctl_grp_id);
    status = nc_get_att(proc_ctl_grp_id, NC_GLOBAL, "l2_flag_names",
            &meta_l3b->flag_names);

    return 0;
}

int32 close_input_hdf(int32 fid, int32 sdfid, int32 vgid, int32 vdata_id[1]) {

    if (vgid != -1) {
        Vdetach(vgid);
        VSdetach(vdata_id[0]);
        if (vdata_id[2] != -1) VSdetach(vdata_id[2]);
        VSdetach(vdata_id[1]);
    }

    SDend(sdfid);
    Vend(fid);
    Hclose(fid);

    return SUCCEED;
}

int get_minmax_rowcol(float32 lat_range[], float32 lon_range[],
        char* proj_type, int32 out_rows,
        int32 *min_out_row, int32 *max_out_row,
        int32 *min_out_col, int32 *max_out_col) {
    int32 n_w, n_e, s_w, s_e, w_0, e_0;

    *min_out_row = (int32) (out_rows * (90 - lat_range[0]) / 180);
    *max_out_row = (int32) (out_rows * (90 - lat_range[1]) / 180 - 1);


    if (strcmp(proj_type, "SIN") == 0) {
        n_w = (int32) ((cos(lat_range[0] * PI / 180) * out_rows) * (lon_range[0] / 180) + out_rows);
        n_e = (int32) ((cos(lat_range[0] * PI / 180) * out_rows) * (lon_range[1] / 180) + out_rows);
        s_w = (int32) ((cos(lat_range[1] * PI / 180) * out_rows) * (lon_range[0] / 180) + out_rows);
        s_e = (int32) ((cos(lat_range[1] * PI / 180) * out_rows) * (lon_range[1] / 180) + out_rows);
        w_0 = (int32) ((out_rows) * (lon_range[0] / 180) + out_rows);
        e_0 = (int32) ((out_rows) * (lon_range[1] / 180) + out_rows);
    }

    if (strcmp(proj_type, "RECT") == 0) {
        n_w = (int32) ((out_rows) * (lon_range[0] / 180) + out_rows);
        n_e = (int32) ((out_rows) * (lon_range[1] / 180) + out_rows);
        s_w = (int32) ((out_rows) * (lon_range[0] / 180) + out_rows);
        s_e = (int32) ((out_rows) * (lon_range[1] / 180) + out_rows);
        w_0 = (int32) ((out_rows) * (lon_range[0] / 180) + out_rows);
        e_0 = (int32) ((out_rows) * (lon_range[1] / 180) + out_rows);
    }

    if ((lat_range[0] * lat_range[1]) < 0) {
        if (n_w <= s_w && n_w <= w_0) *min_out_col = n_w;
        if (s_w <= n_w && s_w <= w_0) *min_out_col = s_w;
        if (w_0 <= s_w && w_0 <= n_w) *min_out_col = w_0;

        if (n_e >= s_e && n_e >= e_0) *max_out_col = n_e - 1;
        if (s_e >= n_e && s_e >= e_0) *max_out_col = s_e - 1;
        if (e_0 >= s_e && e_0 >= n_e) *max_out_col = e_0 - 1;
    } else {
        if (n_w <= s_w) *min_out_col = n_w;
        if (s_w <= n_w) *min_out_col = s_w;

        if (n_e >= s_e) *max_out_col = n_e - 1;
        if (s_e >= n_e) *max_out_col = s_e - 1;
    }

    return 1;
}

