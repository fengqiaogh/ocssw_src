/*
 File: make_L3.c

 Compiler options:
    HP:
    cc -Ae +O2 +Ofastaccess +Onolimit make_L3.c -o make_L3 ./LIB/lib.a -lf -lm -I /usr/local/HDF4.0/include -L/usr/local/HDF4.0/lib -lmfhdf -ldf -ljpeg -lz
    SGI:
    cc -O2 -fullwarn -woff 1209,1210 -64 make_L3.c -o make_L3 ./LIB/lib.a -I./LIB -I$HDFINC -L$HDFLIB -lmfhdf -ldf -lz -ljpeg -lm -lmalloc
 or cc -O2 -fullwarn -woff 1209,1210 -64 make_L3.c -o make_L3 ./LIB/lib.a -I./LIB -I$HDFINC -L/home/jack/libhdf64bit.4.1r3/ -lmfhdf -ldf -lz -lm -ljpeg

  Revision History:
   Version 1.0     May 30, 2000, J. Descloitres
   Version 1.0.1   May 31, 2000, J. Descloitres
   Version 1.0.2  June 14, 2000, J. Descloitres
   Version 1.1    August 8, 2001, J. Gales
   Version 1.11    November 19, 2001, J. Gales

  References and Credits:
   Jacques Descloitres
   University of Maryland - Department of Geography
   NASA / GSFC  Code 922
   Bldg 32 Room S31
   Greenbelt, MD 20771
   Phone: 301 614 5456
   Email: jack@ltpmail.gsfc.nasa.gov
 */

#include <stdint.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>
#include "path.h"
#include "proto.h"
#include "proj_proto.h"
#include "mfhdf.h"
#include <genutils.h>

enum {
    PLATECARREE, SINUSOIDAL, MOLLWEIDE, ROBINSON, HAMMER, NUMBEROFPROJECTIONS
};

enum {
    MINBLUE, MAXNDVI, MINSENZ, LASTIN
};

enum {
    GLOBAL, USERBOX
};

enum {
    BAND1, BAND2, BAND3, BAND4, BAND5, BAND6, BAND7, BAND8, NBANDS
};

enum {
    L1A, L2
};

#define DBL_MAX         1.7976931348623157E+308  /* max decimal value of a "double"*/

#define DEFPROJECTION HAMMER  /* Default output projection */
#define FAST  FALSE
#define ATMCOR  TRUE  /* Do atmospheric correction */
#define UPDATE  TRUE
#define TILTSCANS FALSE
#define MAXMEMDEF 500  /* Maximum memory to allocate (in Mbytes) */
#define STEP  1
#define DEFPIXSZ 10000.  /* Pixel size (in m) for output image files (if compatible) */
#define UNDEF  -999.  /* Undefined value */
#define FILL_UINT8 255  /* Fill value */
#define FILL_INT16 -32767  /* Fill value */
#define MINREFL  0.00  /* Minimum reflectance allowed for min. blue compositing */
#define MAXSENZDEF 55.  /* Maximum view angle */
#define UMASK  002  /* Set permission for created files */
#define BLUEBAND BAND2  /* Blue band for minimum blue compositing */
#define WATER_VAL -16384  /* Water (non-land) pixel value */

#define WAVELENGTH { "412",  "443",  "490",  "510",  "555",  "670",  "765",  "865"} /* Band wavelength */
#define IRRAD  {171.07, 189.91, 194.33, 188.24, 185.93, 152.14, 123.55, 100.00} /* TOA solar irradiance */
#define TAURAY  {0.3139, 0.2341, 0.1561, 0.1324, 0.0942, 0.0439, 0.0257, 0.0155} /* TOA molecular reflectance */
#define OZONE  {.00103, .00400, .02536, .04200, .09338, .04685, .00837, .00485} /* Coef. for ozone absorption */
#define WVAPOR  {0.0000, 0.0000, 0.0000, 0.0000, 0.0008, 0.0043, 0.0007, 0.0057} /* Coef. for water vapor abs. */
#define OXYGEN  {0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000} /* Coef. for oxygen absorption */
#define UOZ  0.270 /* Ozone abundance (cm.atm) */
#define PSURF  1013. /* Sea surface pressure */
#define NLPSURF  2160 /* Number of rows in surface pressure file */
#define NPPSURF  4320 /* Number of columns in surface pressure file */

#define TWO_THIRDS 0.666666666666666666667
#define FOUR_THIRDS 1.333333333333333333333
#define TWO_PI  (OEL_PI * 2)
#define RAD2DEG  OEL_RADEG
#define NpixelsLAC 1285
#define OMF2  0.99330562 /* (1.0 - F) * (1.0 - F) */
#define R  6371007.181 /* Earth radius (m) */
#define MAXLENGTH 100
#define EPSILON  1e-6

#define N_CONVOL 7
#define CONVOL_MAX 30.0 

#define CROSS_PRODUCT(x1, y1, x2, y2, x3, y3)   ( ((x2) - (x1)) * ((y3) - (y1)) - ((y2) - (y1)) * ((x3) - (x1)) )

#define VERSION "1.1.1"

typedef struct {
    int16 *reflTOA[NBANDS], *refl[NBANDS], *solz, *senz, *phi, *smoke, *ndvi, *evi, *gemi, *ifno;
    char projection, projname[MAXLENGTH];
    char filename[MAXLENGTH];
    char *filelist, *cmd_line;
    int Nl, Np, w, h, buflines, bufstep;
    int l1, l2, p1, p2;
    double lon1, lon2, lat1, lat2;
    double xmin, xmax, ymin, ymax;
    char domain;
    double pixsz, lon_center;
    int32 sd_id;
    int startline, startday, Nfiles;
    FILE* input_fp;
    FILE* fileuse_fp;
} OUTPUT;

typedef struct {
    char ndvi, evi, gemi, angles, smoke, gzip, atmcor, TOA;
    char band[NBANDS];
    float maxsenz;
    int maxmem;
    char composite, compositename[MAXLENGTH], path[MAXLENGTH];
    char daystr[14];
    char dayend[14];
} PROCESSINFO;

typedef struct {
    float *lon, *lat, *mus, *muv, *tilt;
    int16 *ndvi, *evi, *senz, *solz, *phi;
    int Nscans, Npixels, startpix, subsampling;
    float *rho[NBANDS], *rhoTOA[NBANDS];
    int doy, doy_prev, year, orbno, step, old_iscan, corruption;
    int32 sd_id;
    double *lon1, *lat1, *lon2, *lat2;
    double *lontop, *lonbot, *lattop, *latbot;
    double *xtop, *xbot, *ytop, *ybot;
    char newfile, intersect, title[100], level;
    int32 *l2_flags;
} INPUT;

int alloc_new_array(void **data, int size, int32 num_type, int maxmem);
void check_process_parameters(PROCESSINFO *process, char **wl);
int get_geolocation_scan(int iscan, INPUT *input, PROCESSINFO process);
int get_lonlat(double xproj, double yproj, double *lon, double *lat, int projtype, double lon_center);
int get_xyproj(double lon, double lat, double *xproj, double *yproj, int projtype, double lon_center);
void init_output_file(OUTPUT *proj, INPUT input, PROCESSINFO process);
void init_output_projection(OUTPUT *proj, PROCESSINFO process);
void memcpy_block(int offset2, int offset1, int Nrows, int Ncols, OUTPUT *proj);
void parse_command_line(int argc, char **argv, OUTPUT *proj, PROCESSINFO *process);
void process_pixel(int ipix, int iline_out, int ipix_out, char **wl, INPUT input, OUTPUT *proj, PROCESSINFO process);
int project_scan(int iscan, char **wl, INPUT *input, OUTPUT *proj, PROCESSINFO process);
int read_file_block(char *filename, char **wl, int offset, int startline, int buflines, OUTPUT *proj);
int read_l2(int iscan, INPUT *input, char **wl);
void resample_geo_scan(int Npixels, float *inlon, float *inlat, double *outlon, double *outlat);
void set_buffer_parameters(OUTPUT *proj, PROCESSINFO *process);
void set_proj_parameters(OUTPUT *proj);
int update_file(char *filename, char **wl, int offset, int startline, int buflines, OUTPUT *proj, char verbose, char gzip);
void update_pixel(INPUT input, int16 ndvi, int ipix, int idx_out, int ifno, OUTPUT *proj);
void write_cmd_line(int32 sd_id, char *cmd_line);
int16 EVI(float blue, float red, float nir);
int16 NDVI(float red, float nir);
int16 GEMI(float red, float nir);
int16 SMOKE(float **rho, int ipix);
unsigned char read_mask(double lon, double lat, int16 close);

float64 sf_refl = 0.0001, sf_angle = 0.01, sf_smoke = 0.0001, sf_vi = 0.0001;

int main(int argc, char **argv) {

    char prefix[MAXLENGTH];
    char filename[MAXLENGTH], inputfile[MAXLENGTH];
    char string[100];
    int badorb;
    int ipix, iscan, j;
    char *ptr;

    OUTPUT proj;
    INPUT input;
    PROCESSINFO process;

    int16 *psurf = NULL;
    int iline_psurf, ipix_psurf, idx_psurf;
    float psurf0 = PSURF;

    char *wl[NBANDS] = WAVELENGTH;
    float taur0[NBANDS] = TAURAY;
    float aoz[NBANDS] = OZONE;
    float awv[NBANDS] = WVAPOR;
    float ao2[NBANDS] = OXYGEN;
    float Es[NBANDS] = IRRAD;
    float taur[NBANDS], rhoray[NBANDS];
    double trup[NBANDS], trdown[NBANDS];
    double coefcalib[NBANDS];
    float m, Ttotrayu, Ttotrayd, rho_pix, dsol, uoz = 0.270;
    int ib;
    float sphalb[5000];

    int32 *scan_year;
    int32 *scan_day;
    int32 *scan_msec;
    int32 start[2];
    int32 edges[2];
    int32 sds_id;
    int hr, min, sec;
    char scantime_str[14];
    int i;
    int k;
    int l;
    float32 convol;
    float32 convol_diff;
    unsigned char noisy_scan;
    float32 *center_lat;

    char* tmpstr;

    FILE *fd;

    setlinebuf(stdout);


    printf("%s %s (%s %s)\n", "LANDBIN", VERSION, __DATE__, __TIME__);

    sphalb[0] = 0;
    for (j = 1; j < 5000; j++) sphalb[j] = csalbr(j / 10000.);

    input.lon = input.lat = input.mus = input.muv = NULL;
    input.senz = input.solz = input.phi = NULL;
    input.lon1 = input.lat1 = input.lon2 = input.lat2 = NULL;
    input.lontop = input.lonbot = input.lattop = input.latbot = NULL;
    input.xtop = input.xbot = input.ytop = input.ybot = NULL;
    input.l2_flags = NULL;

    for (ib = 0; ib < NBANDS; ib++)
        input.rho[ib] = NULL;
    for (ib = 0; ib < NBANDS; ib++)
        input.rhoTOA[ib] = NULL;
    input.ndvi = NULL;
    input.evi = NULL;
    input.step = STEP;

    proj.pixsz = DEFPIXSZ;
    proj.lon_center = 0;
    proj.domain = GLOBAL;
    proj.startday = 0;
    proj.sd_id = -1;
    proj.bufstep = -1;
    proj.lon1 = proj.lon2 = proj.lat1 = proj.lat2 = UNDEF;
    proj.smoke = proj.ndvi = proj.evi = proj.gemi = proj.solz = proj.senz = proj.phi = NULL;
    proj.filelist = NULL;
    proj.Nfiles = -1;
    proj.startline = 0;
    proj.filename[0] = 0;
    proj.input_fp = NULL;
    proj.fileuse_fp = NULL;

    for (ib = 0; ib < NBANDS; ib++)
        process.band[ib] = FALSE;
    process.angles = FALSE;
    process.ndvi = FALSE;
    process.evi = FALSE;
    process.gemi = FALSE;
    process.maxmem = MAXMEMDEF;
    process.maxsenz = MAXSENZDEF;
    process.composite = MINSENZ;
    process.smoke = FALSE;
    process.gzip = FALSE;
    process.TOA = FALSE;
    process.atmcor = ATMCOR;


    strcpy(process.daystr, "0000000000000");
    strcpy(process.dayend, "9999999999999");

    if (argc == 1) exit(0);

    parse_command_line(argc, argv, &proj, &process);

    check_process_parameters(&process, wl);

    set_proj_parameters(&proj);

    printf("Initialization...(120901) \n");

    set_buffer_parameters(&proj, &process);

    init_output_projection(&proj, process);

    umask(UMASK);


    while (1) {


        if (proj.input_fp != NULL) {
            tmpstr = fgets(inputfile, MAXLENGTH, proj.input_fp);

            if (tmpstr == NULL) goto bail;

            inputfile[strlen(inputfile) - 1] = 0;

        } else {

            i = scanf("%s", inputfile);

            if (i == EOF) goto bail;
        }


        /* open input L1 HDF file */
        printf("Opening %s ...\n", inputfile);
        if ((input.sd_id = SDstart(inputfile, DFACC_RDONLY)) == -1) {
            fprintf(stderr, "Cannot open file %s\n", inputfile);
            continue;
        }
        input.newfile = TRUE;

        /* get global attributes */
        if (get_attributes(input.sd_id, prefix, &input.Nscans, &input.Npixels, &input.startpix, &input.subsampling, &input.year, &input.doy, &input.orbno, input.title) == -1) {
            SDend(input.sd_id);
            continue;
        }

        if (strcmp(input.title, "SeaWiFS Level-1A Data") == 0) {
            input.level = L1A;
            if (ATMCOR && !process.atmcor) {
                printf("This file is Level-1A data. Atmospheric correction reactivated.\n");
                process.atmcor = TRUE;
            }
        } else if (strcmp(input.title, "SeaWiFS Level-2 Data") == 0) {
            input.level = L2;
            if (process.atmcor) {
                printf("This file is Level 2 data. Atmospheric correction inhibited.\n");
                process.atmcor = FALSE;
            }
        } else {
            fprintf(stderr, "Product type not supported (%s)\n", input.title);
            SDend(input.sd_id);
            continue;
        }

        if (process.atmcor) {
            printf("Reading terrain elevation...\n");
            alloc_new_array((void *) &psurf, NLPSURF * NPPSURF, DFNT_INT16, process.maxmem);
            sprintf(process.path, "%s", ANC_DIR);
            sprintf(filename, "%s/tbase.img", process.path);
            make_psurf(psurf, filename, NLPSURF, NPPSURF, psurf0);
        }

        /* Truncate granule name */
        *(prefix + 14) = 0;

        if (proj.sd_id < 0) {
            init_output_file(&proj, input, process);
            input.doy_prev = input.doy;
        }

        printf("\n");
        printf("Year                  = %d\n", input.year);
        printf("Day of Year           = %d\n", input.doy);
        printf("Orbit Number          = %d\n", input.orbno);
        printf("Number of scans       = %d\n", input.Nscans);
        printf("Number of pixels      = %d\n", input.Npixels);

        if ((fd = fopen(BADORBLIST, "r")) == NULL)
            fprintf(stderr, "Cannot read file %s\n", BADORBLIST);
        else {
            while ((ptr = fgets(string, 100, fd)) != NULL)
                if (sscanf(string, "%d", &badorb) == 1 &&
                        input.orbno == badorb)
                    break;
            fclose(fd);
            if (ptr) {
                printf("Orbit %d is listed in the bad orbit list. Skipped.\n", input.orbno);
                SDend(input.sd_id);
                continue;
            }
        }

        if (input.doy != proj.startday)
            fprintf(stderr, "\nWARNING !!!!!: Compositing different days.\n\n");

        if (UPDATE &&
                input.doy != input.doy_prev) {

            if (update_file(proj.filename, wl, 0, proj.startline, proj.buflines, &proj, TRUE, process.gzip) == -1)
                exit(-1);
            printf("\n");

        } /* if (input.doy != input.doy_prev) */

        input.doy_prev = input.doy;

        /* Precompute constants */
        dsol = .9856 * (input.doy - 4) * OEL_PI / 180.;
        dsol = 1. / pow(1. - .01673 * cos(dsol), 2);
        for (ib = 0; ib < NBANDS; ib++)
            coefcalib[ib] = OEL_PI / Es[ib] / dsol;

        for (ib = 0; ib < NBANDS; ib++)
            if (process.band[ib])
                input.rho[ib] = (float *) realloc(input.rho[ib], input.Npixels * sizeof (float));
        if (input.level == L2) {
            for (ib = 0; ib < NBANDS; ib++)
                if (process.band[ib] && process.TOA)
                    input.rhoTOA[ib] = (float *) realloc(input.rhoTOA[ib], input.Npixels * sizeof (float));
            if (process.ndvi)
                input.ndvi = (int16 *) realloc(input.ndvi, input.Npixels * sizeof (int16));
            if (process.evi)
                input.evi = (int16 *) realloc(input.evi, input.Npixels * sizeof (int16));
        }

        input.lon = (float *) realloc(input.lon, input.Npixels * sizeof (float));
        input.lat = (float *) realloc(input.lat, input.Npixels * sizeof (float));
        input.mus = (float *) realloc(input.mus, input.Npixels * sizeof (float));
        input.muv = (float *) realloc(input.muv, input.Npixels * sizeof (float));
        input.senz = (int16 *) realloc(input.senz, input.Npixels * sizeof (int16));
        input.solz = (int16 *) realloc(input.solz, input.Npixels * sizeof (int16));
        input.phi = (int16 *) realloc(input.phi, input.Npixels * sizeof (int16));
        input.tilt = (float *) malloc(input.Nscans * sizeof (float));

        input.lon1 = (double *) realloc(input.lon1, (input.Npixels + 1) * sizeof (double));
        input.lat1 = (double *) realloc(input.lat1, (input.Npixels + 1) * sizeof (double));
        input.lon2 = (double *) realloc(input.lon2, (input.Npixels + 1) * sizeof (double));
        input.lat2 = (double *) realloc(input.lat2, (input.Npixels + 1) * sizeof (double));
        input.lontop = (double *) realloc(input.lontop, (input.Npixels + 1) * sizeof (double));
        input.lattop = (double *) realloc(input.lattop, (input.Npixels + 1) * sizeof (double));
        input.lonbot = (double *) realloc(input.lonbot, (input.Npixels + 1) * sizeof (double));
        input.latbot = (double *) realloc(input.latbot, (input.Npixels + 1) * sizeof (double));
        input.xtop = (double *) realloc(input.xtop, (input.Npixels + 1) * sizeof (double));
        input.ytop = (double *) realloc(input.ytop, (input.Npixels + 1) * sizeof (double));
        input.xbot = (double *) realloc(input.xbot, (input.Npixels + 1) * sizeof (double));
        input.ybot = (double *) realloc(input.ybot, (input.Npixels + 1) * sizeof (double));

        input.l2_flags = (int32 *) realloc(input.l2_flags, input.Npixels * sizeof (int32));

        if (input.level == L2)
            printf("Reading file %s ...\n", inputfile);
        else {
            printf("Reading SDS \"l1a_data\" in file %s ...\n", inputfile);
            if (process.atmcor)
                printf("Calibration -> Geolocation -> Computation of surface reflectance (atmospheric correction) -> Projection ...\n");
            else
                printf("Calibration -> Geolocation -> Computation of top-of-the-atmosphere reflectance (no atmospheric correction) -> Projection ...\n");
        }


        input.intersect = FALSE;
        input.old_iscan = -10;
        if (read_sds(input.sd_id, "tilt", input.tilt) == -1) continue;

        if (proj.filelist)
            proj.filelist = (char *) realloc(proj.filelist, (strlen(proj.filelist) + strlen(inputfile) + 8) * sizeof (char));
        else {
            proj.filelist = (char *) malloc((strlen(inputfile) + 6) * sizeof (char)); /* +1 for initial "\n", +3 for "0: ", */
            strcpy(proj.filelist, "\n"); /* +1 for "\n", +1 for NULL */
        }
        proj.Nfiles++;
        sprintf(proj.filelist + strlen(proj.filelist), "%d: %s\n", proj.Nfiles, inputfile);


        /* Read scan year, day, and, msec */
        /* ------------------------------ */
        scan_year = (int32 *) calloc(input.Nscans, sizeof (int32));
        scan_day = (int32 *) calloc(input.Nscans, sizeof (int32));
        scan_msec = (int32 *) calloc(input.Nscans, sizeof (int32));
        center_lat = (float32 *) calloc(input.Nscans, sizeof (float32));

        start[0] = 0;
        edges[0] = input.Nscans;

        sds_id = SDselect(input.sd_id, SDnametoindex(input.sd_id, "year"));
        if (SDreaddata(sds_id, start, NULL, edges, scan_year) < 0) {
            fprintf(stderr, "can't read SDS %s\n", "year");
            return -1;
        }
        SDendaccess(sds_id);

        sds_id = SDselect(input.sd_id, SDnametoindex(input.sd_id, "day"));
        if (SDreaddata(sds_id, start, NULL, edges, scan_day) < 0) {
            fprintf(stderr, "can't read SDS %s\n", "day");
            return -1;
        }
        SDendaccess(sds_id);

        sds_id = SDselect(input.sd_id, SDnametoindex(input.sd_id, "msec"));
        if (SDreaddata(sds_id, start, NULL, edges, scan_msec) < 0) {
            fprintf(stderr, "can't read SDS %s\n", "msec");
            return -1;
        }
        SDendaccess(sds_id);


        start[1] = input.Npixels / 2;
        edges[1] = 1;

        sds_id = SDselect(input.sd_id, SDnametoindex(input.sd_id, "latitude"));
        if (SDreaddata(sds_id, start, NULL, edges, center_lat) < 0) {
            fprintf(stderr, "can't read SDS %s\n", "latitude");
            return -1;
        }
        SDendaccess(sds_id);


        for (iscan = 0; iscan < input.Nscans; iscan += input.step) {

            if (iscan % 500 == 0) printf("Processing scan no. %d...\n", iscan);

            noisy_scan = 0;

            /* Checking that scan is after daystr and before dayend */
            /* ---------------------------------------------------- */
            hr = (int) (scan_msec[iscan] / 1000.0 / 3600.0);
            min = (int) (((scan_msec[iscan] / 1000.0 / 3600.0) - hr) * 60.0);
            sec = (int) (((((scan_msec[iscan] / 1000.0 / 3600.0) - hr) * 60.0) - min) * 60.0);

            sprintf(scantime_str, "%4d%3d%2d%2d%2d", scan_year[iscan], scan_day[iscan],
                    hr, min, sec);
            for (i = 0; i < 13; i++) if (scantime_str[i] == ' ') scantime_str[i] = '0';

            if (strcmp(scantime_str, process.daystr) < 0) {
                /*       printf("%s %s\n", process.daystr, scantime_str);*/
                continue;
            }

            if (strcmp(scantime_str, process.dayend) > 0) {
                /*       printf("%s %s\n", process.dayend, scantime_str);*/
                break;
            }


            if (!TILTSCANS &&
                    (fabs(input.tilt[iscan] - input.tilt[iscan - 1 + 2 * (iscan == 0)]) > 0.1 ||
                    fabs(input.tilt[iscan] - input.tilt[iscan + 1 - 2 * (iscan == input.Nscans - 1)]) > 0.1))
                continue;


            /* Skip "latitude" jumps */
            /* --------------------- */
            if ((fabs(center_lat[iscan] - center_lat[iscan - 1 + 2 * (iscan == 0)]) > 0.1 ||
                    fabs(center_lat[iscan] - center_lat[iscan + 1 - 2 * (iscan == input.Nscans - 1)]) > 0.1))
                continue;


            /* Read and calibrate data */
            if ((input.level == L1A && (input.corruption = compute_l1b(input.sd_id, iscan, input.rho, CALIBTABLE, input.newfile)) == -1) ||
                    (input.level == L2 && read_l2(iscan, &input, wl) == -1))
                continue;

            /* Get navigation info */
            if (input.old_iscan != iscan - input.step) {
                if (get_geolocation_scan(iscan, &input, process) == -1)
                    continue;
                resample_geo_scan(input.Npixels, input.lon, input.lat, input.lon2, input.lat2);
            }

            /* Assumption: bounding box can be directly derived from 1st and last pixels */
            /* To be checked */
            /* Must be refined anyway to not discard scans that intersect 
               the region although no pixel center falls in the region */
            if (proj.lon1 != UNDEF &&
                    (input.lat[input.Npixels - 1] > proj.lat1 ||
                    input.lat[0] < proj.lat2 ||
                    (input.lon[0] < input.lon[input.Npixels - 1] &&
                    (input.lon[0] > proj.lon2 ||
                    input.lon[input.Npixels - 1] < proj.lon1))))
                continue;


            /* L1A processing (if needed) */
            /* -------------------------- */
            if (input.level == L1A &&
                    input.corruption != 1) {

                /* Computation of reflectance & Atmospheric correction */

                for (ipix = 0; ipix < input.Npixels; ipix += input.step) {
                    if (input.solz[ipix] == FILL_INT16) continue;
                    if (process.atmcor) {
                        iline_psurf = (int) ((NLPSURF - 1) / 2. - input.lat[ipix] * (NLPSURF - 1) / OEL_PI + 0.5);
                        ipix_psurf = (int) ((NPPSURF - 1) / 2. + input.lon[ipix] * (NPPSURF - 1) / TWO_PI + 0.5);
                        idx_psurf = iline_psurf * NPPSURF + ipix_psurf;
                        for (ib = 0; ib < NBANDS; ib++)
                            if (process.band[ib]) taur[ib] = taur0[ib] * psurf[idx_psurf] / psurf0;
                        chand(input.phi[ipix] * sf_angle, input.muv[ipix], input.mus[ipix],
                                taur, rhoray, trup, trdown, NBANDS, process.band);
                        m = 1 / input.mus[ipix] + 1 / input.muv[ipix];
                    }
                    for (ib = 0; ib < NBANDS; ib++) {
                        if (process.band[ib]) {
                            if ((rho_pix = input.rho[ib][ipix]) > 0) {
                                rho_pix *= coefcalib[ib] / input.mus[ipix];
                                if (process.atmcor) {
                                    Ttotrayu = ((TWO_THIRDS + input.muv[ipix]) + (TWO_THIRDS - input.muv[ipix]) *
                                            trup[ib]) / (FOUR_THIRDS + taur[ib]);
                                    Ttotrayd = ((TWO_THIRDS + input.mus[ipix]) + (TWO_THIRDS - input.mus[ipix]) *
                                            trdown[ib]) / (FOUR_THIRDS + taur[ib]);
                                    rho_pix = (rho_pix / expf(-m * uoz * aoz[ib]) - rhoray[ib]) / Ttotrayd / Ttotrayu;
                                    if (ao2[ib] > 0) rho_pix /= expf(-m * ao2[ib]);
                                    if (awv[ib] > 0) rho_pix /= expf(-m * awv[ib]);
                                    rho_pix /= (1 + sphalb[(int) (taur[ib] * 10000 + 0.5)] * rho_pix);
                                }
                                input.rho[ib][ipix] = rho_pix;
                            } else
                                input.rho[ib][ipix] = UNDEF;
                        }
                    }
                } /* end atmospheric correction  for (ipix=0 ... */

            } /* if L1A ... */


            /* Find and flag "noisy" scans */
            /* --------------------------- */
            for (ib = 0; ib < NBANDS; ib++) {
                if (input.rho[ib] != 0x0) {
                    convol_diff = 0.0;
                    for (ipix = 0; ipix < input.Npixels; ipix += input.step) {
                        convol = 0.0;
                        for (k = -N_CONVOL / 2; k <= N_CONVOL / 2; k++) {
                            l = ipix + k;
                            if (l < 0) l = 0;
                            if (l >= input.Npixels) l = input.Npixels - 1;

                            convol += input.rho[ib][l];
                        }
                        convol_diff += fabs((convol / N_CONVOL) - input.rho[ib][ipix]);
                    } /* ipix loop */

                    if (convol_diff >= CONVOL_MAX) {
                        noisy_scan = 1;
                        break;
                    }

                } /* if (input.rho[ib] ... */
            } /* ib loop */

            if (noisy_scan) continue;


            /* Project scan and update output projection */
            if (project_scan(iscan, wl, &input, &proj, process) == -1)
                continue;

            input.old_iscan = iscan;
            input.newfile = FALSE;

        } /* for (iscan=0 ... */


        /* Close HDF input file */
        SDend(input.sd_id);


        if (!input.intersect) {
            printf("Note: file %s did not contribute at all to the output projection.\n", inputfile);
        } else {
            printf("Note: file %s actually contributed to the output projection.\n", inputfile);
            if (proj.fileuse_fp != NULL) {
                fprintf(proj.fileuse_fp, "%s\n", inputfile);
            }
        }

        printf("\n");

    } /* while (scanf("%s", inputfile) != EOF) */

bail:

    /* Close input filelist */
    /* -------------------- */
    if (proj.input_fp != NULL) fclose(proj.input_fp);

    /* Close fileuse file */
    /* ------------------ */
    if (proj.fileuse_fp != NULL) fclose(proj.fileuse_fp);


    if (update_file(proj.filename, wl, 0, proj.startline, proj.buflines, &proj, TRUE, process.gzip) == -1)
        exit(-1);


    free(input.lat);
    free(input.lon);
    free(input.senz);
    free(input.solz);
    free(input.phi);
    free(input.mus);
    free(input.muv);
    free(input.tilt);

    free(input.lat1);
    free(input.lon1);
    free(input.lat2);
    free(input.lon2);
    free(input.lattop);
    free(input.lontop);
    free(input.latbot);
    free(input.lonbot);
    free(input.xtop);
    free(input.ytop);
    free(input.xbot);
    free(input.ybot);
    free(input.l2_flags);

    for (ib = 0; ib < NBANDS; ib++) {
        if (input.rho[ib]) free(input.rho[ib]);
        if (input.rhoTOA[ib]) free(input.rhoTOA[ib]);
        if (proj.refl[ib]) free(proj.refl[ib]);
        if (proj.reflTOA[ib]) free(proj.reflTOA[ib]);
    }
    if (proj.smoke) free(proj.smoke);
    if (proj.ndvi) free(proj.ndvi);
    if (proj.evi) free(proj.evi);
    if (proj.gemi) free(proj.gemi);
    if (process.atmcor) free(psurf);
    if (proj.solz) {
        free(proj.solz);
        free(proj.senz);
        free(proj.phi);
    }
    free(proj.ifno);
    /*read_mask((double) 0, (double) 0, (int16) 1); (JMG) */

    free(scan_year);
    free(scan_day);
    free(scan_msec);
    free(center_lat);

    exit(0);

} /* end */







#define IFOV  0.0015835 /* Instantaneous FOV (rad) */

int get_geolocation_scan(int iscan, INPUT *input, PROCESSINFO process) {
    double rmtq[3], geovec[3], up[3], ea[3], no[3];
    double sv, sn, se, sunv, sunn, sune;
    double scan, a, b, c, r, q, Qx, Qz, tmp, uxy, temp, upxy;
    int j, ipix;
    float32 rm[9], pos[3], sun[3], coef[6];
    static double *cosa = NULL, *sina = NULL;
    float sena, senz, sola, solz, phi;

    if (input->level == L2)
        return 0;

    /* Precompute the cosines and sines for all scan angles */
    if (input->newfile) {
        cosa = (double *) realloc(cosa, input->Npixels * sizeof (double));
        sina = (double *) realloc(sina, input->Npixels * sizeof (double));
        for (ipix = 0; ipix < input->Npixels; ipix++) {
            scan = (double) ((input->startpix - 1 + input->subsampling * ipix - (NpixelsLAC - 1) / 2.) * IFOV);
            cosa[ipix] = cos(scan);
            sina[ipix] = sin(scan);
        }
    }

    /* get SDS related to navigation and time */
    if (get_navig_sds_line(input->sd_id, iscan, rm, pos, sun, coef, input->newfile) == -1) {
        fprintf(stderr, "Unable to get geolocation for scan %d\n", iscan);
        return -1;
    }

    for (ipix = 0; ipix < input->Npixels; ipix++) {

        /* Compute the sensor-to-surface vectors for all scan angles */
        a = coef[0] * cosa[ipix] * cosa[ipix] + coef[1] * cosa[ipix] * sina[ipix] + coef[2] * sina[ipix] * sina[ipix];
        b = coef[3] * cosa[ipix] + coef[4] * sina[ipix];
        c = coef[5];

        /* Begin the solution of the quadratic equation */
        r = b * b - 4.0 * a * c;

        /* Check for scan angles past the edge of the Earth */
        if (r < 0) {

            input->lat[ipix] = UNDEF;
            input->lon[ipix] = UNDEF;
            input->senz[ipix] = FILL_INT16;
            input->solz[ipix] = FILL_INT16;
            input->phi[ipix] = FILL_INT16;

        } else {

            /* Solve for the magnitude of the sensor-to-pixel vector and compute its components */
            q = (-b - sqrt(r)) / (2 * a);
            Qx = q * cosa[ipix];
            Qz = q * sina[ipix];

            /* Transform the vector from the sensor frame to the geocentric frame */
            for (j = 0; j < 3; j++) {
                rmtq[j] = Qx * rm[3 * j] + Qz * rm[3 * j + 2];
                geovec[j] = rmtq[j] + pos[j];
            }
            /* Compute the geodetic latitude and longitude for the scan line */
            tmp = sqrt(geovec[0] * geovec[0] + geovec[1] * geovec[1]) * OMF2;
            input->lat[ipix] = atan2(geovec[2], tmp);
            input->lon[ipix] = atan2(geovec[1], geovec[0]);

            /* Compute the local vertical, East, and North unit vectors */
            uxy = geovec[0] * geovec[0] + geovec[1] * geovec[1];
            temp = sqrt(geovec[2] * geovec[2] + OMF2 * OMF2 * uxy);
            up[0] = OMF2 * geovec[0] / temp;
            up[1] = OMF2 * geovec[1] / temp;
            up[2] = geovec[2] / temp;
            upxy = sqrt(up[0] * up[0] + up[1] * up[1]);
            ea[0] = -up[1] / upxy;
            ea[1] = up[0] / upxy;
            ea[2] = 0;
            no[0] = -up[2] * ea[1];
            no[1] = up[2] * ea[0];
            no[2] = up[0] * ea[1] - up[1] * ea[0];

            /* Decompose the spacecraft vector into vertical (up), North (no), and East (ea) components */
            sv = sn = se = 0;
            for (j = 0; j < 3; j++) {
                sv += -rmtq[j] * up[j];
                sn += -rmtq[j] * no[j];
                se += -rmtq[j] * ea[j];
            }
            /* Compute the spacecraft zenith and azimuth angles */
            senz = RAD2DEG * atan2(sqrt(sn * sn + se * se), sv);
            input->senz[ipix] = senz / sf_angle + 0.5;
            if (senz > process.maxsenz) {
                input->solz[ipix] = FILL_INT16;
                continue;
            }
            sena = RAD2DEG * atan2(se, sn);
            if (sena < 0) sena += 360;
            input->muv[ipix] = sv / sqrt(sn * sn + se * se + sv * sv);

            /* Decompose the sun vector into vertical (up), North (no), and East (ea) components */
            sunv = sunn = sune = 0;
            for (j = 0; j < 3; j++) {
                sunv += sun[j] * up[j];
                sunn += sun[j] * no[j];
                sune += sun[j] * ea[j];
            }
            /* Compute the solar zenith and azimuth angles */
            solz = RAD2DEG * atan2(sqrt(sunn * sunn + sune * sune), sunv);
            if (solz >= 90)
                input->solz[ipix] = FILL_INT16;
            else
                input->solz[ipix] = solz / sf_angle + 0.5;
            sola = RAD2DEG * atan2(sune, sunn);
            if (sola < 0) sola += 360;
            phi = sola - sena;
            if (phi > +180) phi -= 360;
            if (phi < -180) phi += 360;
            input->phi[ipix] = floor(phi / sf_angle + 0.5);
            input->mus[ipix] = sunv;

        } /* if (r < 0) ... */

    } /* for (ipix=0 ... */

    return 0;

}

void resample_geo_scan(int Npixels, float *inlon, float *inlat, double *outlon, double *outlat) {
    int ipix;

    outlon[0] = -0.5 * inlon[1] + 1.5 * inlon[0];
    outlat[0] = -0.5 * inlat[1] + 1.5 * inlat[0];
    if (inlon[1] - inlon[0] >= OEL_PI) outlon[0] += 1.5 * TWO_PI;
    if (inlon[0] - inlon[1] >= OEL_PI) outlon[0] -= 0.5 * TWO_PI;
    if (outlon[0] > OEL_PI) outlon[0] -= TWO_PI;
    if (outlon[0] < -OEL_PI) outlon[0] += TWO_PI;

    for (ipix = 0; ipix < Npixels - 1; ipix++) {
        outlon[ipix + 1] = 0.5 * inlon[ipix] + 0.5 * inlon[ipix + 1];
        outlat[ipix + 1] = 0.5 * inlat[ipix] + 0.5 * inlat[ipix + 1];
        if (inlon[ipix] - inlon[ipix + 1] >= OEL_PI) outlon[ipix + 1] += 0.5 * TWO_PI;
        if (inlon[ipix + 1] - inlon[ipix] >= OEL_PI) outlon[ipix + 1] -= 0.5 * TWO_PI;
        if (outlon[ipix + 1] > OEL_PI) outlon[ipix + 1] -= TWO_PI;
        if (outlon[ipix + 1] < -OEL_PI) outlon[ipix + 1] += TWO_PI;
    }

    outlon[Npixels] = 1.5 * inlon[Npixels - 1] - 0.5 * inlon[Npixels - 2];
    outlat[Npixels] = 1.5 * inlat[Npixels - 1] - 0.5 * inlat[Npixels - 2];
    if (inlon[Npixels - 1] - inlon[Npixels - 2] >= OEL_PI) outlon[Npixels] += -0.5 * TWO_PI;
    if (inlon[Npixels - 2] - inlon[Npixels - 1] >= OEL_PI) outlon[Npixels] -= -0.5 * TWO_PI;
    if (outlon[Npixels] > OEL_PI) outlon[Npixels] -= TWO_PI;
    if (outlon[Npixels] < -OEL_PI) outlon[Npixels] += TWO_PI;

}

int read_file_block(char *filename, char **wl, int offset, int startline, int buflines, OUTPUT *proj) {
    int32 sd_id, sds_id;
    int ib;
    char SDSname[H4_MAX_NC_NAME];

    printf("Reading  %s (lines %d to %d)...\n", filename, startline, startline + buflines - 1);
    if ((sd_id = SDstart(filename, DFACC_RDWR)) == -1) {
        fprintf(stderr, "Cannot open file %s\n", filename);
        return -1;
    }
    for (ib = 0; ib < NBANDS; ib++) {
        if (proj->refl[ib]) {
            sprintf(SDSname, "refl_%s", wl[ib]);
            if ((sds_id = open_sds(sd_id, SDSname)) == -1 ||
                    read_sds_block(sds_id, startline, buflines, proj->refl[ib] + offset) == -1)
                return -1;
        }
        if (proj->reflTOA[ib]) {
            sprintf(SDSname, "reflTOA_%s", wl[ib]);
            if ((sds_id = open_sds(sd_id, SDSname)) == -1 ||
                    read_sds_block(sds_id, startline, buflines, proj->reflTOA[ib] + offset) == -1)
                return -1;
        }
    }
    sprintf(SDSname, "input_file");
    if ((sds_id = open_sds(sd_id, SDSname)) == -1 ||
            read_sds_block(sds_id, startline, buflines, proj->ifno + offset) == -1)
        return -1;
    if (proj->senz) {
        if (sprintf(SDSname, "senz") == -1 ||
                (sds_id = open_sds(sd_id, SDSname)) == -1 ||
                read_sds_block(sds_id, startline, buflines, proj->senz + offset) == -1)
            return -1;
    }
    if (proj->solz) {
        if (sprintf(SDSname, "solz") == -1 ||
                (sds_id = open_sds(sd_id, SDSname)) == -1 ||
                read_sds_block(sds_id, startline, buflines, proj->solz + offset) == -1 ||
                sprintf(SDSname, "phi") == -1 ||
                (sds_id = open_sds(sd_id, SDSname)) == -1 ||
                read_sds_block(sds_id, startline, buflines, proj->phi + offset) == -1)
            return -1;
    }
    if (proj->smoke) {
        sprintf(SDSname, "smoke_index");
        if ((sds_id = open_sds(sd_id, SDSname)) == -1 ||
                read_sds_block(sds_id, startline, buflines, proj->smoke + offset) == -1)
            return -1;
    }
    if (proj->ndvi) {
        sprintf(SDSname, "NDVI");
        if ((sds_id = open_sds(sd_id, SDSname)) == -1 ||
                read_sds_block(sds_id, startline, buflines, proj->ndvi + offset) == -1)
            return -1;
    }
    if (proj->evi) {
        sprintf(SDSname, "EVI");
        if ((sds_id = open_sds(sd_id, SDSname)) == -1 ||
                read_sds_block(sds_id, startline, buflines, proj->evi + offset) == -1)
            return -1;
    }
    if (proj->gemi) {
        sprintf(SDSname, "GEMI");
        if ((sds_id = open_sds(sd_id, SDSname)) == -1 ||
                read_sds_block(sds_id, startline, buflines, proj->gemi + offset) == -1)
            return -1;
    }
    SDend(sd_id);

    return 0;

}

int update_file(char *filename, char **wl, int offset, int startline, int buflines, OUTPUT *proj, char verbose, char gzip) {
    int32 sd_id;
    int ib;
    char SDSname[H4_MAX_NC_NAME];
    int16 fill_int16 = FILL_INT16;

    if (startline + buflines > proj->h) buflines = proj->h - startline;
    if (startline >= 0)
        printf("Updating %s (lines %d to %d)...\n", filename, startline, startline + buflines - 1);
    if ((sd_id = SDstart(filename, DFACC_RDWR)) == -1) {
        fprintf(stderr, "Cannot open file %s\n", filename);
        return -1;
    }
    for (ib = 0; ib < NBANDS; ib++) {
        if (proj->refl[ib]) {
            sprintf(SDSname, "refl_%s", wl[ib]);
            write_sds(sd_id, SDSname, proj->refl[ib] + offset, startline, buflines, proj->h, proj->w, DFNT_INT16, &fill_int16, sf_refl, verbose, gzip);
        }
    }
    for (ib = 0; ib < NBANDS; ib++) {
        if (proj->reflTOA[ib]) {
            sprintf(SDSname, "reflTOA_%s", wl[ib]);
            write_sds(sd_id, SDSname, proj->reflTOA[ib] + offset, startline, buflines, proj->h, proj->w, DFNT_INT16, &fill_int16, sf_refl, verbose, gzip);
        }
    }
    sprintf(SDSname, "input_file");
    write_sds(sd_id, SDSname, proj->ifno + offset, startline, buflines, proj->h, proj->w, DFNT_INT16, &fill_int16, 0, verbose, gzip);
    if (proj->senz) {
        sprintf(SDSname, "senz");
        write_sds(sd_id, SDSname, proj->senz + offset, startline, buflines, proj->h, proj->w, DFNT_INT16, &fill_int16, sf_angle, verbose, gzip);
    }
    if (proj->solz) {
        sprintf(SDSname, "solz");
        write_sds(sd_id, SDSname, proj->solz + offset, startline, buflines, proj->h, proj->w, DFNT_INT16, &fill_int16, sf_angle, verbose, gzip);
        sprintf(SDSname, "phi");
        write_sds(sd_id, SDSname, proj->phi + offset, startline, buflines, proj->h, proj->w, DFNT_INT16, &fill_int16, sf_angle, verbose, gzip);
    }
    if (proj->smoke) {
        sprintf(SDSname, "smoke_index");
        write_sds(sd_id, SDSname, proj->smoke + offset, startline, buflines, proj->h, proj->w, DFNT_INT16, &fill_int16, sf_smoke, verbose, gzip);
    }
    if (proj->ndvi) {
        sprintf(SDSname, "NDVI");
        write_sds(sd_id, SDSname, proj->ndvi + offset, startline, buflines, proj->h, proj->w, DFNT_INT16, &fill_int16, sf_vi, verbose, gzip);
    }
    if (proj->evi) {
        sprintf(SDSname, "EVI");
        write_sds(sd_id, SDSname, proj->evi + offset, startline, buflines, proj->h, proj->w, DFNT_INT16, &fill_int16, sf_vi, verbose, gzip);
    }
    if (proj->gemi) {
        sprintf(SDSname, "GEMI");
        write_sds(sd_id, SDSname, proj->gemi + offset, startline, buflines, proj->h, proj->w, DFNT_INT16, &fill_int16, sf_vi, verbose, gzip);
    }

    if (SDsetattr(sd_id, "File list", DFNT_CHAR8, (int) strlen(proj->filelist) + 1, proj->filelist) == -1)
        fprintf(stderr, "Cannot write attribute in file %s\n", filename);

    SDend(sd_id);
    return 0;
}

void memcpy_block(int offset2, int offset1, int Nrows, int Ncols, OUTPUT *proj) {
    int ib, k, irow;
    if (offset1 < 0 || offset2 < 0) return;
    for (k = 0; k < Nrows; k++) {
        if (offset2 < offset1) irow = k;
        else irow = Nrows - 1 - k;
        for (ib = 0; ib < NBANDS; ib++) {
            if (proj->refl[ib])
                memcpy(proj->refl[ib] + (offset2 + irow) * Ncols, proj->refl[ib] + (offset1 + irow) * Ncols, Ncols * sizeof (int16));
            if (proj->reflTOA[ib])
                memcpy(proj->reflTOA[ib] + (offset2 + irow) * Ncols, proj->reflTOA[ib] + (offset1 + irow) * Ncols, Ncols * sizeof (int16));
        }
        memcpy(proj->ifno + (offset2 + irow) * Ncols, proj->ifno + (offset1 + irow) * Ncols, Ncols * sizeof (int16));
        if (proj->senz)
            memcpy(proj->senz + (offset2 + irow) * Ncols, proj->senz + (offset1 + irow) * Ncols, Ncols * sizeof (int16));
        if (proj->solz) {
            memcpy(proj->solz + (offset2 + irow) * Ncols, proj->solz + (offset1 + irow) * Ncols, Ncols * sizeof (int16));
            memcpy(proj->phi + (offset2 + irow) * Ncols, proj->phi + (offset1 + irow) * Ncols, Ncols * sizeof (int16));
        }
        if (proj->smoke)
            memcpy(proj->smoke + (offset2 + irow) * Ncols, proj->smoke + (offset1 + irow) * Ncols, Ncols * sizeof (int16));
        if (proj->ndvi)
            memcpy(proj->ndvi + (offset2 + irow) * Ncols, proj->ndvi + (offset1 + irow) * Ncols, Ncols * sizeof (int16));
        if (proj->evi)
            memcpy(proj->evi + (offset2 + irow) * Ncols, proj->evi + (offset1 + irow) * Ncols, Ncols * sizeof (int16));
        if (proj->gemi)
            memcpy(proj->gemi + (offset2 + irow) * Ncols, proj->gemi + (offset1 + irow) * Ncols, Ncols * sizeof (int16));
    }
}

void process_pixel(int ipix, int iline_out, int ipix_out, char **wl, INPUT input, OUTPUT *proj, PROCESSINFO process) {
    int idx_out, old_startline, offset;
    int16 ndvi;
    char better;

    iline_out -= proj->l1;
    ipix_out -= proj->p1;

    if (iline_out < proj->startline ||
            iline_out > proj->startline + proj->buflines - 1) {
        old_startline = proj->startline;
        if (proj->startline < 0)
            proj->startline = iline_out - iline_out % proj->bufstep;
        while (iline_out < proj->startline)
            proj->startline -= proj->bufstep;
        while (iline_out > proj->startline + proj->buflines - 1)
            proj->startline += proj->bufstep;
        if (proj->startline < 0)
            proj->startline = 0;
        if (proj->startline + proj->buflines > proj->h)
            proj->startline = proj->h - proj->buflines;

        if (proj->startline > old_startline &&
                proj->startline < old_startline + proj->buflines - 1) {
            offset = proj->startline - old_startline;
            if (update_file(proj->filename, wl, 0, old_startline, offset, proj, FALSE, process.gzip) == -1)
                exit(-1);
            memcpy_block(0, offset, proj->buflines - offset, proj->w, proj);
            if (read_file_block(proj->filename, wl, (proj->buflines - offset) * proj->w, proj->startline + proj->buflines - offset,
                    MIN(offset, proj->h - proj->startline - proj->buflines + offset), proj) == -1)
                exit(-1);
        } else if (old_startline > proj->startline &&
                old_startline < proj->startline + proj->buflines - 1) {
            offset = old_startline - proj->startline;
            if (update_file(proj->filename, wl, (proj->buflines - offset) * proj->w, old_startline + proj->buflines - offset,
                    offset, proj, FALSE, process.gzip) == -1)
                exit(-1);
            memcpy_block(offset, 0, proj->buflines - offset, proj->w, proj);
            if (read_file_block(proj->filename, wl, 0, proj->startline,
                    MIN(offset, proj->h - proj->startline), proj) == -1)
                exit(-1);
        } else {
            if (update_file(proj->filename, wl, 0, old_startline, proj->buflines, proj, FALSE, process.gzip) == -1)
                exit(-1);
            if (read_file_block(proj->filename, wl, 0, proj->startline,
                    MIN(proj->buflines, proj->h - proj->startline), proj) == -1)
                exit(-1);
        }
    }

    idx_out = (iline_out - proj->startline) * proj->w + ipix_out;
    better = FALSE;

    if (proj->ifno[idx_out] == FILL_INT16) /* First observation for this pixel */

        better = TRUE;

    else {

        switch (process.composite) {

        case MINBLUE:
            if (input.rho[BLUEBAND][ipix] >= MINREFL &&
                    (input.rho[BLUEBAND][ipix] < proj->refl[BLUEBAND][idx_out] * sf_refl ||
                    proj->refl[BLUEBAND][idx_out] == FILL_INT16))
                better = TRUE;
            break;

        case MAXNDVI:
            if (input.level == L1A)
                ndvi = NDVI(input.rho[BAND6][ipix], input.rho[BAND8][ipix]);
            else
                ndvi = input.ndvi[ipix];
            if (ndvi != FILL_INT16 &&
                    (ndvi > proj->ndvi[idx_out] ||
                    proj->ndvi[idx_out] == FILL_INT16))
                better = TRUE;
            break;

        case MINSENZ:
            if (input.senz[ipix] < proj->senz[idx_out])
                better = TRUE;
            break;

        case LASTIN:
            better = TRUE;
            break;

        }
    }

    if (better) {
        if (proj->ndvi &&
                (proj->ifno[idx_out] == FILL_INT16 ||
                process.composite != MAXNDVI)) {
            if (input.level == L1A)
                ndvi = NDVI(input.rho[BAND6][ipix], input.rho[BAND8][ipix]);
            else
                ndvi = input.ndvi[ipix];
        }
        update_pixel(input, ndvi, ipix, idx_out, proj->Nfiles, proj);
    }

}

int get_xyproj(double lon, double lat, double *xproj, double *yproj, int projtype, double lon_center) {

    lon -= lon_center;
    if (lon < -OEL_PI) lon += TWO_PI;
    if (lon > OEL_PI) lon -= TWO_PI;

    switch (projtype) {
    case PLATECARREE: *xproj = lon * RAD2DEG; /* PLATECARREE was not initialized with RAD2DEG factor */
        *yproj = lat * RAD2DEG;
        break;
    case SINUSOIDAL: *xproj = lon * cos(lat) * RAD2DEG; /* SINUSOIDAL was not initialized with RAD2DEG factor */
        *yproj = lat * RAD2DEG;
        break;
    case MOLLWEIDE: if (molwfor(lon, lat, xproj, yproj) == -1) return -1;
        break;
    case ROBINSON: if (robfor(lon, lat, xproj, yproj) == -1) return -1;
        break;
    case HAMMER: if (hamfor(lon, lat, xproj, yproj) == -1) return -1;
        break;
    default: return -1;
    }
    return 0;
}

int get_lonlat(double xproj, double yproj, double *lon, double *lat, int projtype, double lon_center) {

    switch (projtype) {
    case PLATECARREE: *lon = xproj / RAD2DEG;
        *lat = yproj / RAD2DEG;
        break;
    case SINUSOIDAL: *lat = yproj / RAD2DEG;
        if (fabs(*lat) >= OEL_PI) return -1;
        *lon = xproj / cos(*lat) / RAD2DEG;
        if (fabs(*lon) > OEL_PI) return -1;
        break;
    case MOLLWEIDE: if (molwinv(xproj, yproj, lon, lat) == -1) return -1;
        break;
    case ROBINSON: if (robinv(xproj, yproj, lon, lat) == -1) return -1;
        break;
    case HAMMER: if (haminv(xproj, yproj, lon, lat) == -1) return -1;
        break;
    default: return -1;
    }

    *lon += lon_center;
    if (*lon < -OEL_PI) *lon += TWO_PI;
    if (*lon > OEL_PI) *lon -= TWO_PI;

    return 0;

}

void update_pixel(INPUT input, int16 ndvi, int ipix, int idx_out, int ifno, OUTPUT *proj) {
    int ib;

    for (ib = 0; ib < NBANDS; ib++) {
        if (input.rho[ib]) {
            if (input.rho[ib][ipix] != UNDEF &&
                    input.rho[ib][ipix] / sf_refl <= SHRT_MAX)
                proj->refl[ib][idx_out] = (int16) floor(input.rho[ib][ipix] / sf_refl + 0.5);
            else
                proj->refl[ib][idx_out] = FILL_INT16;
        }

        if (input.rhoTOA[ib]) {
            proj->reflTOA[ib][idx_out] = (int16) floor(input.rhoTOA[ib][ipix] / sf_refl + 0.5);
        }
    }


    if (proj->senz)
        proj->senz[idx_out] = input.senz[ipix];
    if (proj->solz)
        proj->solz[idx_out] = input.solz[ipix];
    if (proj->phi)
        proj->phi[idx_out] = input.phi[ipix];

    /*
    mask_lonlat = read_mask((double) (input.lon[ipix]*RAD2DEG), 
                            (double) (input.lat[ipix]*RAD2DEG), (int16) 0);


    mask_lonlat = (input.l2_flags[ipix] & 0x2) / 0x2;

    if (mask_lonlat == 0 && proj->ifno[idx_out] == FILL_INT16) proj->ifno[idx_out] = WATER_VAL;
    if (mask_lonlat == 1)       proj->ifno[idx_out] = ifno;
     */

    proj->ifno[idx_out] = ifno;

    if (proj->smoke)
        proj->smoke[idx_out] = SMOKE(input.rho, ipix);

    if (proj->ndvi)
        proj->ndvi[idx_out] = ndvi;

    if (proj->evi) {
        if (input.level == L1A)
            proj->evi[idx_out] = EVI(input.rho[BAND2][ipix], input.rho[BAND6][ipix], input.rho[BAND8][ipix]);
        else
            proj->evi[idx_out] = input.evi[ipix];
    }

    if (proj->gemi)
        proj->gemi[idx_out] = GEMI(input.rho[BAND6][ipix], input.rho[BAND8][ipix]);

}

int16 SMOKE(float **rho, int ipix) {
    if (rho[BAND1][ipix] < 0 ||
            rho[BAND6][ipix] > 0.4)
        return FILL_INT16;
    else if (rho[BAND2][ipix] < rho[BAND4][ipix] ||
            rho[BAND5][ipix] < rho[BAND6][ipix])
        return 0;
    else
        return (int16) (rho[BAND1][ipix] * (rho[BAND2][ipix] - rho[BAND4][ipix]) * (rho[BAND5][ipix] - rho[BAND6][ipix]) * 10000);
}

int16 NDVI(float red, float nir) {
    float ndvi;

    if (red != UNDEF &&
            nir != UNDEF &&
            red + nir != 0) {
        ndvi = (nir - red) / (nir + red);
        if (ndvi < -2.) ndvi = -2.;
        if (ndvi > 2.) ndvi = 2.;
        return (int16) floor(ndvi / sf_vi + 0.5);
    } else
        return FILL_INT16;
}

int16 EVI(float blue, float red, float nir) {
    float evi;
    static float L = 1, c1 = 6, c2 = 7.5;
    double val;

    if (red == UNDEF ||
            nir == UNDEF)
        return FILL_INT16;
    else {
        if (blue != UNDEF && /* Most cases - EVI formula */
                (blue <= red || red <= nir)) {
            if ((val = L + nir + c1 * red - c2 * blue) == 0)
                return FILL_INT16;
            else {
                evi = 2.5 * (nir - red) / val;
                if (evi < -2.) evi = -2.;
                if (evi > 2.) evi = 2.;
            }
        } else { /* Backup - SAVI formula */
            if ((val = 0.5 + nir + red) == 0)
                return FILL_INT16;
            else {
                evi = 1.5 * (nir - red) / val;
                if (evi < -2.) evi = -2.;
                if (evi > 2.) evi = 2.;
            }
        }
        return (int16) floor(evi / sf_vi + 0.5);
    }

}

#if 0

int16 EVI(float blue, float red, float nir) {
    float evi;
    static float L = 1, c1 = 6, c2 = 7.5;
    double val;

    if (blue != UNDEF &&
            red != UNDEF &&
            nir != UNDEF &&
            (val = L + nir + c1 * red - c2 * blue) != 0) {
        evi = 2.5 * (nir - red) / val;
        if (evi < -2.) evi = -2.;
        if (evi > 2.) evi = 2.;
        return (int16) floor(evi / sf_vi + 0.5);
    } else
        return FILL_INT16;
}
#endif

int16 GEMI(float red, float nir) {
    float gemi;
    double x;
    if (red != UNDEF &&
            nir != UNDEF &&
            red + nir + 0.5 != 0 &&
            nir != 1.) {
        x = (2 * (nir * nir - red * red) + 1.5 * nir + 0.5 * red) / (red + nir + 0.5);
        gemi = x * (1 - 0.25 * x) - (red - 0.125) / (1 - red);
        if (gemi < -2.) gemi = -2.;
        if (gemi > 2.) gemi = 2.;
        return (int16) floor(gemi / sf_vi + 0.5);
    } else
        return FILL_INT16;
}

int alloc_new_array(void **data, int size, int32 num_type, int maxmem) {
    int sizeoftype;

    if ((*data = (void *) malloc(size * (sizeoftype = DFKNTsize(num_type)))) == NULL) { /* Allocate a new array */
        printf("  Sorry, can't allocate %.1f Mbytes of memory.\n", size * sizeoftype / 1e6);
        exit(-1);
    }
    printf(" (%.1f Mbytes allocated: %d items of %d bytes)\n", size * sizeoftype / 1e6, size, sizeoftype);
    return 0;
}

void write_cmd_line(int32 sd_id, char *cmd_line) {
    char *old_cmd_line, *new_cmd_line;
    int32 attr_index, num_type, count;
    char name[H4_MAX_NC_NAME] = "Command line";

    new_cmd_line = strdup(cmd_line); /* Default value is current command line */

    if ((attr_index = SDfindattr(sd_id, name)) != -1 &&
            SDattrinfo(sd_id, attr_index, name, &num_type, &count) != -1) {

        old_cmd_line = (char *) malloc((count + 1) * sizeof (char));
        if (SDreadattr(sd_id, attr_index, old_cmd_line) != -1) {
            new_cmd_line = (char *) realloc(new_cmd_line, (count + strlen(cmd_line) + 1 + 1) * sizeof (char));
            sprintf(new_cmd_line, "%s%s", old_cmd_line, cmd_line);
        }
        free(old_cmd_line);

    }

    if (SDsetattr(sd_id, name, DFNT_CHAR8, (int) strlen(new_cmd_line) + 1, new_cmd_line) == -1)
        printf("Unable to write command line in file attributes.\n");

    free(new_cmd_line);

}

void set_proj_parameters(OUTPUT *proj) {

    printf("Pixel size for output projection: %g m\n", proj->pixsz);
    proj->pixsz /= R / RAD2DEG;

    switch (proj->projection) {
    case PLATECARREE: sprintf(proj->projname, "Plate Carree");
        proj->Nl = ceil(180. / proj->pixsz);
        proj->Np = ceil(360. / proj->pixsz);
        break;
    case SINUSOIDAL: sprintf(proj->projname, "Sinusoidal");
        proj->Nl = ceil(180. / proj->pixsz);
        proj->Np = ceil(360. / proj->pixsz);
        break;
    case MOLLWEIDE: sprintf(proj->projname, "Mollweide");
        proj->Nl = ceil(2 * sqrt(2.) / OEL_PI * 180. / proj->pixsz);
        proj->Np = ceil(2 * sqrt(2.) / OEL_PI * 360. / proj->pixsz);
        molwforint(RAD2DEG, 0., 0., 0.);
        molwinvint(RAD2DEG, 0., 0., 0.);
        break;
    case ROBINSON: sprintf(proj->projname, "Robinson");
        proj->Nl = ceil(180. / proj->pixsz);
        proj->Np = ceil(0.9858 * 360. / proj->pixsz);
        robforint(RAD2DEG, 0., 0., 0.);
        robinvint(RAD2DEG, 0., 0., 0.);
        break;
    case HAMMER: sprintf(proj->projname, "Hammer-Aitoff");
        proj->Nl = ceil(2 * RAD2DEG * sqrt(2.) / proj->pixsz);
        proj->Np = ceil(4 * RAD2DEG * sqrt(2.) / proj->pixsz);
        hamforint(RAD2DEG, 0., 0., 0.);
        haminvint(RAD2DEG, 0., 0., 0.);
        break;
    default: fprintf(stderr, "Wrong selection for output projection\n");
        exit(-1);
    }

    printf("Output projection: %s\n", proj->projname);

    if (proj->domain == USERBOX &&
            (proj->l1 < 0 ||
            proj->p1 < 0 ||
            proj->l2 > proj->Nl - 1 ||
            proj->p2 > proj->Np - 1)) {
        if (proj->l1 < 0) proj->l1 = 0;
        if (proj->p1 < 0) proj->p1 = 0;
        if (proj->l2 > proj->Nl - 1) proj->l2 = proj->Nl - 1;
        if (proj->p2 > proj->Np - 1) proj->p2 = proj->Np - 1;
        printf("Requested subwindow cannot exceed the projection dimensions.\n");
        printf("Output subwindow has been adjusted: (column %d, row %d) to (column %d, row %d)\n",
                proj->p1, proj->l1, proj->p2, proj->l2);
    }

    if (proj->domain == GLOBAL) { /* No bounding box has been requested */
        printf("No bounding box requested. Output will be a global projection.\n");
        proj->l1 = 0;
        proj->p1 = 0;
        proj->l2 = proj->Nl - 1;
        proj->p2 = proj->Np - 1;
    }

    proj->h = proj->l2 - proj->l1 + 1;
    proj->w = proj->p2 - proj->p1 + 1;
    printf("Output projection size: %dx%d\n", proj->w, proj->h);

    proj->xmin = proj->pixsz * (proj->p1 - (proj->Np - 1) / 2. - 0.5);
    proj->xmax = proj->pixsz * (proj->p2 - (proj->Np - 1) / 2. + 0.5);
    proj->ymin = proj->pixsz * ((proj->Nl - 1) / 2. - proj->l2 - 0.5);
    proj->ymax = proj->pixsz * ((proj->Nl - 1) / 2. - proj->l1 + 0.5);

}

void parse_command_line(int argc, char **argv, OUTPUT *proj, PROCESSINFO *process) {
    char *arg, *ptr;
    int j, ib, projection;

    proj->cmd_line = (char *) malloc((strlen(argv[0]) + 1 + 1) * sizeof (char));
    sprintf(proj->cmd_line, "\n%s", argv[0]);

    for (j = 1; j < argc; j++) {

        arg = argv[j];

        proj->cmd_line = (char *) realloc(proj->cmd_line, (strlen(proj->cmd_line) + strlen(arg) + 1 + 1) * sizeof (char));
        strcat(proj->cmd_line, " ");
        strcat(proj->cmd_line, arg);

        if (strstr(arg, "-geobox=") == arg) {
            if (sscanf(arg + strlen("-geobox="), "%lg,%lg,%lg,%lg",
                    &proj->lon1, &proj->lat1, &proj->lon2, &proj->lat2) == 4) {
                printf("Geographic area selected: (lon,lat)=(%g,%g)->(%g,%g)\n",
                        proj->lon1, proj->lat1, proj->lon2, proj->lat2);
                if (proj->lon1 >= proj->lon2 ||
                        proj->lat1 <= proj->lat2 ||
                        proj->lon1 < -180 ||
                        proj->lon2 > 180 ||
                        proj->lat1 > 90 ||
                        proj->lat2 < -90) {
                    printf("Bad specification for geographic area. Request ignored.\n");
                    proj->lon1 = proj->lon2 = proj->lat1 = proj->lat2 = UNDEF;
                } else {
                    proj->lon1 /= RAD2DEG;
                    proj->lon2 /= RAD2DEG;
                    proj->lat1 /= RAD2DEG;
                    proj->lat2 /= RAD2DEG;
                }
                printf("\n");
            }
            continue;
        }

        if (strstr(arg, "-box=") == arg) {
            if (sscanf(arg + strlen("-box="), "%d,%d,%d,%d",
                    &proj->p1, &proj->l1, &proj->p2, &proj->l2) == 4) {
                printf("Output subwindow selected: (column %d, row %d) to (column %d, row %d)\n",
                        proj->p1, proj->l1, proj->p2, proj->l2);
                if (proj->l1 >= proj->l2 ||
                        proj->p1 >= proj->p2)
                    printf("Bad specification for output subwindow. Request ignored.\n");
                else
                    proj->domain = USERBOX;
                printf("\n");
            }
            continue;
        }

        if (strstr(arg, "-bands=") == arg) {
            arg += strlen("-bands=");
            while ((ptr = strchr(arg, 44))) *ptr = 32;
            ptr = arg;
            while (sscanf(ptr, "%d", &ib) == 1) {
                if (ib >= 1 && ib <= NBANDS) process->band[ib - 1] = TRUE;
                if ((ptr = strchr(ptr, 32))) ptr++;
                else ptr = strchr(arg, 0);
            }
            continue;
        }

        if (strstr(arg, "-center=") == arg) {
            arg += strlen("-center=");
            if (sscanf(arg, "%lg", &proj->lon_center) == 1)
                printf("Center longitude for output projection: %g deg.\n", proj->lon_center);
            proj->lon_center /= RAD2DEG;
            continue;
        }

        if (strstr(arg, "-pixsz=") == arg) {
            arg += strlen("-pixsz=");
            sscanf(arg, "%lg", &proj->pixsz);
            continue;
        }

        if (strstr(arg, "-maxsenz=") == arg) {
            arg += strlen("-maxsenz=");
            sscanf(arg, "%g", &process->maxsenz);
            printf("Maximum sensor zenith angle allowed: %g deg.\n", process->maxsenz);
            continue;
        }

        if (strstr(arg, "-maxmem=") == arg) {
            arg += strlen("-maxmem=");
            sscanf(arg, "%d", &process->maxmem);
            printf("Maximum memory allocation allowed: %d Mbytes\n", process->maxmem);
            continue;
        }

        if (strstr(arg, "-daystr=") == arg) {
            arg += strlen("-daystr=");
            strncpy(&process->daystr[0], arg, 13);
            continue;
        }

        if (strstr(arg, "-dayend=") == arg) {
            arg += strlen("-dayend=");
            strncpy(&process->dayend[0], arg, 13);
            continue;
        }

        if (strcmp(arg, "-smoke") == 0) {
            printf("Smoke index requested\n");
            process->smoke = TRUE;
            continue;
        }

        if (strcmp(arg, "-ndvi") == 0) {
            printf("NDVI requested\n");
            process->ndvi = TRUE;
            continue;
        }

        if (strcmp(arg, "-evi") == 0) {
            printf("EVI requested\n");
            process->evi = TRUE;
            continue;
        }

        if (strcmp(arg, "-gemi") == 0) {
            printf("GEMI requested\n");
            process->gemi = TRUE;
            continue;
        }

        if (strcmp(arg, "-angles") == 0) {
            printf("View and sun angles requested\n");
            process->angles = TRUE;
            continue;
        }

        if (strcmp(arg, "-gzip") == 0) {
            printf("Compression of output SDSs requested\n");
            process->gzip = TRUE;
            continue;
        }

        if (strcmp(arg, "-TOA") == 0) {
            printf("TOA reflectance requested in output\n");
            process->TOA = TRUE;
            continue;
        }

        if (strstr(arg, "-method=") == arg) {
            arg += strlen("-method=");
            if (strcmp(arg, "minblue") == 0)
                process->composite = MINBLUE;
            else if (strcmp(arg, "maxndvi") == 0)
                process->composite = MAXNDVI;
            else if (strcmp(arg, "minsenz") == 0)
                process->composite = MINSENZ;
            else if (strcmp(arg, "lastin") == 0)
                process->composite = LASTIN;
            else
                printf("Compositing method \"%s\" not supported\n", arg);
            continue;
        }

        if (strstr(arg, "-proj=") == arg) {
            arg += strlen("-proj=");
            projection = DEFPROJECTION;
            if (sscanf(arg, "%d", &projection) == 1) {
                if (projection < 0 ||
                        projection > NUMBEROFPROJECTIONS - 1) { /* char is >= 0 anyway */
                    printf("Invalid projection requested (%d). Request ignored.\n", projection);
                    projection = DEFPROJECTION;
                }
                proj->projection = projection;
            }
            continue;
        }

        if (strstr(arg, "-bufstep=") == arg) {
            arg += strlen("-bufstep=");
            if (sscanf(arg, "%d", &proj->bufstep) == 1)
                printf("Buffer step requested: %d lines\n", proj->bufstep);
            continue;
        }

        if (strstr(arg, "-of=") == arg) {
            arg += strlen("-of=");
            strcpy(proj->filename, arg);
            printf("Output filename requested: %s\n", proj->filename);
            continue;
        }


        if (strstr(arg, "-input=") == arg) {
            arg += strlen("-input=");

            proj->input_fp = fopen(arg, "r");
            if (proj->input_fp == NULL) {
                printf("Input listing file: \"%s\" not found.\n", arg);
                exit(1);
            }
            printf("Input filelist filename requested: %s\n", arg);
            continue;
        }

        if (strstr(arg, "-fileuse=") == arg) {
            arg += strlen("-fileuse=");
            proj->fileuse_fp = fopen(arg, "w");
            printf("File Use file requested: %s\n", arg);
            continue;
        }

    }

}

void set_buffer_parameters(OUTPUT *proj, PROCESSINFO *process) {
    int nbytes, ib;

    nbytes = 0;

    for (ib = 0; ib < NBANDS; ib++)
        if (process->band[ib]) {
            nbytes += DFKNTsize(DFNT_INT16); /* Each individual band */
            if (process->TOA) nbytes += DFKNTsize(DFNT_INT16); /* TOA reflectance */
        }

    nbytes += DFKNTsize(DFNT_INT16); /* Input file SDS */

    if (process->angles ||
            process->composite == MINSENZ)
        nbytes += DFKNTsize(DFNT_INT16); /* View zenith angle */

    if (process->angles)
        nbytes += 2 * DFKNTsize(DFNT_INT16); /* Solar zenith angle and relative azimuth angle */

    if (process->smoke)
        nbytes += DFKNTsize(DFNT_INT16); /* Smoke index */

    if (process->ndvi ||
            process->composite == MAXNDVI)
        nbytes += DFKNTsize(DFNT_INT16); /* NDVI */

    if (process->evi)
        nbytes += DFKNTsize(DFNT_INT16); /* EVI */

    if (process->gemi)
        nbytes += DFKNTsize(DFNT_INT16); /* GEMI */

    if (process->atmcor)
        proj->buflines = MIN(proj->h,
            (int) ((1e6 * process->maxmem - NLPSURF * NPPSURF * DFKNTsize(DFNT_INT16)) / proj->w / nbytes));
    else
        proj->buflines = MIN(proj->h,
            (int) (1e6 * process->maxmem / proj->w / nbytes));

    if (proj->buflines < proj->h) {

        printf("Output projection has to be bufferized: %d lines out of %d\n", proj->buflines, proj->h);

        /* bufstep is the shift (in # of lines) to apply to the buffer when pixel is out of the current buffer */
        /* If bufstep is too small the buffer will be updated too often */
        /* If bufstep is too big the buffer will be flushed before all lines are */
        /* processed, and the buffer will go back and forth in the output projection */

        if (proj->bufstep == -1) {
            proj->bufstep = proj->buflines / 5 + 1;
            printf("No buffer step requested. %d will be used.\n", proj->bufstep);
        } else {
            if (proj->bufstep > proj->h - proj->buflines) {
                proj->bufstep = proj->h - proj->buflines;
                printf("Buffer step adjusted to %d lines\n", proj->bufstep);
            }
        }
        if (process->gzip) {
            printf("Compression cannot be performed on bufferized output SDS. Request ignored.\n");
            process->gzip = FALSE;
        }

    }

}

void check_process_parameters(PROCESSINFO *process, char **wl) {
    int ib;

    printf("Requested bands: ");
    for (ib = 0; ib < NBANDS; ib++)
        if (process->band[ib])
            break;
    if (ib == NBANDS) /* No band was specified - Process all bands */
        for (ib = 0; ib < NBANDS; ib++)
            process->band[ib] = TRUE;
    for (ib = 0; ib < NBANDS; ib++)
        if (process->band[ib])
            printf("%s ", wl[ib]);
    printf("\n");

    switch (process->composite) {
    case MINBLUE: sprintf(process->compositename, "minimum band %d (%s nm)", BLUEBAND + 1, wl[BLUEBAND]);
        break;
    case MAXNDVI: strcpy(process->compositename, "maximum NDVI");
        break;
    case MINSENZ: strcpy(process->compositename, "minimum view angle");
        break;
    case LASTIN: strcpy(process->compositename, "last in on top");
        break;
    default: fprintf(stderr, "Invalid composite criterion (%d).\n", process->composite);
        exit(-1);
    }
    printf("Compositing criterion: %s\n", process->compositename);

    if (!process->angles &&
            process->composite == MINSENZ)
        printf("View angle will be computed.\n");

    printf("\n");

    if (process->composite == MINBLUE &&
            !process->band[BLUEBAND]) {
        fprintf(stderr, "Cannot make composite without band %d (%s nm)\n", BLUEBAND + 1, wl[BLUEBAND]);
        exit(-1);
    }

    if (process->smoke) {
        if (!process->band[BAND1] ||
                !process->band[BAND2] ||
                !process->band[BAND4] ||
                !process->band[BAND5] ||
                !process->band[BAND6]) {
            fprintf(stderr, "Cannot compute smoke index without bands %s %s %s %s %s\n",
                    wl[BAND1], wl[BAND2], wl[BAND4], wl[BAND5], wl[BAND6]);
            exit(-1);
        }
    }

    if (process->evi) {
        if (!process->band[BAND2] ||
                !process->band[BAND6] ||
                !process->band[BAND8]) {
            fprintf(stderr, "Cannot compute EVI without bands %s, %s and %s\n", wl[BAND2], wl[BAND6], wl[BAND8]);
            exit(-1);
        }
    }

    if (process->ndvi ||
            process->composite == MAXNDVI) {
        if (!process->band[BAND6] ||
                !process->band[BAND8]) {
            fprintf(stderr, "Cannot compute NDVI without bands %s and %s\n", wl[BAND6], wl[BAND8]);
            exit(-1);
        }
    }

    if (process->gemi) {
        if (!process->band[BAND6] ||
                !process->band[BAND8]) {
            fprintf(stderr, "Cannot compute GEMI without bands %s and %s\n", wl[BAND6], wl[BAND8]);
            exit(-1);
        }
    }

    if (!process->band[BAND8]) {
        fprintf(stderr, "Cannot estimate cloudiness without band %s\n", wl[BAND8]);
        exit(-1);
    }

}

void init_output_projection(OUTPUT *proj, PROCESSINFO process) {
    int ib, j;

    for (ib = 0; ib < NBANDS; ib++) {
        proj->refl[ib] = NULL;
        if (process.band[ib]) {
            alloc_new_array((void *) &proj->refl[ib], proj->buflines * proj->w, DFNT_INT16, process.maxmem);
            for (j = 0; j < proj->buflines * proj->w; j++)
                proj->refl[ib][j] = FILL_INT16;
        }
    }

    for (ib = 0; ib < NBANDS; ib++) {
        proj->reflTOA[ib] = NULL;
        if (process.TOA &&
                process.band[ib]) {
            alloc_new_array((void *) &proj->reflTOA[ib], proj->buflines * proj->w, DFNT_INT16, process.maxmem);
            for (j = 0; j < proj->buflines * proj->w; j++)
                proj->reflTOA[ib][j] = FILL_INT16;
        }
    }

    alloc_new_array((void *) &proj->ifno, proj->buflines * proj->w, DFNT_INT16, process.maxmem);
    for (j = 0; j < proj->buflines * proj->w; j++)
        proj->ifno[j] = FILL_INT16;

    if (process.angles ||
            process.composite == MINSENZ) {
        alloc_new_array((void *) &proj->senz, proj->buflines * proj->w, DFNT_INT16, process.maxmem);
        for (j = 0; j < proj->buflines * proj->w; j++)
            proj->senz[j] = FILL_INT16;
    }

    if (process.angles) {
        alloc_new_array((void *) &proj->solz, proj->buflines * proj->w, DFNT_INT16, process.maxmem);
        alloc_new_array((void *) &proj->phi, proj->buflines * proj->w, DFNT_INT16, process.maxmem);
        for (j = 0; j < proj->buflines * proj->w; j++)
            proj->solz[j] = FILL_INT16;
        for (j = 0; j < proj->buflines * proj->w; j++)
            proj->phi[j] = FILL_INT16;
    }

    if (process.smoke) {
        alloc_new_array((void *) &proj->smoke, proj->buflines * proj->w, DFNT_INT16, process.maxmem);
        for (j = 0; j < proj->buflines * proj->w; j++)
            proj->smoke[j] = FILL_INT16;
    }

    if (process.ndvi ||
            process.composite == MAXNDVI) {
        alloc_new_array((void *) &proj->ndvi, proj->buflines * proj->w, DFNT_INT16, process.maxmem);
        for (j = 0; j < proj->buflines * proj->w; j++)
            proj->ndvi[j] = FILL_INT16;
    }

    if (process.evi) {
        alloc_new_array((void *) &proj->evi, proj->buflines * proj->w, DFNT_INT16, process.maxmem);
        for (j = 0; j < proj->buflines * proj->w; j++)
            proj->evi[j] = FILL_INT16;
    }

    if (process.gemi) {
        alloc_new_array((void *) &proj->gemi, proj->buflines * proj->w, DFNT_INT16, process.maxmem);
        for (j = 0; j < proj->buflines * proj->w; j++)
            proj->gemi[j] = FILL_INT16;
    }

}

void init_output_file(OUTPUT *proj, INPUT input, PROCESSINFO process) {
    char path[MAXLENGTH];
    float64 attr;

    proj->startday = input.doy;

    /* Open output L3 HDF file */
    if (proj->filename[0] == 0) {
        sprintf(path, "%s", OUT_DIR);
        sprintf(proj->filename, "%s/%4.4d%3.3d.hdf", path, input.year, proj->startday);
    }
    printf("Creating %s ...\n", proj->filename);
    if ((proj->sd_id = SDstart(proj->filename, DFACC_CREATE)) == -1) {
        fprintf(stderr, "Cannot create file %s\n", proj->filename);
        exit(-1);
    }

    if (SDsetattr(proj->sd_id, "Projection", DFNT_CHAR8, (int) strlen(proj->projname) + 1, proj->projname) == -1)
        fprintf(stderr, "Cannot write attribute \"%s\" in file %s\n", "Projection", proj->filename);

    attr = proj->pixsz * R / RAD2DEG;
    if (SDsetattr(proj->sd_id, "Pixel size (meters)", DFNT_FLOAT64, 1, &attr) == -1)
        fprintf(stderr, "Cannot write attribute \"%s\" in file %s\n", "Pixel size", proj->filename);

    if (SDsetattr(proj->sd_id, "Maximum sensor zenith angle (degrees)", DFNT_FLOAT32, 1, &process.maxsenz) == -1)
        fprintf(stderr, "Cannot write attribute \"%s\" in file %s\n", "Maximum sensor zenith angle", proj->filename);

    if (SDsetattr(proj->sd_id, "Compositing criterion", DFNT_CHAR8,
            (int) strlen(process.compositename) + 1, process.compositename) == -1)
        fprintf(stderr, "Cannot write attribute \"%s\" in file %s\n", "Compositing criterion", proj->filename);

    if ((proj->p1 != 0) || (proj->l1 != 0)) {
        if (SDsetattr(proj->sd_id, "Subset Scan Range", DFNT_INT32, 2, &proj->l1) == -1)
            fprintf(stderr, "Cannot write attribute \"%s\" in file %s\n", "Subset Scan Range", proj->filename);
        if (SDsetattr(proj->sd_id, "Subset Pixel Range", DFNT_INT32, 2, &proj->p1) == -1)
            fprintf(stderr, "Cannot write attribute \"%s\" in file %s\n", "Subset Pixel Range", proj->filename);
    }

    write_cmd_line(proj->sd_id, proj->cmd_line);


    SDend(proj->sd_id);

}

int project_scan(int iscan, char **wl, INPUT *input, OUTPUT *proj, PROCESSINFO process) {
    int ipix, k, iline1, iline2, ipix1, ipix2;
    int iline_out, ipix_out;
    double xproj, yproj;
    double xarr[6], yarr[6];
    double xarrmin, xarrmax, yarrmin, yarrmax;

    if (FAST) {
        for (ipix = 0; ipix < input->Npixels; ipix += input->step) {
            if (input->solz[ipix] == FILL_INT16) continue;
            if (get_xyproj((double) input->lon[ipix], (double) input->lat[ipix], &xproj, &yproj,
                    proj->projection, proj->lon_center) == -1) continue;
            iline_out = (proj->Nl - 1) / 2. - yproj / proj->pixsz + 0.5;
            ipix_out = (proj->Np - 1) / 2. + xproj / proj->pixsz + 0.5;
            if (iline_out >= proj->l1 &&
                    iline_out <= proj->l2 &&
                    ipix_out >= proj->p1 &&
                    ipix_out <= proj->p2) {
                process_pixel(ipix, iline_out, ipix_out, wl, *input, proj, process);
                input->intersect = TRUE;
            }
        }
    }

    if (get_geolocation_scan(iscan + input->step - 2 * input->step *
            (iscan + input->step > input->Nscans - 1),
            input, process) == -1)
        return -1;

    if (FAST) return -1;

    memcpy(input->lon1, input->lon2, (input->Npixels + 1) * sizeof (double));
    memcpy(input->lat1, input->lat2, (input->Npixels + 1) * sizeof (double));
    resample_geo_scan(input->Npixels, input->lon, input->lat, input->lon2, input->lat2);

    if (input->old_iscan != iscan - input->step) {
        for (ipix = 0; ipix < input->Npixels + 1; ipix += input->step) {
            input->lontop[ipix] = -1.5 * input->lon1[ipix] + 0.5 * input->lon2[ipix];
            if (input->lon1[ipix] - input->lon2[ipix] > OEL_PI) input->lontop[ipix] += 0.5 * TWO_PI;
            if (input->lon2[ipix] - input->lon1[ipix] > OEL_PI) input->lontop[ipix] -= 1.5 * TWO_PI;
            input->lattop[ipix] = -1.5 * input->lat1[ipix] + 0.5 * input->lat2[ipix];
            if (get_xyproj(input->lontop[ipix], input->lattop[ipix], &input->xtop[ipix],
                    &input->ytop[ipix], proj->projection, proj->lon_center) == -1) continue;
        }
    } else {
        memcpy(input->lontop, input->lonbot, (input->Npixels + 1) * sizeof (double));
        memcpy(input->lattop, input->latbot, (input->Npixels + 1) * sizeof (double));
        memcpy(input->xtop, input->xbot, (input->Npixels + 1) * sizeof (double));
        memcpy(input->ytop, input->ybot, (input->Npixels + 1) * sizeof (double));
    }

    if (iscan + input->step > input->Nscans - 1) {
        for (ipix = 0; ipix < input->Npixels + 1; ipix += input->step) {
            input->lonbot[ipix] = 1.5 * input->lon1[ipix] - 0.5 * input->lon2[ipix];
            if (input->lon1[ipix] - input->lon2[ipix] > OEL_PI) input->lonbot[ipix] -= 0.5 * TWO_PI;
            if (input->lon2[ipix] - input->lon1[ipix] > OEL_PI) input->lonbot[ipix] += 1.5 * TWO_PI;
            input->latbot[ipix] = 1.5 * input->lat1[ipix] - 0.5 * input->lat2[ipix];
            if (get_xyproj(input->lonbot[ipix], input->latbot[ipix], &input->xbot[ipix],
                    &input->ybot[ipix], proj->projection, proj->lon_center) == -1) continue;
        }
    } else {
        for (ipix = 0; ipix < input->Npixels + 1; ipix += input->step) {
            input->lonbot[ipix] = 0.5 * input->lon1[ipix] + 0.5 * input->lon2[ipix];
            if (fabs(input->lon1[ipix] - input->lon2[ipix]) > OEL_PI)
                input->lonbot[ipix] += 0.5 * TWO_PI;
            input->latbot[ipix] = 0.5 * input->lat1[ipix] + 0.5 * input->lat2[ipix];
            if (get_xyproj(input->lonbot[ipix], input->latbot[ipix], &input->xbot[ipix],
                    &input->ybot[ipix], proj->projection, proj->lon_center) == -1) continue;
        }
    }

    if (fabs(input->tilt[iscan] - input->tilt[iscan - 1 + 2 * (iscan == 0)]) > 0.1)
        for (ipix = 0; ipix < input->Npixels + 1; ipix += input->step)
            get_xyproj(input->lon1[ipix], input->lat1[ipix], &input->xtop[ipix],
                &input->ytop[ipix], proj->projection, proj->lon_center);

    if (fabs(input->tilt[iscan] - input->tilt[iscan + 1 - 2 * (iscan == input->Nscans - 1)]) > 0.1)
        for (ipix = 0; ipix < input->Npixels + 1; ipix += input->step)
            get_xyproj(input->lon1[ipix], input->lat1[ipix], &input->xbot[ipix],
                &input->ybot[ipix], proj->projection, proj->lon_center);


    /* Update output projection */

    for (ipix = 0; ipix < input->Npixels; ipix += input->step) {
        if (input->solz[ipix] == FILL_INT16) continue;
        if (proj->lon1 != UNDEF &&
                ((input->latbot[ipix] > proj->lat1 &&
                input->latbot[ipix + input->step] > proj->lat1) ||
                (input->lattop[ipix] < proj->lat2 &&
                input->lattop[ipix + input->step] < proj->lat2) ||
                (input->lontop[ipix] < proj->lon1 &&
                input->lontop[ipix + input->step] < proj->lon1 &&
                input->lonbot[ipix] < proj->lon1 &&
                input->lonbot[ipix + input->step] < proj->lon1) ||
                (input->lontop[ipix] > proj->lon2 &&
                input->lontop[ipix + input->step] > proj->lon2 &&
                input->lonbot[ipix] > proj->lon2 &&
                input->lonbot[ipix + input->step] > proj->lon2)))
            continue;
        xarr[0] = input->xtop[ipix];
        yarr[0] = input->ytop[ipix];
        xarr[1] = input->xtop[ipix + input->step];
        yarr[1] = input->ytop[ipix + input->step];
        xarr[2] = input->xbot[ipix + input->step];
        yarr[2] = input->ybot[ipix + input->step];
        xarr[3] = input->xbot[ipix];
        yarr[3] = input->ybot[ipix];
        xarr[4] = xarr[0];
        yarr[4] = yarr[0];
        xarr[5] = xarr[1];
        yarr[5] = yarr[1];

        xarrmin = DBL_MAX;
        xarrmax = -DBL_MAX;
        yarrmin = DBL_MAX;
        yarrmax = -DBL_MAX;
        for (k = 0; k < 4; k++) {
            if (xarr[k] < xarrmin) xarrmin = xarr[k];
            if (xarr[k] > xarrmax) xarrmax = xarr[k];
            if (yarr[k] < yarrmin) yarrmin = yarr[k];
            if (yarr[k] > yarrmax) yarrmax = yarr[k];
        }

        if (xarrmin > proj->xmax ||
                xarrmax < proj->xmin ||
                yarrmin > proj->ymax ||
                yarrmax < proj->ymin)
            continue;

        for (k = 0; k < 4; k++)
            if (CROSS_PRODUCT(xarr[k + 1], yarr[k + 1],
                    xarr[k], yarr[k],
                    xarr[k + 2], yarr[k + 2]) < 0 ||
                    (xarr[k] * xarr[k + 1] < 0 &&
                    fabs(xarr[k] - xarr[k + 1]) > 20 * proj->pixsz))
                break;
        if (k < 4) continue;

        iline1 = (proj->Nl - 1) / 2. - yarrmax / proj->pixsz + 0.5;
        iline2 = (proj->Nl - 1) / 2. - yarrmin / proj->pixsz + 0.5;
        ipix1 = (proj->Np - 1) / 2. + xarrmin / proj->pixsz + 0.5;
        ipix2 = (proj->Np - 1) / 2. + xarrmax / proj->pixsz + 0.5;

        for (iline_out = iline1; iline_out <= iline2; iline_out++)
            for (ipix_out = ipix1; ipix_out <= ipix2; ipix_out++) {
                if (iline_out < proj->l1 ||
                        iline_out > proj->l2 ||
                        ipix_out < proj->p1 ||
                        ipix_out > proj->p2)
                    continue;
                xproj = (ipix_out - (proj->Np - 1) / 2.) * proj->pixsz;
                yproj = ((proj->Nl - 1) / 2. - iline_out) * proj->pixsz;
                for (k = 0; k < 4; k++) {
                    if (CROSS_PRODUCT(xarr[k + 1], yarr[k + 1],
                            xproj, yproj,
                            xarr[k], yarr[k]) > EPSILON ||
                            CROSS_PRODUCT(xarr[k + 1], yarr[k + 1],
                            xproj, yproj,
                            xarr[k + 2], yarr[k + 2]) < -EPSILON) break;
                }
                if (k >= 4) {
                    process_pixel(ipix, iline_out, ipix_out, wl, *input, proj, process);
                    input->intersect = TRUE;
                }
            }
    }

    return 0;

}

typedef struct {
    char sds_name[H4_MAX_NC_NAME];
    int32 sds_id, num_type;
    int32 start[H4_MAX_VAR_DIMS], edges[H4_MAX_VAR_DIMS];
    float32 slope, intercept;
    void *data;
} SDS;

#define L2FIELDS (2 * NBANDS + 9)

int read_l2(int iscan, INPUT *input, char **wl) {
    static int32 sds_index, rank, num_type, nattrs;
    static int32 start[H4_MAX_VAR_DIMS], edges[H4_MAX_VAR_DIMS], dim_sizes[H4_MAX_VAR_DIMS];
    static int16 *sola = NULL, *sena = NULL;
    float phi;
    int ib, j, k;
    static char first_time = TRUE;
    static SDS L2sds[L2FIELDS];

    if (first_time) {

        for (ib = 0; ib < NBANDS; ib++) {
            k = ib;
            sprintf(L2sds[k].sds_name, "rhos_%s", wl[ib]);
            L2sds[k].data = input->rho[ib];
            ;
            L2sds[k].num_type = DFNT_FLOAT32;
        }
        for (ib = 0; ib < NBANDS; ib++) {
            k = ib + NBANDS;
            sprintf(L2sds[k].sds_name, "rhot_%s", wl[ib]);
            L2sds[k].data = input->rhoTOA[ib];
            L2sds[k].num_type = DFNT_FLOAT32;
        }

        k++;
        strcpy(L2sds[k].sds_name, "l2_flags");
        L2sds[k].data = input->l2_flags;
        L2sds[k].num_type = DFNT_INT32;

        k++;
        strcpy(L2sds[k].sds_name, "longitude");
        L2sds[k].data = input->lon;
        L2sds[k].num_type = DFNT_FLOAT32;
        k++;
        strcpy(L2sds[k].sds_name, "latitude");
        L2sds[k].data = input->lat;
        L2sds[k].num_type = DFNT_FLOAT32;
        k++;
        strcpy(L2sds[k].sds_name, "ndvi");
        L2sds[k].data = input->ndvi;
        L2sds[k].num_type = DFNT_INT16;
        k++;
        strcpy(L2sds[k].sds_name, "evi");
        L2sds[k].data = input->evi;
        L2sds[k].num_type = DFNT_INT16;
        k++;
        strcpy(L2sds[k].sds_name, "solz");
        L2sds[k].data = input->solz;
        L2sds[k].num_type = DFNT_INT16;
        k++;
        strcpy(L2sds[k].sds_name, "senz");
        L2sds[k].data = input->senz;
        L2sds[k].num_type = DFNT_INT16;
        k++;
        strcpy(L2sds[k].sds_name, "sola");
        L2sds[k].data = NULL;
        L2sds[k].num_type = DFNT_INT16;
        k++;
        strcpy(L2sds[k].sds_name, "sena");
        L2sds[k].data = NULL;
        L2sds[k].num_type = DFNT_INT16;

        first_time = FALSE;

    }

    if (input->newfile) {

        if (input->phi) {
            sena = (int16 *) realloc(sena, input->Npixels * sizeof (int16));
            sola = (int16 *) realloc(sola, input->Npixels * sizeof (int16));
            L2sds[L2FIELDS - 2].data = sola;
            L2sds[L2FIELDS - 1].data = sena;
        }
        for (k = 0; k < L2FIELDS; k++) {
            if (L2sds[k].data) {
                if ((sds_index = SDnametoindex(input->sd_id, L2sds[k].sds_name)) == -1) {
                    fprintf(stderr, "Cannot locate SDS %s %d\n", L2sds[k].sds_name, k);
                    return -1;
                }
                if ((L2sds[k].sds_id = SDselect(input->sd_id, sds_index)) == -1) {
                    fprintf(stderr, "Cannot select SDS %s\n", L2sds[k].sds_name);
                    return -1;
                }
                if (SDgetinfo(L2sds[k].sds_id, L2sds[k].sds_name, &rank, dim_sizes, &num_type, &nattrs) < 0) {
                    fprintf(stderr, "can't get info for SDS %s\n", L2sds[k].sds_name);
                    return -1;
                }
                if (strcmp(L2sds[k].sds_name, "latitude") == 0 ||
                        strcmp(L2sds[k].sds_name, "longitude") == 0) {
                    L2sds[k].slope = 1 / RAD2DEG;
                    L2sds[k].intercept = 0;
                } else if (strcmp(L2sds[k].sds_name, "l2_flags") == 0) {
                    L2sds[k].slope = 0;
                    L2sds[k].intercept = 0;
                } else if (L2sds[k].num_type == DFNT_FLOAT32) {
                    L2sds[k].slope = 1.;
                    L2sds[k].intercept = 0;
                } else if (read_attr(L2sds[k].sds_id, "slope", &L2sds[k].slope) < 0 ||
                        read_attr(L2sds[k].sds_id, "intercept", &L2sds[k].intercept) < 0) {
                    fprintf(stderr, "can't read scaling parameters of SDS %s\n", L2sds[k].sds_name);
                    return -1;
                }
                if (num_type != L2sds[k].num_type) {
                    fprintf(stderr, "SDS %s has the wrong data type (%d)\n", L2sds[k].sds_name, num_type);
                    return -1;
                }
                if (dim_sizes[0] != input->Nscans ||
                        dim_sizes[1] != input->Npixels) {
                    fprintf(stderr, "SDS %s has the wrong dimensions (%dx%d)\n", L2sds[k].sds_name, input->Npixels, input->Nscans);
                    return -1;
                }
            }
        }

    }

    start[0] = iscan;
    start[1] = 0;
    edges[0] = 1;
    edges[1] = input->Npixels;
    for (k = 0; k < L2FIELDS; k++) {
        if (L2sds[k].data) {
            if (SDreaddata(L2sds[k].sds_id, start, NULL, edges, L2sds[k].data) < 0) {
                fprintf(stderr, "can't read SDS %s\n", L2sds[k].sds_name);
                return -1;
            }
            if (k < NBANDS) /* Surface Reflectance: int16 -> float32 */
                for (j = input->Npixels - 1; j >= 0; j--)
                    ((float32 *) L2sds[k].data)[j] = L2sds[k].slope * ((float32 *) L2sds[k].data)[j] + L2sds[k].intercept;
            else if (k > 2 * NBANDS && k < L2FIELDS - 6) /* TOA refl, longitude and latitude: float32 -> float32 */
                for (j = 0; j < input->Npixels; j++)
                    ((float32 *) L2sds[k].data)[j] = L2sds[k].slope * ((float32 *) L2sds[k].data)[j] + L2sds[k].intercept;
        }
    }

    if (input->ndvi) /* NDVI, EVI, angles: int16 -> int16 */
        for (j = 0; j < input->Npixels; j++)
            input->ndvi[j] = floor((L2sds[L2FIELDS - 6].slope * input->ndvi[j] + L2sds[L2FIELDS - 6].intercept) / sf_vi + 0.5);
    if (input->evi)
        for (j = 0; j < input->Npixels; j++)
            input->evi[j] = floor((L2sds[L2FIELDS - 5].slope * input->evi[j] + L2sds[L2FIELDS - 5].intercept) / sf_vi + 0.5);
    if (input->solz)
        for (j = 0; j < input->Npixels; j++)
            input->solz[j] = (L2sds[L2FIELDS - 4].slope * input->solz[j] + L2sds[L2FIELDS - 4].intercept) / sf_angle + 0.5;
    if (input->senz)
        for (j = 0; j < input->Npixels; j++)
            input->senz[j] = (L2sds[L2FIELDS - 3].slope * input->senz[j] + L2sds[L2FIELDS - 3].intercept) / sf_angle + 0.5;
    if (input->phi)
        for (j = 0; j < input->Npixels; j++) {
            phi = L2sds[L2FIELDS - 2].slope * sola[j] + L2sds[L2FIELDS - 2].intercept
                    - L2sds[L2FIELDS - 1].slope * sena[j] - L2sds[L2FIELDS - 1].intercept;
            if (phi > +180) phi -= 360;
            if (phi < -180) phi += 360;
            input->phi[j] = floor(phi / sf_angle + 0.5);
        }
    return 0;

}
