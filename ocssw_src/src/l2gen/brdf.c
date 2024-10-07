#include "l12_proto.h"
#include "freearray.h"

static float radeg = RADEG;
static float nw = 1.334;

void foq_morel(int foqopt, l2str *l2rec, float wave[], int32_t nwave, float chl,
        float nLw[], float Fo[], float solz, float senz, float phi, float brdf[]);
void dtran_brdf(l2str *l2rec, int32_t ip, float wave[], int32_t nwave, float Fo[], float nLw[], float chl,
        float brdf[]);
void diff_tran_corr_(int *iphase, float *solz, float *senz, float *phi,
        float *chl, float *taua, float *correct);
int32_t getncDimensionLength(int ncid, int dimId);
int getncFQdim(int ncid, char *file, char *sdsname, int sds_id, int nexp, float *bdat);

/* ---------------------------------------------------------------------------- */
/* ocbrdf - bidirectional reflectance correction for one ocean pixel            */
/*                                                                              */
/* Input:                                                                       */
/* sensorID - sensor identification (for chlorophyll algorithm)                 */
/* brdf_opt - brdf option (bitmask, 1=fresnel(sen), 2=fresnel(sol), 4=f/Q)      */
/* wave     - list of sensor wavelengths (no missing-band gaps)                 */
/* nwave    - number of sensor wavelengths to compute brdf                      */
/* solz     - solar zenith (deg)                                                */
/* senz     - sensor zenith (deg)                                               */
/* phi      - relative azimuth (deg, Gordon definition)                         */
/* nLw      - list of normalized water-leaving radiances, band indexed          */
/* Fo       - list of solar irradinaces, band indexed (mW/cm^2/um)              */
/*                                                                              */
/* Output:                                                                      */
/* brdf     - list of brdf corrections to multiply nLw, band indexed            */
/*                                                                              */
/* Note: for brdf_opt=1, 2, or 3, nLw and Fo can be set to NULL pointer         */
/*                                                                              */
/* Written By:                                                                  */
/* B. Franz, OBPG, October 2004                                                 */

/* ---------------------------------------------------------------------------- */

int ocbrdf(l2str *l2rec, int32_t ip, int32_t brdf_opt, float wave[], int32_t nwave,
        float solz, float senz, float phi, float ws, float chl, float nLw[], float Fo[], float brdf[]) {
    static int firstCall = 1;

    float *temp;
    int32_t status = 0;
    int32_t iw;

    if (firstCall == 1) {
        firstCall = 0;
        if (brdf_opt > 0)
            printf("\nApplying ocean BRDF including:\n");
        if ((brdf_opt & FRESNSEN) > 0)
            printf("    Reflection/refraction for upwelling radiance.\n");
        if ((brdf_opt & FRESNSOL) > 0)
            printf("    Reflection/refraction for downwelling radiance.\n");
        if ((brdf_opt & FOQMOREL) > 0)
            printf("    Morel f/Q\n");
        if ((brdf_opt & QMOREL) > 0)
            printf("    Morel Q\n");
        if ((brdf_opt & DTBRDF) > 0)
            printf("    Gordon diffuse transmittance effect.\n");
    }

    if ((temp = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("-E- : Error allocating memory to temp\n");
        exit(FATAL_ERROR);
    }

    for (iw = 0; iw < nwave; iw++)
        brdf[iw] = 1.0;

    /* transmittance of view path through air & sea interface */
    if ((brdf_opt & FRESNSEN) > 0) {
        float tf = fresnel_sen(senz, 0);
        for (iw = 0; iw < nwave; iw++)
            brdf[iw] *= tf;
    }

    /* transmittance of solar path through air & sea interface */
    if ((brdf_opt & FRESNSOL) > 0) {
        fresnel_sol(wave, nwave, solz, ws, temp, 0);
        for (iw = 0; iw < nwave; iw++)
            brdf[iw] *= temp[iw];
    }

    /* Morel f/Q correction */
    if ((brdf_opt & FOQMOREL) > 0) {
        foq_morel(FOQMOREL, l2rec, wave, nwave, chl, nLw, Fo, solz, senz, phi, temp);
        for (iw = 0; iw < nwave; iw++) {
            brdf[iw] *= temp[iw];
        }
    }

    /* Morel f/Q with free sun (Voss) */
    if ((brdf_opt & QMOREL) > 0) {
        foq_morel(QMOREL, l2rec, wave, nwave, chl, nLw, Fo, solz, senz, phi, temp);
        for (iw = 0; iw < nwave; iw++) {
            brdf[iw] *= temp[iw];
        }
    }

    /* Gordon correction of diffuse transmittance */
    if ((brdf_opt & DTBRDF) > 0) {
        dtran_brdf(l2rec, ip, wave, nwave, Fo, nLw, chl, temp);
        for (iw = 0; iw < nwave; iw++) {
            brdf[iw] *= temp[iw];
        }
    }

    free(temp);

    return (status);
}

/* ---------------------------------------------------------------------------- */
/* fresnel_sen() - effects of the air-sea transmittance for sensor view         */
/*                                                                              */
/* Description:                                                                 */
/*   This computes effects of the air-sea transmittance (depending on sensor    */
/*   zenith angle) on the derived normalized water-leaving radiance.            */
/*   Menghua Wang 5/27/02.                                                      */
/*                                                                              */
/* modified to return fresnel transmittance as option, December 2008, BAF       */

/* ---------------------------------------------------------------------------- */

float fresnel_sen(float senz, int return_tf) {
    static float tf0 = 0.9795218;

    float mu = cos(senz / radeg);
    float sq, r2, q1, fres, tf, brdf;

    sq = sqrt(nw * nw - 1. + mu * mu);
    r2 = pow((mu - sq) / (mu + sq), 2);
    q1 = (1. - mu * mu - mu * sq) / (1. - mu * mu + mu * sq);
    fres = r2 * (q1 * q1 + 1.) / 2.0;
    tf = 1. - fres;
    brdf = tf0 / tf;

    if (return_tf != 0)
        return (tf);
    else
        return (brdf);
}


/* ---------------------------------------------------------------------------- */
/* fresnel_sol() - effects of the air-sea transmittance for solar path          */
/*                                                                              */
/* Description:                                                                 */
/*   This computes the correction factor on normalized water-leaving radiance   */
/*   to account for the solar zenith angle effects on the downward irradiance   */
/*   from above the ocean surface to below the surface.                         */
/*   Menghua Wang 9/27/04.                                                      */
/*                                                                              */
/* Added windspeed dependence, December 2004, BAF                               */
/* Modified to return air-sea transmittance as option, December 2008, BAF       */

/* ---------------------------------------------------------------------------- */
void fresnel_sol(float wave[], int32_t nwave, float solz, float ws, float brdf[], int return_tf) {

#define NTWAVE 6
#define NTWIND 5

    static int firstCall = 1;
    static int *tindx;
    static float *tf0;
    static float twave [] = {412., 443., 490., 510., 555., 670.};
    static float tsigma[] = {0.0, 0.1, 0.2, 0.3, 0.4};

    /* M Wang, personal communication, red-nir iterpolated */
    static float tf0_w[] = {412., 443., 490., 510., 555., 670., 765., 865.};
    static float tf0_v[] = {0.965980, 0.968320, 0.971040, 0.971860, 0.973450, 0.977513, 0.980870, 0.984403};

    static float c[NTWIND][4][NTWAVE] = {
        { /* ws=0.0 */
            { -0.0087, -0.0122, -0.0156, -0.0163, -0.0172, -0.0172}, /* a */
            { 0.0638, 0.0415, 0.0188, 0.0133, 0.0048, -0.0003}, /* b */
            { -0.0379, -0.0780, -0.1156, -0.1244, -0.1368, -0.1430}, /* c */
            { -0.0311, -0.0427, -0.0511, -0.0523, -0.0526, -0.0478} /* d */
        },
        { /* ws=1.9 */
            { -0.0011, -0.0037, -0.0068, -0.0077, -0.0090, -0.0106}, /* a */
            { 0.0926, 0.0746, 0.0534, 0.0473, 0.0368, 0.0237}, /* b */
            { -5.3E-4, -0.0371, -0.0762, -0.0869, -0.1048, -0.1260}, /* c */
            { -0.0205, -0.0325, -0.0438, -0.0465, -0.0506, -0.0541} /* d */
        },
        { /* ws=7.5 */
            { 6.8E-5, -0.0018, -0.0011, -0.0012, -0.0015, -0.0013}, /* a */
            { 0.1150, 0.1115, 0.1075, 0.1064, 0.1044, 0.1029}, /* b */
            { 0.0649, 0.0379, 0.0342, 0.0301, 0.0232, 0.0158}, /* c */
            { 0.0065, -0.0039, -0.0036, -0.0047, -0.0062, -0.0072} /* d */
        },
        { /* ws=16.9 */
            { -0.0088, -0.0097, -0.0104, -0.0106, -0.0110, -0.0111}, /* a */
            { 0.0697, 0.0678, 0.0657, 0.0651, 0.0640, 0.0637}, /* b */
            { 0.0424, 0.0328, 0.0233, 0.0208, 0.0166, 0.0125}, /* c */
            { 0.0047, 0.0013, -0.0016, -0.0022, -0.0031, -0.0036} /* d */
        },
        { /* ws=30 */
            { -0.0081, -0.0089, -0.0096, -0.0098, -0.0101, -0.0104}, /* a */
            { 0.0482, 0.0466, 0.0450, 0.0444, 0.0439, 0.0434}, /* b */
            { 0.0290, 0.0220, 0.0150, 0.0131, 0.0103, 0.0070}, /* c */
            { 0.0029, 0.0004, -0.0017, -0.0022, -0.0029, -0.0033} /* d */
        }
    };

    float x, x2, x3, x4;
    int is, is1, is2;
    int iw, i;
    float sigma;
    float slp;
    float brdf1, brdf2;

    /* on first call, find closest table entry to each input wavelength */
    if (firstCall == 1) {
        firstCall = 0;
        if ((tindx = (int *) calloc(nwave, sizeof (int))) == NULL) {
            printf("-E- : Error allocating memory to tindx\n");
            exit(FATAL_ERROR);
        }
        if ((tf0 = (float *) calloc(nwave, sizeof (float))) == NULL) {
            printf("-E- : Error allocating memory to tf0\n");
            exit(FATAL_ERROR);
        }

        for (iw = 0; iw < nwave; iw++) {
            tindx[iw] = windex(wave[iw], twave, NTWAVE);
            tf0 [iw] = linterp(tf0_w, tf0_v, 8, wave[iw]);
        }
    }

    sigma = 0.0731 * sqrt(MAX(ws, 0.0));

    x = log(cos(MIN(solz, 80.0) / radeg));
    x2 = x*x;
    x3 = x*x2;
    x4 = x*x3;

    /* find bracketing table winds */
    for (is = 0; is < NTWIND; is++)
        if (tsigma[is] > sigma)
            break;
    is2 = MIN(is, NTWIND - 1);
    is1 = is2 - 1;
    slp = (sigma - tsigma[is1]) / (tsigma[is2] - tsigma[is1]);

    /* compute at bounding windspeeds and interpolate */
    for (iw = 0; iw < nwave; iw++) {
        i = tindx[iw];
        brdf1 = 1.0 + c[is1][0][i] * x + c[is1][1][i] * x2
                + c[is1][2][i] * x3 + c[is1][3][i] * x4;
        brdf2 = 1.0 + c[is2][0][i] * x + c[is2][1][i] * x2
                + c[is2][2][i] * x3 + c[is2][3][i] * x4;
        brdf[iw] = brdf1 + slp * (brdf2 - brdf1);
        if (return_tf != 0)
            brdf[iw] = tf0[iw] / brdf[iw];
    }
}


/* ---------------------------------------------------------------------------- */
/* gothic_R - returns total effect of transmission through air-sea interface    */

/* ---------------------------------------------------------------------------- */
void gothic_R(float wave[], int32_t nwave, float solz, float senz, float ws, float R[]) {
    int iw;
    float tf = fresnel_sen(senz, 1);

    fresnel_sol(wave, nwave, solz, ws, R, 1);
    for (iw = 0; iw < nwave; iw++)
        R[iw] = R[iw] * tf / nw / nw;

    return;
}


#define N_W  7
#define N_S  6
#define N_C  6
#define N_N 17
#define N_A 13

/* return closest indice of xtab where xtab[i] < x */
int morel_index(float xtab[], int32_t ntab, float x) {
    int i = 0;

    if (x <= xtab[0])
        i = 0;
    else if (x >= xtab[ntab - 1])
        i = ntab - 2;
    else {
        while (x > xtab[i]) i++;
        i--;
    }

    return (i);
}


/* ---------------------------------------------------------------------------- */
/* foqint_morel() - reads and interpolates f/Q tables of Morel & Gentilli       */
/*                                                                              */
/* wave[] - list of input wavelengths (nm)                                      */
/* nwave  - number of input wavelengths                                         */
/* solz   - solar zenith angle of observation (deg)                             */
/* senzp  - view  zenith angle of observation, below surface (deg)              */
/* phi    - relative azimuth of observation, 0=<--> (deg)                       */
/* chl    - chlorophyll-a concentration (mg/m^3)                                */
/* brdf[] - band-indexed array of f/Q corrections per wavelength (f0/Q0)/(f/Q)  */
/*                                                                              */

/* ---------------------------------------------------------------------------- */
void foqint_morel(char *file, float wave[], int32_t nwave, float solz, float senzp,
        float phi, float chl, float brdf[]) {
    static int firstCall = 1;
    static float *foqtab;
    static float *wavetab;
    static float *solztab;
    static float *chltab;
    static float *senztab;
    static float *phitab;
    static float *lchltab;


    float lchl;
    int i, j, k, l, m;
    int iw, js, kc, ln, ma;
    float ds[2], dc[2], dn[2], da[2];

    static int32_t n_a, n_n, n_c, n_s, n_w;

    char sdsname[H4_MAX_NC_NAME];
    int ncid, ndims, nvars, ngatts, unlimdimid;
    int32 sds_id;
    nc_type rh_type; /* variable type */
    int rh_dimids[H4_MAX_VAR_DIMS]; /* dimension IDs */
    int rh_natts; /* number of attributes */
    int rh_ndims; /* number of dims */
    int status, ndx;
    if (firstCall == 1) {

        /* try netCDF first */
        if (nc_open(file, NC_NOWRITE, &ncid) == NC_NOERR) {

            status = nc_inq(ncid, &ndims, &nvars, &ngatts, &unlimdimid);
            if (status != NC_NOERR) {
                fprintf(stderr, "-E- %s line %d:  Error reading %s from %s.\n",
                        __FILE__, __LINE__, sdsname, file);
                exit(1);
            }

            strcpy(sdsname, "foq");

            status = nc_inq_varid(ncid, sdsname, &sds_id);
            if (status != NC_NOERR) {
                fprintf(stderr, "-E- %s line %d:  Error reading %s from %s.\n",
                        __FILE__, __LINE__, sdsname, file);
                exit(1);
            }
            status = nc_inq_var(ncid, sds_id, 0, &rh_type, &rh_ndims, rh_dimids,
                    &rh_natts);

            if (rh_ndims != 5) {
                fprintf(stderr, "-E- %s line %d:  Wrong number of dimensions for %s.  Need 5 got %d.\n",
                        __FILE__, __LINE__, sdsname, rh_ndims);
                exit(1);

            }
            printf("rh_ndims=%d rh_dimids=%d %d %d %d %d \n", rh_ndims, rh_dimids[0], rh_dimids[1], rh_dimids[2], rh_dimids[3], rh_dimids[4]);

            n_w = getncDimensionLength(ncid, rh_dimids[0]);
            n_s = getncDimensionLength(ncid, rh_dimids[1]);
            n_c = getncDimensionLength(ncid, rh_dimids[2]);
            n_n = getncDimensionLength(ncid, rh_dimids[3]);
            n_a = getncDimensionLength(ncid, rh_dimids[4]);

            printf("morel f/q file dimensions n_a=%d n_n=%d n_c=%d n_s=%d n_w=%d \n", n_a, n_n, n_c, n_s, n_w);

            printf("\nReading foq file %s ndims=%d nvars=%d sds_id=%d var=%s\n\n", file, ndims, nvars, sds_id, sdsname);
            /* Read the data. */
            if ((foqtab = (float *) calloc(n_w * n_s * n_c * n_n * n_a, sizeof (float))) == NULL) {
                printf("-E- : Error allocating memory to tindx\n");
                exit(FATAL_ERROR);
            }
            if (nc_get_var(ncid, sds_id, foqtab) != NC_NOERR) {
                fprintf(stderr, "-E- %s line %d:  Error reading %s from %s.\n",
                        __FILE__, __LINE__, sdsname, file);
                exit(1);
            }

            if ((phitab = (float *) calloc(n_a, sizeof (float))) == NULL) {
                printf("-E- : Error allocating memory to tindx\n");
                exit(FATAL_ERROR);
            }
            if ((senztab = (float *) calloc(n_n, sizeof (float))) == NULL) {
                printf("-E- : Error allocating memory to tindx\n");
                exit(FATAL_ERROR);
            }
            if ((solztab = (float *) calloc(n_s, sizeof (float))) == NULL) {
                printf("-E- : Error allocating memory to tindx\n");
                exit(FATAL_ERROR);
            }
            if ((chltab = (float *) calloc(n_c, sizeof (float))) == NULL) {
                printf("-E- : Error allocating memory to tindx\n");
                exit(FATAL_ERROR);
            }

            if ((wavetab = (float *) calloc(n_w, sizeof (float))) == NULL) {
                printf("-E- : Error allocating memory to tindx\n");
                exit(FATAL_ERROR);
            }
            if ((lchltab = (float *) calloc(n_c, sizeof (float))) == NULL) {
                printf("-E- : Error allocating memory to tindx\n");
                exit(FATAL_ERROR);
            }

            strcpy(sdsname, "phi");
            status = getncFQdim(ncid, file, sdsname, sds_id, n_a, phitab);
            strcpy(sdsname, "senz");
            status = getncFQdim(ncid, file, sdsname, sds_id, n_n, senztab);
            strcpy(sdsname, "solz");
            status = getncFQdim(ncid, file, sdsname, sds_id, n_s, solztab);
            strcpy(sdsname, "chl");
            status = getncFQdim(ncid, file, sdsname, sds_id, n_c, chltab);
            strcpy(sdsname, "wave");
            status = getncFQdim(ncid, file, sdsname, sds_id, n_w, wavetab);

            printf("\nClosing foq file %s\n\n", file);

            /* Close the file */
            if (nc_close(ncid) != NC_NOERR) {
                fprintf(stderr, "-E- %s line %d: error closing %s.\n",
                        __FILE__, __LINE__, file);
                exit(1);
            }
        }


        printf("\nMorel f/Q table from file %s\n", file);

        //        for (i=0; i<n_w; i++)
        //          for (j=0; j<n_s; j++)
        //            for (k=0; k<n_c; k++)
        //              for (l=0; l<n_n; l++) {
        //                for (m=0; m<n_a-1; m++)
        //                {
        //                  ndx = i*n_s*n_c*n_n*n_a+j*n_c*n_n*n_a+k*n_n*n_a+l*n_a+m;
        //                  printf("%f",*(foqtab+ndx));
        //                  //printf("%f",foqtab[i][j][k][l][m]);
        //                }
        //                m = n_a-1;
        //                ndx = i*n_s*n_c*n_n*n_a+j*n_c*n_n*n_a+k*n_n*n_a+l*n_a+m;
        //                //printf("%f\n",foqtab[i][j][k][l][m]);
        //                printf("%f\n",*(foqtab+ndx));
        //            }

        /* create table for log of chlorophyll */
        for (kc = 0; kc < n_c; kc++)
            lchltab[kc] = log(chltab[kc]);

        firstCall = 0;
    }


    lchl = log(MAX(chl, 0.01));

    if (senzp < senztab[0]) {
        senzp = senztab[0];
    }

    /* lower bounding indices */
    js = morel_index(solztab, n_s, solz);
    kc = morel_index(lchltab, n_c, lchl);
    ln = morel_index(senztab, n_n, senzp);
    ma = morel_index(phitab, n_a, phi);

    ds[0] = (solztab[js + 1] - solz) / (solztab[js + 1] - solztab[js]);
    ds[1] = (solz - solztab[js]) / (solztab[js + 1] - solztab[js]);

    dc[0] = (lchltab[kc + 1] - lchl) / (lchltab[kc + 1] - lchltab[kc]);
    dc[1] = (lchl - lchltab[kc]) / (lchltab[kc + 1] - lchltab[kc]);

    dn[0] = (senztab[ln + 1] - senzp) / (senztab[ln + 1] - senztab[ln]);
    dn[1] = (senzp - senztab[ln]) / (senztab[ln + 1] - senztab[ln]);

    da[0] = (phitab [ma + 1] - phi) / (phitab [ma + 1] - phitab [ma]);
    da[1] = (phi - phitab [ma]) / (phitab [ma + 1] - phitab [ma]);

    for (iw = 0; iw < nwave; iw++) {

        /* using nearest wavelength (tables are for MERIS bands) */

        i = windex(wave[iw], wavetab, n_w);

        brdf[iw] = 0.0;

        for (j = 0; j <= 1; j++)
            for (k = 0; k <= 1; k++)
                for (l = 0; l <= 1; l++)
                    for (m = 0; m <= 1; m++) {

                        ndx = i * n_s * n_c * n_n * n_a + (js + j) * n_c * n_n * n_a + (kc + k) * n_n * n_a + (ln + l) * n_a + ma + m;
                        brdf[iw] += ds[j] * dc[k] * dn[l] * da[m]*(*(foqtab + ndx));
                    }
    }

    return;
}

/* -----------------------------------------------------------------------*/
/* getDimensionLength() - get the length of the dimension ID              */

/* -----------------------------------------------------------------------*/
int32_t getncDimensionLength(int ncid, int dimId) {
    size_t length;
    int32_t nl;

    if (nc_inq_dimlen(ncid, dimId, &length) != NC_NOERR) {
        char name[H4_MAX_NC_NAME];
        nc_inq_dim(ncid, dimId, name, &length);
        fprintf(stderr,
                "-E- %s line %d: could not get size of dimension \"%s\" in netCDF File.\n",
                __FILE__, __LINE__, name);
        exit(1);
    }
    nl = length;
    return nl;
}

int getncFQdim(int ncid, char *file, char *sdsname, int sds_id, int nexp, float *bdat) {

    int status;
    nc_type rh_type; /* variable type */
    int rh_dimids[H4_MAX_VAR_DIMS]; /* dimension IDs */
    int rh_natts; /* number of attributes */
    int rh_ndims, nsz; /* number of dims */

    status = nc_inq_varid(ncid, sdsname, &sds_id);
    if (status != NC_NOERR) {
        fprintf(stderr, "-E- %s line %d:  Error reading %s from %s.\n",
                __FILE__, __LINE__, sdsname, file);
        exit(1);
    }
    status = nc_inq_var(ncid, sds_id, 0, &rh_type, &rh_ndims, rh_dimids,
            &rh_natts);

    nsz = getncDimensionLength(ncid, rh_dimids[0]);

    if (nsz != nexp) {
        fprintf(stderr, "-E- %s line %d:  Wrong dimensions for %s.  Need %d got %d.\n",
                __FILE__, __LINE__, sdsname, nexp, nsz);
        exit(1);


    }
    if (nc_get_var(ncid, sds_id, bdat) != NC_NOERR) {
        fprintf(stderr, "-E- %s line %d:  Error reading %s from %s.\n",
                __FILE__, __LINE__, sdsname, file);
        exit(1);
    }

    return status;

}
/* ---------------------------------------------------------------------------- */
/* foq_morel() - computes f/Q correction of Morel & Gentilli by iteration       */
/*                                                                              */
/* foqopt - 0=full f/Q, 1=no normalization to sun overhead (fixed f)            */
/* sensorID - MSl12 sensor identification number                                */
/* wave[] - list of input wavelengths (nm)                                      */
/* nwave  - number of input wavelengths                                         */
/* nLw[]  - normalize water-leaving radiances per wave (mW/cm^2/um/sr)          */
/* Fo[]   - solar irradiance per wave (mW/cm^2/um/sr)                           */
/* solz   - solar zenith angle of observation (deg)                             */
/* senzp  - view  zenith angle of observation, below surface (deg)              */
/* phi    - relative azimuth of observation, 0=<--> (deg)                       */
/* brdf[] - band-indexed array of f/Q corrections per wavelength (f0/Q0)/(f/Q)  */

/* ---------------------------------------------------------------------------- */
void foq_morel(int foqopt, l2str *l2rec, float wave[], int32_t nwave, float chl,
        float nLw[], float Fo[], float solz, float senz, float phi, float brdf[]) {
    static int maxiter = 3;
    static int compchl = 1;

    float senzp, phip;
    float *foq0;
    float *foq;
    float *Rrs;
    int32_t iw, iter;
    int numiter;

    if ((foq0 = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("-E- : Error allocating memory to foq0\n");
        exit(FATAL_ERROR);
    }
    if ((foq = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("-E- : Error allocating memory to foq\n");
        exit(FATAL_ERROR);
    }
    if ((Rrs = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("-E- : Error allocating memory to Rrs\n");
        exit(FATAL_ERROR);
    }

    /* Need view zenith and relative azimuth below water. The MSl12 definition of   */
    /* relative azimuth is consistent with the below-water definition. We just need */
    /* to limit from [-180,180] to [0,180].                                         */
    phip = fabs(phi);
    senzp = asin(sin(senz / radeg) / nw) * radeg;

    /* Compute starting chlorophyll (if not supplied) */

    if (chl < 0.0) {
        compchl = 1;
        numiter = maxiter;
        for (iw = 0; iw < nwave; iw++) {
            Rrs[iw] = nLw[iw] / Fo[iw];
        }
        chl = get_default_chl(l2rec, Rrs);
    } else {
        compchl = 0;
        numiter = 1;
    }

    /* If we retrieved a valid chlorophyll, we can compute the correction. */

    if (chl >= 0.0) {

        for (iter = 0; iter < numiter; iter++) {

            if (foqopt == QMOREL)
                foqint_morel(input->fqfile, wave, nwave, solz, 0.0, 0.0, chl, foq0);
            else
                foqint_morel(input->fqfile, wave, nwave, 0.0, 0.0, 0.0, chl, foq0);

            foqint_morel(input->fqfile, wave, nwave, solz, senzp, phip, chl, foq);

            for (iw = 0; iw < nwave; iw++) {
                brdf [iw] = foq0[iw] / foq[iw];
                Rrs [iw] = nLw[iw] * brdf[iw] / Fo[iw];
            }

            if (compchl) {
                chl = get_default_chl(l2rec, Rrs);
            }
        }
    }

    if (chl < 0.0) {
        for (iw = 0; iw < nwave; iw++) {
            brdf [iw] = 1.0;
        }
    }

    free(foq);
    free(foq0);
    free(Rrs);
    return;
}


/* ---------------------------------------------------------------------------- */
/* qint_morel() - reads and interpolates Q tables of Morel & Gentilli           */
/*                                                                              */
/* wave[] - list of input wavelengths (nm)                                      */
/* nwave  - number of input wavelengths                                         */
/* solz   - solar zenith angle of observation (deg)                             */
/*                                                                              */

/* ---------------------------------------------------------------------------- */
void qint_morel(float wave[], int32_t nwave, float solz, float chl, float Qn[]) {
    static int firstCall = 1;

    static float q0tab [N_C][N_W];
    static float sqtab [N_C][N_W];
    static float dq0tab [N_C][N_W];
    static float dsqtab [N_C][N_W];
    static float wavetab[N_W] = {412.5, 442.5, 490.0, 510.0, 560.0, 620.0, 660.0};
    static float chltab [N_C] = {0.03, 0.1, 0.3, 1.0, 3.0, 10.0};
    static float lchltab[N_C];

    int ic, iw;
    float lchl, q1, q2, s1, s2, qn1, qn2;

    if (firstCall == 1) {

        char filename[FILENAME_MAX];
        char *delim = " ";
        FILE *fp;
        char *tmp_str;
        char line [80];
        char buff [80];

        if ((tmp_str = getenv("OCDATAROOT")) == NULL) {
            printf("OCDATAROOT environment variable is not defined.\n");
            exit(1);
        }

        strcpy(filename, tmp_str);
        strcat(filename, "/common/morel_q0.dat");
        fp = fopen(filename, "r");
        if (!fp) {
            printf("-E- %s line %d: error opening (%s) file", __FILE__, __LINE__, filename);
            exit(1);
        }

        printf("\nLoading Morel Q0 table from file %s\n", filename);

        for (ic = 0; ic < N_C; ic++) {
            if (fgets(line, 80, fp) == NULL) {
                fprintf(stderr, "-E- %s line %d: error reading %s at line %d\n", __FILE__, __LINE__, filename, ic);
                exit(1);
            }
            strcpy(buff, line);
            chltab[ic] = atof(strtok(buff, delim));
            for (iw = 0; iw < N_W; iw++) {
                q0tab[ic][iw] = atof(strtok(NULL, delim));
            }
        }

        strcpy(filename, tmp_str);
        strcat(filename, "/common/morel_sq.dat");
        fp = fopen(filename, "r");
        if (!fp) {
            printf("-E- %s line %d: error opening (%s) file", __FILE__, __LINE__, filename);
            exit(1);
        }

        printf("\nLoading Morel Sq table from file %s\n", filename);

        for (ic = 0; ic < N_C; ic++) {
            if (fgets(line, 80, fp) == NULL) {
                fprintf(stderr, "-E- %s line %d: error reading %s at line %d\n", __FILE__, __LINE__, filename, ic);
                exit(1);
            }
            strcpy(buff, line);
            chltab[ic] = atof(strtok(buff, delim));
            for (iw = 0; iw < N_W; iw++) {
                sqtab[ic][iw] = atof(strtok(NULL, delim));
            }
        }

        // create table for log of chlorophyll
        for (ic = 0; ic < N_C; ic++)
            lchltab[ic] = log(chltab[ic]);

        // precompute derivatives for spline interpolation in wavelength

        for (ic = 0; ic < N_C; ic++) {
            spline(wavetab, &q0tab[ic][0], N_W, 1e30, 1e30, &dq0tab[ic][0]);
            spline(wavetab, &sqtab[ic][0], N_W, 1e30, 1e30, &dsqtab[ic][0]);
        }

        firstCall = 0;
    }

    // find bounding chl indices and weights

    lchl = log(MAX(chl, 0.01));
    ic = morel_index(lchltab, N_C, lchl);

    for (iw = 0; iw < nwave; iw++) {

        if (wave[iw] <= wavetab[0]) {
            q1 = q0tab[ic + 0][0];
            q2 = q0tab[ic + 1][0];
            s1 = sqtab[ic + 0][0];
            s2 = sqtab[ic + 1][0];
        } else if (wave[iw] >= wavetab[N_W - 1]) {
            q1 = q0tab[ic + 0][N_W - 1];
            q2 = q0tab[ic + 1][N_W - 1];
            s1 = sqtab[ic + 0][N_W - 1];
            s2 = sqtab[ic + 1][N_W - 1];
        } else {
            splint(wavetab, &q0tab[ic + 0][0], &dq0tab[ic + 0][0], N_W, wave[iw], &q1);
            splint(wavetab, &q0tab[ic + 1][0], &dq0tab[ic + 1][0], N_W, wave[iw], &q2);
            splint(wavetab, &sqtab[ic + 0][0], &dsqtab[ic + 0][0], N_W, wave[iw], &s1);
            splint(wavetab, &sqtab[ic + 1][0], &dsqtab[ic + 1][0], N_W, wave[iw], &s2);
        }

        qn1 = q1 + s1 * (1.0 - cos(solz / radeg));
        qn2 = q2 + s2 * (1.0 - cos(solz / radeg));

        Qn[iw] = qn1 + (lchl - lchltab[ic])*(qn2 - qn1) / (lchltab[ic + 1] - lchltab[ic]);
    }

    return;
}


/* ---------------------------------------------------------------------------- */
/* dtran_brdf() - computes BRDF effect on diffuse transmittance                 */
/*                                                                              */
/* wrapper for H. Gordon's diff_tran_corr() function.                           */
/*                                                                              */
/* l2rec  - pointer to level-2 record (assumes aerosol terms loaded for ip)     */
/* ip     - scan pixel number (from 0)                                          */
/* wave[] - list of input wavelengths (nm)                                      */
/* nwave  - number of input wavelengths                                         */
/* Fo[]   - solar irradiance per wave (mW/cm^2/um/sr)                           */
/* nLw[]  - normalize water-leaving radiances per wave (mW/cm^2/um/sr)          */
/* brdf[] - band-indexed array of correction to nLw                             */
/*                                                                              */
/* B. Franz, September 2006                                                     */

/* ---------------------------------------------------------------------------- */
void dtran_brdf(l2str *l2rec, int32_t ip, float wave[], int32_t nwave, float Fo[], float nLw[], float chl,
        float brdf[]) {
    static int firstCall = 1;
    static int modindex[] = {4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 15, 16};
    static float wavetab [] = {412., 443., 490., 510., 555., 670., 765., 865.};
    static int *waveindex;
    static int iw865;

    int32_t amodmin = l2rec->aermodmin[ip];
    int32_t amodmax = l2rec->aermodmax[ip];
    float amodrat = l2rec->aerratio[ip];
    float solz = l2rec->l1rec->solz[ip];
    float senz = l2rec->l1rec->senz[ip];
    float phi = l2rec->l1rec->delphi[ip];
    float taua865;

    float *Rrs;
    float **cbrdf;
    int iphase;
    float *correct;
    int iw, im;

    if (firstCall == 1) {
        if ((waveindex = (int *) calloc(nwave, sizeof (int))) == NULL) {
            printf("-E- : Error allocating memory to waveindex\n");
            exit(FATAL_ERROR);
        }
        printf("Initializing BRDF correction to diffuse transmittance.\n");
        /* tables are for seawifs wavelength; use nearest */
        for (iw = 0; iw < nwave; iw++) {
            waveindex[iw] = windex(wave[iw], wavetab, 8);
            printf("Using %3.0f entries for sensor wavelength %3.0f nm\n",
                    wavetab[waveindex[iw]], wave[iw]);
        }
        iw865 = waveindex[7];
        firstCall = 0;
    }
    if ((Rrs = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("-E- : Error allocating memory to Rrs\n");
        exit(FATAL_ERROR);
    }
    if ((correct = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("-E- : Error allocating memory to correct\n");
        exit(FATAL_ERROR);
    }
    if ((cbrdf = (float **) calloc(2, sizeof (float *))) == NULL) {
        printf("-E- : Error allocating memory to cbrdf\n");
        exit(FATAL_ERROR);
    }
    for (im = 0; im < 2; im++)
        if ((cbrdf[im] = (float *) calloc(2, sizeof (float))) == NULL) {
            printf("-E- : Error allocating memory to cbrdf\n");
            exit(FATAL_ERROR);
        }

    /* want relative azimuth in radians, other angles in degrees */
    phi /= radeg;

    /* Compute starting chlorophyll (if not supplied) */
    if (chl < 0.0) {
        for (iw = 0; iw < nwave; iw++) {
            Rrs[iw] = nLw[iw] / Fo[iw];
        }
        chl = get_default_chl(l2rec, Rrs);
    }

    taua865 = l2rec->taua[ip * nwave + iw865];

    /* If we retrieved a valid chlorophyll, we can compute the correction. */

    if (chl > 0.0 && taua865 > 0.0) {

        chl = MAX(MIN(chl, 10.0), 0.03);

        for (im = 0; im <= 1; im++) {

            /* 12-model index within 16-model suite (assumes std 12 models) */
            if (im == 0)
                iphase = modindex[amodmin];
            else
                iphase = modindex[amodmax];

            diff_tran_corr_(&iphase, &solz, &senz, &phi, &chl, &taua865, correct);

            /* compute correction per wavelength, per aerosol model */
            for (iw = 0; iw < nwave; iw++) {
                cbrdf[im][iw] = (1 + correct[waveindex[iw]]);
            }
        }

        /* combine corrections from 2 aerosol models */
        for (iw = 0; iw < nwave; iw++) {
            brdf[iw] = 1.0 / ((1 - amodrat) * cbrdf[0][iw] + amodrat * cbrdf[1][iw]);
        }

    } else {
        for (iw = 0; iw < nwave; iw++)
            brdf [iw] = 1.0;
    }

    free(Rrs);
    free(correct);
    freeArrayf(cbrdf, 2);
    return;
}
