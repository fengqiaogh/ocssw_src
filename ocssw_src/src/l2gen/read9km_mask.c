/****************************************************************/
/*								*/
/*    Reads a binary mask file used in excluding 9km bins	*/
/*    from L3 binned data					*/
/*								*/
/*    result = read9km_mask(char *file, char *mask)	        */
/*								*/
/*    result = 0 - success					*/
/*	       1 - failure					*/
/*								*/
/* Programmer     Organization    Date     Description of change*/
/* -------------- ------------  --------   ---------------------*/
/* Ewa Kwiatkowska   SAIC    11 March 2004  Original development*/
/*								*/
/****************************************************************/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define PI  3.141592653589793


static int32_t *numbin_xkm, *basebin_xkm, *numbin_9km, *basebin_9km;
static float *latbin_xkm;

int read9km_mask(char *file, char *mask) {
    char data[742553];
    FILE *fp;
    int32_t len, i, j;
    unsigned char x[8];

    if ((fp = fopen(file, "r")) == NULL) {
        printf(" Cannot open the mask file %s", file);
        return (1);
    }
    if ((len = (int32_t) fread((void *) data, sizeof (char), 1536, fp)) != 1536) {
        printf(" Cannot read the mask file %s", file);
        fclose(fp);
        return (1);
    }
    if ((len = (int32_t) fread((void *) data, sizeof (char), 742553, fp)) != 742553) {
        printf(" Cannot read the mask file %s", file);
        fclose(fp);
        return (1);
    }
    fclose(fp);

    for (j = 0; j < 8; j++) x[j] = (unsigned char) pow(2.0, (double) j);

    for (i = 0; i < 742553; i++) {
        for (j = 0; j < 8; j++) if (data[i] & x[j]) mask[i * 8 + j] = 1;
            else mask[i * 8 + j] = 0;
    }


    numbin_xkm = NULL;
    basebin_xkm = NULL;
    latbin_xkm = NULL;
    numbin_9km = NULL;
    basebin_9km = NULL;

    return ( 0);

}

int is_masked(int32_t bin, char *mask, int32_t nrows) {

    static int first_check = 1;

    int32_t i, rlow, rhi, rmid, row, rows_9km = 2160, xbin, bin_row, bin_col;
    float radfac, lats, lat, lon;


    if (nrows == 2160) {

        if (mask[bin] > 0) return ( 1);
        else return ( 0);

    } else

        if (nrows > 0) {

        radfac = PI / 180.0;

        if (first_check) {
            /* ----------------- */
            /* Set up bin arrays */
            /* ----------------- */
            if ((numbin_xkm = (int32_t *) malloc(nrows * sizeof (int32_t))) == NULL) {
                printf("-E- %s line %d: Cannot allocate memory to numbin", __FILE__, __LINE__);
                return (-1);
            }
            if ((basebin_xkm = (int32_t *) malloc((nrows + 1) * sizeof (int32_t))) == NULL) {
                printf("-E- %s line %d: Cannot allocate memory to basebin", __FILE__, __LINE__);
                return (-1);
            }
            if ((latbin_xkm = (float *) malloc(nrows * sizeof (float))) == NULL) {
                printf("-E- %s line %d: Cannot allocate memory to basebin", __FILE__, __LINE__);
                return (-1);
            }

            for (i = 0; i < nrows; i++) {
                *(latbin_xkm + i) = (i + 0.5) * (180.0 / nrows) - 90.0;
                *(numbin_xkm + i) = (int32_t) (cos(*(latbin_xkm + i) * radfac) * (2.0 * nrows) + 0.5);
            }

            *basebin_xkm = 1;
            for (i = 1; i < nrows; i++) *(basebin_xkm + i) = *(basebin_xkm + i - 1) + *(numbin_xkm + i - 1);
            basebin_xkm[nrows] = *(basebin_xkm + nrows - 1) + *(numbin_xkm + nrows - 1) - 1;


            if ((numbin_9km = (int32_t *) malloc(rows_9km * sizeof (int32_t))) == NULL) {
                printf("-E- %s line %d: Cannot allocate memory to numbin", __FILE__, __LINE__);
                return (-1);
            }
            if ((basebin_9km = (int32_t *) malloc((rows_9km + 1) * sizeof (int32_t))) == NULL) {
                printf("-E- %s line %d: Cannot allocate memory to basebin", __FILE__, __LINE__);
                return (-1);
            }
            for (i = 0; i < rows_9km; i++) {
                lats = (i + 0.5) * (180.0 / rows_9km) - 90.0;
                numbin_9km[i] = (int32_t) (cos(lats * radfac) * (2.0 * rows_9km) + 0.5);
            }

            *basebin_9km = 1;
            for (i = 1; i < rows_9km; i++) *(basebin_9km + i) = *(basebin_9km + i - 1) + *(numbin_9km + i - 1);
            basebin_9km[rows_9km] = *(basebin_9km + rows_9km - 1) + *(numbin_9km + rows_9km - 1) - 1;
        }


        if (bin < 1)
            bin = 1; /* south pole */
        if (bin > basebin_xkm[nrows])
            bin = basebin_xkm[nrows]; /* north pole */

        /* binary search for row in range [1..nrows] */
        rlow = 1; /* 1-relative */
        rhi = nrows; /* 1-relative */
        while (1) {

            rmid = (rlow + rhi - 1) / 2; /* 0-relative */
            if (*(basebin_xkm + rmid) > bin) rhi = rmid;
            else rlow = rmid + 1;

            if (rlow == rhi) {
                row = rlow;
                break;
            }
        }

        lat = *(latbin_xkm + row - 1);
        lon = 360.0 * (bin - *(basebin_xkm + row - 1) + 0.5) / *(numbin_xkm + row - 1);
        lon = lon - 180.0; /* note, lon returned here may be in 0 to 360 */

        bin_row = (int32_t) ((90.0 + lat) * ((double) rows_9km / (double) 180.0));
        bin_col = (int32_t) ((double) numbin_9km[bin_row] * (lon + 180.0) / (double) 360.0);
        xbin = basebin_9km[bin_row] + bin_col;

        first_check = 0;

        if (mask[xbin] > 0) return ( 1);
        else return ( 0);

    } else {

        if (numbin_xkm != NULL) free(numbin_xkm);
        if (basebin_xkm != NULL) free(basebin_xkm);
        if (latbin_xkm != NULL) free(latbin_xkm);
        if (numbin_9km != NULL) free(numbin_9km);
        if (basebin_9km != NULL) free(basebin_9km);

        return ( 0);
    }

}




