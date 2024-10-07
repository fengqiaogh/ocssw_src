/* ============================================================================================== */
/* aph.c - functions to compute phytoplankton-specific absorption                                 */
/* ============================================================================================== */

#include "l12_proto.h"
#include "giop.h"

/* ---------------------------------------------------------------------------------------------- */
/* aph_bricaud() - compute aph for input wavelength and chl using Bricaud                         */

/* ---------------------------------------------------------------------------------------------- */
float aph_bricaud_1995(float wave, float chl) {
    static int firstCall = 1;
    static float *wtab;
    static float *atab, *datab;
    static float *btab, *dbtab;
    static int ntab = 0;

    float a, b, aph;

    if (firstCall) {

        FILE *fp;
        char *filedir;
        char filename[FILENAME_MAX];
        char line [80];
        int32_t itab;
        float a, b, w;

        if ((filedir = getenv("OCDATAROOT")) == NULL) {
            printf("-E- %s: OCDATAROOT env variable undefined.\n", __FILE__);
            exit(1);
        }
        strcpy(filename, filedir);
        strcat(filename, "/common/aph_bricaud_1995.txt");

        if ((fp = fopen(filename, "r")) == NULL) {
            fprintf(stderr,
                    "-E- %s line %d: unable to open %s for reading\n", __FILE__, __LINE__, filename);
            exit(1);
        }

        printf("Reading aph* from %s.\n", filename);

        // number of lines
        ntab = 0;
        while (fgets(line, 80, fp))
            ntab++;
        rewind(fp);

        // allocate space
        if ((wtab = (float *) calloc(ntab, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for Bricaud aph.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((atab = (float *) calloc(ntab, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for Bricaud aph.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((btab = (float *) calloc(ntab, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for Bricaud aph.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((datab = (float *) calloc(ntab, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for Bricaud aph.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((dbtab = (float *) calloc(ntab, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for Bricaud aph.\n",
                    __FILE__, __LINE__);
            exit(1);
        }

        // read into arrays, skipping header or comment lines
        itab = 0;
        while (itab < ntab) {
            if (fgets(line, 80, fp) == NULL) {
                fprintf(stderr, "-E- %s line %d: error reading %s at line %d\n", __FILE__, __LINE__, filename, itab);
                exit(1);
            }
            if (line[0] != '/' && line[0] != '!') {
                sscanf(line, "%f,%f,%f", &w, &a, &b);
                wtab [itab] = w;
                atab [itab] = a;
                btab [itab] = b;
                itab++;
            }
        }
        ntab = itab;

        // precompute derivatives for spline interpolation
        spline(wtab, atab, ntab, 1e30, 1e30, datab);
        spline(wtab, btab, ntab, 1e30, 1e30, dbtab);

        firstCall = 0;
    }

    if (wave > wtab[ntab - 1]) wave = wtab[ntab - 1];

    // interpolate coefficients to wavelength
    splint(wtab, atab, datab, ntab, wave, &a);
    splint(wtab, btab, dbtab, ntab, wave, &b);
    aph = a * pow(chl, -b);

    return (aph);
}

/* ---------------------------------------------------------------------------------------------- */
/* aph_bricaud() - compute aph for input wavelength and chl using Bricaud                         */

/* ---------------------------------------------------------------------------------------------- */
float aph_bricaud_1998(float wave, float chl) {
    static int firstCall = 1;
    static float *table;
    static float *wtab;
    static float *aphitab, *daphitab;
    static float *ephitab, *dephitab;
    static int ntab = 0;

    float aph, aphi, ephi;

    if (firstCall) {

        char *filedir;
        char filename[FILENAME_MAX];
        int ncol, nrow;

        if ((filedir = getenv("OCDATAROOT")) == NULL) {
            printf("-E- %s: OCDATAROOT env variable undefined.\n", __FILE__);
            exit(1);
        }
        strcpy(filename, filedir);
        strcat(filename, "/common/aph_bricaud_1998.txt");
        printf("Reading aph* from %s.\n", filename);

        nrow = table_row_count(filename);
        if (nrow <= 0) {
            printf("-E- %s line %d: error opening (%s) file", __FILE__, __LINE__, filename);
            exit(1);
        }

        ncol = table_column_count(filename);
        table = table_read_r4(filename, ncol, nrow);
        ntab = nrow;

        wtab = &table[nrow * 0];
        aphitab = &table[nrow * 3];
        ephitab = &table[nrow * 4];

        // precompute derivatives for spline interpolation
        if ((daphitab = (float *) calloc(ntab, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for Bricaud aph.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((dephitab = (float *) calloc(ntab, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for Bricaud aph.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        spline(wtab, aphitab, ntab, 1e30, 1e30, daphitab);
        spline(wtab, ephitab, ntab, 1e30, 1e30, dephitab);

        firstCall = 0;
    }

    if (wave > wtab[ntab - 1]) wave = wtab[ntab - 1];

    // interpolate coefficients to wavelength
    splint(wtab, aphitab, daphitab, ntab, wave, &aphi);
    splint(wtab, ephitab, dephitab, ntab, wave, &ephi);
    aph = aphi * pow(chl, ephi - 1);

    return (aph);
}

float aph_bricaud(float wave, float chl) {
    return (aph_bricaud_1998(wave, chl));
}

/* ---------------------------------------------------------------------------------------------- */
/* aph_bricaud() - compute aph for input wavelength and chl using Bricaud                         */

/* ---------------------------------------------------------------------------------------------- */
float pderiv_aph_bricaud_1998(float wave, float chl) {
    static int firstCall = 1;
    static float *table;
    static float *wtab;
    static float *aphitab, *daphitab;
    static float *ephitab, *dephitab;
    static int ntab = 0;

    float aphi, ephi, daphdchl;

    if (firstCall) {

        char *filedir;
        char filename[FILENAME_MAX];
        int ncol, nrow;

        if ((filedir = getenv("OCDATAROOT")) == NULL) {
            printf("-E- %s: OCDATAROOT env variable undefined.\n", __FILE__);
            exit(1);
        }
        strcpy(filename, filedir);
        strcat(filename, "/common/aph_bricaud_1998.txt");
        //printf("Reading aph* from %s.\n", filename);

        nrow = table_row_count(filename);
        if (nrow <= 0) {
            printf("-E- %s line %d: error opening (%s) file", __FILE__, __LINE__, filename);
            exit(1);
        }

        ncol = table_column_count(filename);
        table = table_read_r4(filename, ncol, nrow);
        ntab = nrow;

        wtab = &table[nrow * 0];
        aphitab = &table[nrow * 3];
        ephitab = &table[nrow * 4];

        // precompute derivatives for spline interpolation
        if ((daphitab = (float *) calloc(ntab, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for Bricaud aph.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((dephitab = (float *) calloc(ntab, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for Bricaud aph.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        spline(wtab, aphitab, ntab, 1e30, 1e30, daphitab);
        spline(wtab, ephitab, ntab, 1e30, 1e30, dephitab);

        firstCall = 0;
    }

    if (wave > wtab[ntab - 1]) wave = wtab[ntab - 1];

    // interpolate coefficients to wavelength
    splint(wtab, aphitab, daphitab, ntab, wave, &aphi);
    splint(wtab, ephitab, dephitab, ntab, wave, &ephi);
    //aph = aphi * pow(chl, ephi - 1);

    daphdchl = aphi*(ephi-1.)* pow(chl,ephi-2.);

    return (daphdchl);
}

float aph_bricaud_pderiv(float wave, float chl) {
    return (pderiv_aph_bricaud_1998(wave, chl));
}



/* ---------------------------------------------------------------------------------------------- */
/* aph_ciotti() - compute aph for input wavelength and size fraction using Cioti                  */

/* ---------------------------------------------------------------------------------------------- */
float aph_ciotti(float wave, float sf) {
    static int firstCall = 1;
    static float *wtab;
    static float *ptab, *dptab; // pico
    static float *mtab, *dmtab; // micro
    static int ntab = 0;

    float pico, micro, aph;

    if (firstCall) {

        FILE *fp;
        char *filedir;
        char filename[FILENAME_MAX];
        char line [80];
        int32_t itab;
        float wv;

        if ((filedir = getenv("OCDATAROOT")) == NULL) {
            printf("-E- %s: OCDATAROOT env variable undefined.\n", __FILE__);
            exit(1);
        }
        strcpy(filename, filedir);
        strcat(filename, "/common/aph_ciotti_2002_2006.txt");

        if ((fp = fopen(filename, "r")) == NULL) {
            fprintf(stderr,
                    "-E- %s line %d: unable to open %s for reading\n", __FILE__, __LINE__, filename);
            exit(1);
        }

        printf("Reading aph* from %s.\n", filename);

        // number of lines
        ntab = 0;
        while (fgets(line, 80, fp))
            ntab++;
        rewind(fp);

        // allocate space
        if ((wtab = (float *) calloc(ntab, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for Ciotti aph.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((ptab = (float *) calloc(ntab, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for Ciotti aph.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((mtab = (float *) calloc(ntab, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for Ciotti aph.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((dptab = (float *) calloc(ntab, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for Ciotti aph.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((dmtab = (float *) calloc(ntab, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for Ciotti aph.\n",
                    __FILE__, __LINE__);
            exit(1);
        }

        // read into arrays, skipping header or comment lines
        itab = 0;
        while (itab < ntab) {
            if (fgets(line, 80, fp) == NULL) {
                fprintf(stderr, "-E- %s line %d: error reading %s at line %d\n", __FILE__, __LINE__, filename, itab);
                exit(1);
            }
            if (line[0] != '/' && line[0] != '!') {
                sscanf(line, "%f %f %f", &wv, &pico, &micro);
                wtab [itab] = wv;
                ptab [itab] = pico * 0.0230 / 0.891;
                mtab [itab] = micro * 0.0086 / 1.249;
                itab++;
            }
        }
        ntab = itab;

        // precompute derivatives for spline interpolation
        spline(wtab, ptab, ntab, 1e30, 1e30, dptab);
        spline(wtab, mtab, ntab, 1e30, 1e30, dmtab);

        firstCall = 0;
    }

    if (wave > wtab[ntab - 1]) wave = wtab[ntab - 1];

    // interpolate coefficients to wavelength
    splint(wtab, ptab, dptab, ntab, wave, &pico);
    splint(wtab, mtab, dmtab, ntab, wave, &micro);

    aph = sf * pico + (1.0 - sf) * micro;

    return (aph);
}


/* ---------------------------------------------------------------------------------------------- */
/* get_aphstar() - compute aph* for center wavelength, wavelength width, function type, and proxy */

/* ---------------------------------------------------------------------------------------------- */
float get_aphstar(float wave, int dwave, int ftype, float proxy) {
    int32_t npts;
    float wave1, wave2, wv;
    float aph = 0.0;

    if (dwave > 1) {
        npts = dwave + 1;
        wave1 = wave - dwave / 2;
        wave2 = wave + dwave / 2;
    } else {
        npts = 1;
        wave1 = wave;
        wave2 = wave;
    }

    for (wv = wave1; wv <= wave2; wv += 1.0) {
        switch (ftype) {
        case APHBRICAUD:
            aph += aph_bricaud(wv, proxy);
            break;
        case APHCIOTTI:
            aph += aph_ciotti(wv, proxy);
            break;
        }
    }

    aph /= npts;

    return (aph);
}

/* ---------------------------------------------------------------------------------------------- */
/* get_uaphstar() -  compute aph partial derivative wrt chl for center wavelength,                */
/*                    wavelength width type, and proxy                                            */

/* ---------------------------------------------------------------------------------------------- */
float get_aphstar_pderiv(float wave, int dwave, int ftype, float proxy) {
    int32_t npts;
    float wave1, wave2, wv;
    float daphdchl = 0.0;

    if (dwave > 1) {
        npts = dwave + 1;
        wave1 = wave - dwave / 2;
        wave2 = wave + dwave / 2;
    } else {
        npts = 1;
        wave1 = wave;
        wave2 = wave;
    }

    for (wv = wave1; wv <= wave2; wv += 1.0) {
        switch (ftype) {
        case APHBRICAUD:
            daphdchl += aph_bricaud_pderiv(wv, proxy);
            break;
        case APHCIOTTI:
            daphdchl += 0;  //set to zero at present
            break;
        }
    }

    daphdchl /= npts;

    return (daphdchl);
}
