/* =========================================================== */
/* Module read_target.c                                        */
/*                                                             */
/* Functions to open and read a recalibration target file.     */
/*                                                             */
/* Written By:                                                 */
/*                                                             */
/*     B. A. Franz                                             */
/*     SAIC General Sciences Corp.                             */
/*     NASA/SIMBIOS Project                                    */
/*     April 1998                                              */
/*                                                             */
/* =========================================================== */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "l12_proto.h"
#include "filehdr_struc.h"
/*
#include "target_struc.h"
#include "filehandle.h"
 */
#include "read_l3bin.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_spline.h>

static FILE *fp = NULL;
static char current_file[FILENAME_MAX];


/* ----------------------------------------------------------- */
/* close_target() - close recal target file                    */

/* ----------------------------------------------------------- */
void close_target(void) {
    fclose(fp);
}


/* ----------------------------------------------------------- */
/* read_targethdr() - loads header info                        */

/* ----------------------------------------------------------- */
int read_targethdr(filehandle *file) {
    filehdr hdr;

    fseek(fp, 0, SEEK_SET);

    if (fread(&hdr, 1, sizeof (hdr), fp) != sizeof (hdr)) {
        printf("-E- %s: File read error.\n", __FILE__);
        return (1);
    }

    if (hdr.length < 0 || hdr.length > 1000000 ||
            hdr.npix < 0 || hdr.npix > 10000) {
        printf("-E- %s: Header values out of range.\n", __FILE__);
        printf("  Record length    = %d\n", hdr.length);
        printf("  Pixels per scan  = %d\n", hdr.npix);
        return (1);
    }

    file->sensorID = hdr.sensorID;
    file->length = hdr.length;
    file->npix = hdr.npix;
    file->format = hdr.format;
    file->nscan = hdr.nscan;
    file->mode = READ;

    return (0);
}


/* ----------------------------------------------------------- */
/* open_target() - opens file if not already opened            */

/* ----------------------------------------------------------- */
int open_target(filehandle *file) {
    if (fp == NULL || strcmp(file->name, current_file) != 0) {
        if (fp != NULL) close_target();
        if ((fp = fopen(file->name, "r")) == NULL) {
            printf("-E- %s: Error opening %s for reading.\n",
                    __FILE__, file->name);
            return (1);
        }
        strcpy(current_file, file->name);
        if (read_targethdr(file) != 0) {
            printf("-E- %s: Error reading header for %s.\n",
                    __FILE__, file->name);
            return (1);
        }
    }

    return (0);
}


/* ----------------------------------------------------------- */
/* read_target() - reads one recal target record               */
/*                                                             */
/* B. A. Franz, GSC, SIMBIOS Project, March 1998               */

/* ----------------------------------------------------------- */
int read_target(filehandle *file, int32_t recnum, tgstr *tgrec) {
    /*                                                         */
    /* Open the input file if it is not already open           */
    /*                                                         */
    if (open_target(file) != 0) {
        printf("-E- %s: File open error.\n", __FILE__);
        return (1);
    }

    if (feof(fp)) {
        printf("-I- %s: End of target file %s reached.",
                __FILE__, file->name);
        return (1);
    }

    if (fseek(fp, (recnum + 1) * file->length, SEEK_SET) != 0) {
        printf("-E- %s: Error seeking record %d in %s.",
                __FILE__, recnum, file->name);
        return (1);
    }

    if (fread(tgrec->data, 1, file->length, fp) != file->length) {
        return (1);
    }

    tgrec->sensorID = file->sensorID;
    tgrec->length = file->length;
    tgrec->npix = file->npix;

    return (0);
}

/**
 * bin_match - a binary search routine to find the nearest bin in a list of bin numbers.
 * TODO: replace this function with code Don adding to libbin++
 *
 * @param nbins - number of bins
 * @param bins - list of bin numbers
 * @param bin_num - bin number to seach
 * @return
 */
int32_t bin_match(int32_t nbins, int32_t *bins, int32_t bin_num) {
    int32_t jl = -1, ju = nbins, jm = 0;
    int32_t ascnd;

    ascnd = (bins[nbins - 1] >= bins[0]);
    while (ju - jl > 1) {
        jm = (ju + jl) / 2;
        if (ascnd == (bin_num >= bins[jm]))
            jl = jm;
        else
            ju = jm;
    }

    if (bin_num == bins[jl]) return (jl);
    if (jl + 1 < nbins && bin_num == bins[jl + 1]) return (jl + 1);
    if (jl > 0 && bin_num == bins[jl - 1]) return (jl - 1);
    if (bin_num == bins[0]) return (0);
    if (bin_num == bins[nbins - 1]) return (nbins - 1);

    return (-1);
}

/**
 * lonlat2bin - returns a L3 bin number given a longitude and latitude
 * TODO: replace this function with code Don adding to libbin++
 * @param l3bin - l3bin structure
 * @param lon - longitude
 * @param lat - latitude
 * @return
 */
int32_t lonlat2bin(l3binstr *l3bin, float lon, float lat) {
    int32_t row;
    int32_t col;
    int32_t bin;

    lat = MAX(MIN(lat, 90), -90);
    lon = MAX(MIN(lon, 180), -180);

    row = MIN(((90 + lat) * l3bin->nrows / 180.0), l3bin->nrows - 1);
    col = MIN(((lon + 180.0) * l3bin->numbin[row] / 360.0), l3bin->numbin[row]);

    bin = l3bin->basebin[row] + col;

    return (bin);
}


/* ----------------------------------------------------------- */
/* read_target_l3() - loads one recal target record            */
/*                                                             */
/* B. A. Franz, GSC, SIMBIOS Project, March 1998               */
/* ----------------------------------------------------------- */

/**
 * read_target_l3 - loads a vicarous calibration target record,
 * applying band shifting if necessary
 *
 * @param file - filehandle structure for L3 bin target file
 * @param l1rec
 * @param nbands
 * @param tgrec
 * @return
 */
int read_target_l3(filehandle *file, l1str *l1rec, int32_t nbands, tgstr *tgrec) {
    static int firstCall = 1;
    static l3binstr l3bin;
    static int l3nwaves;
    static float *l3wvls;
    static double *l3wvlinterp;
    static int *wvlmatch;
    static int nvisbands;

    int32_t ip, ib, ipb, band;
    int32_t bin, idx, nobs, nscenes;
    float chl, aot;
    static int needBandShift = 0;
    int l3nbands;
    static float *l3vals;
    static double *l3valsinterp;
    static gsl_interp_accel *acc;
    static gsl_spline *spline_steffen;
    float interpband; // for linear and spline interp

    if (firstCall) {
        l3nbands = rdsensorinfo(tgrec->sensorID, 0, "fwave", (void **) &l3wvls);

        read_l3bin(file->name, &l3bin, l3nbands);
        l3nwaves = l3bin.nwave;
        if ((l3bin.sensorID != l1rec->l1file->sensorID) && (input->band_shift_opt != 2)) {
            needBandShift = 1;
        }
        wvlmatch = (int *) malloc(l1rec->l1file->nbands * sizeof (int));

        for (ib = 0; ib < l1rec->l1file->nbands; ib++) {
            wvlmatch[ib] = windex(l1rec->l1file->fwave[ib], l3bin.wavelengths, l3nwaves);
            if (l1rec->l1file->iwave[ib] < 700) {
                nvisbands++;
            }
        }
        l3wvlinterp = (double *) calloc(l3nwaves, sizeof (double));
        for (band = 0; band < l3nwaves; band++) {
            l3wvlinterp[band] = l3bin.wavelengths[band];
        }
        l3valsinterp = (double *) calloc(l3nwaves, sizeof (double));
        l3vals = (float *) calloc(l3nwaves, sizeof (float));

        acc = gsl_interp_accel_alloc();
        spline_steffen = gsl_spline_alloc(gsl_interp_steffen, l3nwaves);
        firstCall = 0;
    }

    for (ip = 0; ip < l1rec->npix; ip++) {
        // initialize target rec
        for (ib = 0; ib < nbands; ib++) {
            ipb = ip * nbands + ib;
            tgrec->nLw[ipb] = BAD_FLT;
            tgrec->Lw[ipb] = BAD_FLT;
        }
        tgrec->solz[ip] = BAD_FLT;

        // check L1 masking
        if (input->vcal_depth < 0) {
            if (l1rec->dem[ip] > input->vcal_depth)
                continue;
        }
        // locate bin 
        bin = lonlat2bin(&l3bin, l1rec->lon[ip], l1rec->lat[ip]);
        idx = bin_match(l3bin.nbins, l3bin.bins, bin);

        // check masking and copy values 
        if (idx >= 0) {
            chl = l3bin.chl[idx];
            aot = l3bin.tau[idx];
            nobs = l3bin.nobs[idx];
            nscenes = l3bin.nscenes[idx];
            for (band = 0; band < l3nwaves; band++) {
                l3valsinterp[band] = l3bin.data[idx][band];
                l3vals[band] = l3bin.data[idx][band];
                if (!l3bin.hasRrs) {
                    l3vals[band] /= l1rec->Fo[ib];
                    l3valsinterp[band] /= l1rec->Fo[ib];
                } // Need Rrs for bandshift

            }

            gsl_spline_init(spline_steffen, l3wvlinterp, l3valsinterp, l3nwaves);

            if (chl <= input->chlthreshold && aot <= input->aotthreshold
                    && nobs >= input->vcal_min_nbin && nscenes >= input->vcal_min_nscene) {

                for (ib = 0; ib < nvisbands; ib++) {
                    ipb = ip * nbands + ib;
                    // Only bandshift if L2 wvl is different from L3 wvl even if sensors differ
                    if (needBandShift && (l3bin.wavelengths[wvlmatch[ib]] != l1rec->l1file->fwave[ib])) {
                        // option 0: linear interpolation
                        if (input->band_shift_opt == 0 ||  l1rec->l1file->fwave[ib] > 700) {
                            if (l1rec->l1file->fwave[ib] > l3bin.wavelengths[0] &&
                                l1rec->l1file->fwave[ib] < l3bin.wavelengths[l3nwaves-1]) {
                            // interpband = linterp(l3wvls, l3vals, l3nwaves, l1rec->l1file->fwave[ib]);
                                interpband = gsl_spline_eval(spline_steffen, l1rec->l1file->fwave[ib], acc);
                            } else {
                                // don't extrapolate, just use the l3 value as if we didn't need to shift
                                interpband = l3bin.data[idx][wvlmatch[ib]];
                            }
                        // option 1: bio-optical band shift
                        } else if (input->band_shift_opt == 1 && l1rec->l1file->fwave[ib] <= 700) {
                            interpband = bioBandShift(l3bin.wavelengths, l3vals, chl, l3nwaves, l1rec->l1file->fwave[ib]);
                        } else {
                            printf("-E- %s:%d: band_shift_opt error.\n", __FILE__, __LINE__);
                            exit(EXIT_FAILURE);
                        }

                        interpband *= l1rec->Fo[ib];
                        tgrec->nLw[ipb] = interpband;

                        //the nLws in the target files are now in W/m2/um/sr, so a factor of 10 too big
                        if (l3bin.hasRrs == 0)
                            tgrec->nLw[ipb] /= 10.;
                    } else {
                        if (l3bin.hasRrs == 1)
                            tgrec->nLw[ipb] = l3bin.data[idx][wvlmatch[ib]] * l1rec->Fo[ib];
                        else
                            tgrec->nLw[ipb] = l3bin.data[idx][wvlmatch[ib]] / 10.;
                    }
                }
                for (ib = nvisbands; ib < nbands; ib++) {
                    ipb = ip * nbands + ib;
                    tgrec->nLw[ipb] = 0.0;
                }
            }
        }
    }

    tgrec->sensorID = l1rec->l1file->sensorID;
    tgrec->npix = l1rec->npix;

    return (0);
}
