/* ---------------------------------------------------------- */
/* filter.c - Filtering control module for MSl12              */
/*                                                            */
/* Written By: B. Franz, SAIC GSC, NASA/SIMBIOS, October 1998 */
/* ---------------------------------------------------------- */

#include "l12_proto.h"
#include <ctype.h>
#include <stdint.h>

static float badval = BAD_FLT;

typedef struct fnode_str {
    int32_t i;
    int32_t j;
    float v;
} fnode;

void fctl_init(fctlstr *fctl) {
    fctl->nfilt = 0;
    fctl->nscan = 1;
    fctl->npix = 1;
}

void get_kernel(filstr *f) {
    int type = 0;
    int32_t i;
    int32_t j;
    int32_t nnx;
    int32_t nx = f->nx;
    int32_t ny = f->ny;

    if (nx == -1) {

        type = 1;
        f->nx = nx = ny;
    }

    if (nx == -2) {
        type = 2;
        f->nx = nx = 1;
    }

    if ((f->kernel = (char *) malloc(nx * ny * sizeof (char))) == NULL) {
        printf("-E- %s %d: Unable to allocate workspace of %d bytes for filter kernel\n",
                __FILE__, __LINE__, nx * ny);
        exit(1);
    }

    for (j = 0; j < ny; j++)
        for (i = 0; i < nx; i++)
            f->kernel[j * nx + i] = 0;

    if (type == 0) {

        /* square */

        for (j = 0; j < ny; j++)
            for (i = 0; i < nx; i++)
                f->kernel[j * nx + i] = 1;
        if (f->minfill <= 0)
            f->minfill = nx * ny / 2 + 1;
        else
            f->minfill = MIN(MAX(f->minfill, 1), nx * ny / 2 + 1);

    } else if (type == 2) {

        /* along-track */

        for (j = 0; j < ny; j++)
            for (i = 0; i < nx; i++)
                f->kernel[j * nx + i] = 1;
        if (f->minfill <= 0)
            f->minfill = nx * ny / 2 + 1;

    } else {

        /* diamond */

        for (j = 0; j < ny; j++) {
            nnx = nx - 2 * abs(ny / 2 - j);
            for (i = nx / 2 - nnx / 2; i <= nx / 2 + nnx / 2; i++)
                f->kernel[j * nx + i] = 1;
        }
        if (f->minfill <= 0)
            f->minfill = nx * ny / 4 + 1;
        else
            f->minfill = MIN(MAX(f->minfill, 1), nx * ny / 4 + 1);
    }
}

int fctl_set(fctlstr *fctl, int32_t npix, char *fname,
        int32_t band, int32_t nx, int32_t ny, int32_t minfill, int32_t nbands) {
    int i = fctl->nfilt;
    int32_t j;
    int32_t func = -1;

    /* filter width can't exceed number of pixels per scan */
    nx = MIN(nx, npix);

    /* force filter window to odd number */
    if (nx >= 0)
        nx = MAX(1, (nx / 2)*2 + 1);
    ny = MAX(1, (ny / 2)*2 + 1);


    if (strstr(fname, "dilate"))
        func = FDILATE;
    else if (strstr(fname, "stlight"))
        func = FSTLIGHT;
    else if (strstr(fname, "clean"))
        func = FCLEAN;
    else if (strstr(fname, "ltrmed"))
        func = FLTRMED;
    else if (strstr(fname, "ltrmean"))
        func = FLTRMEAN;
    else if (strstr(fname, "ltriqmean"))
        func = FLTRIQMEAN;
    else if (strstr(fname, "ltmed"))
        func = FLTMED;
    else if (strstr(fname, "ltmean"))
        func = FLTMEAN;
    else if (strstr(fname, "btdetavg"))
        func = FBTDETAVG;
    else if (strstr(fname, "test"))
        func = FTEST;
    else if (strstr(fname, "epsmean"))
        func = FEPSMEAN;
    else if (strstr(fname, "ltrreject"))
        func = FLTRREJECT;



    switch (func) {
    case FDILATE:
        if (band < 1 || band > L1_NFLAGS) {
            if (want_verbose)
                printf("-E- %s %d: bogus band number %d for filter %s\n",
                    __FILE__, __LINE__, band, fname);
            return (0);
        }
        if (want_verbose)
            printf("Setting %d x %d dilation filter on %s mask\n",
                nx, ny, l2_flag_lname[band - 1]);
        band--;
        break;
    case FSTLIGHT:
        if (band < 1 || band > L1_NFLAGS) {
            if (want_verbose)
                printf("-E- %s %d: bogus band number %d for filter %s\n",
                    __FILE__, __LINE__, band, fname);
            return (0);
        }
        if (want_verbose)
            printf("Setting %d x %d straylight filter on %s mask\n",
                nx, ny, l2_flag_lname[band - 1]);
        band--;
        break;
    case FCLEAN:
        if (band < 1 || band > L1_NFLAGS) {
            if (want_verbose)
                printf("-E- %s %d: bogus band number %d for filter %s\n",
                    __FILE__, __LINE__, band, fname);
            return (0);
        }
        if (want_verbose)
            printf("Setting %d x %d cleaning filter on %s mask\n",
                nx, ny, l2_flag_lname[band - 1]);
        band--;
        break;
    case FLTMEAN:
        if (want_verbose)
            printf("Setting %d x %d averaging filter on Lt(%d)\n",
                nx, ny, band);
        if (band < 1 || band > nbands) {
            if (want_verbose)
                printf("-E- %s %d: bogus band number %d for filter %s\n",
                    __FILE__, __LINE__, band, fname);
            return (0);
        }
        band--;
        break;
    case FLTRMEAN:
        if (want_verbose)
            printf("Setting %d x %d averaging filter on Lt(%d)-Lr(%d)\n",
                nx, ny, band, band);
        if (band < 1 || band > nbands) {
            if (want_verbose)
                printf("-E- %s %d: bogus band number %d for filter %s\n",
                    __FILE__, __LINE__, band, fname);
            return (0);
        }
        band--;
        break;
    case FLTRREJECT:
        if (want_verbose)
            printf("Setting %d x %d rejection filter on Lt(%d)-Lr(%d)\n",
                nx, ny, band, band);
        if (band < 1 || band > nbands) {
            if (want_verbose)
                printf("-E- %s %d: bogus band number %d for filter %s\n",
                    __FILE__, __LINE__, band, fname);
            return (0);
        }
        band--;
        break;
    case FLTRIQMEAN:
        if (want_verbose)
            printf("Setting %d x %d interquartile averaging filter on Lt(%d)-Lr(%d)\n",
                nx, ny, band, band);
        if (band < 1 || band > nbands) {
            if (want_verbose)
                printf("-E- %s %d: bogus band number %d for filter %s\n",
                    __FILE__, __LINE__, band, fname);
            return (0);
        }
        band--;
        break;
    case FLTMED:
        if (want_verbose)
            printf("Setting %d x %d median filter on Lt(%d)\n",
                nx, ny, band);
        if (band < 1 || band > nbands) {
            if (want_verbose)
                printf("-E- %s %d: bogus band number %d for filter %s\n",
                    __FILE__, __LINE__, band, fname);
            return (0);
        }
        band--;
        break;
    case FLTRMED:
        if (want_verbose)
            printf("Setting %d x %d median filter on Lt(%d)-Lr(%d)\n",
                nx, ny, band, band);
        if (band < 1 || band > nbands) {
            if (want_verbose)
                printf("-E- %s %d: bogus band number %d for filter %s\n",
                    __FILE__, __LINE__, band, fname);
            return (0);
        }
        band--;
        break;
    case FEPSMEAN:
        if (want_verbose)
            printf("Setting %d x %d smoothing filter on epsilon(%d)\n",
                nx, ny, band);
        if (band < 1 || band + 1 > nbands) {
            if (want_verbose)
                printf("-E- %s %d: bogus band number %d for filter %s\n",
                    __FILE__, __LINE__, band, fname);
            return (0);
        }
        band--;
        break;
    case FBTDETAVG:
        if (want_verbose)
            printf("Setting %d x %d filter on IR band %d for detector %d\n",
                nx, ny, band, minfill);
        band--; /* band number    (0-based) */
        minfill--; /* detector number (0-based)*/
        break;
    case FTEST:
        if (want_verbose)
            printf("Setting %d x %d test filter on %d\n",
                nx, ny, band);
        band--;
        break;
    default:
        printf("-E- %s %d: unknown filter function %s\n",
                __FILE__, __LINE__, fname);
        return (0);
    }

    fctl->nfilt++;
    ;
    fctl->npix = MAX(nx, fctl->npix);
    fctl->nscan = MAX(ny, fctl->nscan);

    fctl->f[i].func = func;
    fctl->f[i].band = band;
    fctl->f[i].nx = nx;
    fctl->f[i].ny = ny;
    fctl->f[i].minfill = minfill;

    get_kernel(&fctl->f[i]);

    if (want_verbose) {
        printf("\nFilter Kernel");
        for (j = 0; j < fctl->f[i].nx * fctl->f[i].ny; j++) {
            if ((j % fctl->f[i].nx) == 0) printf("\n");
            printf("%d ", fctl->f[i].kernel[j]);
        }
        printf("\n\nMinimum fill set to %d pixels\n", fctl->f[i].minfill);
        printf("\n\n");
    }

    return (1);
}

int rdfilter(char *file, fctlstr *fctl, int32_t nbands) {
    FILE *fp;
    char line[80];
    char *p1, *p2;
    char fname[80];
    int band;
    int nscan;
    int32_t npix;
    int32_t minfill;
    int i;

    if (want_verbose)
        printf("Opening filter file %s\n", file);

    if ((fp = fopen(file, "r")) == NULL) {
        printf("The specified filter file (%s) was not found.\n", file);
        exit(1);
    }

    while (fgets(line, 80, fp)) {
        if (line[0] == '#' || line[0] == '\n')
            continue;

        p1 = line;
        if (!(p2 = strchr(p1, ','))) {
            printf("-E- %s %d: filter parsing error on %s. "
                    "Expecting comma-separated list.\n",
                    __FILE__, __LINE__, file);
            exit(1);
        }

        memset(fname, '\0', sizeof (fname));
        strncpy(fname, p1, p2 - p1);
        for (i = 0; i < (p2 - p1); i++)
            fname[i] = tolower(fname[i]);

        p1 = p2 + 1;
        if (!(p2 = strchr(p1, ',')))
            continue;
        band = atoi(p1);

        p1 = p2 + 1;
        if (!(p2 = strchr(p1, ',')))
            continue;
        npix = atoi(p1);

        p1 = p2 + 1;
        if (!(p2 = strchr(p1, ',')))
            continue;
        nscan = atoi(p1);

        p1 = p2 + 1;
        minfill = atoi(p1);

        fctl_set(fctl, npix, fname, band, npix, nscan, minfill, nbands);

    }

    return (0);
}

void fdilate(l1qstr *l1que, int32_t nx, int32_t ny, int flag, char kernel[], l1str *l1rec) {
    int32_t nscan = l1que->nq;
    int32_t npix = l1rec->npix;
    int32_t ip, ip1, ip2;
    int32_t is, is1, is2;
    int32_t i, j, k;
    char *pflag[NQMAX];
    char *pf;
    uint32_t flagID = pow(2L, flag);

    /* Get pointer to desired l1rec flag */
    switch (flagID) {
    case HILT: pf = l1rec->hilt;
        break;
    case LAND: pf = l1rec->land;
        break;
    case CLOUD: pf = l1rec->cloud;
        break;
    case HIGLINT: pf = l1rec->glint;
        break;
    default:
        fprintf(stderr,
                "-E- %s line %d: The specified flag, %d (%d), is not supported\n",
                __FILE__, __LINE__, flag + 1, flagID);
        fprintf(stderr, "for dilation filtering.\n");
        exit(FATAL_ERROR);
    }

    /* Build array of pointers to equivalent queue flag */
    for (is = 0; is < nscan; is++)
        switch (flagID) {
        case HILT: pflag[is] = l1que->r[is].hilt;
            break;
        case LAND: pflag[is] = l1que->r[is].land;
            break;
        case CLOUD: pflag[is] = l1que->r[is].cloud;
            break;
        case HIGLINT: pflag[is] = l1que->r[is].glint;
            break;
        }

    /* Compute queue scan limits for the ROI */
    is = nscan / 2;
    is1 = MIN(MAX(0, is - ny / 2), nscan - 1);
    is2 = MAX(MIN(nscan - 1, is + ny / 2), 0);

    /* Loop over each pixel of the central scan */

    for (ip = 0; ip < l1rec->npix; ip++) {

        /* compute pixel neighbor limits for the ROI */
        ip1 = MIN(MAX(0, ip - nx / 2), npix - 1);
        ip2 = MAX(MIN(npix - 1, ip + nx / 2), 0);

        /* If the pixel is already flagged, just mask it and go on. */
        /* Also mask the surrounding pixels in the queue records.   */
        /* This ensures that the dilation masking will be visible   */
        /* to later smoothing processes.                            */

        if (*(pf + ip)) {
            l1rec->mask[ip] = ON;
            for (j = is1; j <= is2; j++) for (i = ip1; i <= ip2; i++) {
                    k = (j - is1) * nx + (i - ip1);
                    if (kernel[k])
                        l1que->r[j].mask[i] = ON;
                }
            continue;
        }

        /* Search the ROI for the desired flag */
        for (j = is1; j <= is2; j++) for (i = ip1; i <= ip2; i++) {
                k = (j - is1) * nx + (i - ip1);
                if (kernel[k] && *(pflag[j] + i)) {
                    *(pf + ip) = ON; /* set flag in l1rec   */
                    l1rec->mask[ip] = ON; /* set mask in l1rec   */
                    goto next_pixel;
                }
            }
    next_pixel:
        ;
    }
}

void fstlight(l1qstr *l1que, int32_t nx, int32_t ny, int32_t dscan, int flag, char kernel[], l1str *l1rec)
/*
 *  W. Robinson, SAIC, 16 Aug 2011 - modify for VIIRS aggregation zones 
 */ {
    int32_t nscan = l1que->nq;
    int32_t npix = l1rec->npix;
    int32_t ip, ip1, ip2;
    int32_t is, is1, is2;
    int32_t i, j, k;
    int32_t ip_scan, filt_dist, ipmod;
    char *pflag[NQMAX];
    char *pf;
    uint32_t flagID = pow(2L, flag);

    /* Get pointer to desired l1rec flag */
    switch (flagID) {
    case HILT: pf = l1rec->hilt;
        break;
    case LAND: pf = l1rec->land;
        break;
    case CLOUD: pf = l1rec->cloud;
        break;
    case HIGLINT: pf = l1rec->glint;
        break;
    default:
        fprintf(stderr,
                "-E- %s line %d: The specified flag, %d (%d), is not supported\n",
                __FILE__, __LINE__, flag + 1, flagID);
        fprintf(stderr, "for dilation filtering.\n");
        exit(FATAL_ERROR);
    }

    /* Build array of pointers to equivalent queue flag */
    for (is = 0; is < nscan / dscan; is++)
        switch (flagID) {
        case HILT: pflag[is] = l1que->r[is].hilt;
            break;
        case LAND: pflag[is] = l1que->r[is].land;
            break;
        case CLOUD: pflag[is] = l1que->r[is].cloud;
            break;
        case HIGLINT: pflag[is] = l1que->r[is].glint;
            break;
        }

    /* Compute queue scan limits for the ROI */
    is = nscan / 2 / dscan;
    is1 = MIN(MAX(0, is - ny / dscan / 2), nscan / dscan - 1);
    is2 = MAX(MIN(nscan / dscan - 1, is + ny / dscan / 2), 0);
    /* Loop over each pixel of the central scan */

    for (ip = 0; ip < l1rec->npix; ip++) {

        /* compute pixel neighbor limits for the ROI */
        /* for aggregated VIIRS, limits are aggregation zone dependent */
        if (((l1rec->l1file->sensorID == VIIRSN) ||
	     (l1rec->l1file->sensorID == VIIRSJ1) ||
	     (l1rec->l1file->sensorID == VIIRSJ2)) && 
            (l1rec->scn_fmt == 0)) {
            ip_scan = l1rec->pixnum[ip]; /* scan pixel */
            int spix = ip_scan - ip;
            filt_dist = -nx / 2;
            viirs_pxcvt_agdel(ip_scan, filt_dist, &ipmod);
            ipmod -= spix;
            ip1 = MIN(MAX(0, ipmod), npix - 1);

            filt_dist = nx / 2;
            viirs_pxcvt_agdel(ip_scan, filt_dist, &ipmod);
            ipmod -= spix;
            ip2 = MAX(MIN(npix - 1, ipmod), 0);
        } else {
            ip1 = MIN(MAX(0, ip - nx / 2), npix - 1);
            ip2 = MAX(MIN(npix - 1, ip + nx / 2), 0);
        }

        /* If the pixel is already flagged, skip it. */
        if (*(pf + ip))
            continue;

        /* Search the ROI for the desired flag */
        for (j = is1; j <= is2; j++) for (i = ip1; i <= ip2; i++) {
                k = (j - is1) * nx + (i - ip1);
                if (kernel[k] && *(pflag[j] + i)) {
                    l1rec->stlight[ip] = ON; /* set straylight flag in l1rec   */
                    goto next_pixel;
                }
            }
    next_pixel:
        ;
    }
}

void fclean(l1qstr *l1que, int32_t nx, int32_t ny, int flag, char kernel[], l1str *l1rec) {
    int32_t nscan = l1que->nq;
    int32_t npix = l1rec->npix;
    int32_t ip, ip1, ip2;
    int32_t is, is1, is2;
    int32_t i, j, k;
    char *pflag[NQMAX];
    char *pf;
    uint32_t flagID = pow(2L, flag);
    int32_t cnt = 0;
    int32_t tot = 0;

    /* Get pointer to desired l1rec flag */
    switch (flagID) {
    case HILT: pf = l1rec->hilt;
        break;
    case LAND: pf = l1rec->land;
        break;
    case CLOUD: pf = l1rec->cloud;
        break;
    case HIGLINT: pf = l1rec->glint;
        break;
    default:
        fprintf(stderr,
                "-E- %s line %d: The specified flag, %d (%d), is not supported\n",
                __FILE__, __LINE__, flag + 1, flagID);
        fprintf(stderr, "for dilation filtering.\n");
        exit(FATAL_ERROR);
    }

    /* Build array of pointers to equivalent queue flag */
    for (is = 0; is < nscan; is++)
        switch (flagID) {
        case HILT: pflag[is] = l1que->r[is].hilt;
            break;
        case LAND: pflag[is] = l1que->r[is].land;
            break;
        case CLOUD: pflag[is] = l1que->r[is].cloud;
            break;
        case HIGLINT: pflag[is] = l1que->r[is].glint;
            break;
        }

    /* Compute queue scan limits for the ROI */
    is = nscan / 2;
    is1 = MIN(MAX(0, is - ny / 2), nscan - 1);
    is2 = MAX(MIN(nscan - 1, is + ny / 2), 0);

    /* Loop over each pixel of the central scan */

    for (ip = 0; ip < l1rec->npix; ip++) {

        /* compute pixel neighbor limits for the ROI */
        ip1 = MIN(MAX(0, ip - nx / 2), npix - 1);
        ip2 = MAX(MIN(npix - 1, ip + nx / 2), 0);

        /* If the pixel is already flagged, just mask it and go on. */
        if (*(pf + ip)) {
            l1rec->mask[ip] = ON;
            continue;
        }

        /* Search the ROI for the desired flag */
        for (j = is1; j <= is2; j++) for (i = ip1; i <= ip2; i++) {
                k = (j - is1) * nx + (i - ip1);
                if (kernel[k]) {
                    tot++;
                    if (*(pflag[j] + i))
                        cnt++;
                }
            }

        /* If every surrounding pixel is flagged, then flag the center pixel */
        if (tot - cnt == 1) {
            *(pf + ip) = ON; /* set flag in l1rec   */
            l1rec->mask[ip] = ON; /* set mask in l1rec   */
        }
    }
}

void fLTmean(l1qstr *l1que, int32_t nx, int32_t ny, int ib, int32_t minfill, char kernel[],
        l1str *l1rec) {
    int32_t nscan = l1que->nq;
    int32_t npix = l1rec->npix;
    int32_t ip, ip1, ip2;
    int32_t is, is1, is2;
    int32_t ii, iib;
    int32_t i, j, k;
    float x;
    int32_t cnt;

    /* define range of scan numbers within queue */
    is = nscan / 2;
    is1 = MIN(MAX(0, is - ny / 2), nscan - 1);
    is2 = MAX(MIN(nscan - 1, is + ny / 2), 0);

    /*                                           */
    /* Loop through each pixel in the scan       */
    /*                                           */
    for (ip = 0; ip < l1rec->npix; ip++) {

        iib = ip * l1rec->l1file->nbands + ib;

        /* if the pixel is already masked, skip it */
        if (l1rec->mask[ip] || l1rec->Lt[iib] == badval)
            continue;

        /* define pixel range around the central pixel */
        ip1 = MIN(MAX(0, ip - nx / 2), npix - 1);
        ip2 = MAX(MIN(npix - 1, ip + nx / 2), 0);


        /*                                                */
        /* Compute the mean over the pixel and scan range */
        /* and replace the central pixel value with it.   */
        /*                                                */
        x = 0.0;
        cnt = 0;
        for (i = ip1; i <= ip2; i++) {
            for (j = is1; j <= is2; j++) {
                ii = i * l1rec->l1file->nbands + ib;
                k = (j - is1) * nx + (i - ip1);
                if (kernel[k] && !l1que->r[j].mask[i] && l1que->r[j].Lt[ii] != badval) {
                    x += l1que->r[j].Lt[ii];
                    cnt++;
                }
            }
        }
        if (cnt >= minfill)
            l1rec->Lt[iib] = x / cnt;
        else {
            l1rec->filter[ip] = ON;
            l1rec->mask[ip] = ON;
        }
    }
}

void fLTRmean(l1qstr *l1que, int32_t nx, int32_t ny, int ib, int32_t minfill, char kernel[],
        l1str *l1rec) {
    int32_t nscan = l1que->nq;
    int32_t npix = l1rec->npix;
    int32_t ip, ip1, ip2;
    int32_t is, is1, is2;
    int32_t ii, iib;
    int32_t i, j, k;
    float x;
    int32_t cnt;

    /* define range of scan numbers within queue */
    is = nscan / 2;
    is1 = MIN(MAX(0, is - ny / 2), nscan - 1);
    is2 = MAX(MIN(nscan - 1, is + ny / 2), 0);

    /*                                           */
    /* Loop through each pixel in the scan       */
    /*                                           */
    for (ip = 0; ip < l1rec->npix; ip++) {

        iib = ip * l1rec->l1file->nbands + ib;

        /* if the pixel is already masked, skip it */
        if (l1rec->mask[ip] || l1rec->Lt[iib] == badval)
            continue;

        /* define pixel range around the central pixel */
        ip1 = MIN(MAX(0, ip - nx / 2), npix - 1);
        ip2 = MAX(MIN(npix - 1, ip + nx / 2), 0);


        /*                                                */
        /* Compute the mean over the pixel and scan range */
        /* and replace the central pixel value with it.   */
        /*                                                */
        x = 0.0;
        cnt = 0;
        for (i = ip1; i <= ip2; i++) {
            for (j = is1; j <= is2; j++) {
                ii = i * l1rec->l1file->nbands + ib;
                k = (j - is1) * nx + (i - ip1);
                if (kernel[k] && !l1que->r[j].mask[i] && l1que->r[j].Lt[ii] != badval) {
                    x += (l1que->r[j].Lt[ii] - l1que->r[j].Lr[ii]);
                    cnt++;
                }
            }
        }
        if (cnt >= minfill)
            l1rec->Lt[iib] = x / cnt + l1rec->Lr[iib];
        else {
            l1rec->filter[ip] = ON;
            l1rec->mask[ip] = ON;
        }
    }
}

int compfloat(float *x, float *y) {
    if (*x < *y)
        return (-1);
    else
        return ( 1);
}

int compfnode(fnode *x, fnode *y) {
    if (x->v < y->v)
        return (-1);
    else
        return ( 1);
}

void fLTRreject(l1qstr *l1que, int32_t nx, int32_t ny, int ib, int32_t minfill, char kernel[],
        l1str *l1rec) {
    static fnode *medx = NULL;
    static fnode *mad = NULL;
    static int32_t len = 0;

    int32_t nscan = l1que->nq;
    int32_t npix = l1rec->npix;
    int32_t ip, ip1, ip2;
    int32_t is, is1, is2;
    int32_t ii, iib;
    int32_t i, j, k;
    float x;
    int32_t cnt;
    float median;
    float madev;

    /* define range of scan numbers within queue */
    is = nscan / 2;
    is1 = MIN(MAX(0, is - ny / 2), nscan - 1);
    is2 = MAX(MIN(nscan - 1, is + ny / 2), 0);

    /* allocate sufficient workspace for the filter size */
    if (medx == NULL || len < nx * ny) {
        len = nx*ny;
        if (medx != NULL) free(medx);
        if ((medx = (fnode *) malloc(len * sizeof (fnode))) == NULL) {
            printf("-E- %s %d: Unable to allocate workspace for median\n",
                    __FILE__, __LINE__);
            return;
        }
    }
    if (mad == NULL || len < nx * ny) {
        len = nx*ny;
        if (mad != NULL) free(mad);
        if ((mad = (fnode *) malloc(len * sizeof (fnode))) == NULL) {
            printf("-E- %s %d: Unable to allocate workspace for median absolute deviation\n",
                    __FILE__, __LINE__);
            return;
        }
    }

    /*                                           */
    /* Loop through each pixel in the scan       */
    /*                                           */
    for (ip = 0; ip < l1rec->npix; ip++) {

        iib = ip * l1rec->l1file->nbands + ib;

        /* if the pixel is already masked, skip it */
        if (l1rec->mask[ip] || l1rec->Lt[iib] == badval)
            continue;

        /* define pixel range around the central pixel */
        ip1 = MIN(MAX(0, ip - nx / 2), npix - 1);
        ip2 = MAX(MIN(npix - 1, ip + nx / 2), 0);


        /*                                                */
        /* Compute the mean over the pixel and scan range */
        /* and replace the central pixel value with it.   */
        /*                                                */
        cnt = 0;
        for (i = ip1; i <= ip2; i++) {
            for (j = is1; j <= is2; j++) {
                ii = i * l1rec->l1file->nbands + ib;
                k = (j - is1) * nx + (i - ip1);
                if (kernel[k] && !l1que->r[j].mask[i] && l1que->r[j].Lt[ii] != badval) {
                    medx[cnt].v = l1que->r[j].Lt[ii] - l1que->r[j].Lr[ii];
                    medx[cnt].i = i;
                    medx[cnt].j = j;
                    cnt++;
                }
            }
        }
        if (cnt >= minfill && cnt >= 2) {
            qsort(medx, cnt, sizeof (fnode),
                    (int (*)(const void *, const void *)) compfnode);

            median = medx[cnt / 2].v;

            cnt = 0;
            for (i = ip1; i <= ip2; i++) {
                for (j = is1; j <= is2; j++) {
                    ii = i * l1rec->l1file->nbands + ib;
                    k = (j - is1) * nx + (i - ip1);
                    if (kernel[k] && !l1que->r[j].mask[i] && l1que->r[j].Lt[ii] != badval) {
                        x = (l1que->r[j].Lt[ii] - l1que->r[j].Lr[ii]) - median;
                        mad[cnt].v = fabs(x);
                        mad[cnt].i = i;
                        mad[cnt].j = j;
                        cnt++;
                    }
                }
            }
            qsort(mad, cnt, sizeof (fnode),
                    (int (*)(const void *, const void *)) compfnode);

            madev = 1.4826 * mad[cnt / 2].v; // const. 1.4826 makes MAD a consistent
            // esitmator - equiv. to stdev

            if (fabs((l1rec->Lt[iib] - l1rec->Lr[iib] - median) / madev) >= 3) {
                l1rec->filter[ip] = ON;
                //l1rec->mask[ip]   = ON;
            }
        }
    }
}

void fEPSmean(l1qstr *l1que, int32_t nx, int32_t ny, int ib1, int32_t minfill, char kernel[],
        l1str *l1rec) {
    int32_t nscan = l1que->nq;
    int32_t npix = l1rec->npix;
    int32_t ip, ip1, ip2;
    int32_t is, is1, is2;
    int32_t ii, iib;
    int32_t i, j, k;
    float x1;
    float x2;
    int32_t cnt;
    float r_bar;
    float La1;
    float La2;
    float La;
    int ib2 = l1rec->l1file->nbands - 1;

    /* define range of scan numbers within queue */
    is = nscan / 2;
    is1 = MIN(MAX(0, is - ny / 2), nscan - 1);
    is2 = MAX(MIN(nscan - 1, is + ny / 2), 0);

    /*                                           */
    /* Loop through each pixel in the scan       */
    /*                                           */
    for (ip = 0; ip < l1rec->npix; ip++) {

        iib = ip * l1rec->l1file->nbands;

        /* if the pixel is already masked, skip it */
        if (l1rec->mask[ip] ||
                l1rec->Lt[iib + ib1] == badval ||
                l1rec->Lt[iib + ib2] == badval)
            continue;

        /* define pixel range around the central pixel */
        ip1 = MIN(MAX(0, ip - nx / 2), npix - 1);
        ip2 = MAX(MIN(npix - 1, ip + nx / 2), 0);

        /*                                                */
        /* Compute the mean over the pixel and scan range */
        /* and replace the central pixel value with it.   */
        /*                                                */
        x1 = 0.0;
        x2 = 0.0;
        cnt = 0;
        La = 0.0;
        for (i = ip1; i <= ip2; i++) {
            for (j = is1; j <= is2; j++) {
                k = (j - is1) * nx + (i - ip1);
                if (kernel[k] && !l1que->r[j].mask[i]) {
                    La1 = 0.0;
                    ii = i * l1rec->l1file->nbands + ib1;
                    if (l1que->r[j].Lt[ii] != badval) {
                        La1 = (l1que->r[j].Lt[ii]
                                / l1que->r[j].tg_sol[ii]
                                / l1que->r[j].tg_sen[ii]
                                - l1que->r[j].tLf[ii]
                                - l1que->r[j].Lr[ii])
                                / l1que->r[j].t_o2[ii];
                    }
                    La2 = 0.0;
                    ii = i * l1rec->l1file->nbands + ib2;
                    if (l1que->r[j].Lt[ii] != badval) {
                        La2 = (l1que->r[j].Lt[ii]
                                / l1que->r[j].tg_sol[ii]
                                / l1que->r[j].tg_sen[ii]
                                - l1que->r[j].tLf[ii]
                                - l1que->r[j].Lr[ii])
                                / l1que->r[j].t_o2[ii];
                    }

                    /* Note: we may be averaging negative values (dark pixels) */
                    x1 += La1;
                    x2 += La2;
                    cnt++;

                    /* Center-pixel long-wavelength value */
                    if (i == ip && j == is) {
                        La = La2;
                    }
                }
            }
        }

        if (cnt >= minfill) {
            ii = ip * l1rec->l1file->nbands + ib1;
            /* if region is anomalously dark, we'll just leave it alone */
            /* here and handle it elsewhere                             */
            if (x1 > 0 && x2 > 0) {
                r_bar = x1 / x2; /* local average band ratio */
                La1 = La*r_bar; /* modify La1 at center_pixel */
                l1rec->Lt[ii] = (La1 * l1rec->t_o2[ii] + l1rec->Lr[ii] + l1rec->tLf[ii])
                        * l1rec->tg_sol[ii]
                        * l1rec->tg_sen[ii];
            }
        } else {
            l1rec->filter[ip] = ON;
            l1rec->mask[ip] = ON;
        }
    }
}

void fBTdetavg(l1qstr *l1que, int32_t nx, int32_t ny, int ib, int32_t id, char kernel[],
        l1str *l1rec) {
    int32_t ndet = 10; /* need to generalize this later */
    int32_t nscan = l1que->nq;
    int32_t npix = l1rec->npix;
    int32_t detnum = l1rec->detnum;
    int32_t scanum = l1rec->iscan;
    int32_t ip, ip1, ip2;
    int32_t is, is1, is2;
    int32_t ii, ipb;
    int32_t i, j, k;
    float x;
    int32_t cnt;

    /* only apply to a specific detector number */
    if (detnum != id) return;

    /* define range of scan numbers within queue */
    is = nscan / 2;
    is1 = MIN(MAX(0, is - ny / 2), nscan - 1);
    is2 = MAX(MIN(nscan - 1, is + ny / 2), 0);

    /* limit to detector boundaries */
    if (detnum == 0) is1 = is + 1;
    if (detnum == ndet - 1) is2 = is - 1;

    printf("%d %d %d %d\n", ib, id, is1, is2);

    /*                                           */
    /* Loop through each pixel in the scan       */
    /*                                           */
    for (ip = 0; ip < l1rec->npix; ip++) {

        ipb = ip * NBANDSIR + ib;

        /* define pixel range around the central pixel */
        ip1 = MIN(MAX(0, ip - nx / 2), npix - 1);
        ip2 = MAX(MIN(npix - 1, ip + nx / 2), 0);

        /*                                                */
        /* Compute the mean over the pixel and scan range */
        /* and replace the central pixel value with it.   */
        /*                                                */
        x = 0.0;
        cnt = 0;
        for (i = ip1; i <= ip2; i++) {
            for (j = is1; j <= is2; j++) {
                if (j == scanum) continue;
                ii = i * NBANDSIR + ib;
                k = (j - is1) * nx + (i - ip1);
                if (kernel[k] && !l1que->r[j].mask[i] && l1que->r[j].Ltir[ii] != badval) {
                    x += l1que->r[j].Bt[ii];
                    cnt++;
                }
            }
        }
        if (cnt > 0) {
            l1rec->Bt[ipb] = x / cnt;
            l1rec->filter[ip] = ON;
        } else {
            l1rec->filter[ip] = ON;
            l1rec->mask[ip] = ON;
        }
    }
}

void fLTmed(l1qstr *l1que, int32_t nx, int32_t ny, int ib, int32_t minfill, char kernel[],
        l1str *l1rec) {
    static float *x = NULL;
    static int32_t len = 0;

    int32_t nscan = l1que->nq;
    int32_t npix = l1rec->npix;
    int32_t ip, ip1, ip2;
    int32_t is, is1, is2;
    int32_t ii, iib;
    int32_t i, j, k;
    int32_t cnt;

    /* allocate sufficient workspace for the filter size */
    if (x == NULL || len < nx * ny) {
        len = nx*ny;
        if (x != NULL) free(x);
        if ((x = (float *) malloc(len * sizeof (float))) == NULL) {
            printf("-E- %s %d: Unable to allocate workspace for median\n",
                    __FILE__, __LINE__);
            return;
        }
    }

    /* define range of scan numbers within queue */
    is = nscan / 2;
    is1 = MIN(MAX(0, is - ny / 2), nscan - 1);
    is2 = MAX(MIN(nscan - 1, is + ny / 2), 0);

    /*                                           */
    /* Loop through each pixel in the scan       */
    /*                                           */
    for (ip = 0; ip < l1rec->npix; ip++) {

        iib = ip * l1rec->l1file->nbands + ib;

        /* if the pixel is already masked, skip it */
        if (l1rec->mask[ip] || l1rec->Lt[iib] == badval)
            continue;

        /* define pixel range around the central pixel */
        ip1 = MIN(MAX(0, ip - nx / 2), npix - 1);
        ip2 = MAX(MIN(npix - 1, ip + nx / 2), 0);


        /*                                                */
        /* Compute the mean over the pixel and scan range */
        /* and replace the central pixel value with it.   */
        /*                                                */
        cnt = 0;
        for (i = ip1; i <= ip2; i++) {
            for (j = is1; j <= is2; j++) {
                ii = i * l1rec->l1file->nbands + ib;
                k = (j - is1) * nx + (i - ip1);
                if (kernel[k] && !l1que->r[j].mask[i] && l1que->r[j].Lt[ii] != badval) {
                    x[cnt] = l1que->r[j].Lt[ii];
                    cnt++;
                }
            }
        }
        if (cnt >= minfill) {
            qsort(x, cnt, sizeof (float),
                    (int (*)(const void *, const void *)) compfloat);
            l1rec->Lt[iib] = x[cnt / 2];
        } else {
            l1rec->filter[ip] = ON;
            l1rec->mask[ip] = ON;
        }
    }
}

void fLTRmed(l1qstr *l1que, int32_t nx, int32_t ny, int ib, int32_t minfill, char kernel[],
        l1str *l1rec) {
    static fnode *x = NULL;
    static int32_t len = 0;

    int32_t nscan = l1que->nq;
    int32_t npix = l1rec->npix;
    int32_t ip, ip1, ip2;
    int32_t is, is1, is2;
    int32_t ii, iib;
    int32_t i, j, k;
    int32_t cnt;

    /* allocate sufficient workspace for the filter size */
    if (x == NULL || len < nx * ny) {
        len = nx*ny;
        if (x != NULL) free(x);
        if ((x = (fnode *) malloc(len * sizeof (fnode))) == NULL) {
            printf("-E- %s %d: Unable to allocate workspace for median\n",
                    __FILE__, __LINE__);
            return;
        }
    }

    /* define range of scan numbers within queue */
    is = nscan / 2;
    is1 = MIN(MAX(0, is - ny / 2), nscan - 1);
    is2 = MAX(MIN(nscan - 1, is + ny / 2), 0);

    /*                                           */
    /* Loop through each pixel in the scan       */
    /*                                           */
    for (ip = 0; ip < l1rec->npix; ip++) {

        iib = ip * l1rec->l1file->nbands + ib;

        /* if the pixel is already masked, skip it */
        if (l1rec->mask[ip] || l1rec->Lt[iib] == badval)
            continue;

        /* define pixel range around the central pixel */
        ip1 = MIN(MAX(0, ip - nx / 2), npix - 1);
        ip2 = MAX(MIN(npix - 1, ip + nx / 2), 0);

        /*                                                      */
        /* Accumulate non-masked values over the pixel and scan */
        /*                                                      */
        cnt = 0;
        for (i = ip1; i <= ip2; i++) {
            for (j = is1; j <= is2; j++) {
                ii = i * l1rec->l1file->nbands + ib;
                k = (j - is1) * nx + (i - ip1);
                if (kernel[k] && !l1que->r[j].mask[i] && l1que->r[j].Lt[ii] != badval) {
                    x[cnt].v = l1que->r[j].Lt[ii] - l1que->r[j].Lr[ii];
                    x[cnt].i = i;
                    x[cnt].j = j;
                    cnt++;
                }
            }
        }
        if (cnt >= minfill) {
            qsort(x, cnt, sizeof (fnode),
                    (int (*)(const void *, const void *)) compfnode);

            l1rec->Lt[iib] = x[cnt / 2].v + l1rec->Lr[iib];

        } else {
            l1rec->filter[ip] = ON;
            l1rec->mask[ip] = ON;
        }
    }
}

void fLTRiqmean(l1qstr *l1que, int32_t nx, int32_t ny, int ib, int32_t minfill, char kernel[],
        l1str *l1rec) {
    static fnode *x = NULL;
    static int32_t len = 0;

    int32_t nscan = l1que->nq;
    int32_t npix = l1rec->npix;
    int32_t ip, ip1, ip2;
    int32_t is, is1, is2;
    int32_t ii, iib;
    int32_t i, j, k;
    int32_t cnt, num;
    float val;

    /* allocate sufficient workspace for the filter size */
    if (x == NULL || len < nx * ny) {
        len = nx*ny;
        if (x != NULL) free(x);
        if ((x = (fnode *) malloc(len * sizeof (fnode))) == NULL) {
            printf("-E- %s %d: Unable to allocate workspace for median\n",
                    __FILE__, __LINE__);
            return;
        }
    }

    /* define range of scan numbers within queue */
    is = nscan / 2;
    is1 = MIN(MAX(0, is - ny / 2), nscan - 1);
    is2 = MAX(MIN(nscan - 1, is + ny / 2), 0);

    /*                                           */
    /* Loop through each pixel in the scan       */
    /*                                           */
    for (ip = 0; ip < l1rec->npix; ip++) {

        iib = ip * l1rec->l1file->nbands + ib;

        /* if the pixel is already masked, skip it */
        if (l1rec->mask[ip] || l1rec->Lt[iib] == badval)
            continue;

        /* define pixel range around the central pixel */
        ip1 = MIN(MAX(0, ip - nx / 2), npix - 1);
        ip2 = MAX(MIN(npix - 1, ip + nx / 2), 0);

        /*                                                      */
        /* Accumulate non-masked values over the pixel and scan */
        /*                                                      */
        cnt = 0;
        for (i = ip1; i <= ip2; i++) {
            for (j = is1; j <= is2; j++) {
                ii = i * l1rec->l1file->nbands + ib;
                k = (j - is1) * nx + (i - ip1);
                if (kernel[k] && !l1que->r[j].mask[i] && l1que->r[j].Lt[ii] != badval) {
                    x[cnt].v = l1que->r[j].Lt[ii] - l1que->r[j].Lr[ii];
                    x[cnt].i = i;
                    x[cnt].j = j;
                    cnt++;
                }
            }
        }
        if (cnt >= minfill) {
            qsort(x, cnt, sizeof (fnode),
                    (int (*)(const void *, const void *)) compfnode);

            num = 0;
            val = 0.0;

            for (i = cnt / 4; i < 3 * cnt / 4; i++) {
                num++;
                val += x[i].v;
            }
            val /= num;

            l1rec->Lt[iib] = val + l1rec->Lr[iib];

        } else {
            l1rec->filter[ip] = ON;
            l1rec->mask[ip] = ON;
        }
    }
}

void fEPSiqmean(l1qstr *l1que, int32_t nx, int32_t ny, int ib1, int32_t minfill, char kernel[],
        l1str *l1rec) {
    static fnode *x1 = NULL;
    static fnode *x2 = NULL;
    static int32_t len = 0;

    int32_t nscan = l1que->nq;
    int32_t npix = l1rec->npix;
    int32_t ip, ip1, ip2;
    int32_t is, is1, is2;
    int32_t ii, iib;
    int32_t i, j, k;
    int32_t cnt, num;
    float r_bar;
    float La1;
    float La2;
    float La;
    int ib2 = l1rec->l1file->nbands - 1;

    /* allocate sufficient workspace for the filter size */
    if (x1 == NULL || len < nx * ny) {
        len = nx*ny;
        if (x1 != NULL) free(x1);
        if ((x1 = (fnode *) malloc(len * sizeof (fnode))) == NULL) {
            printf("-E- %s %d: Unable to allocate workspace for median\n",
                    __FILE__, __LINE__);
            return;
        }
    }
    if (x2 == NULL || len < nx * ny) {
        len = nx*ny;
        if (x2 != NULL) free(x2);
        if ((x2 = (fnode *) malloc(len * sizeof (fnode))) == NULL) {
            printf("-E- %s %d: Unable to allocate workspace for median\n",
                    __FILE__, __LINE__);
            return;
        }
    }

    /* define range of scan numbers within queue */
    is = nscan / 2;
    is1 = MIN(MAX(0, is - ny / 2), nscan - 1);
    is2 = MAX(MIN(nscan - 1, is + ny / 2), 0);

    /*                                           */
    /* Loop through each pixel in the scan       */
    /*                                           */
    for (ip = 0; ip < l1rec->npix; ip++) {

        iib = ip * l1rec->l1file->nbands;

        /* if the pixel is already masked, skip it */
        if (l1rec->mask[ip] ||
                l1rec->Lt[iib + ib1] <= 0.0 ||
                l1rec->Lt[iib + ib2] <= 0.0)
            continue;

        /* define pixel range around the central pixel */
        ip1 = MIN(MAX(0, ip - nx / 2), npix - 1);
        ip2 = MAX(MIN(npix - 1, ip + nx / 2), 0);


        /*                                                      */
        /* Accumulate non-masked values over the pixel and scan */
        /*                                                      */
        cnt = 0;
        La = 0.0;
        for (i = ip1; i <= ip2; i++) {
            for (j = is1; j <= is2; j++) {

                k = (j - is1) * nx + (i - ip1);
                if (kernel[k] &&
                        !l1que->r[j].mask[i] &&
                        l1que->r[j].Lt[i * l1rec->l1file->nbands + ib1] > 0.0 &&
                        l1que->r[j].Lt[i * l1rec->l1file->nbands + ib2] > 0.0) {

                    La1 = 0.0;
                    ii = i * l1rec->l1file->nbands + ib1;
                    La1 = (l1que->r[j].Lt[ii]
                            / l1que->r[j].tg_sol[ii]
                            / l1que->r[j].tg_sen[ii]
                            - l1que->r[j].tLf[ii]
                            - l1que->r[j].Lr[ii])
                            / l1que->r[j].t_o2[ii];

                    La2 = 0.0;
                    ii = i * l1rec->l1file->nbands + ib2;
                    La2 = (l1que->r[j].Lt[ii]
                            / l1que->r[j].tg_sol[ii]
                            / l1que->r[j].tg_sen[ii]
                            - l1que->r[j].tLf[ii]
                            - l1que->r[j].Lr[ii])
                            / l1que->r[j].t_o2[ii];

                    /* Center-pixel long-wavelength value */
                    if (i == ip && j == is)
                        La = La2;

                    x1[cnt].v = La1;
                    x1[cnt].i = i;
                    x1[cnt].j = j;

                    x2[cnt].v = La2;
                    x2[cnt].i = i;
                    x2[cnt].j = j;
                    cnt++;

                }
            }
        }

        if (cnt >= minfill) {

            qsort(x1, cnt, sizeof (fnode),
                    (int (*)(const void *, const void *)) compfnode);
            qsort(x2, cnt, sizeof (fnode),
                    (int (*)(const void *, const void *)) compfnode);

            num = 0;
            La1 = 0.0;
            La2 = 0.0;

            for (i = cnt / 4; i <= 3 * cnt / 4; i++) {
                num++;
                La1 += x1[i].v;
                La2 += x2[i].v;
            }
            La1 /= num;
            La2 /= num;

            /* if region is anomalously dark, we'll just leave it alone */
            /* here and handle it elsewhere                             */

            if (La1 > 0.0 && La2 > 0.0 && La1 / La2 < 2.0) {
                ii = ip * l1rec->l1file->nbands + ib1;
                r_bar = La1 / La2; /* local average band ratio */
                La1 = La*r_bar; /* modify La1 at center_pixel */
                l1rec->Lt[ii] = (La1 * l1rec->t_o2[ii] + l1rec->Lr[ii] + l1rec->tLf[ii])
                        * l1rec->tg_sol[ii]
                        * l1rec->tg_sen[ii];
            }

        } else {
            l1rec->filter[ip] = ON;
            l1rec->mask[ip] = ON;
        }
    }
}

void filter(fctlstr *fctl, l1qstr *l1que, l1str *l1rec, int32_t dscan) {
    int i;
    int32_t func, band, nx, ny, minfill;
    char *kernel;

    for (i = 0; i < fctl->nfilt; i++) {

        func = fctl->f[i].func;
        band = fctl->f[i].band;
        nx = fctl->f[i].nx;
        ny = fctl->f[i].ny;
        minfill = fctl->f[i].minfill;
        kernel = fctl->f[i].kernel;

        switch (func) {
        case FDILATE: fdilate(l1que, nx, ny, band, kernel, l1rec);
            break;
        case FSTLIGHT: fstlight(l1que, nx, ny, dscan, band, kernel, l1rec);
            break;
        case FCLEAN: fclean(l1que, nx, ny, band, kernel, l1rec);
            break;
        case FLTMEAN: fLTmean(l1que, nx, ny, band, minfill, kernel, l1rec);
            break;
        case FBTDETAVG: fBTdetavg(l1que, nx, ny, band, minfill, kernel, l1rec);
            break;
        case FLTRMEAN: fLTRmean(l1que, nx, ny, band, minfill, kernel, l1rec);
            break;
        case FLTRIQMEAN: fLTRiqmean(l1que, nx, ny, band, minfill, kernel, l1rec);
            break;
        case FLTMED: fLTmed(l1que, nx, ny, band, minfill, kernel, l1rec);
            break;
        case FLTRMED: fLTRmed(l1que, nx, ny, band, minfill, kernel, l1rec);
            break;
        case FTEST: fLTRiqmean(l1que, nx, ny, band, minfill, kernel, l1rec);
            break;
        case FEPSMEAN: fEPSiqmean(l1que, nx, ny, band, minfill, kernel, l1rec);
            break;
        case FLTRREJECT: fLTRreject(l1que, nx, ny, band, minfill, kernel, l1rec);
            break;
        default:
            printf("-E- %s %d: unknown filter function %d\n",
                    __FILE__, __LINE__, func);
            return;
        }

    }

    return;
}


