/* =========================================================== */
/* Module aer_io.c                                             */
/*                                                             */
/* Functions to open and read an aerosol specification file.   */
/*                                                             */
/* Written By:                                                 */
/*                                                             */
/*     B. A. Franz                                             */
/*     SAIC General Sciences Corp.                             */
/*     NASA/SIMBIOS Project                                    */
/*     March 2001                                              */
/*                                                             */
/* =========================================================== */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "filehdr_struc.h"
#include "filehandle.h"
#include "aer_struc.h"

static FILE *fp = NULL;
static char current_file[FILENAME_MAX];


/* ----------------------------------------------------------- */
/* close_aer() - close aerosol file                            */

/* ----------------------------------------------------------- */
void close_aer(void) {
    fclose(fp);
}


/* ----------------------------------------------------------- */
/* read_aerhdr() - loads header info                           */

/* ----------------------------------------------------------- */
int read_aerhdr(filehandle *file) {
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
/* open_aer() - opens file if not already opened               */

/* ----------------------------------------------------------- */
int open_aer(filehandle *file) {
    if (fp == NULL || strcmp(file->name, current_file) != 0) {
        if (fp != NULL) close_aer();
        if ((fp = fopen(file->name, "r")) == NULL) {
            printf("-E- %s: Error opening %s for reading.\n",
                    __FILE__, file->name);
            return (1);
        }
        strcpy(current_file, file->name);
        if (read_aerhdr(file) != 0) {
            printf("-E- %s: Error reading header for %s.\n",
                    __FILE__, file->name);
            return (1);
        }
    }

    return (0);
}


/* ----------------------------------------------------------- */
/* read_aer() - reads one aerosol record                       */
/*                                                             */
/* B. A. Franz, GSC, SIMBIOS Project, March 2001               */

/* ----------------------------------------------------------- */
int read_aer(filehandle *file, int32_t recnum, aestr *aerec) {
    /*                                                         */
    /* Open the input file if it is not already open           */
    /*                                                         */
    if (open_aer(file) != 0) {
        printf("-E- %s: File open error.\n", __FILE__);
        return (1);
    }

    if (feof(fp)) {
        printf("-I- %s: End of aer file %s reached.",
                __FILE__, file->name);
        return (1);
    }

    if (fseek(fp, (recnum + 1) * file->length, SEEK_SET) != 0) {
        printf("-E- %s: Error seeking record %d in %s.",
                __FILE__, recnum, file->name);
        return (1);
    }

    if (fread(aerec->data, 1, file->length, fp) != file->length) {
        return (1);
    }

    aerec->length = file->length;
    aerec->npix = file->npix;

    return (0);
}





