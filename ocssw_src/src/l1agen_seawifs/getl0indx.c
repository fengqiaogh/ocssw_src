/* ============================================================== */
/* Module:  getl0indx.c                                           */
/* Purpose: generate content and quality index for L0 file        */
/* Author:  B.A. Franz, General Scences Corp., 9/97               */
/* ============================================================== */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <time.h>

#include "swl0_proto.h"

INT32 getl0indx(char *filename, swl0ctl *l0ctl, swl0indx *indx) {
    swl0hdr hdr;
    BYTE mnf[L0LEN];
    INT16 scid[2];
    INT16 ttag[4];
    INT16 mnftype;
    INT16 mnfnum;
    FLOAT64 startSec;
    FLOAT64 stopSec;
    FLOAT64 utime;
    FLOAT64 lastTime;
    INT32 irec;
    INT32 srec;
    INT32 erec;
    FILE *fp;
    INT32 nrecs;
    INT32 numbits;
    INT32 numerrs;
    INT32 totnumbits = 0;
    INT32 totnumerrs = 0;
    INT16 timeValid;

    /*									*/
    /* Open file for reading						*/
    /*									*/
    if ((fp = fopen(filename, "r")) == NULL) {
        fprintf(stderr,
                "-E- %s line %d: unable to open %s for reading\n",
                __FILE__, __LINE__, filename);
        return (FATAL_ERROR);
    }

    /*									*/
    /* Read Level-0 header                                       	*/
    /*									*/
    if (fread(&hdr, sizeof (hdr), 1, fp) != 1) {
        fprintf(stderr,
                "-E- %s line %d: error reading %s\n",
                __FILE__, __LINE__, filename);
        return (FATAL_ERROR);
    }
    if (endianess() == 1) {
        swapc_bytes((char *) &hdr.numbits, 4, 1);
        swapc_bytes((char *) &hdr.errbits, 4, 1);
        swapc_bytes((char *) &hdr.startTime, 4, 1);
        swapc_bytes((char *) &hdr.stopTime, 4, 1);
        swapc_bytes((char *) &hdr.numlines, 4, 1);
    }

    /* Check header magic number */
    if (strncmp((char *) hdr.id, "CWIF", 4) != 0) {
        fprintf(stderr,
                "-E- %s line %d: error in header format of %s\n",
                __FILE__, __LINE__, filename);
        fprintf(stderr, "  Header should start with CWIF string.\n");
        fprintf(stderr, "  This may not be level-0 format.\n");
    }

    /* Invoke filetype over-ride switch, if set */
    if (l0ctl->fileType >= 0)
        hdr.type = l0ctl->fileType;

    /* Header AOS and LOS should be greater than 09/01/1997 */
    if (hdr.stopTime <= hdr.startTime || hdr.startTime < 873072000) {
        if (hdr.type != HRPT) {
            fprintf(stderr,
                    "-E- %s line %d: error in header format of %s\n",
                    __FILE__, __LINE__, filename);
            fprintf(stderr, "  AOS and/or LOS times make no sense.\n");
        }
        hdr.startTime = 873072000;
        hdr.stopTime = (INT32) time(NULL);
    }

    /* Check validity of HRPT/GAC flag in header                        */
    if (hdr.type != HRPT && hdr.type != GAC) {
        fprintf(stderr,
                "-E- %s line %d: error in header format of %s\n",
                __FILE__, __LINE__, filename);
        fprintf(stderr, "  HRPT/Stored-GAC flag invalid, assuming HRPT.\n");
        hdr.type = HRPT;
    }

    /*									*/
    /* Determine total number of frames and allocate frame quality index*/
    /*									*/
    nrecs = (filesize(filename) - sizeof (hdr)) / sizeof (mnf);
    if (nrecs < 1) {
        fprintf(stderr,
                "-E- %s line %d: input file %s is empty\n",
                __FILE__, __LINE__, filename);
        return (FATAL_ERROR);
    }
    if ((indx->rec = (frameqc *) calloc(nrecs, sizeof (frameqc))) == NULL) {
        fprintf(stderr,
                "-E- %s line %d: error allocating memory for index array\n",
                __FILE__, __LINE__);
        return (FATAL_ERROR);
    }

    /*									*/
    /* Load or initialize level-0 indx structure			*/
    /*									*/
    strncpy(indx->l0file, filename, FILENAME_MAX);
    indx->timeAOS = hdr.startTime;
    indx->timeLOS = hdr.stopTime;
    indx->type = hdr.type;
    indx->nrecs = nrecs;
    indx->ngac = 0;
    indx->nlac = 0;
    indx->nlun = 0;
    indx->nsol = 0;
    indx->nigc = 0;
    indx->ntdi = 0;
    indx->nund = 0;
    indx->timeError = 0;
    indx->scidError = 0;
    indx->sohError = 0;
    indx->bitRateRep = ((double) hdr.errbits) / MAX(hdr.numbits, 1);
    srec = 0;
    erec = nrecs - 1;


    /*									*/
    /* Get bounding times.  For HRPT, all frame timetags should be      */
    /* between AOS and LOS times. GAC/LAC should be between LOS - 1 day */
    /* and LOS time.                                                    */
    /*									*/
    if (hdr.type == HRPT) {
        startSec = hdr.startTime;
        stopSec = hdr.stopTime;
    } else {
        startSec = hdr.startTime - 86400.;
        stopSec = hdr.stopTime + 60.;
    }


    /*									*/
    /* Read through end of input file, load frame index		     	*/
    /*									*/
    printf("\nChecking telemetry for general corruption.\n");
    for (irec = 0; irec < nrecs; irec++) {

        if (fread(mnf, L0LEN, 1, fp) != 1) {
            fprintf(stderr,
                    "-E- %s line %d: error reading %s at rec %d\n",
                    __FILE__, __LINE__, indx->l0file, irec);
            return (FATAL_ERROR);
        }

        /*                                                              */
        /* Initialize error flags for this frame                        */
        /*                                                              */
        indx->rec[irec].bitError = 0;
        indx->rec[irec].scidError = 0;
        indx->rec[irec].tRanError = 0;
        indx->rec[irec].tSeqError = 0;
        indx->rec[irec].tDifError = 0;
        indx->rec[irec].sgaError = 0;
        indx->rec[irec].sacError = 0;
        indx->rec[irec].saaError = 0;

        /*                                                              */
        /* Get frame type, number, and time                             */
        /*                                                              */
        memcpy(scid, &mnf[O_SCID], sizeof (scid));
        memcpy(ttag, &mnf[O_TIME], sizeof (ttag));
        if (endianess() == 1) {
            swapc_bytes((char *) scid, 2, 2);
            swapc_bytes((char *) ttag, 2, 4);
        }
        mnfnum = scid2mnfnum(scid);
        mnftype = scid2mnftype(scid);
        utime = ttag2unix(ttag);

        indx->rec[irec].time = utime;
        indx->rec[irec].mnfnum = mnfnum;
        indx->rec[irec].mnftype = mnftype;

        /* Check for Spacecraft ID field error */
        if (scid[0] != 403 && scid[0] != 147 && scid[0] != 275) {
            indx->rec[irec].scidError = 1;
            printf("Invalid S/C ID, %d, at frame %d.\n", scid[0], irec);
        }
        if (hdr.type == HRPT && mnftype != 0) {
            indx->rec[irec].scidError = 1;
            printf("Invalid frame type for HRPT at frame %d.\n", irec);
        }
        if (hdr.type == GAC && ((mnftype > 4 && mnftype < 15) || mnftype > 15)) {
            indx->rec[irec].scidError = 1;
            printf("Invalid frame type for GAC/LAC at frame %d.\n", irec);
        }

        /*                                                              */
        /* Check for bit errors                                         */
        /*                                                              */
        indx->rec[irec].bitError = bitError(mnf, &numbits, &numerrs);
        indx->rec[irec].numBits = numbits;
        indx->rec[irec].errBits = numerrs;
        totnumbits += numbits;
        totnumerrs += numerrs;
        if (mnf[2] > 0)
            indx->rec[irec].bitError = 1;
        if (indx->rec[irec].bitError == 1)
            printf("Frame %5d, FF bit errors: %d, L0 bit errors: %d\n",
                irec, mnf[2], numerrs);
        if (indx->rec[irec].errBits > l0ctl->maxBitErrors) {
            printf("Max number of bit errors exceeded for this frame.\n");
            indx->rec[irec].maxErrBits = 1;
        }

        /*                                                              */
        /* Check for errors in SOH headers                              */
        /*                                                              */
        if (mnfnum == 1) {
            if (sohHdrError(&mnf[O_SGA]))
                indx->rec[irec].sgaError = 1;
            if (sohHdrError(&mnf[O_SAC]))
                indx->rec[irec].sacError = 1;
            if (sohHdrError(&mnf[O_SAA]))
                indx->rec[irec].saaError = 1;
        }

        /*                                                              */
        /* Compute pixel-to-pixel variance                              */
        /*                                                              */
        indx->rec[irec].pixVariance = pixVariance(mnf);
    }

    /*                                                                   */
    /* If this is HRPT data, the scene should be defined by contiguous   */
    /* frames with a standard time increment, so we can redefine the     */
    /* start and end times based on timetag continuity.                  */
    /*                                                                   */
    if (hdr.type == HRPT) {

        printf("\nComputing AOS and LOS from frame timetags.\n");

        /* Find valid start time */
        timeValid = 0;
        irec = 0;
        while (irec < nrecs && !(timeValid = timeContiguous(indx, irec)))
            irec++;

        if (timeValid == 1) {
            srec = MAX(irec - 1, 0);
            startSec = indx->rec[srec].time;
            indx->timeAOS = startSec;
            printf("Contiguous timetags at frame %d\n", srec);
            for (irec = srec - 1; irec >= 0; irec--)
                indx->rec[irec].tRanError = 1;
        } else {
            fprintf(stderr,
                    "-E- %s line %d: no contiguous timetags found in %s\n",
                    __FILE__, __LINE__, filename);
            exit(FATAL_ERROR);
        }


        /* Find valid stop time */
        timeValid = 0;
        irec = nrecs - 1;
        while ((irec > 0 &&
                (!(timeValid = timeContiguous(indx, irec)) ||
                indx->rec[irec].time > startSec + l0ctl->stopTimeDelta)) ||
                indx->rec[irec].time < startSec)
            irec--;

        if (timeValid) {
            erec = MIN(irec + 1, nrecs - 1);
            stopSec = indx->rec[erec].time;
            indx->timeLOS = stopSec;
            printf("Contiguous timetags at frame %d\n", erec);
            for (irec = erec + 1; irec < nrecs; irec++)
                indx->rec[irec].tRanError = 1;
        } else {
            fprintf(stderr,
                    "-E- %s line %d: no contiguous timetags found in %s\n",
                    __FILE__, __LINE__, filename);
            exit(FATAL_ERROR);
        }
    }
    indx->srec = srec;
    indx->erec = erec;


    /*                                                                   */
    /* If this is GAC data, the scene should be defined by contiguous    */
    /* GAC frames at the start.  We can use that knowledge to reset the  */
    /* time limits for GAC, but we only do that on request of the user.  */
    /*                                                                   */
    if (hdr.type == GAC && l0ctl->timerangeFactor > 0) {

        printf("\nComputing valid timerange from frame timetags.\n\n");

        /* Find valid start time */
        timeValid = 0;
        irec = 0;
        while (irec < nrecs && (timeValid = timeContiguous(indx, irec)) == 0)
            irec++;
        if (timeValid == 1) {
            startSec = indx->rec[irec].time;
            printf("Contiguous timetags at frame %d\n", irec);
        } else {
            fprintf(stderr,
                    "-E- %s line %d: no contiguous timetags found in %s\n",
                    __FILE__, __LINE__, filename);
            exit(FATAL_ERROR);
        }

        stopSec = startSec + l0ctl->timerangeFactor * 86400.0;
        startSec = startSec - l0ctl->timerangeFactor * 86400.0;
        indx->timeAOS = startSec;
        indx->timeLOS = stopSec;
    }


    /*                                                                  */
    /* Now, loop through index and Set frame quality flags              */
    /*                                                                  */
    printf("\nChecking timetag differences, limits, and sequence.\n");
    for (irec = srec; irec <= erec; irec++) {

        mnfnum = indx->rec[irec].mnfnum;
        mnftype = indx->rec[irec].mnftype;

        /* Time range check */
        if (indx->rec[irec].time < startSec ||
                indx->rec[irec].time > stopSec) {
            printf("Time Range Error at Frame %d\n", irec);
            indx->rec[irec].tRanError = 1;
        }

        /* Time difference check */
        if (!timeConsistent(indx, irec)) {
            indx->rec[irec].tDifError = 1;
        }

        /* Time sequence check */
        indx->rec[irec].tSeqError = timeSeqError(indx, irec);
        if (indx->rec[irec].tSeqError == 1)
            printf("Time Sequence Error at Frame %d\n", irec);

        /* Time shift check */
        if (mnftype == GACTYPE)
            indx->rec[irec].tShfError =
                timeShifted(indx, irec, &(indx->rec[irec].timeShift));
        else {
            indx->rec[irec].tShfError = 0;
            indx->rec[irec].timeShift = 0.0;
        }

        /*                                                              */
        /* Update accumulation statistics                               */
        /*                                                              */
        switch (mnftype) {
        case LACTYPE: indx->nlac += 1;
            break;
        case LUNTYPE: indx->nlun += 1;
            break;
        case SOLTYPE: indx->nsol += 1;
            break;
        case IGCTYPE: indx->nigc += 1;
            break;
        case TDITYPE: indx->ntdi += 1;
            break;
        case GACTYPE: indx->ngac += 1;
            break;
        default: indx->nund += 1;
            break;
        }
        indx->scidError += indx->rec[irec].scidError;
        indx->sohError += indx->rec[irec].sgaError ||
                indx->rec[irec].sacError ||
                indx->rec[irec].saaError;
        indx->timeError += indx->rec[irec].tRanError ||
                indx->rec[irec].tSeqError ||
                indx->rec[irec].tDifError;
    }

    /*                                                              */
    /* Compute total bit errors from level-0 data                   */
    /*                                                              */
    indx->bitRateComp = (float) totnumerrs / (float) totnumbits;

    /*                                                               */
    /* Now, make another pass to repeat the time sequence test. This */
    /* test uses only unflagged timetags, so it may catch additional */
    /* problems on a second pass. In addition, we impose a strict    */
    /* requirement that HRPT frames be monotonically increasing.     */
    /*                                                               */
    printf("\nPerforming pass 2 time sequence test.\n");
    lastTime = startSec;
    for (irec = srec; irec <= erec; irec++) {

        if (!timeError(indx, irec)) {
            /* Time sequence check */
            indx->rec[irec].tSeqError = timeSeqError(indx, irec);
            if (indx->rec[irec].tSeqError == 1) {
                printf("Time Sequence Error at Frame %d\n", irec);
                indx->rec[irec].tSeqError = 1;
                indx->timeError++;
            }
        }

        if (indx->type == HRPT && !timeError(indx, irec)) {
            /* Time increasing check */
            if (indx->rec[irec].time < lastTime) {
                printf("Decreasing timetag at frame %d\n", irec);
                indx->rec[irec].tSeqError = 1;
                indx->timeError++;
            } else
                lastTime = indx->rec[irec].time;
        }
    }

    fclose(fp);

    indx->timeFirst = (double) indx->timeLOS;
    indx->timeLast = (double) indx->timeAOS;

    for (irec = srec; irec <= erec; irec++) {
        if (!timeError(indx, irec)) {
            indx->timeFirst = MIN(indx->rec[irec].time, indx->timeFirst);
            indx->timeLast = MAX(indx->rec[irec].time, indx->timeLast);
        }
    }

    return (nrecs);
}
