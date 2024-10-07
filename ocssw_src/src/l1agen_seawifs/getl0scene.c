/* ============================================================== */
/* Module:  getl0scene.c                                          */
/* Purpose: Generate index of contiguous scenes with the L0 file  */
/* Author:  B.A. Franz, General Scences Corp., 9/97               */
/* ============================================================== */

#include <string.h>
#include "swl0_proto.h"

INT16 goodFrame(swl0indx *indx, INT16 irec);
INT16 goodTLM(swl0indx *indx, INT16 irec);
INT16 newScene(INT16 thisRec, INT16 lastRec, swl0indx *indx);

INT16 getl0scene(swl0indx *indx, swl0scene scene[]) {
    INT16 nscenes = 0;
    INT16 irec;
    INT16 wantFrame;
    INT16 lastRec;
    INT16 is;
    INT16 firstInstFound;

    /*									*/
    /* HRPT files will be held to one scene				*/
    /*									*/

    if (indx->type == HRPT) {

        nscenes = 1;
        scene[0].nrec = 0;
        scene[0].type = 1;
        firstInstFound = 0;

        for (irec = 0; irec < indx->nrecs; irec++) {

            wantFrame = (indx->rec[irec].mnftype == LACTYPE) &&
                    goodFrame(indx, irec);

            if (wantFrame && (indx->rec[irec].mnfnum != 2))
                firstInstFound = 1;

            /* Delay start rec for first valid inst telemetry */
            if (!firstInstFound)
                wantFrame = 0;

            if (wantFrame) {
                scene[0].indx[scene[0].nrec] = irec;
                if (goodTLM(indx, irec))
                    scene[0].qual[scene[0].nrec] = 0;
                else
                    scene[0].qual[scene[0].nrec] = 1;
                scene[0].nrec++;
            }
        }

        if (scene[0].nrec < MINFRAMES)
            nscenes = 0;

    } else {


        nscenes = -1;
        lastRec = -1;

        for (irec = 0; irec < indx->nrecs; irec++) {

            if (goodFrame(indx, irec)) {

                if (newScene(irec, lastRec, indx)) {
                    nscenes++;
                    if (nscenes > 0 && scene[nscenes - 1].nrec < MINFRAMES) {
                        nscenes--;
                        printf("getl0scene: scene %d of size %d is too small.\n",
                                nscenes, scene[nscenes].nrec);
                    }
                    if (nscenes >= MAXSCENES) {
                        printf("getl0scene: Warning: Max Scene Limit Exceeded\n");
                        break;
                    }
                    scene[nscenes].nrec = 0;
                    firstInstFound = 0;
                    lastRec = irec;
                }

                /* If LAC, delay start of scene for first inst telemetry */
                if ((indx->rec[irec].mnfnum != 2) ||
                        (indx->rec[irec].mnftype == GACTYPE))
                    firstInstFound = 1;

                if (firstInstFound) {
                    scene[nscenes].indx[scene[nscenes].nrec] = irec;
                    if (goodTLM(indx, irec))
                        scene[nscenes].qual[scene[nscenes].nrec] = 0;
                    else
                        scene[nscenes].qual[scene[nscenes].nrec] = 1;
                    scene[nscenes].nrec++;
                    lastRec = irec;
                }
            }
        }

        nscenes++;
        if (nscenes > 0 && scene[nscenes - 1].nrec < MINFRAMES) {
            nscenes--;
            fprintf(stderr,
                    "-W- %s line %d: scene %d of size %d is too small.\n",
                    __FILE__, __LINE__, nscenes, scene[nscenes].nrec);
        }

    }


    for (is = 0; is < nscenes; is++) {
        strncpy(scene[is].l0file, indx->l0file, FILENAME_MAX);
        scene[is].srec = scene[is].indx[ 0 ];
        scene[is].erec = scene[is].indx[ scene[is].nrec - 1 ];
        scene[is].crec = scene[is].indx[ scene[is].nrec / 2 ];
        scene[is].stime = indx->rec[ scene[is].srec ].time;
        scene[is].etime = indx->rec[ scene[is].erec ].time;
        scene[is].ctime = indx->rec[ scene[is].crec ].time;
        scene[is].mnftype = indx->rec[ scene[is].srec ].mnftype;
        scene[is].type = indx->type;
        scene[is].nscan = scene[is].nrec;
        if (scene[is].mnftype == GACTYPE) { /* 5 GAC scans/frame    */
            scene[is].nscan *= 5;
            scene[is].etime += (DTGAC / 5 * 4);
        }
        scene[is].orbnum = getorbnum(scene[is].ctime);
    }

    return nscenes;
}


/* ----------------------------------------------------------------- */
/* goodFrame() - returns 1 if critical frame info is good.           */

/* ----------------------------------------------------------------- */
INT16 goodFrame(swl0indx *indx, INT16 irec) {
    INT16 status = 0;

    /* Don't want any timing or scid errors */
    status = !(indx->rec[irec].tRanError ||
            indx->rec[irec].tSeqError ||
            indx->rec[irec].tDifError ||
            indx->rec[irec].scidError ||
            indx->rec[irec].maxErrBits);

    return (status);
}


/* ----------------------------------------------------------------- */
/* goodTLM() - returns 1 if SOH telemetry looks OK.                  */

/* ----------------------------------------------------------------- */
INT16 goodTLM(swl0indx *indx, INT16 irec) {
    INT16 status = 1;

    /* 
     * Only care about mnf #1, since that is where we find the
     * critical navigation telemetry.
     */
    if (indx->rec[irec].mnfnum == 1)
        status = !(indx->rec[irec].bitError ||
            indx->rec[irec].sgaError ||
            indx->rec[irec].saaError ||
            indx->rec[irec].sacError);

    return (status);
}


/* ----------------------------------------------------------------- */
/* newScene() - returns 1 if the current record should start a new   */
/*              scene.                                               */

/* ----------------------------------------------------------------- */
INT16 newScene(INT16 thisRec, INT16 lastRec, swl0indx *indx) {
    FLOAT64 delTime;

    if (lastRec < 0) {
        printf("First scene located.\n");
        return (1);
    }

    if (indx->rec[thisRec].mnftype != indx->rec[lastRec].mnftype) {
        printf("Scene break for frame type change: %3d %3d\n",
                indx->rec[thisRec].mnftype, indx->rec[lastRec].mnftype);
        return (1);
    }

    delTime = indx->rec[thisRec].time - indx->rec[lastRec].time;

    if ((indx->rec[thisRec].mnftype == GACTYPE) &&
            (delTime < 0 || delTime > DELSCENEGAC)) {
        printf("Scene break for GAC time difference: %lf\n", delTime);
        return (1);
    }


    if ((indx->rec[thisRec].mnftype != GACTYPE) &&
            (delTime < 0 || delTime > DELSCENELAC)) {
        printf("Scene break for LAC time difference: %lf\n", delTime);
        return (1);
    }

    return (0);
}


