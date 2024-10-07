/* -------------------------------------------------------------------- */
/* getl1rec() - combines L0 data with preprocessed navigation and       */
/*              converted engineering telemetry to produce one or five  */
/*              L1 records.                                             */
/*                                                                      */
/* Synopsis:                                                            */
/*                                                                      */
/*   nrecs = getl1rec(file,sceneFrameNum,scene,navinp,navblk,tiltblk,   */
/*                    l1rec)                                            */
/*                                                                      */
/*   nrecs          INT32        number of output L1 recs               */
/*   sceneFrameNum  INT32        frame number relative to scene record  */
/*   scene          swl0scene*   pointer to scene record                */
/*   navinp         input_sType* array of navigation input records      */
/*   navblk         input_sType* array of processed navigation blocks   */
/*   tiltblk         tilt_states_sType  tilt state information          */
/*   l1rec          swl1rec*     array of output L1 records             */
/*                                                                      */
/* Description:                                                         */
/*                                                                      */
/* getl1rec returns -1 if it can't read the file, otherwise it returns  */
/* the number of output L1 records.  For GAC data, nrec=5, as there are */
/* 5 scans in a GAC frame.  All other types return nrec=1. The input    */
/* scene record contains, among other things, the list of frame numbers */
/* within the L0 file which are associated with this scene.  So, the    */
/* scene record, incombination with the sceneFrameNum serve as a pointer*/
/* to the L0 frame of the file.  The navinp array contains one record   */
/* for each frame of the scene, but the navblk contains one record for  */
/* each scan of the scene.  So, for GAC frames the number of records in */
/* the navblk will be 5x the number of records in the navinp array. The */
/* primary purpose of this function is to coordinate the combination of */
/* these various input sources.                                         */
/*                                                                      */
/* Modification History:                                                */
/*                                                                      */
/* 22 Oct 1997, B. A. Franz, General Sciences Corp.                     */
/*     - initial development                                            */
/*                                                                      */
/* ---------------------------------------------------------------------*/

#include <stdio.h>
#include <math.h>
#include "swl0_proto.h"

INT32 getl1rec(INT16 sceneFrameNum,
        swl0scene *scene,
        swl0ctl *l0ctl,
        input_sType navinp[],
        navblk_sType navblk[],
        tilt_states_sType *tiltblk,
        swl1rec l1rec[]) {
    static FILE *fp = NULL;

    INT16 fileFrameNum; /* frame number relative to file start   */
    INT16 sceneScanNum; /* scan number relative to scene start   */
    INT32 fileBytePos; /* file position of desired L0 frame     */
    INT16 irec; /* output l1 record index                */
    INT16 nrecs = 1; /* number of output l1 records           */
    BYTE mnf[L0LEN]; /* raw minor frame                       */

    INT32 npix; /* number of pixels per scan             */
    INT32 spix; /* first pixel in LAC scan               */
    INT32 ipix; /* pixel spacing in LAC scan             */

    INT32 o_inst; /* byte offset to first inst tlm block   */
    INT32 o_data; /* byte offset to first image data block */
    INT32 o_gtdipix; /* offset to tdi pixel of image block    */
    INT32 o_startpix; /* offset to start pixel of image block  */
    INT32 o_darkpix; /* offset to first dark restore pixel    */
    INT32 o_image; /* offset to first image pixel           */
    INT32 o_stoppix; /* offset to stop pixel of image block   */
    INT32 s_data; /* size of image data block in bytes     */
    INT32 s_image; /* size of image scan in bytes           */
    INT32 s_inst; /* size of instrument telemetry block    */

    static INT16 lastScanTemp[8] = {0, 0, 0, 0, 0, 0, 0, 0};
    static INT16 lastMirrorSide = 0;
    static INT16 lastTdi[8] = {0, 0, 0, 0, 0, 0, 0, 0};

    INT16 scid[2];
    INT16 mnfnum;
    INT16 mnftype;
    INT16 inst[44];
    INT32 nbits;
    INT32 nbiterrs;

    int endian = 0;
    int i, j;

    if (endianess() == 1)
        endian = 1;

    /*                                                               */
    /* Define any data-type specific offset and size parameters      */
    /*                                                               */
    if (scene->mnftype == GACTYPE) {
        nrecs = 5;
        spix = SPIXGAC;
        ipix = IPIXGAC;
        npix = NPIXGAC;
        o_inst = O_INSTGAC;
        o_data = O_DATAGAC;
        s_data = S_DATAGAC;
    } else {
        nrecs = 1;
        spix = SPIXLAC;
        ipix = IPIXLAC;
        npix = NPIXLAC;
        o_inst = O_INSTLAC;
        o_data = O_DATALAC;
        s_data = S_DATALAC;
    }
    s_image = npix * sizeof (l1rec[0].data[0][0]) * NBANDS;
    s_inst = sizeof (l1rec[0].inst);


    /*                                                                */
    /* Open the L0 file for reading, if it is not already opened.     */
    /* Note:  this assumes only one L0 file per process.              */
    /*                                                                */
    if (fp == NULL)
        if ((fp = fopen(scene->l0file, "r")) == NULL) {
            fprintf(stderr,
                    "-E- %s line %d: unable to open %s for reading\n",
                    __FILE__, __LINE__, scene->l0file);
            return (1);
        }

    /*                                                                */
    /* Read minor frame from L0 file                                  */
    /*                                                                */
    fileFrameNum = scene->indx[sceneFrameNum];
    fileBytePos = sizeof (swl0hdr) + fileFrameNum * (L0LEN);
    if (fseek(fp, fileBytePos, SEEK_SET) != 0) {
        fprintf(stderr,
                "-E- %s line %d: error reading %s\n",
                __FILE__, __LINE__, scene->l0file);
        return (1);
    }
    if (fread(mnf, L0LEN, 1, fp) != 1) {
        fprintf(stderr,
                "-E- %s line %d: error reading %s\n",
                __FILE__, __LINE__, scene->l0file);
        return (1);
    }

    /* Get basic frame attributes */
    memcpy(scid, &mnf[O_SCID], sizeof (scid));
    if (endian)
        swapc_bytes((char *) scid, 2, 2);
    mnfnum = scid2mnfnum(scid);
    mnftype = scid2mnftype(scid);
    if (mnftype != scene->mnftype) {
        fprintf(stderr,
                "-E- %s line %d: frame type mismatch\n",
                __FILE__, __LINE__);
        return (1);
    }

    /*                                                                */
    /* For GAC data, there will be five l1 records generated.  For    */
    /* LAC and LAC-like data, only one l1 record is generated.  Load  */
    /* each output record in sequence.                                */
    /*                                                                */
    for (irec = 0; irec < nrecs; irec++) {

        /* First, clear the record so there's no confusion */
        memset(&l1rec[irec], 0, sizeof (swl1rec));

        l1rec[irec].type = scene->mnftype;
        l1rec[irec].npix = npix;

        /*                                                            */
        /* Grab raw telemetry from minor frame. Note that the ttag,   */
        /* scid, and soh data are only recorded once per frame, so    */
        /* they are replicated for each GAC scan.                     */
        /*                                                            */
        memcpy(l1rec[irec].ttag, &mnf[O_TIME], sizeof (l1rec[0].ttag));
        memcpy(l1rec[irec].scid, &mnf[O_SCID], sizeof (l1rec[0].scid));
        memcpy(l1rec[irec].soh, &mnf[O_SOH], sizeof (l1rec[0].soh));
        memcpy(l1rec[irec].inst, &mnf[o_inst + irec * s_inst], s_inst);

        if (endian) {
            swapc_bytes((char *) l1rec[irec].ttag, 2, 4);
            swapc_bytes((char *) l1rec[irec].scid, 2, 2);
            swapc_bytes((char *) l1rec[irec].inst, 2, 44);
        }

        /*                                                            */
        /* Grab raw image data from minor frame.                      */
        /*                                                            */
        o_gtdipix = o_data + irec * s_data + 0;
        o_startpix = o_data + irec * s_data + 16;
        o_darkpix = o_data + irec * s_data + 32;
        o_image = o_data + irec * s_data + 48;
        o_stoppix = o_data + irec * s_data + s_data - 16;

        memcpy(l1rec[irec].startpix, &mnf[o_startpix],
                sizeof (l1rec[0].startpix));
        memcpy(l1rec[irec].darkpix, &mnf[o_darkpix],
                sizeof (l1rec[0].darkpix));
        memcpy(&l1rec[irec].data[0][0], &mnf[o_image],
                s_image);
        memcpy(l1rec[irec].stoppix, &mnf[o_stoppix],
                sizeof (l1rec[0].stoppix));

        if (endian) {
            swapc_bytes((char *) l1rec[irec].startpix, 2, 8);
            swapc_bytes((char *) l1rec[irec].darkpix, 2, 8);
            swapc_bytes((char *) l1rec[irec].data, 2, s_image / 2);
            swapc_bytes((char *) l1rec[irec].stoppix, 2, 8);
        }

        /*                                                            */
        /* Set gain and tdi.  For HRPT we just force it to nominal    */
        /* values, to avoid problems with noise.                      */
        /*                                                            */
        if (scene->type == HRPT && l0ctl->gainSetting >= 0) {
            for (i = 0; i < NBANDS; i++) {
                l1rec[irec].gain[i] = l0ctl->gainSetting;
                l1rec[irec].tdi[i] = 0;
            }
        } else {
            if (scene->mnftype == TDITYPE || scene->mnftype == IGCTYPE) {
                /* Apply previous TDI setting to current scan */
                memcpy(l1rec[irec].tdi, lastTdi, sizeof (lastTdi));
                memcpy(lastTdi, &mnf[o_gtdipix], sizeof (lastTdi));
            } else {
                memcpy(l1rec[irec].tdi, &mnf[o_gtdipix], sizeof (l1rec[0].tdi));
            }
            if (endian)
                swapc_bytes((char *) l1rec[irec].tdi, 2, 8);

            /* Need to separate the GAIN from the TDI */
            for (i = 0; i < NBANDS; i++) {
                l1rec[irec].gain[i] = l1rec[irec].tdi[i] >> 8;
                l1rec[irec].tdi[i] = l1rec[irec].tdi[i] & (INT16) 255;
            }
        }

        /*                                                            */
        /* Need to determine half-angle mirror side from start pixel. */
        /* If the field is not all 0 or all 1, it is corrupted, so we */
        /* assume the mirror side based on last mirror side.          */
        /*                                                            */
        if (l1rec[irec].startpix[7] == 0)
            l1rec[irec].side = 0;
        else if (l1rec[irec].startpix[7] == 1023)
            l1rec[irec].side = 1;
        else {
            /*
              fprintf(stderr,
                      "-W- %s: start pix corrupted in frame %d\n",
                      __FILE__,fileFrameNum);
             */
            if (scene->mnftype == GACTYPE)
                l1rec[irec].side = lastMirrorSide;
            else
                l1rec[irec].side = (lastMirrorSide + 1) % 2;
        }
        lastMirrorSide = l1rec[irec].side;


        /*                                                            */
        /* Grab converted spacecrft and instrument tlm from nav input.*/
        /* The navigation input is passed in with 1 record per frame, */
        /* so for GAC we will be replicating the sc_ana and sc_dis.   */
        /*                                                            */
        memcpy(l1rec[irec].sc_ana, navinp[sceneFrameNum].sc_ana,
                sizeof (l1rec[0].sc_ana));
        memcpy(l1rec[irec].sc_dis, navinp[sceneFrameNum].sc_dis,
                sizeof (l1rec[0].sc_dis));
        memcpy(l1rec[irec].inst_ana,
                &(navinp[sceneFrameNum].inst_ana[irec][0]),
                sizeof (l1rec[0].inst_ana));
        memcpy(l1rec[irec].inst_dis,
                &(navinp[sceneFrameNum].inst_dis[irec][0]),
                sizeof (l1rec[0].inst_dis));


        /* For GAC scans beyond the first, need to add to msecs       */
        l1rec[irec].msec = (INT32) (navinp[sceneFrameNum].msec +
                (irec * 4 * DTLAC * 1000)) % 86400000;

        /*                                                            */
        /* The navigation block is passed in with 1 record per scan,  */
        /* so for GAC it has already been expanded. Same for tiltblk. */
        /*                                                            */
        sceneScanNum = sceneFrameNum * nrecs + irec;
        memcpy(l1rec[irec].orb_vec, navblk[sceneScanNum].orb_vec,
                sizeof (l1rec[0].orb_vec));
        memcpy(l1rec[irec].l_vert, navblk[sceneScanNum].l_vert,
                sizeof (l1rec[0].l_vert));
        memcpy(l1rec[irec].sun_ref, navblk[sceneScanNum].sun_ref,
                sizeof (l1rec[0].sun_ref));
        memcpy(l1rec[irec].att_ang, navblk[sceneScanNum].att_ang,
                sizeof (l1rec[0].att_ang));
        memcpy(l1rec[irec].sen_mat, navblk[sceneScanNum].sen_mat,
                sizeof (l1rec[0].sen_mat));
        memcpy(l1rec[irec].scan_ell, navblk[sceneScanNum].scan_ell,
                sizeof (l1rec[0].scan_ell));
        memcpy(l1rec[irec].nflag, navblk[sceneScanNum].nflag,
                sizeof (l1rec[0].nflag));
        l1rec[irec].tilt = tiltblk->tilt[sceneScanNum];

        /*                                                            */
        /* If the navigation quality is OK, geolocate the scanline    */
        /* and load the geolocation arrays.                           */
        /*                                                            */
        if (navblk[sceneScanNum].nflag[0] == 0) {
            geonav_(navblk[sceneScanNum].orb_vec,
                    navblk[sceneScanNum].sen_mat,
                    navblk[sceneScanNum].scan_ell,
                    navblk[sceneScanNum].sun_ref,
                    &spix, &ipix, &npix,
                    l1rec[irec].lat,
                    l1rec[irec].lon,
                    l1rec[irec].solz,
                    l1rec[irec].sola,
                    l1rec[irec].senz,
                    l1rec[irec].sena);
        }

        l1rec[irec].slat = MIN(l1rec[irec].lat[0], 90.0);
        l1rec[irec].slon = MIN(l1rec[irec].lon[0], 180.0);
        l1rec[irec].clat = MIN(l1rec[irec].lat[npix / 2], 90.0);
        l1rec[irec].clon = MIN(l1rec[irec].lon[npix / 2], 180.0);
        l1rec[irec].elat = MIN(l1rec[irec].lat[npix - 1], 90.0);
        l1rec[irec].elon = MIN(l1rec[irec].lon[npix - 1], 180.0);
        l1rec[irec].csol_z = l1rec[irec].solz[npix / 2];


        /*                                                            */
        /* Need to copy focal-plane temps from instrument telemetry,  */
        /* but not every scan has valid instrument telemetry, so in   */
        /* that case we use the nearest temperature data available.   */
        /*                                                            */
        if (valid_instlm(mnftype, mnfnum, irec)) {

            l1rec[irec].scan_temp[0] = l1rec[irec].inst[ 7] & 255;
            l1rec[irec].scan_temp[1] = l1rec[irec].inst[ 7] & 255;
            l1rec[irec].scan_temp[2] = l1rec[irec].inst[ 8] & 255;
            l1rec[irec].scan_temp[3] = l1rec[irec].inst[ 8] & 255;
            l1rec[irec].scan_temp[4] = l1rec[irec].inst[ 9] & 255;
            l1rec[irec].scan_temp[5] = l1rec[irec].inst[ 9] & 255;
            l1rec[irec].scan_temp[6] = l1rec[irec].inst[10] & 255;
            l1rec[irec].scan_temp[7] = l1rec[irec].inst[10] & 255;
            for (i = 0; i < NBANDS; i++)
                lastScanTemp[i] = l1rec[irec].scan_temp[i];

            /* set engineering quality flags if the inst tlm is valid.*/
            getEngQual(l1rec[irec].inst_ana, l1rec[irec].eng_qual);

        } else {

            if ((scene->mnftype == GACTYPE) && (irec == 0)) {
                /*                                                    */
                /* If the first scan in a GAC frame does not have     */
                /* valid instrument telemetry, we can use the insttlm */
                /* from the next scan to get scan temp. Otherwise, we */
                /* can rely on the lastScanTemp array.                */
                /*                                                    */
                memcpy(inst, &mnf[o_inst + s_inst], s_inst);
                if (endian)
                    swapc_bytes((char *) inst, 2, 44);


                l1rec[irec].scan_temp[0] = inst[ 7] & 255;
                l1rec[irec].scan_temp[1] = inst[ 7] & 255;
                l1rec[irec].scan_temp[2] = inst[ 8] & 255;
                l1rec[irec].scan_temp[3] = inst[ 8] & 255;
                l1rec[irec].scan_temp[4] = inst[ 9] & 255;
                l1rec[irec].scan_temp[5] = inst[ 9] & 255;
                l1rec[irec].scan_temp[6] = inst[10] & 255;
                l1rec[irec].scan_temp[7] = inst[10] & 255;

            } else
                /*                                                    */
                /* In most cases, the use of the previous scan temps  */
                /* will be perfectly accurate.  If, however, this is  */
                /* the first frame of a LAC-type scene, the last scan */
                /* temp may not be set or may not apply.  This issue  */
                /* is avoided if LAC-type scenes are forced to start  */
                /* on a valid instrument frame.                       */
                /*                                                    */
                for (i = 0; i < NBANDS; i++)
                    l1rec[irec].scan_temp[i] = lastScanTemp[i];
        }


        /*                                                              */
        /* Set the scan-line flags                                      */
        /* 0 : # bit errors                                             */
        /* 3 : # bits tested                                            */
        /* 1 : navigation tlm quality flag                              */
        /* 2 : scan number within frame (1-15)                          */
        /*                                                              */

        bitError(mnf, &nbits, &nbiterrs);
        nbits += (256 * mnf[0] + mnf[1]);
        nbiterrs += mnf[2];
        l1rec[irec].s_flags[0] = (BYTE) MIN(nbiterrs, 255);
        l1rec[irec].s_flags[3] = (BYTE) MIN(nbits / 5, 255);

        l1rec[irec].s_flags[1] = scene->qual[sceneFrameNum];

        if (scene->mnftype == GACTYPE)
            l1rec[irec].s_flags[2] = (mnfnum - 1) * 5 + irec + 1;
        else
            l1rec[irec].s_flags[2] = 1;


        /*                                                              */
        /* Get saturated and zero-pixel counts                          */
        /*                                                              */
        for (i = 0; i < npix; i++)
            for (j = 0; j < NBANDS; j++) {
                if (l1rec[irec].data[i][j] >= 1023)
                    l1rec[irec].s_satp[j]++;
                if (l1rec[irec].data[i][j] - l1rec[irec].darkpix[j] <= 0) {
                    l1rec[irec].s_zerop[j]++;
                }
            }

    }


    return (nrecs);
}




