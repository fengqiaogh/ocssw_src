/* -------------------------------------------------------------------- */
/* getnavdata() - loads telemetry structure with raw navigation inputs  */
/* -------------------------------------------------------------------- */

#include <stdio.h>
#include "swl0_proto.h"

INT32 getnavdata(swl0scene *scene, input_sType navinp[]) {
    BYTE mnf[L0LEN];
    INT16 ttag[4];
    INT16 scid[2];
    INT16 mnftype;
    INT16 mnfnum;
    FILE *fp;
    INT32 irec;
    INT32 i;
    int j;
    int ninst; /* number of instrument telemetry block    */
    int oinst; /* byte offset to first inst tlm block     */
    int sinst = 88; /* byte length of one inst tlm block       */
    timeTab *trec; /* temporal anomaly table entry            */
    int endian = 0;

    if (endianess() == 1)
        endian = 1;

    if (scene->mnftype == GACTYPE) {
        ninst = 5;
        oinst = O_INSTGAC;
    } else {
        ninst = 1;
        oinst = O_INSTLAC;
    }

    memset(navinp, 0, sizeof (input_sType) * scene->nrec);

    /*									*/
    /* Open file for reading						*/
    /*									*/
    if ((fp = fopen(scene->l0file, "r")) == NULL) {
        fprintf(stderr,
                "-E- %s line %d: unable to open %s for reading\n",
                __FILE__, __LINE__, scene->l0file);
        return (1);
    }


    /*									*/
    /* Loop through scene record index and load	nav input structure     */
    /*									*/
    for (i = 0; i < scene->nrec; i++) {

        irec = scene->indx[i];

        if (fseek(fp, sizeof (swl0hdr) + irec * (L0LEN), SEEK_SET) != 0) {
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

        /* Set telemetry quality flag */
        navinp[i].flag = scene->qual[i];

        /* Get basic frame attributes */
        memcpy(scid, &mnf[O_SCID], sizeof (scid));
        if (endian)
            swapc_bytes((char *) scid, 2, 2);
        mnfnum = scid2mnfnum(scid);
        mnftype = scid2mnftype(scid);

        /* This should never happen if scene index is OK */
        if (mnftype != scene->mnftype) {
            fprintf(stderr,
                    "-E- %s line %d: frame type mismatch\n",
                    __FILE__, __LINE__);
            return (1);
        }

        memcpy(navinp[i].sc_id, &mnf[O_SCID], 4);
        if (endian)
            swapc_bytes((char *) navinp[i].sc_id, 2, 2);

        memcpy(ttag, &mnf[O_TIME], 8);
        if (endian)
            swapc_bytes((char *) ttag, 2, 4);
        ttag2ydmsec(ttag,
                &navinp[i].iyear, &navinp[i].iday, (INT32 *) & navinp[i].msec);

        /* Correct for time shifts */
        trec = temporal_anomaly(mnftype, navinp[i].iyear, navinp[i].iday, navinp[i].msec);
        if (trec != NULL) {
            if (trec->action == 1) {
                printf("Correcting frame %d of time %d %d %d by %d msecs\n",
                        irec, navinp[i].iyear, navinp[i].iday, navinp[i].msec, trec->delmsec);
                navinp[i].msec -= trec->delmsec;
                navinp[i].nflag[5] = 1;
                if (i == 0)
                    scene->stime -= trec->delmsec / 1000.0;
                else if (i == scene->nrec - 1)
                    scene->etime -= trec->delmsec / 1000.0;
            } else if (trec->action == 2) {
                printf("Failing frame %d of time %d %d %d for time correction ambiguity\n",
                        irec, navinp[i].iyear, navinp[i].iday, navinp[i].msec);
                navinp[i].nflag[5] = 1;
                navinp[i].nflag[0] = 1;
            }
        }

        /* Convert SOH telemetry */
        conv_soh_(&mnf[O_SOH], navinp[i].sc_ana, navinp[i].sc_dis);

        /* Convert instrument tlm, if it is valid */
        for (j = 0; j < ninst; j++)
            if (valid_instlm(mnftype, mnfnum, j))
                conv_ins_(&mnf[oinst + sinst * j ],
                    &(navinp[i].inst_ana[j][0]),
                    &(navinp[i].inst_dis[j][0]));
    }

    fclose(fp);

    return (i);

}




