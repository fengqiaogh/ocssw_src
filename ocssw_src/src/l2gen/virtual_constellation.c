#include "l12_proto.h"
#include "chl.h"
#include "fann.h"

#define NBANDS_VC  7
static int wave_vc[] = {412, 443, 490, 510, 531, 555, 670};
static float badval = BAD_FLT;
static float Rrs_norm = 0.025;
static float ratio_531_norm = 2.;

float Rrs_vc_seawifs(l2str *l2rec, int32_t ip, int32_t ibvc) {
    static int firstCall = 1;

    int ib;
    int ipb = ip * l2rec->l1rec->l1file->nbands;
    static int iwave_nn[] = {412, 443, 490, 510, 555, 670};
    float *Rrs = &l2rec->Rrs[ipb];
    fann_type Rrs_in[4];
    float Rrs_vc = badval;

    static struct fann *nntab;

    if (firstCall) {
        // load tables
        char nnFile [FILENAME_MAX] = "";
        sprintf(nnFile, "%s%s%d%s", input->vcnnfile, "_", wave_vc[4], ".dat");
        nntab = fann_create_from_file(nnFile);
        firstCall = 0;
    }

    // build NN input Rrs
    for (ib = 1; ib < 5; ib++) {
        ipb = ip * l2rec->l1rec->l1file->nbands + bindex_get(iwave_nn[ib]);
        Rrs_in[ib - 1] = l2rec->Rrs[ipb] / Rrs_norm;
    }

    // call conversion function

    switch (wave_vc[ibvc]) {
    case 412:
        Rrs_vc = Rrs[0];
        break;
    case 443:
        Rrs_vc = Rrs[1];
        break;
    case 490:
        Rrs_vc = Rrs[2];
        break;
    case 510:
        Rrs_vc = Rrs[3];
        break;
    case 531:
        Rrs_vc = (float) *fann_run(nntab, Rrs_in) * ratio_531_norm * Rrs[1];
        break;
    case 555:
        Rrs_vc = Rrs[4];
        break;
    case 670:
        Rrs_vc = Rrs[5];
        break;
    }

    return (Rrs_vc);
}

float Rrs_vc_modisa(l2str *l2rec, int32_t ip, int32_t ibvc) {
    static int firstCall = 1;
    static int nband_nn = 7;
    static int iwave_nn[] = {412, 443, 488, 531, 547, 667};
    static struct fann * nntab[NBANDS_VC];

    int ib;
    int ipb;
    fann_type *Rrs_in;
    float Rrs_vc = badval;

    if ((Rrs_in = (fann_type *) calloc(l2rec->l1rec->l1file->nbands, sizeof (fann_type))) == NULL) {
        printf("-E- : Error allocating memory to Rrs_vc_modisa\n");
        exit(FATAL_ERROR);
    }

    if (firstCall) {
        // load tables
        char nnFile [FILENAME_MAX] = "";
        for (ib = 0; ib < nband_nn; ib++) {
            sprintf(nnFile, "%s%s%d%s", input->vcnnfile, "_", wave_vc[ib], ".dat");
            nntab[ib] = fann_create_from_file(nnFile);
        }
        firstCall = 0;
    }

    // build NN input Rrs
    for (ib = 0; ib < nband_nn - 1; ib++) {
        ipb = ip * l2rec->l1rec->l1file->nbands + bindex_get(iwave_nn[ib]);
        Rrs_in[ib] = l2rec->Rrs[ipb] / Rrs_norm;
        ;
    }

    // call conversion function

    switch (wave_vc[ibvc]) {
    case 412:
    case 443:
    case 490:
    case 510:
    case 555:
    case 670:
        Rrs_vc = (float) *fann_run(nntab[ibvc], Rrs_in);
        Rrs_vc *= Rrs_norm;
        break;
    case 531:
        Rrs_vc = Rrs_in[3] / Rrs_in[1];
        Rrs_vc *= (float) *fann_run(nntab[ibvc], Rrs_in);
        Rrs_vc *= Rrs_norm;
        break;
    }
    free(Rrs_in);
    return (Rrs_vc);
}

float Rrs_vc_modist(l2str *l2rec, int32_t ip, int32_t ibvc) {
    static int firstCall = 1;
    static int nband_nn = 7;
    static int iwave_nn[] = {412, 443, 488, 531, 547, 667};
    static struct fann * nntab[NBANDS_VC];

    int ib;
    int ipb;
    fann_type *Rrs_in;
    float Rrs_vc = badval;

    if (firstCall) {
        // load tables
        char nnFile [FILENAME_MAX] = "";
        for (ib = 0; ib < nband_nn; ib++) {
            sprintf(nnFile, "%s%s%d%s", input->vcnnfile, "_", wave_vc[ib], ".dat");
            nntab[ib] = fann_create_from_file(nnFile);
        }
        firstCall = 0;
    }
    if ((Rrs_in = (fann_type *) calloc(l2rec->l1rec->l1file->nbands, sizeof (fann_type))) == NULL) {
        printf("-E- %s line %d : error allocating memory for Rrs_in in Rrs_vc_modist.\n",
                __FILE__, __LINE__);
        exit(1);
    }

    // build NN input Rrs
    for (ib = 0; ib < nband_nn - 1; ib++) {
        ipb = ip * l2rec->l1rec->l1file->nbands + bindex_get(iwave_nn[ib]);
        Rrs_in[ib] = l2rec->Rrs[ipb] / Rrs_norm;
        ;
    }

    // call conversion function

    switch (wave_vc[ibvc]) {
    case 412:
    case 443:
    case 490:
    case 510:
    case 555:
    case 670:
        Rrs_vc = (float) *fann_run(nntab[ibvc], Rrs_in);
        Rrs_vc *= Rrs_norm;
        break;
    case 531:
        Rrs_vc = Rrs_in[3] / Rrs_in[1];
        Rrs_vc *= (float) *fann_run(nntab[ibvc], Rrs_in);
        Rrs_vc *= Rrs_norm;
        break;
    }

    free(Rrs_in);
    return (Rrs_vc);
}

float Rrs_vc_viirs(l2str *l2rec, int32_t ip, int32_t ibvc) {
    static int firstCall = 1;
    static int nband_nn = 7;
    static int iwave_nn[] = {410, 443, 486, 551, 671};
    static struct fann * nntab[NBANDS_VC];

    int ib;
    int ipb;
    fann_type *Rrs_in;
    fann_type Rrs_531_in[4];

    float Rrs_vc = badval;

    if (firstCall) {
        // load tables
        char nnFile [FILENAME_MAX] = "";
        for (ib = 0; ib < nband_nn; ib++) {
            sprintf(nnFile, "%s%s%d%s", input->vcnnfile, "_", wave_vc[ib], ".dat");
            nntab[ib] = fann_create_from_file(nnFile);
        }
        firstCall = 0;
    }
    if ((Rrs_in = (fann_type *) calloc(l2rec->l1rec->l1file->nbands, sizeof (fann_type))) == NULL) {
        printf("-E- %s line %d : error allocating memory for Rrs_in in Rrs_vc_modist.\n",
                __FILE__, __LINE__);
        exit(1);
    }

    // build NN input Rrs
    for (ib = 0; ib < nband_nn - 2; ib++) {
        ipb = ip * l2rec->l1rec->l1file->nbands + bindex_get(iwave_nn[ib]);
        Rrs_in[ib] = l2rec->Rrs[ipb] / Rrs_norm;
        if (ib > 0 && ib < 5)
            Rrs_531_in[ib - 1] = l2rec->Rrs[ipb] / Rrs_norm;
    }

    // call conversion function

    switch (wave_vc[ibvc]) {
    case 412:
    case 443:
    case 490:
    case 510:
    case 555:
    case 670:
        Rrs_vc = (float) *fann_run(nntab[ibvc], Rrs_in) * Rrs_norm;
        break;
    case 531:
        Rrs_vc = (float) *fann_run(nntab[1], Rrs_in) * Rrs_norm;
        Rrs_vc *= (float) *fann_run(nntab[ibvc], Rrs_531_in) * ratio_531_norm;
        break;
    }
    free(Rrs_in);
    return (Rrs_vc);
}

float Rrs_vc_convert(l2str *l2rec, int ip, int ib) {
    switch (l2rec->l1rec->l1file->sensorID) {
    case SEAWIFS:
        return (Rrs_vc_seawifs(l2rec, ip, ib));
        break;
        break;
    case MODISA:
        return (Rrs_vc_modisa(l2rec, ip, ib));
        break;
        break;
    case MODIST:
        return (Rrs_vc_modist(l2rec, ip, ib));
        break;
        break;
    case VIIRSN:
    case VIIRSJ1:
    case VIIRSJ2:
        return (Rrs_vc_viirs(l2rec, ip, ib));
        break;
        break;
    default:
        printf("%s Line %d: VC Rrs conversion not defined for this sensor.\n", __FILE__, __LINE__);
        exit(1);
        break;
    }
}

float chl_vc(l2str *l2rec, int32_t ip) {
    static int firstCall = 1;
    static float a[] = {0.2515, -2.3798, 1.5823, -0.6372, -0.5692};
    static int ib1 = 1; // 443
    static int ib2 = 2; // 490
    static int ib3 = 5; // 555

    float rat, minRrs;
    float Rrs1, Rrs2, Rrs3;
    float chl = chlbad;

    if (firstCall) {
        firstCall = 0;
    }

    /*                                                      */
    /* Compute desired products at each pixel               */
    /*                                                      */
    Rrs1 = Rrs_vc_convert(l2rec, ip, ib1);
    Rrs2 = Rrs_vc_convert(l2rec, ip, ib2);
    Rrs3 = Rrs_vc_convert(l2rec, ip, ib3);

    minRrs = MIN(Rrs1, Rrs2);

    if (Rrs3 > 0.0 && Rrs2 > 0.0 && minRrs > -0.001) {
        rat = MAX(Rrs1, Rrs2) / Rrs3;
        if (rat > minrat && rat < maxrat) {
            rat = log10(rat);
            chl = (float) pow(10.0, (a[0] + rat * (a[1] + rat * (a[2] + rat * (a[3] + rat * a[4])))));
            chl = (chl > chlmin ? chl : chlmin);
            chl = (chl < chlmax ? chl : chlmax);
        }
    }

    return (chl);
}

void virtual_constellation(l2str *l2rec, l2prodstr *p, float prod[]) {
    int32_t prodnum = p->cat_ix;
    int32_t ibvc = p->prod_ix;
    int32_t ip;

    /*                                                      */
    /* Compute desired products at each pixel               */
    /*                                                      */
    for (ip = 0; ip < l2rec->l1rec->npix - 1; ip++) {
        switch (prodnum) {
        case CAT_Rrs_vc:
            prod[ip] = Rrs_vc_convert(l2rec, ip, ibvc);
            break;
        case CAT_chl_vc:
            prod[ip] = chl_vc(l2rec, ip);
            break;
        }
    }

}
