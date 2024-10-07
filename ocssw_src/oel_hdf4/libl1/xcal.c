#include "xcal.h"

#include <stdio.h>
#include <stdlib.h>

#include <mfhdf.h>

static int *xcal_wave;
static int nband = 0;
static int nwave;
static double pixwt = 1.0;
static int detfac = 1;

int xcal_band(int wave) {
    int iband = -1;
    int iw;
    static int firstCall = 1;

    if (firstCall) {
        if ((xcal_wave = (int *) calloc(nwave, sizeof (int))) == NULL) {
            printf("-E- : Error allocating memory to run_soa_sma\n");
            exit(EXIT_FAILURE);
        }
        firstCall = 0;
    }
    for (iw = 0; iw < nband; iw++)
        if (xcal_wave[iw] == wave)
            iband = iw;

    return (iband);
}

double *get_xcal(l1str *l1rec, int type, int wave) {
    typedef double m_array[XTNMSIDE][XTNDET][MAXPIX];
    static m_array *rvs, *m12, *m13;

    int bandnum;

    static int firstCall = 1;
    int32_t ic, ip, im, id, it, it1, it2;
    int32_t ncoef = 0, ntime = 0, ndet = 0, nmside = 0;

    if (firstCall == 1) {
        nwave = l1rec->l1file->nbands;
        if ((m12 = (m_array *) calloc(nwave, sizeof (m_array))) == NULL) {
            printf("-E- : Error allocating memory to dray_for_i\n");
            exit(EXIT_FAILURE);
        }
        if ((m13 = (m_array *) calloc(nwave, sizeof (m_array))) == NULL) {
            printf("-E- : Error allocating memory to dray_for_i\n");
            exit(EXIT_FAILURE);
        }
        if ((rvs = (m_array *) calloc(nwave, sizeof (m_array))) == NULL) {
            printf("-E- : Error allocating memory to dray_for_i\n");
            exit(EXIT_FAILURE);
        }
        firstCall = 0;
    }
    if ((bandnum = xcal_band(wave)) < 0) {

        double rvstab[XTNORDER][XTNDET][XTNMSIDE][XTNTIME];
        double m12tab[XTNORDER][XTNDET][XTNMSIDE][XTNTIME];
        double m13tab[XTNORDER][XTNDET][XTNMSIDE][XTNTIME];
        int16_t yrtab [XTNTIME];
        int16_t dytab [XTNTIME];
        double sectab[XTNTIME];

        int32_t sensorID = l1rec->l1file->sensorID;
        double rvsc1[XTNORDER];
        double rvsc2[XTNORDER];
        double m12c1[XTNORDER];
        double m12c2[XTNORDER];
        double m13c1[XTNORDER];
        double m13c2[XTNORDER];
        double ap, ap2, ap3, ap4, ap5, x1, x2;
        double sec, wt;

        char *file;
        file = (char *) calloc(FILENAME_MAX, sizeof (char)); //[FILENAME_MAX] = "";
        char *name;
        name = (char *) calloc(H4_MAX_NC_NAME, sizeof (char)); //[H4_MAX_NC_NAME]  = "";
        char *sdsname;
        sdsname = (char *) calloc(H4_MAX_NC_NAME, sizeof (char)); //[H4_MAX_NC_NAME]  = "";
        int32_t sd_id;
        int32_t rank;
        int32_t nt;
        int32_t dims[H4_MAX_VAR_DIMS];
        int32_t nattrs;
        int32_t sds_id;
        int32_t start[5] = {0, 0, 0, 0, 0};
        int32_t end [5] = {1, 1, 1, 1, 1};
        int32_t status;

        bandnum = nband++;
        xcal_wave[bandnum] = wave;

        sprintf(file, "%s%s%d%s", l1_input->xcal_file, "_", wave, ".hdf");

        printf("Loading XCAL rvs and polarization sensitivities from %s\n", file);

        sd_id = SDstart(file, DFACC_RDONLY);
        if (sd_id == FAIL) {
            fprintf(stderr, "-E- %s line %d: SDstart(%s, %d) failed.\n",
                    __FILE__, __LINE__, file, DFACC_RDONLY);
            exit(1);
        }

        status = SDreadattr(sd_id, SDfindattr(sd_id, "Number of Coefficients"), &ncoef);
        if (status != 0) {
            printf("-E- %s Line %d:  Error reading global attribute %s from %s.\n",
                    __FILE__, __LINE__, "Number of Coefficients", file);
            exit(1);
        }
        if (ncoef > XTNORDER) {
            printf("-E- %s Line %d:  mismatch in ncoef\n", __FILE__, __LINE__);
            exit(1);
        }

        status = SDreadattr(sd_id, SDfindattr(sd_id, "Number of Detectors"), &ndet);
        if (status != 0) {
            printf("-E- %s Line %d:  Error reading global attribute %s from %s.\n",
                    __FILE__, __LINE__, "Number of Detectors", file);
            exit(1);
        }
        if (ndet > XTNDET) {
            printf("-E- %s Line %d:  mismatch in ndet\n", __FILE__, __LINE__);
            exit(1);
        }

        status = SDreadattr(sd_id, SDfindattr(sd_id, "Number of Mirror Sides"), &nmside);
        if (status != 0) {
            printf("-E- %s Line %d:  Error reading global attribute %s from %s.\n",
                    __FILE__, __LINE__, "Number of Mirror Sides", file);
            exit(1);
        }
        if (nmside > XTNMSIDE) {
            printf("-E- %s Line %d:  mismatch in nmside\n", __FILE__, __LINE__);
            exit(1);
        }

        status = SDreadattr(sd_id, SDfindattr(sd_id, "Length of Time Series ndates"), &ntime);
        if (status != 0) {
            printf("-E- %s Line %d:  Error reading global attribute %s from %s.\n",
                    __FILE__, __LINE__, "Length of Time Series ndates", file);
            exit(1);
        }
        if (ntime > XTNTIME) {
            printf("-E- %s Line %d:  mismatch in ntime\n", __FILE__, __LINE__);
            exit(1);
        }

        strcpy(sdsname, "year");
        start[0] = 0;
        sds_id = SDselect(sd_id, SDnametoindex(sd_id, sdsname));
        status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
        status = SDreaddata(sds_id, start, NULL, dims, (VOIDP) yrtab);
        if (status != 0) {
            printf("-E- %s Line %d:  Error reading SDS %s from %s.\n",
                    __FILE__, __LINE__, sdsname, file);
            exit(1);
        } else {
            status = SDendaccess(sds_id);
        }

        strcpy(sdsname, "day");
        start[0] = 0;
        sds_id = SDselect(sd_id, SDnametoindex(sd_id, sdsname));
        status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
        status = SDreaddata(sds_id, start, NULL, dims, (VOIDP) dytab);
        if (status != 0) {
            printf("-E- %s Line %d:  Error reading SDS %s from %s.\n",
                    __FILE__, __LINE__, sdsname, file);
            exit(1);
        } else {
            status = SDendaccess(sds_id);
        }

        strcpy(sdsname, "M11");
        sds_id = SDselect(sd_id, SDnametoindex(sd_id, sdsname));
        for (ic = 0; ic < ncoef; ic++) for (id = 0; id < ndet; id++) for (im = 0; im < nmside; im++) for (it = 0; it < ntime; it++) {
                        start[0] = ic;
                        start[1] = id;
                        start[2] = im;
                        start[3] = it;
                        status = SDreaddata(sds_id, start, NULL, end, (VOIDP) & rvstab[ic][id][im][it]);
                        if (status != 0) {
                            printf("-E- %s Line %d:  Error reading SDS %s from %s.\n",
                                    __FILE__, __LINE__, sdsname, file);
                            exit(1);
                        }
                    }
        status = SDendaccess(sds_id);

        strcpy(sdsname, "m12");
        sds_id = SDselect(sd_id, SDnametoindex(sd_id, sdsname));
        for (ic = 0; ic < ncoef; ic++) for (id = 0; id < ndet; id++) for (im = 0; im < nmside; im++) for (it = 0; it < ntime; it++) {
                        start[0] = ic;
                        start[1] = id;
                        start[2] = im;
                        start[3] = it;
                        status = SDreaddata(sds_id, start, NULL, end, (VOIDP) & m12tab[ic][id][im][it]);
                        if (status != 0) {
                            printf("-E- %s Line %d:  Error reading SDS %s from %s.\n",
                                    __FILE__, __LINE__, sdsname, file);
                            exit(1);
                        }
                    }
        status = SDendaccess(sds_id);

        strcpy(sdsname, "m13");
        sds_id = SDselect(sd_id, SDnametoindex(sd_id, sdsname));
        for (ic = 0; ic < ncoef; ic++) for (id = 0; id < ndet; id++) for (im = 0; im < nmside; im++) for (it = 0; it < ntime; it++) {
                        start[0] = ic;
                        start[1] = id;
                        start[2] = im;
                        start[3] = it;
                        status = SDreaddata(sds_id, start, NULL, end, (VOIDP) & m13tab[ic][id][im][it]);
                        if (status != 0) {
                            printf("-E- %s Line %d:  Error reading SDS %s from %s.\n",
                                    __FILE__, __LINE__, sdsname, file);
                            exit(1);
                        }
                    }
        status = SDendaccess(sds_id);

        status = SDend(sd_id);


        for (it = 0; it < ntime; it++) {
            if (yrtab[it] > 0)
                sectab[it] = yds2unix(yrtab[it], dytab[it], 0.0);
            else
                break;
        }
        if (it != ntime) {
            printf("-E- %s line %d: time field mismatch %d %d\n", __FILE__, __LINE__, it, ntime);
            exit(1);
        }

        sec = l1rec->scantime;
        for (it = 0; it < ntime; it++) {
            if (sectab[it] > sec) break;
        }
        if (it == 0)
            it2 = it1 = it;
        else if (it == ntime)
            it2 = it1 = ntime - 1;
        else {
            it2 = MAX(MIN(it, ntime - 1), 1);
            it1 = it2 - 1;
        }

        for (im = 0; im < nmside; im++) for (id = 0; id < ndet; id++) {

                for (ic = 0; ic < ncoef; ic++) {
                    rvsc1[ic] = rvstab[ic][id][im][it1];
                    m12c1[ic] = m12tab[ic][id][im][it1];
                    m13c1[ic] = m13tab[ic][id][im][it1];
                    rvsc2[ic] = rvstab[ic][id][im][it2];
                    m12c2[ic] = m12tab[ic][id][im][it2];
                    m13c2[ic] = m13tab[ic][id][im][it2];
                }
                for (ic = ncoef; ic < XTNORDER; ic++) {
                    rvsc1[ic] = 0;
                    m12c1[ic] = 0;
                    m13c1[ic] = 0;
                    rvsc2[ic] = 0;
                    m12c2[ic] = 0;
                    m13c2[ic] = 0;
                }

                if (it2 == it1)
                    wt = 0.0;
                else
                    wt = (sec - sectab[it1]) / (sectab[it2] - sectab[it1]);

                if (sensorID == MODISA || sensorID == MODIST) {
                    switch (l1_input->resolution) {
                    case 250:
                        pixwt = 0.25;
                        detfac = 4;
                        break;
                    case 500:
                        pixwt = 0.5;
                        detfac = 2;
                        break;
                    }
                }

                for (ip = 0; ip < l1rec->npix; ip++) {

                    ap = (double) l1rec->pixnum[ip] * pixwt;
                    ap2 = ap*ap;
                    ap3 = ap2*ap;
                    ap4 = ap3*ap;
                    ap5 = ap4*ap;

                    x1 = rvsc1[0] + ap * rvsc1[1] + ap2 * rvsc1[2] + ap3 * rvsc1[3] + ap4 * rvsc1[4] + ap5 * rvsc1[5];
                    x2 = rvsc2[0] + ap * rvsc2[1] + ap2 * rvsc2[2] + ap3 * rvsc2[3] + ap4 * rvsc2[4] + ap5 * rvsc2[5];
                    rvs[bandnum][im][id][ip] = x1 + wt * (x2 - x1);

                    x1 = m12c1[0] + ap * m12c1[1] + ap2 * m12c1[2] + ap3 * m12c1[3] + ap4 * m12c1[4] + ap5 * m12c1[5];
                    x2 = m12c2[0] + ap * m12c2[1] + ap2 * m12c2[2] + ap3 * m12c2[3] + ap4 * m12c2[4] + ap5 * m12c2[5];
                    m12[bandnum][im][id][ip] = x1 + wt * (x2 - x1);

                    x1 = m13c1[0] + ap * m13c1[1] + ap2 * m13c1[2] + ap3 * m13c1[3] + ap4 * m13c1[4] + ap5 * m13c1[5];
                    x2 = m13c2[0] + ap * m13c2[1] + ap2 * m13c2[2] + ap3 * m13c2[3] + ap4 * m13c2[4] + ap5 * m13c2[5];
                    m13[bandnum][im][id][ip] = x1 + wt * (x2 - x1);
                }
            }
        free(file);
        free(name);
        free(sdsname);
    }

    if (l1rec->mside > 1)
        return 0x0;

    switch (type) {
    case XRVS: return (&rvs[bandnum][l1rec->mside][l1rec->detnum / detfac][0]);
        break;
    case XM12: return (&m12[bandnum][l1rec->mside][l1rec->detnum / detfac][0]);
        break;
    case XM13: return (&m13[bandnum][l1rec->mside][l1rec->detnum / detfac][0]);
        break;
    default:
        printf("-E- %s line %d: bad xcal type specified\n", __FILE__, __LINE__);
        exit(1);
    }


}

double *get_fpm_xcal(char *fpm_file)
{
    FILE *fp;
    char *line = NULL;
    char *ptr = NULL;
    double *fpm_gains = NULL;
    size_t bufsz = 100;
    ssize_t llen;
    int nbands, nsca;
    int linenum, isca, ib, isb;

    if ((fp = fopen(fpm_file,"r")) == NULL)
    {
        printf("-E- %s line %d: Cannot open fpm xcal file %s\n",__FILE__,__LINE__, fpm_file);
        exit(1);
    }
    if ((line = calloc(bufsz,sizeof(char))) == NULL)
    {
        printf("-E- %s line %d: Cannot allocate memory for line buffer\n",__FILE__,__LINE__);
        exit(1);
    }

    linenum = 0;
    isca = 0;
    while((llen = getline(&line, &bufsz, fp)) != -1)
    {
        if (linenum == 0)
        {
            if (strstr(line,"FPM") == NULL)
            {
                printf("-E- %s line %d: Missing FPM numbers in the first line of fpm xcal file %s\n",__FILE__,__LINE__, fpm_file);
                exit(1);
            }
            else
            {
                sscanf(line,"%*[^=]= %d", &nsca);
                linenum++;
                continue;
            }
        }

        if (linenum == 1)
        {
            if (strstr(line,"BANDS") == NULL)
            {
                printf("-E- %s line %d: Missing FPM numbers in the first line of fpm xcal file %s\n",__FILE__,__LINE__, fpm_file);
                exit(1);
            }
            else
            {
                sscanf(line,"%*[^=]= %d", &nbands);
                fpm_gains = (double *)calloc(nbands, nsca*sizeof(double) );
                if ( fpm_gains == NULL )
                {
                    printf("-E- %s line %d: Cannot allocate memory for fpm gain buffer\n",__FILE__,__LINE__);
                    exit(1);
                }
                linenum++;
                continue;
            }
        }

        //printf("nbands = %d, nsca = %d\n", nbands, nsca);
        if ((ptr = strtok(line,",")) != NULL)
        {
            for (ib = 0; ib< nbands; ib++)
            {
                isb = isca*nbands + ib;
                fpm_gains[isb] = atof(ptr);
                ptr = strtok(NULL,",");
            }
        }
        linenum++;
        isca++;
    }

    /* for (isca = 0; isca < nsca; isca++)
        for (ib = 0; ib < nbands; ib++)
        {
            isb = isca*nbands + ib;
            printf("fpm_gains = %lf\n", fpm_gains[isb]);
        } */

    free( line );
    fclose(fp);

    return fpm_gains;
}
