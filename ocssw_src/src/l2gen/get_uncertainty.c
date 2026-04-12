//#include <stdio.h>
#include "l12_proto.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <sys/time.h>

#include <hdf5.h>
//--------------------------//

#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}

unsigned long int random_seed() {
    /* Seed generator for gsl. */
    struct timeval tv;
    gettimeofday(&tv, 0);
    return (tv.tv_sec + tv.tv_usec);
}

float make_noise(float sigma) {
    unsigned long randSeed = random_seed();
    float noise;
    gsl_rng *rng;
    rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng, randSeed);
    noise = gsl_ran_gaussian(rng, sigma);
    gsl_rng_free(rng);
    return (noise);
}

/* get the No. and  corresponding index of bands that used for calculating rhownir */
/* please note chl is used for calculating rhownir, so the bands used to calculate chla should also be included */
/*  it only works for chl_oci right now*/
void get_bands_rhownir(int *band_index, l1str* l1rec){
    
    uncertainty_t *uncertainty=l1rec->uncertainty;
    int *nbands=&uncertainty->nbands_rhownir;
    int num_oc=0,i,ib2,ib5,ib6;
    int nwave=l1rec->l1file->nbands;
    float *wave=l1rec->l1file->fwave;
    float *temp_index=(float *)malloc(nwave*sizeof(float));

    /* bands used for calculating rhownir*/
    ib6 = windex(670, wave, nwave);
    if (fabs(670 - wave[ib6]) > 50) {
        printf("%s line %d: can't find reasonable red band\n", __FILE__, __LINE__);
        printf("looking for 670, found %f\n", wave[ib6]);
        exit(EXIT_FAILURE);
    }

    ib5 = windex(555, wave, nwave);
    if (fabs(555 - wave[ib5]) > 15) {
        printf("%s line %d: can't find reasonable green band\n", __FILE__, __LINE__);
        printf("looking for 555, found %f\n", wave[ib5]);
        exit(EXIT_FAILURE);
    }
    ib2 = windex(443, wave, nwave);
    if (fabs(443 - wave[ib2]) > 5) {
        printf("%s line %d: can't find reasonable blue band\n", __FILE__, __LINE__);
        printf("looking for 443, found %f\n", wave[ib2]);
        exit(EXIT_FAILURE);
    }

    band_index[0]=ib2;
    band_index[1]=ib5;
    band_index[2]=ib6;

    /* bands used for calculating chla*/

     switch (l1rec->l1file->sensorID) {
    case SEAWIFS:
    case OCTS:
    case OCM1:
    case OCM2:
    case MERIS:
    case HICO:
    case OCI:
    case HAWKEYE:
    case AVIRIS:
    case OLCIS3A:
    case OLCIS3B:
        num_oc=4;
        for(i=0;i<num_oc;i++)
            band_index[3+i]=windex(input->chloc4w[i]*1.0,wave, nwave);
        break;
    case MODIST:
    case MODISA:
    case CZCS:
    case VIIRSN:
    case VIIRSJ1:
    case VIIRSJ2:
    case OCRVC:
    case GOCI:
    case OLIL8:
    case OLIL9:
    case PRISM:
    case SGLI:
    case MSIS2A:
    case MSIS2B:
         num_oc=3;
        for(i=0;i<num_oc;i++)
            band_index[3+i]=windex(input->chloc3w[i]*1.0,wave, nwave);
        break;
    case L5TM:
    case L7ETMP:
    case MISR:
         num_oc=2;
        for(i=0;i<num_oc;i++)
            band_index[3+i]=windex(input->chloc2w[i]*1.0,wave, nwave);
        break;
    default:
        printf("%s Line %d: need a default chlorophyll algorithm for this sensor\n",
                __FILE__,__LINE__);
        exit(1);
        break;
    }

    // put the index in order
    int j=0,ifexist=0,total=0;
    for(j=0;j<num_oc;j++){
        ifexist=0;
        for (i = 0; i < 3; i++) {
            if (band_index[3 + j] == band_index[i])
                ifexist = 1;
        }
        if(!ifexist){
            band_index[3+total]=band_index[3+j];
            total++;
        }    
    }

    *nbands=3+total;

    for (i = 0; i < *nbands - 1; i++)
        for (j = i + 1; j < *nbands; j++) {
            if (band_index[j] < band_index[i]) {
                total = band_index[i];
                band_index[i] = band_index[j];
                band_index[j] = total;
            }
        }

    free(temp_index);

}

/*
if uncertainty products are produced, calculate the uncertainty of the uncertainty sources, return NULL.
Otherwise, calculate the sensor noise, that is used by the mbac AC algorithm, return sensor_noise.

*/
float *get_uncertainty(l1str *l1rec) {

#define N_anc 7

    int ip,iw,ib,ipb,iw_table;
    int nbands=l1rec->l1file->nbands;
    uncertainty_t *uncertainty=l1rec->uncertainty;
    int npix=l1rec->npix;
    float *wave=l1rec->l1file->fwave;
    static int firstcall=1;
    static float *SNR_scale, *noise_coef, *corr_coef_rhot, *rel_unc_vc,*temp_corr_coef;
    int polynomial_order=4, iorder;
    float tmp_poly;
    static float *noise_temp=NULL;
    //static char model_type[20]="\0";
    static int nbands_table=0, *bandindex;
    static float  *wave_table;
    static int sensorID;
    static int *bindx_rhownir, nbands_rhownir;

    float scaled_lt;
    static int uncertainty_lut=1;

    if(firstcall){
        firstcall=0;
        bindx_rhownir=(int *)malloc(nbands*sizeof(int));

        if(uncertainty){
            corr_coef_rhot=uncertainty->corr_coef_rhot;
            get_bands_rhownir(bindx_rhownir,l1rec);
            nbands_rhownir=uncertainty->nbands_rhownir;
        }
            

        if(input->uncertaintyfile[0]=='\0'){
            uncertainty_lut=0;
            return NULL;
        }

        sensorID=l1rec->l1file->sensorID;

        bandindex=(int *)malloc(nbands*sizeof(int));

        for(iw=0;iw<nbands;iw++)
            bandindex[iw]=iw;

        noise_temp = (float *)malloc(nbands * npix * sizeof(float));

        char filename[FILENAME_MAX];
        sprintf(filename, "%s",input->uncertaintyfile);
        printf("Reading uncertainty from: %s\n", filename);

        int ncid;
        int32 sds_id,group_id;
        nc_type rh_type; /* variable type */
        int rh_dimids[H4_MAX_VAR_DIMS]; /* dimension IDs */
        int rh_natts; /* number of attributes */
        int rh_ndims; /* number of dims */
        int status;

        if (nc_open(filename, NC_NOWRITE, &ncid) == NC_NOERR) {

            status = nc_inq_varid(ncid, "wavelength", &sds_id);
            if (status != NC_NOERR) {
                fprintf(stderr, "-E- %s line %d:  Error reading %s from %s.\n",
                        __FILE__, __LINE__, "wavelength", filename);
                exit(1);
            }
            status = nc_inq_var(ncid, sds_id, 0, &rh_type, &rh_ndims, rh_dimids,
                    &rh_natts);

            size_t tmpSize;
            DPTB( nc_inq_dimlen( ncid, rh_dimids[0], &tmpSize ) );
            nbands_table = tmpSize;

            wave_table=(float *)malloc(nbands_table*sizeof(float));
            SNR_scale=(float *)malloc(nbands_table*sizeof(float));
            noise_coef=(float *)malloc(nbands_table*5*sizeof(float));

            temp_corr_coef=(float *)malloc(nbands_table*nbands_table*sizeof(float));
            rel_unc_vc= (float *)malloc (nbands_table*sizeof(float));

            if (nc_get_var(ncid, sds_id, wave_table) != NC_NOERR) {
                fprintf(stderr, "-E- %s line %d:  Error reading %s from %s.\n",
                        __FILE__, __LINE__, "wave", filename);
                exit(1);
            }
            if(nbands_table!=nbands){
                printf("The band No. in uncertainty.nc is not equal to the band No. of sensor,interpolation is used\n");

                for(iw=0;iw<nbands;iw++) {

                    for(iorder=0;iorder<nbands_table;iorder++){
                        if(wave[iw]<=wave_table[iorder])
                            break;
                    }
                    if(iorder==0)
                        iw_table=0;
                    else if(iorder==nbands_table)
                        iw_table=nbands_table-1;
                    else{
                        iw_table=iorder;
                        if(fabs(wave[iw]-wave_table[iorder])>fabs(wave[iw]-wave_table[iorder-1]))
                            iw_table=iorder-1;
                    }
                    bandindex[iw]=iw_table;
                }
            }

            status = nc_inq_ncid(ncid,"random_noise",&group_id);
            status = nc_inq_varid(group_id, "SNR_scale", &sds_id);
            if (nc_get_var(group_id, sds_id, SNR_scale) != NC_NOERR) {
                fprintf(stderr, "-E- %s line %d:  Error reading %s from %s.\n",
                        __FILE__, __LINE__, "SNR_scale", filename);
                exit(1);
            }
            status = nc_inq_varid(group_id, "sensor_noise", &sds_id);
           /* if(nc_get_att_text(group_id,sds_id,"model_type",model_type)!=NC_NOERR){
                fprintf(stderr, "-E- %s line %d:  Error reading model_type attribute from %s.\n",
                        __FILE__, __LINE__,filename);
                exit(1);
            }*/
            if (nc_get_var(group_id, sds_id, noise_coef) != NC_NOERR) {
                fprintf(stderr, "-E- %s line %d:  Error reading %s from %s.\n",
                        __FILE__, __LINE__, "sensor_noise", filename);
                exit(1);
            }

            status = nc_inq_ncid(ncid,"vicarious_cal",&group_id);
            status = nc_inq_varid(group_id, "correlation_coefficient", &sds_id);
            if (nc_get_var(group_id, sds_id, temp_corr_coef) != NC_NOERR) {
                fprintf(stderr, "-E- %s line %d:  Error reading %s from %s.\n",
                        __FILE__, __LINE__, "correlation_coefficient", filename);
                exit(1);
            }
            status = nc_inq_varid(group_id, "relative_uncertainty", &sds_id);
            if (nc_get_var(group_id, sds_id, rel_unc_vc) != NC_NOERR) {
                fprintf(stderr, "-E- %s line %d:  Error reading %s from %s.\n",
                        __FILE__, __LINE__, "relative_uncertainty", filename);
                exit(1);
            }

            if (nc_close(ncid) != NC_NOERR) {
                fprintf(stderr, "-E- %s line %d: error closing %s.\n",
                        __FILE__, __LINE__, filename);
                exit(1);
            }
        }

    }

    if(!uncertainty){
        if(!uncertainty_lut)
            return NULL;
        for (ip=0;ip<npix;ip++)
            for(iw=0;iw<nbands;iw++) {
                tmp_poly=0.;
                ipb=ip*nbands+iw;
                scaled_lt=l1rec->Lt[ipb]*10.;

                iw_table=bandindex[iw];

                //if (strstr(model_type, "polynomial")) {
                    float pow_scaled = 1;
                    for (iorder = 0; iorder <= polynomial_order; iorder++) {
                        tmp_poly += noise_coef[iw_table * (polynomial_order + 1) + iorder] * pow_scaled;
                        pow_scaled *= scaled_lt;
                    }
              //  }

                noise_temp[ipb]=tmp_poly/SNR_scale[iw_table]/10.;
                if(sensorID==OCI)
                    noise_temp[ipb]=sqrt(tmp_poly)/SNR_scale[iw_table]/10.;

               noise_temp[ipb] =
                    sqrt(l1rec->Lt[ipb] * rel_unc_vc[iw_table] * l1rec->Lt[ipb] * rel_unc_vc[iw_table] +
                         noise_temp[ipb] * noise_temp[ipb]);

            }
        return noise_temp;
    }

    if(uncertainty){
        for(iw=0;iw<nbands;iw++){
            iw_table=bandindex[iw];

            for(ib=0;ib<nbands;ib++){
                corr_coef_rhot[iw * nbands + ib] = temp_corr_coef[iw_table * nbands_table + bandindex[ib]];
                if (ib < iw)
                    corr_coef_rhot[iw * nbands + ib] =
                        temp_corr_coef[bandindex[ib] * nbands_table +iw_table];
            }
        }

        for (iw = 0; iw < nbands_rhownir; iw++)
            uncertainty->bindx_rhownir[iw] = bindx_rhownir[iw];
    }

    // calculate uncertainty of uncertainty sources
    for (ip=0;ip<npix;ip++) {
        for(iw=0;iw<nbands;iw++){
            tmp_poly=0.;
            ipb=ip*nbands+iw;
            scaled_lt=l1rec->Lt[ipb]*10.;

            iw_table=bandindex[iw];
           // if (strstr(model_type, "polynomial")) {
                float pow_scaled = 1;
                for (iorder = 0; iorder <= polynomial_order; iorder++) {
                    tmp_poly += noise_coef[iw_table * (polynomial_order + 1) + iorder] * pow_scaled;
                    pow_scaled *= scaled_lt;
                }
           // }

            uncertainty->dsensor[ipb]=tmp_poly/SNR_scale[iw_table]/10.;

            if(sensorID==OCI)
                uncertainty->dsensor[ipb]=sqrt(tmp_poly)/SNR_scale[iw_table]/10.;

            uncertainty->dvc[ipb]=rel_unc_vc[iw_table]*l1rec->Lt[ipb];

            noise_temp[ipb] =sqrt(uncertainty->dvc[ipb] * uncertainty->dvc[ipb] +uncertainty->dsensor[ipb]*uncertainty->dsensor[ipb]);
        }
    }

    /*  error in gas transmittance */

    gas_trans_uncertainty(l1rec);

    return noise_temp;
}


