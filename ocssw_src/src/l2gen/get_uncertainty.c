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
    static float Lr_unc=0.0;
    int polynomial_order=4, iorder;
    float tmp_poly;
    static float *noise_temp=NULL;
    static char model_type[20]="\0";
    static int nbands_table=0, *bandindex;
    static float  *wave_table;
    static int sensorID;

    float scaled_lt;
    static int uncertainty_lut=1;

    if(firstcall){
        firstcall=0;

        if(uncertainty)
            corr_coef_rhot=uncertainty->corr_coef_rhot;

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

            status = nc_inq_varid(ncid, "wave", &sds_id);
            if (status != NC_NOERR) {
                fprintf(stderr, "-E- %s line %d:  Error reading %s from %s.\n",
                        __FILE__, __LINE__, "wave", filename);
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
            if(nc_get_att_text(group_id,sds_id,"model_type",model_type)!=NC_NOERR){
                fprintf(stderr, "-E- %s line %d:  Error reading model_type attribute from %s.\n",
                        __FILE__, __LINE__,filename);
                exit(1);
            }
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
            for(iw=0;iw<nbands;iw++){
                tmp_poly=0.;
                ipb=ip*nbands+iw;
                scaled_lt=l1rec->Lt[ipb]*10.;

                iw_table=bandindex[iw];

                if(strstr(model_type,"polynomial"))
                    for(iorder=0;iorder<=polynomial_order;iorder++)
                        tmp_poly+=noise_coef[iw_table*(polynomial_order+1)+iorder]*pow(scaled_lt,iorder);

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

            for(ib=0;ib<nbands;ib++)
                corr_coef_rhot[iw*nbands+ib]=temp_corr_coef[iw_table*nbands_table+bandindex[ib]];

            if(sensorID==OCI){
                for(ib=0;ib<nbands;ib++){
                    corr_coef_rhot[iw*nbands+ib]=0.;
                    if(ib==iw)
                        corr_coef_rhot[iw*nbands+ib]=1.;
                }
            }
        }
    }

    // calculate uncertainty of uncertainty sources
    for (ip=0;ip<npix;ip++) {
        for(iw=0;iw<nbands;iw++){
            tmp_poly=0.;
            ipb=ip*nbands+iw;
            scaled_lt=l1rec->Lt[ipb]*10.;

            iw_table=bandindex[iw];
            if(strstr(model_type,"polynomial"))
                for(iorder=0;iorder<=polynomial_order;iorder++)
                    tmp_poly+=noise_coef[iw_table*(polynomial_order+1)+iorder]*pow(scaled_lt,iorder);

            uncertainty->dsensor[ipb]=tmp_poly/SNR_scale[iw_table]/10.;

            if(sensorID==OCI)
                uncertainty->dsensor[ipb]=sqrt(tmp_poly)/SNR_scale[iw_table]/10.;

            uncertainty->dvc[ipb]=rel_unc_vc[iw_table]*l1rec->Lt[ipb];
            uncertainty->dLr[ipb]=l1rec->Lr[ipb]*Lr_unc;

            noise_temp[ipb] =sqrt(uncertainty->dvc[ipb] * uncertainty->dvc[ipb] +uncertainty->dsensor[ipb]*uncertainty->dsensor[ipb]);
        }
    }

    /*  error in gas transmittance */

    gas_trans_uncertainty(l1rec);

    return noise_temp;
}


