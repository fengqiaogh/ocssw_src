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

void noise_model_ocis(l1str *l1rec, float *noise) {

	static int firstcall=1;
	static size_t nwave,ncoef;
    static float *Lnoise_c;
    static int16 *wave;

	int nbands=l1rec->l1file->nbands;
	int npix=l1rec->npix;
	float *wave_oci=l1rec->l1file->fwave;
    static int *wvindx;
	int i,j;

	if(firstcall){
		firstcall=0;

		wvindx=(int *)malloc(nbands*sizeof(int));

        char *filedir;
        char filename[FILENAME_MAX];
        if ((filedir = getenv("OCDATAROOT")) == NULL) {
            printf("-E- %s: OCDATAROOT env variable undefined.\n", __FILE__);
            exit(EXIT_FAILURE);
        }
        strcpy(filename, filedir);
        strcat(filename, "/");
        strcat(filename, "ocis");
        strcat(filename, "/");
        strcat(filename, "oci_snr.h5");
        printf("OCI SNR FILE: %s\n", filename);

        int ncid, varid1, varid2;
        /* error handling. */
        int retval;

        /* Open the file. NC_NOWRITE tells netCDF we want read-only access to the file.*/
        if (retval = nc_open(filename, NC_NOWRITE, &ncid))
            ERR(retval);

        /* Get the varid of the data variable, based on its name. */
		 if ((retval = nc_inq_varid(ncid, "LNOISE_C", &varid1)) )
            ERR(retval);

        int dimids[2];
		 if(retval = nc_inq_vardimid(ncid, varid1, dimids))
            ERR(retval);

		 if(retval = nc_inq_dimlen(ncid, dimids[1], &nwave))
            ERR(retval);
		 if(retval = nc_inq_dimlen(ncid, dimids[0], &ncoef))
            ERR(retval);

		 Lnoise_c=(float *)malloc(ncoef*nwave*sizeof(float));

        /* Read the data. */
		 if ((retval = nc_get_var_float(ncid, varid1, Lnoise_c)) )
            ERR(retval);

		 wave=(int16_t *)malloc(nwave*sizeof(int16_t));
		 if ((retval = nc_inq_varid(ncid, "WAVE", &varid2)) )
            ERR(retval);
		 if ((retval = nc_get_var_short(ncid, varid2, wave)) )
            ERR(retval);

        /* Close the file, freeing all resources. */
        if ((retval = nc_close(ncid)))
            ERR(retval);

        float wavedif;
		 for (i=0;i<nbands;i++) {
			 wavedif=fabs(wave_oci[i]-(float)wave[0]);
			 wvindx[i]=0;
			 for(j=1;j<nwave;j++){
				 if(fabs(wave_oci[i]-(float)wave[j])<wavedif){
					 wavedif=fabs(wave_oci[i]-(float)wave[j]);
					 wvindx[i]=j;
                }
            }
        }
        free(wave);
    }

    float scaled_lt;
    for (i=0;i<npix;i++) {
    	for(j=0;j<nbands;j++){
    		scaled_lt=l1rec->Lt[i*nbands+j]*10;  //covert from mw.cm-2.sr-1.um-1 to w.m-2.sr-1.um-1.
    		noise[i*nbands+j]=sqrt(Lnoise_c[wvindx[j]]+Lnoise_c[nwave+wvindx[j]]*scaled_lt );// /snr_scale[wvindx[j]];
        }
    }
}

///derive the index of OCI wavelength in MODIS wavelength domain
int noise_index_oci(float wave)
{
    int i, index;
	float wave_modis[16]={412,443,469,488,531,547,555,645,667,678,748,859,869,1240,1640,2130};

	if(wave<wave_modis[0])
		index=0;
	else if(wave>wave_modis[15])
		index=15;
	else{
		for(i=0;i<16;i++){
			if(wave<wave_modis[i])
                break;
        }
		index=i;
		if(fabs(wave-wave_modis[i-1])< fabs(wave-wave_modis[i]))
			index=i-1;
    }
    return index;
}

void set_error_input(l1str *l1rec, float *sensor_noise)
{
    filehandle *l1file=l1rec->l1file;
    uncertainty_t *uncertainty=l1rec->uncertainty;
    int32_t sensorID=l1file->sensorID;
    int32_t nbands=l1file->nbands;
    float *wave=l1file->fwave;
    int32_t npix=l1rec->npix;
    int32_t ip,ib,ipb;
    float mu,mu0;
    //float *pr=l1rec->pr;
    float *Tau_r=l1file->Tau_r;
    static float p0 = STDPR;
    //int ib869=windex(869,l1file->fwave,nbands);
    //float *vcgain=l1_input->gain;
    //float *offset=l1_input->offset;
    static float *dvc_rel;
    static float *noise;
    static int firstcall=1;

    //read the relative uncertainty coefficients in VC from a file?
    float dvc_rel_temp[16]={0.0077, 0.0083,1.0,0.0093,0.0105,0.0116,0.0121,1.0,0.0236,0.0249,0.035,1.0,0.05,1.0,1.0,1.0};//new version 869nm=5%, Moby 5%

    //                     412    443   469  488     531    547    555    645   667    678    748    859     869  1240 1640 2130
    //float corr_nir_s[16]={0.2851,0.3336,0.0,0.3993,0.4906,0.5589,0.5870,0.0000,0.9443,0.9585,1.0000,0.0000,0.9724,1.0,1.0,1.0};
    //float corr_nir_l[16]={0.2172,0.2116,0.0,0.2177,0.2973,0.3710,0.4002,0.0000,0.8534,0.8784,0.9724,0.0000,1.0000,1.0,1.0,1.0};

    //                     412    443   469  488     531    547    555  645  667    678   748 859 869 1240 1640 2130
    //float corr_nir_s[16]={0.2851,0.3336,0.0,0.3993,0.4906,0.5589,0.5870,0.0,1.000,1.000,1.0,0.0,1.0,0.0,0.0,0.0};
    //float corr_nir_l[16]={0.2172,0.2116,0.0,0.2177,0.2973,0.3710,0.4002,0.0,1.000,1.000,1.0,0.0,1.0,0.0,0.0,0.0};


    if(firstcall){
        firstcall=0;
        dvc_rel=(float *)malloc(nbands*sizeof(float));
        noise=(float *)malloc(npix*nbands*sizeof(float));

        for(ib=0;ib<nbands;ib++)
            dvc_rel[ib]=0.;

        if(uncertainty){
            switch(sensorID){
                case OCI:
                case OCIS:
                for(ib=0;ib<nbands;ib++){

                    if(wave[ib]<400)
                        dvc_rel[ib]=0.0047;
                    else if(wave[ib]<=865)
                        dvc_rel[ib]=0.004;

                    }
                ib=bindex_get(820.0);
                dvc_rel[ib]=0.0075;

                //for glimr
                for(ib=0;ib<nbands;ib++){
                        // uncertainty->corr_nir_s[ib]=0.59;
                    //uncertainty->corr_nir_l[ib]=0.53;
                        // dvc_rel[ib]=0.005;
                    }
                //uncertainty->corr_nir_s[222]=0.97;
                    // uncertainty->corr_nir_l[176]=0.97;
                    break;
                case MODISA:
                //float dvc_rel[16]={0.0074, 0.0078,1.0,0.0082,0.0078,0.0082,0.0084,1.0,0.0144,0.0151,0.021,1.0,0.03,1.0,1.0,1.0}; //new version 869nm=3%, Moby 5%
                //float dvc_rel_temp[16]={0.0077, 0.0083,1.0,0.0093,0.0105,0.0116,0.0121,1.0,0.0236,0.0249,0.035,1.0,0.05,1.0,1.0,1.0};//new version 869nm=5%, Moby 5%
                //float dvc_rel[16]={0.0056, 0.0063,1.0,0.0075,0.0101,0.0114,0.0118,1.0,0.0236,0.0249,0.035,1.0,0.05,1.0,1.0,1.0};//new version 869nm=5%, Moby 3%
                for(ib=0;ib<nbands;ib++){
                    dvc_rel[ib]=dvc_rel_temp[ib];
                    }
                    break;
                default:
                    break;
            }
        }


    }

    switch (sensorID) {
        case OCI:
        case OCIS:
        noise_model_ocis(l1rec,noise);
            break;

        case SEAWIFS:
        case MERIS:
        case OCTS:
        case OCM1:
        case OCM2:
        case MOS:
        case HICO:
        case CZCS:
        case OSMI:
        case VIIRSN:
        case VIIRSJ1:
        case VIIRSJ2:
        case OCRVC:
        case GOCI:
        printf("%s Line %d: need a noise mode for this sensor\n",
                __FILE__, __LINE__);
            exit(1);
            break;
        default:
        printf("%s Line %d: need a noise mode for this sensor\n",
                __FILE__, __LINE__);
            exit(1);
            break;
    }

     //lt_agregat_ocis(l1rec);

    if(!uncertainty){
        for(ip=0;ip<npix; ip++)
            for(ib=0;ib<nbands;ib++){
                ipb=ip*nbands+ib;
                sensor_noise[ipb]=noise[ipb]/10.;
            }
        return;
    }

    for(ip=0;ip<npix; ip++){

        for(ib=0;ib<nbands;ib++){
            ipb=ip*nbands+ib;

            //indx=noise_index_oci(wave[ib]);

            // random noise for sensor
            //noise[ipb]=0.;
            uncertainty->dsensor[ipb]=noise[ipb]/10.;   //covert w.m-2.sr-1.um-1 to mw.cm-2.sr-1.um-1 .

            //uncertainty->dsensor[ipb]*=vcgain[ib];

            uncertainty->dvc[ipb]=dvc_rel[ib]*l1rec->Lt[ipb];
            //uncertainty->dvc[ipb]=0.0;

            //if(ib!=ib869)
            //  uncertainty->dsensor[ipb]+=0.005*l1rec->Lt[ipb];

            // Rayleigh systematic error based on Wang et al.(2005)
            // uncertainty->dLr[ipb]=l1rec->Lr[ipb]*0.001;
            //uncertainty->dLr[ipb]=0.0;
        }
    }

    // error in ancillary data
    float rel_error=0.0;
    for(ip=0;ip<npix; ip++)
    {
        uncertainty->drh[ip]=rel_error*l1rec->rh[ip];
        uncertainty->doz[ip]=rel_error*l1rec->oz[ip];
        uncertainty->dpr[ip]=rel_error*l1rec->pr[ip];
        uncertainty->dwv[ip]=rel_error*l1rec->wv[ip];
        uncertainty->dws[ip]=rel_error*l1rec->ws[ip];

        uncertainty->dno2_tropo[ip]=rel_error*l1rec->no2_tropo[ip];
        uncertainty->dno2_strat[ip]=rel_error*l1rec->no2_strat[ip];
    }

    //error in diffuse transmittance resulted from the ancillary data
    for(ip=0;ip<npix;ip++)
    {
        mu0=l1rec->csolz[ip];
        mu= l1rec->csenz[ip];
        for(ib=0;ib<nbands;ib++)
        {
            ipb=ip*nbands+ib;
            uncertainty->dt_sol[ipb]=l1rec->t_sol [ipb]*0.5 / p0 * Tau_r[ib] / mu0*uncertainty->dpr[ip];
            uncertainty->dt_sen[ipb]=l1rec->t_sen [ipb]*0.5 / p0 * Tau_r[ib] / mu*uncertainty->dpr[ip];
        }
    }

    /*  error in gas transmittance */

    gas_trans_uncertainty(l1rec);

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

        /* if(l1rec->l1file->sensorID==OCI ||l1rec->l1file->sensorID==OCIS  ){
             if(firstcall){
                 firstcall=0;

                 if(!uncertainty)
                     noise_temp=(float *)malloc(nbands*npix*sizeof(float));
             }
             set_error_input(l1rec,noise_temp);
             return noise_temp;
         }*/


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
                if(sensorID==OCI || sensorID==OCIS)
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

            if(sensorID==OCI || sensorID==OCIS){
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

            if(sensorID==OCI || sensorID==OCIS)
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


