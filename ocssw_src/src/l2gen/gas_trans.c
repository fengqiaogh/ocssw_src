/* =========================================================== */
/* Module gaseous_transmittance.c                              */
/*                                                             */
/* Computes sensor-specific transmittance through various      */
/* atmospheric gases.                                          */
/*                                                             */
/* B. Franz, NASA/OBPG, July 2006                              */
/*                                                             */
/* adding the capability to read LUT with air mass factor      */
/* as one additional dimension                                 */
/* M. Zhang, SAIC/NASA, Nov. 2024                              */
/* =========================================================== */

#include "l12_proto.h"
#include <allocate3d.h>
#include "atrem_corl1.h"

float get_wv_band_ratio(l1str *l1rec,int32_t ip,float window1, float absorp_band,float window2);
float get_airmass_oxygen(l1str *l1rec,int32_t ip,float window1, float absorp_band,float window2);
int32_t get_index_lowerbound(float *table_val, int32_t num_val, float val);
int32_t get_index_uppperbound(float *table_val, int32_t num_val, float val);

static int32_t model = 5;
static int32_t amf;
static size_t num_models;
static size_t num_wavelengths;
static size_t num_water_vapors;
static size_t num_airmass;
static float *wvtbl = NULL;
static float *t_co2 = NULL;
static float *t_o2 = NULL;
static float *t_co = NULL;
static float *t_ch4 = NULL;
static float *t_n2o = NULL;
static float *cwv_all = NULL;

static float *amf_mixed = NULL;
static float *amf_wv = NULL;

static int32_t index_amf_solz;
static int32_t index_amf_total;
//interpolation ratio for amf at direction of solz and two-way
static float ratio_solz;
static float ratio_total;
static float amf_solz;
static float amf_senz;
static float amf_total;

void load_gas_tables(l1str *l1rec) {
    /*
        netcdf file that has water vapor transmittance as a function of wavelength and
        cwv (6 profiles x number of bands x 220 water vapor value)
        filename = /OCDATAROOT/sensor[/subsensor]/<sensorName>_gas_trans.nc
    */

    /* This will be the netCDF ID for the file and data variable. */
    int32_t ncid, varid;
    int32_t num_water_vapors_id, num_models_id, num_wavelengths_id, num_airmass_id;
    amf = 0;
    num_airmass=1;

    /* Open the file */
    if ((nc_open(input->gas_transmittance_file, NC_NOWRITE, &ncid)) != NC_NOERR) {
        printf("-E- %s: Failed to open %s\n", __FILE__, input->gas_transmittance_file);
        exit(EXIT_FAILURE);
    } else {
        printf("Reading gas transmittance file: %s.\n", input->gas_transmittance_file);
    }
    if ((nc_inq_dimid(ncid, "n_air_mass_factor", &num_airmass_id)) == NC_NOERR) {
        amf = 1;
        if ((nc_inq_dimlen(ncid, num_airmass_id, &num_airmass)) != NC_NOERR) {
            printf("-E- %s: Failed to read dimension n_water_vapor\n", __FILE__);
            exit(EXIT_FAILURE);
        }
    }
    if ((nc_inq_dimid(ncid, "n_water_vapor", &num_water_vapors_id)) == NC_NOERR) {
        if ((nc_inq_dimlen(ncid, num_water_vapors_id, &num_water_vapors)) != NC_NOERR) {
            printf("-E- %s: Failed to read dimension n_water_vapor\n", __FILE__);
            exit(EXIT_FAILURE);
        }
    } else {
        printf("-E- %s: Failed to find dimension n_water_vapor\n", __FILE__);
        exit(EXIT_FAILURE);
    }
    if((nc_inq_dimid(ncid, "nmodels", &num_models_id)) == NC_NOERR){
        if((nc_inq_dimlen(ncid, num_models_id, &num_models)) != NC_NOERR){
            printf("-E- %s: Failed to read dimension nmodels\n", __FILE__);
            exit(EXIT_FAILURE);
        }
    } else {
        printf("-E- %s: Failed to find dimension nmodels\n", __FILE__);
        exit(EXIT_FAILURE);
    }

    if((nc_inq_dimid(ncid, "nwavelengths", &num_wavelengths_id)) == NC_NOERR){
        if((nc_inq_dimlen(ncid, num_wavelengths_id, &num_wavelengths)) != NC_NOERR){
            printf("-E- %s: Failed to read dimension nwavelengths\n", __FILE__);
            exit(EXIT_FAILURE);
        }
    } else {
        printf("-E- %s: Failed to find dimension nwavelengths\n", __FILE__);
        exit(EXIT_FAILURE);
    }

    /* Read the water vapor transmittance */
    if ((wvtbl = (float *)malloc(num_models*num_wavelengths*num_airmass*num_water_vapors*sizeof(float))) == NULL) {
        printf("-E- %s line %d : error allocating memory for water vapor transmittance table.\n",
               __FILE__, __LINE__);
        exit(1);
    }
    if ((nc_inq_varid(ncid, "water_vapor_transmittance", &varid)) == NC_NOERR) {
        if ((nc_get_var_float(ncid, varid,wvtbl)) != NC_NOERR){
            printf("-E- %s: failed to read water_vapor_transmittance from %s\n", __FILE__,
                   input->gas_transmittance_file);
            exit(EXIT_FAILURE);
        }
    } else {
        printf("-E- %s: '%s' does not have water_vapor_transmittance.\n",
                    __FILE__, input->gas_transmittance_file);
        exit(EXIT_FAILURE);
    }

    if (amf) {
        if ((amf_mixed = (float *)malloc(num_airmass * sizeof(float))) == NULL) {
            printf("Error: allocating memory for air mass factor mixed\n");
            exit(EXIT_FAILURE);
        }
        if ((nc_inq_varid(ncid, "air_mass_factor_mixed", &varid)) == NC_NOERR) {
            if ((nc_get_var_float(ncid, varid, amf_mixed)) != NC_NOERR) {
                printf("-E- %s: failed to read air mass factor mixed from %s\n", __FILE__,
                       input->gas_transmittance_file);
                exit(EXIT_FAILURE);
            }
        } else {
            printf("-E- %s: '%s' does not have air mass factor mixed\n", __FILE__,
                   input->gas_transmittance_file);
            exit(EXIT_FAILURE);
        }

        if ((amf_wv = (float *)malloc(num_airmass * sizeof(float))) == NULL) {
            printf("Error: allocating memory for air mass factor for water vapor\n");
            exit(EXIT_FAILURE);
        }
        if ((nc_inq_varid(ncid, "air_mass_factor_wv", &varid)) == NC_NOERR) {
            if ((nc_get_var_float(ncid, varid, amf_wv)) != NC_NOERR) {
                printf("-E- %s: failed to read air mass factor for water vapor from %s\n", __FILE__,
                       input->gas_transmittance_file);
                exit(EXIT_FAILURE);
            }
        } else {
            printf("-E- %s: '%s' does not have air mass factor for water vapor\n", __FILE__,
                   input->gas_transmittance_file);
            exit(EXIT_FAILURE);
        }
    }

    /* Read the water vapor table */
    if ((cwv_all = (float *)malloc(num_water_vapors*sizeof(float))) == NULL) {
        printf("Error: allocating memory for water vapor table\n");
        exit(EXIT_FAILURE);
    }
    if ((nc_inq_varid(ncid, "water_vapor", &varid)) == NC_NOERR) {
        if ((nc_get_var_float(ncid, varid, cwv_all)) != NC_NOERR) {
            printf("-E- %s: failed to read water_vapor from %s\n", __FILE__, input->gas_transmittance_file);
            exit(EXIT_FAILURE);
        }
    } else {
        printf("-E- %s: '%s' does not have water_vapor\n",
                    __FILE__, input->gas_transmittance_file);
        exit(EXIT_FAILURE);
    }

    /* Read the carbon monoxide transmittance */
    if ((t_co = (float *)malloc(num_wavelengths*num_airmass*sizeof(float))) == NULL) {
        printf("-E- %s line %d : error allocating memory for carbon monoxide transmittance table.\n",
               __FILE__, __LINE__);
        exit(1);
    }    
    if ((nc_inq_varid(ncid, "carbon_monoxide_transmittance", &varid)) == NC_NOERR) {
        if ((nc_get_var_float(ncid, varid, t_co)) != NC_NOERR) {
            printf("-E- %s: failed to read carbon_monoxide_transmittance from %s\n", __FILE__,
                   input->gas_transmittance_file);
            exit(EXIT_FAILURE);
        }
    } else {
        printf("-E- %s: '%s' does not have carbon_monoxide_transmittance\n",
                    __FILE__, input->gas_transmittance_file);
        exit(EXIT_FAILURE);
    }

    /* Read the carbon dioxide transmittance */
    if ((t_co2 = (float *)malloc(num_wavelengths*num_airmass*sizeof(float))) == NULL) {
        printf("-E- %s line %d : error allocating memory for carbon dioxide transmittance table.\n", __FILE__,
               __LINE__);
        exit(1);
    }
    if ((nc_inq_varid(ncid, "carbon_dioxide_transmittance", &varid)) == NC_NOERR) {
        if ((nc_get_var_float(ncid, varid, t_co2)) != NC_NOERR) {
            printf("-E- %s: failed to read carbon_dioxide_transmittance from %s\n", __FILE__,
                   input->gas_transmittance_file);
            exit(EXIT_FAILURE);
        }
    } else {
        printf("-E- %s: '%s' does not have carbon_dioxide_transmittance\n",
                    __FILE__, input->gas_transmittance_file);
        exit(EXIT_FAILURE);
    }

    /* Read the oxygen transmittance */
    if ((t_o2 = (float *)malloc(num_wavelengths*num_airmass*sizeof(float))) == NULL) {
        printf("-E- %s line %d : error allocating memory for oxygen transmittance table.\n", __FILE__,
               __LINE__);
        exit(1);
    }
    if ((nc_inq_varid(ncid, "oxygen_transmittance", &varid)) == NC_NOERR) {
        if ((nc_get_var_float(ncid, varid, t_o2)) != NC_NOERR) {
            printf("-E- %s: failed to read oxygen_transmittance from %s\n", __FILE__,
                   input->gas_transmittance_file);
            exit(EXIT_FAILURE);
        }
    } else {
        printf("-E- %s: '%s' does not have oxygen_transmittance\n",
                    __FILE__, input->gas_transmittance_file);
        exit(EXIT_FAILURE);
    }

    /* Read the nitrous oxide transmittance */
    if ((t_n2o = (float *)malloc(num_wavelengths*num_airmass*sizeof(float))) == NULL) {
        printf("-E- %s line %d : error allocating memory for nitrous oxide transmittance table.\n", __FILE__,
               __LINE__);
        exit(1);
    }
    if ((nc_inq_varid(ncid, "nitrous_oxide_transmittance", &varid)) == NC_NOERR) {
        if ((nc_get_var_float(ncid, varid, t_n2o)) != NC_NOERR) {
            printf("-E- %s: failed to read nitrous_oxide_transmittance from %s\n", __FILE__,
                   input->gas_transmittance_file);
            exit(EXIT_FAILURE);
        }
    } else {
        printf("-E- %s: '%s' does not have nitrous_oxide_transmittance\n",
                    __FILE__, input->gas_transmittance_file);
        exit(EXIT_FAILURE);
    }

    /* Read the methane transmittance */
    if ((t_ch4 = (float *)malloc(num_wavelengths*num_airmass*sizeof(float))) == NULL) {
        printf("-E- %s line %d : error allocating memory for methane transmittance table.\n", __FILE__,
               __LINE__);
        exit(1);
    }
    if ((nc_inq_varid(ncid, "methane_transmittance", &varid)) == NC_NOERR) {
        if ((nc_get_var_float(ncid, varid, t_ch4)) != NC_NOERR) {
            printf("-E- %s: failed to read methane_transmittance from %s\n", __FILE__,
                   input->gas_transmittance_file);
            exit(EXIT_FAILURE);
        }
    } else {
        printf("-E- %s: '%s' does not have methane_transmittance\n",
                    __FILE__, input->gas_transmittance_file);
        exit(EXIT_FAILURE);
    }

    /* Read the Nitrogen dioxide */
    if ((nc_inq_varid(ncid, "k_no2", &varid)) == NC_NOERR) {
        if ((nc_get_var_float(ncid, varid, l1rec->l1file->k_no2)) != NC_NOERR) {
            printf("-E- %s: failed to read k_no2 from %s\n", __FILE__,
                   input->gas_transmittance_file);
            exit(EXIT_FAILURE);
        }
    } 
    /* Read the ozone cross-section  */
    if ((nc_inq_varid(ncid, "k_oz", &varid)) == NC_NOERR) {
        if ((nc_get_var_float(ncid, varid, l1rec->l1file->k_oz)) != NC_NOERR) {
            printf("-E- %s: failed to read k_oz from %s\n", __FILE__,
                   input->gas_transmittance_file);
            exit(EXIT_FAILURE);
        }
    } 

    /* Close the file */
    if ((nc_close(ncid)) != NC_NOERR){
        printf("-E- %s: failed to close %s\n",
            __FILE__, input->gas_transmittance_file);
        exit(EXIT_FAILURE);
    }
}

int32_t get_index_lowerbound(float *table_val, int32_t num_val, float val) {
    int32_t index;
    int32_t i;

    for (i = 0; i < num_val; i++)
        if (val < table_val[i])
            break;
    index=MAX(i-1,0);
    index=MIN(index,num_val-2);
    return index;
}

int32_t get_index_upperbound(float *table_val, int32_t num_val, float val) {
    int32_t index;

    for (index = 0; index < num_val; index++)
        if (val >= table_val[index])
            break;

    return index;
}


void ozone_transmittance(l1str *l1rec, int32_t ip) {
    float tau_oz;
    int32_t iw;

    filehandle *l1file = l1rec->l1file;
    int32_t nwave = l1file->nbands;
    int32_t ipb = ip*nwave;

    for (iw = 0; iw < nwave; iw++) {
        tau_oz = l1rec->oz[ip] * l1file->k_oz[iw];
        l1rec->tg_sol[ipb + iw] *= exp(-(tau_oz / l1rec->csolz[ip]));
        if(amf) {
            l1rec->tg[ipb + iw]*=exp(-tau_oz* (1.0/ l1rec->csolz[ip]+1.0/ l1rec->csenz[ip]));
        } else {
            l1rec->tg_sen[ipb + iw] *= exp(-(tau_oz / l1rec->csenz[ip]));
        }

    }
}

void co2_transmittance(l1str *l1rec, int32_t ip) {
    int32_t iw;

    filehandle *l1file = l1rec->l1file;
    int32_t nwave = l1file->nbands;
    int32_t ipb = ip*nwave;

    if ((input->gas_opt & GAS_TRANS_TBL_BIT) == 0) {
        if (t_co2 == NULL) {
        rdsensorinfo(l1file->sensorID, l1_input->evalmask, "t_co2", (void **) &t_co2);
        }
    }

    for (iw = 0; iw < nwave; iw++) {
        
        if(amf){
            int32_t index=iw*num_airmass;
            float t_co2_interp;

            t_co2_interp=t_co2[index+index_amf_solz]*(1-ratio_solz)+t_co2[index+index_amf_solz+1]*ratio_solz;
            l1rec->tg_sol[ipb + iw] *= t_co2_interp;

            t_co2_interp=t_co2[index+index_amf_total]*(1-ratio_total)+t_co2[index+index_amf_total+1]*ratio_total;
            l1rec->tg[ipb + iw] *= t_co2_interp;
        } else {
            l1rec->tg_sol[ipb + iw] *= pow(t_co2[iw], amf_solz);
            l1rec->tg_sen[ipb + iw] *= pow(t_co2[iw], amf_senz);
        }
    }

}

void co_transmittance(l1str *l1rec, int32_t ip) {
    int32_t iw;

    filehandle *l1file = l1rec->l1file;
    int32_t nwave = l1file->nbands;
    int32_t ipb = ip*nwave;

    /* Only compute if gas transmittance table was requested */
    if ((input->gas_opt & GAS_TRANS_TBL_BIT) != 0) {
        for (iw = 0; iw < nwave; iw++) {
            if (amf) {
                int32_t index=iw*num_airmass;
                
                float t_co_interp=t_co[index+index_amf_solz]*(1-ratio_solz)+t_co[index+index_amf_solz+1]*ratio_solz;
                l1rec->tg_sol[ipb + iw] *= t_co_interp;
                
                t_co_interp=t_co[index+index_amf_total]*(1-ratio_total)+t_co[index+index_amf_total+1]*ratio_total;
                l1rec->tg[ipb + iw] *= t_co_interp;
            } else {
                l1rec->tg_sol[ipb + iw] *= pow(t_co[iw], amf_solz);
                l1rec->tg_sen[ipb + iw] *= pow(t_co[iw], amf_senz);
            }
        }
    } else {
        printf("-E- carbon monoxide transmittance can only be computed if gas_opt includes bit  %d\n", GAS_TRANS_TBL_BIT);
        exit(EXIT_FAILURE);

    }
}

void ch4_transmittance(l1str *l1rec, int32_t ip) {
    int32_t iw;

    filehandle *l1file = l1rec->l1file;
    int32_t nwave = l1file->nbands;
    int32_t ipb = ip*nwave;

    /* Only compute if gas transmittance table was requested */
    if ((input->gas_opt & GAS_TRANS_TBL_BIT) != 0) {
        for (iw = 0; iw < nwave; iw++) {
            if (amf) {
                int32_t index=iw*num_airmass;
                
                float t_ch4_interp=t_ch4[index+index_amf_solz]*(1-ratio_solz)+t_ch4[index+index_amf_solz+1]*ratio_solz;
                l1rec->tg_sol[ipb + iw] *= t_ch4_interp;
                
                t_ch4_interp=t_ch4[index+index_amf_total]*(1-ratio_total)+t_ch4[index+index_amf_total+1]*ratio_total;
                l1rec->tg[ipb + iw] *= t_ch4_interp;
            } else {
                l1rec->tg_sol[ipb + iw] *= pow(t_ch4[iw], amf_solz);
                l1rec->tg_sen[ipb + iw] *= pow(t_ch4[iw], amf_senz);
            }
        }
    } else {
        printf("-E- methane transmittance can only be computed if gas_opt includes bit  %d\n", GAS_TRANS_TBL_BIT);
        exit(EXIT_FAILURE);
    }
}

void o2_transmittance(l1str *l1rec, int32_t ip) {
    int32_t iw;

    filehandle *l1file = l1rec->l1file;
    int32_t nwave = l1file->nbands;
    int32_t ipb = ip*nwave;

    /* Only compute if gas transmittance table was requested */
    if ((input->gas_opt & GAS_TRANS_TBL_BIT) != 0) {
        int index_amf_solz_o2;
        float ratio_solz_o2;
        float ratio_total_o2;
        int index_amf_total_o2;

        if (amf) {
            float amf_total_o2 = get_airmass_oxygen(l1rec, ip, 753.0221, 761.7891, 776.81335);
            float scaling_factor = amf_total_o2 / amf_total;

            index_amf_solz_o2 = get_index_lowerbound(amf_mixed, num_airmass, amf_solz * scaling_factor);
            index_amf_total_o2 = get_index_lowerbound(amf_mixed, num_airmass, amf_total * scaling_factor);

            ratio_solz_o2 = (amf_solz * scaling_factor - amf_mixed[index_amf_solz_o2]) /
                            (amf_mixed[index_amf_solz_o2 + 1] - amf_mixed[index_amf_solz_o2]);
            ratio_total_o2 = (amf_total * scaling_factor - amf_mixed[index_amf_total_o2]) /
                             (amf_mixed[index_amf_total_o2 + 1] - amf_mixed[index_amf_total_o2]);
        }

        for (iw = 0; iw < nwave; iw++) {
            if (amf) {
                int32_t index=iw*num_airmass;
                float t_o2_interp;

                if (l1rec->l1file->sensorID != OCI) {
                    t_o2_interp = t_o2[index + index_amf_solz] * (1 - ratio_solz) +
                                  t_o2[index + index_amf_solz + 1] * ratio_solz;
                    l1rec->tg_sol[ipb + iw] *= t_o2_interp;

                    t_o2_interp = t_o2[index + index_amf_total] * (1 - ratio_total) +
                                  t_o2[index + index_amf_total + 1] * ratio_total;
                    l1rec->tg[ipb + iw] *= t_o2_interp;
                } else {
                    t_o2_interp = t_o2[index + index_amf_solz_o2] * (1 - ratio_solz_o2) +
                                  t_o2[index + index_amf_solz_o2 + 1] * ratio_solz_o2;
                    l1rec->tg_sol[ipb + iw] *= t_o2_interp;

                    t_o2_interp = t_o2[index + index_amf_total_o2] * (1 - ratio_total_o2) +
                                  t_o2[index + index_amf_total_o2 + 1] * ratio_total_o2;
                    l1rec->tg[ipb + iw] *= t_o2_interp;
                }
            } else {
                l1rec->tg_sol[ipb + iw] *= pow(t_o2[iw], amf_solz);
                l1rec->tg_sen[ipb + iw] *= pow(t_o2[iw], amf_senz);
            }
        }
    }
}

void n2o_transmittance(l1str *l1rec, int32_t ip) {
    int32_t iw;

    filehandle *l1file = l1rec->l1file;
    int32_t nwave = l1file->nbands;
    int32_t ipb = ip*nwave;

    /* Only compute if gas transmittance table was requested */
    if ((input->gas_opt & GAS_TRANS_TBL_BIT) != 0) {
        for (iw = 0; iw < nwave; iw++) {
            if (amf) {
                int32_t index=iw*num_airmass;
                
                float t_n2o_interp=t_n2o[index+index_amf_solz]*(1-ratio_solz)+t_n2o[index+index_amf_solz+1]*ratio_solz;
                l1rec->tg_sol[ipb + iw] *= t_n2o_interp;
                
                t_n2o_interp=t_n2o[index+index_amf_total]*(1-ratio_total)+t_n2o[index+index_amf_total+1]*ratio_total;
                l1rec->tg[ipb + iw] *= t_n2o_interp;
            } else {
                l1rec->tg_sol[ipb + iw] *= pow(t_n2o[iw], amf_solz);
                l1rec->tg_sen[ipb + iw] *= pow(t_n2o[iw], amf_senz);
            }
        }
    } else {
        printf("-E- nitrous oxide transmittance can only be computed if gas_opt includes bit %d\n", GAS_TRANS_TBL_BIT);
        exit(EXIT_FAILURE);
    }
}

void no2_transmittance(l1str *l1rec, int32_t ip) {

    float a_285, a_225;
    float tau_to200;
    float no2_tr200;
    int32_t iw;

    filehandle *l1file = l1rec->l1file;
    int32_t nwave = l1file->nbands;
    int32_t ipb = ip*nwave;

    float sec0 = 1.0 / l1rec->csolz[ip];
    float sec = 1.0 / l1rec->csenz[ip];

    if (l1rec->no2_tropo[ip] > 0.0)
        /* compute tropo no2 above 200m (Z.Ahmad)
        no2_tr200 = exp(12.6615 + 0.61676*log(no2_tropo));
           new, location-dependent method */
        no2_tr200 = l1rec->no2_frac[ip] * l1rec->no2_tropo[ip];
    else
        no2_tr200 = 0.0;


    for (iw = 0; iw < nwave; iw++) {

        if (l1file->k_no2[iw] > 0.0) {

            a_285 = l1file->k_no2[iw] * (1.0 - 0.003 * (285.0 - 294.0));
            a_225 = l1file->k_no2[iw] * (1.0 - 0.003 * (225.0 - 294.0));

            tau_to200 = a_285 * no2_tr200 + a_225 * l1rec->no2_strat[ip];

            l1rec->tg_sol[ipb + iw] *= exp(-(tau_to200 * sec0));
            if(amf) {
                l1rec->tg    [ipb + iw] *= exp(-(tau_to200 * (sec+sec0)));
            } else {
                l1rec->tg_sen[ipb + iw] *= exp(-(tau_to200 * sec));
            }

        }
    }
}

void h2o_transmittance(l1str *l1rec, int32_t ip) {
    static float *a_h2o = NULL;
    static float *b_h2o = NULL;
    static float *c_h2o = NULL;
    static float *d_h2o = NULL;
    static float *e_h2o = NULL;
    static float *f_h2o = NULL;
    static float *g_h2o = NULL;
    int32_t index;

    float t_h2o;
    int32_t iw;

    filehandle *l1file = l1rec->l1file;
    int32_t nwave = l1file->nbands;
    int32_t ipb = ip*nwave;
    float wv = l1rec->wv[ip];

    if (amf && input->watervapor_bands) {
        wv = 0;
        for (iw = 0; iw < input->nbands_watervapor;) {
            wv += get_wv_band_ratio(l1rec, ip, input->watervapor_bands[iw], input->watervapor_bands[iw + 1],
                                    input->watervapor_bands[iw + 2]);
            iw += 3;
        }
        wv /= (input->nbands_watervapor / 3);
        l1rec->wv[ip]=wv;
    }

    // Apply water vapor transmittance only for the MSE and multi-band AC from the netcdf
    // if (input->aer_opt == AERRHMSEPS || input->aer_opt == AERRHMSEPS_lin || input->aer_opt == AERRHSM) {
    if ((input->gas_opt & GAS_TRANS_TBL_BIT) != 0) {
        int32_t ja = 0,ja_sen=0,ja_sol=0;
        float ratio_wv;
        float f00,f11,f01,f10;
        float ratio_amf_solz,ratio_amf_total,tempratio;
        int32_t index_amf_wv_solz,index_amf_wv_total;

        if (amf) {
            index_amf_wv_solz = get_index_lowerbound(amf_wv, num_airmass, amf_solz);
            index_amf_wv_total = get_index_lowerbound(amf_wv, num_airmass, amf_total);

            ratio_amf_solz = (amf_solz - amf_wv[index_amf_wv_solz]) /
                             (amf_wv[index_amf_wv_solz + 1] - amf_wv[index_amf_wv_solz]);
            ratio_amf_total = (amf_total - amf_wv[index_amf_wv_total]) /
                              (amf_wv[index_amf_wv_total + 1] - amf_wv[index_amf_wv_total]);
        }
        ja = get_index_lowerbound(cwv_all, num_water_vapors, wv );
        ja_sen = get_index_lowerbound(cwv_all, num_water_vapors, wv*amf_senz );
        ja_sol = get_index_lowerbound(cwv_all, num_water_vapors, wv*amf_solz );

        ratio_wv=(wv -cwv_all[ja])/(cwv_all[ja+1]-cwv_all[ja]);
        

        for (iw = 0; iw < nwave; iw++) {
            if (amf) {
                index=model*num_wavelengths*num_airmass*num_water_vapors+iw*num_airmass*num_water_vapors;

                f00 = wvtbl[index+index_amf_wv_solz*num_water_vapors+ja];
                f10 = wvtbl[index+(index_amf_wv_solz+1)*num_water_vapors+ja];
                f01 = wvtbl[index+index_amf_wv_solz*num_water_vapors+ja+1];
                f11 = wvtbl[index+(index_amf_wv_solz+1)*num_water_vapors+ja+1];

                t_h2o = (1. - ratio_amf_solz)*(1. - ratio_wv) * f00 + ratio_amf_solz * ratio_wv * f11 + ratio_amf_solz * (1. - ratio_wv) * f10 + ratio_wv * (1. - ratio_amf_solz) * f01;
                l1rec->tg_sol[ipb + iw] *= t_h2o;

                f00 = wvtbl[index+index_amf_wv_total*num_water_vapors+ja];
                f10 = wvtbl[index+(index_amf_wv_total+1)*num_water_vapors+ja];
                f01 = wvtbl[index+index_amf_wv_total*num_water_vapors+ja+1];
                f11 = wvtbl[index+(index_amf_wv_total+1)*num_water_vapors+ja+1];

                t_h2o = (1. - ratio_amf_total)*(1. - ratio_wv) * f00 + ratio_amf_total * ratio_wv * f11 + ratio_amf_total * (1. - ratio_wv) * f10 + ratio_wv * (1. - ratio_amf_total) * f01;
                l1rec->tg[ipb + iw] *= t_h2o;
            } else {
                index=model*num_wavelengths*num_water_vapors+iw*num_water_vapors;

                tempratio=(wv*amf_solz -cwv_all[ja_sol])/(cwv_all[ja_sol+1]-cwv_all[ja_sol]);
                t_h2o=wvtbl[index+ja_sol]*(1-tempratio)+wvtbl[index+ja_sol+1]*tempratio;
                l1rec->tg_sol[ipb + iw] *= t_h2o;

                tempratio=(wv*amf_senz -cwv_all[ja_sen])/(cwv_all[ja_sen+1]-cwv_all[ja_sen]);
                t_h2o=wvtbl[index+ja_sen]*(1-tempratio)+wvtbl[index+ja_sen+1]*tempratio;
                l1rec->tg_sen[ipb + iw] *= t_h2o;
            }
        }
    }        // otherwise apply Zia's tabel from Bo-cai
    else {
        if (a_h2o == NULL) {
            rdsensorinfo(l1file->sensorID, l1_input->evalmask, "a_h2o", (void **) &a_h2o);
            rdsensorinfo(l1file->sensorID, l1_input->evalmask, "b_h2o", (void **) &b_h2o);
            rdsensorinfo(l1file->sensorID, l1_input->evalmask, "c_h2o", (void **) &c_h2o);
            rdsensorinfo(l1file->sensorID, l1_input->evalmask, "d_h2o", (void **) &d_h2o);
            rdsensorinfo(l1file->sensorID, l1_input->evalmask, "e_h2o", (void **) &e_h2o);
            rdsensorinfo(l1file->sensorID, l1_input->evalmask, "f_h2o", (void **) &f_h2o);
            rdsensorinfo(l1file->sensorID, l1_input->evalmask, "g_h2o", (void **) &g_h2o);
        }

        for (iw = 0; iw < nwave; iw++) {
            t_h2o = a_h2o[iw] + wv * (b_h2o[iw] + wv * (c_h2o[iw] + wv * (d_h2o[iw]
                    + wv * (e_h2o[iw] + wv * (f_h2o[iw] + wv * g_h2o[iw])))));
            l1rec->tg_sol[ipb + iw] *= pow(t_h2o, 1.0 / l1rec->csolz[ip]);
            if(amf) {
                l1rec->tg    [ipb + iw] *= pow(t_h2o, 1.0 / l1rec->csenz[ip]+1.0 / l1rec->csolz[ip]);
            } else {
                l1rec->tg_sen[ipb + iw] *= pow(t_h2o, 1.0 / l1rec->csenz[ip]);
            }
        }
    }
    return;
}

void gaseous_transmittance(l1str *l1rec, int32_t ip) {
    static int32_t firstRun = TRUE;
    int ib, ipb;
    int nwave = l1rec->l1file->nbands;

    if ((input->gas_opt & ATREM_BIT) != 0) {
        if (input->oxaband_opt == 1 && ((input->atrem_opt & ATREM_O2) != 0)){
            printf("ATREM O2 correction is incompatible with the Ding and Gordon approach\n");
            printf("Either unset the atrem_opt bit %d or set oxaband_opt=0\n",ATREM_O2);
            exit(EXIT_FAILURE);
        }
        static float *rhot, *tg_tot;
        float airmass, A;
        
        if (firstRun) {
            if ((rhot = (float *) calloc(nwave, sizeof (float))) == NULL) {
                printf("-E- : Error allocating memory to rhot\n");
                exit(EXIT_FAILURE);
            }
            if ((tg_tot = (float *) calloc(nwave, sizeof (float))) == NULL) {
                printf("-E- : Error allocating memory to tg_tot\n");
                exit(EXIT_FAILURE);
            }
            firstRun = FALSE;
        }
        /* Copy surface reflectance to temp. var for atrem calculation */
        for (ib = 0; ib < nwave; ib++) {
            ipb = ip * nwave + ib;
            rhot [ib] = M_PI * l1rec->Lt[ipb] / l1rec->Fo[ib] / l1rec->csolz[ip];
        }

        get_atrem_cor(l1rec, ip, rhot, tg_tot, &l1rec->tg_sol[ip * nwave], &l1rec->tg_sen[ip * nwave]);
        if (input->atrem_splitpaths == 0) {
            airmass = 1.0 / l1rec->csolz[ip] + 1.0 / l1rec->csenz[ip];

            for (ib = 0; ib < nwave; ib++) {
                ipb = ip * nwave + ib;
                // ATREM didn't do the path splitting, do it here
                A = -log(tg_tot[ib]) / airmass; //effective optical depth
                l1rec->tg_sol[ipb] = exp(-A / l1rec->csolz[ip]);
                l1rec->tg_sen[ipb] = exp(-A / l1rec->csenz[ip]);
                l1rec->tg[ipb] = l1rec->tg_sen[ipb] * l1rec->tg_sol[ipb];
            }
        }
    } else {
        if ((input->gas_opt & GAS_TRANS_TBL_BIT) != 0) {
            if (firstRun) {
                load_gas_tables(l1rec);
                firstRun = FALSE;
            }
            /*
                o2 transmittance is special
                If the oxaband_opt is set to use the gas_transmittance tables, oblige
            */
            amf_solz=1.0/l1rec->csolz[ip];
            amf_senz=1.0/l1rec->csenz[ip];
            amf_total=amf_solz+amf_senz;

            if (amf) {
                index_amf_solz = get_index_lowerbound(amf_mixed, num_airmass, amf_solz);
                index_amf_total = get_index_lowerbound(amf_mixed, num_airmass, amf_total);

                ratio_solz = (amf_solz - amf_mixed[index_amf_solz]) /
                             (amf_mixed[index_amf_solz + 1] - amf_mixed[index_amf_solz]);
                ratio_total = (amf_total - amf_mixed[index_amf_total]) /
                              (amf_mixed[index_amf_total + 1] - amf_mixed[index_amf_total]);
            }

            if (input->oxaband_opt == 2) {
                o2_transmittance(l1rec, ip);
            }
        }

        if ((input->gas_opt & O3_BIT) != 0) {
            ozone_transmittance(l1rec, ip);
        }

        if ((input->gas_opt & CO2_BIT) != 0) {
            co2_transmittance(l1rec, ip);
        }

        if ((input->gas_opt & NO2_BIT) != 0) {
            no2_transmittance(l1rec, ip);
        }

        if ((input->gas_opt & H2O_BIT) != 0) {
            h2o_transmittance(l1rec, ip);
        }

        if ((input->gas_opt & CO_BIT) != 0) {
            co_transmittance(l1rec, ip);
        }
        if ((input->gas_opt & CH4_BIT) != 0) {
            ch4_transmittance(l1rec, ip);
        }
        if ((input->gas_opt & N2O_BIT) != 0) {
            n2o_transmittance(l1rec, ip);
        }
        if (amf) {
            for (ib = 0; ib < nwave; ib++) {
                ipb = ip * nwave + ib;
                l1rec->tg_sen[ipb] = l1rec->tg[ipb] / l1rec->tg_sol[ipb];
            }
        } else {
            for (ib = 0; ib < nwave; ib++) {
                ipb = ip * nwave + ib;
                l1rec->tg[ipb] = l1rec->tg_sen[ipb] * l1rec->tg_sol[ipb];
            }

        }
    }
}

void gas_trans_uncertainty(l1str *l1rec) {
    filehandle *l1file = l1rec->l1file;
    uncertainty_t *uncertainty=l1rec->uncertainty;
    int32_t nwave = l1file->nbands;
    int npix = l1rec->npix;
    int32_t ipb, ip, ib;
    float tg_oz, tg_co2, tg_no2, dt_oz, dt_co2, dt_no2;

    float mu0, mu, tau_oz;

    float a_285, a_225;
    float tau_to200, dtau_to200;
    float no2_tr200, dno2_tr200;

    if(!t_co2)
        load_gas_tables(l1rec);


    for (ip = 0; ip < npix; ip++) {
        mu0 = l1rec->csolz[ip];
        mu = l1rec->csenz[ip];

        if (l1rec->no2_tropo[ip] > 0.0) {
            no2_tr200 = l1rec->no2_frac[ip] * l1rec->no2_tropo[ip];
            dno2_tr200 = l1rec->no2_frac[ip] * uncertainty->dno2_tropo[ip];
        } else {
            no2_tr200 = 0.0;
            dno2_tr200 = 0.0;
        }

        for (ib = 0; ib < nwave; ib++) {
            ipb = ip * nwave + ib;

            tau_oz = l1rec->oz[ip] * l1file->k_oz[ib];

            /* calculate error in tg_sol */
            tg_oz = exp(-tau_oz / mu0);
            dt_oz = tg_oz / mu0 * l1file->k_oz[ib] * uncertainty->doz[ip];

            tg_co2 = pow(t_co2[ib], 1.0 / mu0);
            dt_co2 = 0.0;

            tg_no2 = 1.0;
            dt_no2 = 0.0;

            if (l1file->k_no2[ib] > 0) {
                a_285 = l1file->k_no2[ib] * (1.0 - 0.003 * (285.0 - 294.0));
                a_225 = l1file->k_no2[ib] * (1.0 - 0.003 * (225.0 - 294.0));

                tau_to200 = a_285 * no2_tr200 + a_225 * l1rec->no2_strat[ip];
                dtau_to200 = sqrt(pow(a_285*dno2_tr200, 2) + pow(a_225 * uncertainty->dno2_strat[ip], 2));

                tg_no2 = exp(-tau_to200 / mu0);
                dt_no2 = tg_no2 / mu0*dtau_to200;
            }
            uncertainty->dtg_sol[ipb] = sqrt(pow(tg_no2 * tg_co2*dt_oz, 2) + pow(tg_oz * tg_co2*dt_no2, 2) + pow(tg_no2 * tg_oz*dt_co2, 2));


            /* calculate error in tg_sen */
            tg_oz = exp(-tau_oz / mu);
            dt_oz = tg_oz / mu * l1file->k_oz[ib] * uncertainty->doz[ip];

            tg_co2 = pow(t_co2[ib], 1.0 / mu);
            dt_co2 = 0.0;

            tg_no2 = 1.0;
            dt_no2 = 0.0;
            if (l1file->k_no2[ib] > 0) {
                a_285 = l1file->k_no2[ib] * (1.0 - 0.003 * (285.0 - 294.0));
                a_225 = l1file->k_no2[ib] * (1.0 - 0.003 * (225.0 - 294.0));

                tau_to200 = a_285 * no2_tr200 + a_225 * l1rec->no2_strat[ip];

                dtau_to200 = sqrt(pow(a_285*dno2_tr200, 2) + pow(a_225 * uncertainty->dno2_strat[ip], 2));

                tg_no2 = exp(-tau_to200 / mu);
                dt_no2 = tg_no2 / mu*dtau_to200;
            }
            uncertainty->dtg_sen[ipb] = sqrt(pow(tg_no2 * tg_co2*dt_oz, 2) + pow(tg_oz * tg_co2*dt_no2, 2) + pow(tg_no2 * tg_oz*dt_co2, 2));
        }
    }
}

// TODO: Modify this function to work with non-AMF table
// ...*maybe*, if we push all sensors to use AMF, this would not be necessary
float get_wv_band_ratio(l1str *l1rec,int32_t ip,float window1, float absorp_band,float window2){

    static int32_t firstcall=1;
    static float *tran_interp;

    int32_t i;
    int32_t index;
    int32_t amf_index;
    int32_t band1;
    int32_t band2;
    int32_t band_absorp;
    float wv;
    float rhot_interp;
    float trans_wv_true;
    float amf_ratio;
    float rhot[3];

    filehandle *l1file = l1rec->l1file;

    if(firstcall){
        firstcall=0;
        tran_interp=(float *)malloc(num_water_vapors * sizeof(float));
    }

    // derive a transmittance using a line height (or in this case, depth) approach
    band1 = windex(window1, l1file->fwave, l1file->nbands);
    rhot[0] = M_PI * l1rec->Lt[ip * l1file->nbands + band1] / l1rec->Fo[band1] / l1rec->csolz[ip];

    band2 = windex(window2, l1file->fwave, l1file->nbands);
    rhot[1] = M_PI * l1rec->Lt[ip * l1file->nbands + band2] / l1rec->Fo[band2] / l1rec->csolz[ip];

    band_absorp = windex(absorp_band, l1file->fwave, l1file->nbands);
    rhot[2] = M_PI * l1rec->Lt[ip * l1file->nbands + band_absorp] / l1rec->Fo[band_absorp] / l1rec->csolz[ip];

    rhot_interp = rhot[0] + ((absorp_band - window1) / (window2 - window1)) * (rhot[1] - rhot[0]);
    trans_wv_true = rhot[2] / rhot_interp;

    // For the given absorption band and pixel air mass factor, interpolate the water vapor transmittance table for the
    // each tabular water vapor
    amf_index = get_index_lowerbound(amf_wv, num_airmass, amf_total);
    amf_ratio = (amf_total - amf_wv[amf_index]) / (amf_wv[amf_index + 1] - amf_wv[amf_index]);
    index = (model * num_wavelengths * num_airmass * num_water_vapors) +
             (band_absorp * num_airmass * num_water_vapors) +
             (amf_index * num_water_vapors);

    for (i = 0; i < num_water_vapors; i++) {
        tran_interp[i] = wvtbl[index + i] * (1 - amf_ratio) +
                            wvtbl[index + num_water_vapors + i] *  amf_ratio;
    }

    // Find the bounding transmittance index matching the "true" (computed) transmittance
    index = get_index_upperbound(tran_interp, num_water_vapors, trans_wv_true);
    // retrieve water vapor by interpolating the tabular column water vapor assocaited with the "true" transmittance
    wv = cwv_all[index] + (trans_wv_true - tran_interp[index]) * (cwv_all[index] - cwv_all[index - 1]) /
                            (tran_interp[index] - tran_interp[index - 1]);

    return (wv);
}

float get_airmass_oxygen(l1str *l1rec,int32_t ip,float window1, float absorp_band,float window2){

    float amf_interp;
    int32_t i;

    filehandle *l1file = l1rec->l1file;
    int32_t nwave = l1file->nbands;
    int32_t ipb = ip*nwave;
    float * wave=l1file->fwave;
    float *Lt=&l1rec->Lt[ipb];
    float u0=l1rec->csolz[ip];
    float *F0=l1rec->Fo;
    float rhot_interp,trans_o2_true;
    float rhot[3];
    
    int band1,band2,band_absorp;

    band1=windex(window1,wave,nwave);
    rhot[0]=M_PI*Lt[band1]/F0[band1]/u0;

    band2=windex(window2,wave,nwave);
    rhot[1]=M_PI*Lt[band2]/F0[band2]/u0;

    band_absorp=windex(absorp_band,wave,nwave);
    rhot[2]=M_PI*Lt[band_absorp]/F0[band_absorp]/u0;

    rhot_interp=rhot[0]+(absorp_band-window1)*(rhot[1]-rhot[0])/(window2-window1);

    trans_o2_true=rhot[2]/rhot_interp;

    ipb=band_absorp*num_airmass;
    for (i = 0; i < num_airmass; i++) {
        if (trans_o2_true >= t_o2[ipb+i])
            break;
    }
    if (i == 0)
        i = 1;
    if (i == num_airmass)
        i = num_airmass - 1;

    amf_interp = amf_mixed[i] + (trans_o2_true - t_o2[ipb+i]) * (amf_mixed[i] - amf_mixed[i - 1]) /
                          (t_o2[ipb+i] - t_o2[ipb+i - 1]);

    return (amf_interp);
}
