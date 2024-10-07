/* =========================================================== */
/* Module gaseous_transmittance.c                              */
/*                                                             */
/* Computes sensor-specific transmittance through various      */
/* atmospheric gases.                                          */
/*                                                             */
/* B. Franz, NASA/OBPG, July 2006                              */
/* =========================================================== */

#include "l12_proto.h"
#include <allocate3d.h>
#include "atrem_corl1.h"

size_t num_models, num_wavelengths, num_water_vapors;
static float ***wvtbl = NULL;
static float *t_co2 = NULL;
static float *t_o2 = NULL;
static float *t_co = NULL;
static float *t_ch4 = NULL;
static float *t_n2o = NULL;
static float *cwv_all = NULL;

void load_gas_tables(l1str *l1rec) {
    /*
        netcdf file that has water vapor transmittance as a function of wavelength and
        cwv (6 profiles x number of bands x 220 water vapor value)
        filename = /OCDATAROOT/sensor[/subsensor]/<sensorName>_gas_trans.nc
    */
    char *filedir;
    char filename[FILENAME_MAX];
    if ((filedir = getenv("OCDATAROOT")) == NULL) {
        printf("-E- %s: OCDATAROOT env variable undefined.\n", __FILE__);
        return;
    }
    strcpy(filename, filedir);

    strcat(filename, "/");
    strcat(filename, sensorId2SensorDir(l1rec->l1file->sensorID));

    // SeaWiFS has a subsensorID, but the gas transmittance table isn't GAC/LAC specific
    if ((l1rec->l1file->sensorID != SEAWIFS) && (l1rec->l1file->subsensorID != -1)) {
        strcat(filename, "/");
        strcat(filename, subsensorId2SubsensorDir(l1rec->l1file->subsensorID));
    }

    char *sensorName = strdup(sensorId2SensorName(l1rec->l1file->sensorID));
    lowcase(sensorName);
    strcat(filename, "/");
    strcat(filename, sensorName);
    strcat(filename, "_gas_transmittance.nc");
    free(sensorName);

    /* This will be the netCDF ID for the file and data variable. */
    int32_t ncid, varid;
    int32_t num_water_vapors_id, num_models_id, num_wavelengths_id;

    /* Open the file */
    if ((nc_open(filename, NC_NOWRITE, &ncid)) != NC_NOERR) {
        printf("-E- %s: Failed to open %s\n", __FILE__, filename);
        exit(EXIT_FAILURE);
    }

    if((nc_inq_dimid(ncid, "n_water_vapor", &num_water_vapors_id)) == NC_NOERR){
        if((nc_inq_dimlen(ncid, num_water_vapors_id, &num_water_vapors)) != NC_NOERR){
            printf("-E- %s: Failed to read dimension n_water_vapor\n", __FILE__);
            exit(EXIT_FAILURE);
        }
    } else{
        printf("-E- %s: Failed to find dimension n_water_vapor\n", __FILE__);
        exit(EXIT_FAILURE);
    }
    if((nc_inq_dimid(ncid, "nmodels", &num_models_id)) == NC_NOERR){
        if((nc_inq_dimlen(ncid, num_models_id, &num_models)) != NC_NOERR){
            printf("-E- %s: Failed to read dimension nmodels\n", __FILE__);
            exit(EXIT_FAILURE);
        }
    } else{
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
    wvtbl = allocate3d_float(num_models, num_wavelengths, num_water_vapors);
    if(!wvtbl) {
        printf("Error: allocating memory for water vapor transmittance tables\n");
        exit(EXIT_FAILURE);
    }

    if ((nc_inq_varid(ncid, "water_vapor_transmittance", &varid)) == NC_NOERR) {
        if ((nc_get_var_float(ncid, varid, &wvtbl[0][0][0])) != NC_NOERR){
            printf("-E- %s: failed to read water_vapor_transmittace from %s\n",
                        __FILE__, filename);
            exit(EXIT_FAILURE);
        }
    } else {
        printf("-E- %s: '%s' does not have water_vapor_transmittace.\n",
                    __FILE__, filename);
        exit(EXIT_FAILURE);
    }

    /* Read the water vapor table */
    if ((cwv_all = (float *) calloc(num_water_vapors, sizeof(float))) == NULL) {
        printf("Error: allocating memory for water vapor table\n");
        exit(EXIT_FAILURE);
    }
    if ((nc_inq_varid(ncid, "water_vapor", &varid)) == NC_NOERR) {
        if ((nc_get_var_float(ncid, varid, &cwv_all[0])) != NC_NOERR){
            printf("-E- %s: failed to read water_vapor from %s\n",
                        __FILE__, filename);
            exit(EXIT_FAILURE);
        }
    } else {
        printf("-E- %s: '%s' does not have water_vapor\n",
                    __FILE__, filename);
        exit(EXIT_FAILURE);
    }

    /* Read the carbon monoxide transmittance */
    if ((t_co = (float *) calloc(num_wavelengths, sizeof(float))) == NULL) {
        printf("-E- %s line %d : error allocating memory for carbon monoxide transmitance table.\n",
                __FILE__, __LINE__);
        exit(1);
    }
    if ((nc_inq_varid(ncid, "carbon_monoxide_transmittance", &varid)) == NC_NOERR) {
        if ((nc_get_var_float(ncid, varid, &t_co[0])) != NC_NOERR){
            printf("-E- %s: failed to read carbon_monoxide_transmittance from %s\n",
                        __FILE__, filename);
            exit(EXIT_FAILURE);
        }
    } else {
        printf("-E- %s: '%s' does not have carbon_monoxide_transmittance\n",
                    __FILE__, filename);
        exit(EXIT_FAILURE);
    }

    /* Read the carbon dioxide transmittance */
    if ((t_co2 = (float *) calloc(num_wavelengths, sizeof(float))) == NULL) {
        printf("-E- %s line %d : error allocating memory for carbon dioxide transmitance table.\n",
                __FILE__, __LINE__);
        exit(1);
    }
    if ((nc_inq_varid(ncid, "carbon_dioxide_transmittance", &varid)) == NC_NOERR) {
        if ((nc_get_var_float(ncid, varid, &t_co2[0])) != NC_NOERR){
            printf("-E- %s: failed to read carbon_dioxide_transmittance from %s\n",
                        __FILE__, filename);
            exit(EXIT_FAILURE);
        }
    } else {
        printf("-E- %s: '%s' does not have carbon_dioxide_transmittance\n",
                    __FILE__, filename);
        exit(EXIT_FAILURE);
    }

    /* Read the oxygen transmittance */
    if ((t_o2 = (float *) calloc(num_wavelengths, sizeof(float))) == NULL) {
        printf("-E- %s line %d : error allocating memory for oxygen transmitance table.\n",
                __FILE__, __LINE__);
        exit(1);
    }
    if ((nc_inq_varid(ncid, "oxygen_transmittance", &varid)) == NC_NOERR) {
        if ((nc_get_var_float(ncid, varid, &t_o2[0])) != NC_NOERR){
            printf("-E- %s: failed to read oxygen_transmittance from %s\n",
                        __FILE__, filename);
            exit(EXIT_FAILURE);
        }
    } else {
        printf("-E- %s: '%s' does not have oxygen_transmittance\n",
                    __FILE__, filename);
        exit(EXIT_FAILURE);
    }

    /* Read the nitrous oxide transmittance */
    if ((t_n2o = (float *) calloc(num_wavelengths, sizeof(float))) == NULL) {
        printf("-E- %s line %d : error allocating memory for nitrous oxide transmitance table.\n",
                __FILE__, __LINE__);
        exit(1);
    }
    if ((nc_inq_varid(ncid, "nitrous_oxide_transmittance", &varid)) == NC_NOERR) {
        if ((nc_get_var_float(ncid, varid, &t_n2o[0])) != NC_NOERR){
            printf("-E- %s: failed to read nitrous_oxide_transmittance from %s\n",
                        __FILE__, filename);
            exit(EXIT_FAILURE);
        }
    } else {
        printf("-E- %s: '%s' does not have nitrous_oxide_transmittance\n",
                    __FILE__, filename);
        exit(EXIT_FAILURE);
    }

    /* Read the methane transmittance */
    if ((t_ch4 = (float *) calloc(num_wavelengths, sizeof(float))) == NULL) {
        printf("-E- %s line %d : error allocating memory for methane transmitance table.\n",
                __FILE__, __LINE__);
        exit(1);
    }
    if ((nc_inq_varid(ncid, "methane_transmittance", &varid)) == NC_NOERR) {
        if ((nc_get_var_float(ncid, varid, &t_ch4[0])) != NC_NOERR){
            printf("-E- %s: failed to read methane_transmittance from %s\n",
                        __FILE__, filename);
            exit(EXIT_FAILURE);
        }
    } else {
        printf("-E- %s: '%s' does not have methane_transmittance\n",
                    __FILE__, filename);
        exit(EXIT_FAILURE);
    }

    /* Close the file */
    if ((nc_close(ncid)) != NC_NOERR){
        printf("-E- %s: failed to close %s\n",
            __FILE__, filename);
        exit(EXIT_FAILURE);
    }
}

int32_t get_wvindex(float *wvtable,int32_t nwv, double wv)
{
	int32_t index;
	int32_t i;

	if(wv>wvtable[nwv-1])
		i=nwv-1;
	else
	{
		for(i=0;i<nwv;i++)
			if(wv<wvtable[i])
				break;
	}

	index=i;
	if( (wv-wvtable[i-1]) < (wvtable[i]-wv) )
		index=i-1;

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
        l1rec->tg_sen[ipb + iw] *= exp(-(tau_oz / l1rec->csenz[ip]));
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
        if(l1rec->l1file->sensorID == OCI) {
            l1rec->tg_sol[ipb + iw] *= pow(pow(t_co2[iw],0.5), 1.0 / l1rec->csolz[ip]);
            l1rec->tg_sen[ipb + iw] *= pow(pow(t_co2[iw],0.5), 1.0 / l1rec->csenz[ip]);
        } else {
            l1rec->tg_sol[ipb + iw] *= pow(t_co2[iw], 1.0 / l1rec->csolz[ip]);
            l1rec->tg_sen[ipb + iw] *= pow(t_co2[iw], 1.0 / l1rec->csenz[ip]);
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
            if(l1rec->l1file->sensorID == OCI) {
                l1rec->tg_sol[ipb + iw] *= pow(pow(t_co[iw],0.5), 1.0 / l1rec->csolz[ip]);
                l1rec->tg_sen[ipb + iw] *= pow(pow(t_co[iw],0.5), 1.0 / l1rec->csenz[ip]);
            } else {
                l1rec->tg_sol[ipb + iw] *= pow(t_co[iw], 1.0 / l1rec->csolz[ip]);
                l1rec->tg_sen[ipb + iw] *= pow(t_co[iw], 1.0 / l1rec->csenz[ip]);
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
            if(l1rec->l1file->sensorID == OCI) {
                l1rec->tg_sol[ipb + iw] *= pow(pow(t_ch4[iw],0.5), 1.0 / l1rec->csolz[ip]);
                l1rec->tg_sen[ipb + iw] *= pow(pow(t_ch4[iw],0.5), 1.0 / l1rec->csenz[ip]);
            } else {
                l1rec->tg_sol[ipb + iw] *= pow(t_ch4[iw], 1.0 / l1rec->csolz[ip]);
                l1rec->tg_sen[ipb + iw] *= pow(t_ch4[iw], 1.0 / l1rec->csenz[ip]);
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
        for (iw = 0; iw < nwave; iw++) {
            if(l1rec->l1file->sensorID == OCI) {
                l1rec->tg_sol[ipb + iw] *= pow(pow(t_o2[iw],0.5), 1.0 / l1rec->csolz[ip]);
                l1rec->tg_sen[ipb + iw] *= pow(pow(t_o2[iw],0.5), 1.0 / l1rec->csenz[ip]);
            } else {
                l1rec->tg_sol[ipb + iw] *= pow(t_o2[iw], 1.0 / l1rec->csolz[ip]);
                l1rec->tg_sen[ipb + iw] *= pow(t_o2[iw], 1.0 / l1rec->csenz[ip]);
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
            if(l1rec->l1file->sensorID == OCI) {
                l1rec->tg_sol[ipb + iw] *= pow(pow(t_n2o[iw],0.5), 1.0 / l1rec->csolz[ip]);
                l1rec->tg_sen[ipb + iw] *= pow(pow(t_n2o[iw],0.5), 1.0 / l1rec->csenz[ip]);
            } else {
                l1rec->tg_sol[ipb + iw] *= pow(t_n2o[iw], 1.0 / l1rec->csolz[ip]);
                l1rec->tg_sen[ipb + iw] *= pow(t_n2o[iw], 1.0 / l1rec->csenz[ip]);
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
            l1rec->tg_sen[ipb + iw] *= exp(-(tau_to200 * sec));

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
    int32_t model = 5;

    float t_h2o;
    int32_t iw;

    filehandle *l1file = l1rec->l1file;
    int32_t nwave = l1file->nbands;
    int32_t ipb = ip*nwave;
    float wv = l1rec->wv[ip];

    // Apply water vapor transmittance only for the MSE and multi-band AC from the netcdf
    // if (input->aer_opt == AERRHMSEPS || input->aer_opt == AERRHMSEPS_lin || input->aer_opt == AERRHSM) {
        if ((input->gas_opt & GAS_TRANS_TBL_BIT) != 0) {
        int32_t ja_sol = 0;
        int32_t ja_sen = 0;

        ja_sol = get_wvindex(cwv_all,num_water_vapors,wv / l1rec->csolz[ip]);
        ja_sen = get_wvindex(cwv_all,num_water_vapors,wv / l1rec->csenz[ip]);

        for (iw = 0; iw < nwave; iw++) {
            if(l1rec->l1file->sensorID == OCI) {
                l1rec->tg_sol[ipb + iw] *= pow(wvtbl[model][iw][ja_sol],0.5);
                l1rec->tg_sen[ipb + iw] *= pow(wvtbl[model][iw][ja_sen],0.5);
            } else {
                l1rec->tg_sol[ipb + iw] *= wvtbl[model][iw][ja_sol];
                l1rec->tg_sen[ipb + iw] *= wvtbl[model][iw][ja_sen];
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
            l1rec->tg_sen[ipb + iw] *= pow(t_h2o, 1.0 / l1rec->csenz[ip]);
        }
    }
    return;
}

void gaseous_transmittance(l1str *l1rec, int32_t ip) {
    static int32_t firstRun = TRUE;

    if ((input->gas_opt & ATREM_BIT) != 0) {
        if (input->oxaband_opt == 1 && ((input->atrem_opt & ATREM_O2) != 0)){
            printf("ATREM O2 correction is incompatible with the Ding and Gordon approach\n");
            printf("Either unset the atrem_opt bit %d or set oxaband_opt=0\n",ATREM_O2);
            exit(EXIT_FAILURE);
        }
        static float *rhot, *tg_tot;
        float airmass, A;
        int ib, ipb;
        int nwave = l1rec->l1file->nbands;
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
