/* =================================================================== */
/* module raman.c - Raman scattering correction for Rrs                */
/*                                                                     */
/* This module contains the functions to correct sensor-observed       */
/* above-water remote sensing reflectances, Rrs. for Raman Scattering  */
/* effects.                                                            */
/*                                                                     */
/* References:                                                         */
/* Lee et al (1994) Model for the interpretation of hyperspectral      */
/* remote-sensing reflectance, Applied Optics, 33(24), 5721-5732       */
/* doi:10.1364/AO.33.005721                                            */
/*                                                                     */
/* Lee et al (2013) Penetration of UV-visible solar radiation in the   */
/* global oceans: Insights from ocean color remote sensing, Journal of */
/* Geophysical Research Oceans, 118, 4241-4255, doi:10.1002/jgrc.20308 */
/*                                                                     */
/* McKinna et al. (2016) Implementation of an analytical Raman         */
/* scattering correction for satellite ocean-color processing,         */
/* Opt. Express, 24(14), A1123-A1137, doi: 10.1364/OE.24.0A1123        */
/*                                                                     */
/* Mobley. C.D. (2010) Hydrolight Technical Note 10: Interpretation    */
/* of Raman Scattering Computations, Sequoia Scientific.               */
/*                                                                     */
/* Sathyendranath and Platt (1998) Ocean-color model incorporating     */
/* transspectral processes, Applied Optics, 37(12), 2216-2227,         */
/* doi:10.1364/AO.37.002216                                            */
/*                                                                     */
/* Walrafen (1967) Raman spectral studies of the effects of temperature*/
/* on water structure, Journal of Chemical Physics, 47, 118-121,       */
/* doi:10.1063/1.711834                                                */
/*                                                                     */
/* Westberry et al. (2013) Influence if Raman scattering on ocean color*/
/* inversion models, Applied Optics, 52(22), 5552-5561,                */
/* doi:10.1364/AO.52.005552                                            */
/*                                                                     */
/*                                                                     */
/* Implementation:                                                     */
/* L. McKinna NASA/OBPG/SAIC, April 2015                               */
/*                                                                     */
/* Notes:                                                              */
/* Feb 2019: Downwelling irradiance model updated to be consistent     */
/* with atmospheric correction model. Added support for OLCIA, OLCIB,  */
/* SGLI, and VIIRS-J1.                                                 */
/*                                                                     */ 
/* May 2021: replaced hard-coded coefficients with netCDF files.       */
/* Added support for PACE OCI.                                         */
/* =================================================================== */

#define NORAMAN 0
#define LEE2013 1
#define WESTBERRY2013 2
#define LEE1994 3
#define WAVEMIN 380
#define WAVEMAX 700

#include <stdio.h>
#include <math.h>
#include <netcdf.h>
#include <gsl/gsl_fit.h>
#include "l12_proto.h"
#include "sensorDefs.h"

static int32_t nbandVis; //Number of visible bands
static float *lambda; //Sensor wavelengths

//Wavelength indices
static int idx412;
static int idx443;
static int idx488;
static int idx550;
static int idx670;
static int bandNum550; //Number of bands shorter than 500 nm;
static int limitVisFlag;

//static float angstEst; //Estimate Angstrom exponent (Taua power law slope)

static float nWater = 1.344; //Refractive index of seawater

//Bio-optical model parameters
static float *aw; //water absorption coefficient at sensor bands
static float *bbw; //water backscattering coefficient at sensor bands
static float *awR; //water absorption coeficient at Raman excitation bands
static float *bbwR; //water backscattering coefficient at Raman excitation bandss
static float *bbtot; //total backscattering modelled at sensor bands
static float *atot; //total absorption modelled at sensor bands
static float *atotR; //total absorption coefficient at Raman excitation bands
static float *bbtotR; //total backscattering coefficient at Raman excitation bands
static float *aphStar; //Normalised phytoplankton absorption spectral shape at sensor bands
static float *aphStarR; //Normalised phytoplankton absorption spectral shape at Raman bands
static float *kdSen; //Diffuse attenuation coefficient at sensor bands
static float *kdRam; //Diffuse attenuation coefficient at raman bands
static float *kappaSen; 
static float *kappaRam;
// at sensor bands
static float aph443; //QAA-derived phytoplankton absorption at 443 nm      
static float adg443; //QAA-derived Absorption of CDOM and colored detrital matter at 443 nm
static float bbp443; //QAA-derived backscattering coefficient at 443 nm
static float Ybbp; //Power exponent of bbp
static float Sdg; //Exponential slope of adg

static float *aph1Sen; //Bricaud aph0 coefficient at sensor bands
static float *aph2Sen; //Bricaud aph1 coefficient at sensor bands
static float *aph1Ram; //Bricaud aph0 coefficient at Raman bands
static float *aph2Ram; //Bricaud aph1 coefficient at Raman bands


//Radtran atmospheric model parameters
static float *eDRam; //Down-welling irradiance at Raman excitation bands
static float *eDSen; //Down-welling irradiance at sensor bands

static float *t_sol_ram; //Down-welling irradiance at Raman excitation bands

//static float *tauaRam; //Aerosol optical thickness at Raman bands

static float *ramLam; //Raman excitation wavelengths
static float *ramBandW; //Raman band widths
static float *bRex; //Raman scattering coefficient

static float *f0BarRam; //TOA solar irradiance at Raman bands
static float *rayRam; //Rayleigh transmittance at Raman bands
static float *ozRam; //Ozone absorption at Raman bands
static float *wvRam; //Water vapor absorption at Raman bands
static float *oxyRam; //Oxygen absorption at Raman bands
static float *no2Ram; //Oxygen absorption at Raman bands

static float *Rrs_ram; //Rrs due to Raman scattering

//Lee et al. 2013 Raman correction 1 coefficients
//Pointers for interpolation/extrapolation
static float *alphaCor;
static float *betaCor0;
static float *betaCor1;

//Lee et al 2013 Empirical coefficients
static float lamCor1[6] = {412.0, 443.0, 488.0, 531.0, 551.0, 667.0};
static float alpha[6] = {0.003, 0.004, 0.011, 0.015, 0.017, 0.018};
static float beta0[6] = {0.014, 0.015, 0.010, 0.010, 0.010, 0.010};
static float beta1[6] = {-0.022, -0.023, -0.051, -0.070, -0.080, -0.081};


/* ------------------------------------------------------------------- */
/* getRamanCoeffs() get spectral coefficients from netCDF file         */
/* ------------------------------------------------------------------- */

void getRamanCoeff(int grpid, char* varName, float* variable) {
    
    int varid;
    int dimid;
    size_t varLen;
    
    //get the variable id
    if (nc_inq_varid(grpid, varName, &varid) != NC_NOERR) {
        fprintf(stderr, "-E- %s line %d: could not find netCDF variable \"%s\" in netCDF File.\n",
                __FILE__, __LINE__, varName);
        exit(1);
    }
    
    //get the dimension id for the variable
    if (nc_inq_dimid(grpid,varName, &dimid) != NC_NOERR ) {
        fprintf(stderr, "-E- %s line %d: could not find %s in netCDF File.\n",
                __FILE__, __LINE__, varName);
    }
    
    //get the dimension length of the variable
    if (nc_inq_dimlen(grpid, dimid, &varLen) != NC_NOERR ) {
        fprintf(stderr, "-E- %s line %d: could not find netCDF variable dimensions in netCDF File.\n",
                __FILE__, __LINE__);
        exit(1);   
    }
      
    //Check to see if the dimensions of the coefficient arrays are consistent with
    //visible bands of the sensor
    if (varLen != nbandVis && limitVisFlag!=1) {
        fprintf(stderr, "-E- %s line %d: number of Raman coefficients %ld does not equal sensors visible"
                " bands %d.\n",
                __FILE__, __LINE__,varLen,nbandVis);
        exit(1);   
    }
    
    //get the spectral coefficient from the file
    if (nc_get_var_float (grpid, varid, variable) != NC_NOERR ) { 
        fprintf(stderr, "-E- %s line %d: could not read netCDF variable %s in netCDF File.\n",
                __FILE__, __LINE__,varName);
        exit(1);   
    }
         	
}

/* ---------------------------------------------------------------------------*/
/* raman_pixel_alloc() Allocate pointer memory for Raman computations         */
/* ---------------------------------------------------------------------------*/
void raman_pixel_alloc(l2str *l2rec) {

    //Allocate static pointer memory
    lambda = (float*) allocateMemory(nbandVis * sizeof (float), "lambda");

    aw = (float*) allocateMemory(nbandVis * sizeof (float), "aw");
    bbw = (float*) allocateMemory(nbandVis * sizeof (float), "bbw");
    awR = (float*) allocateMemory(nbandVis * sizeof (float), "awR");
    bbwR = (float*) allocateMemory(nbandVis * sizeof (float), "bbwR");
    aphStar = (float*) allocateMemory(nbandVis * sizeof (float), "aphStar");
    aphStarR = (float*) allocateMemory(nbandVis * sizeof (float), "aphStarR");
    bbtot = (float*) allocateMemory(nbandVis * sizeof (float), "bbtot");
    atot = (float*) allocateMemory(nbandVis * sizeof (float), "atot");
    bbtotR = (float*) allocateMemory(nbandVis * sizeof (float), "bbtotR");
    atotR = (float*) allocateMemory(nbandVis * sizeof (float), "atotR");

    kdSen = (float*) allocateMemory(nbandVis * sizeof (float), "kdSen");
    kdRam = (float*) allocateMemory(nbandVis * sizeof (float), "kdRam");
    kappaSen = (float*) allocateMemory(nbandVis * sizeof (float), "kappaSen");
    kappaRam = (float*) allocateMemory(nbandVis * sizeof (float), "kappaRam");

    betaCor0 = (float*) allocateMemory(nbandVis * sizeof (float), "betaCor0");
    betaCor1 = (float*) allocateMemory(nbandVis * sizeof (float), "betaCor1");
    alphaCor = (float*) allocateMemory(nbandVis * sizeof (float), "alphaCor");

    eDRam = (float*) allocateMemory(nbandVis * sizeof (float), "eDRam");
    eDSen = (float*) allocateMemory(nbandVis * sizeof (float), "eDSen");
    
    t_sol_ram = (float*) allocateMemory(nbandVis * sizeof (float), "t_sol_ram");
    
    Rrs_ram = (float*) allocateMemory(l2rec->l1rec->npix *
            l2rec->l1rec->l1file->nbands * sizeof (float), "Rrs_ram");
    
    ramLam =  (float*) allocateMemory(nbandVis * sizeof (float), "ramLam");
    ramBandW =  (float*) allocateMemory(nbandVis * sizeof (float), "ramBandW");
    bRex =  (float*) allocateMemory(nbandVis * sizeof (float), "bRex");   
    f0BarRam =  (float*) allocateMemory(nbandVis * sizeof (float), "foBarRam");
    rayRam =  (float*) allocateMemory(nbandVis * sizeof (float), "rayRam");
    ozRam =  (float*) allocateMemory(nbandVis * sizeof (float), "ozRam");
    wvRam =  (float*) allocateMemory(nbandVis * sizeof (float), "wvRam");
    oxyRam =  (float*) allocateMemory(nbandVis * sizeof (float), "oxyRam");
    no2Ram =  (float*) allocateMemory(nbandVis * sizeof (float), "no2Ram");
    aph1Ram =  (float*) allocateMemory(nbandVis * sizeof (float), "aph1Ram");
    aph2Ram =  (float*) allocateMemory(nbandVis * sizeof (float), "aph2Ram");
    
    aph1Sen =  (float*) allocateMemory(nbandVis * sizeof (float), "aph1Sen");
    aph2Sen =  (float*) allocateMemory(nbandVis * sizeof (float), "aph2Sen");

}

/* ---------------------------------------------------------------------------*/
/* get_raman_coeffs() Allocate pointer memory for Raman computations         */

/* ---------------------------------------------------------------------------*/
void get_raman_coeffs(l2str *l2rec) {
    
    char filename[FILENAME_MAX];
    char *filedir;
    
    if ((filedir = getenv("OCDATAROOT")) == NULL) {
        printf("-E- %s: OCDATAROOT env variable undefined.\n", __FILE__);
        exit(1);
    }
    
    filehandle *l1file = l2rec->l1rec->l1file;
    strcpy(filename, filedir);
    strcat(filename, "/");
    strcat(filename, sensorId2SensorDir(l1file->sensorID));
    if(l1file->subsensorID > SEAWIFS_LAC) {  // ignore seawifs subsensor ID
        strcat(filename, "/");
        strcat(filename, subsensorId2SubsensorDir(l1file->subsensorID));
    }
    strcat(filename, "/raman.nc");
   
    if(l1file->sensorID == OCIS) {
        limitVisFlag = 1;
    } else if(l1file->sensorID == OCI) {
        limitVisFlag = 1;
    } else {
        limitVisFlag = 0;
    }
    
    printf("Loading Raman coefficients from: %s.\n", filename);
    
    int ncid;
    int grpid;
    //Open netCDF file
    if (nc_open(filename, NC_NOWRITE, &ncid) != NC_NOERR) {
        fprintf(stderr, "-E- %s line %d: could not open netCDF File \"%s\".\n",
                __FILE__, __LINE__, filename);
        exit(1);
    }
    
    //get the group id for raman_coeffs
    if (nc_inq_grp_ncid(ncid, "raman_coeffs", &grpid) != NC_NOERR) {
        fprintf(stderr, "-E- %s line %d: could not find netCDF group \"raman_coeffs\".\n",
                __FILE__, __LINE__);
        exit(1);
    }
    
    getRamanCoeff(grpid, "ramLam",ramLam);
    getRamanCoeff(grpid, "ramBandW",ramBandW);
    getRamanCoeff(grpid, "bRex",bRex);
    getRamanCoeff(grpid, "f0BarRam",f0BarRam);
    getRamanCoeff(grpid, "rayRam",rayRam);
    getRamanCoeff(grpid, "ozRam",ozRam);
    getRamanCoeff(grpid, "wvRam",wvRam);
    getRamanCoeff(grpid, "oxyRam",oxyRam);
    getRamanCoeff(grpid, "no2Ram",no2Ram);
    getRamanCoeff(grpid, "aph1Ram",aph1Ram);
    getRamanCoeff(grpid, "aph2Ram",aph2Ram);

    //get the group id for sensor_coeffs
    if (nc_inq_grp_ncid(ncid, "sensor_coeffs", &grpid) != NC_NOERR) {
        fprintf(stderr, "-E- %s line %d: could not find netCDF group \"sensor_coeffs\".\n",
                __FILE__, __LINE__);
        exit(1);
    }
    
    getRamanCoeff(grpid, "aph1Sen",aph1Sen);
    getRamanCoeff(grpid, "aph2Sen",aph2Sen);

    //close the netcdf
    if (nc_close(ncid) != NC_NOERR){
        fprintf(stderr, "-E- %s line %d: could not close netCDF \"%s\".\n",
                __FILE__, __LINE__,filename);
        exit(1);
    }

}
/* -----------------------------------------------------------------------*/
/* set_raman_aph_uv() - Extrapolates aph* into the UV (below 400 nm).     */
/*                      This methods assumes that aph can be linearly     */
/*                      extrapolated into the UV                          */

/* -----------------------------------------------------------------------*/
void set_raman_aph_uv(l2str *l2rec, int ip) {

    int iw;
    float aphM, aphC;
    float aphStar443;

    //Calculate aphStar at 443 nm using Bricaud et al 1998 
    aphStar443 = aph1Sen[idx443] * pow(l2rec->chl[ip], (aph2Sen[idx443]));

    //Using sensor-specific coefficients, calculate aphstar at all bands and
    //normalize to 1.0 at 443 nm
    for (iw = 0; iw < nbandVis; iw++) {
        aphStar[iw] = (aph1Sen[iw] * pow(l2rec->chl[ip], (aph2Sen[iw])))
                / aphStar443;
        aphStarR[iw] = (aph1Ram[iw] * pow(l2rec->chl[ip], (aph2Ram[iw])))
                / aphStar443;

    }

    //Linear model of aph betweeen 412 and 443 nm
    aphM = (aphStar[idx443] - aphStar[idx412]) / (lambda[idx443] - lambda[idx412]);
    aphC = aphStar[idx443] - aphM * lambda[idx443];

    //If outside the Bricaud in UV, extrapolate as linear function.        
    for (iw = 0; iw < nbandVis; iw++) {
        if (ramLam[iw] <= 400) {
            aphStarR[iw] = (aphM * ramLam[iw] + aphC);
        }
    }

}

/* -----------------------------------------------------------------------*/
/* fit_raman_tsol() - Interpolates/extrapolates tramsittanc.e from sun to  */
/*                      surface (t_sol) at Raman excitation bands. The    */
/*                      assumption is that t_sol follows a power law.     */
/* -----------------------------------------------------------------------*/
void fit_raman_tsol(l2str *l2rec, int ip) {

    int iw;

    double logTsol[nbandVis];
    double waveRatio[nbandVis];

    double c0, c1, cov00, cov01, cov11, sumsq;
    int nbands = l2rec->l1rec->l1file->nbands;

    //Log transform data so function has a linear form
    for (iw = 0; iw < nbandVis; iw++) {
        logTsol[iw] = log(l2rec->l1rec->t_sol[ip * nbands + iw] / l2rec->l1rec->t_sol[ip * nbands + idx488]);
        waveRatio[iw] = log(lambda[iw] / lambda[idx488]);
    }

    //Use a linear fit on long-transformed data to compute power law exponent (range 412-490 nm)
    gsl_fit_linear(waveRatio, 1, logTsol, 1, nbandVis, &c0, &c1, &cov00, &cov01, &cov11, &sumsq);


    //Extrapolate/interpolate using power law model
    for (iw = 0; iw < nbandVis; iw++) {
        t_sol_ram[iw] = l2rec->l1rec->t_sol[ip * nbands + idx488] * pow((ramLam[iw] / lambda[idx488]), c1);
    }

}

/* -----------------------------------------------------------------------*/
/* raman_qaa() - computes a first order estimate of the constituent IOPs  */
/* -----------------------------------------------------------------------*/
void raman_qaa(l2str *l2rec, int ip) {

    int iw, idxref;

    float aref, bbpref, numer, denom, rho;
    float rat, zeta, xi, term1, term2;

    float rrs_a[nbandVis];
    float rrs_s[nbandVis];
    float u[nbandVis];
    float at[nbandVis];
    float bbt[nbandVis];

    float g0 = 0.08945;
    float g1 = 0.1247;
    float acoefs[3] = {-1.146, -1.366, -0.469};

    //Calculate the sub-surface remote sensing reflectance
    for (iw = 0; iw < nbandVis; iw++) {
        rrs_a[iw] = l2rec->Rrs[ip * l2rec->l1rec->l1file->nbands + iw];
        rrs_s[iw] = rrs_a[iw] / (0.52 + 1.7 * rrs_a[iw]);
    }

    /*pre-test of Rrs550*/
    if (rrs_a[idx550] <= 0.0)
        rrs_a[idx550] = 0.001;

    /* pre-test 670 */
    if ((rrs_a[idx670] > 20.0 * pow(rrs_a[idx550], 1.5)) || (rrs_a[idx670] < 0.9 * pow(rrs_a[idx550], 1.7))) {
        rrs_a[idx670] = 1.27 * pow(rrs_a[idx550], 1.47) + 0.00018 * pow(rrs_a[idx488] / rrs_a[idx550], -3.19);
    }

    /*Quadratic formula to get the quantity u*/
    for (iw = 0; iw < nbandVis; iw++) {
        u[iw] = (sqrt(g0 * g0 + 4.0 * g1 * rrs_s[iw]) - g0) / (2.0 * g1);
    }

    /*Determine whether to use 550 or 670 nm as reference then compute total*/
    /*absorption at the reference wavelength, aref*/
    if (rrs_s[idx670] >= 0.0015) {
        aref = aw[idx670] + 0.07 * pow(rrs_a[idx670] / rrs_s[idx443], 1.1);
        idxref = idx670;

    } else {
        numer = rrs_a[idx443] + rrs_a[idx488];
        denom = rrs_a[idx550] + 5 * rrs_a[idx670]*(rrs_a[idx670] / rrs_a[idx488]);
        rho = log10(numer / denom);
        rho = acoefs[0] + acoefs[1] * rho + acoefs[2] * rho*rho;
        aref = aw[idx550] + pow(10.0, rho);
        idxref = idx550;
    }

    /*Calculate the backscattering coefficient at the reference wavelength*/
    bbpref = ((u[idxref] * aref) / (1.0 - u[idxref])) - bbw[idxref];

    /*Steps to compute power exponent of bbp coefficient*/
    rat = rrs_s[idx443] / rrs_s[idx550];
    Ybbp = 2.0 * (1.0 - 1.2 * exp(-0.9 * rat));

    /*Compute backscattering coefficient at 443nm, this is needed later for */
    /*calculating bbp at the Raman excitation bands.s                        */
    bbp443 = bbpref * pow((lambda[idxref] / lambda[idx443]), Ybbp);

    /*Compute the total backscattering coefficient for sensor bands*/
    for (iw = 0; iw < nbandVis; iw++) {
        bbt[iw] = bbpref * pow((lambda[idxref] / lambda[iw]), Ybbp) + bbw[iw];
    }

    /*Calculate the total absorption coefficient at sensor bands*/
    for (iw = 0; iw < nbandVis; iw++) {
        at[iw] = ((1.0 - u[iw]) * bbt[iw]) / u[iw];
    }

    /*Calculate the exponential slope coefficient of absorption by colored */
    /*dissolved and detrital matter */
    Sdg = 0.015 + 0.002 / (0.6 + rrs_s[idx443] / rrs_s[idxref]);

    /*Compute the absorption coefficient of colored dissolved and detrital */
    /*matter at 443 nm, adg443*/
    zeta = 0.74 + 0.2 / (0.8 + rrs_s[idx443] / rrs_s[idxref]);
    xi = exp(Sdg * (lambda[idx443] - lambda[idx412]));
    term1 = (at[idx412] - zeta * at[idx443]) - (aw[idx412] - zeta * aw[idx443]);
    term2 = xi - zeta;
    adg443 = term1 / term2;

    /*Calculate the absorption of phytoplankton at 443 nm*/
    aph443 = at[idx443] - adg443 - aw[idx443];

}

/* -----------------------------------------------------------------------*/
/* raman_iops() - computes IOPs at raman excitation bands using a         */    
/* -----------------------------------------------------------------------*/
void raman_and_sensor_iops() {

    int iw;
    float bbpR, adgR, aphR;
    float bbp, adg, aph;

    /*loop over sensor and raman bands*/
    for (iw = 0; iw < nbandVis; iw++) {

        /*Compute IOPs at Raman bands*/
        bbpR = bbp443 * pow((lambda[idx443] / ramLam[iw]), Ybbp);
        adgR = adg443 * exp(-Sdg * (ramLam[iw] - lambda[idx443]));
        aphR = aph443 * aphStarR[iw];

        /*Compute IOPs at sensor bands*/
        bbp = bbp443 * pow((lambda[idx443] / lambda[iw]), Ybbp);
        adg = adg443 * exp(-Sdg * (lambda[iw] - lambda[idx443]));
        aph = aph443 * aphStar[iw];

        /*total absorption and backscattering coefficients at Raman bands*/
        atotR[iw] = awR[iw] + adgR + aphR;
        bbtotR[iw] = bbwR[iw] + bbpR;

        /*total absorption and backscattering coefficients at sensor bands*/
        atot[iw] = aw[iw] + adg + aph;
        bbtot[iw] = bbw[iw] + bbp;
    }

}

/* -----------------------------------------------------------------------*/
/* k_functions() - computes k functions from IOP bio-optical models       */
/* -----------------------------------------------------------------------*/
void raman_k_func(l2str *l2rec, int ip) {

    int iw;
    float solzRad, solzRadSw, muD, muU;

    //Convert solar and sensor zenith angles from degrees to radianss
    solzRad = l2rec->l1rec->solz[ip]*(M_PI / 180.0);
    solzRadSw = asin(sin(solzRad) / nWater); //subsurface solar zenith angle
    
    //Means cosines
    muD = cos(solzRadSw);
    muU = 0.5;

    for (iw = 0; iw < nbandVis; iw++) {

        //Diffuse attenuation coefficient
        kdSen[iw] = (atot[iw] + bbtot[iw]) / muD;
        kdRam[iw] = (atotR[iw] + bbtotR[iw]) / muD;

        kappaSen[iw] = (atot[iw] + bbtot[iw]) / muU;
        kappaRam[iw] = (atotR[iw] + bbtotR[iw]) / muU;

    }

}

/* -----------------------------------------------------------------------*/
/* radtran_raman_ed() - calculate total downwelling irradiance            */
/* -----------------------------------------------------------------------*/
void raman_radtran_ed(l2str *l2rec, int ip) {

    int32_t iw, ipb;
    
    l1str *l1rec = l2rec->l1rec;

    //Constants
    float p0 = 29.92 * 33.8639; //standard pressure (convert mmHg -> HPa);
    float radCon = M_PI / 180.0; //degrees -> radians conversion
    float ws = l2rec->l1rec->ws[ip]; //wind speed
    float nw2= pow(nWater,2.);  //refractive index of water squared

    //Define model variables
    float sin2Solz, rf0, a0,a1,a2,a3, rho_fres;
    float solz, solzRad, cosSolz, mPres, mAtm, f0Ex;
    float th2oEx, to2Ex, tno2Ex, to3Ex, tgEx;
    float a_285R, a_225R,tau_to200R;
    float no2_tr200R, eDRam_above;

    //Load per-pixel variables from L2 record
    solz = l1rec->solz[ip]; //solar zenith

    //Estimate the transmittance from sun to surface at Raman excitation bands
    fit_raman_tsol(l2rec, ip);

    //Convert solar zenith in degrees to radians
    solzRad = l1rec->solz[ip]*radCon;

    //Calculate cosine of the solar zenith angle
    cosSolz = cos(solzRad);
    
    //Claculate sin2(solz)
    sin2Solz = pow(sin(solzRad),2.);
    
    //Compute Fresnel reflection coefficient for a wind-blown surface following
    //Haltrin 2002
    //Fresnel specular reflection coefficient
    rf0 = 0.5*( pow( (cosSolz - pow(nw2 - sin2Solz ,0.5))/(cosSolz + pow(nw2 + sin2Solz ,0.5)) ,2.) 
            +   pow( (nw2*cosSolz - pow(nw2 - sin2Solz ,0.5))/(nw2*cosSolz + pow(nw2 + sin2Solz ,0.5)) ,2.) );
    
    //Wind-specific coefficients
    a0 = 0.001*(6.944831 - 1.912076*ws  + 0.03654833*ws*ws);
    a1 = 0.7431368 + 0.0679787*ws - 0.0007171*ws*ws;
    a2 = 0.5650262 + 0.0061502*ws - 0.0239810*ws*ws + 0.0010695*ws*ws*ws;
    a3 = -0.4128083 - 0.1271037*ws + 0.0283907*ws*ws - 0.0011706*ws*ws*ws;
    
    //Fresnel reflectance as a function of wind speed
    rho_fres = a0 + rf0*(a1 + rf0*(a2 + a3*rf0));
    

    //---------------------    
    //Calculate the atmospheric pathlengths
    //Use the Kasten abd Young (1989) coefficients
    // mAtm = 1.0 / ( cosSolz + 0.50572*(96.07995 - solz)**(-1.6364) );
    mAtm = 1.0 / (cosSolz + 0.50572 * pow((86.07995 - solz), -1.6364));

    //Greg and Carder (1990) values **NOT USED**
    //mAtm = 1.0 / ( np.cos(solzRad) + 0.15* pow((93.885 - solz),-1.253) );

    //Pathlength corrected for non-standard atmoshperic pressure
    mPres = mAtm * (l1rec->pr[ip] / p0);

    
    //compute no2 above 200m per get_trans.c
    if (l2rec->l1rec->no2_tropo[ip] > 0.0)
    //compute tropo no2 above 200m (Z.Ahmad)    
        no2_tr200R = l1rec->no2_frac[ip] *l1rec->no2_tropo[ip];
    else
        no2_tr200R = 0.0;
    
    //---------------------//
    //Spectral calculations - loop over sensor
    for (iw = 0; iw < nbandVis; iw++) {
        
        ipb =  ip*l1rec->l1file->nbands  + iw;

        //Correct for sun-earth distance
        f0Ex = f0BarRam[iw] * l1rec->fsol; //TOA solar irradiance at sensor bands

        //Compute ozone transmittance
        to3Ex = exp(-(ozRam[iw] * l1rec->oz[ip] / l1rec->csolz[ip]) );
        
        //compute no2 absorption transmittance
        a_285R = no2Ram[iw] * (1.0 - 0.003 * (285.0 - 294.0));
        a_225R = no2Ram[iw] * (1.0 - 0.003 * (225.0 - 294.0));
        tau_to200R = a_285R * no2_tr200R + a_225R * l1rec->no2_strat[ip];
        tno2Ex = exp(-(tau_to200R / l1rec->csolz[ip]) );
        
        //Compute water vapour transmittance
        th2oEx = exp((-0.2385 * wvRam[iw] * l1rec->wv[ip] * mAtm) /
                (pow((1.0 + 20.07 * wvRam[iw] * l1rec->wv[ip] * mAtm), 0.45)));
               

        //Gas Transmittance due to oxygen/gases
        to2Ex = exp((-1.41 * oxyRam[iw] * mPres) / (pow((1.0 + 118.3 * oxyRam[iw] * mPres), 0.45)));
        
        tgEx = to3Ex*tno2Ex*th2oEx*to2Ex;
        
        //      
        eDRam_above = f0Ex * tgEx * t_sol_ram[iw]* cos(solzRad) * 10.0;
        eDSen[iw] = l1rec->Fo[iw] * l1rec->tg_sol[ipb] * l1rec->t_sol[ipb]* cos(solzRad) * 10.0;
        
        eDRam[iw] = 1.03*(1.0 - rho_fres)*eDRam_above;
 
    }

}

/* -------------------------------------------------------------------------*/
/* raman_cor_1() - computes the Rrs corrected for Raman scattering reported */
/*                  by Lee et al, 2013                                      */
/* -------------------------------------------------------------------------*/
void raman_cor_lee1(l2str *l2rec, int ip) {

    int iw;
    float Rrs443, Rrs550;
    float rFactor;
    float rrs_a;
    int nbands = l2rec->l1rec->l1file->nbands;
    //If pixel is already masked, return and do not attempt the correction.
    if (l2rec->l1rec->mask[ip]) {
        return;
    }

    //Get Rrs443 and Rrs550 from the L2 record
    Rrs443 = l2rec->Rrs[ip * nbands + idx443];
    Rrs550 = l2rec->Rrs[ip * nbands + idx550];

    //Loop over vis bands
    for (iw = 0; iw < nbandVis; iw++) {
        
        if (lambda[iw] >= WAVEMIN && lambda[iw] <= WAVEMAX) {

        //Raman correction factor - Eq. 11 in Lee et al 2013
        rFactor = alphaCor[iw] * Rrs443 / Rrs550 + betaCor0[iw] * pow(Rrs550, betaCor1[iw]);

        rrs_a = l2rec->Rrs[ip * nbands + iw];

        Rrs_ram[ip * nbands + iw] = rrs_a - rrs_a / (1.0 + rFactor);
        
        } else { 
            Rrs_ram[ip * nbands + iw] = 0.0;
        }

    }

}

/* -----------------------------------------------------------------------*/
/* raman_cor_westberry() - computes the Rrs corrected for Raman scattering*/
/*                         following the method of Westberry et al 2013   */
/*                                                                        */
/*                         Note: only the first-order scattering term is  */ 
/*                               used as the higher order term often to   */
/*                               resulted in: Rrs_raman > Rrs.            */
/* -----------------------------------------------------------------------*/
void raman_cor_westberry(l2str *l2rec, int ip) {
    
    //We are only considering the first-order Raman scattering

    int iw;
    int nbands = l2rec->l1rec->l1file->nbands;
    float term0, term1, numer1, denom1;
    float Rrs_rc;


    //If pixel is already masked, return and do not attempt the correction.
    if (l2rec->l1rec->mask[ip]) {
        return;
    }

    //Extrapolate aphStar into the UV Raman excitation bands
    set_raman_aph_uv(l2rec, ip);

    //Run QAA to estimate IOP magnitudes and spectral slopes
    raman_qaa(l2rec, ip);

    //Get total absorption and scattering coefficients at both sensor
    //and Raman excitation bands.
    raman_and_sensor_iops();

    //Get the k-functions
    raman_k_func(l2rec, ip);

    //Get the down-welling irradiances
    raman_radtran_ed(l2rec, ip);

    //From Mobley 2010 - assume isotropic emission
    term0 = 1.0 / (4.0 * M_PI * nWater * nWater);
        

    //Loop over visible bands    
    for (iw = 0; iw < nbandVis; iw++) {
        
        //Only calculate correction for visible bands
        if (lambda[iw] >= WAVEMIN && lambda[iw] <= WAVEMAX) {

            //First term in in Westberry et al. 2013, Eq 7.
            numer1 = bRex[iw] * eDRam[iw];
            denom1 = (kdRam[iw]+kappaSen[iw])*eDSen[iw];
            term1 = numer1 / denom1;

            //Raman corrected Rrs (first order scattering term only)
            Rrs_rc = term0 * term1;

            //Fill Rrs_raman array
            Rrs_ram[ip * nbands + iw] = Rrs_rc;
        
        } else {
            Rrs_ram[ip * nbands + iw] = 0.0;
        }

    }

}

/* -----------------------------------------------------------------------*/
/* raman_cor_lee2() - computes the Rrs corrected for Raman scattering     */
/*                          following the method of Lee et al. 1994       */

/* -----------------------------------------------------------------------*/
void raman_cor_lee2(l2str *l2rec, int ip) {

    int iw;
    int nbands = l2rec->l1rec->l1file->nbands;
    float Rrs_rc;

    //If pixel is already masked, return and do not attempt the correction.
    if (l2rec->l1rec->mask[ip]) {
        return;
    }

    //Extrapolate aphStar into the UV Raman excitation bands
    set_raman_aph_uv(l2rec, ip);

    //Run QAA to estimate IOP magnitudes and spectral slopes
    raman_qaa(l2rec, ip);

    //Get total absorption and scattering coefficients at both sensor
    //and Raman excitation bands.
    raman_and_sensor_iops();

    //Get the down-welling irradiances
    raman_radtran_ed(l2rec, ip);

    //Loop over vis bands
    for (iw = 0; iw < nbandVis; iw++) {
        
         if (lambda[iw] >= WAVEMIN && lambda[iw] <= WAVEMAX) {
        //Raman-corrected Rrs 
        Rrs_rc = 0.072 * (bRex[iw] * eDRam[iw]) / ((2.0 * atot[iw] + atotR[iw]) * eDSen[iw]);

        //Fill Rrs_raman array
        Rrs_ram[ip * nbands + iw] = Rrs_rc;
        
         } else {
             Rrs_ram[ip * nbands + iw] = 0.0;
         }

    }

}

/* -----------------------------------------------------------------------*/
/* run_raman_cor() - run the Raman scattering correction algorithm        */

/* -----------------------------------------------------------------------*/
void run_raman_cor(l2str *l2rec, int ip) {

    int iw;
    static int firstCall = 1;
    static int sensorSupported;
    //int32_t ocSensorID;

    //Number of visible bands
    nbandVis = rdsensorinfo(l2rec->l1rec->l1file->sensorID, l1_input->evalmask, "NbandsVIS", NULL);

    //First call 
    if (firstCall) {

        //Supported sensors list
        switch (l2rec->l1rec->l1file->sensorID){
        case MODISA:
        case MODIST:
        case SEAWIFS:
        case VIIRSN:
        case VIIRSJ1:
        case VIIRSJ2:
        case MERIS:
        case OCTS:
        case OLCIS3A:
        case OLCIS3B:
        case SGLI:
        case OCIS:
        case OCI:
            sensorSupported = 1;
            break;
        default:
            //Force l2gen input to raman_opt=0 
            input->raman_opt = NORAMAN;
            //Set flag to 1 - supported sensor
            sensorSupported = 0;
        }
        

        //No Raman correction selected
        if (input->raman_opt == 0) {

            printf("\n");
            if (!sensorSupported) {
                printf("Raman correction unsupported for this sensor. \n");
            }

            printf("No Raman scattering correction calculated for Rrs. \n");
            printf("\n");

            //Raman correction selected
        } else if (input->raman_opt > 0 && input->raman_opt <= 3) {

            printf("\n");
            printf("Compute Raman scattering correction #%d. \n", input->raman_opt);
            printf("\n");

            //Invalid Raman correction selected
        } else {
            printf("-E- %s line %d : '%d' is not a valid Raman correction option.\n",
                    __FILE__, __LINE__, input->raman_opt);
            exit(1);
        }


        /*Allocate memory according to number of sensor bands*/
        raman_pixel_alloc(l2rec);

        //If sensor is supported for Raman correction, coefficients are allocated
        if (sensorSupported) {

            /*Get the sensor-specific coefficients*/
            get_raman_coeffs(l2rec);

            /*Loop over visible bands*/
            for (iw = 0; iw < nbandVis; iw++) {

                //Get the visible band centers
                lambda[iw] = l2rec->l1rec->l1file->fwave[iw];

                //Interpolate the Lee et al. (2013) correction coefficients to sensor resolution
                alphaCor[iw] = linterp(lamCor1, alpha, 6, lambda[iw]);
                betaCor0[iw] = linterp(lamCor1, beta0, 6, lambda[iw]);
                betaCor1[iw] = linterp(lamCor1, beta1, 6, lambda[iw]);

                //Populate the pure water absorption data//
                awR[iw] = aw_spectra(ramLam[iw], ramBandW[iw]);
                bbwR[iw] = bbw_spectra(ramLam[iw], ramBandW[iw]);
                aw[iw] = l2rec->l1rec->l1file->aw[iw];
                bbw[iw] = l2rec->l1rec->l1file->bbw[iw];
                

            }

            //Find the sensor wavelength indices, do this once.
            idx412 = windex(412.0, lambda, nbandVis);
            idx443 = windex(443.0, lambda, nbandVis);
            idx488 = windex(488.0, lambda, nbandVis);
            idx550 = windex(547.0, lambda, nbandVis);
            idx670 = windex(667.0, lambda, nbandVis);
            

            //How many bands between 412 and 490?
            for (iw = 0; iw < nbandVis; iw++) {
                if (lambda[iw] > 550.0) {
                    bandNum550 = iw - 1;
                    break;
                }
            }

        }
        firstCall = 0;
    }
          

    /*No raman correction is applied, raman_opt=0*/
    if (input->raman_opt == NORAMAN) {
        /*Fill Rrs_raman with zeros*/
        for (iw = 0; iw < l2rec->l1rec->l1file->nbands; iw++) {
            Rrs_ram[ip * l2rec->l1rec->l1file->nbands + iw] = 0.0;
        }

        /*Apply the Lee et al 2013 correction, raman_opt=1*/
    } else if (input->raman_opt == LEE2013) {
        raman_cor_lee1(l2rec, ip);

        /*Apply the Westberry et al 2013 correction, raman_opt=2*/
    } else if (input->raman_opt == WESTBERRY2013) {
        raman_cor_westberry(l2rec, ip);

        /*Apply the Lee et al 1994 Raman correction, raman_opt=3*/
    } else if (input->raman_opt == LEE1994) {
        raman_cor_lee2(l2rec, ip);
    }


    l2rec->Rrs_raman = Rrs_ram;

}

/* ---------------------------------------------------------------------------*/
/* ---------------------------------------------------------------------------*/
