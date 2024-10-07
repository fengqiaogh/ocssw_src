/* =================================================================== */
/* module get_bpar.c - compute PAR at depth                            */
/*                                                                     */
/* Experimental code to estimate photosynthetically active radiation   */
/* (PAR) at the seafloor or other user-defined depth. Code assumes     */                  
/* clear sky conditions only.                                          */
/*                                                                     */
/* References:                                                         */
/* Magno-Canto et al. 2019. Model for deriving benthic irradiance in   */ 
/* the Great Barrier Reef from MODIS satellite imagery, Optics Express,*/
/* 27(20), A1350-A1371.                                                */  
/*                                                                     */
/* Implementation:                                                     */
/* L. McKinna  Go2Q,     Jan 2020                                      */
/* M. Canto    JCU/AIMS, Jan 2020                                      */
/*                                                                     */
/* User inputs:                                                        */
/* bpar_validate_opt (default  = 0)                                    */
/*      0: compute using solar zenith at local solar noon              */
/*      1: compute using solar zentih at overpass time                 */
/*                                                                     */
/* bpar_elev_opt (default=0)                                           */
/*      0: use bathymetry data                                         */
/*      1: use scene-wide user-defined depth value                     */
/*                                                                     */
/* bpar_elev_value (default=30.0) {only used when bpar_elev_opt=1}     */
/*      float: input positive value for scene-wide depth value         */
/*                                                                     */
/* Notes:                                                              */
/* Daily integrated values do not consider variability in atmospheric  */
/* conditions and assume a cloud free sky.                             */
/* =================================================================== */

#include <stdio.h>
#include <math.h>
#include "l12_proto.h"

#define KD_MAX  6.4
#define KD_MIN  0.016
        
static float badval = BAD_FLT;
static int32_t LastRecNum = -1;
int32_t timeStepNum = 50;

static float *iparb; //PAR at seafloor
static float *parb; //Daily PAR at seafloor
static double *ed; //above-water down-welling irradiance
static double *edt; //above-water down-welling irradiance
static double *edsubq; //sub-surface down-welling irradiance
static double *kdtot; //diffuse attenuation coefficient
static double *edzq; //quantum irradiance at seafloor in micromol/m2/s
static double *solzArr; //Array of solar zenith angles from sunrise to midday
static double *solTimeArr;
static double *partemp;

static int32_t nbandVis;         //number of visible bands for a given sensor
static double solStepSec;        //integral time step in seconds

/* ------------------------------------------------------------------------*/
/* bpar_ran() - determin in bpar has been run for scanline                 */
/* ------------------------------------------------------------------------*/
int bpar_ran(int recnum)
{ 
    if ( recnum == LastRecNum )
        return 1;
    else
        return 0; 
}

/* ---------------------------------------------------------------------- */
/* bpar_init()  - allocates memory for benthic par algorithm             */
/* ---------------------------------------------------------------------- */
void bpar_init(l2str *l2rec) {
    
    l1str *l1rec = l2rec->l1rec;
    int32_t nbands = l1rec->l1file->nbands;
    
    ed = (double*) allocateMemory(l1rec->npix *nbands * sizeof(double), "ed"); 
    edsubq = (double*) allocateMemory(l1rec->npix *nbands * sizeof(double), "edsubq");
    kdtot = (double*) allocateMemory(l1rec->npix *nbands * sizeof(double), "ktot");
    edzq = (double*) allocateMemory(l1rec->npix *nbands * sizeof(double),  "edzq");
    iparb = (float*) allocateMemory(l1rec->npix * sizeof(float), "iparb");
    
    edt = (double*) allocateMemory( timeStepNum*nbands * sizeof(double), "edt");
    solTimeArr = (double*) allocateMemory(timeStepNum * sizeof(double), "solTimeArr");
    solzArr = (double*) allocateMemory(timeStepNum* sizeof(double), "solzArr");
    partemp = (double*) allocateMemory(timeStepNum * sizeof(double), "partemp");
    parb = (float*) allocateMemory(l1rec->npix * sizeof(float), "iparb");
 
}

/* ---------------------------------------------------------------------- */
/* calc_isleap() - calculates if year is leap or not                      */
/* ---------------------------------------------------------------------- */
int cacl_isleap(int year) {
    
    int isleap;
    
    //test if it is a leap year
    if (year%4  == 0 ) {
        if (year%100 == 0) { 
            if (year%400 == 0) {
                isleap = 1;
            } else {
                isleap = 0;
            }
         } else {
            isleap= 1;
         }
    } else {
        isleap = 0;
    }
    
    return isleap;
}

/* ---------------------------------------------------------------------- */
/* calc_solz_t() - calculates solar zenith angle for given solar time     */
/* ---------------------------------------------------------------------- */
double calc_solz_t(l2str *l2rec, int ip, double soltime) {
    
    l1str *l1rec = l2rec->l1rec;
    float lat_rad = l1rec->lat[ip] * M_PI/180.;  //latitude in radians
    
    float date_angle;       //date angle
    float year_len;         //number of days in a year
    double soldecl;         //solar declination angle
    double solz_t;          //solar zenith angle for given solar time
    
    int16_t year, day;
    double sec;
    //Get year, day, and second from scantime
    unix2yds(l2rec->l1rec->scantime, &year, &day, &sec);
    
        
    if (cacl_isleap(year)){
        year_len=366.;
    } else {
        year_len=365.;
    }

    //compute date angle for overpass
    date_angle = 2.*M_PI * (day-1) / year_len;
       
    //Compute solar declination angle
    soldecl =  0.006918 - 0.399912*cos(date_angle) + 0.070257*sin(date_angle)
            - 0.006758*cos(2.*date_angle) + 0.000907*sin(2.*date_angle)
            -0.002697*cos(3.*date_angle) +0.00148*sin(3.*date_angle);
    
    solz_t = acos( sin(lat_rad)*sin(soldecl) + cos(lat_rad)*cos(soldecl)*cos(soltime) );
    
    return solz_t;
            
}

/* ------------------------------------------------------------------------ */
/* calc_solz_t_array() - calculates solar zenith angle for given solar time */
/* ------------------------------------------------------------------------ */
void calc_solz_t_array(l2str *l2rec, int ip) {
    
    l1str *l1rec = l2rec->l1rec;
    float lat_rad = l1rec->lat[ip] * M_PI/180.;  //latitude in radians
    
    int ix;
    float date_angle;       //date angle
    float year_len;         //number of days in a year
    double soldecl;         //solar declination angle
    double sunrise;
    double solTimeAngleStep;
    double solTimeTemp;
    int16_t year, day;
    double sec;
    
    //Get year, day, and second from scantime
    unix2yds(l2rec->l1rec->scantime, &year, &day, &sec);
        
    if (cacl_isleap(year)){
        year_len=366.;
    } else {
        year_len=365.;
    }

    //compute date angle for overpass
    date_angle = 2.*M_PI * (day-1.) / year_len;
       
    //Compute solar declination angle
    soldecl =  0.006918 - 0.399912*cos(date_angle) + 0.070257*sin(date_angle)
            - 0.006758*cos(2.*date_angle) + 0.000907*sin(2.*date_angle)
            -0.002697*cos(3.*date_angle) +0.00148*sin(3.*date_angle);
    
    //compute sunrise time angle 
    sunrise = -1.*acos(-1.*tan(soldecl)*tan(lat_rad)); //convention midday at 0
    
    //Create angle array from sunrise to solar noon
    solTimeAngleStep = (0.0 - sunrise)/(timeStepNum) ;
    
    //What is the solar angle step in seconds, this is needed for daily PAR integration
    solStepSec = (24.*60.*60.)*solTimeAngleStep/(2.*M_PI);
    
    //Create an array of solar time angles from sunrise to midday
    solTimeArr[0] = 0.0;
    solzArr[0]= acos( sin(lat_rad)*sin(soldecl) + cos(lat_rad)*cos(soldecl)*cos(solTimeArr[0]) );

    for (ix = 1; ix < timeStepNum; ix++) {
        
        solTimeTemp = solTimeArr[ix-1] - solTimeAngleStep;
        solTimeArr[ix] = solTimeTemp;
        
        //Calculate solar zenith angle from sunrise to solar noon
        solzArr[ix]= acos( sin(lat_rad)*sin(soldecl) + cos(lat_rad)*cos(soldecl)*cos(solTimeArr[ix]) );
        
    }
}

/* ---------------------------------------------------------------------- */
/* calc_benthic_irrad() - calculates instantaneous benthic PAR            */
/* ---------------------------------------------------------------------- */
double calc_benthic_ipar(l2str *l2rec, int ip, double z, double solza) {
    
    l1str *l1rec = l2rec->l1rec;
    int32_t nbands = l1rec->l1file->nbands;

    float *wave = l1rec->l1file->fwave;  //wavelength array
    
    int qcFlag = -1;
    
    float iparbtemp;

    int iw, ipb, ipb1, ipb2;
    float rho_dir;          //direct beam reflectance  
    float rho_dif = 0.066;  //diffuse reflectance for uniform overcast sky    
    float nw = 1.341;       //refractive index of water
    float solzw;            //Sub-surface solar zenith angle
    float trans;            //transmittance coefficient
    float cloud = 0;        //proportion of cloud cover

    //Values taken from Reinart et al (1998), JGR, 103(C4), 7749-7752
    double h = 6.6255E-34;         //Planck's constant J s
    double c = 2.9979E8;           //Speed of light m /s
    //double umole = 6.022E17;        //micromole photons / (m^2 s)
    double mole = 6.022E23;        //mole photons / (m^2 s)
    double nmtom = 1E-9;            //convert nm to m
    double quanta = nmtom/(mole*h*c); //quanta energy conversion factor umol / (s W nm)
    
    
    //compute reflectance of the direct (solar) beam
    if ( solza != 0.0 ) {
        solzw = asin(sin(solza)/nw);
        rho_dir = 0.5* ( pow( (sin(solza - solzw)/sin(solza + solzw)), 2) + 
                         pow( (tan(solza - solzw)/tan(solza + solzw)), 2)  );
    } else {
        rho_dir = pow( ((nw - 1.)/(nw + 1.)), 2);
    }

    //Assume a clear, non-overcast sky & calculate transmittance coefficient.
    //However, I have included diffuse reflectance in case we modify this in
    //future.
    trans = (1.0 - rho_dir)*( 1.0 - cloud ) + (1.0 - rho_dif)*cloud;


    for (iw = 0; iw < nbandVis; iw++) {

        ipb = ip*nbands+iw;

        //Compute incident above-water down-welling irradiance W /m^2/nm
        //**Normalize by cos(solz) for Ed where sun is directly overhead**
        //Note 1: the factor of 0.01 converts from uW/cm^2/nm to W/m^2/nm
        //Note 2: to scale to Ed at specific time, you could multiply by cos(solz)
        //If sun is directly overhead, cos(0) = 1.
        ed[ipb] = l1rec->Fo[iw]* l1rec->tg_sol[ipb] * l1rec->t_sol[ipb] *cos(solza)* 0.01 ;

        //Compute sub-surface down-welling irradiance
        //Convert to photon flux units of umol photos / (m^2  s)
        edsubq[ipb] = quanta*wave[iw]*trans*ed[ipb];
        
        //Compute spectral diffuse attenuation coefficient following 
        //Lee et al. 2005
        if (l1rec->mask[ip] || l2rec->a[ipb] <= 0.0 || l2rec->bb[ipb] <= 0.0 ) {
            kdtot[ipb] = badval;
            edzq[ipb] = badval;

        } else {
            //Compute diffuse attenuation coefficient using Lee et al (2005)
            kdtot[ipb] = (1.0 + 0.005 * (solza*180./M_PI)) 
                        * l2rec->a[ipb] + 4.18 
                            * (1.0 - 0.52* exp(-10.18* l2rec->a[ipb]))
                                * l2rec->bb[ipb];

            //Check range of kd values
            if (kdtot[ipb] > KD_MAX) {
                kdtot[ipb] = KD_MAX;
                l1rec->flags[ip] |= PRODWARN;
            } else
            if (kdtot[ipb] < KD_MIN) {
                kdtot[ipb] = KD_MIN;
                l1rec->flags[ip] |= PRODWARN;
            }

            //Spectral photon flux the seafloor in unit of mol photon / (m^2 um)
            edzq[ipb] = edsubq[ipb]*exp(-z*kdtot[ipb]);

        }
    }
        
    //Integrate over band using trapezoid rule
    iparbtemp = 0.0;
     
    for (iw = 0; iw < nbandVis-1 ; iw++) {
        
        ipb1 = ip*nbands+iw;
        ipb2 = ip*nbands + (iw+1);
        
        if (l1rec->mask[ip] || edzq[ipb1] < 0 || edzq[ipb2] < 0 ) {
            qcFlag = 1;
            l1rec->flags[ip] |= PRODFAIL;
            break;
            
        } else {
           qcFlag = 0;
           iparbtemp += 0.5*(wave[iw+1] - wave[iw])*(edzq[ipb2]+edzq[ipb1]);     
        }            
    }
    
    if (qcFlag) {
        return BAD_FLT;
    } else {
        return iparbtemp;
        
    }
   
}

/* ---------------------------------------------------------------------- */
/* calc_benthic_par() - calculates daily benthic PAR                      */
/* ---------------------------------------------------------------------- */
double calc_benthic_par(l2str *l2rec, int ip, double z) {
    
    l1str *l1rec = l2rec->l1rec;
    
    double parbtemp;
    double ipartemp;
    double solztemp;
    int qcFlag ,ix;
 
    //calculate solar zenith array
    calc_solz_t_array(l2rec, ip);
    
    //Calculate instantaneous benthic par from sunrise to solar noon
    for (ix = 0; ix < timeStepNum; ix++) {
        solztemp = solzArr[ix];
        ipartemp =  calc_benthic_ipar(l2rec, ip, z, solztemp);
        partemp[ix]= ipartemp;
    }

    qcFlag = 0; 
    parbtemp = 0.0;
    for (ix = 0; ix < timeStepNum-1 ; ix++) {   
        
        //Check for bad values and flags
        if (l1rec->mask[ip] || partemp[ix] < 0 ) {
            qcFlag = 1;
            l1rec->flags[ip] |= PRODFAIL;
            break;
            
        } else { 
            qcFlag = 0;
            parbtemp += 0.5*solStepSec*(partemp[ix+1]+partemp[ix]);  
        }   
    }
    
    if (qcFlag) {
        return BAD_FLT;        
    } else {
        //Multiply by two as we integrated over half the day only.  
        //We assume the solar elevation is symmetric about noon.
        return 2.*parbtemp;
    }
}

/* ---------------------------------------------------------------------- */
/* run_bpar() - calculates run benthic par algorithm                      */
/* ---------------------------------------------------------------------- */
void run_bpar(l2str *l2rec) {
    
    l1str *l1rec = l2rec->l1rec;
    filehandle *l1file = l1rec->l1file;
    static int32_t bparScanNum = 1;
    
    double depth;   //water-column depth for a given pixel
    double solz_iparb; //zenith angle at pixel's solar noon
    int ip;
    

    if ( bparScanNum) {

        nbandVis = rdsensorinfo(l1file->sensorID, l1_input->evalmask,
                "NbandsVIS", NULL);

        /*Allocate memory*/
        bpar_init(l2rec);
        
        if (input->iop_opt == 0) {
            printf("\n");
            printf("bPAR model requires iop model selection. Using SWIM model (iop_opt=8).\n");
            input->iop_opt=IOPSWIM;
        }
       
        if (input->bpar_validate_opt == 0) {
            printf("bPAR model validation option off. Compute iparb for solar noon. \n");
        } else if (input->bpar_validate_opt == 1)  {
            printf("bPAR model option validation on. Compute ipbar using overpass solar zenith.\n");
        } else {
            fprintf(stderr,
                "-E- %s line %d: invalid input value for bpar_validate_opt: \"%d\".\n",
                __FILE__, __LINE__,input->bpar_validate_opt);
            exit(1);
        }
        
        if (input->bpar_elev_opt == 0) {
            printf("bPAR model using bathymetry raster file. \n\n");
        } if (input->bpar_elev_opt == 1)  {
            printf("bPAR model using user-supplied value of %f m for water depth. \n\n", input->bpar_elev_value);
        } else if (input->bpar_elev_opt !=0 && input->bpar_elev_opt !=1) {
            fprintf(stderr,
                "-E- %s line %d: invalid input value for bpar_elev_opt: \"%d\".\n\n",
                __FILE__, __LINE__,input->bpar_elev_opt);
            exit(1);
        }
        bparScanNum = 0;
        
    }
        
    //Get iops - recommend use of shallow water IOP model for opticaly shallow
    //waters. SWIM is iop_opt=8
    get_iops(l2rec,input->iop_opt);
        
    /*iterate over pixels*/
    for (ip = 0; ip < l1rec->npix; ip++) {
        

        if (input->bpar_elev_opt==0) {
            depth = 0.0 - l1rec->dem[ip];
        } else if (input->bpar_elev_opt==1) {
            depth = input->bpar_elev_value;
        }

        iparb[ip]  = BAD_FLT;
        parb[ip]  = BAD_FLT;

        //Check if we are on the continental shelf (300m conservative)
        //Ignore deep water pixels
        if (depth > 300.) {

            //If water deeper than 300m set as bad value
            iparb[ip] = BAD_FLT;;
            parb[ip] = BAD_FLT;

        } else {
            
            if (input->bpar_validate_opt==0) {
                solz_iparb = calc_solz_t(l2rec, ip, 0.0);
            }
            
            if (input->bpar_validate_opt==1) {
                solz_iparb =  l1rec->solz[ip] * M_PI/180;
            }
            
            iparb[ip]  = calc_benthic_ipar(l2rec, ip, depth, solz_iparb);
            parb[ip] = calc_benthic_par(l2rec,ip, depth);

        }

    }
    
    LastRecNum = l1rec->iscan;
      
}

/* -------------------------------------------------------------------- ------*/
/* get_ibpar() - returns requested benthic irradiance product for full scanline*/
/* ---------------------------------------------------------------------------*/
void get_bpar(l2str *l2rec,  l2prodstr *p, float prod[]) { 
    
    l1str *l1rec = l2rec->l1rec;
    int ip;
    
    if (!bpar_ran(l1rec->iscan))
        run_bpar(l2rec);
    
    for (ip = 0; ip < l2rec->l1rec->npix; ip++) {
        
        switch (p->cat_ix) {
        case CAT_iparb:
            prod[ip] = (float) iparb[ip];
            break;
         case CAT_parb:
            prod[ip] = (float) parb[ip];
            break;
        default:
            printf("Error: %s : Unknown product specifier: %d\n", __FILE__, p->cat_ix);
            exit(FATAL_ERROR);
            break;
        }
    }

}
