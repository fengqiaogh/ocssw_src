/*=========================================================================*/
/*                             get_sdp.c                                   */
/*                                                                         */
/* Description:                                                            */ 
/* This algorithm calculates multiple phytoplankton pigment concentrations */
/* based on spectral derivative of hyperspectral radiometry                */
/*                                                                         */
/*                                                                         */
/* Notes:                                                                  */
/*                                                                         */
/* References:                                                             */
/* Kramer, S., Siegel, D.A., Maritorena, S. and D. Catlett. (2022) Modeling*/
/* surface ocean phytoplankton pigments from hyperspectral remote sensing  */
/* reflectance on global scales, Remote Sens. Environ., 112879             */
/* doi: 10.1016/j.rse.2021.112879.                                         */
/*                                                                         */
/* Implementation:                                                         */
/* L. McKinna (Go2Q), 19 Dec 2022                                          */
/*=========================================================================*/

#include <math.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_filter.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_statistics.h>
#include <levmar.h>
#include "l12_proto.h"
#include "get_sdp.h"

static int32_t LastRecNum = -1;

static int npig = 13;       //Number of pigments in model
static int nAcoeff = 299;   //Number of model coefficents
static int nreps = 100;     //Number of replicates in model coefficients
static char *pignames[13] = {"tchla","zea","dvchla","butfuco","hexfuco","allo",
                        "mvchlb","neo","viola","fuco","chlc12","chlc3","perid"};


static int16  *sdpflags;   //Binary quality control flags  
static float *pigArr;  //Pigment array pointer (this is where we'll store the solutions)
static float *pigArr_unc; //Pigment uncertinty array pointer (this is where we'll store the 
                       //uncertainty in solutions)
static float *mRrs_out; //Array holding forward-modelled Rrs
static float *Rrs_diff_out; //Array holdig difference between modeled and observed Rrs.
static float *d2Rrs_out; //Array holding 2nd Derivative of RrsDiff


//static int nfitbands = 6;
//static float fitbands[6] = {412,443,490,530,555,670};

//static int nfitbands = 12;
//static float fitbands[12] = {400, 425, 450, 475, 500, 525, 550, 575, 600, 625, 650, 675};

static int nfitbands = 56;
static float fitbands[56] = {400,405,410,415,420,425,430,435,440,445,450,455,460,465,470,475,
                             480,485,490,495,500,505,510,515,520,525,530,535,540,545,550,555,
                             560,565,570,575,580,585,590,595,600,605,610,615,620,625,630,635,
                             640,645,650,655,660,665,670,675};


void sdp_qssa( double *p, double *mrrs_s, int m, int n, void *data);

/* ------------------------------------------------------------------- */
/*sdp_ran() - record if run_sdp has been executed for scan line       */

/* ------------------------------------------------------------------- */
int sdp_ran(int recnum) {
    if (recnum == LastRecNum)
        return 1;
    else
        return 0;
}

/* -----------------------------------------------------------------*/
/* aph_kramer() - get A and B model coefficients for modelling aph  */
/*                 as a function of chl                             */

/****************Note: may moved into aph.c after testing************/

/* -----------------------------------------------------------------*/
int aph_kramer_2022(float wave, float *A, float *B) {
    
    static int firstCall = 1;
    
    static float *table;
    static float *wtab;
    static float *Atab, *Btab;
    static float *dAtab, *dBtab;
    static int ntab = 0;

    
    if (firstCall) {

        char *filedir;
        char filename[FILENAME_MAX];
        int ncol, nrow;

        if ((filedir = getenv("OCDATAROOT")) == NULL) {
            printf("-E- %s: OCDATAROOT env variable undefined.\n", __FILE__);
            exit(1);
        }
        strcpy(filename, filedir);
        strcat(filename, "/common/aph_kramer_2022.txt");
        printf("Reading aph* for sdp model from %s.\n", filename);

        nrow = table_row_count(filename);
        
        if (nrow <= 0) {
            printf("-E- %s line %d: error opening (%s) file", __FILE__, __LINE__, filename);
            exit(1);
        }

        ncol = table_column_count(filename);
        table = table_read_r4(filename, ncol, nrow);
        ntab = nrow;

        wtab = &table[nrow * 0]; //wavelength 
        Atab = &table[nrow * 1]; //Spectral A coefficient
        Btab = &table[nrow * 2]; //Spectral B coefficent

        // precompute derivatives for spline interpolation
        if ((dAtab = (float*) calloc(ntab, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for Kramer aph.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((dBtab = (float*) calloc(ntab, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for Kramer aph.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        spline(wtab, Atab, ntab, 1e30, 1e30, dAtab);
        spline(wtab, Btab, ntab, 1e30, 1e30, dBtab);

        firstCall = 0;
    }

    if (wave > wtab[ntab - 1]){
        wave = wtab[ntab - 1];
    }

    // interpolate coefficients to wavelength to get A and B
    splint(wtab, Atab, dAtab, ntab, wave, A);
    splint(wtab, Btab, dBtab, ntab, wave, B);
    
    return 0;

}

/* ------------------------------------------------------------------*/
/* get_sdp_aw() - read netcf file and return purewater absorption    */
/*                    coefficients                                   */

int get_sdp_aw(float wave, float *A)  {

    static int firstCall = 1;

    int32_t status, netcdf_input;
    size_t nlam;
    int ncid, grpid, varid1, varid2, dimid ;

    static float *wtab;
    static float *atab;
    static float *dAtab;
    static int ntab = 0;

    if (firstCall) {

        char *filedir;
        char filename[FILENAME_MAX];      
        
        if ((filedir = getenv("OCDATAROOT")) == NULL) {
            printf("-E- %s: OCDATAROOT env variable undefined.\n", __FILE__);
            exit(1);
        }
        
        strcpy(filename, filedir);
        strcat(filename, "/common/sdp_pure_water.nc");
        printf("Reading SDP pure water coefficients from %s.\n", filename);
        
        /* Does the file exist? */
        if (access(filename, F_OK) || access(filename, R_OK)) {
            printf("-E- %s: SDP pure water coefficient file '%s' does not exist or cannot open.\n",
                    __FILE__, filename);
            exit(EXIT_FAILURE);
        }
        
        /* test for NetCDF input file */
        status = nc_open(filename, NC_NOWRITE, &ncid);
        netcdf_input = (status == NC_NOERR);

        if (netcdf_input) {
            
            /* Get the group ids of our two groups. */
            if ((nc_inq_ncid(ncid, "water_absorption", &grpid)) == NC_NOERR) {

                //Get the dimensions of the variable
                nc_inq_dimid(grpid,"n_wavelength", &dimid);
                nc_inq_dimlen(grpid, dimid, &nlam);
                ntab = (int) nlam;

                //Allocate memory
                atab = (float*) allocateMemory(ntab * sizeof (float), "atab");
                wtab = (float*) allocateMemory(ntab * sizeof (float), "wtab");
                dAtab = (float*) allocateMemory(ntab * sizeof (float), "dAtab");
                
                nc_inq_varid(grpid, "aw", &varid1);
                nc_get_var_float(grpid, varid1, atab);

                nc_inq_varid(grpid, "wavelength", &varid2);
                nc_get_var_float(grpid, varid2, wtab);                
            
            }
        

        spline(wtab, atab, nlam, 1e30, 1e30, dAtab);

        firstCall = 0;

        }
    }

    if (wave > wtab[ntab - 1]){
        wave = wtab[ntab - 1];
    }

    // interpolate coefficients to wavelength to get A
    splint(wtab, atab, dAtab, ntab, wave, A);

    return 0;
    
}


/* -----------------------------------------------------------------*/
/* get_sdp_coeff() - read netcf file and return the A and C model   */
/*                    coefficients                                  */
/* input:                                                           */
/*       pigment name (char)                                        */
/* output:                                                          */
/*       A (float array)                                            */
/*       C (float)                                                  */

/* -----------------------------------------------------------------*/
int get_sdp_coeff(double ** A, double **C, char **pigments, int n, int m) {
    
    static int firstCall = 1;
    
    int32_t i, status, netcdf_input;
    int ncid, varidA, varidC, grpidA, grpidC;
    
    char varNameA[16]; 
    char varNameC[16];
    
    if (firstCall) {
        
        char *filedir;
        char filename[FILENAME_MAX];      
        
        if ((filedir = getenv("OCDATAROOT")) == NULL) {
            printf("-E- %s: OCDATAROOT env variable undefined.\n", __FILE__);
            exit(1);
        }
        
        strcpy(filename, filedir);
        strcat(filename, "/common/sdp_AC_coeffs_all.nc");
        printf("Reading SDP model coefficients from %s.\n", filename);
        
        /* Does the file exist? */
        if (access(filename, F_OK) || access(filename, R_OK)) {
            printf("-E- %s: SDP coefficient file '%s' does not exist or cannot open.\n",
                    __FILE__, filename);
            exit(EXIT_FAILURE);
        }
        
        /* test for NetCDF input file */
        status = nc_open(filename, NC_NOWRITE, &ncid);
        netcdf_input = (status == NC_NOERR);
        
        if (netcdf_input) {
            
            for (i=0; i<npig; i++) {
                
                strncpy(varNameA, "A",16);
                strncat(varNameA,pigments[i],16);
                strncpy(varNameC, "C",16);
                strncat(varNameC,pigments[i],16);
            
                /* Get the group ids of our two groups. */
                if ((nc_inq_ncid(ncid, "A_coeff", &grpidA)) == NC_NOERR) {
                    if ((nc_inq_varid(grpidA, varNameA, &varidA)) == NC_NOERR) {
                        /* Read the coefficeint  data. */
                        nc_get_var_double(grpidA, varidA, A[i]);
                    }
                }
                
                if (( nc_inq_ncid(ncid, "C_coeff", &grpidC)) == NC_NOERR) {	
                    if ((nc_inq_varid(grpidC, varNameC, &varidC)) == NC_NOERR) {
                        nc_get_var_double(grpidC, varidC, C[i]);
                    }        
                }
            }
        }
        
        firstCall = 0;
    }

    return 0;  
}


/* ------------------------------------------------------------------- */
/* sdp_interp() - linear 1D interpolation function                     */

/* Inputs                                                               */
/*       x: initial x array                                            */
/*       y: intial y array                                              */
/*       xinterp: new x array                                           */
/*       n: length or x and y                                           */
/*       m: length of xinterp and yinterp                               */
/*                                                                      */
/* Output                                                               */
/*       yinterp: new y array                                           */

/* ------------------------------------------------------------------- */

int sdp_interp(float *x, float *y, float *xinterp, float *yinterp,
                int m, int n) {
        int iw;
        double *xd;
        double *yd;
        
        xd = (double*)calloc(n, sizeof(double));
        yd = (double*)calloc(n, sizeof(double));
        
        gsl_interp *interp_work;
        
        for (iw=0; iw < n; iw++) {
            xd[iw] = x[iw];
            yd[iw] = y[iw];
        }
        
        
        interp_work = gsl_interp_alloc(gsl_interp_linear, m);
        gsl_interp_init(interp_work,  xd, yd, m);

        for (iw=0; iw < n; iw++) {
            yinterp[iw] = (float) gsl_interp_eval(interp_work, xd, yd, (double) xinterp[iw], NULL);
        }

    gsl_interp_free(interp_work);
    free(xd);
    free(yd);
    return 0;
}

/* ------------------------------------------------------------------- */
/* sdp_mean_filter() - simple mean moving average filter               */

/* ------------------------------------------------------------------- */

int sdp_mean_filter(float *y, float *yfilt, int n, int window) {

    int iw, ih, halfwin;
    float temp;

    if(window % 2 == 0) {
            printf("-E- %s line %d : mean moving average filter window must be odd value.\n",
            __FILE__, __LINE__);
            exit(1);
    }

    halfwin = (window-1)/2;

    for (iw=0; iw< n; iw++) {
        
        if (iw <= halfwin) {
            temp = y[iw];
        } else if (iw >= halfwin) {
            temp = y[iw];
        } else {
            for (ih=0;ih<window;ih++) {
                temp += y[iw+ih];
            }
            temp /= window;
        }
        yfilt[iw] = temp;
    }
    return 0;
}

/* ------------------------------------------------------------------- */
/* sdp_filter() - simple median filter                                 */

/* ------------------------------------------------------------------- */
int sdp_filter(float *y, float *yfilt, int n, int window) {
    
    int iw;
    int success;
    float temp;
    
    gsl_filter_median_workspace *medfilt;
    medfilt = gsl_filter_median_alloc(window);
    
    gsl_vector *yin = gsl_vector_alloc(n);
    gsl_vector *yout = gsl_vector_alloc(n);
    
    for (iw=0; iw<n; iw++) {
        gsl_vector_set(yin, iw, y[iw]);
        gsl_vector_set(yout, iw, 0.0);
    }
    
    success= gsl_filter_median(GSL_FILTER_END_PADVALUE, yin, yout, medfilt);
    
   for (iw=0; iw<n; iw++) {
        temp = gsl_vector_get(yout, iw) ;
        yfilt[iw] = temp;
    }
    
    //cleanup
    gsl_filter_median_free(medfilt);
    gsl_vector_free(yin);
    gsl_vector_free(yout);
    return success;
}

/* ------------------------------------------------------------------- */
/* sdp_deriv2() - numerical second derivative                          */

/* ------------------------------------------------------------------- */

int sdp_deriv2(float *x, float *y, float *deriv2, int n) {
    
    int i;
    float d1, d2;
    
    //zero pad first and last element of output
    deriv2[0] =0.0;
    deriv2[n] = 0.0;
    
    for (i=1; i<=n-1; i++) {
        
        //1st derative at two points
        d1 = (y[i] - y[i-1]) /(x[i]-x[i-1]);   //1st Deriv at point 1
        d2 = (y[i+1] - y[i])/(x[i+1]-x[i]);    //1st Deriv at point 2
        
        //Compute 2nd deriv from d1 and d2
        deriv2[i] = (d2 - d1)/ ((x[i+1] - x[i-1])/2.);


    }
    
    return 0;
}

/* ------------------------------------------------------------------- */


/* ------------------------------------------------------------------- */
/* sdp_qssa() - Quasi-single scattering approximation (QSSA) forward   */
/*              reflectance model                                      */

/* ------------------------------------------------------------------- */


void sdp_qssa( double *p, double *mrrs_s, int m, int n, void *data){
    
    sdpstr *s = ((sdpstr *)data);
    
    int iw;
    double atot, bbtot, u;
    double aph, adg, bbp;
    
    for (iw = 0; iw < n; iw++) {

        //Calculate the total IOPs 
        aph = s->aphA[iw]*pow(p[0],s->aphB[iw]);
        adg = p[1]*s->adgstar[iw];
        bbp = p[2]*s->bbpstar[iw];
        atot = s->aw[iw] + aph + adg;
        bbtot = s->bbw[iw] + bbp;

        
        u = bbtot/(atot + bbtot);
        mrrs_s[iw] = (s->g[0] + s->g[1]*u)*u; 

    } 
} 

/* ------------------------------------------------------------------- */
/* sdp_qssaf() - Quasi-single scattering approximation (QSSA) forward   */
/*              reflectance model at fine resolution 1nm                */

/* ------------------------------------------------------------------- */


void sdp_qssaf( double *p, double *mrrs_s, int m, int n, void *data){
    
    sdpstr *s = ((sdpstr *)data);
    
    int iw;
    double atot, bbtot, u;
    double aph, adg, bbp;
    
    for (iw = 0; iw < n; iw++) {

        //Calculate the total IOPs 
        aph = s->faphA[iw]*pow(p[0],s->faphB[iw]);
        adg = p[1]*s->fadgstar[iw];
        bbp = p[2]*s->fbbpstar[iw];
        atot = s->faw[iw] + aph + adg;
        bbtot = s->fbbw[iw] + bbp;
        
        u = bbtot/(atot + bbtot);
        mrrs_s[iw] = (s->g[0] + s->g[1]*u)*u; 
    } 
} 

/* ------------------------------------------------------------------- */
/* cal_adstar() - compute shape of adgstar                             */ 


/*two options, 0 for spectral fit model, 1 for 1nm forward model       */
/* ------------------------------------------------------------------- */
void calc_adgstar(sdpstr *s, int opt) {
    
    int iw;
   
    switch (opt) {
    case 0:
        for (iw = 0; iw < s->nfitbands; iw++) {
            /*Colored dissolved and detrital matter absorption shape*/       
            s->adgstar[iw] = exp((s->bands[iw] - s->bands[s->i440])*s->adg_s);
        }
        break;
    case 1:
        for (iw = 0; iw < s->nfbands; iw++) {
            /*Colored dissolved and detrital matter absorption shape*/       
            s->fadgstar[iw] = exp((s->fbands[iw] - s->fbands[s->i440f])*s->adg_s);
        }
        break;
    default:
        printf(
                "-E- %s, %d: Spectral Derivative Pigment adgstar "
                "model invalid input option value\n",
                __FILE__, __LINE__);
        exit(1);
        break;
    }
}
        

/* ------------------------------------------------------------------- */
/* cal_bbpstar() - compute shape of adgstar                             */

/*two options, 0 for spectral fit model, 1 for 1nm forward model       */
/* ------------------------------------------------------------------- */
void calc_bbpstar(sdpstr *s, int opt) {
    
    int iw;
   
    switch (opt) {
    case 0:
        for (iw = 0; iw < s->nfitbands; iw++) {
            /*Particle backscattering shape*/       
            s->bbpstar[iw] = pow((s->bands[iw]/s->bands[s->i440] ),s->bbp_s);
        }
        break;
    case 1:
        for (iw = 0; iw < s->nfbands; iw++) {
            /*Particle backscattering shape*/       
            s->fbbpstar[iw] = pow(( s->fbands[iw]/s->fbands[s->i440f] ),s->bbp_s);
            }
        break;
    default:
        printf(
                "-E- %s, %d: Spectral Derivative Pigment bbpstar "
                "model invalid input option value\n",
                __FILE__, __LINE__);
        exit(1);
        break;
    }
}


/* ------------------------------------------------------------------- */
/* sdp_jacqssa() - Compute Jacobian of cost function                   */
/*                 e    = y - y(x)                                     */
/*                 dedx =   - dy(x)/dx                                 */

/* dimension of jac is: npar x nbands                                  */

/* ------------------------------------------------------------------- */
void sdp_jacqssa(double *p, double *jac, int n, int m, void* data) {
    
    sdpstr *s = ((sdpstr *)data);

    int i, j;
    double atot, bbtot, aph, adg, bbp, u;
    double A, B;
    double drdchl;
    double drdadg;
    double drdbbp;
    
    for (i=j=0; i< m; i++) {
        
        aph = s->aphA[i]*pow(p[0], s->aphB[i]);
        adg = p[1] * s->adgstar[i];
        bbp = p[2] * s->bbpstar[i];
        atot = s->aw[i] + aph + adg;
        bbtot = s->bbw[i] + bbp;
        
        u = bbtot/(atot + bbtot);
        A = s->aphA[i];
        B = s->aphB[i];
            
        drdchl =  -(s->g[0] + 2.*s->g[1]*u)*(bbtot*A*B*pow(p[0],(B-1.))) / pow( bbtot + atot,2.);
        drdadg =  -(s->g[0] + 2.*s->g[1]*u)*(bbtot*s->adgstar[i])/ pow( bbtot + atot, 2.);
        drdbbp =   (s->g[0] + 2.*s->g[1]*u)*(atot*s->bbpstar[i]) / pow( bbtot + atot, 2.);
        
        jac[j++] = drdchl;
        jac[j++] = drdadg;
        jac[j++] = drdbbp;     
 
    }
}


/* ------------------------------------------------------------------- */
/* rund sdp_saa() - runs sdp semianalytical algorithm to get to mRrs    */

/* ------------------------------------------------------------------- */
int run_sdp_saa(sdpstr *s) {
  
    int iw;
    int statusLM;
    int maxit = 500;

    double p[3] = {0.01, 0.1, 0.0005};
    double info[LM_INFO_SZ]; 
    double opts[LM_OPTS_SZ];
  
    double *rrs_s;
    rrs_s = (double *) calloc(s->nfitbands, sizeof (double)); 
    
    //Intialize output solution array
    for (iw=0;iw<s->nfree;iw++){
        s->popt[iw] = BAD_FLT;
    }

    //Compute subsurface remote-sensing reflectance
    for (iw = 0; iw < s->nfitbands; iw++) {
        rrs_s[iw] = s->Rrs_a[iw]/(0.52 + 1.7*s->Rrs_a[iw]);    
    }

    //model slope of bbp
    s->bbp_s =  2.0 * (1.0 - 1.2 * expf(-0.9 * (rrs_s[s->i440]/rrs_s[s->i555])));

    //model slope of adg
    s->adg_s = -(0.01447 + 0.00033*s->Rrs_a[s->i490]/s->Rrs_a[s->i555]);

    //Compute the spectral shapes adstar and bbpstar
    calc_adgstar(s, 0);
    calc_bbpstar(s, 0);
   
    /*Optimisation control parameters*/
    opts[0] = 1E-3; //1E-3 LM_INIT_MU; /* scale factor for initial mu */
    opts[1] = 1E-10; //1E-10;  convergence in gradient
    opts[2] = 1E-5; //1E-5 relative error desired in the approximate solution
    opts[3] = 1E-7; //1E-7 relative error desired in the sum of squares

    /*Run L-M optimization with analytical jacobian*/    
    statusLM = dlevmar_der(sdp_qssa, sdp_jacqssa, p,  rrs_s, s->nfree, 
    s->nfitbands, maxit, opts, info, NULL, NULL, (void *) s);

    free(rrs_s);
    
    //return solution if one exists
    if (statusLM <= 0) { 
        return 1;
    } else {
        for (iw=0;iw<s->nfree;iw++){
            s->popt[iw] = p[iw];
        }
        return 0;
    }
}

/* ------------------------------------------------------------------- */
/* check_l2_flags() - checks to see if l2 flags are set that indicate  */
/*                    poor quality Rrs                                 */    

/* Flags checked are: ATMFAIL(1), LAND (2), HIGLINT(4), PRODFAIL (31) */
/*                    STRAYLIGHT (9), CLOUD (10)                      */


int check_l2_flags(l2str *l2rec, int ip) {
    
    int i;
    int nflag = 6;
    int bit_position[6] = {1,2,4,9,10,31};
    
    int status = 0;
    
    //Check bit position in l2flags
    for (i=0;i<nflag;i++) {
        if (l2rec->l1rec->flags[ip] & (1<<bit_position[i])) {
            status = 1;
            break;
        } else {
            status = 0;
        }
    } 
    
    return status;
}

/* ------------------------------------------------------------------- */
/* set_sdp_flags() - checks range of derived pigments and sets quality */
/*                    flags                                            */

int set_sdp_flags(float pigvalue, int pigid, int ip, int16  *flag) {
    
    if (pigvalue < pigmin || pigvalue > pigmax) {
        
        if (pigid == 0) {
            *flag |= BADTCHL;
        }
        if (pigid == 1) {
            *flag |= BADZEA;
        }
        if (pigid == 2) {
            *flag |= BADDVCHLA;
        }
        if (pigid == 3) {
            *flag |= BADBUTOFUCO;
        }
        if (pigid == 4) {
            *flag |= BADHEXOFUCO;
        }
        if (pigid == 5) {
            *flag |= BADALLO;
        }
        if (pigid == 6) {
            *flag |= BADMVCHLB;
        }
        if (pigid == 7) {
            *flag |= BADNEO;
        }
        if (pigid == 8) {
            *flag |= BADVIOLA;
        }
        if (pigid == 9) {
            *flag |= BADFUCO;
        }
        if (pigid == 10) {
            *flag |= BADCHLC12 ;
        }
        if (pigid == 11) {
            *flag |= BASCHLC3;
        }
        if (pigid == 12) {
            *flag |= BADPERID;
        }
        
        pigvalue = BAD_FLT;
        
        return 1;
        
    } else { 
    
        return 0;
    }
}


/*----------------------------------------------------------------------*/
/* extract_band_3d() - utility function to extract selected band subset */

/*----------------------------------------------------------------------*/
static void extract_band_3d(float *in_buf, float *out_buf, int numPixels, int numBands) {
    float * out_ptr = out_buf;
    for (int pix = 0; pix < numPixels; pix++) {
        float *in_ptr = in_buf + pix * numBands;
        for (int band_3d = 0; band_3d < input->nwavelengths_3d; band_3d++) {
            int band = input->wavelength_3d_index[band_3d];
            *out_ptr = in_ptr[band];
            out_ptr++; 
        }
    }
}


/* ------------------------------------------------------------------- */
/* run_sdp() - computes the pigment concentrations from 2nd derivativ  */
/*             of the model residuals                                  */

/* ------------------------------------------------------------------- */

int run_sdp(l2str *l2rec) {

    static int firstCall = 1;
    
    static sdpstr sdpdata;
    static sdpstr* s = &sdpdata;

    int iw, iwx, ig, ip, irep, ipb, iwxf;
    int i400, i700, badrrs;
        
    float finewave = 1.0;
    float band_temp;
    double pigtemp_med;
    double pigtemp_std;
    double pigtemp;

    l1str *l1rec = l2rec->l1rec;
    float *bands = l1rec->l1file->fwave;
    int32_t npix = l1rec->npix;
    int32_t nbands = l1rec->l1file->nbands;
    
    static float *vbands;
    static double *fmrrs_s;
    static float *fmRrs_a;
    static float *RrsDiff;
    static float *deriv2;
    static float *fRrs_a;
    static float *fsmthRrs_a;
    static float *Rrs_a;
    static double *pigtempArr;

    if (firstCall) {

        firstCall = 0;

        //Check sensor compatability. Must be hyperspectral sensor.
        switch (l1rec->l1file->sensorID) {
        case OCI:
            break;
        case HICO:
            break;
        case PRISM:
            break;
        default:
            printf(
                    "-E- %s, %d: Spectral Derivative Pigment model only applicable to hypersectral data\n",
                    __FILE__, __LINE__);
            exit(1);
            break;
        }
      
        s->nfree = 3;
        i400 = windex(400, bands, nbands);
        i700 = windex(700, bands, nbands);
        s->nvbands = i700 - i400 + 1; //Number of sensor visible bands
        s->nfitbands = nfitbands; //Number of fit bands
        s->nfbands = (705 - 395 + 1)/finewave; //Number of fine resolution bands

        //Allocate output arrays
        if ((pigArr =(float *)  calloc(npix * npig, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for SDP.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((pigArr_unc =(float *)  calloc(npix * npig, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for SDP.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        //Allocate output arrays
        if ((sdpflags =(int16 *)  calloc(npix, sizeof (int16))) == NULL) {
            printf("-E- %s line %d : error allocating memory for SDP.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((vbands = (float *) calloc(s->nvbands, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for SDP.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((s->g = (float *) calloc(2, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory in init_iop_flag.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((s->bindx = (int *) calloc(s->nvbands, sizeof (int))) == NULL) {
            printf("-E- %s line %d : error allocating memory in init_iop_flag.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((s->fbindx = (int *) calloc(nfitbands, sizeof (int))) == NULL) {
            printf("-E- %s line %d : error allocating memory in init_iop_flag.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((s->bands = (float *) calloc(s->nvbands, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for SDP.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((s->fbands = (float *) calloc(s->nfbands, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for SDP.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((s->aw = (float *) calloc(s->nvbands, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for SDP.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((s->bbw = (float *) calloc(s->nvbands, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for SDP.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((s->aphA = (float *) calloc(s->nvbands, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for SDP.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((s->aphB = (float *)  calloc(s->nvbands, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for SDP.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((s->faw = (float *) calloc(s->nfbands, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for SDP.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((s->fbbw = (float *) calloc(s->nfbands, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for SDP.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((s->faphA = (float *) calloc(s->nfbands, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for SDP.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((s->faphB = (float *)  calloc(s->nfbands, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for SDP.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((s->adgstar = (float *)  calloc(s->nvbands, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for SDP.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((s->bbpstar = (float *)  calloc(s->nvbands, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for SDP.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((s->fadgstar = (float *)  calloc(s->nfbands, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for SDP.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((s->fbbpstar = (float *)  calloc(s->nfbands, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for SDP.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((s->Rrs_a = (float *)  calloc(s->nvbands, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for SDP.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((s->popt = (double*)  calloc(s->nfree, sizeof (double))) == NULL) {
            printf("-E- %s line %d : error allocating memory for SDP.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((s->popt = (double*)  calloc(s->nfree, sizeof (double))) == NULL) {
            printf("-E- %s line %d : error allocating memory for SDP.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((pigtempArr = (double*)  calloc(nreps, sizeof (double))) == NULL) {
            printf("-E- %s line %d : error allocating memory for SDP.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((mRrs_out = calloc(npix * nbands, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for SDP.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((Rrs_diff_out = calloc(npix * nbands, sizeof (float))) == NULL) {
            printf("-E- %s line %d : errorn allocating memory for SDP.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((d2Rrs_out = calloc(npix * nbands, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for SDP.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        
        
        //2D array of PCA coefficients (13 x 29000). Replicates are interleved in each row.
        s->Acoeff = (double **) allocate2d_double(npig, nreps*nAcoeff);
        
        //2D rray of PCA intercept coefficients (13 x 100)
        s->Ccoeff = (double **) allocate2d_double(npig, nreps);


        //Get the indices for the model fit bands
        for (iw = 0; iw < s->nfitbands; iw++) {
            iwx = windex(fitbands[iw], bands, nbands); 
            s->bindx[iw] = iwx;
        }


         //Populate model spectral constants at model fit bands
        for (iw=0; iw<s->nfitbands; iw++) {            
            
            //Subset of the model fit bands
            s->bands[iw] = bands[s->bindx[iw]];
            
            //Get phytoplankton absorption model coefficients
            aph_kramer_2022(s->bands[iw], &s->aphA[iw], &s->aphB[iw]);

            //get pure water absorption coefficients at fit bands
            get_sdp_aw(s->bands[iw], &s->aw[iw]) ;
            
        }
         
        //Compute reflectance model coefficients at fine resolution (400 - 700 nm, 1nm res)
        for (iw=0;iw<s->nfbands;iw++) {
            
            if (!iw){
                s->fbands[0] = 395;
            } else {
                s->fbands[iw] += s->fbands[iw -1] + 1;
            }
            
            //Get phytoplankton absorption model coefficients at 1nm
            aph_kramer_2022(s->fbands[iw], &s->faphA[iw], &s->faphB[iw]);
          
            //Get the water absorption coefficient at 1 nm
            get_sdp_aw(s->fbands[iw], &s->faw[iw]) ;
        }

        //Get bands at which model fits are calculated
        for (iw = 0; iw < s->nfitbands; iw++) {
            iwxf = windex(fitbands[iw], s->fbands, s->nfbands);
            s->fbindx[iw] = iwxf;
        }

        //Find array index for required wavelengths
        s->i440 = windex(440, s->bands, s->nvbands);
        s->i490 = windex(490, s->bands, s->nvbands);
        s->i555 = windex(550, s->bands, s->nvbands);
        s->i440f = windex(440, s->fbands, s->nfbands);
        
        //Get the SPD PCA model coefficients
        get_sdp_coeff(s->Acoeff, s->Ccoeff, pignames, npig, nAcoeff);

        //set Gordon rrs model coefficients
        s->g[0] = 0.0949;
        s->g[1] = 0.0794;
        
    }
    
    
    if ((fmrrs_s = (double *)  calloc(s->nfbands, sizeof (double))) == NULL) {
        printf("-E- %s line %d : error allocating memory for SDP.\n",
                __FILE__, __LINE__);
        exit(1);
    }
    if ((Rrs_a = (float *)  calloc(nbands, sizeof (float))) == NULL) {
        printf("-E- %s line %d : error allocating memory for SDP.\n",
                __FILE__, __LINE__);
        exit(1);
    }
    if ((fRrs_a = (float *)  calloc(s->nfbands, sizeof (float))) == NULL) {
        printf("-E- %s line %d : error allocating memory for SDP.\n",
                __FILE__, __LINE__);
        exit(1);
    }
    if ((fsmthRrs_a = (float *)  calloc(s->nfbands, sizeof (float))) == NULL) {
        printf("-E- %s line %d : error allocating memory for SDP.\n",
                __FILE__, __LINE__);
        exit(1);
    }
    if ((fmRrs_a = (float *)  calloc(s->nfbands, sizeof (float))) == NULL) {
        printf("-E- %s line %d : error allocating memory for SDP.\n",
                __FILE__, __LINE__);
        exit(1);
    }
    if ((RrsDiff= (float *)  calloc(s->nfbands, sizeof (float))) == NULL) {
        printf("-E- %s line %d : error allocating memory for SDP.\n",
                __FILE__, __LINE__);
        exit(1);
    }
    if ((deriv2= (float *)  calloc(s->nfbands, sizeof (float))) == NULL) {
        printf("-E- %s line %d : error allocating memory for SDP.\n",
                __FILE__, __LINE__);
        exit(1);
    }
    
    
    
    //Initialize outputs for scanline
    for (ip=0; ip < npix; ip++) {
        
        sdpflags[ip] = 0;
        
        for (ig=0; ig < npig; ig++) {
          pigArr[ip*npig + ig] = BAD_FLT;
          pigArr_unc[ip*npig + ig] = BAD_FLT;
        }

        for (iw=0;iw<nbands;iw++){

            mRrs_out[ip*nbands + iw] = BAD_FLT; 
            Rrs_diff_out[ip*nbands + iw] = BAD_FLT;
            d2Rrs_out[ip*nbands + iw] = BAD_FLT; 

        }
   
    }
   
    //Process pixels across a scan line
    for (ip = 0; ip<l1rec->npix; ip++) {
           
        ipb = ip*l1rec->l1file->nbands;
        
        for (iw=0;iw<nbands;iw++){
            Rrs_a[iw] = l2rec->Rrs[ipb + iw];
        }

        //Interpolate the observed Rrs to 1 nm resolution
        sdp_interp(bands, Rrs_a, s->fbands, fRrs_a, nbands, s->nfbands);
        
        //Smooth the interpolated input Rrs data using 5nm median filter
        //sdp_filter(fRrs_a, fsmthRrs_a, s->nfbands, 5);

        sdp_mean_filter(fRrs_a, fsmthRrs_a, s->nfbands, 5);
        
        //Iniitalize bad rrs flag (helps us flag invalid Rrs_a values)
        badrrs=0;

        //Initialize the seawater absorption coefficient at 1nm resolution
        for (iw = 0; iw < s->nfitbands; iw++) {

            //Load seawater absorption and backscattering coefficients
            //at model fit bands
            //s->aw[iw] = l1rec->sw_a[ipb + s->bindx[iw]];

            get_sdp_aw(bands[s->bindx[iw]], &s->aw[iw]) ;

            s->bbw[iw] = l1rec->sw_bb[ipb + s->bindx[iw]]; 

            //Get the smoothed Rrs values at model fit bands
            s->Rrs_a[iw] = fsmthRrs_a[s->fbindx[iw]];

            //Check for bad Rrs
            if (s->Rrs_a[iw] == BAD_FLT){
                badrrs += 1;
            }
        }
        
        //Calculate bbw at fine 1nm resolution
        for (iw = 0; iw < s->nfbands; iw++) {
            s->fbbw[iw] = seawater_bb(s->fbands[iw], l1rec->sstref[ip], 
                    l1rec->sssref[ip], 0.039);
        }
           
            
        //Run the inverse semianalytical model
        if(run_sdp_saa(s)) {
            l2rec->l1rec->flags[ip] |= PRODFAIL;
        };


        //Compute the spectral shapes adsgtar and bbpstar at 1nm resolution
        calc_adgstar(s, 1);
        calc_bbpstar(s, 1);

        //Run forward model to generatd Rrs at fine resolution (e.g., 1nm)
        sdp_qssaf(s->popt, fmrrs_s, s->nfree, s->nfbands, s);
        

        //Compute the difference between observed and modelled above water Rrs
        //at 1 nm resolution
        for (iw=0;iw<s->nfbands;iw++) {
            fmRrs_a[iw] = (0.52 * fmrrs_s[iw]) / (1.0 - 1.7 * fmrrs_s[iw]);
            RrsDiff[iw] = fsmthRrs_a[iw] - fmRrs_a[iw]; 
        }

        //compute the 2nd derivative
        sdp_deriv2(s->fbands, RrsDiff, deriv2, s->nfbands);

        //Use all 100 x 13 A pigment coefficients. Use all 100 C coefficients.
        //Compute the pigment concentrations for each 100 coefficient replicates
        for (ig=0;ig<npig;ig++){ 

            pigtemp_med = 0.0;
            pigtemp_std = 0.0;

            for (irep=0;irep<nreps;irep++) {

                //reset temporary variables
                pigtemp=0.0;
                pigtempArr[irep] = 0.0;
                
                for (iw=0;iw<nAcoeff;iw++) {
                    //Compute PCA sum 401 - 699 nm
                    //Note: PCA model coefficients start at 401nm and derivative Rrs starts at 399nm.
                    //So we need to start deriv2 multiplication at +6nm.
                    pigtemp += deriv2[iw+6] * s->Acoeff[ig][irep*nAcoeff + iw];
                }
                
                //Add the intercept for the ith replicate pigment    
                pigtempArr[irep] = pigtemp + s->Ccoeff[ig][irep];       

            }
            
            //Compute median and standard deviation of the 100 replicate pigment computations
            pigtemp_med = gsl_stats_median(pigtempArr,1, nreps);
            pigtemp_std = sqrt(gsl_stats_variance(pigtempArr, 1, nreps));
          
            //check range of pigments and and flag s necessary
            set_sdp_flags(pigtemp_med, ig, ip, &sdpflags[ip]);  
            
            pigArr[ip*npig + ig] = pigtemp_med;
            pigArr_unc[ip*npig + ig] = pigtemp_std;
        }

        //Extract the modelled Rrs, Rrs_diff, and 2nd derivative for output
        //at bands closest to OCI
        for (iw = 0; iw < nbands; iw++) {
            
            band_temp = bands[iw];
            iwx = windex(band_temp, s->fbands, s->nfbands);

            if (band_temp > 400 && band_temp < 700) {
                mRrs_out[ip*nbands + iw] = fmRrs_a[iwx];
                Rrs_diff_out[ip*nbands + iw]  = RrsDiff[iwx];
                d2Rrs_out[ip*nbands + iw] = deriv2[iwx];

            } else {
                mRrs_out[ip*nbands + iw] = BAD_FLT;
                Rrs_diff_out[ip*nbands + iw] = BAD_FLT;
                d2Rrs_out[ip*nbands + iw] = BAD_FLT;
            }
        
        }

    }

    LastRecNum = l1rec->iscan;
    
    free(Rrs_a);
    free(fmrrs_s);
    free(fRrs_a);
    free(fmRrs_a);
    free(RrsDiff);
    free(deriv2);
    return 0;
}


/* ------------------------------------------------------------------- */
/* get_sdp() - returns requested spectral derivative pigment products  */

/* ------------------------------------------------------------------- */
void get_sdp(l2str *l2rec, l2prodstr *p, float prod[]) {
    
    l1str *l1rec = l2rec->l1rec;
    int ib = p->prod_ix;
    
    int32_t ip, ipg, ipb;
    int32_t npix = l2rec->l1rec->npix;

    //Check if calculations have been computed already
    if (!sdp_ran(l1rec->iscan)) {
        run_sdp(l2rec);
    }

    if(p->rank == 3) {
        
        switch (p->cat_ix) {

        case CAT_mRrs_sdp:
            extract_band_3d( mRrs_out, prod, l2rec->l1rec->npix, l2rec->l1rec->l1file->nbands);
            break;
        case CAT_mRrs_diff_sdp:
            extract_band_3d( Rrs_diff_out, prod, l2rec->l1rec->npix, l2rec->l1rec->l1file->nbands);
            break;
        case CAT_d2Rrs_diff_sdp:
            extract_band_3d( d2Rrs_out, prod, l2rec->l1rec->npix, l2rec->l1rec->l1file->nbands);
            break;

        }

    } else {
    
        for (ip = 0; ip < npix; ip++) {
            
            ipg = ip* npig;
            ipb = ip * l1rec->l1file->nbands + ib;
        
            switch (p->cat_ix) {
                case CAT_tchl_sdp:
                    prod[ip] = (float) pigArr[ipg ];
                    break;
                case CAT_zea_sdp:
                    prod[ip] = (float) pigArr[ipg + 1];
                    break;
                case CAT_dvchla_sdp:
                    prod[ip] = (float)  pigArr[ipg + 2];
                    break;
                case CAT_butfuco_sdp:
                    prod[ip] = (float) pigArr[ipg + 3];
                    break;
                case CAT_hexfuco_sdp:
                    prod[ip] = (float) pigArr[ipg + 4];
                    break;
                case CAT_allo_sdp:
                    prod[ip] = (float) pigArr[ipg+ 5];
                    break;
                case CAT_mvchlb_sdp:
                    prod[ip] = (float) pigArr[ipg + 6];
                    break;
                case CAT_neo_sdp:
                    prod[ip] = (float) pigArr[ipg + 7];
                    break;
                case CAT_viola_sdp:
                    prod[ip] = (float) pigArr[ipg+ 8];
                    break;
                case CAT_fuco_sdp:
                    prod[ip] = (float) pigArr[ipg+ 9];
                    break;
                case CAT_chlc12_sdp:
                    prod[ip] = (float) pigArr[ipg + 10];
                    break;
                case CAT_chlc3_sdp:
                    prod[ip] = (float) pigArr[ipg + 11];
                    break;
                case CAT_perid_sdp:
                    prod[ip] = (float) pigArr[ipg + 12];
                    break;
                case CAT_flags_sdp:
                    prod[ip] = sdpflags[ip];
                    break;
                case CAT_mRrs_sdp:
                    prod[ip] = (float) mRrs_out[ipb];
                    break;
                case CAT_mRrs_diff_sdp:
                    prod[ip] = (float) Rrs_diff_out[ipb];
                    break;       
                case CAT_d2Rrs_diff_sdp:
                    prod[ip] = (float) d2Rrs_out[ipb];
                   break;    
                default:
                    printf("Error: %s : Unknown product specifier: %d\n", __FILE__, p->cat_ix);
                    exit(FATAL_ERROR);
                    break;
                }
        }
    }

}