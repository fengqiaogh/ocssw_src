/*=========================================================================*/
/*                             get_gpig.c                                   */
/*                                                                         */
/* Description:                                                            */ 
/* This algorithm calculates multiple phytoplankton pigment concentrations */
/* based on a bio-optical model that parameterized the phytoplankton       */
/* absorption coefficient using Guassian bands.                            */
/*                                                                         */
/*                                                                         */
/* Notes:                                                                  */
/*                                                                         */
/* References:                                                             */
/* Chase, A.P., Boss, E., Cetinic, I. and W.H. Slade (2017) Estimation of  */
/* Phytoplankton Accessory Pigments From Hyperspectral Reflectance Spectra:*/
/* Toward a Global Algorithm, JGR Oceans,122(12)9725-9743,                 */
/* doi: 10.1002/2017JC012859                                               */
/*                                                                         */
/* Implementation:                                                         */
/* L. McKinna (GO2Q/NASA GSFC), Dec 2022                                   */
/*=========================================================================*/

#include <math.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_filter.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_statistics.h>
#include <levmar.h>
#include "l12_proto.h"
#include "get_gpig.h"
#include "trustopt.h"

static int32_t LastRecNum = -1;

static int npig = 4;       //Number of pigments in model
static int ngauss = 8;
//static int ngauss = 12;
static int nfree = 15;
//static int nfree = 19;
static int n_mc = 10000;

static int16  *gpigflags;   //Binary quality control flags  
static float *pigArr;  //Pigment array pointer (this is where we'll store the solutions)
static float *pigArr_unc; //Pigment uncertinty array pointer (this is where we'll store the 
                       //uncertainty in solutions)
static float *mRrs_out; //Array holding forward-modelled Rrs
static double *pig_mc;

// Define the center peak locations (nm) and widths (nm) of the Gaussian functions
// sig = sigma, where FWHM = sigma*2.355 and FWHM = full width at half max
static float gauss_center[8] = {384,413,435,461,464,490,532,583};
//static float gauss_center[12] = {384,413,435,461,464,490,532,583,623,644,655,676};
static float gauss_fwhm[8] = {23,9,14,11,19,19,20,20};
//static float gauss_fwhm[12] = {23,9,14,11,19,19,20,20,15,12,12,9};


//Model coefficients for deriving four pigments
//order of coefficients: Ctchla, chlc12, tchlb, ppc
//static char *pignames[4] = {"tchla","chlc12","tchlb","ppc"}; //Derived pigments
static float pigwave[4] = {435,461,464,490}; //Gaussian peaks asscociated with pigments
static float pigCoeffsA[4] = {0.048,0.043,0.033, 0.079};
static float pigCoeffsB[4] = {0.643,0.561,0.327,0.823};
static float pigCoeffsA_unc[4] = {0.008,0.009,0.013,0.024};
static float pigCoeffsB_unc[4] = {0.068,0.059,0.074,0.105}; 


void gpig_qssa( double *p, double *mrrs_s, int m, int n, void *data);

/* ------------------------------------------------------------------- */
/*gpig_ran() - record if gpig_ran has been executed for scan line       */

/* ------------------------------------------------------------------- */
int gpig_ran(int recnum) {
    if (recnum == LastRecNum)
        return 1;
    else
        return 0;
}

/* ------------------------------------------------------------------- */
/*gpig_randn() - generate random  number that follow normal distribution */
/* ------------------------------------------------------------------- */
double gpig_randn (double mu, double sigma) {
    
    double U1, U2, W, mult;
    static double X1, X2;
    static int call = 0;

    if (call == 1) {
        call = !call;
        return (mu + sigma * (double) X2);
    }
    
    W = 0;
    while ((W >= 1) || (W == 0)) {
      U1 = -1 + ((double) rand () / RAND_MAX) * 2;
      U2 = -1 + ((double) rand () / RAND_MAX) * 2;
      W = pow (U1, 2) + pow (U2, 2);
    }

    mult = sqrt ((-2 * log (W)) / W);
    X1 = U1 * mult;
    X2 = U2 * mult;

    call = !call;

    return (mu + sigma * (double) X1);
}

/* ------------------------------------------------------------------*/
/* get_gpig_aw() - read netcf file and return purewater absorption    */
/*                    coefficients                                   */

int get_gpig_aw(float wave, float *A)  {

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
        strcat(filename, "/common/gpig_pure_water.nc");
        printf("Reading GPIG pure water coefficients from %s.\n", filename);
        
        /* Does the file exist? */
        if (access(filename, F_OK) || access(filename, R_OK)) {
            printf("-E- %s: GPIG pure water coefficient file '%s' does not exist or cannot open.\n",
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

/* ------------------------------------------------------------------*/
/* get_gpig_aw_ts_cor() - read netcf file and return seawater T-S    */
/*                        correction coefficients                    */

int get_gpig_aw_ts_cor(float wave, float *psiT, float *psiS)  {

    static int firstCall = 1;

    int32_t status, netcdf_input;
    size_t nlam;
    int ncid, grpid, varid1, varid2, varid3, dimid ;

    static float *wtab;
    static float *psiTtab;
    static float *psiStab;
    static float *dStab;
    static float *dTtab;
    static int ntab = 0;

    if (firstCall) {

        char *filedir;
        char filename[FILENAME_MAX];      
        
        if ((filedir = getenv("OCDATAROOT")) == NULL) {
            printf("-E- %s: OCDATAROOT env variable undefined.\n", __FILE__);
            exit(1);
        }
        
        strcpy(filename, filedir);
        strcat(filename, "/common/gpig_aw_ts_cor.nc");
        printf("Reading GPIG pure water coefficients from %s.\n", filename);
        
        /* Does the file exist? */
        if (access(filename, F_OK) || access(filename, R_OK)) {
            printf("-E- %s: GPIG pure water coefficient file '%s' does not exist or cannot open.\n",
                    __FILE__, filename);
            exit(EXIT_FAILURE);
        }
        
        /* test for NetCDF input file */
        status = nc_open(filename, NC_NOWRITE, &ncid);
        netcdf_input = (status == NC_NOERR);

        if (netcdf_input) {
            
            /* Get the group ids of our two groups. */
            if ((nc_inq_ncid(ncid, "water_abs_ts_corr", &grpid)) == NC_NOERR) {

                //Get the dimensions of the variable
                nc_inq_dimid(grpid,"n_wavelength", &dimid);
                nc_inq_dimlen(grpid, dimid, &nlam);
                ntab = (int) nlam;

                //Allocate memory
                psiTtab = (float*) allocateMemory(ntab * sizeof (float), "psiTtab");
                psiStab = (float*) allocateMemory(ntab * sizeof (float), "psiStab");
                wtab = (float*) allocateMemory(ntab * sizeof (float), "wtab");
                dTtab = (float*) allocateMemory(ntab * sizeof (float), "dTtab");
                dStab = (float*) allocateMemory(ntab * sizeof (float), "dStab");
                
                nc_inq_varid(grpid, "psiT", &varid1);
                nc_get_var_float(grpid, varid1, psiTtab);

                nc_inq_varid(grpid, "psiT", &varid2);
                nc_get_var_float(grpid, varid2, psiTtab);

                nc_inq_varid(grpid, "wavelength", &varid3);
                nc_get_var_float(grpid, varid3, wtab);                
            
            }
        
            spline(wtab, psiTtab, nlam, 1e30, 1e30, dTtab);
            spline(wtab, psiStab, nlam, 1e30, 1e30, dStab);

            firstCall = 0;

        }
    }

    if (wave > wtab[ntab - 1]){
        wave = wtab[ntab - 1];
    }

    // interpolate coefficients to wavelength to get psiT and psiS
    // coefficients
    splint(wtab, psiTtab, dStab, ntab, wave, psiT);
    splint(wtab, psiStab, dStab, ntab, wave, psiS);

    return 0;
    
}


/* ------------------------------------------------------------------- */
/* gpig_interp() - linear 1D interpolation function                     */

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

int gpig_interp(float *x, float *y, float *xinterp, float *yinterp,
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
/* gpig_qssa() - Quasi-single scattering approximation (QSSA) forward   */
/*              reflectance model                                      */

/* ------------------------------------------------------------------- */

void gpig_qssa( double *p, double *umod, int m, int n, void *data){

    
    gpigstr *gp = ((gpigstr *)data);
    
    int iw, ig;
    double atot, bbtot;
    double aph, ap, ag, anap;
    double cp, bbp, bp;
    double numer, denom;
    
    for (iw = 0; iw < n; iw++) {

        aph = 0;
        //Calculate the total IOPs
        for (ig=0; ig<ngauss; ig++){
            aph += p[ig]*gp->aphGauss[ig][iw];
        }  
        
        anap = p[8]*exp(-p[9]*(gp->fitbands[iw] - gp->fitbands[gp->i400f]));
        ag = p[10]*exp(-p[11]*(gp->fitbands[iw] - gp->fitbands[gp->i400f]));

        cp = p[12]*pow(gp->fitbands[iw]/gp->fitbands[gp->i400f],-p[13]);
        ap = aph + anap;
        bp = cp - ap;
        bbp = p[14]*bp;


        atot = gp->aw[iw] + aph + ag + anap;
        bbtot = gp->bbw[iw] + bbp;
        numer = bbtot ;
        denom = bbtot + atot;

        umod[iw] = numer/denom;

    } 
} 

/* ------------------------------------------------------------------- */
/* gpig_qssa2() - Quasi-single scattering approximation (QSSA) forward   */
/*              reflectance model                                      */

/* ------------------------------------------------------------------- */

int gpig_qssa2( int m, int n, double *p, double *umod, double **dydpar, void *data){

    
    gpigstr *gp = ((gpigstr *)data);
    
    int iw, ig;
    double atot, bbtot;
    double aph, ap, ag, anap;
    double cp, bbp, bp;
    double numer, denom;
    double dcpd12;
    double dcpd13;
    double ddenomdp, dnumerdp;
    double temp = 0;

    double *dadp;
    double *dbbdp;
    dadp = (double *) calloc(n, sizeof (double)); 
    dbbdp = (double *) calloc(n, sizeof (double)); 

    for (iw = 0; iw < m; iw++) {

        aph=0;
        anap=0;
        ag=0;
        ap=0;
        cp=0;
        
        //Calculate the total IOPs
        for (ig=0; ig < ngauss; ig++){
            dadp[ig] = gp->aphGauss[ig][iw];   
            aph += p[ig]*gp->aphGauss[ig][iw];
        }
        
        anap = p[8]*exp(-p[9]*(gp->fitbands[iw] - gp->fitbands[gp->i400f]));
        ag = p[10]*exp(-p[11]*(gp->fitbands[iw] - gp->fitbands[gp->i400f]));

        dadp[8] =  exp(-p[9]*(gp->fitbands[iw] - gp->fitbands[gp->i400f]));
        dadp[9] = -(gp->fitbands[iw] - gp->fitbands[gp->i400f])*anap;

        dadp[10] = exp(-p[11]*(gp->fitbands[iw] - gp->fitbands[gp->i400f]));
        dadp[11] = -(gp->fitbands[iw] - gp->fitbands[gp->i400f])*ag;

        dadp[12] = 0.0;
        dadp[13] = 0.0;
        dadp[14] = 0.0;

        cp = p[12]*pow(gp->fitbands[iw]/gp->fitbands[gp->i400f],-p[13]);

        dcpd12 =   pow(gp->fitbands[iw]/gp->fitbands[gp->i400f],-p[13]);
        dcpd13 = -cp*log(gp->fitbands[iw]/gp->fitbands[gp->i400f]);

        ap = aph + anap;
        bp = cp - ap;
        bbp = p[14]*bp;

        for (ig=0;ig<n;ig++) {
            dbbdp[ig] = -p[14]*dadp[ig] ;   
        }

        dbbdp[8] =  -p[14]*dadp[8];
        dbbdp[9] =  -p[14]*dadp[9];
        dbbdp[10] = 0.0;
        dbbdp[11] = 0.0;
        dbbdp[12] = p[14]*dcpd12;
        dbbdp[13] = p[14]*dcpd13;
        dbbdp[14] = bp;


        atot = gp->aw[iw] + aph + ag + anap;
        bbtot = gp->bbw[iw] + bbp;
        numer = bbtot ;
        denom = bbtot + atot;

        umod[iw] = numer/denom;

        for (ig=0;ig<n;ig++){

            temp = 0.0;
            dnumerdp = dbbdp[ig];
            ddenomdp = dbbdp[ig] + dadp[ig];
            temp = (dnumerdp*denom - ddenomdp*numer )/ pow(denom,2.);
            dydpar[iw][ig] = temp;

        }    
    } 

    free(dadp);
    free(dbbdp);
    return 0;
} 

/* ------------------------------------------------------------------- */
/* gpig_jacqssa() - Compute Jacobian of cost function                   */
/*                 e    = y - y(x)                                     */
/*                 dedx =   - dy(x)/dx                                 */

/* dimension of jac is: npar x nbands                                  */

/* ------------------------------------------------------------------- */


void gpig_jacqssa(double *p, double *jac, int n, int m, void* data) {
    
    gpigstr *gp = ((gpigstr *)data);

    int i, j, ig;
    double aph, ap, ag, anap;
    double cp, bbp, bp;
    double atot, bbtot;
    double numer, denom;
    double dcpd12;
    double dcpd13;
    double temp;

    double *dadp;
    double *dbbdp;

    dadp = (double *) calloc(n, sizeof (double)); 
    dbbdp = (double *) calloc(n, sizeof (double)); 
    
    double dnumerdp,ddenomdp;
      
     //                 p0    p1   02   p3   p4   p5   p6   p7   p8    p9     p10   p11     p12   p13    p14
     //free parameters: mg1, mg2, mg3, mg4, mg5, mg6, mg7, mg8, mnap,  s_nap, mag,  s_ag,   mcp, s_cp,  bbpbp   
    j = 0;
    for (i=0; i< m; i++) {

        aph=0;
        anap=0;
        ag=0;
        ap=0;
        cp=0;

        for (ig=0;ig<ngauss;ig++) {
            dadp[ig] = gp->aphGauss[ig][i];   
            aph += p[ig]*gp->aphGauss[ig][i];    
        }

        anap = p[8]*exp(-p[9]*(gp->  fitbands[i] - gp->fitbands[gp->i400f]));
        dadp[8] =  exp(-p[9]*(gp->fitbands[i] - gp->fitbands[gp->i400f]));
        dadp[9] = -(gp->fitbands[i] - gp->fitbands[gp->i400f])*anap;
        
        ag = p[10]*exp(-p[11]*(gp->fitbands[i] - gp->fitbands[gp->i400f]));
        dadp[10] = exp(-p[11]*(gp->fitbands[i] - gp->fitbands[gp->i400f]));
        dadp[11] = -(gp->fitbands[i] - gp->fitbands[gp->i400f])*ag;

        dadp[12] = 0.0;
        dadp[13] = 0.0;
        dadp[14] = 0.0;

        ap = aph + anap;
        
        cp = p[12]*pow(gp->fitbands[i]/gp->fitbands[gp->i400f],-p[13]);
        dcpd12 =   pow(gp->fitbands[i]/gp->fitbands[gp->i400f],-p[13]);
        dcpd13 = -cp*log(gp->fitbands[i]/gp->fitbands[gp->i400f]);
        
        bp = cp - ap;
        bbp = p[14]*bp;

        for (ig=0;ig<ngauss;ig++) {
            dbbdp[ig] = -p[14]*dadp[ig] ;   
        }

        dbbdp[8] =  -p[14]*dadp[8];
        dbbdp[9] =  -p[14]*dadp[9];
        dbbdp[10] = 0.0;
        dbbdp[11] = 0.0;
        dbbdp[12] = p[14]*dcpd12;
        dbbdp[13] = p[14]*dcpd13;
        dbbdp[14] = bp;
        

        atot = gp->aw[i] + aph +anap + ag;
        bbtot = gp->bbw[i] + bbp;
        
        numer = bbtot;
        denom = bbtot + atot;
        //u = numer/denom
        //dudp = (dnumerdp*denom - ddenomdp*numer) / (denom)^2.

        for (ig=0;ig<n;ig++){

            temp = 0.0;

            dnumerdp = dbbdp[ig];
            ddenomdp = dbbdp[ig] + dadp[ig];

            temp = (dnumerdp*denom - ddenomdp*numer )/ pow(denom,2.);
            jac[j++] = temp;

        }
 
    }

    free(dadp);
    free(dbbdp);

}

/* ------------------------------------------------------------------- */
/* run gpig_saa() - runs gpig semianalytical algorithm to get to mRrs    */

/* ------------------------------------------------------------------- */
int run_gpig_saa(gpigstr *gp) {
  
    int iw;
    int statusLM;
    int maxit = 1500;

    double rrs_s, urrs_s;
    double temp;
    //coefficients: mg1,   mg2,   mg3,   mg4,   mg5,   mg6,   mg7,   mg8,   mnap,  s_nap, mag,  s_ag,   mcp,  s_cp,  bbpbp   
    double pl[15] ={0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.005, 0.0, 0.005,  0.0, 0.0,   0.005};
    double p[15] = {0.01,  0.01,  0.01,  0.01,  0.01,  0.01,  0.01,  0.01,  0.005, 0.011, 0.1,  0.0185, 0.1,  1.0,   0.01};
    double pu[15] ={0.5,   0.5,   0.5,   0.5,   0.5,   0.5,   0.5,   0.5,   0.05,  0.016, 0.8,  0.02,   1.0,  1.3,   0.1};

    double info[LM_INFO_SZ]; 
    double opts[LM_OPTS_SZ];
  
    double *uterm;
    double *dscl;
    uterm = (double *) calloc(gp->nfitbands, sizeof (double)); 
    dscl = (double *) calloc(gp->nfitbands, sizeof (double)); 
    
    //Intialize output solution array
    for (iw=0;iw<gp->nfree;iw++){
        gp->popt[iw] = BAD_FLT;
    } 

    //Compute subsurface remote-sensing reflectance
    for (iw = 0; iw < gp->nfitbands; iw++) {
        rrs_s = gp->Rrs_a[iw]/(0.52 + 1.7*gp->Rrs_a[iw]);
        uterm[iw] =  (-1.0*gp->g[0] + sqrt(pow(gp->g[0],2.) + 4.*gp->g[1]*rrs_s  ))/(2*gp->g[1]); 
        //printf("%f,%f,%f\n",gp->Rrs_a[iw],rrs_s,uterm[iw]);
        
        //Compute input uncertainties in uterm
        //deriv1 = 0.52/pow(0.52 + 1.7*gp->Rrs_a[iw],2);
        //urrs_s = sqrt(pow(deriv1*gp->Rrs_unc[iw],2));
        //deriv2 = 1./sqrt(pow(gp->g[0],2.) + 4.*gp->g[1]*rrs_s );

        //Uncertainties go in dscl array for scaling cost function.
        //dscl[iw] = sqrt(pow(deriv2*urrs_s,2)); 

        //Ali's code
        urrs_s = gp->Rrs_unc[iw]/(0.52 + 1.7 * gp->Rrs_unc[iw]);
        temp = (-1.0*gp->g[0] + sqrt(pow(gp->g[0],2.) + 4.*gp->g[1]*urrs_s  ))/(2*gp->g[1]); 
        dscl[iw] = pow(temp,2.);


    }
    //printf("\n");

   
    /*Optimisation control parameters*/
    opts[0] = 1E-3; //1E-3 LM_INIT_MU; /* scale factor for initial mu */
    opts[1] = 1E-12; //1E-10;  convergence in gradient
    opts[2] = 1E-12; //1E-5 relative error desired in the approximate solution
    opts[3] = 1E-12; //1E-7 relative error desired in the sum of squares
    opts[4] = 1E-6; //1E-6  LM_DIFF_DELTA;  /*step used in difference approximation to the Jacobian*/


    /*Run L-M optimization with analytical jacobian*/    
    //statusLM = dlevmar_der(gpig_qssa, gpig_jacqssa, p,  uterm, gp->nfree, 
    //gp->nfitbands, maxit, opts, info, NULL, NULL, (void *) gp);
    
    //statusLM = dlevmar_dif(gpig_qssa, p, uterm, gp->nfree, gp->nfitbands,
    //        maxit, opts, info, NULL, NULL, (void *) gp);
    

    //Constrained box constraints LM
    //statusLM = dlevmar_bc_dif(gpig_qssa, p, uterm, gp->nfree, gp->nfitbands,
    //       pl, pu, dscl, maxit, opts, info, NULL,
    //        NULL, (void *) gp);

    statusLM = dlevmar_bc_der(gpig_qssa,gpig_jacqssa, p, uterm, gp->nfree, gp->nfitbands,
           pl, pu, dscl, maxit, opts, info, NULL,
            NULL, (void *) gp);


    free(uterm);

    
    //return solution if one exists
    if (statusLM <= 0) { 
        return 1;
    } else {
        for (iw=0;iw<gp->nfree;iw++){
            gp->popt[iw] = p[iw];
        }
        return 0;
    }
}


/* ------------------------------------------------------------------- */
/* run gpig_saa2() - runs gpig semianalytical algorithm to get to mRrs    */

/* ------------------------------------------------------------------- */
int run_gpig_saa2(gpigstr *gp) {
  
    int i, iw;
    int status;

    trustopt_result result;
    trustopt_par pars[15];
   
    double rrs_s, urrs_s, temp;
    
   //coefficients: mg1,   mg2,   mg3,   mg4,   mg5,   mg6,   mg7,   mg8,    mnap,  s_nap, mag,  s_ag,   mcp,  s_cp,  bbpbp   
    double pl[15] ={0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.005, 0.01, 0.005,   0.01, 0.0,   0.005};
    double p[15] = {0.01,  0.01,  0.01,  0.01,  0.01,  0.01,  0.01,  0.01,  0.005, 0.011, 0.1,  0.0185, 0.1,  1.0,   0.01};
    double pu[15] ={0.5,   0.5,   0.5,   0.5,   0.5,   0.5,   0.5,   0.5,   0.05,  0.016, 0.8,  0.02,   1.0,  1.3,   0.015};
    
    
    //Set up the parameter struct
    for (i=0;i<15;i++) {
        //set all parameters to free flag
        pars[i].fixed = 0;
        //set lower and upper limits flag
        pars[i].limited[0] = 1;  
        pars[i].limited[1] = 1;  

        //set lower limit values
        pars[i].limits[0] = pl[i];
        //set lower limit vlalues
        pars[i].limits[1] = pu[i];

        //Use analytical jacobian
        pars[i].side=3;
        pars[i].step = 1E-6;
        pars[i].relstep = 1E-6;

    } 

    double *uterm;
    double **covy;
    uterm = (double *) calloc(gp->nfitbands, sizeof (double)); 
    covy = (double **) calloc(gp->nfitbands, sizeof (double *));

    for (i = 0; i < gp->nfitbands; i++) {
        covy[i] = (double *) calloc(gp->nfitbands , sizeof (double));
    }
    
    //Intialize output solution array
    for (iw=0;iw<gp->nfree;iw++){
        gp->popt[iw] = BAD_FLT;
    } 

    //Compute subsurface remote-sensing reflectance
    for (iw = 0; iw < gp->nfitbands; iw++) {
        rrs_s = gp->Rrs_a[iw]/(0.52 + 1.7*gp->Rrs_a[iw]);
        uterm[iw] =  (-1.0*gp->g[0] + sqrt(pow(gp->g[0],2.) + 4.*gp->g[1]*rrs_s  ))/(2*gp->g[1]);  
        
        //Compute input uncertainties in uterm
        //deriv1 = 0.52/pow(0.52 + 1.7*gp->Rrs_a[iw],2);
        //urrs_s = sqrt(pow(deriv1*gp->Rrs_unc[iw],2));
        //deriv2 = 1./sqrt(pow(gp->g[0],2.) + 4.*gp->g[1]*rrs_s );

        //Set diaganoal elements of covariablce matrix,
        //off-diagonal intialised to zero.
        //covy[iw][iw] = pow(deriv2*urrs_s,2.);

        //Ali's code
        urrs_s = gp->Rrs_unc[iw]/(0.52 + 1.7 * gp->Rrs_unc[iw]);
        temp = (-1.0*gp->g[0] + sqrt(pow(gp->g[0],2.) + 4.*gp->g[1]*urrs_s  ))/(2*gp->g[1]); 
        covy[iw][iw] = pow(temp,2.);

    }

    status = trustopt(gpig_qssa2, gp->nfitbands, gp->nfree, 1, 
                    p, uterm, covy, pars, NULL,  
                    (void *) gp, &result);


    gpig_qssa2( gp->nfitbands,gp->nfree, p, uterm, covy, (void *) gp);
    
    free(uterm);
    
    //return solution if one exists
    if (status < 0) { 
        return 1;
    } else {
        for (iw=0;iw<gp->nfree;iw++){
            gp->popt[iw] = p[iw];
        }
        return 0;
    }
    
}


/* ------------------------------------------------------------------- */
/* check_l2_flags_gpig() - checks to see if l2 flags are set that indicate  */
/*                    poor quality Rrs                                 */    

/* Flags checked are: ATMFAIL(1), LAND (2), HIGLINT(4), PRODFAIL (31) */
/*                    STRAYLIGHT (9), CLOUD (10)                      */


int check_l2_flags_gpig(l2str *l2rec, int ip) {
    
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
/* set_gpig_flags() - checks range of derived pigments and sets quality */
/*                    flags                                            */

int set_gpig_flags(float pigvalue, int pigid, int ip, int16  *flag) {
    
    if (pigvalue < pigmin || pigvalue > pigmax) {
        
        if (pigid == 0) {
            *flag |= BADTCHL;
        }
        if (pigid == 1) {
            *flag |= BADCHLC12;
        }
        if (pigid == 2) {
            *flag |= BADCHLB;
        }
        if (pigid == 3) {
            *flag |= BADPPC;
        }
        
        pigvalue = BAD_FLT;
        
        return 1;
        
    } else { 
    
        return 0;
    }
}


/* ------------------------------------------------------------------- */
/* run_gpig() - computes the pigment concentrations from Guassian bands*/

/* ------------------------------------------------------------------- */

int run_gpig(l2str *l2rec) {

    static int firstCall = 1;
    
    static gpigstr gpigdata;
    static gpigstr* gp = &gpigdata;

    int iw, iwx, ig, ip, ipb, iwxf;
    int i600, i700, badrrs;
        
    double pigtemp;
    double t_pope = 22.;

    l1str *l1rec = l2rec->l1rec;
    float *bands = l1rec->l1file->fwave;
    int32_t npix = l1rec->npix;
    int32_t nbands = l1rec->l1file->nbands;

    uncertainty_t *uncertainty=l1rec->uncertainty;
    
    float randn1;
    float randn2;
    float randn_A;
    float randn_B;

    

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
      
        gp->nfree = nfree;
        gp->i400 = windex(400, bands, nbands); //Find index of 400nm in sensor
        i600 = windex(600, bands, nbands); //Find index of 600nm in sensor
        i700 = windex(700,bands, nbands);//Find index of 700nm in sensor

        gp->nvbands = i700 - gp->i400 + 1; //Number of sensor visible bands
        gp->nfitbands =i600 - gp->i400 + 1;; //Number of fit bands

        //Allocate output arrays
        if ((pigArr =(float *)  calloc(npix * npig, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for GPIG.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((pig_mc =(double *)  calloc(n_mc, sizeof (double))) == NULL) {
            printf("-E- %s line %d : error allocating memory for GPIG.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((pigArr_unc =(float *)  calloc(npix * npig, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for GPIG.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((gpigflags =(int16 *)  calloc(npix, sizeof (int16))) == NULL) {
            printf("-E- %s line %d : error allocating memory for GPIG.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((gp->g = (float *) calloc(2, sizeof (float))) == NULL) {
            printf("-E- %s line %d :  error allocating memory for GPIG\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((gp->bindx = (int *) calloc(gp->nfitbands, sizeof (int))) == NULL) {
            printf("-E- %s line %d : error allocating memory for GPIG.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((gp->gindx = (int *) calloc(npig, sizeof (int))) == NULL) {
            printf("-E- %s line %d : error allocating memory for GPIG\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((gp->aw = (float *) calloc(gp->nfitbands, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for GPIG.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((gp->psiT = (float *) calloc(gp->nfitbands, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for GPIG.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((gp->psiS = (float *) calloc(gp->nfitbands, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for GPIG.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((gp->bbw = (float *) calloc(gp->nfitbands, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for GPIG\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((gp->Rrs_a = (float *)  calloc(gp->nfitbands, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for GPIG.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((gp->uterm = (float *)  calloc(gp->nfitbands, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for GPIG.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((gp->popt = (double*)  calloc(gp->nfree, sizeof (double))) == NULL) {
            printf("-E- %s line %d : error allocating memory for GPIG.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((gp->popt = (double*)  calloc(gp->nfree, sizeof (double))) == NULL) {
            printf("-E- %s line %d : error allocating memory for GPIG.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((gp->fitbands = (float*)  calloc(gp->nfitbands, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for GPIG.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((gp->Rrs_a = (float*) calloc(gp->nfitbands, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for GPIG.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((gp->Rrs_unc = (float*) calloc(gp->nfitbands, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for GPIG.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((mRrs_out = calloc(npix * nbands, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for GPIG.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
           
        //2D array of Gaussian bands
        gp->aphGauss = (float **) allocate2d_float(ngauss, gp->nfitbands);

        //Get the OCI band indices for the model(400 - 600 nm)
        for (iw = 0; iw < gp->nfitbands; iw++) {
            iwxf = gp->i400 + iw;
            gp->fitbands[iw] = bands[iwxf];
            iwx = windex(gp->fitbands[iw], bands, nbands); 
            gp->bindx[iw] = iwx;
        }

        //Find index of 400nm in truncated set of fit bands
        gp->i400f = windex(400, gp->fitbands, gp->nfitbands);

        //Get the indices of the Gaussains band peaks assciociated with pigment models
        for (ig=0;ig<npig;ig++) {
            iwx = windex( pigwave[ig], gauss_center, ngauss); 
            gp->gindx[ig] = iwx;
        } 

        //Generate Gaussian bands
        for (ig=0;ig<ngauss;ig++) {
            for (iw = 0; iw< gp->nfitbands;iw++){
                gp->aphGauss[ig][iw] =  exp(-0.5*pow((gp->fitbands[iw] - gauss_center[ig])/gauss_fwhm[ig],2 ));
            }
        }

        //set Gordon rrs model coefficients
        gp->g[0] = 0.0949;
        gp->g[1] = 0.0794;
        
    }
       
    //Initialize outputs for scanline
    for (ip=0; ip < npix; ip++) {

        gpigflags[ip] = 0;
        
        for (ig=0; ig < npig; ig++) {
          pigArr[ip*npig + ig] = BAD_FLT;
          pigArr_unc[ip*npig + ig] = BAD_FLT;
        }

        for (iw=0;iw<nbands;iw++){
            mRrs_out[ip*nbands + iw] = BAD_FLT;  
        }
   
    }
   
    //Process pixels across a scan line
    for (ip = 0; ip<l1rec->npix; ip++) {
           
        ipb = ip*l1rec->l1file->nbands;
        
        badrrs=0;

        //Get per-pixel spectral values
        for (iw = 0; iw < gp->nfitbands; iw++) {

            gp->Rrs_a[iw] = l2rec->Rrs[ipb + gp->bindx[iw]];

            if (uncertainty) {
                gp->Rrs_unc[iw] = l2rec->Rrs_unc[ipb + gp->bindx[iw]];
            } else {
               gp->Rrs_unc[iw] = 0.05*l2rec->Rrs[ipb + gp->bindx[iw]]; 
            }

            //Load seawater absorption and backscattering coefficients
            //at model fit bands
            get_gpig_aw(gp->fitbands[iw], &gp->aw[iw]) ;
            get_gpig_aw_ts_cor(gp->fitbands[iw], &gp->psiT[iw], &gp->psiS[iw]);

            gp->aw[iw] += gp->psiT[iw]*(l1rec->sstref[ip] - t_pope) + gp->psiS[iw]* l1rec->sssref[ip];

            //for testing use fixed SST and SSS for calculating bbw
            gp->bbw[iw] = seawater_bb(gp->fitbands[iw], l1rec->sstref[ip], l1rec->sssref[ip], 0.039);
            
            //Check for bad Rrs
            if (gp->Rrs_a[iw] == BAD_FLT){
                badrrs += 1;
            }
            
        }                 
            
        //Run the inverse semianalytical model
        if(run_gpig_saa(gp)) {
            l2rec->l1rec->flags[ip] |= PRODFAIL;
        };

        ///Run the inverse semianalytical model
        //if(run_gpig_saa2(gp)) {
        //    //l2rec->l1rec->flags[ip] |= PRODFAIL;
        //};


        for (ig=0;ig<npig;ig++){

            for (iw=0; iw<n_mc; iw++) {
                randn1 = (-1+2*((float)rand())/RAND_MAX)*pigCoeffsA_unc[ig];
                randn2 = (-1+2*((float)rand())/RAND_MAX)*pigCoeffsB_unc[ig];
                randn_A = randn1 + pigCoeffsA[ig];
                randn_B = randn2 + pigCoeffsB[ig];
                pig_mc[iw] = pow((gp->popt[gp->gindx[ig]] / randn_A ), (1./randn_B));  
            }

            pigtemp = (float) gsl_stats_median(pig_mc,0,n_mc);
            //printf("%d, %f\n",ig,pigtemp);
            pigArr[ip*npig + ig] = pigtemp;
        }

    }

    LastRecNum = l1rec->iscan;
    
    return 0;
}


/* ------------------------------------------------------------------- */
/* get_gpig() - returns requested spectral derivative pigment products  */

/* ------------------------------------------------------------------- */
void get_gpig(l2str *l2rec, l2prodstr *p, float prod[]) {
    
    l1str *l1rec = l2rec->l1rec;   
    int32_t ip, ipg;
    int32_t npix = l2rec->l1rec->npix;

    //Check if calculations have been computed already
    if (!gpig_ran(l1rec->iscan)) {
        run_gpig(l2rec);
    }
   
    for (ip = 0; ip < npix; ip++) {
        
        ipg = ip* npig;
        
        switch (p->cat_ix) {
            case CAT_tchl_gpig:
                prod[ip] = (float) pigArr[ipg ];
                break;
            case CAT_chlc12_gpig:
                prod[ip] = (float) pigArr[ipg + 1];
                break;
            case CAT_tchlb_gpig:
                prod[ip] = (float) pigArr[ipg + 2];
                break;
            case CAT_ppc_gpig:
                prod[ip] = (float) pigArr[ipg+ 3];
                break;
            case CAT_flags_gpig:
                prod[ip] = gpigflags[ip];
                break;   
            default:
                printf("Error: %s : Unknown product specifier: %d\n", __FILE__, p->cat_ix);
                exit(FATAL_ERROR);
                break;
            }
    }


}