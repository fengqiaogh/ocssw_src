/* =================================================================== */
/* trustopt.c  (trustopt)                                              */
/*                                                                     */ 
/*This program contains trust region optimization for least sqayares   */
/*with reflective boundary constraints. The approach is modified from  */
/* theoretical formulations devloped by Coleman and Li (1996).         */
/*                                                                     */ 
/* Dependencies: Gnu Scientific Library (GSL)                          */
/*                                                                     */ 
/* Version 1.0, July 2024                                              */
/*                                                                     */ 
/* Implementation: L. McKinna (Go2Q)                                   */ 
/*                                                                     */ 
/* References:                                                         */
/*                                                                     */ 
/* Coleman, T.F. and Y. Li (1996) An Interior Trust Region Approach    */
/* for Nonlinear Minimization Subject to Bounds, SIAM Journal on       */
/* Optimization, 6(2),418-445, doi: 10.1137/0806023.                   */
/*                                                                     */ 
/* Nocedal, J. and S. Wright (2006) Numerical Optimization (2nd Ed),   */
/* Springer-Verlag, New York City, USA, ISBN-10: 0-387-30303-0         */
/*                                                                     */
/*                                                                     */
/*Pauli Virtanen, Ralf Gommers, Travis E. Oliphant, Matt Haberland,  
Tyler Reddy, David Cournapeau, Evgeni Burovski, Pearu Peterson, 
Warren Weckesser, Jonathan Bright, Stéfan J. van der Walt, 
Matthew Brett, Joshua Wilson, K. Jarrod Millman, Nikolay Mayorov, 
Andrew R. J. Nelson, Eric Jones, Robert Kern, Eric Larson, CJ Carey, 
İlhan Polat, Yu Feng, Eric W. Moore, Jake VanderPlas, Denis Laxalde, 
Josef Perktold, Robert Cimrman, Ian Henriksen, E.A. Quintero, 
Charles R Harris, Anne M. Archibald, Antônio H. Ribeiro, 
Fabian Pedregosa, Paul van Mulbregt, and SciPy 1.0 Contributors. (2020) 
SciPy 1.0: Fundamental Algorithms for Scientific Computing in Python. 
Nature Methods, 17(3), 261-272. DOI: 10.1038/s41592-019-0686-2.
*/

/* =================================================================== */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include "trustopt.h"

#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sort_vector.h>

static int tropt_jac(trustopt_func funct,
	      int m, int nfree, int *ifree, int npar, double *x,
	      gsl_matrix *K, double epsfcn,
	      void *priv, int *nfev,
	      double *step, double *dstep, int *dside,
	      double **dydpar);

static int calc_tropt_gradf(gsl_matrix *Wy, gsl_vector *deltay, 
        gsl_vector *temp1, gsl_matrix *K, gsl_vector *gradf);

static int calc_hessian(gsl_matrix *Wy, gsl_matrix *W2, 
            gsl_matrix *K, gsl_matrix *H);

static int tropt_cost_actual(gsl_vector *deltay, gsl_vector *temp1, 
        gsl_matrix *Wy, int covyflag, double *fcost);

static int calc_cost_model(gsl_vector *xstep, gsl_vector *g,
         gsl_matrix *B, gsl_vector *temp2, double *mcost);

static  int set_affine_scaling_matrix(int *isllimit, int *isulimit, 
        double *llimit, double *ulimit, double *xall, gsl_vector *g,
        int nfree, int *ifree, gsl_vector *v, gsl_vector *dv ) ;

static int calc_mhat(gsl_vector *gradf, gsl_matrix *invD, gsl_matrix *H,
                int n, gsl_matrix *Mhat, gsl_matrix *Chat);

// static int calc_C(gsl_vector *g, gsl_matrix *B, gsl_matrix *D,
//         gsl_matrix *invD, int n, gsl_matrix *C);

// static int calc_ghat(gsl_vector *g, gsl_matrix *invD,
//         gsl_vector *ghat) ;

//static int scale_K(gsl_matrix *K, gsl_matrix *invD,
//        gsl_matrix *Khat);

static int tropt_matrix_invert(gsl_matrix *A, gsl_matrix *invA,
        int diagflag);

static int calc_quadratic_1d(gsl_matrix *J, gsl_vector *g, 
            gsl_vector *s, gsl_vector *s0,
            gsl_matrix *DIAG, int n, int m,
            double *a, double *b, double *c);

static int  eval_quadratic_1d(double a, double b, double c,
            double lb, double ub, double *tout, double *yout); 

// static int posdef_check(gsl_matrix *A, int n);

static int solve_subproblem(gsl_vector *g, gsl_matrix *H, 
        gsl_vector *p, double *delta, int n);

static int set_scaling_matrix(gsl_vector *g, gsl_matrix *D,
            gsl_matrix *invD, int n, int iter);

static int select_step(double *xall, double *xnew, double *llimit, double *ulimit,
            int *ifree, int *isllimit, int *isulimit, int nfree, int m,
            gsl_matrix *K_hat, gsl_matrix *Hhat, gsl_matrix *Chat, 
            gsl_vector *gradfhat, gsl_vector *xstep, gsl_vector *xstep_hat, 
            gsl_matrix *invD, gsl_matrix *D, 
            double *delta, double theta, double *mcost);

static int check_is_bounded(gsl_vector *s, double *x, int *ifree,
                  int nfree, double *ulimit, double *llimit,
                  double *xnew, gsl_vector *hit_boundary);

static double calc_alpha(double *xnew, double *ulimit, double *llimit, 
                        int *isulimit, int *isllimit, int *ifree, 
                        int nfree, gsl_vector *xstep, 
                        gsl_vector *hit_boundary);

double calc_lnsph_inter(gsl_vector *xhat, gsl_vector *shat, 
                double *delta);

static int update_trust_region(double mu, double eta, 
            double theta, double gamma1, double gamma2,
            double *actual_reduction, double *pred_reduction,
            double norm_xstep, double*rho, double *delta);
//Utility functions.
static int get_sign(double a) ;
static double tropt_min(double a, double b);
static double tropt_max(double a, double b); 
static double tropt_norm2_xall(double *x, int *ifree, int n);
// static double vector_abs_max(gsl_vector *v);
static int array_min(double *v, int n, double *min, int *idx);
static int tropt_free2D(double **Arr, int m);


#define tropt_call(funct, m, n, x, fvec, dydpar, priv) (*(funct))(m,n,x,fvec,dydpar,priv);

/***********************************************************************/
/*                       trustopt                                      */
/*                                                                     */
/*                                                                     */

int trustopt(trustopt_func funct, int m, int npar, int covydiag, 
     double *xall, double *yall, double **covy, 
     trustopt_par *pars, trustopt_config *config, 
            void *private_data, trustopt_result *result) {
    
    static trustopt_config conf; //configuration structure conf

    int i,j;
    int iter = 0;         // outer loop iteration counter
    int status=0;         //integer code for success/fail/error
    int nfree;            //integer number of free parameters in model
    int covyflag;         //Covariance of y flag       (0, none provided)
                          //                           (1, provided)
    int iflag;            //Function evaluation output flag
    int nfev = 0;         //Number of function evaluations
    int convflag=0;       //binary flag for solution convergence
    int tempint;          //temporary integer
    int limset=0;         //flag to indicate a boundary constraint exists
    
    double *fcost;        //Actual cost function
    double *mcost;        //Modeled cost function
    double *fcost_old;    //Old actual cost function
    double *mcost_old;    //Old modeled cost function

    double *actual_reduction;    //Actual cost fucntion reduction
    double *pred_reduction;      //Model cost function reduction

    double *rhoc_numer;
    double *rhoc_denom;
    double *rho_f;   //ratio of actual and model reduction

    double *delta; //trust region.

    //Allocate cost function values
    fcost = calloc(1,sizeof(double));
    fcost_old = calloc(1,sizeof(double));
    mcost = calloc(1,sizeof(double));
    mcost_old = calloc(1,sizeof(double));
    delta = calloc(1,sizeof(double));

    pred_reduction = calloc(1,sizeof(double));
    actual_reduction = calloc(1,sizeof(double));
    rhoc_denom = calloc(1,sizeof(double));
    rhoc_numer = calloc(1,sizeof(double));
    rho_f = calloc(1,sizeof(double));

    *fcost = DBL_MAX;

    double tempd0;   //temporary double value;
    double tempd1;   //temporary double value;
    double tempd2;   //temporary double value
    
    double *yhat;       //function evaluation guess
    double *yold;       //previous function evaluation
    double *xnew;       //array holding updated parameter values

    int *side=0;
    double *step = 0;
    double *dstep = 0;
 
    static double **dydpar;    //analytical Jacobian (if supplied)
   
    int *ifree;      //Index array of free parameters
    int *pfixed;     //index of fixed parameters
    int *isllimit;   //Binary array of lower limit set
    int *isulimit;   //Binary array of upper limit set
    double *ulimit;  //Array of lower limit values
    double *llimit;  //Array of upper limit values
    
     
    //Set the default values in the config structure values
    conf.maxiter = 1500;
    conf.eps1 = 1E-8;    //gradient tolerance (gtol)
    conf.eps2 = 1E-8;    //parameter tolerance (xtol)  
    conf.eps3 = 1E-8;    //cost function (ftol)

    //Trust region parameters
    conf.delta0 = 1.0;     //Initial trust region radius
    conf.delta_min = 1E-6; //Min trust region radius
    conf.delta_max = 100.; //Max trust region radius
    conf.gamma1 = 0.25;    //scale factor to decrease trust region
    conf.gamma2 = 2.0;     //scale factor to increase trust region
    conf.mu = 0.25;        //Trust region ratio acceptance threshold
    conf.eta = 0.75;       //trust region ratio increase threshold
    conf.theta = 0.95;     //Threhold for stepback from boundary

    conf.prnt = 0;
    conf.step_size = DBL_EPSILON;
    
    /*Determine if the user has input configuration parameters and, if so,*/
    /*update the conf structure*/
    if (config) {
        if (config->maxiter > 0) {
            conf.maxiter = config->maxiter;
        }
        if (config->eps2 > 0) {
            conf.eps1 = config->eps1;
        }
        if (config->eps2 > 0) {
            conf.eps2 = config->eps2;
        }
        if (config->eps3 > 0) {
            conf.eps3 = config->eps3;
        }
        if (config->prnt> 0) {
            conf.prnt = config->prnt;
        }
        if (config->step_size> 0) {
            conf.step_size = config->step_size;
        }
        if (config->delta0> 0) {
            conf.delta0 = config->delta0;
        }
        if (config->delta_min> 0) {
            conf.delta_min = config->delta_min;
        }
        if (config->delta_max> 0) {
            conf.delta_max = config->delta_max;
        }
        if (config->gamma1> 0) {
            conf.gamma1 = config->gamma1;
        }
        if (config->gamma2> 0) {
            conf.gamma2 = config->gamma2;
        }
        if (config->mu> 0) {
            conf.mu = config->mu;
        }
        if (config->eta> 0) {
            conf.eta = config->eta;
        }
        if (config->theta> 0) {
            conf.theta = config->theta;
        }
    }
    
    if (result) {
        result->niter = -999;
        result->nfev = -999;
        result->cost =  -999;
        result->costa = -999;
        result->costm = -999;
        result->status = -999;
        result->redX2 = -999;
    }

    /*---------------------*/
    /* Identify any fixed parameters and set flags in ifree array*/  
    ifree = (int *) calloc(npar, sizeof (int));
    pfixed = (int *) calloc (npar, sizeof (int));
    nfree = 0;    

    if (pars) for (i=0; i<npar; i++) {
        pfixed[i] = (pars[i].fixed)?1:0;
    }

    if (pars) {
        for (i = 0, j=0; i < npar; i++) {
                if (pfixed[i] == 0) {      
                    ifree[j++] = i;
                    nfree++; 
                }
            }    
    } else { 
        nfree = npar;
        for (i = 0; i < npar; i++) {
            ifree[i] = i;
        }
    }

    if (nfree == 0){
        status = TROPT_ERR_NFREE;
        fprintf(stderr, "%s", "-E- Number of free parameters set to zero.\n");
        goto cleanup_main;
    }

    //Check if upper/lower limit constraints and set
    isllimit = (int *)    calloc(nfree, sizeof (int));
    isulimit = (int *)    calloc(nfree, sizeof (int));
    llimit   = (double *) calloc(nfree, sizeof (double));
    ulimit   = (double *) calloc(nfree, sizeof (double));
    
    if (pars) {
        for (i=0; i<nfree;i++) {
            isllimit[i] = 0;
            isulimit[i] = 0;
            if (pars[ifree[i]].limited[0])  {
                isllimit[i] = 1;
                limset=1;
                if (isfinite(pars[ifree[i]].limits[0])) {
                    //user-supplied lower limit
                    llimit[i] = pars[ifree[i]].limits[0];  
                } else {
                    //If no user supplied limit, lower bound
                   printf("\n%s%s%d\n", __FILE__ ,": -E- No lower-limit set for bounded par #:", i);
                   exit(TROPT_NO_LLIMSET);
                }
            }             
            if (pars[ifree[i]].limited[1] )  {
                isulimit[i] = 1;
                limset=1;
                if (isfinite(pars[ifree[i]].limits[1])) {
                    //user-supplied upper limit
                    ulimit[i] = pars[ifree[i]].limits[1];
                } else {
                    //If no user supplied limit, upper bound
                    printf("\n%s%s%d\n", __FILE__ ,": -E- No upper-limit set for bounded par #:", i);
                    exit(TROPT_NO_ULIMSET);
                   
                }
            }
        }
    }
    
    //Check for Jacobian derivative parameters    
    dstep = (double *) calloc(nfree, sizeof (double));
    step =  (double *) calloc(nfree, sizeof (double));
    side =  (int *)    calloc(nfree, sizeof (int));

    if (pars) {
        for (i = 0; i < npar; i++) {
            step[ifree[i]] = pars[ifree[i]].step;
            dstep[i] = pars[ifree[i]].relstep;
            side[i] = pars[ifree[i]].side;
        }
    }
    
    /*Allocate GSL matrices and vectors*/
    gsl_matrix *covym = gsl_matrix_calloc(m,m);

    //Covariance matrix of just the free parameters
    gsl_matrix *Shat = gsl_matrix_calloc(nfree, nfree);
    gsl_vector *delta_y = gsl_vector_calloc(m);
    gsl_vector *xstep = gsl_vector_calloc(nfree);
    gsl_vector *xstep_hat = gsl_vector_calloc(nfree);

    //Allocate Jacobian matrix of the function (df/dx) and set elements to zero
    gsl_matrix *K = gsl_matrix_calloc(m,nfree);
    gsl_matrix *K_hat = gsl_matrix_calloc(m,nfree);

    //Allocate the gradient of the cost function
    gsl_vector *gradf = gsl_vector_calloc(nfree);
    gsl_vector *gradfhat = gsl_vector_calloc(nfree);

    //Allocate Hessian matrix for the cost function
    gsl_matrix *H = gsl_matrix_calloc(nfree,nfree);
    gsl_matrix *H_hat = gsl_matrix_calloc(nfree,nfree);

    //Define the affine scaling matrix, D, and inverse scaling matrix
    //Allocate later if boundary constraints set (limset =1)
    gsl_matrix *invDc = gsl_matrix_calloc(nfree,nfree);
    gsl_matrix *Dc = gsl_matrix_calloc(nfree,nfree);
    gsl_matrix *invD = gsl_matrix_calloc(nfree,nfree);
    gsl_matrix *D = gsl_matrix_calloc(nfree,nfree);
    gsl_matrix *invDaff = gsl_matrix_calloc(nfree,nfree);
    gsl_matrix *Daff = gsl_matrix_calloc(nfree,nfree);
    gsl_matrix *Chat = gsl_matrix_calloc(nfree,nfree);
    gsl_matrix *C = gsl_matrix_calloc(nfree,nfree);
    gsl_matrix *Mhat = gsl_matrix_calloc(nfree,nfree);

    //Inverse covariance matrices
    gsl_matrix *Wy = gsl_matrix_calloc(m,m);
    //Scaled inverse covariance matrix
    gsl_matrix *What;
    
    //Working matrices
    gsl_matrix *W0 = gsl_matrix_calloc(nfree,nfree);
    gsl_matrix *W1 = gsl_matrix_calloc(nfree,nfree);
    //gsl_matrix *W2 = gsl_matrix_calloc(nfree,m);
    gsl_matrix *W2 = gsl_matrix_calloc(m,nfree);

    //Working vectors
    gsl_vector *tempy = gsl_vector_calloc(m);
    gsl_vector *tempv1 = gsl_vector_calloc(m);
    gsl_vector *tempv2 = gsl_vector_calloc(nfree);
    gsl_vector *tempv3 = gsl_vector_calloc(nfree);
    
    if (nfree <= m) {
        What = gsl_matrix_calloc(nfree,nfree);
    } else {
        What = gsl_matrix_calloc(m,m); 
    }
    
    //Forward model parameters
    yhat =   (double *) calloc(m, sizeof (double));
    yold =   (double *) calloc(m, sizeof (double));
    xnew =   (double *) calloc(npar, sizeof (double));

    //Forward model Jacobian matrix as 1-D array
    dydpar = (double **) calloc(m, sizeof (double *));
    for (i = 0; i < m; i++) {
        dydpar[i] = (double *) calloc(npar, sizeof (double));
    }


    /*==============*/
    /*Begin checks for initialization errors*/

    if (covydiag != 0) {
        if(covydiag !=1) {
            status = TROPT_COVYFLAG_ERR;
            goto cleanup_main;
        }
    }   

    //Check if input function has been provided
    if (funct == 0) {
        status = TROPT_ERR_NOFUNC;
        goto cleanup_main;
    }
    //Check if number of parameters is greater than zero
    if (npar <= 0) {
        status = TROPT_ERR_NPAR;
        goto cleanup_main;
    }
    //Check if number of free parameters is greater than zero
    if (nfree <= 0) {
        status = TROPT_ERR_NFREE;
        goto cleanup_main;
    }
    //Check adequate degrees of freedom exist
    if (m < nfree) {
        status = TROPT_ERR_DOF;
        goto cleanup_main;
    }
    // Sanity checking on input configuration
    if ((conf.maxiter <= 0) || 
        (conf.eps1 <= 0) || 
        (conf.eps2 <= 0) ||
        (conf.eps3 <= 0) || 
        (conf.delta0 <= 0) ||
        (conf.mu <=0) ||
        (conf.theta <=0) ||
        (conf.eta <= 0) ||
        (conf.gamma1 <= 0) ||
        (conf.gamma2 <= 0) ||
        (conf.prnt < 0))
        {
        status = TROPT_ERR_PARAM;
        goto cleanup_main;
    }
    //End initialization error checks
    /*===============*/


    //Populate gsl covym covariance matrices if C array provided
    if (covy != NULL){
        covyflag = 1;
        for (i = 0; i < m; i++) {
            for (j = 0; j < m; j++) {        
                gsl_matrix_set(covym, i, j, covy[i][j]);
            }
        }       
    } else {
        //No measurement uncertainty covariance provided
        covyflag = 0;
    }
    
    //Compute weight matrix for cost function. This is the inverse of the 
    //covariance matrix. However, if covy not supplied set it to and Identity matrix.
    if (covyflag) {
        if (covydiag == 1) {
            iflag = tropt_matrix_invert(covym, Wy, 1);
        } else {
            iflag = tropt_matrix_invert(covym, Wy, 0);  
        }
    } else {
        //Set Weight matrix to identidy matrix
        for (i = 0; i < m; i++) {
            for (j = 0; j < m; j++) {
                if (i == j) {
                    gsl_matrix_set(Wy, i, j, 1.0);
                } else {
                    gsl_matrix_set(Wy, i, j, 0.0);
                }
            }
        }
    }
    
    
    /*===============*/
    //First function evaluation
    //Initialize function evalution counter;

    convflag = 0;

    // Evaluate user function with initial parameter values
    iflag = tropt_call(funct, m, npar, xall, yhat, dydpar, private_data);
    nfev += 1;

    //Initial Jacobian matrix
    iflag = tropt_jac(funct, m, nfree, ifree, npar, xall, K,
                conf.step_size, private_data, &nfev,
                step, dstep, side, dydpar); 
    
    /* Make a new copy */
    //set initial residul vector: deltay = y - yhat
    for (i = 0; i < m; i++) {
        gsl_vector_set(delta_y,i, yall[i] - yhat[i]); 
    }

    //First gradient of the cost function
    calc_tropt_gradf(Wy, delta_y, tempv1, K, gradf);

    //Set ellipsoid scaling matrix (for conditioning)
    set_scaling_matrix(gradf, D, invD, nfree, 0);

    if (limset) {
        //Calculate the affine scaling matrix (for boundary constraints)
        set_affine_scaling_matrix(isllimit, isulimit, 
        llimit, ulimit, xall, gradf, nfree, ifree, tempv2,tempv3 );

        for (i=0;i<nfree;i++) {
            if (gsl_vector_get(tempv3,i) !=0 ) {
                tempd0 = sqrt(fabs(gsl_vector_get(tempv2,i)));
                gsl_matrix_set(Daff,i,i,1./tempd0);
                gsl_matrix_set(invDaff,i,i,tempd0);
            }
        }

        //Compute initial scaling matrix invDc 
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, invDaff, invD, 0.0, W1);
        gsl_matrix_memcpy(invDc,W1);

        //Compute initial combined scaling matrix Dc
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Daff, D, 0.0, W1);
        gsl_matrix_memcpy(Dc,W1);

    }

    //Compute initial value of the cost function
    iflag =  tropt_cost_actual(delta_y, tempv1, Wy, covyflag, fcost_old);  
    
    //Initialize xnew array, set to initial guess
    for (i = 0; i<nfree; i++) {
        xnew[ifree[i]] = xall[ifree[i]];
    }

    //Initialise actual vs predicted reduction ratio
    *rho_f  = 1.;   

    //Test if the maximum value of the gradient of the cost function is less 
    //than tolerance eps3.
    //If so, the initial guess x is very close, return x unchanged as the solution.
    if ( gsl_blas_dnrm2(gradf) < conf.eps1) {
        status = TROPT_SUCCESS_CLOSE;
        convflag = 1;
        goto cleanup_main;
    }
    
    //Set the initial trust region radius
    *delta = conf.delta0;

    /*===============*/
    //Begin Trust region parameter and evaluation counter  
    iter = 0;

    //solve uncostrained trust region problem

    if (!limset) {
        
        while ((!convflag) && (iter <= conf.maxiter) ) {

            //Compute first scaled jacobian
            gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, K, invD, 0.0, K_hat);
        
            //First scaled Hessian
            calc_hessian(Wy, W2, K_hat, H);

            //First scaled gradient of the cost function
            calc_tropt_gradf(Wy, delta_y, tempv1, K_hat, gradf);

            //Solve the sub-problem in hat space, get xstep_hat
            solve_subproblem(gradf, H, xstep_hat, delta, nfree);

            //Scale xstep_hat -> xstep
            //xstep = D xstep_hat
            gsl_blas_dgemv(CblasNoTrans, 1.0, D, xstep_hat, 0.0, xstep);

            //Check step size is finite
            tempint = 0;
            for (i = 0; i < nfree; i++) {
                tempint = isfinite(gsl_vector_get(xstep,i));
                //take action if non-finite step-size detected, nudge it.
                if (tempint <= 0) {
                    gsl_vector_set(xstep,i,DBL_EPSILON);
                }           
            }

            //Update the state vector with step
            for (i = 0; i < nfree; i++) {
                xnew[ifree[i]] = xall[ifree[i]] + gsl_vector_get(xstep,i);
            }

            //Compute the forward model at xnew 
            iflag = tropt_call(funct, m, npar, xnew, yhat, dydpar, private_data);
            nfev +=1;


            //Compute difference between observation of forward model
            for (i = 0; i < m; i++) {
                gsl_vector_set(delta_y,i, yall[i] - yhat[i]);    
            }

            //First gradient of the cost function
            iflag = calc_tropt_gradf(Wy, delta_y, tempv1, K_hat, gradf);

            //Compute the actual cost function
            iflag = tropt_cost_actual(delta_y, tempv1, Wy, covyflag, fcost);
            
            //Compute estimated cost function
            iflag = calc_cost_model(xstep_hat, gradfhat, H, tempv2, pred_reduction);

            //rhof_numer = f(x + xstep) - f(x)
            *actual_reduction = *fcost_old - *fcost;

            //ratio of actual reduction to predicted reduction
            *rho_f = *actual_reduction / -(*pred_reduction); 

            if (isnan(*fcost)) {
                status = TROPT_COSTNAN_ERR; 
                convflag = -1;
                goto cleanup_main;
            };

            //Norm of the step length
            tempd0 = gsl_blas_dnrm2(gradf);
            tempd1 = gsl_blas_dnrm2(xstep);
            tempd2 =  tropt_norm2_xall(xnew,ifree,nfree);

            //Evaluate and updte the trust region
            iflag = update_trust_region(conf.mu, conf.eta, conf.theta, 
                            conf.gamma1, conf.gamma2, 
                            actual_reduction, pred_reduction, 
                            tempd1, rho_f, delta);

            //test step tolerance xtol
            //printf("norm_xstep = %f, xtol thresh = %.3e\n",tempd1,conf.eps2*(conf.eps2 + tempd2));
            if (tempd1 < conf.eps2*(conf.eps2 + tempd2)) {
                status = TROPT_SUCCESS_XTOL;
                convflag = 1;
                //goto cleanup_main;
            }

            //Check cost function tolerance (ftol)
            //printf("df = %f, f = %.3e\n",*actual_reduction,1E-8 * (*fcost));
            if (*actual_reduction < conf.eps3* (*fcost)) {
                status = TROPT_SUCCESS_FTOL;
                convflag = 1;
                //goto cleanup_main;
            } 

            //Cost function gradient tolerance (gtol)
            //printf("grad_tol = %.3e, gtol_thres = %.3e\n",tempd0, conf.eps1);
            if ((tempd0 < conf.eps1) ) {
                status = TROPT_SUCCESS_JGRAD;
                convflag = 1;
                //goto cleanup_main;
            }

            //Test for convergence/termination
            //Check to see if max iterations were reached
            if (iter == conf.maxiter){
                status = TROPT_MAXITER;
                convflag = 1;
                //goto cleanup_main;
            }

            if (*rho_f > 0.25) {

                if(convflag) {
                    goto cleanup_main;
                }

                for (i = 0; i < nfree; i++) {
                    //accept step and update state vector
                    xall[ifree[i]] += gsl_vector_get(xstep,i);  
                } 

                *fcost_old = *fcost; 

                iflag = tropt_jac(funct, m, nfree, ifree, npar, xall, K,
                conf.step_size, private_data, &nfev,
                step, dstep, side, dydpar);

                //Update the caling matrix D for iteration > 0
                set_scaling_matrix(gradf, D, invD, nfree, iter);

            }
            //update iteration counter
            iter++;



        } //end while
    } //end if not limset


    //If any boundary constraints set
    if (limset) {
        
        //Begin main iteration loop
        //while ( (!convflag) && (iter <= conf.maxiter) ) {
        while ( iter <= conf.maxiter ) {

            //K_hat = K invDc, combained parameter and affine scaling
            gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, K, invDc, 0.0, K_hat);
            
            //Scale Hessian H -> H_hat (update in place)
            calc_hessian(Wy, W2, K_hat, H);

            //Gradient of the cost function gradf (in place) 
            calc_tropt_gradf(Wy, delta_y, tempv1, K, gradf);

            //Scaled gradient of the cost function gradfhat (in place) 
            calc_tropt_gradf(Wy, delta_y, tempv1, K_hat, gradfhat);

            //Calculate matrix Mhat and Chat
            calc_mhat(gradf, invD, H, nfree, Mhat, Chat);

            //If limset, subproblem is solved in hat space
            solve_subproblem(gradfhat, Mhat, xstep_hat, delta, nfree);

            //Scale xstep_hat -> xstep
            //xstep = Dc xstep_hat, where Dc = Daff D
            gsl_blas_dgemv(CblasNoTrans, 1.0, Dc, xstep_hat, 0.0, xstep);
            

            //Check if the step size if finite
            tempint = 0;
            for (i = 0; i < nfree; i++) {
                tempint = isfinite(gsl_vector_get(xstep,i));

                //take action if non-finite step-size detected, nudge it.
                if (!tempint) {
                    gsl_vector_set(xstep,i,DBL_EPSILON);
                }           
            }

            select_step(xall, xnew, llimit, ulimit, ifree, 
                        isllimit, isulimit, nfree, m, 
                        K_hat, Mhat, Chat, gradfhat, 
                        xstep, xstep_hat, invDc, Dc, 
                        delta, conf.theta, pred_reduction);

            //Make step strictly feasible
            /*
            for (i=0;i<nfree;i++){

                xnew[ifree[i]] = xall[ifree[i]] + gsl_vector_get(xstep,i); 

                if (isllimit[ifree[i]]) {
                    if (xnew[ifree[i]] < llimit[ifree[i]] ) {
                        xnew[ifree[i]] = 0.95*llimit[ifree[i]]; 
                    }
                }
                if (isulimit[ifree[i]]) {
                    if (xnew[ifree[i]] > ulimit[ifree[i]] ) {
                        xnew[ifree[i]] = 0.95*ulimit[ifree[i]] ;

                    }
                }

            }*/
            
             
                            
            //Call forward model with updated x parameters 
            iflag = tropt_call(funct, m, npar, xnew, yhat, dydpar, private_data);
                        
            //Compute estimated cost function
            iflag = calc_cost_model(xstep_hat, gradfhat, H, tempv2, pred_reduction);

            for (i = 0; i < m; i++) {
                gsl_vector_set(delta_y,i, yall[i] - yhat[i]); 
            }

            iflag = tropt_cost_actual(delta_y, tempv1, Wy, covyflag, fcost);

            *actual_reduction =  *fcost -*fcost_old ;
            *rho_f = *actual_reduction / -(*pred_reduction);

            //if (*pred_reduction < 0) {
            //    *rho_f = *actual_reduction / -(*pred_reduction);
            //} else if (*pred_reduction == *actual_reduction) {
            //    *rho_f = 1;
            //} else {
            //    *rho_f = 0;
            //}

            //printf("actual_reduct = %.3e\n",*actual_reduction);
            //printf("pred_reduct = %.3e\n",*pred_reduction);
            //printf("rho_f = %.3e\n",*rho_f);
            //printf("iter = %d\n",iter);
            
            if (isnan(*fcost)) {
                status = TROPT_COSTNAN_ERR; 
                convflag = -1;
                goto cleanup_main;
            };

            //Norm of the step length
            tempd0 = gsl_blas_dnrm2(gradf);
            tempd1 = gsl_blas_dnrm2(xstep);
            tempd2 =  tropt_norm2_xall(xnew,ifree,nfree);

            //Evaluate and updte the trust region
            iflag = update_trust_region(conf.mu, conf.eta, conf.theta, 
                            conf.gamma1, conf.gamma2, 
                            actual_reduction, pred_reduction, 
                            tempd1, rho_f, delta);

            //Convergence tests here. 
            //test step tolerance xtol
            //if (tempd1 < conf.eps2*(conf.eps2 + tempd2)) {
            if (gsl_blas_dnrm2(xstep) < conf.eps2) {
                status = TROPT_SUCCESS_XTOL;
                convflag = 1;
                //goto cleanup_main;
            }

            //Check cost function tolerance (ftol)
            //if (fabs(*actual_reduction) < conf.eps3 * (*fcost)) {
            if (*fcost < conf.eps3 ) {
            //if (fabs(*actual_reduction) < conf.eps3 ) {
                status = TROPT_SUCCESS_FTOL;
                convflag = 1;
                //goto cleanup_main;
            } 

            //Cost function gradient tolerance (gtol)
            if (tempd0 < conf.eps1) {
                status = TROPT_SUCCESS_JGRAD;
                convflag = 1;
                //goto cleanup_main;
            }

            //Test for convergence/termination
            //Check to see if max iterations were reached
            if (iter == conf.maxiter){
                status = TROPT_MAXITER;
                convflag = 1;
                //goto cleanup_main;
            }


            if (actual_reduction > 0) {
                for (i = 0; i < nfree; i++) {
                    //accept step and update state vector
                    //printf("xstep = %.3e\n",gsl_vector_get(xstep,i));
                    xall[ifree[i]] = xnew[ifree[i]]; 
                    //printf("proposed step = %.3e\n",gsl_vector_get(xstep,i));
                } 
                //printf("\n");
                if((convflag) && (*rho_f > 0.25)) {
                    //printf("OUT STATUS = %d\n",status);
                    goto cleanup_main;
                }

                *fcost_old = *fcost; 

            }

            iflag = tropt_jac(funct, m, nfree, ifree, npar, xnew, K,
            conf.step_size, private_data, &nfev,
            step, dstep, side, dydpar);

            //First gradient of the cost function
            calc_tropt_gradf(Wy, delta_y, tempv1, K, gradf);

            //Set ellipsoid scaling matrix (for conditioning)
            set_scaling_matrix(gradf, D, invD, nfree, iter);

            //Calculate the affine scaling matrix (for boundary constraints)
            set_affine_scaling_matrix(isllimit, isulimit, 
            llimit, ulimit, xall, gradf, nfree, ifree, tempv2,tempv3 );

            for (i=0;i<nfree;i++) {
                if (gsl_vector_get(tempv3,i) !=0 ) {
                    tempd0 = sqrt(fabs(gsl_vector_get(tempv2,i)));
                    gsl_matrix_set(Daff,i,i,1./tempd0);
                    gsl_matrix_set(invDaff,i,i,tempd0);
                }
            }

            //Compute initial scaling matrix invDc 
            gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, invDaff, invD, 0.0, W1);
            gsl_matrix_memcpy(invDc,W1);

            //Compute initial combined scaling matrix Dc
            gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Daff, D, 0.0, W1);
            gsl_matrix_memcpy(Dc,W1);
            
            //update iteration counter
            iter++;


        } //End limset outer while
 
        
    } //end if limset
    

   
   /*===============*/
   //Begin cleanup and report back info integer flag;
    
    cleanup_main:  
        //If error encountered, return 1
        if (convflag < 0) {
            iflag = 1;
        } else {
            iflag = 0;
        }

        //print outputs to info struct
        result->niter = iter;
        result->nfev = nfev;
        result->cost = *fcost;
        result->costm = *mcost;
        result->status = status;
        result->redX2 = *fcost/(double)(m-nfree+1);
        result->status = convflag;

        if (ifree) free(ifree);
        if (pfixed) free(pfixed);
        if (llimit) free(llimit);
        if (ulimit) free (ulimit);
        if (isulimit) free(isulimit);
        if (isllimit) free(isllimit);
        if (side) free(side);
        if (step) free(step);
        if (dstep) free(dstep);

        if(yhat) free(yhat);
        if(yold) free(yold); 
        if(xnew) free(xnew);  

        if(dydpar) tropt_free2D(dydpar,m);

        if(fcost) free(fcost);
        if(mcost) free(mcost);
        if(fcost_old) free(fcost_old); 
        if(mcost_old) free(mcost_old); 
        if(pred_reduction) free(pred_reduction); 
        if(actual_reduction) free(actual_reduction); 
        if(delta) free(delta);

        if (rho_f) free (rho_f);
        if (rhoc_numer) free(rhoc_numer);
        if (rhoc_denom) free(rhoc_denom);

        if(covym) gsl_matrix_free(covym);
        if (K)  gsl_matrix_free(K);
        if (K_hat) gsl_matrix_free(K_hat);
        if (H) gsl_matrix_free(H);
        if (H_hat) gsl_matrix_free(H_hat);
        if (Wy) gsl_matrix_free(Wy);
        if (W0) gsl_matrix_free(W0);
        if (W1) gsl_matrix_free(W1);
        if (W2) gsl_matrix_free(W2);
        if(What) gsl_matrix_free(What);
        if(Shat) gsl_matrix_free(Shat);
        if(C) gsl_matrix_free(C);
        if(Chat) gsl_matrix_free(Chat);
        if(Mhat) gsl_matrix_free(Mhat);
        
        if(D) gsl_matrix_free(D);
        if(invD) gsl_matrix_free(invD);
        if(Dc) gsl_matrix_free(Dc);
        if(invDc) gsl_matrix_free(invDc);
        if(Daff) gsl_matrix_free(Daff);
        if(invDaff) gsl_matrix_free(invDaff);

        if (tempy) gsl_vector_free(tempy);
        if (tempv1) gsl_vector_free(tempv1);
        if (tempv2) gsl_vector_free(tempv2);
        if (tempv3) gsl_vector_free(tempv3);
        if (delta_y) gsl_vector_free(delta_y);
        if (xstep) gsl_vector_free(xstep);
        if (xstep_hat) gsl_vector_free(xstep_hat);
        if (gradf) gsl_vector_free(gradf);
        if (gradfhat) gsl_vector_free(gradfhat);

        return iflag;
}


/**************************************************************/
/*                     select_step                            */
/*                                                            */
/* Select subproblem step (s), reflected step (r), or gradient*/
/*  step  (g)                                                 */
/*                                                            */
/* Method follows that used in Python Scipy library           */
/* */


static int select_step(double *xall, double *xnew, double *llimit, double *ulimit,
            int *ifree, int *isllimit, int *isulimit, int nfree, int m,
            gsl_matrix *K_hat, gsl_matrix *Hhat, gsl_matrix *C, 
            gsl_vector *gradfhat,gsl_vector *xstep, gsl_vector *xstep_hat, 
            gsl_matrix *invD, gsl_matrix *D, 
            double *delta, double theta, double *mcost) {

    int i, status, outflag;
    double step_to_bound, step_to_delta;
    double ag_stride, r_stride, p_stride;
    double p_value, r_value, ag_value;
    double tempd0,tempd1,tempd2,tempd3;
    double r_llimit_stride, r_ulimit_stride;

    gsl_vector *tempv2 =     gsl_vector_alloc(nfree);
    gsl_vector *tempv3 =     gsl_vector_alloc(nfree);
    gsl_vector *rxstep_hat = gsl_vector_alloc(nfree);
    gsl_vector *rxstep =     gsl_vector_alloc(nfree);
    gsl_vector *x_on_bound = gsl_vector_alloc(nfree);
    gsl_vector *gxstep_hat = gsl_vector_alloc(nfree);
    gsl_vector *gxstep =     gsl_vector_alloc(nfree);

    status = check_is_bounded(xstep, xall, ifree, nfree, 
                            ulimit, llimit, xnew, tempv2);

    if (status) {

        calc_cost_model(xstep_hat, gradfhat, Hhat, tempv2, mcost);

        gsl_vector_free(tempv2);
        gsl_vector_free(tempv3);
        gsl_vector_free(rxstep);
        gsl_vector_free(rxstep_hat);
        gsl_vector_free(gxstep);
        gsl_vector_free(gxstep_hat);
        gsl_vector_free(x_on_bound);
        
        return 0;

    } else { 

        //What is the step to the boundary
        p_stride = calc_alpha(xall, ulimit, llimit, isulimit, 
                isllimit, ifree, nfree, xstep, tempv2);
        
        //If we hit the boundary, change sign of the 
        //gradient of the cost function
        for (i = 0; i<nfree; i++) {      
            tempd0 = gsl_vector_get(xstep_hat,i);
            if(gsl_vector_get(tempv2,i)){ 
                gsl_vector_set(rxstep_hat,i,-1*tempd0);
            } else { 
                gsl_vector_set(rxstep_hat,i,tempd0);
            }
        }

        //Scale by D to get reflected step in parameter space
        gsl_blas_dgemv(CblasNoTrans, 1.0, D, rxstep_hat, 0.0, rxstep);

        //Scale step to hit bound
        gsl_vector_scale(xstep, p_stride);
        gsl_vector_scale(xstep_hat, p_stride);

        //Temporary solution that is on boundary
        for (i = 0; i<nfree; i++) {      
            //Populate the the x_on_bound vector
            tempd0 = xall[ifree[i]] + gsl_vector_get(xstep,i);
            xnew[ifree[i]] = xall[ifree[i]] + gsl_vector_get(xstep,i);
            gsl_vector_set(x_on_bound,i,tempd0);
        }
        
        //Compute reflected scaling - which is better?
        //(i) To upper/lower limit, or 
        //(ii) or to the trust region boundary?
        //Compute distance along reflected path to trust region
        step_to_delta = calc_lnsph_inter(xstep_hat, rxstep_hat, delta);

        //Compute distance path to trust boundary
        step_to_bound = calc_alpha(xnew, ulimit, llimit, isulimit, 
                isllimit, ifree, nfree, rxstep, tempv2);

        //Compute the smallest scaling step
        r_stride = tropt_min(step_to_delta,step_to_bound);
        
        
        if (r_stride > 0) {
            r_llimit_stride = (1. - theta) * p_stride/r_stride; 
            if (r_stride == step_to_bound) {
                r_ulimit_stride = theta*step_to_bound; 
            } else {
                r_ulimit_stride = step_to_delta;
            }
        } else {
            r_llimit_stride= 0;
            r_ulimit_stride = -1;
        }

        //reflected scaling
        if (r_llimit_stride <= r_ulimit_stride) {

            calc_quadratic_1d(Hhat, gradfhat, rxstep_hat, xstep_hat, NULL, nfree, m,
                        &tempd0, &tempd1, &tempd2);

            eval_quadratic_1d(tempd0, tempd1, tempd2, r_llimit_stride, r_ulimit_stride, 
                                &r_stride, &r_value);

            gsl_vector_scale(rxstep_hat, r_stride);
            gsl_vector_add(rxstep_hat,xstep_hat);
            gsl_blas_dgemv(CblasNoTrans, 1.0, D, rxstep_hat, 0.0, rxstep);

        } else {
            r_value = DBL_MAX;
        }

        //Scale by theta to make step interior
        gsl_vector_scale(xstep,theta);
        gsl_vector_scale(xstep_hat,theta);

        //Compute cost for interior xstep
        calc_cost_model(xstep_hat, gradfhat, Hhat, tempv2, &p_value); 

        //Compute reflected gradient rgrad
        gsl_vector_memcpy(gxstep_hat,gradfhat);
        gsl_vector_scale(gxstep_hat,-1.0);
        gsl_blas_dgemv(CblasNoTrans, 1.0, D, gxstep_hat, 0.0, gxstep);

        step_to_delta= *delta/gsl_blas_dnrm2(gxstep_hat);

        step_to_bound = calc_alpha(xall, ulimit, llimit, isulimit, 
                isllimit, ifree, nfree, gxstep, tempv2);

        if (step_to_bound < step_to_delta) {
                ag_stride = theta*step_to_bound;
            } else {
                ag_stride = step_to_delta;
        }

        //Set gxstep_hat and gxstep
        gsl_vector_scale(gxstep_hat,ag_stride);
        calc_cost_model(gxstep_hat, gradfhat, Hhat, tempv2, &ag_value); 
        calc_quadratic_1d(Hhat, gradfhat, gxstep_hat, NULL, NULL, nfree, m,
                    &tempd0, &tempd1, &tempd2);

        //update ag_stride
        eval_quadratic_1d(tempd0, tempd1, 0, 0, ag_stride, 
                        &tempd3, &ag_value);

        //If here is zero reduction, set to max
        gsl_vector_scale(gxstep,ag_stride);
        gsl_vector_scale(gxstep_hat,ag_stride);

        calc_cost_model(gxstep_hat, gradfhat, Hhat, tempv2, &ag_value);      

        if ((p_value < r_value) && (p_value < ag_value)) {
            *mcost = p_value;
            outflag = 1;
        } else if ((r_value < p_value) && (r_value < ag_value)) { 
            gsl_vector_memcpy(xstep,rxstep);
            gsl_vector_memcpy(xstep_hat , rxstep_hat);
            *mcost = r_value;
            outflag=2;
        } else {
            gsl_vector_memcpy(xstep,gxstep);
            gsl_vector_memcpy(xstep_hat , gxstep_hat);
            *mcost = ag_value;
            outflag=3;
        }    

        gsl_vector_free(tempv2);
        gsl_vector_free(tempv3);
        gsl_vector_free(rxstep);
        gsl_vector_free(rxstep_hat);
        gsl_vector_free(gxstep);
        gsl_vector_free(gxstep_hat);
        gsl_vector_free(x_on_bound);

        return outflag;
    }     
}


/**************************************************************/
/*                     calc_s                                 */
/*                                                            */
/* Calculate step, s, as:                                     */
/*                                                            */
/*                     -1                                     */
/* s = -( H + lambda I)   gradf                                */
/*                                                            */
/* We use eigevalue decompositon to express the Hessian matrix*/
/* as:                                                        */
/*                    T                                      */
/*  H = -Q diag(egval) Q                                       */
/*                                                            */
/* as rexpress s as:                                          */
/*                                  -1   T                    */
/* s =  -Q (diag(E) + diag(lambda))   Q     gradf             */
/*                                                            */
/*                                                            */
/*Input:                                                      */
/* evec: gsl_matrix (n,n), eigenvectors of Hessian            */
/* eval: gsl_vectpor (n), eigenvalues of Hessian              */
/* Qtg: gsl_vetor (n), precomputed QTg                        */
/* n: int, number free parameters                             */
/* lambda: double, L-M scaling parameter                      */
/*                                                            */
/* Output:                                                    */
/* outflag: integer                                           */

static int calc_s(gsl_matrix *evec, gsl_vector *eval, gsl_vector *Qtg, 
                int n, double *lambda, gsl_vector *s) {


    int i,j;
    double diag;
    
    //Allocate memory
    gsl_matrix *invD = gsl_matrix_calloc(n,n); 
    gsl_vector *temp = gsl_vector_calloc(n);

    //Create the inverse diagonal matrix invD
    for (i=0;i<n;i++){
        for (j=0;j<n;j++) { 
            if (j == i) {
                if (*lambda == gsl_vector_get(eval,i)) {
                    diag = gsl_vector_get(eval,i);
                } else {
                    diag = gsl_vector_get(eval,i) + (*lambda);
                }
                gsl_matrix_set(invD,i,j,1./diag);
            } else {
                gsl_matrix_set(invD,i,j,0.0);  
            }
        }
    }

    //                                   -1  T
    //Compute -Q (diag(E) + diag(lambda))   Q g
    gsl_blas_dgemv(CblasNoTrans, 1.0, invD, Qtg, 0.0, temp); 
    gsl_blas_dgemv(CblasNoTrans,-1.0, evec, temp, 0.0, s);

    //Cleanup
    gsl_matrix_free(invD);
    gsl_vector_free(temp);

    return 0;
}


/**************************************************************/
/*                        calc_ds                             */
/*                                                            */
/* Calculate the derivative of the step with respect to the   */
/* scaling parameter lambda                                   */
/*                                                            */
/*                                                            */
/*                                  -2   T                    */
/* ds =  Q (diag(E) + diag(lambda))   Q     gradf             */
/*                                                            */ 
/*Input:                                                      */
/* evec: gsl_matrix (n,n), eigenvectors of Hessian            */
/* eval: gsl_vectpor (n), eigenvalues of Hessian              */
/* Qtg: gsl_vetor (n), precomputed QTg                        */
/* n: int, number free parameters                             */
/* lambda: double, L-M scaling parameter                      */
/*                                                            */
/* Output:                                                    */
/* outflag: integer                                           */

static int calc_ds(gsl_matrix *evec, gsl_vector *eval, gsl_vector *Qtg, 
                int n, double *lambda, gsl_vector *ds) {

    int i,j;
    double diag;
    
    //Allocate memory
    gsl_matrix *invD2 = gsl_matrix_calloc(n,n);
    gsl_vector *temp = gsl_vector_calloc(n);

    //Create the inverse diagonal matrix
    for (i=0;i<n;i++){
        for (j=0;j<n;j++) { 
            if (j == i) {
                if (*lambda == gsl_vector_get(eval,i)) {
                    diag = pow( gsl_vector_get(eval,i) , 2.);
                } else {
                    diag = pow( gsl_vector_get(eval,i) + *lambda , 2.);
                }
                gsl_matrix_set(invD2,i,j,1./diag);
            } else {
               gsl_matrix_set(invD2,i,j,0.0); 
            }
        }
    }

    //           -2     T
    //Compute   Q  D   Q g
    gsl_blas_dgemv(CblasNoTrans,  1.0, invD2, Qtg, 0.0, temp);
    gsl_blas_dgemv(CblasNoTrans, 1.0, evec, temp, 0.0, ds);
    
    gsl_matrix_free(invD2);
    gsl_vector_free(temp);

    return 0;
}

/**************************************************************/
/*                     calc_phi                               */
/*                                                            */
/*         phi = ||s|| - delta                                */
/*                                                            */
/* Difference between step and trust region (inverse space)   */
/* phi = 0: step is on the boundary                           */
/* phi > 0: step is outside boundary                          */
/* phi < 0: step is inside boundary                           */
/*                                                            */
/*                                                            */
/*Input:                                                      */
/* evec: gsl_matrix (n,n), eigenvectors of Hessian            */
/* eval: gsl_vectpor (n), eigenvalues of Hessian              */
/* Qtg: gsl_vetor (n), precomputed QTg                        */
/* n: int, number free parameters                             */
/* delta: double, trust region radius                         */
/* lambda: double, L-M scaling parameter                      */
/*                                                            */
/* Output:                                                    */
/* phi: double, distance from step to trust region boundary   */
/*                                                            */

static double calc_phi(gsl_matrix *evec, gsl_vector *eval, 
                gsl_vector *Qtg, int n, double *delta, 
                double *lambda) {

    double norm_s, phi;

    gsl_vector *s  = gsl_vector_alloc(n);

    calc_s(evec, eval, Qtg, n, lambda, s) ;
    norm_s = gsl_blas_dnrm2(s);

    phi = norm_s - *delta;

    gsl_vector_free(s);

    if (norm_s > 0) {
        return phi;   
    } else {
        return INFINITY;
    }
}

/**************************************************************/
/*                        calc_dphi                           */
/*                                                            */
/* Derivative of the difference between step and trust region */
/* with respect to lambda                                     */
/*                                                            */
/*          T                                                 */
/* dphi =  s ds / ||s||                                       */
/*                                                            */
/*Input:                                                      */
/* evec: gsl_matrix (n,n), eigenvectors of Hessian            */
/* eval: gsl_vectpor (n), eigenvalues of Hessian              */
/* Qtg: gsl_vetor (n), precomputed QTg                        */
/* n: int, number free parameters                             */
/* delta: double, trust region radius                         */
/* lambda: double, L-M scaling parameter                      */
/*                                                            */
/* Output:                                                    */
/* dphi: double, dphi/dlambda                                 */
/*                                                            */

static double calc_dphi(gsl_matrix *evec, gsl_vector *eval, 
                gsl_vector *Qtg, int n, double *lambda) {

    double norm_s, dphi; 
    double temp;

    gsl_vector *s  = gsl_vector_alloc(n);
    gsl_vector *ds = gsl_vector_alloc(n);

    calc_s (evec, eval, Qtg, n, lambda, s) ;
    calc_ds(evec, eval, Qtg, n, lambda, ds) ;

    norm_s = gsl_blas_dnrm2(s);
    gsl_blas_ddot(s,ds,&temp);

    dphi = temp/norm_s;

    gsl_vector_free(s);
    gsl_vector_free(ds);

    if (norm_s > 0) {
        return dphi;
    } else {
        return INFINITY;
    }
}

/**************************************************************/
/*                      calc_root_lambda                      */
/*                                                            */

static int calc_lambda( gsl_matrix *evec, gsl_vector *eval,
                    gsl_vector *Qtg, gsl_vector *s, int n,
                    double *delta, double *lambda) {

    double uq,lq;
    double phi,dphi;
    double norm_qtg;
    double lambda0;
    int iter = 0;

    norm_qtg = gsl_blas_dnrm2(Qtg);
    phi =   calc_phi(evec, eval, Qtg, n, delta, lambda);
    dphi =  calc_dphi(evec, eval, Qtg, n, lambda);

    uq = norm_qtg / (*delta);
    lq = -phi/dphi;

    lambda0 = tropt_max(0.001*uq,pow(uq*lq,2));

    while ((fabs(phi) < 0.01*(*delta)) || iter > 10) {

        phi =  calc_phi(evec, eval, Qtg, n, delta, lambda);
        dphi = calc_dphi(evec, eval, Qtg, n, lambda);
        *lambda = lambda0 - ((1./phi + *delta)/(*delta)) * (phi/dphi); 
        
        if ((*lambda >= lq) && (*lambda <= uq)) {
            *lambda = lambda0 - ((phi + *delta)/(*delta)) * (phi/dphi); 
        }  else {
            *lambda = fmax(0.001*uq, sqrt(lq*uq));
        }

        lq = fmax(lq, lambda0 - (phi/dphi));

        if (phi < 0) {
            uq = lambda0;
        } 

        lambda0 = *lambda;
        iter++;
    }
    return 0;
}


/*************************************************************/
/*                    calc_alpha                             */

static double calc_alpha(double *xnew, double *ulimit, 
                    double *llimit, int *isulimit, int 
                    *isllimit, int *ifree, int nfree, 
                    gsl_vector *xstep, 
                    gsl_vector *hit_boundary) {

    int i;
    double alpha, tempd0, tempd1, tempd2;
    gsl_vector *tempv = gsl_vector_alloc(nfree);

    for (i = 0;i<nfree;i++) {

        //Flag if we hit upper (1), lower (-1), or no (0) limit
        if (xnew[ifree[i]] >= ulimit[ifree[i]]) {
            gsl_vector_set(hit_boundary,i,1);
        } else if (xnew[ifree[i]] <= ulimit[ifree[i]]) {
            gsl_vector_set(hit_boundary,i,-1);
        } else {
            gsl_vector_set(hit_boundary,i,0);
        }

        //first check if step is non-zero
        if (gsl_vector_get(xstep,i)!=0) {

            //Calculate distances to nearest breakpoints

            if (isulimit[i]) {
                //Compute normalized step to upper boundary if it set
                tempd0 = (ulimit[i] - xnew[ifree[i]] ) / gsl_vector_get(xstep,i);
            } else {
                //No upper boundary
                tempd0 = 1.0;
            }

            if(isllimit[i]) {
                //Compute normalized step to lower boundary if set 
                tempd1 = (llimit[i] - xnew[ifree[i]]  )/ gsl_vector_get(xstep,i);
            } else {
                //No lower boundary
                tempd1 = 1.0;
            } 

        } else {
            tempd1 = 1.0;
        }

        //Find max of distance to breakpoint for free parameter
        tempd2 = tropt_max(tempd0,tempd1);

        gsl_vector_set(tempv,i,tempd2);
    }

    //Get distance to nearest breakpoint along xstep
    
    tempd0 = gsl_vector_min(tempv);
    alpha = tropt_min(1.,tempd0);
    
    //Free gsl vector
    gsl_vector_free(tempv);

    return alpha;
}

/**************************************************************/
/*                        check_is_bounded                    */
/* Does the step result in a solution that is within the     */
/*predefined boundaruy*/

static int check_is_bounded(gsl_vector *s, double *x, int *ifree,
                  int nfree, double *ulimit, double *llimit,
                  double *xnew, gsl_vector *hit_boundary) {
    
    int i, counter;
    double temp;

    counter = 0;
    //Check if all candidate steps result in  is within the  bounnds
    for (i=0;i<nfree;i++){

        temp = x[ifree[i]] + gsl_vector_get(s,i);

        if ((temp <= ulimit[ifree[i]] ) && (temp >= llimit[ifree[i]])) {
            counter++;
        }

        //Update vector indicating if (and where) the boundary is crossed'
        //0 - not crossed boundary
        //1 - crossed upper boundary
        //-1 - cross lower boundary
        if (temp >= ulimit[ifree[i]]) {
            gsl_vector_set(hit_boundary,i,1);
        } else if (temp <= llimit[ifree[i]]) {
            gsl_vector_set(hit_boundary,i,-1);
        } else {
            gsl_vector_set(hit_boundary,i,0);
        }
    }

    //update the xnew vector
    for (i=0;i<nfree;i++){
        xnew[ifree[i]] = x[ifree[i]] + gsl_vector_get(s,i);
    }

    //If all values within the boundary, then step is updated
    if (counter == nfree) {
        return 1; // is in bounds
    } else {
        return 0; // is not in bounds
    }

}

/*************************************************************/
/*              make_feasibl                                 */
/*                                                           */
/* Make proposed xnew x + step feasible based on upper and   */
/* lower bounds                                              */
/*                                                           */


// static void make_feasible(double *xall, double *xnew,
//          gsl_vector *xstep, 
//         int *isllimit, int *isulimit, 
//         double *llimit, double *ulimit,  
//         int *ifree, int nfree) {

//         int i;

        
//         for(i=0;i<nfree;i++) {

//             xnew[ifree[i]] = xall[ifree[i]] + gsl_vector_get(xstep,i);

//             if (isllimit[ifree[i]]) {
//                 if (xnew[ifree[i]] <= llimit[ifree[i]] ) {
//                     xnew[ifree[i]] =  + DBL_EPSILON;
                    
//                 }
//             }
//             if (isulimit[ifree[i]]) {
//                 if (xnew[ifree[i]] >= ulimit[ifree[i]] ) {
//                     xnew[ifree[i]] = ulimit[ifree[i]] - DBL_EPSILON;
//                 }
//             }

//             if ((xnew[ifree[i]] < llimit[ifree[i]]) || (xnew[ifree[i]] > ulimit[ifree[i]])) {
//                 xnew[ifree[i]] = 0.5*(llimit[ifree[i]]+ulimit[ifree[i]]);            
//             }

//         }
// }

/**************************************************************/
/*             calc_lnsph_inter                               */
/*                                                            */
/* Calculate the intersection of line between xnew and xold   */
/* with the boundary of the trust region's n-sphere.          */
/*                                                            */
/* Given the radius of the trust regionm, r, and the vectors: */
/* x0 = (x0_1, x0_2, ..., x0_n)                               */
/* x  = (x_1,   x_2, ...,  x_n)                               */
/*                                                            */
/* We can write the equuatio of the n-sphere as:              */
/*                                                            */
/*             2               2                   2    2     */
/* (x_1 - x0_1)  + (x_1 - x0_1) + ... + (x_n - x0_n) = r      */
/*                                                            */
/*                                                            */
/* We can respxress as a quadtratic                           */


double calc_lnsph_inter(gsl_vector *xhat, gsl_vector *shat, 
                double *delta) {

    double a, b, c, q, t1, t2;
    
    //     T         T        T         2
    //a = s s,  b = x s, c = x x - delta

    gsl_blas_ddot(shat,shat,&a);
    gsl_blas_ddot(xhat,shat,&b);
    gsl_blas_ddot(xhat,xhat,&c);
    c -= pow(*delta,2);

    q = -(b + get_sign(b)*sqrt(b*b - a*c));
    t1 = q/a;
    t2 = c/q;

    //return positive root
    if (t1 < t2) {
        return t2;
    } else {
        return t1;
    }
}

/*************************************************************/
/*                    calc_quadratic_1d                      */
/*                                                           */
/*                       T                   T               */
/* f(t) = 0.5 *(s_0 + s*t) Mhat (s0 + s*t) + g (s0 + s*t)    */
/*                                                           */
/*               T        2       T                 T         */
/* f(t) = 0.5*[ s  Mhat s t  +  2 s_0 Mhat s t  + s_0 Mhat S_0]*/
/*              T      T                                      */
/*           + g s0 + g st                                    */

static int calc_quadratic_1d(gsl_matrix *H, gsl_vector *g, 
                            gsl_vector *s, gsl_vector *s0,
                            gsl_matrix *DIAG, int n, int m,
                            double *a, double *b, double *c) {

    double tempd0, tempd2;

    gsl_vector *tempv0 = gsl_vector_calloc(n);
    gsl_vector *tempv1 = gsl_vector_calloc(n);
    
    gsl_blas_dgemv(CblasNoTrans,1.0, H, s, 0.0, tempv0);
    gsl_blas_ddot(s,tempv0,a);
    *a *=0.5;
    gsl_blas_ddot(g,s,b);
    *c = 0;

    if (s0) {
        gsl_blas_dgemv(CblasNoTrans,1.0, H, s0, 0.0, tempv1);
        gsl_blas_ddot(tempv1,s,&tempd2);
        *b += tempd2;

        gsl_blas_ddot(s0,tempv1,&tempd0);
        gsl_blas_ddot(g,s0,&tempd2);
        *c = (0.5*tempd0 + tempd2);
    }

    //cleanup
    gsl_vector_free(tempv0);
    gsl_vector_free(tempv1);

    return 0;
}


/*******************************************************************/
/*                    eval_quadratic_1d                           */

static int  eval_quadratic_1d(double a, double b, double c,
                             double lb, double ub, 
                             double *tout, double *yout) {

    int i, idx,idx_max;
    double extrm, tempd1, tempd2;
    double *t = (double*) calloc(3,sizeof(double));
    double *y = (double*) calloc(3,sizeof(double));

    t[0] = lb;
    t[1] = ub;
    t[2] = DBL_MAX;
    y[2] = DBL_MAX;

    idx_max = 2;

    if (a) {
        extrm = -0.5 * b / a;
        if ((extrm > lb) && (extrm < ub)) {
            t[2] = extrm;
            idx_max = 3;
        }
    }

    for (i=0;i<idx_max;i++) {
        y[i] = t[i] * (a * t[i] + b) + c;
    }

    //find min y value and corresponding idx
    array_min(y, 3, yout, &idx);

    //value of t corresponding to min y
    tempd1 = t[idx];
    *tout = tempd1;

    tempd2 = y[idx];
    *yout = tempd2;

    //cleanup
    free(t);
    free(y);
    return 0;
}

/*******************************************************************/
/*                    eval_quadratic                              */

// static double eval_quadratic(gsl_matrix *K, gsl_vector *gradf, 
//                         gsl_vector *s, gsl_matrix *C, int n ) {

//     double tempd0,tempd1;
//     gsl_vector *tempv = gsl_vector_calloc(n);
    
//     gsl_blas_dgemv(CblasNoTrans, 1.0, K, s, 0.0, tempv);
//     gsl_blas_ddot(tempv,tempv, &tempd0);

//     gsl_blas_ddot(s,gradf, &tempd1);

//     tempd0 *= 0.5;
//     tempd0 += tempd1;
    
//     gsl_vector_free(tempv);
    
//     return tempd0;
// }


/*******************************************************************/
/*                      solve_subproblem                           */
/*                                                                 */
/*Input:                                                           */
/* g: vector (n), Gradient of the cost functon                     */
/* B: matrix (n,n),  Hessian  of cost function                     */
/* W0: matrix (n,n), working matrix                                */
/* delta: double, trust region radius                              */
/* n: int, dimension of square Hessian matrix                      */
/*                                                                 */
/*                                                                 */
/*Output:                                                          */
/* s: vector (n),  solution step                                   */
/* outflag: integer                                                */
/*         -1: solution not found                                  */
/*          0: Unconstrained solution                              */
/*          1: Constrained solution (Easy)                         */
/*          2: Constrained solution (Hard 1)                       */
/*          3: Constrained solution (Hard 2)                       */
/*                                                                 */

static int solve_subproblem(gsl_vector *g, gsl_matrix *H, 
             gsl_vector *s, double *delta, int n) {
            
    int i, posdef, outflag, evalmin_idx;
    double temp, norm_s, delta_eps, evalmin;

    double *lambda;
    lambda = calloc(1,sizeof(double));

    delta_eps = (*delta) + sqrt(DBL_EPSILON);

    //Allocate memory for eigenvalue decomposition
    gsl_vector *eval = gsl_vector_calloc (n);
    gsl_vector *eval_sorted = gsl_vector_calloc (n);
    gsl_matrix *evec = gsl_matrix_calloc (n,n);
    gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (n);

    gsl_vector *Qtg  = gsl_vector_calloc(n);
    gsl_matrix *WORK  = gsl_matrix_calloc(n,n);
   
    //Copy matrix H in to W
    gsl_matrix_memcpy (WORK, H);

    //Compute eigenvalues. 
    gsl_eigen_symmv(WORK, eval, evec, w);

    //Sort eigenvalues into ascending order
    gsl_vector_memcpy (eval_sorted, eval);
    
    //Sorted eigenvalues
    //gsl_sort_vector(eval_sorted);
    evalmin = gsl_vector_min(eval);
    evalmin_idx = gsl_vector_min_index (eval);

    //         T
    //Compute Q g, where Q is the matrix of eigenvectors
    // and g is the gradient of the cost function
    gsl_blas_dgemv(CblasTrans, -1.0, evec, g, 0, Qtg);

    //Loop over the eigenvalues
    for (i=0;i<n;i++){
        //Check if any eigenvalues are negative
        if (gsl_vector_get(eval,i) < 0.0) {
            //Set positive definite flag false (0) and exit loop
            //as soon as a negative eigenvalue is encountered
            posdef = 0;
            break; 
        } else {
            posdef = 1;
        }
    }

    //Copy matrix H in to W again
    gsl_matrix_memcpy (WORK, H);

    //printf("is posdef = %d\n",posdef);

    if (posdef) {
        //If Hessian is positive definite, use Cholesky decomposition
        // and solve the Newton step 

        //Perform the cholesky decompositon of H, H replaced with 
        gsl_linalg_cholesky_decomp1(WORK);
        
        //Compute the inverse matrix of H from Cholesky decomposition
        gsl_linalg_cholesky_invert(WORK);

        //Solve  g = -inv(B) s to
        gsl_blas_dgemv(CblasNoTrans, -1.0, WORK, g, 0.0, s);

        //Compute euclidean norm
        norm_s = gsl_blas_dnrm2(s);

        //test if norm of direction is within the trust region
        if (norm_s <= delta_eps) {

            //CONSTRAINED SOLUTION FOUND.
            //Step lies inside trust region ||p|| <= delta
            //Return p as solution.
            //Output flag of 0 means the "Easy" case of the subproblem
            outflag = 0;
            goto cleanup_sub;

        } else {

            //UNCONSTRAINED POSDEF SOLUTION CASE (EASY)
            // Step lies outside trust region ||p|| > delta
            //Go on to iterative solution with initial
            //lambda set very small
            *lambda = DBL_EPSILON;
            outflag = 1;
        }

    } else {
        //UNCONSTRAINED SOLUTION CASE (HARD 1 or 2)
        //Matrix H is not positive definite
        //Move onto iterative solution and set 
        //the initial guess to the the smallest eigenvalue
        *lambda = -1*evalmin;
        outflag = 2;
    
    }

    //Calculate  = ||s|| - delta
    //Set problem as (A + lamda I ) s = g, solve for s.
    //If outflag > 0, then solution for p outstide trust region
    if (outflag > 0) {

        //Calculate L-M factor such that solution is found within
        //the trust region
        calc_lambda( evec, eval, Qtg, s, n,delta, lambda);

        //Calculate step
        calc_s(evec, eval, Qtg, n, lambda, s);

        //Calculate norm of p
        norm_s = gsl_blas_dnrm2(s);

        if (norm_s <= delta_eps) {
            goto cleanup_sub;
        }

    }

    //Solve HARD CASE 2
    //Eigenvector correspoding to eigmin is orthogonal to function gradient (qtg = 0)
    for (i=0;i<n;i++){
        temp = gsl_vector_get(eval,i) - *lambda;
        if (fabs(temp) < 0.0 + DBL_EPSILON) {
            gsl_vector_set(Qtg,i,0.0);
        }
    }

    calc_s(evec, eval, Qtg, n, lambda, s);

    temp = sqrt(fmax( pow(*delta,2) - pow(gsl_blas_dnrm2(s),2), 0.0)); 
    
    for (i=0;i<n;i++){
        temp = gsl_vector_get(s,i) + temp*gsl_matrix_get(evec,i,evalmin_idx);
        gsl_vector_set(s,i,temp);
    }

    //Output flag is 3 indicating the hard case was solved
    outflag = 3;
    goto cleanup_sub;

    //Free allocated memory
    cleanup_sub:
        if (lambda) free (lambda);
        if (w) gsl_eigen_symmv_free (w);
        if (eval) gsl_vector_free (eval);
        if (eval_sorted) gsl_vector_free (eval_sorted);
        if (Qtg) gsl_vector_free(Qtg);
        if (evec) gsl_matrix_free (evec);  
        if (WORK) gsl_matrix_free(WORK);
        
        return outflag;
}

/**************************************************************/
/*                      calc_cost_actual                      */
/*                                                            */
/* Computes actual least squares cost function                */
/*                                                            */
/*               T                                            */
/* fcost = 0.5 ey  W ey                                       */
/*                                                            */
/*Input:                                                      */
/* ey: vector (m), which is (y - F(x))                        */
/* temp1: vector (m), working temporary vector                */
/* Wy: matrix (m,m), weights matrix                           */
/* int: covyflag, is weight matrix populated                  */
/*                                                            */
/*Output:                                                     */
/* fcost: double, least squares cost function                 */

static int tropt_cost_actual(gsl_vector *deltay, gsl_vector *temp1, 
        gsl_matrix *Wy, int covyflag, double *fcost) {
        
    //compute actual cost function fcost
    if (covyflag) {
        gsl_blas_dgemv(CblasNoTrans, 1.0, Wy, deltay, 0.0, temp1);
        gsl_blas_ddot(deltay, temp1, fcost);
    } else if (!covyflag) {
        gsl_blas_ddot(deltay, deltay, fcost);
    }

    *fcost *= 0.5;    
    return 0; 
}

/**************************************************************/
/*                      calc_cost_model                       */
/*                                                            */
/* Computes cost function approximation model at x + xstep    */
/* via Taylor series expansion:                               */
/*                                                            */
/*           T            T                                   */
/* mcost =  g s +  0.5 * s H s                                */
/*                                                            */
/*                                                            */
/*                                                            */
/*Input:                                                      */
/* xnew: vector (n), xold + xstep                             */
/* temp1: vector(n), temporary working vector                 */
/* g: vector (n), Jacobian of the cost function at x          */
/* H: matrix (n,n), Hessian matrix of the cost function at x  */
/*                                                            */
/*Output:                                                     */
/* mcost: double, modeled cost function                       */

static int calc_cost_model(gsl_vector *xstep, 
        gsl_vector *g, gsl_matrix *H, gsl_vector *temp2, 
        double *mcost) {

    double term2; 
    double term3;

    //second term in Taylor series
    gsl_blas_ddot(g, xstep, &term2);

    //compute actual cost function fcost
    gsl_blas_dgemv(CblasNoTrans, 1.0, H, xstep, 0.0, temp2);
    gsl_blas_ddot(xstep, temp2, &term3);
    *mcost = term2 + 0.5* (term3);
    
    return 0;

}

/**************************************************************/
/*                    calc_tropt_gradf                        */
/*                                                            */
/* Calculate the gradient,gradf, and Hessian, H, of the least */
/* squares cost function Fx:                                  */
/*                                                            */
/*                T                                           */
/* Fx = [y - F(x)] Wy [y - F(x)]                              */
/*                                                            */
/* The gradient of the cost function:                         */
/*                                                            */
/*          T                                                 */
/* gradf = K Wy [y - F(x)]                                    */
/*                                                            */
/*Input:                                                      */
/* Wy: matrix (m x m), matrix of weights                      */
/* deltay: vector (m), deltay = y - F(x)                      */
/* temp1: vector (m), temporary working vector                */
/* K: matrix (m,n), Jacobian of the model function (dF/dx)    */
/*                                                            */
/*Output:                                                     */
/* gradf: vector (n), Jacobian of the cost function at x      */

static int calc_tropt_gradf(gsl_matrix *Wy, gsl_vector *deltay, 
        gsl_vector *temp1, gsl_matrix *K, gsl_vector *gradf) {
    
    //Weigthing matrix (inverse of covariance) provided

    //compute gradf
    gsl_blas_dgemv(CblasNoTrans, 1.0, Wy, deltay, 0.0, temp1);
    gsl_blas_dgemv(CblasTrans, 1.0, K, temp1, 0.0, gradf);

    return 0;
}

/**************************************************************/
/*                    calc_hessian                            */
/*                                                            */
/* Approximate Hessian, H, of the least square coost function */
/* Fx                                                         */
/*                                                            */
/*                T                                           */
/* Fx = [y - F(x)] Wy [y - F(x)]                              */
/*                                                            */
/* as:                                                        */
/*                                                            */
/*      T                                                     */
/* H = K Wy K                                                 */
/*                                                            */
/*Input:                                                      */
/* Wy: matrix (m x m), matrix of weights                      */
/* W2:, matrix (n,m), temporary working matrix                */
/* K: matrix (m,n), Jacobian of the model function (dF/dx)    */
/*                                                            */
/*Output:                                                     */
/* H: matrix (n,n), Hessian matrix of the cost function at x  */

static int calc_hessian(gsl_matrix *Wy, gsl_matrix *W2, 
            gsl_matrix *K, gsl_matrix *H) {
    
    //compute H
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Wy, K, 0.0, W2);
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, K, W2, 0.0, H);    

    return 0;
}

/**************************************************************/
/*                   tropt_matrix_invert                      */
/*                                                            */
/* Function uses LU decomposition to compute inverse of square*/
/* matrix                                                     */
/*                                                            */
/*Input:                                                      */
/* A: gsl matrix (mxm), square matrix                         */
/* diagflag: int, is matrix diagaonal: 0 (no) or 1 (yes).     */
/*                                                            */
/*Output:                                                     */
/* invA: gsl mmatrix (mxm), inverse of input square matrix    */
/*                                                            */

static int tropt_matrix_invert(gsl_matrix *A, gsl_matrix *invA, int diagflag){
    
    unsigned int n = A->size1;
    unsigned int m = A->size2;
    int s, i, j;
    double temp;
    gsl_permutation *p = gsl_permutation_alloc(m);
    
    if (m <= 0 || n <= 0) {
        printf("-E- %s line %d : input matrix dimensions less than equal to zero\n",
                __FILE__, __LINE__);
        exit(1);
    } else if (m != n) {
        printf("-E- %s line %d : input matrix dimensions not equal\n",
                __FILE__, __LINE__);
        exit(1);
    } else {
        //if not a diagonal matrix
        if (diagflag == 0) {
            // use a LU decomposition for the square matrix
            // Compute LU decomposition of A
            gsl_linalg_LU_decomp(A, p, &s);    
            gsl_linalg_LU_invert(A, p, invA);
        //if a diagonal matrix, set just the inverse of diagonal elements    
        } else {
            for (i = 0; i < A->size1; i++) {
                for (j = 0; j < A->size2; j++) {
                    if (j == i) {
                        temp = 1./gsl_matrix_get(A,i,j);  
                        gsl_matrix_set(invA,i,j,temp);
                    } else {
                        gsl_matrix_set(invA,i,j,0.0);
                    }
                }
            }
        }
    }
    //Cleanup
    gsl_permutation_free(p);
    return 0;
}

/**************************************************************/
/*                      calc_C                                */
/*                                                            */
/* Calculte the scaled matrix C:                              */
/*                                                            */
/* C = D diag(g) Jv D                                         */
/* Jv = diag(sign(g))                                         */
/*                                                            */
/* where,  D is a scaling matrix and g is the gradient of the */
/* cost function                                              */   
/*                                                            */
/*Input:                                                      */
/* g: vector (m,n),  gradient of the cost function (dF/dx)    */
/* B: matrix (n,n) Hessian matrix of cost function            */
/* invD: matrix (nxn),  inverse diagonal scaling matrix invD  */
/*                                                            */
/*Output:                                                     */
/* C: matrix (n,n)                                            */

// static int calc_C(gsl_vector *g, gsl_matrix *B, gsl_matrix *D,
//         gsl_matrix *invD, int n, gsl_matrix *C) {
    
//     int i, j, sign;
//     double temp;

//     gsl_matrix *diagG = gsl_matrix_calloc(n,n);

//     //Jv: Diagonal matrix whose ith element is the sign of the 
//     //ith element of the gradient of the cost function
//     gsl_matrix *Jv = gsl_matrix_calloc(n,n);

//     for (i=0;i<n;i++) {
//         for (j=0;j<n;j++) {
//             if (i == j) {
//                 temp = gsl_vector_get(g,i);
//                 sign = get_sign(temp);
//                 gsl_matrix_set(diagG,i,j,temp);
//                 gsl_matrix_set(Jv,i,j,sign);
//             }
//         }
//     }

//     //Calculate C
//     gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Jv, D, 0.0, C);
//     gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, diagG, C, 0.0, C);
//     gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, D, C, 0.0, C);

//     gsl_matrix_free(diagG);
//     gsl_matrix_free(Jv);

//     return 0;
// }

/**************************************************************/
/*                      calc_Mhat                             */
/*                                                            */
/* Calculte the scaled matrix Mhat                            */
/*         -1   -1                                            */
/* Mhat = D  M D                                              */
/*                                                            */
/* where,  D is a scaling matrix and                          */
/*                                                            */
/* M = H + C                                                  */
/*                                                            */
/* B is the Hessian matrix of the cost function and           */
/*                                                            */
/* C = D diag(g) Jv D                                         */
/* Jv = diag(sign(g))                                         */
/*                                                            */   
/*         -1  -1                                             */
/* Mhat = D  M  D + diag(g)diag(sign(g))                      */
/*                                                            */
/*Input:                                                      */
/* g: vector (m,n),  gradient of the cost function (dF/dx)    */
/* H: matrix (n,n) Hessian matrix of cost function            */
/* invD: matrix (nxn),  inverse diagonal scaling matrix invD  */
/*                                                            */
/*Output:                                                     */
/* Mhat: matrix (n,n), scaled matrix                          */
/* C: matrix (n,n)                                            */


static int calc_mhat(gsl_vector *gradf, gsl_matrix *invD, gsl_matrix *H, 
                int n, gsl_matrix *Mhat, gsl_matrix *Chat) {
    
    int i, j, sign;
    double temp;
    
    gsl_matrix *diagG = gsl_matrix_calloc(n,n);
    gsl_vector *gradfhat = gsl_vector_calloc(n);

    gsl_blas_dgemv(CblasNoTrans, 1.0, invD, gradf, 0.0, gradfhat);

    //Jv: Diagonal matrix whose ith element is the sign of the 
    //ith element of the gradient of the cost function
    gsl_matrix *Jv = gsl_matrix_calloc(n,n);

    for (i=0;i<n;i++) {
        for (j=0;j<n;j++) {
            if (i == j) {
                temp = gsl_vector_get(gradfhat,i);
                sign = get_sign(temp);
                gsl_matrix_set(diagG,i,j,temp);
                gsl_matrix_set(Jv,i,j,sign);
            } else {
                gsl_matrix_set(diagG,i,j,0.0);
                gsl_matrix_set(Jv,i,j,0.0);
            }
        }
    }

    //Calculate Chat
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, diagG, Jv, 0.0, Chat);

    //Calculate Mhat
    gsl_matrix_memcpy(Mhat,H);
    gsl_matrix_add (Mhat, Chat);

    gsl_matrix_free(diagG);
    gsl_matrix_free(Jv);
    gsl_vector_free(gradfhat)
;
    return 0;
}

/*******************************************************************/
/*                      set_affine_scaling_matrix                  */
/*                                                                 */
/* Generate the inverse of scaling matrix, invD. The matrix  D is  */
/* diagonal with elements determined from a vector, v. The values  */
/* of the matrix D are related to whether boundary constraints and */
/* their are set and,values.                                       */
/*                                                                 */
/* The sign of ith element the gradient of the cost function, g,   */
/* is taken into account as well as the ith  upper or lower limits */
/* if they exist.                                                  */
/*                                                                 */
/* if g < 0  and u < inf,  then v = x - u                          */
/* if g >= 0 and l > -inf, then v = x - l                          */
/* if g < 0  and u = inf,  then v = -1                             */
/* if g >= 0 and l = -inf, then v = 1                              */
/*                                                                 */
/*Input:                                                           */
/* isllimit: vector (n),  binary vector defining if lower limit    */
/* isulimit: vector (n),  binary vector defining if upper limit    */
/* llimit: vector (n),  lower limit if exists                      */
/* llimit: vector (n),  upper limit if exists                      */
/* xall: gsl vector(n), free parameter values                      */
/* g: gsl vector(n), gradient of the cost function at xall         */
/* n: int, number of free parameters                               */
/*                                                                 */
/*Output:                                                          */
/* d: gsl_vector(n,n), diagonal scaling matrix elements            */
/* dv: gsl_vector(n),  derivative of scaling matrix elements       */

static  int set_affine_scaling_matrix(int *isllimit, int *isulimit, 
        double *llimit, double *ulimit, double *xall, gsl_vector *g,
        int nfree, int *ifree, gsl_vector *v, gsl_vector *dv ) {

    int i;
    double vtemp = 0.0;

    for (i=0;i<nfree;i++) {
        
        if (isulimit[i] == 1) {
            if (gsl_vector_get(g,i) < 0) {
                vtemp = ulimit[ifree[i]] - xall[ifree[i]];
                gsl_vector_set(v,i,vtemp);
                gsl_vector_set(dv,i,-1.0);
            } 
        } else {
            gsl_vector_set(v,i,1.0);
            gsl_vector_set(dv,i,0.0);
        }

        if (isllimit[i] == 1) {
            if (gsl_vector_get(g,i) > 0)  {
                vtemp = xall[ifree[i]] - llimit[ifree[i]];
                gsl_vector_set(v,i,vtemp);
                gsl_vector_set(dv,i,1.0);
            }
        } else {
            gsl_vector_set(v,i,1.0);
            gsl_vector_set(dv,i,0.0);
        }
    }

    return 0;

}

/*******************************************************************/
/*                      set_scaling_matrix                         */
/*                                                                 */
/* Generate the inverse of scaling matrix, invD. The matrix  D is  */
/* diagonal with elements determined from a vector, v. The values  */
/* of the matrix D are related to whether boundary constraints and */
/* their are set and,values.                                       */
/*                                                                 */
static int set_scaling_matrix(gsl_vector *g, gsl_matrix *D,
            gsl_matrix *invD, int n, int iter) {

    int i,j;
    double tempd, tempd_new, tempd_old;

    if (iter == 0) {
        //If iter=0, set intial scaling matrix to identity
        //set diagonal elements only, others are zero.
        for (i=0;i<n;i++) {
            for (j=0;j<n;j++) {    
                if (j == i) {
                    gsl_matrix_set(D,i,j, 1.0);
                    gsl_matrix_set(invD,i,j, 1.0);
                } else {
                    gsl_matrix_set(D,i,j, 0.0);
                    gsl_matrix_set(invD,i,j, 0.0);
                }
            }
        } 
    } else {
        //if iter > 0, then update scaling
        for (i=0;i<n;i++) {
            
            tempd_old = gsl_matrix_get(D,i,i);
            tempd_new = fabs(gsl_vector_get(g,i));
            tempd = tropt_max(tempd_old,tempd_new);
            tempd = sqrt(tempd);

            for (j=0;j<n;j++) {    
                if (j == i) {
                    gsl_matrix_set(D,i,j, 1./tempd);
                    gsl_matrix_set(invD,i,j,tempd);
                } else {
                    gsl_matrix_set(D,i,j, 0.0);
                    gsl_matrix_set(invD,i,j, 0.0);
                }
            }
        } 
    }
    return 0;
}

/**************************************************************/
/*                      calc_ghat                             */
/*                                                            */
/* Scale the gradient of the cost function, g, using a scaling*/
/* matrix D.                                                  */
/*          -1                                                */
/* g_hat = D  g                                               */
/*                                                            */
/*Input:                                                      */
/* g: vector (n),  gradient of the cost function              */
/* invD: matrix (nxn),  inverse of scaling matrix D           */
/*                                                            */
/*Output:                                                     */
/* ghat: matrix (m,n),Jacobian matrix of function (dF/dx) */

// static int calc_ghat(gsl_vector *g, gsl_matrix *invD,
//         gsl_vector *ghat) {

//     gsl_blas_dgemv(CblasNoTrans, 1.0, invD, g, 0.0, ghat);

//     return 0;
// }

/**************************************************************/
/*                      posdef_check                          */
/*                                                            */
/* Check if Hessian matrix, H, is positive definite by        */
/* examing the eigenvalues to ensure they are all positive.   */
/*                                                            */
/*                                                            */
/*Input:                                                      */
/* H: matrix (n,n),  Hessian  matrix                          */
/* n: int, dimension of square Hessian matrix                 */
/*                                                            */
/*Output:                                                     */
/* flag: int, 1 is posdef, gt 0 is non-posdef                 */

// static int posdef_check(gsl_matrix *A, int n) {

//     int i, flag;

//     //Allocate memory for eigenvalue decomposition
//     gsl_vector *s = gsl_vector_calloc (n);
//     gsl_matrix *V = gsl_matrix_calloc (n, n);
//     gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (n);

//     gsl_eigen_symmv(A, s, V, w);

//     flag = 1;
//     for (i=0;i<n;i++){
//         //If any eigenvalues are negative, flag as not posdef
//         if (gsl_vector_get(s,i) < 0) {
//             flag = 0;
//         }
//     }

//     //Free allocated memory
//     gsl_eigen_symmv_free (w);
//     gsl_vector_free (s);
//     gsl_matrix_free (V);

//     return flag;
// }



/*******************************************************************/
/*                      update_trust_region                        */
/*                                                                 */
/* Update the trust region radius, delta per Noecedal & Wright     */
/*                                                                 */                                                  
/*                                                                 */  
/*                                                                 */
/*Input:                                                           */
/* mu, double,  threshold 1 (default 0.25                          */
/* eta, double, threshold 2 (default 0.75                          */
/* theta, double, threshold 3 (default 0.95)                       */
/* gamma1, double, update scaling     (default 0.25)               */
/* gamma1, double, update scaling    (default 2.0)                 */
/* actual_reduction, double, actual reduction in cost function     */
/* prec_reduction, double, predicted reduction in cost function    */
/* xstep: vector (n), candidate step                               */
/* delta: double, trust region radius                              */
/*                                                                 */
/*Output:                                                          */
/* rho: double, ratio of actual to predicted reduction             */
/* delta: double, trust region radius  (update in place)           */


static int update_trust_region(double mu, double eta, 
            double theta, double gamma1, double gamma2,
            double *actual_reduction, double *pred_reduction,
            double norm_xstep, double *rho, double *delta) {

    //*rho = *actual_reduction / -(*pred_reduction);
    
    //Check that pred_reduction less than zero;
    if (*pred_reduction < 0) {
        *rho = *actual_reduction / -(*pred_reduction);
    } else if (*pred_reduction == *actual_reduction) {
        *rho = 1;
    } else {
        *rho = 0;
    }

    

    //Following Noceal and Wright
    if (*rho < 0.25)  {
        *delta = 0.25*norm_xstep;
    
    }  else {

        if ((*rho > 0.75)  && (norm_xstep >= theta * (*delta))) {
            *delta = tropt_min(2.*(*delta), 100.);
        }
    }

    return 0;

}


/***********************************************************************/
/*                    trustopt_jac                                     */
/*                                                                     */
/* Function that computes and estimate of the function jacobian matrix */
/* using finite difference approximation                               */

static int tropt_jac(trustopt_func funct,
        int m, int nfree, int *ifree, int npar, double *x, 
        gsl_matrix *K, double epsfcn,
        void *priv, int *nfev,
        double *step, double *dstep, int *dside,
        double **dydpar) {

    int i,j;
    int iflag = 0;

    int dsidei;
    
    double eps,h,temp;
    double K_temp;
    
    double *yhat0;
    double *yhat1;
    double *yhat2; 

    yhat0 = (double *) calloc(m, sizeof (double));
    yhat1 = (double *) calloc(m, sizeof (double));
    yhat2 = (double *) calloc(m, sizeof (double));

    temp = fmax(epsfcn,DBL_EPSILON);
    eps = sqrt(temp);

    //Reset Jacobian matrix to zero
    gsl_matrix_scale(K,0.0); 

    //Check which parameters have analytical derivatives and which
    //need numerical ones 
    for (j=0; j<nfree; j++) {  /* Loop through free parameters only */

        if (dside[ifree[j]] == 3) {
            /* If there are any parameters requiring analytical derivatives,
            then compute them first. */
            for (i=0; i<m; i++) {
                K_temp = dydpar[i][ifree[j]];    
                gsl_matrix_set(K,i,j,K_temp);  
            }
            
        }
        
        if (dside[ifree[j]] != 3) {
            /*Compute numerical derivatives*/ 
            if (dside){
                dsidei = dside[ifree[j]];
            } else {
                dsidei = 0;
            }
        
            iflag = tropt_call(funct, m, npar, x, yhat0, dydpar, priv);
            
            if (iflag < 0){
                goto cleanup_jac;
            }
                  
            //Set default step size for numerical derivative
            temp = x[ifree[j]];
            h = eps * fabs(temp); 
            
            //Check for absolute step size provided by user
            if (step  &&  step[ifree[j]] > 0){
                h = step[ifree[j]]; 
            }
            //Check for relative step size provided by uers
            if (dstep && dstep[ifree[j]] > 0) {
                h = fabs(dstep[ifree[j]]*temp);
            }
            if (h == 0.0) {
                h = eps;
            }

            /* If negative step requested*/
            if ((dside && dsidei == -1) || (dside && dsidei == 0 )) {
                h = -h;
            }

            x[ifree[j]] = temp + h;

            iflag = tropt_call(funct, m, npar, x, yhat1, dydpar, priv);

            if (nfev) {
                *nfev = *nfev + 1;
            } else {
                *nfev = 1;
            }
            if (iflag < 0 ) {
                goto cleanup_jac;
            }

            x[ifree[j]] = temp;

            if (dsidei <= 1) {
            /* One sided derivative*/
                for (i=0; i<m; i++) {
                    K_temp = (yhat1[i] - yhat0[i])/h; /* fjac[i+m*j] */
                    //Check of free or fixed parameter
                    gsl_matrix_set(K,i,ifree[j],K_temp);
                }
            } else {  /* dside > 2 */
                /* Two sided derivative*/

                x[ifree[j]] = temp - h;
                iflag = tropt_call(funct, m, npar, x, yhat2, dydpar, priv);

                if (nfev) {
                    *nfev = *nfev + 1;
                } else {
                    *nfev = 1;
                }
                
                if (iflag < 0 ) {
                    goto cleanup_jac;
                }
                
                x[ifree[j]] = temp;

                for (i=0; i<m; i++) {
                    K_temp = (yhat1[i] - yhat2[i])/(2*h);
                    //Check of free or fixed parameter 
                    gsl_matrix_set(K,i,j,K_temp);
                }

            } 
        } 
    }
    
  cleanup_jac:
    if (yhat0) free(yhat0);
    if (yhat1) free(yhat1);
    if (yhat2) free(yhat2);

    if (iflag < 0){
        return iflag;
    } else {
        return 0;
    }
    
}


/***********************************************************************/
/*                    Utility functions */

//Get the sign of a double (+ve = 1, zero = 0, -ve, = -1)
static int get_sign(double a) {
    int sign;
    if (a < 0) {
        sign = -1;
    } else if (a > 0) {
        sign = 1;
    } else if (a == 0) {
        sign = 0;
    }
    return sign;
}

//Find the maximum of two doubles
static double tropt_max(double a, double b) {
    if (a >= b) {
        return(a);
    } else {
        return(b);
    }
}

//Find the minimum of two doubles
static double tropt_min(double a, double b) {
    if (a <= b) {
      return(a);
    } else {
      return(b);
    }
}


//Find the 2-norm of 1-D C array
static double tropt_norm2_xall(double *x, int *ifree, int n) {

    int i;
    double temp = 0;

    for (i=0;i<n;i++) {
        temp += pow(x[ifree[i]],2.);
    }

    return sqrt(temp);
}

//Find the maximum absolute value in gsl_vector v
// static double vector_abs_max(gsl_vector *v) {

//     int i;
//     double max = 0;
//     double temp;

//     for (i=0;i<v->size;i++) {
//         temp = gsl_vector_get(v,i);
//         if (fabs(temp) > max) {
//             max = fabs(temp);
//         }
//     }
//     return max;

// }

//Find the minumum value in 1D c array and index of min
//array_min(y, 3, yout, &idx);
static int array_min(double *v, int n, double *min, int *idx) {
    int i;
    *min = DBL_MAX;
    double temp;

    for (i=0;i<n;i++) {
        temp = v[i];
        if (temp < *min) {
            *min = temp;
            *idx = i;
        }
    }
    return 0;
}

//Free a 2D C array
static int tropt_free2D(double **Arr, int m){
    int i;
    for (i = 0; i < m; ++i) {
        free(Arr[i]);
    }
    free(Arr);
    return 0;
}
