#ifndef TRUSTOPT_H
#define TRUSTOPT_H

/* This is a C library.  Allow compilation with a C++ compiler */
#ifdef __cplusplus
extern "C" {
#endif


struct trustopt_par_struct {
  int fixed;        /* 1 = fixed; 0 = free */
  int limited[2];   /* 1 = low/upper limit; 0 = no limit */
  double limits[2]; /* lower/upper limit boundary value */

  double step;      /* Step size for finite difference */
  double relstep;   /* Relative step size for finite difference */
  int side;         /* Sidedness of finite difference derivative 
		        0 - one-sided derivative computed automatically
		        1 - one-sided derivative (f(x+h) - f(x)  )/h
		       -1 - one-sided derivative (f(x)   - f(x-h))/h
		        2 - two-sided derivative (f(x+h) - f(x-h))/(2*h) 
			      3 - user-computed analytical derivatives*/
};

// Definition of configuration structure for trustopt
struct trustopt_config_struct {

    double maxiter;   //maximum iterations
    double eps1;      //gradient tolerance (gtol)
    double eps2;      //parameters tolerance  (xtol)
    double eps3;      //chi-squared tolerance (ftol)
    double step_size; //finite differnce step size;

    double delta0;     //Trust region radius initial
    double delta_min;  //Trust region radius min
    double delta_max;  //Trust region radium max
    double gamma1;     //Trust region scale 1
    double gamma2;     //Trust region scale 2
    double mu;         //Trust region acceptance tolerance
    double eta;        //Trust region increase tolerance
    double theta;      //stepback parameter
    int prnt;          //print outputs

} ;


// Definition of output structure for trustopt
struct trustopt_result_struct {
  double redX2;      /*reduced Chi-squared error term*/
  int niter;         /*Number of iterations*/
  int nfev;          /*Number of function evaluations */
  double cost;       /*cost function value*/
  double costm;      /*measurement contribution to cost function*/
  double costa;      /*aprior contribution to cost function*/
  int status;        /* Model status code (see below for defintions) */

}; 


//Struct containing the configuration paramters
typedef struct trustopt_par_struct trustopt_par;
typedef struct trustopt_config_struct trustopt_config;
typedef struct trustopt_result_struct trustopt_result;

/* Enforce type of fitting function */
typedef int (*trustopt_func)(int m, /* Number of functions (elts of fvec) */
		       int n, /* Number of variables (elts of x) */
		       double *x,      /* I - Parameters */
		       double *fvec,   /* O - function values */
		       double **dydpar,  /* O - function derivatives (optional)*/
		       void *private_data); /* I/O - function private data*/


extern int trustopt(trustopt_func funct, int m, int npar, int covydiag, 
        double *xall, double *yall, double **covy, 
        trustopt_par *pars, trustopt_config *config, 
        void *private_data, trustopt_result *result); 

//Define error codes here
#define TROPT_ERR_NOFUNC     -1  //No input function was supplied
#define TROPT_ERR_FUNCCALL   -2  //Error calling the model function
#define TROPT_ERR_NPAR       -3  //Number of paramters <= zero
#define TROPT_ERR_NFREE      -4  //Number of free paramters <= zero
#define TROPT_ERR_PARAM      -5  //Negative configuration paramters provided
#define TROPT_ERR_DOF        -6  // Degrees of freedom inadequate
#define TROPT_COVYFLAG_ERR   -7  //Y covariance diagaonal flag set incorrectly. 
                              //Must be 0 or 1.
#define TROPT_COVXFLAG_ERR   -8  //X covariance diagaonal flag set incorrectly. 
                               //Must be 0 or 1.
#define TROPT_COSTNAN_ERR    -9  //Cost function is NaN
#define TROPT_NO_LLIMSET     -10 //no lower limit set for free parameter
#define TROPT_NO_ULIMSET     -11 //no upper limit set for free parameter


/* Potential success status codes */
#define TROPT_SUCCESS_FTOL    1   // Convergence based on cost function tolerance
#define TROPT_SUCCESS_XTOL    2   // Convergence based on cost function tolerance
#define TROPT_SUCCESS_JGRAD   3   // Convergence in gradient of cost function
#define TROPT_SUCCESS_CLOSE   4   //Initial guess very close to solution
#define TROPT_SUCCESS_STEP    5   // Change in parameter step is below threshold
#define TROPT_MAXITER         6   // Maximum number of iterations reached


#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* TRUSTOPT_H */



