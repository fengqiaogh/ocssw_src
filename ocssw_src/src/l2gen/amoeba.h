typedef struct {
    short int nfunc;
    short int npnts;
    double merit;
    double *y;
    double *wgt;
    double *fitfunc;
    double *yfit;
    int niter;
    void *meta;
} FITSTRUCT;

short amoeba(double *, FITSTRUCT *,
        double (*) (FITSTRUCT *, double []), float);
