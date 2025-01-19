#ifndef _GET_GIP_H
#define _GET_GPIG_H

static const float pigmin = 0.001;
static const float pigmax = 1000.0;

#define GPIGBADRRS   0x0001;
#define GPIGSAAFAIL  0x0002;
#define BADTCHL      0x0004;
#define BADCHLC12    0x0008;
#define BADCHLB      0x0010;
#define BADPPC       0x0020;


typedef struct gpig_data_str {
    
    int nvbands;
    int nfitbands;
    int nfbands;
    int nfree;
    int i400;
    int i400f;
       
    int *bindx;
    int *gindx;
    
    float *fitbands;
    
    float *g;
    float *aw;
    float *psiT;
    float *psiS;
    float *bbw;
    float *Rrs_unc;
    float **aphGauss;
  
    double *popt;
    float *Rrs_a;
    float *uterm;

} gpigstr;

#endif