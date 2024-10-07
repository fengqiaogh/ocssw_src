#ifndef _GET_SDP_H
#define _GET_SDP_H

static const float pigmin = 0.001;
static const float pigmax = 1000.0;

#define SDPBADRRS   0x0001;
#define SDPSAAFAIL  0x0002;
#define BADTCHL     0x0004;
#define BADZEA      0x0008;
#define BADDVCHLA   0x0010;
#define BADBUTOFUCO 0x0020;
#define BADHEXOFUCO 0x0040;
#define BADALLO     0x0080;
#define BADMVCHLB   0x0100
#define BADNEO      0x0200;
#define BADVIOLA    0x0400;
#define BADFUCO     0x0800;
#define BADCHLC12   0x1000;
#define BASCHLC3    0x2000;
#define BADPERID    0x4000;

typedef struct sdp_data_str {
    
    int nvbands;
    int nfitbands;
    int nfbands;
    int nfree;
    int i440;
    int i440f;
    int i490;
    int i555;
    
    int *bindx;
    int *fbindx;
    
    float *bands;
    float *fbands;
    
    float adg_s;
    float bbp_s;
    
    float *g;
    float *aw;
    float *bbw;
    float *adgstar;
    float *bbpstar;
    float *aphA;
    float *aphB;
  
    float *faw;
    float *fbbw;
    float *fadgstar;
    float *fbbpstar;
    float *faphA;
    float *faphB;
    double **Acoeff;
    double **Ccoeff;
    
    double *popt;
    float *Rrs_a;

} sdpstr;

#endif