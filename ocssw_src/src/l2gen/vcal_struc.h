#ifndef _VCAL_STRUC_H
#define _VCAL_STRUC_H

typedef struct vcal_struct {
    int32_t sensorID;
    int32_t length;
    int32_t npix;
    int32_t nbands;
    char *data;
    float *vLt;
    float *Lw;
    float *nLw;
    float *tLw;
    float *brdfsat;
    float *brdftgt;
} vcstr;

#endif



