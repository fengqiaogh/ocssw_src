#ifndef _TARGET_STRUC_H
#define _TARGET_STRUC_H

typedef struct target_struct {
    int32_t sensorID;
    int32_t subsensorID;
    int32_t length;
    int32_t npix;
    int32_t mode;
    char *data;
    int32_t *year;
    int32_t *day;
    int32_t *msec;
    float *solz;
    float *Lw;
    float *nLw;
} tgstr;

#endif



