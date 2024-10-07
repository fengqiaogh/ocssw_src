#ifndef _FILEHDR_STRUC_H
#define _FILEHDR_STRUC_H

typedef struct filehdr_struct {
    int32_t sensorID;
    int32_t length;
    int32_t npix;
    int32_t format;
    int32_t nscan;
} filehdr;

#endif



