#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <dfutils.h>
#include <genutils.h>
#include <productInfo.h>

#include <hdf.h>
#include <mfhdf.h>

/**
 * getDataTypeInt - returns the HDF definitions for the productInfo dataype
 * @param p_info
 * @return 
 */
int16_t getDataTypeInt(productInfo_t *p_info) {
    int16_t nt;

    if (strcmp(p_info->dataType, "byte") == 0)
        nt = DFNT_INT8;
    else if (strcmp(p_info->dataType, "ubyte") == 0)
        nt = DFNT_UINT8;
    else if (strcmp(p_info->dataType, "short") == 0)
        nt = DFNT_INT16;
    else if (strcmp(p_info->dataType, "ushort") == 0)
        nt = DFNT_UINT16;
    else if (strcmp(p_info->dataType, "int") == 0)
        nt = DFNT_INT32;
    else if (strcmp(p_info->dataType, "uint") == 0)
        nt = DFNT_UINT32;
    else if (strcmp(p_info->dataType, "float") == 0)
        nt = DFNT_FLOAT32;
    else if (strcmp(p_info->dataType, "double") == 0)
        nt = DFNT_FLOAT64;
    else {
        printf("-E- %s %d: datatype %s is not valid\n", __FILE__, __LINE__, p_info->dataType);
        exit(1);
    }

    return nt;
}
/* -------------------------------------------------------- */

/* -------------------------------------------------------- */
int16_t *float2int16(float32 fbuf[], int32_t spix, int32_t npix, int incr,
        float slope, float offset) {
    static int32_t npix_alloc = 0;
    static int16 *ibuf = NULL;
    float32 fval;
    int32_t i;

    double maxval = slope * 32767 + offset;
    double minval = slope * (-32766) + offset;

    /*
     * -32766 is used for minval to allow -32767 to remain the sentinal bad value
     * Yes, this does leave one digit hanging out down low, but meh....
     */

    if (npix > npix_alloc) {
        npix_alloc = npix;
        if (ibuf)
            free(ibuf);
        if ((ibuf = calloc(npix, sizeof (int16))) == NULL) {
            fprintf(stderr,
                    "-E- %s line %d: Unable to allocate buffer space.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
    }

    for (i = 0; i < npix; i++) {
        fval = fbuf[spix + i * incr];
        if (fval == BAD_FLT)
            ibuf[i] = BAD_INT;
        else if (fval >= maxval)
            ibuf[i] = 32767;
        else if (fval <= minval)
            ibuf[i] = -32766;
        else
            ibuf[i] = round((fval - offset) / slope);
    }

    return (ibuf);
}

/* -------------------------------------------------------- */

/* -------------------------------------------------------- */
uint16_t *float2uint16(float32 fbuf[], int32_t spix, int32_t npix, int incr,
        float slope, float offset) {
    static int32_t npix_alloc = 0;
    static uint16 *ibuf = NULL;
    float32 fval;
    int32_t i;

    double maxval = slope * 65534 + offset;
    double minval = offset;

    if (npix > npix_alloc) {
        npix_alloc = npix;
        if (ibuf)
            free(ibuf);
        if ((ibuf = calloc(npix, sizeof (uint16))) == NULL) {
            fprintf(stderr,
                    "-E- %s line %d: Unable to allocate buffer space.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
    }

    for (i = 0; i < npix; i++) {
        fval = fbuf[spix + i * incr];
        if (fval == BAD_FLT)
            ibuf[i] = BAD_UINT;
        else if (fval >= maxval)
            ibuf[i] = 65534;
        else if (fval <= minval)
            ibuf[i] = 0;
        else
            ibuf[i] = round((fval - offset) / slope);
    }

    return (ibuf);
}
/* -------------------------------------------------------- */

/* -------------------------------------------------------- */
uint8_t *float2uint8(float32 fbuf[], int32_t spix, int32_t npix, int incr,
        float slope, float offset) {
    static int32_t npix_alloc = 0;
    static uint8_t *ibuf = NULL;
    float32 fval;
    int32_t i;

    double maxval = slope * 254 + offset;
    double minval = slope * 0 + offset;

    if (npix > npix_alloc) {
        npix_alloc = npix;
        if (ibuf)
            free(ibuf);
        if ((ibuf = calloc(npix, sizeof (int8))) == NULL) {
            fprintf(stderr,
                    "-E- %s line %d: Unable to allocate buffer space.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
    }

    for (i = 0; i < npix; i++) {
        fval = fbuf[spix + i * incr];
        if (fval == BAD_FLT)
            ibuf[i] = BAD_UBYTE;
        if (fval >= maxval)
            ibuf[i] = 254;
        else if (fval <= minval)
            ibuf[i] = 0;
        else
            ibuf[i] = round((fval - offset) / slope);
    }

    return (ibuf);
}

/* -------------------------------------------------------- */

/* -------------------------------------------------------- */
int8_t *float2int8(float32 fbuf[], int32_t spix, int32_t npix, int incr,
        float slope, float offset) {
    static int32_t npix_alloc = 0;
    static int8_t *ibuf = NULL;
    float32 fval;
    int32_t i;

    double maxval = slope * 127 + offset;
    double minval = slope * -127 + offset;

    if (npix > npix_alloc) {
        npix_alloc = npix;
        if (ibuf)
            free(ibuf);
        if ((ibuf = calloc(npix, sizeof (int8))) == NULL) {
            fprintf(stderr,
                    "-E- %s line %d: Unable to allocate buffer space.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
    }

    for (i = 0; i < npix; i++) {
        fval = fbuf[spix + i * incr];
        if (fval == BAD_FLT)
            ibuf[i] = BAD_BYTE;
        if (fval >= maxval)
            ibuf[i] = 127;
        else if (fval <= minval)
            ibuf[i] = -127;
        else
            ibuf[i] = round((fval - offset) / slope);
    }

    return (ibuf);
}


/* -------------------------------------------------------- */

/* -------------------------------------------------------- */
void *scale_sds(float *data, productInfo_t *p, int npix) {
    void *pbuf;

    switch (getDataTypeInt(p)) {
    case DFNT_INT16:
        pbuf = (void *) float2int16(data, 0, npix, 1, p->scaleFactor, p->addOffset);
        break;
    case DFNT_UINT16:
        pbuf = (void *) float2uint16(data, 0, npix, 1, p->scaleFactor, p->addOffset);
        break;
    case DFNT_INT8:
        pbuf = (void *) float2int8(data, 0, npix, 1, p->scaleFactor, p->addOffset);
        break;
    case DFNT_UINT8:
        pbuf = (void *) float2uint8(data, 0, npix, 1, p->scaleFactor, p->addOffset);
        break;
    case DFNT_FLOAT32:
        fprintf(stderr,
                "-W- %s Line %d: Stubbornly refusing to scale floating point data: \n\t%s\n", __FILE__, __LINE__, p->productName);
        pbuf = data;
        break;
    default:
        fprintf(stderr,
                "-W- %s Line %d: Unknown data type %s.\n", __FILE__, __LINE__, p->dataType);
        pbuf = data;
        break;
    }

    return (pbuf);
}


/* -------------------------------------------------------- */

/* -------------------------------------------------------- */
float *unscale_sds(void *data, productInfo_t *p, int32_t spix, int32_t npix, int incr) {
    static int32_t npix_alloc = 0;
    static float32 *fbuf = NULL;

    float32 fval, *fptr;
    int16 ival, *iptr;
    uint8_t bval, *bptr;
    int32_t i;

    if (npix > npix_alloc) {
        npix_alloc = npix;
        if (fbuf)
            free(fbuf);
        if ((fbuf = calloc(npix, sizeof (float32))) == NULL) {
            fprintf(stderr,
                    "-E- %s line %d: Unable to allocate buffer space.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
    }

    switch (getDataTypeInt(p)) {
    case DFNT_INT8:
    case DFNT_UINT8:
        bptr = (uint8_t *) data;
        for (i = 0; i < npix && i < npix; i++) {
            bval = bptr[spix + i * incr];
            fbuf[i] = bval * p->scaleFactor + p->addOffset;
        }
        break;
    case DFNT_INT16:
    case DFNT_UINT16:
        iptr = (int16 *) data;
        for (i = 0; i < npix && i < npix; i++) {
            ival = iptr[spix + i * incr];
            fbuf[i] = ival * p->scaleFactor + p->addOffset;
        }
        break;
    case DFNT_FLOAT32:
        fptr = (float32 *) data;
        for (i = 0; i < npix && i < npix; i++) {
            fval = fptr[spix + i * incr];
            fbuf[i] = fval * p->scaleFactor + p->addOffset;
        }
        break;
    default:
        fprintf(stderr, "-W- %s Line %d: Unknown data type %s product %d.\n",
                __FILE__, __LINE__, p->dataType, p->cat_ix);
        break;
    }

    return (fbuf);
}


