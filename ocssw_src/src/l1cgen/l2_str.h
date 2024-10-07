//
//  l2_str.h
//
//
//  Created by Martin Montes on 10/26/20.
//

#ifndef L2_STR_H
#define L2_STR_H

#include "l1c_filehandle.h"
#include <stdbool.h>
#include <stdint.h>
#include "l1c_input.h"

namespace l1c {

class l2_str {
   protected:
   public:
    // methods
    l2_str();
    virtual ~l2_str();

    // Open, read, close
    virtual int32_t openl2_ocis_l1c(L1C_input *l1cinput, l2_str *l2str, l1c_filehandle *l1cfile,
                                    int16_t *file_id);
    virtual int32_t readl2_ocis_l1c(l2_str *l2str, l1c_filehandle *l1cfile, int16_t *file_id, int32_t recnum);
    virtual int32_t closel2_ocis_l1c(l2_str *l2str, l1c_filehandle *l1cfile);

    // global attributes
    size_t npix;
    size_t iscan;
    size_t nscan;
    // scan line attributes--
    double *ev_mid_time;
    size_t *scan_quality_flag;
    size_t spix;
    size_t epix;
    size_t dpix;

    // structure--pointers to data arrays
    // structure--pointers to data arrays
    float **Lt;       // dim depends between sensors--OCI is bands x pixels --
    float **Lt_blue;  //[num_views][num_pol][num_blue_bands][num_pixels]
    float **Lt_red;
    float **Lt_SWIR;

    float **l2prod;  // # l2 products x number of pixels
    size_t nl2prod;
    float *slopeprod;
    float *offsetprod;
    float *tilt;  // sensor tilt

    // sensor/sun geometry
    float *senz;
    float *sena;
    float *solz;
    float *sola;
    float *delphi;
    float *scattang;

    // navigation--
    float att_ang[3];
    float orb_pos[3];
    float orb_vel[3];

    // geolocation--
    float *senazpix;
    float *latpix;
    float *lonpix;
    float *latpix2;  // lat2/lon are the sline +1
    float *lonpix2;

    l1c_filehandle *l1cfile;  // not sure why this thing is here?
};

// prototypes------
//  Open, read, close
int32_t openl2_ocis_l1c(L1C_input *l1cinput, l2_str *l2str, l1c_filehandle *l1cfile, int16_t *file_id);
int32_t readl2_ocis_l1c(l2_str *l2str, l1c_filehandle *l1cfile, int16_t *file_id, int32_t recnum);
int32_t closel2_ocis_l1c(l2_str *l2str, l1c_filehandle *l1cfile);

}  // namespace l1c

#endif /* L2_STR_H */
