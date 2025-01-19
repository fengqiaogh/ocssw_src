//  l1c_str.h
//  Created by Martin Montes on 8/15/2022
//

#ifndef L1C_STR_H
#define L1C_STR_H

#include "l1c_filehandle.h"
#include <stdbool.h>
#include <stdint.h>
#include "l1c_input.h"
#include <filehandle.h>
#include <netcdf>

#include "l1c_latlongrid.h"

namespace l1c {

class l1c_str {
   protected:
   public:
    // methods
    l1c_str();
    virtual ~l1c_str();

    // Open, read, close
    virtual int32_t openl1b_oci_l1c(l1c_str *l1cstr, l1c_filehandle *l1cfile);
    virtual int32_t closel1b_oci_l1c(l1c_str *l1cstr, l1c_filehandle *l1cfile);

    virtual int32_t openl1b_spex_l1c(l1c_str *l1cstr, l1c_filehandle *l1cfile);
    virtual int32_t readl1b_spex_l1c(l1c_str *l1cstr, l1c_filehandle *l1cfile, int32_t recnum);
    virtual int32_t closel1b_spex_l1c(l1c_str *l1cstr, l1c_filehandle *l1cfile);

    virtual int32_t openl1b_harp2_l1c(l1c_str *l1cstr, l1c_filehandle *l1cfile);
    virtual int32_t readl1b_harp2_l1c(l1c_str *l1cstr, l1c_filehandle *l1cfile, int32_t recnum);
    virtual int32_t closel1b_harp2_l1c(l1c_str *l1cstr, l1c_filehandle *l1cfile);

    virtual int32_t openl1b_misr_l1c(l1c_str *l1cstr, l1c_filehandle *l1cfile, filehandle *l1file);
    virtual int32_t readl1b_misr_l1c(l1c_str *l1cstr, l1c_filehandle *l1cfile, filehandle *l1file,
                                     int32_t recnum);
    virtual int32_t closel1b_misr_l1c(l1c_str *l1cstr, l1c_filehandle *l1cfile, filehandle *l1file);

    // global attributes
    size_t npix;
    size_t iscan;
    size_t nscan;
    size_t nviews;
    size_t nbands;
    // scan line attributes--
    double *ev_mid_time;
    size_t *scan_quality_flag;
    size_t spix;
    size_t epix;
    size_t dpix;

    // structure--pointers to data arrays
    float *Fobar;
    float *Lt_tot;  // dim depends between sensors--OCI is bands x pixels --
    float **Lt;
    float **Lt_blue;  //[num_views][num_pol][num_blue_bands][num_pixels]
    float **Lt_red;
    float **Lt_SWIR;
    float *blue_lambdas;
    float *red_lambdas;
    float *SWIR_lambdas;

    // spex
    float **I;
    float **I_polsample;
    float **I_lambdas;
    float **pol_lambdas;
    uint8_t *viewport;

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
    double *timepix;
    float *tilt;
    float *senazpix;
    float *latpix;
    float *lonpix;
    float *latpix2;  // lat2/lon are the sline +1
    float *lonpix2;
    float *latpix3;  // lat2/lon are the sline +1
    float *lonpix3;

    // harp2
    float **senazpix_3d;
    float **latpix_3d;
    float **lonpix_3d;
    float **latpix2_3d;  // lat2/lon are the sline +1
    float **lonpix2_3d;

    float *latnad;  // nadir latitude
    float *lonnad;
    float *lonershift;  // earth rotation shift
    float *terr_height;
    float *cloud_height;

    l1c_filehandle *l1cfile;
    filehandle *l1file;
};

// prototypes------
//  Open, read, close
int32_t openl1b_oci_l1c(l1c_str *l1cstr, l1c_filehandle *l1cfile);
int32_t closel1b_oci_l1c(l1c_str *l1cstr, l1c_filehandle *l1cfile);

int32_t openl1b_spex_l1c(l1c_str *l1cstr, l1c_filehandle *l1cfile);
int32_t readl1b_spex_l1c(l1c_str *l1cstr, l1c_filehandle *l1cfile, int32_t recnum);
int32_t closel1b_spex_l1c(l1c_str *l1cstr, l1c_filehandle *l1cfile);

int32_t openl1b_harp2_l1c(l1c_str *l1cstr, l1c_filehandle *l1cfile);
int32_t readl1b_harp2_l1c(l1c_str *l1cstr, l1c_filehandle *l1cfile, int32_t recnum);
int32_t closel1b_harp2_l1c(l1c_str *l1cstr, l1c_filehandle *l1cfile);

int32_t openl1b_misr_l1c(l1c_str *l1cstr, l1c_filehandle *l1cfile, filehandle *l1file);
int32_t readl1b_misr_l1c(l1c_str *l1cstr, l1c_filehandle *l1cfile, filehandle *l1file, int32_t recnum);
int32_t closel1b_misr_l1c(l1c_str *l1cstr, l1c_filehandle *l1cfile, filehandle *l1file);

}  // namespace l1c

#endif /* L1C_STR_H */
