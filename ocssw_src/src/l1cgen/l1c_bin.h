//  l1c_bin.h
//  Created by Martin Montes on 1/30/2023
//

#ifndef L1C_BIN_H
#define L1C_BIN_H

#include "l1c_filehandle.h"
#include <stdbool.h>
#include <stdint.h>
#include "l1c_input.h"
#include <filehandle.h>

namespace l1c {

class bin_L1C {
   protected:
   public:
    // methods
    bin_L1C();
    virtual ~bin_L1C();

    // Open, read, close
    virtual int32_t open_binvars(bin_L1C *binstr, l1c_filehandle *l1cfile);
    virtual int32_t read_binvars(bin_L1C *binstr, l1c_filehandle *l1cfile, int32_t recnum);
    virtual int32_t close_binvars_l1c(bin_L1C *binstr, l1c_filehandle *l1cfile);

    // global attributes
    int16_t num_gridlines;
    int16_t nbinx;
    size_t **nrec_2D;   // row/col
    size_t ***nrec_3D;  // row/col/view

    // OCIS binned
    float **diff_h2;          // row/col
    float ***sca_3D;          // row/col/view
    float ****QC_bitwise_4D;  // row/col/view/bands
    float ***QC_3D;           // row/col/view
    float ****I_4D;           // row/col/view/bands
    float ****I_noise_4D;     // #pixels

    // OCIS line by line
    short *obs_per_view;
    float *QC_bitwise;
    float *QC;
    float *I;
    float *I_noise;

    filehandle *l1file;
};

// prototypes------
//  Open, read, close
int32_t open_binvars(bin_L1C *binstr, l1c_filehandle *l1cfile);
int32_t read_binvars(bin_L1C *binstr, l1c_filehandle *l1cfile, int32_t recnum);
int32_t close_binvars_l1c(bin_L1C *binstr, l1c_filehandle *l1cfile);

}  // namespace l1c

#endif /* L1C_BIN_H */
