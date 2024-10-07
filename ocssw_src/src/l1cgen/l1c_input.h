//
//  l1c_input.h
//
//
//  Created by Martin Montes on 10/21/20.
//

#ifndef L1C_INPUT_h
#define L1C_INPUT_h
#include <stdio.h>
#include "clo.h"
#include "l1c_filehandle.h"

#define DEF_FLAG                                                                                   \
    "ATMFAIL,LAND,HILT,HISATZEN,STRAYLIGHT,CLDICE,COCCOLITH,LOWLW,CHLFAIL,CHLWARN,NAVWARN,ABSAER," \
    "MAXAERITER,ATMWARN,HISOLZEN,NAVFAIL,FILTER"

#ifdef __cplusplus
extern "C" {
#endif

namespace l1c {

class L1C_input {
   protected:
   public:
    char infile[FILENAME_MAX];
    char ofile[FILENAME_MAX];
    char outlist[FILENAME_MAX];
    char l1c_grid[FILENAME_MAX];
    char l1c_anc[FILENAME_MAX];
    char start_time[50];
    char end_time[50];
    char demfile[100];
    char history[200];
    char pversion[256];
    char doi[256];
    char oformat[20];        // output file type
    char oformat_depth[20];  // output file color depth l1brsgen only

    char l2prod[2048];
    bool verbose;

    std::vector<std::string> files;  // this is c++ object?
    std::vector<std::string> files_l1c;

    float south;
    float north;
    float west;
    float east;

    // telemetry stauff----L1A--HKT
    int gransize;       // granule size in minutes
    int grantype;       // granule size in minutes
    int bintype;        // binning type 0:discrete, 1:area-wighting
    int32_t swath_num;  // swath id in one orbit, 1 or 2, ascending or descending
    int start_timeflag;

    //    char all_l1cprod[MAXPRODl1c][50];//string length rather than number of strings
    int selgran[10];    // selected granules id (not indexes) for L1C processing, up to 10 files
    int32_t ix_l1cprod[3];  // 3x1 array with selected l1c products, 1: selected
    int32_t l1c_pflag;      // l1c processing flag, 0: no, 1: yes
    float grid_resolution;  // grid resolution in km
    int32_t sensor;
    int32_t selday;
    int32_t selmon;
    int32_t selyear;
  

    int32_t projection;
    int32_t sort_method;
    int32_t cloud_height;  // 0 is ground or DEM, 1 is CTH cloud top height registration
    bool terrain_correct;  // terrain distortion correction , 1: yes
    int cloud_correct;     // cloud distortion correction , 1: yes
    int cloud_type;
    int demcloud_flag;
                        // multi  attributes (view, pol, bands)
    bool overlap_vflag;  // tells if we want merged views
    bool overlap_pflag;  // tells if we want merged polarizations
    bool overlap_bflag;  // tells if we want merged spectral bands
                         // uncertainty params l1c merged products
    int32_t unc_meth;   // uncertainity calculation method
    float unc_thres_v;  // uncertainity threshold of angular merged products as %
    float unc_thres_p;  // same but for polarization
    float unc_thres_b;  // same but for multispectral products, same view and polarization

    // constr/destruc
    L1C_input();
    virtual ~L1C_input();
    // methods
    virtual int32_t l1c_inputmain(int argc, char **argv, L1C_input *l1cinput, l1c_filehandle *l1cfile,
                                  const char *prog, const char *version);
    virtual int32_t l1c_init_options(clo_optionList_t *list, const char *prog, const char *version);
    virtual int32_t l1c_load_input(clo_optionList_t *list, L1C_input *l1cinput);
    virtual int32_t l1c_usage(const char *prog, const char *version);
    virtual int32_t l1c_input_init(L1C_input *l1cinput);
};

// prototypes l1c input
int32_t l1c_inputmain(int argc, char **argv, L1C_input *l1cinput, l1c_filehandle *l1cfile, const char *prog,
                      const char *version);
int32_t l1c_init_options(clo_optionList_t *list, const char *prog, const char *version);
int32_t l1c_load_input(clo_optionList_t *list, L1C_input *l1cinput);
int32_t l1c_usage(const char *prog, const char *version);
int32_t l1c_input_init(L1C_input *l1cinput);

#ifdef __cplusplus
}
#endif

}  // end namespace
#endif /* l1c_input_h */
