#include <iostream>
#include <stdint.h>
#include <stdlib.h>
#include <cmath>
#include <vector>
#include <netcdf>
#include <sstream>

// Chunking the arrays
#define CHUNK_CACHE_SIZE  256 * 1024 * 1024  // 256MiB of cache memory.
#define CHUNK_CACHE_NELEMS 1033
#define CHUNK_CACHE_PREEMPTION .75
#define VARCHUNK_CACHE_SIZE  4 * 1024 * 1024  // 4Mib of cache memory.
#define VARCHUNK_CACHE_NELEMS 1033
#define VARCHUNK_CACHE_PREEMPTION .75
#define CHUNKBANDS 40
#define CHUNKPIXELS 2000
#define CHUNKLINES 16
// Do something like this
//   vector<size_t> chunkVec{CHUNKBANDS,CHUNKPIXELS,CHUNKLINES};
// for dimensions less than these defaults, fix it to the dimension, e.g.:
//      if (dimBands < CHUNKBANDS)
//        chunkVec[0] = dimBands;

typedef struct {
  double master_clock;
  double MCE_clock;

  double sc_to_tilt[3][3];
  double tilt_axis[3];
  double tilt_angles[2];
  double tilt_home;
  double tilt_to_oci_mech[3][3];
  double oci_mech_to_oci_opt[3][3];
  double rta_axis[3];
  double ham_axis[3];
  double ham_at_angles[2];
  double ham_ct_angles[2];
  double rta_enc_scale;
  double ham_enc_scale;

  int32_t rta_nadir[2];

  double as_planarity[5];
  double at_planarity[5];
} geo_struct;

#define SWAP_2(x) ( (((x) & 0xff) << 8) | ((unsigned short)(x) >> 8) )

#define SWAP_4(x) ( ((x) << 24) | \
         (((x) << 8) & 0x00ff0000) | \
         (((x) >> 8) & 0x0000ff00) | \
         ((x) >> 24) )

int geolocate_oci( std::string l1a_filename, std::string geo_lut_filename,
                   geo_struct& geoLUT,
                   std::string l1b_filename, std::string dem_file,
                   bool radiance, std::string doi, const std::string ephFile, std::string pversion);

int read_mce_tlm( netCDF::NcFile *l1afile, geo_struct& geo_lut,
                  netCDF::NcGroup egid,
                  uint32_t nmcescan, uint32_t nenc, int32_t& ppr_off,
                  double& revpsec, double&secpline, int16_t& board_id,
                  int32_t *mspin, int32_t *ot_10us,
                  uint8_t *enc_count, float **hamenc, float **rtaenc,
                  int16_t &iret);

int get_ev( double secpline, int16_t *dtype, int16_t *lines, int16_t *iagg,
            uint16_t& pcdim, uint16_t& psdim, double& ev_toff,
            float *clines, float *slines, double *deltc, double *delts,
            bool dark, int16_t& iret);

int get_oci_vecs( uint32_t nscan, uint16_t pdim, double as_planarity[5],
                  double at_planarity[5], int32_t *rta_nadir,
                  double ham_ct_angles[2], double ev_toff, int32_t spin,
                  uint8_t hside, float *clines,
                  double *delt, double revpsec, int32_t ppr_off,
                  int16_t board_id, uint32_t nmcescan, int32_t *mspin,
                  uint8_t *enc_count, float **hamenc, float **rtaenc,
                  float **pview, double *theta, int16_t& iret);

int createField( netCDF::NcGroup &ncGrp, const char *sname, const char *lname, 
                   const char *standard_name, const char *units,
                   const char *description,
                   void *fill_value, const char *flag_values,
                   const char *flag_meanings,
                   double low, double high, 
                   double scale, double offset, 
                   int nt, std::vector<netCDF::NcDim>& varVec,
                   std::string coordinates);

int get_nib_nbb( uint32_t ntaps, size_t *ia, uint32_t ntb[16], int16_t jagg[16],
                 uint32_t& nib, uint32_t& nbb);

int get_agg_mat( size_t *ia, int16_t jagg[16], uint32_t ntb[16],
                 uint32_t nib, uint32_t nbb, float **amat, float **gmat);

int check_scan_times( uint32_t nscan, double *sstime, short *sfl);
