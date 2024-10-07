#ifndef _GEOLOCATE_OCI_H_
#define _GEOLOCATE_OCI_H_

#include "common.h"

#include <netcdf>

using namespace netCDF;
using namespace netCDF::exceptions;
/*
#define SWAP_2(x) ( (((x) & 0xff) << 8) | ((unsigned short)(x) >> 8) )

#define SWAP_4(x) ( ((x) << 24) | \
         (((x) << 8) & 0x00ff0000) | \
         (((x) >> 8) & 0x0000ff00) | \
         ((x) >> 24) )
*/
typedef float quat_array[4];
typedef float orb_array[3];

double constexpr PI = 3.14159265358979323846;
double constexpr RADEG = 180 / PI;
double constexpr DTOR = PI / 180;

float constexpr focal_length = 45.184;


class geoFile {
  //  netCDF::NcFile *geofile;

  std::string fileName;
  
  int ngrps;
  int ndims;
  
  netCDF::NcDim ncDims[1000];
  
 public:
  geoFile();
  ~geoFile();

  netCDF::NcFile *geofile;
  netCDF::NcGroup ncGrps[10];

  int createFile( const char* filename, const char* cdlfile,
                  size_t sdim, int *ncid, int *gid,
                  uint32_t bbb, uint32_t rbb,
                  int16_t pcdim, int16_t psdim, size_t nswband, int32_t *rta_nadir,
                  bool radiance);

  int parseDims( std::string dimString, std::vector<netCDF::NcDim>& varDims);

  int close();
};

int get_nbands( uint32_t ntaps, int16_t jagg[16], uint32_t& nbb);

int j2000_to_ecr(int32_t iyr, int32_t idy, double sec, double ecmat[3][3]);
int j2000_to_mod(int32_t iyr, int32_t idy, double sec, double j2mod[3][3]);
int get_nut(int32_t iyr, int32_t idy, double xnut[3][3]);
int get_ut1(int32_t iyr, int32_t idy, double &ut1utc);
int ephparms(double t, double &xls, double &gs, double &xlm, double &omega);
int nutate(double t, double xls, double gs, double xlm, double omega,
           double &dpsi, double &eps, double &epsm);
int gha2000(int32_t iyr, double day, double &gha);
int mtoq(double rm[3][3], double q[4]);
int qprod(double q1[4], float q2[4], double q3[4]);
int qprod(float q1[4], float q2[4], float q3[4]);
int orb_interp(size_t n_SC_rec, size_t sdim,
               double *torb, orb_array *p, orb_array *v,
               double *time, orb_array *posi, orb_array *veli);
int q_interp(size_t n_SC_rec, size_t sdim, double *tq, quat_array *q,
             double *time, quat_array *qi);
int tilt_interp(size_t n_tilts, size_t sdim, double *ttilt, float *tiltin,
                double *time, float *tilt);
int l_sun(size_t sdim, int32_t iyr, int32_t iday,
          double *sec, orb_array *sunr, double *rs);
int sun2000(size_t sdim, int32_t iyr, int32_t idy,
            double *sec, orb_array *sun, double *rs);
int qtom(float quat[4], double rm[3][3]);
int scan_ell(float p[3], double sm[3][3], double coef[10]);

int oci_geonav(const char* dem_file, float pos[3], float vel[3],
               double smat[3][3], double coef[10],
               float sunr[3], float **pview, size_t npix, double *delt,
               float *xlat, float *xlon, short *solz, short *sola,
               short *senz, short *sena, short *height,
               float& lonwest, float& loneast,
               float& latsouth, float& latnorth);

int mat2rpy(float pos[3], float vel[3], double smat[3][3], float rpy[3]);



inline
int expandEnvVar( std::string *sValue) {
  if ( (*sValue).find_first_of( "$" ) == std::string::npos) return 0;
  std::string::size_type posEndIdx = (*sValue).find_first_of( "/" );
  if ( posEndIdx == std::string::npos) return 0;
  const std::string envVar = sValue->substr (1, posEndIdx - 1);
  char *envVar_str = getenv(envVar.c_str());
  if (envVar_str == 0x0) {
    printf("Environment variable: %s not defined.\n", sValue->c_str());
    exit(1);
  }
  *sValue = envVar_str + (*sValue).substr( posEndIdx);

  return 0;
}


geoFile::geoFile() {

}


geoFile::~geoFile() {
  
}


#endif  // _GEOLOCATE_OCI_H_
