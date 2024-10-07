#include <fstream>
#include <iostream>
#include <sstream>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix_double.h>

#include "geolocate_hawkeye.h"
#include "hawkeyeUtil.h"
#include "netcdf.h"

#define VERSION "0.772"

//    Modification history:
//  Programmer     Organization   Date     Ver   Description of change
//  ----------     ------------   ----     ---   ---------------------
//  Joel Gales     FutureTech     05/21/18 0.10  Original development
//  Joel Gales     FutureTech     10/18/18 0.20  Change offset in alp
//  Joel Gales     FutureTech     10/23/18 0.30  Add support for rpy
//  Joel Gales     SAIC           11/02/18 0.40  Initialize output arrays
//                                               to fill values.  Skip rpy
//                                               processing if stime = -999.
//
//  Joel Gales     SAIC           11/06/18 0.50  Set bit 2 (=4) of qfl when
//                                               stime = -999
//  Gwyn Fireman   SAIC           02/21/20 0.55  Allow different att, orb dims
//  Joel Gales     SAIC           03/21/20 0.56  Pass sdim to createFile() to
//                                               set number_of_scans
//  Liang Hong     SAIC           11/02/20 0.561 quatr memory allocation correction
//  Liang Hong     SAIC           02/24/21 0.661 added time and roll offset
//  Liang Hong     SAIC           07/29/21 0.761 corrected center-pixel offset and
//                                               roll offset, added static pitch offset
//  Liang Hong     SAIC           01/12/22 0.771 added scale and offset attributes
//  Liang Hong     SAIC           04/06/22 0.772 default output file name if not specified

using namespace std;

int main(int argc, char *argv[]) {
  string l1a_name;
  string geo_name;
  
  cout << "geolocate_hawkeye " << VERSION << " ("
       << __DATE__ << " " << __TIME__ << ")" << endl;

  if (argc == 1) {
    cout << endl
         << "geolocate_hawkeye input_l1a_filename output_geo_filename" << endl;

    return 0;
  } else if (argc == 2) {
    l1a_name.assign(argv[1]);
    string str_L1A("L1A");
    size_t L1A_found = l1a_name.find(str_L1A);
    if ((L1A_found == string::npos)) {
		geo_name.assign(l1a_name+".GEO");
    } else {
    	geo_name = l1a_name;
    	geo_name.replace(L1A_found,str_L1A.length(),"GEO");
    }
  } else {
  	l1a_name.assign(argv[1]);
  	geo_name.assign(argv[2]);
  }
  
  int status;
  int grpid, varid, dimid;
  int dimids[NC_MAX_VAR_DIMS];

  
  int l1a_ncid;
  int l1a_ngrps;
  int l1a_gid[10];
  size_t n_att_rec, n_orb_rec, sdim, pdim, t_len; 
  nc_type t_type;
  float rolloff, timeoff;

  enum l1a_grps { scan_attr,
                  telem,
                  navigation,
                  earth_view };
  enum geo_grps { scan_attr_geo,
                  nav_geo,
                  geolocation };

  status = nc_open(l1a_name.c_str(), NC_NOWRITE, &l1a_ncid);
  check_err(status, __LINE__, __FILE__);

  status = nc_inq_grps(l1a_ncid, &l1a_ngrps, l1a_gid);
  check_err(status, __LINE__, __FILE__);

  // number of pixels
  status = nc_inq_dimid(l1a_ncid, "number_of_pixels", &dimid);
  check_err(status, __LINE__, __FILE__);
  nc_inq_dimlen(l1a_ncid, dimid, &pdim);

  grpid = l1a_gid[navigation];  // navigation_data group

  // att time
  status = nc_inq_varid(grpid, "att_time", &varid);
  check_err(status, __LINE__, __FILE__);
  nc_inq_vardimid(grpid, varid, dimids);
  nc_inq_dimlen(l1a_ncid, dimids[0], &n_att_rec);

  double *att_time = new double[n_att_rec]();
  status = nc_get_var_double(grpid, varid, att_time);
  check_err(status, __LINE__, __FILE__);

  // att quat
  status = nc_inq_varid(grpid, "att_quat", &varid);
  check_err(status, __LINE__, __FILE__);

  quat_array *att_quat = new quat_array[n_att_rec]();
  status = nc_get_var_float(grpid, varid, &att_quat[0][0]);
  check_err(status, __LINE__, __FILE__);

  // orb time
  status = nc_inq_varid(grpid, "orb_time", &varid);
  check_err(status, __LINE__, __FILE__);
  nc_inq_vardimid(grpid, varid, dimids);
  nc_inq_dimlen(l1a_ncid, dimids[0], &n_orb_rec);

  double *orb_time = new double[n_orb_rec]();
  status = nc_get_var_double(grpid, varid, orb_time);
  check_err(status, __LINE__, __FILE__);

  // orb pos
  status = nc_inq_varid(grpid, "orb_pos", &varid);
  check_err(status, __LINE__, __FILE__);

  orb_array *orb_pos = new orb_array[n_orb_rec]();
  status = nc_get_var_float(grpid, varid, &orb_pos[0][0]);
  check_err(status, __LINE__, __FILE__);

  // orb vel
  status = nc_inq_varid(grpid, "orb_vel", &varid);
  check_err(status, __LINE__, __FILE__);

  orb_array *orb_vel = new orb_array[n_orb_rec]();
  status = nc_get_var_float(grpid, varid, &orb_vel[0][0]);
  check_err(status, __LINE__, __FILE__);
  
  // Check for time and roll offset attributes
  status = nc_inq_att(grpid, NC_GLOBAL, "roll_offset", &t_type, &t_len);
  if (status != NC_NOERR) {
      rolloff = 0.0;
  } else {
      status = nc_get_att_float(grpid, NC_GLOBAL, "roll_offset", &rolloff);  
  } 
  cout << "roll_offset=" <<rolloff << endl;
  
  status = nc_inq_att(grpid, NC_GLOBAL, "time_offset", &t_type, &t_len);
  if (status != NC_NOERR) {
      timeoff = 0.0;
  } else {
      status = nc_get_att_float(grpid, NC_GLOBAL, "time_offset", &timeoff);  
  } 
  cout << "time_offset=" <<timeoff << endl;

  // scan_time
  grpid = l1a_gid[scan_attr];  // scan_line_attributes group
  status = nc_inq_varid(grpid, "scan_time", &varid);
  check_err(status, __LINE__, __FILE__);
  nc_inq_vardimid(grpid, varid, dimids);
  nc_inq_dimlen(l1a_ncid, dimids[0], &sdim);

  double *stime = new double[sdim]();
  status = nc_get_var_double(grpid, varid, stime);
  check_err(status, __LINE__, __FILE__);

  // Get L1A start date
  char tstart[25] = "";
  char tend[25] = "";
  status = nc_get_att_text(l1a_ncid, NC_GLOBAL, "time_coverage_start", tstart);
  check_err(status, __LINE__, __FILE__);

  status = nc_get_att_text(l1a_ncid, NC_GLOBAL, "time_coverage_end", tend);
  check_err(status, __LINE__, __FILE__);

  int32_t iyr, iday, msec;
  isodate2ydmsec(tstart, &iyr, &iday, &msec);

  // Transform orbit and attitude from J2000 to ECR
  quat_array *quatr = new quat_array[n_att_rec]();
  orb_array *posr = new orb_array[n_orb_rec]();
  orb_array *velr = new orb_array[n_orb_rec]();

  double omegae = 7.29211585494e-5;
  double ecmat[3][3];

  // Orbit
  for (size_t i = 0; i < n_orb_rec; i++) {
    j2000_to_ecr(iyr, iday, orb_time[i], ecmat);

    for (size_t j = 0; j < 3; j++) {
      posr[i][j] = ecmat[j][0] * orb_pos[i][0] +
                   ecmat[j][1] * orb_pos[i][1] +
                   ecmat[j][2] * orb_pos[i][2];
      velr[i][j] = ecmat[j][0] * orb_vel[i][0] +
                   ecmat[j][1] * orb_vel[i][1] +
                   ecmat[j][2] * orb_vel[i][2];
    }
    velr[i][0] += posr[i][1] * omegae;
    velr[i][1] -= posr[i][0] * omegae;
  }  // i loop

  // Attitude
  for (size_t i = 0; i < n_att_rec; i++) {
    double ecq[4], qt2[4];
    float qt1[4];
    j2000_to_ecr(iyr, iday, att_time[i], ecmat);
    mtoq(ecmat, ecq);

    memcpy(qt1, &att_quat[i][1], 3 * sizeof(float));
    qt1[3] = att_quat[i][0];

    qprod(ecq, qt1, qt2);
    for (size_t j = 0; j < 4; j++) quatr[i][j] = qt2[j];
  }  // i loop

  // Interpolate orbit and attitude to scan times
  quat_array *quati = new quat_array[sdim]();
  q_interp(n_att_rec, sdim, att_time, quatr, stime, quati);
  
  // apply offsets to orbit times
  for (size_t i = 0; i < sdim; i++) {
      stime[i] = stime[i] + timeoff;
  }
  
  orb_array *posi = new orb_array[sdim]();
  orb_array *veli = new orb_array[sdim]();
  orb_array *rpy = new orb_array[sdim]();
  orb_array *ang = new orb_array[sdim]();
  orb_interp(n_orb_rec, sdim, orb_time, posr, velr, stime, posi, veli);



  float *xlon = new float[sdim * pdim]();
  float *xlat = new float[sdim * pdim]();
  short *solz = new short[sdim * pdim]();
  short *sola = new short[sdim * pdim]();
  short *senz = new short[sdim * pdim]();
  short *sena = new short[sdim * pdim]();
  short *range = new short[sdim * pdim]();
  uint8_t *qfl = new uint8_t[sdim * pdim]();

  // Initialize output arrays
  for (size_t i = 0; i < sdim * pdim; i++) {
    xlon[i] = -999.9;
    xlat[i] = -999.9;

    solz[i] = -32768;
    sola[i] = -32768;
    senz[i] = -32768;
    sena[i] = -32768;
    range[i] = -32768;
  }

  // Generate pointing vector array
  orb_array *pview = new orb_array[pdim]();
  for (size_t i = 0; i < pdim; i++) {
    double alp = atan(((double)pdim / 2 - i + 40) * 0.01 / focal_length);   // Ver 0.761
    pview[i][0] = 0.0;
    pview[i][1] = sin(alp);
    pview[i][2] = cos(alp);
  }

  // Get Sun vectors
  orb_array *sunr = new orb_array[sdim]();
  l_sun(sdim, iyr, iday, stime, sunr);

  // Get S/C-to-sensor matrix (placeholder for now)
  double sc_to_hwk[3][3];
  for (size_t i = 0; i < 3; i++)
    for (size_t j = 0; j < 3; j++) sc_to_hwk[i][j] = 0.0;
  // double sr = sin(rolloff/RADEG);
  // double cr = cos(rolloff/RADEG);
  sc_to_hwk[0][2] = 1.0;
  // sc_to_hwk[1][0] = -sr;
  sc_to_hwk[1][1] = 1.0; //cr;
  sc_to_hwk[2][0] = -1.0; //-cr;
  // sc_to_hwk[2][1] = -sr;
  float pitchoff = -0.8;
  float rpyoff[3] = {rolloff,pitchoff,0.0};  // Ver 0.761

  // Geolocate each scan line
  gsl_matrix *smat = gsl_matrix_alloc(3, 3);

  for (size_t i = 0; i < sdim; i++) {
    if (stime[i] > -900.0) {  // replaced if (stime[i] != -999.0) due to stime + timeoff shift
      // Convert quaternion to matrix
      double om[3][3];
      double rpym[3][3];
      double qmat[3][3];
      qtom(quati[i], qmat);
    
      gsl_matrix_view A = gsl_matrix_view_array(&sc_to_hwk[0][0], 3, 3);
      gsl_matrix_view B = gsl_matrix_view_array(&qmat[0][0], 3, 3);

      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0,
                     &A.matrix, &B.matrix, 0.0, smat);

	  // Compute and adjust attitude angles
	  double *ptr_C = gsl_matrix_ptr(smat, 0, 0);
	  mat2rpy(posi[i], veli[i], (double(*)[3])ptr_C, rpy[i], om);

	  for (size_t j = 0; j < 3; j++) ang[i][j] = rpy[i][j] + rpyoff[j]; 

	  // Get sensor orientation matrix and scan ellipse coefficients
	  euler(ang[i],rpym);

	  //smat = rpym#om
      gsl_matrix_view C = gsl_matrix_view_array(&rpym[0][0], 3, 3);
      gsl_matrix_view D = gsl_matrix_view_array(&om[0][0], 3, 3);

      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0,
                     &C.matrix, &D.matrix, 0.0, smat);
      
      double coef[10];
      ptr_C = gsl_matrix_ptr(smat, 0, 0);
      scan_ell(posi[i], (double(*)[3])ptr_C, coef);	

      // Geolocate pixels
      size_t ip = i * pdim;
      uni_geonav(posi[i], veli[i], (double(*)[3])ptr_C, coef,
                 sunr[i], pview, pdim, &xlat[ip], &xlon[ip],
                 &solz[ip], &sola[ip], &senz[ip], &sena[ip], &range[ip]);
      // mat2rpy(posi[i], veli[i], (double(*)[3])ptr_C, rpy[i]);
    } else {
      qfl[i] |= 4;
    }
  }  // scan loop
  gsl_matrix_free(smat);

  int geo_ncid;
  int geo_gid[10];
  createFile(geo_name.c_str(),
             "$OCDATAROOT/hawkeye/Hawkeye_GEO_Data_Structure.cdl",
             sdim, &geo_ncid, geo_gid);

  char buf[32];
  strcpy(buf, unix2isodate(now(), 'G'));
  nc_put_att_text(geo_ncid, NC_GLOBAL, "date_created", strlen(buf), buf);

  string varname;

  varname.assign("scan_time");
  status = nc_inq_varid(geo_gid[scan_attr_geo], varname.c_str(), &varid);
  status = nc_put_var_double(geo_gid[scan_attr_geo], varid, stime);
  check_err(status, __LINE__, __FILE__);

  status = nc_put_att_text(geo_ncid, NC_GLOBAL,
                           "time_coverage_start", strlen(tstart), tstart);

  status = nc_put_att_text(geo_ncid, NC_GLOBAL,
                           "time_coverage_end", strlen(tend), tend);

  varname.assign("att_quat");
  status = nc_inq_varid(geo_gid[nav_geo], varname.c_str(), &varid);
  status = nc_put_var_float(geo_gid[nav_geo], varid, (float *)quati);
  check_err(status, __LINE__, __FILE__);

  varname.assign("att_ang");
  status = nc_inq_varid(geo_gid[nav_geo], varname.c_str(), &varid);
  status = nc_put_var_float(geo_gid[nav_geo], varid, (float *)ang);
  check_err(status, __LINE__, __FILE__);

  varname.assign("orb_pos");
  status = nc_inq_varid(geo_gid[nav_geo], varname.c_str(), &varid);
  status = nc_put_var_float(geo_gid[nav_geo], varid, (float *)posi);
  check_err(status, __LINE__, __FILE__);

  varname.assign("orb_vel");
  status = nc_inq_varid(geo_gid[nav_geo], varname.c_str(), &varid);
  status = nc_put_var_float(geo_gid[nav_geo], varid, (float *)veli);
  check_err(status, __LINE__, __FILE__);

  varname.assign("sun_ref");
  status = nc_inq_varid(geo_gid[nav_geo], varname.c_str(), &varid);
  status = nc_put_var_float(geo_gid[nav_geo], varid, (float *)sunr);
  check_err(status, __LINE__, __FILE__);

  varname.assign("sc_to_hawkeye");
  status = nc_inq_varid(geo_gid[nav_geo], varname.c_str(), &varid);
  status = nc_put_var_double(geo_gid[nav_geo], varid, (double *)sc_to_hwk);
  check_err(status, __LINE__, __FILE__);

  varname.assign("latitude");
  status = nc_inq_varid(geo_gid[geolocation], varname.c_str(), &varid);
  status = nc_put_var_float(geo_gid[geolocation], varid, xlat);
  check_err(status, __LINE__, __FILE__);

  varname.assign("longitude");
  status = nc_inq_varid(geo_gid[geolocation], varname.c_str(), &varid);
  status = nc_put_var_float(geo_gid[geolocation], varid, xlon);
  check_err(status, __LINE__, __FILE__);

  varname.assign("range");
  status = nc_inq_varid(geo_gid[geolocation], varname.c_str(), &varid);
  status = nc_put_var_short(geo_gid[geolocation], varid, range);
  check_err(status, __LINE__, __FILE__);

  varname.assign("quality_flag");
  status = nc_inq_varid(geo_gid[geolocation], varname.c_str(), &varid);
  status = nc_put_var_ubyte(geo_gid[geolocation], varid, qfl);
  check_err(status, __LINE__, __FILE__);

  varname.assign("sensor_azimuth");
  status = nc_inq_varid(geo_gid[geolocation], varname.c_str(), &varid);
  status = nc_put_var_short(geo_gid[geolocation], varid, sena);
  check_err(status, __LINE__, __FILE__);

  varname.assign("sensor_zenith");
  status = nc_inq_varid(geo_gid[geolocation], varname.c_str(), &varid);
  status = nc_put_var_short(geo_gid[geolocation], varid, senz);
  check_err(status, __LINE__, __FILE__);

  varname.assign("solar_azimuth");
  status = nc_inq_varid(geo_gid[geolocation], varname.c_str(), &varid);
  status = nc_put_var_short(geo_gid[geolocation], varid, sola);
  check_err(status, __LINE__, __FILE__);

  varname.assign("solar_zenith");
  status = nc_inq_varid(geo_gid[geolocation], varname.c_str(), &varid);
  status = nc_put_var_short(geo_gid[geolocation], varid, solz);
  check_err(status, __LINE__, __FILE__);

  delete[](att_time);
  delete[](orb_time);
  delete[](att_quat);
  delete[](orb_pos);
  delete[](orb_vel);
  delete[](stime);
  delete[](quatr);
  delete[](posr);
  delete[](velr);
  delete[](posi);
  delete[](veli);
  delete[](rpy);
  delete[](quati);
  delete[](xlon);
  delete[](xlat);
  delete[](solz);
  delete[](sola);
  delete[](senz);
  delete[](sena);
  delete[](range);
  delete[](qfl);
  delete[](pview);
  delete[](sunr);

  return 0;
}

int j2000_to_ecr(int32_t iyr, int32_t idy, double sec, double ecmat[3][3]) {
  // Get J2000 to ECEF transformation matrix

  // Arguments:

  // Name               Type    I/O     Description
  // --------------------------------------------------------
  // iyr         I       I      Year, four digits
  // idy         I       I      Day of year
  // sec        R*8      I      Seconds of day
  // ecmat(3,3)  R       O      J2000 to ECEF matrix

  // Get transformation from J2000 to mean-of-date inertial
  double j2mod[3][3];
  j2000_to_mod(iyr, idy, sec, j2mod);

  // Get nutation and UT1-UTC (once per run)
  double xnut[3][3], ut1utc;
  get_nut(iyr, idy, xnut);
  get_ut1(iyr, idy, ut1utc);

  // Compute Greenwich hour angle for time of day
  double day = idy + (sec + ut1utc) / 86400;
  double gha, gham[3][3];
  gha2000(iyr, day, gha);

  gham[0][0] = cos(gha / RADEG);
  gham[1][1] = cos(gha / RADEG);
  gham[2][2] = 1;
  gham[0][1] = sin(gha / RADEG);
  gham[1][0] = -sin(gha / RADEG);

  gham[0][2] = 0;
  gham[2][0] = 0;
  gham[1][2] = 0;
  gham[2][1] = 0;

  // Combine all transformations
  gsl_matrix_view A = gsl_matrix_view_array(&gham[0][0], 3, 3);  // gham
  gsl_matrix_view B = gsl_matrix_view_array(&xnut[0][0], 3, 3);  // xnut
  gsl_matrix *C = gsl_matrix_alloc(3, 3);

  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, &A.matrix, &B.matrix, 0.0, C);

  gsl_matrix_view D = gsl_matrix_view_array(&j2mod[0][0], 3, 3);  // j2mod
  gsl_matrix *E = gsl_matrix_alloc(3, 3);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, C, &D.matrix, 0.0, E);
  double *ptr_E = gsl_matrix_ptr(E, 0, 0);

  memcpy(ecmat, ptr_E, 9 * sizeof(double));

  gsl_matrix_free(C);
  gsl_matrix_free(E);

  return 0;
}

int j2000_to_mod(int32_t iyr, int32_t idy, double sec, double j2mod[3][3]) {
  // Get J2000 to MOD (precession) transformation

  // Arguments:

  // Name               Type    I/O     Description
  // --------------------------------------------------------
  // iyr         I       I      Year, four digits
  // idy         I       I      Day of year
  // sec        R*8      I      Seconds of day
  // j2mod(3,3)  R       O      J2000 to MOD matrix

  int16_t iyear = iyr;
  int16_t iday = idy;

  double t = jday(iyear, 1, iday) - (double)2451545.5 + sec / 86400;
  t /= 36525;

  double zeta0 = t * (2306.2181 + 0.302 * t + 0.018 * t * t) / RADEG / 3600;
  double thetap = t * (2004.3109 - 0.4266 * t - 0.04160 * t * t) / RADEG / 3600;
  double xip = t * (2306.2181 + 1.095 * t + 0.018 * t * t) / RADEG / 3600;

  j2mod[0][0] = -sin(zeta0) * sin(xip) + cos(zeta0) * cos(xip) * cos(thetap);
  j2mod[0][1] = -cos(zeta0) * sin(xip) - sin(zeta0) * cos(xip) * cos(thetap);
  j2mod[0][2] = -cos(xip) * sin(thetap);
  j2mod[1][0] = sin(zeta0) * cos(xip) + cos(zeta0) * sin(xip) * cos(thetap);
  j2mod[1][1] = cos(zeta0) * cos(xip) - sin(zeta0) * sin(xip) * cos(thetap);
  j2mod[1][2] = -sin(xip) * sin(thetap);
  j2mod[2][0] = cos(zeta0) * sin(thetap);
  j2mod[2][1] = -sin(zeta0) * sin(thetap);
  j2mod[2][2] = cos(thetap);

  return 0;
}

int get_nut(int32_t iyr, int32_t idy, double xnut[3][3]) {
  int16_t iyear = iyr;
  int16_t iday = idy;

  double t = jday(iyear, 1, iday) - (double)2451545.5;

  double xls, gs, xlm, omega;
  double dpsi, eps, epsm;
  ephparms(t, xls, gs, xlm, omega);
  nutate(t, xls, gs, xlm, omega, dpsi, eps, epsm);

  xnut[0][0] = cos(dpsi / RADEG);
  xnut[1][0] = -sin(dpsi / RADEG) * cos(epsm / RADEG);
  xnut[2][0] = -sin(dpsi / RADEG) * sin(epsm / RADEG);
  xnut[0][1] = sin(dpsi / RADEG) * cos(eps / RADEG);
  xnut[1][1] = cos(dpsi / RADEG) * cos(eps / RADEG) * cos(epsm / RADEG) +
               sin(eps / RADEG) * sin(epsm / RADEG);
  xnut[2][1] = cos(dpsi / RADEG) * cos(eps / RADEG) * sin(epsm / RADEG) -
               sin(eps / RADEG) * cos(epsm / RADEG);
  xnut[0][2] = sin(dpsi / RADEG) * sin(eps / RADEG);
  xnut[1][2] = cos(dpsi / RADEG) * sin(eps / RADEG) * cos(epsm / RADEG) -
               cos(eps / RADEG) * sin(epsm / RADEG);
  xnut[2][2] = cos(dpsi / RADEG) * sin(eps / RADEG) * sin(epsm / RADEG) +
               cos(eps / RADEG) * cos(epsm / RADEG);

  return 0;
}

int ephparms(double t, double &xls, double &gs, double &xlm, double &omega) {
  //  This subroutine computes ephemeris parameters used by other Mission
  //  Operations routines:  the solar mean longitude and mean anomaly, and
  //  the lunar mean longitude and mean ascending node.  It uses the model
  //  referenced in The Astronomical Almanac for 1984, Section S
  //  (Supplement) and documented for the SeaWiFS Project in "Constants
  //  and Parameters for SeaWiFS Mission Operations", in TBD.  These
  //  parameters are used to compute the solar longitude and the nutation
  //  in longitude and obliquity.

  //  Sun Mean Longitude
  xls = (double)280.46592 + ((double)0.9856473516) * t;
  xls = fmod(xls, (double)360);

  //  Sun Mean Anomaly
  gs = (double)357.52772 + ((double)0.9856002831) * t;
  gs = fmod(gs, (double)360);

  //  Moon Mean Longitude
  xlm = (double)218.31643 + ((double)13.17639648) * t;
  xlm = fmod(xlm, (double)360);

  //  Ascending Node of Moon's Mean Orbit
  omega = (double)125.04452 - ((double)0.0529537648) * t;
  omega = fmod(omega, (double)360);

  return 0;
}

int nutate(double t, double xls, double gs, double xlm, double omega,
           double &dpsi, double &eps, double &epsm) {
  //  This subroutine computes the nutation in longitude and the obliquity
  //  of the ecliptic corrected for nutation.  It uses the model referenced
  //  in The Astronomical Almanac for 1984, Section S (Supplement) and
  //  documented for the SeaWiFS Project in "Constants and Parameters for
  //  SeaWiFS Mission Operations", in TBD.  These parameters are used to
  //  compute the apparent time correction to the Greenwich Hour Angle and
  //  for the calculation of the geocentric Sun vector.  The input
  //  ephemeris parameters are computed using subroutine ephparms.  Terms
  //  are included to 0.1 arcsecond.

  //  Nutation in Longitude
  dpsi = -17.1996 * sin(omega / RADEG) + 0.2062 * sin(2. * omega / RADEG) -
         1.3187 * sin(2. * xls / RADEG) + 0.1426 * sin(gs / RADEG) -
         0.2274 * sin(2. * xlm / RADEG);

  //  Mean Obliquity of the Ecliptic
  epsm = (double)23.439291 - ((double)3.560e-7) * t;

  //  Nutation in Obliquity
  double deps = 9.2025 * cos(omega / RADEG) + 0.5736 * cos(2. * xls / RADEG);

  //  True Obliquity of the Ecliptic
  eps = epsm + deps / 3600;

  dpsi = dpsi / 3600;

  return 0;
}

int get_ut1(int32_t iyr, int32_t idy, double &ut1utc) {
  int16_t iyear = iyr;
  int16_t iday = idy;

  static int32_t ijd[25000];
  static double ut1[25000];
  string utcpole = "$OCVARROOT/modis/utcpole.dat";
  static bool first = true;
  int k = 0;
  if (first) {
    string line;
    expandEnvVar(&utcpole);
    istringstream istr;

    ifstream utcpfile(utcpole.c_str());
    if (utcpfile.is_open()) {
      getline(utcpfile, line);
      getline(utcpfile, line);
      while (getline(utcpfile, line)) {
        //        cout << line << '\n';
        istr.clear();
        istr.str(line.substr(0, 5));
        istr >> ijd[k];
        istr.clear();
        istr.str(line.substr(44, 9));
        istr >> ut1[k];
        k++;
      }
      ijd[k] = 0;
      utcpfile.close();
      first = false;
    } else {
      cout << utcpole.c_str() << " not found" << endl;
      exit(1);
    }
  }  // if (first)

  k = 0;
  int mjd = jday(iyear, 1, iday) - 2400000;
  while (ijd[k] > 0) {
    if (mjd == ijd[k]) {
      ut1utc = ut1[k];
      break;
    }
    k++;
  }

  return 0;
}

int gha2000(int32_t iyr, double day, double &gha) {
  //  This subroutine computes the Greenwich hour angle in degrees for the
  //  input time.  It uses the model referenced in The Astronomical Almanac
  //  for 1984, Section S (Supplement) and documented for the SeaWiFS
  //  Project in "Constants and Parameters for SeaWiFS Mission Operations",
  //  in TBD.  It includes the correction to mean sidereal time for nutation
  //  as well as precession.
  //

  //  Compute days since J2000
  int16_t iday = day;
  double fday = day - iday;
  int jd = jday(iyr, 1, iday);
  double t = jd - (double)2451545.5 + fday;

  //  Compute Greenwich Mean Sidereal Time      (degrees)
  double gmst = (double)100.4606184 + (double)0.9856473663 * t +
                (double)2.908e-13 * t * t;

  //  Check if need to compute nutation correction for this day
  double xls, gs, xlm, omega;
  double dpsi, eps, epsm;
  ephparms(t, xls, gs, xlm, omega);
  nutate(t, xls, gs, xlm, omega, dpsi, eps, epsm);

  //  Include apparent time correction and time-of-day
  gha = gmst + dpsi * cos(eps / RADEG) + fday * 360;
  gha = fmod(gha, (double)360);

  return 0;
}

int euler(float a[3], double xm[3][3]) {
	//  Computes coordinate transformation matrix corresponding to Euler 
	//   sequence; assumes order of rotations is X, Y, Z
	//
	//  Reference:  Wertz, Appendix E

	double xm1[3][3];
	double xm2[3][3];
	double xm3[3][3];
	double xmm[3][3];
	
	for(size_t i=0;i<3;i++) {
		for (size_t j=0;j<3;j++) {
			xm1[i][j]=0;
			xm2[i][j]=0;
			xm3[i][j]=0;
			xmm[i][j]=0;
			xm[i][j]=0;
		}
	}

	double c1=cos(a[0]/RADEG);
	double s1=sin(a[0]/RADEG);
	double c2=cos(a[1]/RADEG);
	double s2=sin(a[1]/RADEG);
	double c3=cos(a[2]/RADEG);
	double s3=sin(a[2]/RADEG);

	// Convert individual rotations to matrices
	xm1[0][0]=1.0;
	xm1[1][1]=c1;
	xm1[2][2]=c1;
	xm1[1][2]=s1;
	xm1[2][1]=-s1;
	xm2[1][1]=1.0;
	xm2[0][0]=c2;
	xm2[2][2]=c2;
	xm2[2][0]=s2;
	xm2[0][2]=-s2;
	xm3[2][2]=1.0;
	xm3[1][1]=c3;
	xm3[0][0]=c3;
	xm3[0][1]=s3;
	xm3[1][0]=-s3;

	// Compute total rotation as xm3*xm2*xm1
	//xmm=xm2#xm1	
	gsl_matrix_view A = gsl_matrix_view_array(&xm2[0][0], 3, 3);
	gsl_matrix_view B = gsl_matrix_view_array(&xm1[0][0], 3, 3);
	gsl_matrix *E = gsl_matrix_alloc(3, 3);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0,
                     &A.matrix, &B.matrix, 0.0, E);
    double *ptr_E = gsl_matrix_ptr(E, 0, 0);
    memcpy(xmm, ptr_E, 9 * sizeof(double));

	//xm=xm3#xmm
	gsl_matrix_view C = gsl_matrix_view_array(&xm3[0][0], 3, 3);
	gsl_matrix_view D = gsl_matrix_view_array(&xmm[0][0], 3, 3);
	gsl_matrix *F = gsl_matrix_alloc(3, 3);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0,
                     &C.matrix, &D.matrix, 0.0, F);
    double *ptr_F = gsl_matrix_ptr(F, 0, 0);
    memcpy(xm, ptr_F, 9 * sizeof(double));                     

    gsl_matrix_free(E);
    gsl_matrix_free(F);
	return 0;
}

int mtoq(double rm[3][3], double q[4]) {
  //  Convert direction cosine matrix to equivalent quaternion

  double e[3];

  //  Compute Euler angle
  double phi;
  double cphi = (rm[0][0] + rm[1][1] + rm[2][2] - 1) / 2;
  if (fabs(cphi) < 0.98) {
    phi = acos(cphi);
  } else {
    double ssphi = (pow(rm[0][1] - rm[1][0], 2) +
                    pow(rm[2][0] - rm[0][2], 2) +
                    pow(rm[1][2] - rm[2][1], 2)) /
                   4;
    phi = asin(sqrt(ssphi));
    if (cphi < 0) phi = PI - phi;
  }

  //  Compute Euler axis
  e[0] = (rm[2][1] - rm[1][2]) / (sin(phi) * 2);
  e[1] = (rm[0][2] - rm[2][0]) / (sin(phi) * 2);
  e[2] = (rm[1][0] - rm[0][1]) / (sin(phi) * 2);
  double norm = sqrt(e[0] * e[0] + e[1] * e[1] + e[2] * e[2]);
  e[0] = e[0] / norm;
  e[1] = e[1] / norm;
  e[2] = e[2] / norm;

  //  Compute quaternion
  q[0] = e[0] * sin(phi / 2);
  q[1] = e[1] * sin(phi / 2);
  q[2] = e[2] * sin(phi / 2);
  q[3] = cos(phi / 2);

  return 0;
}

int qprod(double q1[4], float q2[4], double q3[4]) {
  // Compute the product of two quaternions q3 = q1*q2

  q3[0] = q1[0] * q2[3] + q1[1] * q2[2] - q1[2] * q2[1] + q1[3] * q2[0];
  q3[1] = -q1[0] * q2[2] + q1[1] * q2[3] + q1[2] * q2[0] + q1[3] * q2[1];
  q3[2] = q1[0] * q2[1] - q1[1] * q2[0] + q1[2] * q2[3] + q1[3] * q2[2];
  q3[3] = -q1[0] * q2[0] - q1[1] * q2[1] - q1[2] * q2[2] + q1[3] * q2[3];

  return 0;
}

int qprod(float q1[4], float q2[4], float q3[4]) {
  // Compute the product of two quaternions q3 = q1*q2

  q3[0] = q1[0] * q2[3] + q1[1] * q2[2] - q1[2] * q2[1] + q1[3] * q2[0];
  q3[1] = -q1[0] * q2[2] + q1[1] * q2[3] + q1[2] * q2[0] + q1[3] * q2[1];
  q3[2] = q1[0] * q2[1] - q1[1] * q2[0] + q1[2] * q2[3] + q1[3] * q2[2];
  q3[3] = -q1[0] * q2[0] - q1[1] * q2[1] - q1[2] * q2[2] + q1[3] * q2[3];

  return 0;
}

int orb_interp(size_t n_orb_rec, size_t sdim, double *torb, orb_array *p,
               orb_array *v, double *time, orb_array *posi, orb_array *veli) {
  //  Purpose: Interpolate orbit position and velocity vectors to scan line
  //  times
  //
  //
  //  Calling Arguments:
  //
  //  Name         Type    I/O     Description
  //  --------     ----    ---     -----------
  //  torb(*)     double   I      Array of orbit vector times in seconds of day
  //  p(3,*)       float   I      Array of orbit position vectors for
  //                                each time torb
  //  v(3,*)       float   I      Array of orbit velocity vectors for
  //                                each time torb
  //  time(*)     double   I      Array of time in seconds of day
  //                               for every scan line
  //  pi(3,*)     float    O      Array of interpolated positions
  //  vi(3,*)     float    O      Array of interpolated velocities
  //
  //
  //  By: Frederick S. Patt, GSC, August 10, 1993
  //
  //  Notes:  Method uses cubic polynomial to match positions and velocities
  //   at input data points.
  //
  //  Modification History:
  //
  //  Created IDL version from FORTRAN code.  F.S. Patt, SAIC, November 29, 2006
  //

  double a0[3], a1[3], a2[3], a3[3];

  /*
  //  Make sure that first orbit vector precedes first scan line
  k = where (time lt torb(0))
  if (k(0) ne -1) then begin
     posi(*,k) = 0.0
     veli(*,k) = 0.0
     orbfl(k) = 1
     print, 'Scan line times before available orbit data'
     i1 = max(k) + 1
  endif
  */

  size_t ind = 0;
  for (size_t i = 0; i < sdim; i++) {
    //  Find input orbit vectors bracketing scan
    for (size_t j = ind; j < n_orb_rec; j++) {
      if (time[i] > torb[j] && time[i] <= torb[j + 1]) {
        ind = j;
        break;
      }
    }

    //  Set up cubic interpolation
    double dt = torb[ind + 1] - torb[ind];
    for (size_t j = 0; j < 3; j++) {
      a0[j] = p[ind][j];
      a1[j] = v[ind][j] * dt;
      a2[j] = 3 * p[ind + 1][j] - 3 * p[ind][j] - 2 * v[ind][j] * dt -
              v[ind + 1][j] * dt;
      a3[j] = 2 * p[ind][j] - 2 * p[ind + 1][j] + v[ind][j] * dt +
              v[ind + 1][j] * dt;
    }

    //  Interpolate orbit position and velocity components to scan line time
    double x = (time[i] - torb[ind]) / dt;
    double x2 = x * x;
    double x3 = x2 * x;
    for (size_t j = 0; j < 3; j++) {
      posi[i][j] = a0[j] + a1[j] * x + a2[j] * x2 + a3[j] * x3;
      veli[i][j] = (a1[j] + 2 * a2[j] * x + 3 * a3[j] * x2) / dt;
    }
  }  // i-loop

  return 0;
}

int q_interp(size_t n_att_rec, size_t sdim, double *tq, quat_array *q,
             double *time, quat_array *qi) {
  //  Purpose: Interpolate quaternions to scan line times
  //
  //

  size_t ind = 0;
  for (size_t i = 0; i < sdim; i++) {
    //  Find input attitude vectors bracketing scan
    for (size_t j = ind; j < n_att_rec; j++) {
      if (time[i] > tq[j] && time[i] <= tq[j + 1]) {
        ind = j;
        break;
      }
    }

    //  Set up quaternion interpolation
    double dt = tq[ind + 1] - tq[ind];
    double qin[4];
    qin[0] = -q[ind][0];
    qin[1] = -q[ind][1];
    qin[2] = -q[ind][2];
    qin[3] = q[ind][3];

    double e[3], qr[4];
    qprod(qin, q[ind + 1], qr);
    memcpy(e, qr, 3 * sizeof(double));
    double sto2 = sqrt(e[0] * e[0] + e[1] * e[1] + e[2] * e[2]);
    for (size_t j = 0; j < 3; j++) e[j] /= sto2;

    // Interpolate quaternion to scan times
    double x = (time[i] - tq[ind]) / dt;
    float qri[4], qp[4];
    for (size_t j = 0; j < 3; j++) qri[j] = e[j] * sto2 * x;
    qri[3] = 1.0;
    qprod(q[ind], qri, qp);
    memcpy(qi[i], qp, 4 * sizeof(float));
  }

  return 0;
}

int l_sun(size_t sdim, int32_t iyr, int32_t iday, double *sec,
          orb_array *sunr) {
  //  Computes unit Sun vector in geocentric rotating coordinates, using
  //  subprograms to compute inertial Sun vector and Greenwich hour angle

  //  Get unit Sun vector in geocentric inertial coordinates
  sun2000(sdim, iyr, iday, sec, sunr);

  //  Get Greenwich mean sidereal angle
  for (size_t i = 0; i < sdim; i++) {
    double day = iday + sec[i] / 86400;
    double gha;
    gha2000(iyr, day, gha);
    double ghar = gha / RADEG;

    //  Transform Sun vector into geocentric rotating frame
    float temp0 = sunr[i][0] * cos(ghar) + sunr[i][1] * sin(ghar);
    float temp1 = sunr[i][1] * cos(ghar) - sunr[i][0] * sin(ghar);
    sunr[i][0] = temp0;
    sunr[i][1] = temp1;
  }

  return 0;
}

int sun2000(size_t sdim, int32_t iyr, int32_t idy, double *sec,
            orb_array *sun) {
  //  This subroutine computes the Sun vector in geocentric inertial
  //  (equatorial) coordinates.  It uses the model referenced in The
  //  Astronomical Almanac for 1984, Section S (Supplement) and documented
  //  for the SeaWiFS Project in "Constants and Parameters for SeaWiFS
  //  Mission Operations", in TBD.  The accuracy of the Sun vector is
  //  approximately 0.1 arcminute.

  float xk = 0.0056932;  // Constant of aberration

  //   Compute floating point days since Jan 1.5, 2000
  //    Note that the Julian day starts at noon on the specified date
  int16_t iyear = iyr;
  int16_t iday = idy;

  for (size_t i = 0; i < sdim; i++) {
    double t =
        jday(iyear, 1, iday) - (double)2451545.0 + (sec[i] - 43200) / 86400;

    double xls, gs, xlm, omega;
    double dpsi, eps, epsm;
    //  Compute solar ephemeris parameters
    ephparms(t, xls, gs, xlm, omega);

    // Compute nutation corrections for this day
    nutate(t, xls, gs, xlm, omega, dpsi, eps, epsm);

    //  Compute planet mean anomalies
    //   Venus Mean Anomaly
    double g2 = 50.40828 + 1.60213022 * t;
    g2 = fmod(g2, (double)360);

    //   Mars Mean Anomaly
    double g4 = 19.38816 + 0.52402078 * t;
    g4 = fmod(g4, (double)360);

    //  Jupiter Mean Anomaly
    double g5 = 20.35116 + 0.08309121 * t;
    g5 = fmod(g5, (double)360);

    //  Compute solar distance (AU)
    double rs =
        1.00014 - 0.01671 * cos(gs / RADEG) - 0.00014 * cos(2. * gs / RADEG);

    //  Compute Geometric Solar Longitude
    double dls = (6893. - 4.6543463e-4 * t) * sin(gs / RADEG) +
                 72. * sin(2. * gs / RADEG) - 7. * cos((gs - g5) / RADEG) +
                 6. * sin((xlm - xls) / RADEG) +
                 5. * sin((4. * gs - 8. * g4 + 3. * g5) / RADEG) -
                 5. * cos((2. * gs - 2. * g2) / RADEG) -
                 4. * sin((gs - g2) / RADEG) +
                 4. * cos((4. * gs - 8. * g4 + 3. * g5) / RADEG) +
                 3. * sin((2. * gs - 2. * g2) / RADEG) - 3. * sin(g5 / RADEG) -
                 3. * sin((2. * gs - 2. * g5) / RADEG);

    double xlsg = xls + dls / 3600;

    //  Compute Apparent Solar Longitude// includes corrections for nutation
    //  in longitude and velocity aberration
    double xlsa = xlsg + dpsi - xk / rs;

    //   Compute unit Sun vector
    sun[i][0] = cos(xlsa / RADEG);
    sun[i][1] = sin(xlsa / RADEG) * cos(eps / RADEG);
    sun[i][2] = sin(xlsa / RADEG) * sin(eps / RADEG);
  }  // i-loop

  return 0;
}

int qtom(float q[4], double rm[3][3]) {
  // Convert quaternion to equivalent direction cosine matrix

  rm[0][0] = q[0] * q[0] - q[1] * q[1] - q[2] * q[2] + q[3] * q[3];
  rm[0][1] = 2 * (q[0] * q[1] + q[2] * q[3]);
  rm[0][2] = 2 * (q[0] * q[2] - q[1] * q[3]);
  rm[1][0] = 2 * (q[0] * q[1] - q[2] * q[3]);
  rm[1][1] = -q[0] * q[0] + q[1] * q[1] - q[2] * q[2] + q[3] * q[3];
  rm[1][2] = 2 * (q[1] * q[2] + q[0] * q[3]);
  rm[2][0] = 2 * (q[0] * q[2] + q[1] * q[3]);
  rm[2][1] = 2 * (q[1] * q[2] - q[0] * q[3]);
  rm[2][2] = -q[0] * q[0] - q[1] * q[1] + q[2] * q[2] + q[3] * q[3];

  return 0;
}

int scan_ell(float p[3], double sm[3][3], double coef[10]) {
  //  This program calculates the coefficients which
  //  represent the Earth scan track in the sensor frame.

  //  The reference ellipsoid uses an equatorial radius of 6378.137 km and
  //  a flattening factor of 1/298.257 (WGS 1984).

  //  Calling Arguments
  //
  //  Name              Type    I/O     Description
  //
  //  pos(3)    R*4      I      ECR Orbit Position Vector (km)
  //  smat(3,3) R*4      I      Sensor Orientation Matrix
  //  coef(10)  R*4      O      Scan path coefficients

  double re = 6378.137;
  double f = 1 / 298.257;
  double omf2 = (1 - f) * (1 - f);

  //  Compute constants for navigation model using Earth radius values
  double rd = 1 / omf2;

  //  Compute coefficients of intersection ellipse in scan plane
  coef[0] = 1 + (rd - 1) * sm[0][2] * sm[0][2];
  coef[1] = 1 + (rd - 1) * sm[1][2] * sm[1][2];
  coef[2] = 1 + (rd - 1) * sm[2][2] * sm[2][2];
  coef[3] = (rd - 1) * sm[0][2] * sm[1][2] * 2;
  coef[4] = (rd - 1) * sm[0][2] * sm[2][2] * 2;
  coef[5] = (rd - 1) * sm[1][2] * sm[2][2] * 2;
  coef[6] = (sm[0][0] * p[0] + sm[0][1] * p[1] + sm[0][2] * p[2] * rd) * 2;
  coef[7] = (sm[1][0] * p[0] + sm[1][1] * p[1] + sm[1][2] * p[2] * rd) * 2;
  coef[8] = (sm[2][0] * p[0] + sm[2][1] * p[1] + sm[2][2] * p[2] * rd) * 2;
  coef[9] = p[0] * p[0] + p[1] * p[1] + p[2] * p[2] * rd - re * re;

  return 0;
}

int uni_geonav(float pos[3], float vel[3], double smat[3][3], double coef[10],
               float sunr[3], orb_array *pview, size_t npix, float *xlat,
               float *xlon, short *solz, short *sola, short *senz, short *sena,
               short *range) {
  //  This subroutine performs navigation of a scanning sensor on the
  //  surface of an ellipsoid based on an input orbit position vector and
  //  spacecraft orientation matrix.  It uses a closed-form algorithm for
  //  determining points on the ellipsoidal surface which involves
  //  determining the intersection of the scan plan with the ellipsoid.
  //  The sensor view vectors in the sensor frame are passed in as a 3xN array.

  //  The reference ellipsoid is set according to the scan
  //  intersection coefficients in the calling sequence// an equatorial
  //  radius of 6378.137 km. and a flattening factor of 1/298.257 are
  //  used by both the Geodetic Reference System (GRS) 1980 and the
  //  World Geodetic System (WGS) 1984.

  //  It then computes geometric parameters using the pixel locations on
  //  the Earth, the spacecraft position vector and the unit Sun vector in
  //  the geocentric rotating reference frame.  The outputs are arrays of
  //  geodetic latitude and longitude, solar zenith and azimuth and sensor
  //  zenith and azimuth.  The azimuth angles are measured from local
  //  North toward East.  Flag values of 999. are returned for any pixels
  //  whose scan angle is past the Earth's horizon.

  //  Reference: "Exact closed-form geolocation geolocation algorithm for
  //  Earth survey sensors", F. S. Patt and W. W. Gregg, IJRS, Vol. 15
  //  No. 18, 1994.

  //  Calling Arguments

  //  Name      Type    I/O     Description
  //
  //  pos(3)    R*4      I      ECR Orbit Position Vector (km) at scan
  //                                 mid-time
  //  vel(3)      R*4      I      ECR Orbit Velocity Vector (km/sec)
  //  smat(3,3) R*4      I      Sensor Orientation Matrix
  //  coef(10)  R*4      I      Scan path coefficients
  //  sunr(3)   R*4      I      Sun unit vector in geocentric rotating
  //                             reference frame
  //  pview(3,*)  R*4      I      Array of sensor view vectors
  //  npix        R*4      I      Number of pixels to geolocate
  //  xlat(*)   R*4      O      Pixel geodetic latitudes
  //  xlon(*)   R*4      O      Pixel geodetic longitudes
  //  solz(*)   I*2      O      Pixel solar zenith angles
  //  sola(*)   I*2      O      Pixel solar azimuth angles
  //  senz(*)   I*2      O      Pixel sensor zenith angles
  //  sena(*)   I*2      O      Pixel sensor azimuth angles
  //

  //    Program written by:     Frederick S. Patt
  //                            General Sciences Corporation
  //                            October 20, 1992
  //
  //    Modification History:

  //       Created universal version of the SeaWiFS geolocation algorithm
  //       by specifying view vectors as an input.  F. S. Patt, SAIC, 11/17/09

  // Earth ellipsoid parameters
  float f = 1 / 298.257;
  float omf2 = (1 - f) * (1 - f);
  gsl_vector *C = gsl_vector_alloc(3);
  
  for (size_t i = 0; i < npix; i++) {
    // Compute sensor-to-surface vectors for all scan angles
    // Compute terms for quadratic equation
    float o = coef[0] * pview[i][0] * pview[i][0] +
              coef[1] * pview[i][1] * pview[i][1] +
              coef[2] * pview[i][2] * pview[i][2] +
              coef[3] * pview[i][0] * pview[i][1] +
              coef[4] * pview[i][0] * pview[i][2] +
              coef[5] * pview[i][1] * pview[i][2];

    float p =
        coef[6] * pview[i][0] + coef[7] * pview[i][1] + coef[8] * pview[i][2];

    float q = coef[9];

    float r = p * p - 4 * q * o;

    xlat[i] = -999;
    xlon[i] = -999;

    solz[i] = -999;
    sola[i] = -999;
    senz[i] = -999;
    sena[i] = -999;
    range[i] = -999;

    //  Check for scan past edge of Earth
    if (r >= 0) {
      //  Solve for magnitude of sensor-to-pixel vector and compute components
      float d = (-p - sqrt(r)) / (2 * o);
      double x1[3];
      for (size_t j = 0; j < 3; j++) x1[j] = d * pview[i][j]; 

      //  Transform vector from sensor to geocentric frame
      
      gsl_matrix_view A = gsl_matrix_view_array((double *)smat, 3, 3);
      gsl_vector_view B = gsl_vector_view_array(x1, 3);

      gsl_blas_dgemv(CblasTrans, 1.0, &A.matrix, &B.vector, 0.0, C);

      float rh[3], geovec[3];
      double *ptr_C = gsl_vector_ptr(C, 0);
      for (size_t j = 0; j < 3; j++) {
        rh[j] = ptr_C[j];
        geovec[j] = pos[j] + rh[j];
      }  

      // Compute the local vertical, East and North unit vectors
      float uxy = geovec[0] * geovec[0] + geovec[1] * geovec[1];
      float temp = sqrt(geovec[2] * geovec[2] + omf2 * omf2 * uxy);

      float up[3];
      up[0] = omf2 * geovec[0] / temp;
      up[1] = omf2 * geovec[1] / temp;
      up[2] = geovec[2] / temp;
      float upxy = sqrt(up[0] * up[0] + up[1] * up[1]);

      float ea[3];
      ea[0] = -up[1] / upxy;
      ea[1] = up[0] / upxy;
      ea[2] = 0.0;

      // no = crossp(up,ea)
      float no[3];
      no[0] = -up[2] * ea[1];
      no[1] = up[2] * ea[0];
      no[2] = up[0] * ea[1] - up[1] * ea[0];

      //  Compute geodetic latitude and longitude
      xlat[i] = RADEG * asin(up[2]);	
      xlon[i] = RADEG * atan2(up[1], up[0]);
      range[i] = (short)((d - 400) * 10);

      // Transform the pixel-to-spacecraft and Sun vectors into local frame
      float rl[3], sl[3];
      rl[0] = -ea[0] * rh[0] - ea[1] * rh[1] - ea[2] * rh[2];
      rl[1] = -no[0] * rh[0] - no[1] * rh[1] - no[2] * rh[2];
      rl[2] = -up[0] * rh[0] - up[1] * rh[1] - up[2] * rh[2];

      sl[0] = sunr[0] * ea[0] + sunr[1] * ea[1] + sunr[2] * ea[2];
      sl[1] = sunr[0] * no[0] + sunr[1] * no[1] + sunr[2] * no[2];
      sl[2] = sunr[0] * up[0] + sunr[1] * up[1] + sunr[2] * up[2];

      //  Compute the solar zenith and azimuth
      solz[i] = (short)(100 * RADEG *
                        atan2(sqrt(sl[0] * sl[0] + sl[1] * sl[1]), sl[2]));

      // Check for zenith close to zero
      if (solz[i] > 0.01) sola[i] = (short)(100 * RADEG * atan2(sl[0], sl[1]));

      // Compute the sensor zenith and azimuth
      senz[i] = (short)(100 * RADEG *
                        atan2(sqrt(rl[0] * rl[0] + rl[1] * rl[1]), rl[2]));

      // Check for zenith close to zero
      if (senz[i] > 0.01) sena[i] = (short)(100 * RADEG * atan2(rl[0], rl[1]));

    }  // if on-earth
  }    // pixel loop

  gsl_vector_free(C);

  return 0;
}

int mat2rpy(float pos[3], float vel[3], double smat[3][3], float rpy[3], double om[3][3]) {
  //  This program calculates the attitude angles from the ECEF orbit vector and
  //  attitude matrix.  The rotation order is (1,2,3).

  //  The reference ellipsoid uses an equatorial radius of 6378.137 km and
  //  a flattening factor of 1/298.257 (WGS 1984).

  //  Calling Arguments

  //  Name                Type    I/O     Description
  //
  //  pos(3)      R*4      I      Orbit position vector (ECEF)
  //  vel(3)      R*4      I      Orbit velocity vector (ECEF)
  //  smat(3,3)   R*8      I      Sensor attitude matrix (ECEF to sensor)
  //  rpy(3)      R*4      O      Attitude angles (roll, pitch, yaw)
  //  om(3,3)	R*4	 O	ECEF-to-orbital frame matrix

  double rem = 6371;
  double f = 1 / (double)298.257;
  double omegae = 7.29211585494e-5;
  double omf2 = (1 - f) * (1 - f);

  // Determine local vertical reference axes
  double p[3], v[3];
  for (size_t j = 0; j < 3; j++) {
    p[j] = (double)pos[j];
    v[j] = (double)vel[j];
  }
  v[0] -= p[1] * omegae;
  v[1] += p[0] * omegae;

  //  Compute Z axis as local nadir vector
  double pm = sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
  double omf2p = (omf2 * rem + pm - rem) / pm;
  double pxy = p[0] * p[0] + p[1] * p[1];
  double temp = sqrt(p[2] * p[2] + omf2p * omf2p * pxy);

  double z[3];
  z[0] = -omf2p * p[0] / temp;
  z[1] = -omf2p * p[1] / temp;
  z[2] = -p[2] / temp;

  // Compute Y axis along negative orbit normal
  double on[3];
  on[0] = v[1] * z[2] - v[2] * z[1];
  on[1] = v[2] * z[0] - v[0] * z[2];
  on[2] = v[0] * z[1] - v[1] * z[0];

  double onm = sqrt(on[0] * on[0] + on[1] * on[1] + on[2] * on[2]);

  double y[3];
  for (size_t j = 0; j < 3; j++) y[j] = -on[j] / onm;

  // Compute X axis to complete orthonormal triad (velocity direction)
  double x[3];
  x[0] = y[1] * z[2] - y[2] * z[1];
  x[1] = y[2] * z[0] - y[0] * z[2];
  x[2] = y[0] * z[1] - y[1] * z[0];
    
  // Store local vertical reference vectors in matrix
  // double om[3][3];
  memcpy(&om[0][0], &x, 3 * sizeof(double));
  memcpy(&om[1][0], &y, 3 * sizeof(double));
  memcpy(&om[2][0], &z, 3 * sizeof(double));

  // Compute orbital-to-spacecraft matrix
  double rm[3][3];
  gsl_matrix_view rm_mat = gsl_matrix_view_array(&rm[0][0], 3, 3);

  int s;

  gsl_permutation *perm = gsl_permutation_alloc(3);

  // Compute the LU decomposition of this matrix

  double B_mat[3][3];
  memcpy(&B_mat[0][0], &om[0][0], 9 * sizeof(double));
  gsl_matrix_view B = gsl_matrix_view_array(&B_mat[0][0], 3, 3);

  gsl_linalg_LU_decomp(&B.matrix, perm, &s);

  // Compute the  inverse of the LU decomposition
  double inv[3][3];
  gsl_matrix_view inv_mat = gsl_matrix_view_array(&inv[0][0], 3, 3);

  gsl_linalg_LU_invert(&B.matrix, perm, &inv_mat.matrix);

  gsl_matrix_view A = gsl_matrix_view_array(&smat[0][0], 3, 3);

  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &A.matrix, &inv_mat.matrix,
                 0.0, &rm_mat.matrix);

  gsl_permutation_free(perm);

  // Compute attitude angles
  rpy[0] = RADEG * atan(-rm[2][1] / rm[2][2]);
  // double cosp = sqrt(rm[2][1] * rm[2][1] + rm[2][2] * rm[2][2]);
  // if (rm[2][2] < 0) cosp *= -1;
  // rpy[1] = RADEG * atan2(rm[2][0], cosp);
  rpy[1] = RADEG * asin(rm[2][0]);
  rpy[2] = RADEG * atan(-rm[1][0] / rm[0][0]);

  return 0;
}
