#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix_double.h>
#include <boost/algorithm/string/trim.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include "geolocate_oci.h"

#include <netcdf>

#include "nc4utils.h"
#include "timeutils.h"
#include <clo.h>
#include "global_attrs.h"
#include <terrain.h>
#include <genutils.h>
#define ERROR_EXIT(dim) {if(dim==0){printf("--Error--: the dimension %s is zero. Exiting. See line %d in file %s",#dim,__LINE__, __FILE__);exit(EXIT_FAILURE);} }

//    Modification history:
//  Programmer     Organization   Date     Ver   Description of change
//  ----------     ------------   ----     ---   ---------------------
//  Joel Gales     SAIC           09/13/21 0.01  Original development
//  Joel Gales     SAIC           06/14/22 0.03  Add support for tilt change
//  Gwyn Fireman   SAIC           07/10/23 0.07  Read global metadata from json file
//  Joel Gales     SAIC           08/08/23 0.08  Check scan times
//  Joel Gales     SAIC           08/18/23 0.80  Add earth_sun distance
//                                               correction attribute, support
//                                               for reflectance output
//  Joel Gales     SAIC           08/29/23 0.81  Change fill values
//  Joel Gales     SAIC           10/03/23 0.85  Change tilt output name to
//                                               tilt_angle
//                                               Change to callable routine
//                                               rather than standalone program

using namespace std;

/**
 * @brief Read a predicted PACE ephemeris file and convert vectors to ECR.
 * @param ephfile Name of the ephemeris file. Expected to be non-null and not empty.
 * @param otime Out-parameter indicating time of each record.
 * @param posj In/out-parameter indicating position of each record.
 * @param velj In/out-parameter indicating velocity of each record.
 * @param l1aEpoch Unix time of the start of the L1A file's first day
*/
void read_def_eph(string ephfile, vector<double>& otime, 
    vector<vector<double>>& posj, vector<vector<double>>& velj, double l1aEpoch) {

  ifstream inStream(ephfile, ios::in);
  if (!inStream.is_open()) {
    cerr << "-E- Unable to open definitive ephemeris file\n";
    exit(EXIT_FAILURE);
  }

  string tmpBuf;

  //size_t numRecords = 1921; // Number of records per file (32 hours)
  vector<double> posVec(3);
  vector<double> velVec(3);
  string ephStr;
  vector<string> tokens; // Datetime, 3d position, 3d velocity.
  vector<double> secondsVec(1921);
  double oldTime = 0;
  
  while (inStream.peek() != EOF) {
    getline(inStream, ephStr, '\n');

    istringstream lineStream(ephStr);
    string token; // Each individual token in this line

    tokens.clear();
    while (lineStream >> token) {
      tokens.push_back(token);
    }

    if(tokens.size() != 7 || tokens[0].find("T") == string::npos) {
      continue;
    }

    // skip duplicate records
    double theTime = isodate2unix(tokens[0].c_str());
    if(oldTime == theTime)
      continue;
    oldTime = theTime;

    posVec.clear();
    posVec.push_back(stod(tokens[1]));
    posVec.push_back(stod(tokens[2]));
    posVec.push_back(stod(tokens[3]));

    velVec.clear();
    velVec.push_back(stod(tokens[4]));
    velVec.push_back(stod(tokens[5]));
    velVec.push_back(stod(tokens[6]));

    posj.push_back(posVec);
    velj.push_back(velVec);
    otime.push_back(theTime - l1aEpoch);

    ephStr.clear();
    lineStream.clear();
  }
  inStream.close();

  if (otime.size() == 0) {
    cerr << "-E- No records found in definitive ephemeris file\n";
    exit(EXIT_FAILURE);
  }

}

int geolocate_oci( string l1a_filename, string geo_lut_filename,
                   geo_struct& geoLUT,
                   string l1b_filename, string dem_file, bool radiance,
                   string doi, const string ephFile, string pversion) { // add ephemeris file

  ofstream tempOut;

  // Open and read data from L1Afile
  NcFile *l1afile = new NcFile( l1a_filename, NcFile::read);
   
  NcGroup ncGrps[6];

  ncGrps[0] = l1afile->getGroup( "scan_line_attributes");
  ncGrps[1] = l1afile->getGroup( "spatial_spectral_modes");
  ncGrps[2] = l1afile->getGroup( "engineering_data");
  ncGrps[3] = l1afile->getGroup( "navigation_data");
  ncGrps[4] = l1afile->getGroup( "onboard_calibration_data");
  ncGrps[5] = l1afile->getGroup( "science_data");

  NcGroupAtt att;
  NcVar var;
  
  // Append call sequence to existing history
  string history = get_history(l1afile);
  //  history.append(call_sequence(argc, argv));

  // Get date (change this when year and day are added to time field)
  string tstart, tend;
  att = l1afile->getAtt("time_coverage_start");
  att.getValues(tstart);

  att = l1afile->getAtt("time_coverage_end");
  att.getValues(tend);

  size_t pos = tstart.find("T");
  if(pos == string::npos) {
    cout << "-E- L1A file time_coverage_start, bad format\n";
    exit(EXIT_FAILURE);
  }
  string granule_date = tstart.substr(0,pos);

  // Vector of the positions of OCI for 32 hours starting at tstart with a resolution of 1 minute
  vector<vector<double>> ephPosVec;
  // Vector of the velocities of OCI for 32 hours starting at tstart with a resolution of 1 minute
  vector<vector<double>> ephVelVec;
  vector<double> ephOTime;

  if (!ephFile.empty()) {
    double l1a_start = isodate2unix(tstart.c_str());
    double l1a_end = isodate2unix(tend.c_str());

    int16_t yr,mo,dy;
    double secs;
    unix2ymds(l1a_start, &yr, &mo, &dy, &secs);
    double l1a_epoch = ymds2unix(yr, mo, dy, 0.0);

    read_def_eph(ephFile, ephOTime, ephPosVec, ephVelVec, l1a_epoch);

    // Ensure ephem file covers L1A file time
    if((l1a_start < (ephOTime.front() + l1a_epoch)) ||
        (l1a_end > (ephOTime.back() + l1a_epoch))) {
      cerr << "-E- Definitive ephemeris file does not cover L1A file time" << endl;
      exit(EXIT_FAILURE);
    }
  }

  // Scan time, spin ID and HAM side
  NcDim nscan_dim = l1afile->getDim("number_of_scans");
  uint32_t nscan = nscan_dim.getSize();

  // Number of swir bands

  uint32_t nswband = l1afile->getDim("SWIR_bands").getSize();

  NcDim mcescan_dim = l1afile->getDim("number_of_mce_scans");
  uint32_t nmcescan = mcescan_dim.getSize();

  double *sstime = new double[nscan];
  var = ncGrps[0].getVar( "scan_start_time");
  var.getVar( sstime);

  int32_t *spin = new int32_t[nscan];
  var = ncGrps[0].getVar( "spin_ID");
  var.getVar( spin);

  uint8_t *hside = new uint8_t[nscan];
  var = ncGrps[0].getVar( "HAM_side");
  var.getVar( hside);
  
  uint32_t sdim=0;
  for (size_t i=0; i<nscan; i++) {
    if ( spin[i] > 0) {
      sstime[sdim] = sstime[i];
      hside[sdim] = hside[i];
      sdim++;
    }
  }

  // Check for and fill in missing scan times
  short *sfl = new short[sdim]();
  check_scan_times( sdim, sstime, sfl);
  
  NcDim natt_dim = l1afile->getDim("att_records");
  uint32_t n_att_rec = natt_dim.getSize();
  ERROR_EXIT(n_att_rec);

  double *att_time = new double[n_att_rec];
  var = ncGrps[3].getVar( "att_time");
  var.getVar( att_time);

  NcDim nquat_dim = l1afile->getDim("quaternion_elements");
  uint32_t n_quat_elems = nquat_dim.getSize();
  ERROR_EXIT(n_quat_elems);

  float **att_quat = new float *[n_att_rec];
  att_quat[0] = new float[n_quat_elems*n_att_rec];
  for (size_t i=1; i<n_att_rec; i++)
    att_quat[i] = att_quat[i-1] + n_quat_elems;

  var = ncGrps[3].getVar( "att_quat");
  var.getVar( &att_quat[0][0]);

  NcDim norb_dim = l1afile->getDim("orb_records");
  uint32_t n_orb_rec = norb_dim.getSize();
  ERROR_EXIT(n_orb_rec);

  vector<double> orb_time(n_orb_rec);
  if (ephFile.empty()) {
    var = ncGrps[3].getVar( "orb_time");
    var.getVar( orb_time.data());
  }

  float **orb_pos = new float *[n_orb_rec];
  orb_pos[0] = new float[3*n_orb_rec];
  for (size_t i=1; i<n_orb_rec; i++) orb_pos[i] = orb_pos[i-1] + 3;
  var = ncGrps[3].getVar( "orb_pos");
  var.getVar( &orb_pos[0][0]);

  float **orb_vel = new float *[n_orb_rec];
  orb_vel[0] = new float[3*n_orb_rec];
  for (size_t i=1; i<n_orb_rec; i++) orb_vel[i] = orb_vel[i-1] + 3;
  var = ncGrps[3].getVar( "orb_vel");
  var.getVar( &orb_vel[0][0]);

  NcDim ntilt_dim = l1afile->getDim("tilt_samples");
  uint32_t n_tilts = ntilt_dim.getSize();
  ERROR_EXIT(n_tilts);

  float *tiltin = new float[n_tilts];
  var = ncGrps[3].getVar( "tilt");
  var.getVar( tiltin);

  double *ttime = new double[n_tilts];
  var = ncGrps[3].getVar( "tilt_time");
  var.getVar( ttime);


  // **************************** //
  // *** Read geolocation LUT *** //
  // **************************** //

  NcFile *geoLUTfile = new NcFile( geo_lut_filename, NcFile::read);
  //  geo_struct geoLUT;

  NcGroup gidTime, gidCT, gidRTA_HAM, gidPlanarity;
  
  gidTime = geoLUTfile->getGroup( "time_params");
  var = gidTime.getVar( "master_clock");
  var.getVar( &geoLUT.master_clock);
  var = gidTime.getVar( "MCE_clock");
  var.getVar( &geoLUT.MCE_clock);

  gidCT = geoLUTfile->getGroup( "coord_trans");
  var = gidCT.getVar( "sc_to_tilt");
  var.getVar( &geoLUT.sc_to_tilt);
  var = gidCT.getVar( "tilt_axis");
  var.getVar( &geoLUT.tilt_axis);
  var = gidCT.getVar( "tilt_angles");
  var.getVar( &geoLUT.tilt_angles);
  var = gidCT.getVar( "tilt_home");
  var.getVar( &geoLUT.tilt_home);
  var = gidCT.getVar( "tilt_to_oci_mech");
  var.getVar( &geoLUT.tilt_to_oci_mech);
  var = gidCT.getVar( "oci_mech_to_oci_opt");
  var.getVar( &geoLUT.oci_mech_to_oci_opt);

  gidRTA_HAM = geoLUTfile->getGroup( "RTA_HAM_params");
  var = gidRTA_HAM.getVar( "RTA_axis");
  var.getVar( &geoLUT.rta_axis);
  var = gidRTA_HAM.getVar( "HAM_axis");
  var.getVar( &geoLUT.ham_axis);
  var = gidRTA_HAM.getVar( "HAM_AT_angles");
  var.getVar( &geoLUT.ham_at_angles);
  var = gidRTA_HAM.getVar( "HAM_CT_angles");
  var.getVar( &geoLUT.ham_ct_angles);
  var = gidRTA_HAM.getVar( "RTA_enc_scale");
  var.getVar( &geoLUT.rta_enc_scale);
  var = gidRTA_HAM.getVar( "HAM_enc_scale");
  var.getVar( &geoLUT.ham_enc_scale);
  var = gidRTA_HAM.getVar( "RTA_nadir");
  var.getVar( geoLUT.rta_nadir);

  gidPlanarity = geoLUTfile->getGroup( "planarity");
  var = gidPlanarity.getVar( "along_scan_planarity");
  var.getVar( &geoLUT.as_planarity);
  var = gidPlanarity.getVar( "along_track_planarity");
  var.getVar( &geoLUT.at_planarity);
  
  geoLUTfile->close();

  // MCE telemetry
  int32_t ppr_off;
  double revpsec, secpline;

  int32_t *mspin = new int32_t[nmcescan];
  int32_t *ot_10us = new int32_t[nmcescan];
  uint8_t *enc_count = new uint8_t[nmcescan];
  int16_t board_id;
  
  NcDim nenc_dim = l1afile->getDim("encoder_samples");
  uint32_t nenc = nenc_dim.getSize();

  float **hamenc = new float *[nmcescan];
  hamenc[0] = new float[nenc*nmcescan];
  for (size_t i=1; i<nmcescan; i++) hamenc[i] = hamenc[i-1] + nenc;

  float **rtaenc = new float *[nmcescan];
  rtaenc[0] = new float[nenc*nmcescan];
  for (size_t i=1; i<nmcescan; i++) rtaenc[i] = rtaenc[i-1] + nenc;

  int16_t iret;
  read_mce_tlm( l1afile, geoLUT, ncGrps[2], nmcescan, nenc, ppr_off, revpsec,
                secpline, board_id, mspin, ot_10us, enc_count, hamenc, rtaenc,
                iret);

  int32_t iyr, iday, msec;
  isodate2ydmsec((char *) tstart.c_str(), &iyr, &iday, &msec);

  double ecmat[3][3];
  quat_array *quatr = new quat_array[n_att_rec];
  orb_array *posr; // Array of 1x3 vectors
  orb_array *velr; // Array of 1x3 vectors

  // Transform orbit (if needed) and attitude from J2000 to ECR
  if (!ephFile.empty()) {

    size_t ephOTimeLBIndex = -1;
    size_t ephOTimeUBIndex = -1;

    for (size_t i = 0; i < ephOTime.size(); i++) {
      double value = ephOTime[i];

      if (value <= sstime[0]) {
        ephOTimeLBIndex = i;
      }
      if (value >= sstime[nscan - 1]) {
        ephOTimeUBIndex = i;
        break;
      }
    }


    double omega_e = 7.29211585494e-5;
    size_t num_eph_recs = ephOTimeUBIndex - ephOTimeLBIndex + 1;

    n_orb_rec = num_eph_recs;
    posr = new orb_array[n_orb_rec]; // Array of 1x3 vectors
    velr = new orb_array[n_orb_rec]; // Array of 1x3 vectors

    bzero(posr, n_orb_rec * sizeof(orb_array));
    bzero(velr, n_orb_rec * sizeof(orb_array));
    
    orb_time.resize(n_orb_rec);

    for (size_t idxRec = 0; idxRec < num_eph_recs; idxRec++) {
      size_t ephIndex = idxRec + ephOTimeLBIndex;
      orb_time[idxRec] = ephOTime[ephIndex];

      j2000_to_ecr(iyr, iday, ephOTime[ephIndex], ecmat);
      
      for (size_t idxLine = 0; idxLine < 3; idxLine++) {
        for (size_t idxValue = 0; idxValue < 3; idxValue++) {
          posr[idxRec][idxLine] += ecmat[idxLine][idxValue] * ephPosVec[ephIndex][idxValue];
          velr[idxRec][idxLine] += ecmat[idxLine][idxValue] * ephVelVec[ephIndex][idxValue];
        }
      }

      velr[idxRec][0] = velr[idxRec][0] + posr[idxRec][1] * omega_e;
      velr[idxRec][1] = velr[idxRec][1] - posr[idxRec][0] * omega_e;
    }
  } else {
    posr = new orb_array[n_orb_rec]; // Array of 1x3 vectors
    velr = new orb_array[n_orb_rec]; // Array of 1x3 vectors
    // Orbit
    for (size_t i = 0; i < n_orb_rec; i++) {
      //  j2000_to_ecr(iyr, iday, orb_time[i], ecmat);

      for (size_t j = 0; j < 3; j++) {
        posr[i][j] = orb_pos[i][j] / 1000;
        velr[i][j] = orb_vel[i][j] / 1000;
      }
    }  // i loop
  }


  // Attitude
  for (size_t i = 0; i < n_att_rec; i++) {
    double ecq[4], qt2[4];
    float qt1[4];
    j2000_to_ecr(iyr, iday, att_time[i], ecmat);
    mtoq(ecmat, ecq);

    memcpy(qt1, &att_quat[i][0], 3 * sizeof(float));
    qt1[3] = att_quat[i][3];

    qprod(ecq, qt1, qt2);
    for (size_t j = 0; j < 4; j++) quatr[i][j] = qt2[j];
  }  // i loop

  // ******************************************** //
  // *** Get spatial and spectral aggregation *** //
  // ******************************************** //
  NcDim spatzn_dim = l1afile->getDim("spatial_zones");
  uint32_t spatzn = spatzn_dim.getSize();

  int16_t *dtype = new int16_t[spatzn];
  var = ncGrps[1].getVar( "spatial_zone_data_type");
  var.getVar( dtype);

  int16_t *lines = new int16_t[spatzn];
  var = ncGrps[1].getVar( "spatial_zone_lines");
  var.getVar( lines);

  int16_t *iagg = new int16_t[spatzn];
  var = ncGrps[1].getVar( "spatial_aggregation");
  var.getVar( iagg);

  float clines[32400], slines[4050];
  uint16_t pcdim, psdim;
  double ev_toff, deltc[32400], delts[4050];
  bool dark = false;
  get_ev( secpline, dtype, lines, iagg, pcdim, psdim, ev_toff,
          clines, slines, deltc, delts, dark, iret);
  if ( iret < 0) {
    cout << "No Earth view in file: " << l1a_filename << endl;
    l1afile->close();
    return 1;
  }

  double *evtime = new double[sdim];
  for (size_t i = 0; i < sdim; i++)
    if (sstime[i] == -32767) evtime[i] = sstime[i];
    else
      evtime[i] = sstime[i] + ev_toff;

  // Interpolate orbit, attitude and tilt to scan times
  orb_array *posi = new orb_array[sdim]();
  orb_array *veli = new orb_array[sdim]();
  orb_array *rpy = new orb_array[sdim]();
  orb_interp(n_orb_rec, sdim, orb_time.data(), posr, velr, evtime, posi, veli);

  quat_array *quati = new quat_array[sdim]();
  q_interp(n_att_rec, sdim, att_time, quatr, evtime, quati);

  float *tilt = new float[sdim];
  tilt_interp(n_tilts, sdim, ttime, tiltin, evtime, tilt);
  for (size_t i = 0; i < sdim; i++)
       sfl[i] |= 1;  //  set tilt flag for all scans
  for (size_t i = 0; i < sdim; i++) {
      tilt[i] += geoLUT.tilt_home; // Add tilt home position to angles
      tilt[i] = tilt[i] < geoLUT.tilt_angles[0] ? geoLUT.tilt_angles[0] : tilt[i];
      tilt[i] = tilt[i] > geoLUT.tilt_angles[1] ? geoLUT.tilt_angles[1] : tilt[i];
      // Fred's tilt detection 
      // for any 2 adjecent and = tilts clear tilt flag (ie, bit AND with 110b)
      if( ( i > 0 ) && ( tilt[i-1] == tilt[i] ) ) {
          sfl[i-1] &= 6;
          sfl[i] &= 6;
      }
  }
  float *xlon = new float[sdim * pcdim]();
  float *xlat = new float[sdim * pcdim]();
  short *solz = new short[sdim * pcdim]();
  short *sola = new short[sdim * pcdim]();
  short *senz = new short[sdim * pcdim]();
  short *sena = new short[sdim * pcdim]();
  
  uint8_t *qfl = new uint8_t[sdim * pcdim]();
  // short *range = new short[sdim * psdim]();
  short *height = new short[sdim * pcdim]();

  // Initialize output arrays
  for (size_t i = 0; i < sdim * pcdim; i++) {
    xlon[i] = -32767;
    xlat[i] = -32767;

    solz[i] = -32767;
    sola[i] = -32767;
    senz[i] = -32767;
    sena[i] = -32767;
  }

  // Generate pointing vector array
  float **pview = new float *[pcdim];
  pview[0] = new float[3*pcdim];
  for (size_t i=1; i<pcdim; i++) pview[i] = pview[i-1] + 3;

  // Get Sun vectors
  orb_array *sunr = new orb_array[sdim]();
  
  double *rs = new double[sdim];
  l_sun(sdim, iyr, iday, evtime, sunr, rs);
  double distcorr = pow(rs[sdim/2],2);

  // S/C matrices
  gsl_matrix *sc2tiltp = gsl_matrix_alloc(3, 3);
  gsl_matrix *sc2ocim = gsl_matrix_alloc(3, 3);
  gsl_matrix *sc_to_oci = gsl_matrix_alloc(3, 3);
  gsl_matrix *smat = gsl_matrix_alloc(3, 3);

  double *thetap = new double[pcdim]();
  
  //////////////////////////////
  // Geolocate each scan line //
  //////////////////////////////
  float lonwest = +180;
  float loneast = -180;
  float latsouth = +90;
  float latnorth = -90;
  for (size_t iscn = 0; iscn < sdim; iscn++) {
    
    if (sstime[iscn] != -999.0) {

      gsl_matrix_view A;
      gsl_matrix_view B;
      double *ptr_C;
      
      // Model tilt rotation using a quaternion
      float qt[4];
      qt[0] = geoLUT.tilt_axis[0]*sin(tilt[iscn]/2/RADEG);
      qt[1] = geoLUT.tilt_axis[1]*sin(tilt[iscn]/2/RADEG);
      qt[2] = geoLUT.tilt_axis[2]*sin(tilt[iscn]/2/RADEG);
      qt[3] = cos(tilt[iscn]/2/RADEG);

      double tiltm[3][3];
      qtom(qt, tiltm);

      // Combine tilt and alignments

      // sc2tiltp = tiltm#geo_lut.sc_to_tilt
      A = gsl_matrix_view_array( (double *) geoLUT.sc_to_tilt, 3, 3);
      B = gsl_matrix_view_array( (double *) tiltm, 3, 3);
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0,
                     &A.matrix, &B.matrix, 0.0, sc2tiltp);
      ptr_C = gsl_matrix_ptr(sc2tiltp, 0, 0);

      // sc2ocim = geo_lut.tilt_to_oci_mech#sc2tiltp
      B = gsl_matrix_view_array( (double *)  geoLUT.tilt_to_oci_mech, 3, 3);
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0,
                     sc2tiltp, &B.matrix, 0.0, sc2ocim);
      ptr_C = gsl_matrix_ptr(sc2ocim, 0, 0);

      // sc_to_oci = geo_lut.oci_mech_to_oci_opt#sc2ocim
      B = gsl_matrix_view_array( (double *)  geoLUT.oci_mech_to_oci_opt, 3, 3);
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0,
                     sc2ocim, &B.matrix, 0.0, sc_to_oci);
      ptr_C = gsl_matrix_ptr(sc_to_oci, 0, 0);

      
      // Convert quaternion to matrix
      double qmat[3][3];
      qtom(quati[iscn], qmat);

      // smat = sc_to_oci#qmat
      B = gsl_matrix_view_array(&qmat[0][0], 3, 3);
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0,
                     sc_to_oci, &B.matrix, 0.0, smat);

      // Compute attitude angles (informational only)
      mat2rpy(posi[iscn], veli[iscn], qmat, rpy[iscn]);

      // Get scan ellipse coefficients
      ptr_C = gsl_matrix_ptr(smat, 0, 0);
      double coef[10];
      scan_ell(posi[iscn], (double(*)[3])ptr_C, coef);

      // Generate pointing vector and relative time arrays in instrument frame
      get_oci_vecs( nscan, pcdim, geoLUT.as_planarity, geoLUT.at_planarity,
                    geoLUT.rta_nadir, geoLUT.ham_ct_angles,
                    ev_toff, spin[iscn], hside[iscn], clines, deltc,
                    revpsec, ppr_off, board_id,  nmcescan, mspin, enc_count,
                    &hamenc[0], &rtaenc[0], pview, thetap, iret);

      sfl[iscn] |= 4*iret;

      // Geolocate pixels
      size_t ip = iscn * pcdim;
      oci_geonav(dem_file.c_str(), posi[iscn], veli[iscn], (double(*)[3])ptr_C,
                 coef, sunr[iscn], pview, pcdim, delts, &xlat[ip], &xlon[ip],
                 &solz[ip], &sola[ip], &senz[ip], &sena[ip], &height[ip],
                 lonwest, loneast, latsouth, latnorth);
    } else {
      qfl[iscn] |= 4;
    }
  }  // scan loop
  gsl_matrix_free(smat);
  gsl_matrix_free(sc2tiltp);
  gsl_matrix_free(sc2ocim);
  gsl_matrix_free(sc_to_oci);

  int geo_ncid;
  int geo_gid[10];

  // Get number of bands
  NcDim ntaps_dim = l1afile->getDim("number_of_taps");
  uint32_t ntaps = ntaps_dim.getSize();

  int16_t *bagg = new int16_t[ntaps];
  var = ncGrps[1].getVar( "blue_spectral_mode");
  var.getVar( bagg);

  int16_t *ragg = new int16_t[ntaps];
  var = ncGrps[1].getVar( "red_spectral_mode");
  var.getVar( ragg);

  size_t *ia = new size_t[ntaps];
  uint32_t ntb[16];
  uint32_t bib, rib;
  uint32_t bbb, rbb;
  get_nib_nbb( ntaps, ia, ntb, bagg, bib, bbb);
  get_nib_nbb( ntaps, ia, ntb, ragg, rib, rbb);
  
  static geoFile outfile;
  
  outfile.createFile(l1b_filename.c_str(),
                     //  "$OCDATAROOT/oci/OCI_GEO_Data_Structure.cdl",
                     "$OCDATAROOT/oci/OCI_GEO_Level-1B_Data_Structure.cdl",
                     sdim, &geo_ncid, geo_gid, bbb, rbb,
                     pcdim, psdim, nswband, geoLUT.rta_nadir, radiance);

  string varname;

  /*
  char buf[32];
  strcpy(buf, unix2isodate(now(), 'G'));
  nc_put_att_text(geo_ncid, NC_GLOBAL, "date_created", strlen(buf), buf);

  varname.assign("scan_time");
  status = nc_inq_varid(geo_gid[scan_attr_geo], varname.c_str(), &varid);
  status = nc_put_var_double(geo_gid[scan_attr_geo], varid, stime);
  check_err(status, __LINE__, __FILE__);

  status = nc_put_att_text(geo_ncid, NC_GLOBAL,
                           "time_coverage_start", strlen(tstart), tstart);

  status = nc_put_att_text(geo_ncid, NC_GLOBAL,
                           "time_coverage_end", strlen(tend), tend);
  */

  
  vector<size_t> start, count;
  
  start.clear();
  start.push_back(0);
  start.push_back(0);
  start.push_back(0);

  count.push_back(sdim);

  // Need to increment ncGrps index by 1 for new geo_l1b cdl file
  
  var = outfile.ncGrps[1].getVar( "time");
  var.putVar( start, count, evtime);
  var.putAtt("units", "seconds since " + granule_date);

  var = outfile.ncGrps[1].getVar( "HAM_side");
  var.putVar( start, count, hside);

  var = outfile.ncGrps[1].getVar( "scan_quality_flags");
  var.putVar( start, count, sfl);

  var = outfile.ncGrps[3].getVar( "tilt_angle");
  var.putVar( start, count, tilt);
  
  count.push_back(4);

  var = outfile.ncGrps[3].getVar( "att_quat");
  var.putVar( start, count, quati);

  count.pop_back();
  count.push_back(3);
  
  var = outfile.ncGrps[3].getVar( "att_ang");
  var.putVar( start, count, rpy);

  var = outfile.ncGrps[3].getVar( "orb_pos");
  var.putVar( start, count, posi);

  var = outfile.ncGrps[3].getVar( "orb_vel");
  var.putVar( start, count, veli);

  var = outfile.ncGrps[3].getVar( "sun_ref");
  var.putVar( start, count, sunr);

  count.pop_back();
  count.push_back(pcdim);
  
  var = outfile.ncGrps[2].getVar( "latitude");
  var.putVar( start, count, xlat);

  var = outfile.ncGrps[2].getVar( "longitude");
  var.putVar( start, count, xlon);

  var = outfile.ncGrps[2].getVar( "quality_flag");
  var.putVar( start, count, qfl);

  var = outfile.ncGrps[2].getVar( "sensor_azimuth");
  var.putVar( start, count, sena);

  var = outfile.ncGrps[2].getVar( "sensor_zenith");
  var.putVar( start, count, senz);

  var = outfile.ncGrps[2].getVar( "solar_azimuth");
  var.putVar( start, count, sola);

  var = outfile.ncGrps[2].getVar( "solar_zenith");
  var.putVar( start, count, solz);

  //  var = outfile.ncGrps[2].getVar( "range");
  //var.putVar( start, count, range);

  var = outfile.ncGrps[2].getVar( "height");
  var.putVar( start, count, height);

  // write global attributes, including history and date_created
  set_global_attrs(outfile.geofile, history, doi, pversion);

  outfile.geofile->putAtt("earth_sun_distance_correction",
                          NC_DOUBLE, 1, &distcorr);

  outfile.geofile->putAtt( "cdm_data_type", "swath");

  outfile.geofile->putAtt("geospatial_lat_min", NC_FLOAT, 1, &latsouth);
  outfile.geofile->putAtt("geospatial_lat_max", NC_FLOAT, 1, &latnorth);
  outfile.geofile->putAtt("geospatial_lon_min", NC_FLOAT, 1, &lonwest);
  outfile.geofile->putAtt("geospatial_lon_max", NC_FLOAT, 1, &loneast);

  outfile.geofile->putAtt( "geospatial_bounds_crs", "EPSG:4326");
  /*
geospatial_bounds = "POLYGON ((85.129562 -78.832314,81.489693 49.646988...))"
   */
  outfile.close();

  delete[](att_time);
  // delete[](orb_time);
  delete[](quatr);
  delete[](posr);
  delete[](velr);
  delete[](posi);
  delete[](veli);
  delete[](rpy);
  delete[](quati);
  delete[](tilt);
  delete[](xlon);
  delete[](xlat);
  delete[](solz);
  delete[](sola);
  delete[](senz);
  delete[](sena);
  //  delete[](range);
  delete[](height);
  delete[](tiltin);
  delete[](ttime);
  delete[](qfl);
  delete[]pview[0]; delete[](pview);
  delete[](sunr);

  delete [] hamenc[0]; delete [] hamenc;
  delete [] rtaenc[0]; delete [] rtaenc;
   
  delete [](sstime);
  delete [](evtime);
  delete [](thetap);
  delete [](spin);
  delete [](mspin);
  delete [](ot_10us);
  delete []att_quat[0]; delete [](att_quat);
  delete []orb_vel[0]; delete [](orb_vel);
  delete []orb_pos[0]; delete[](orb_pos);
  delete [](sfl);
  delete [](hside);
  delete [](enc_count);
  delete [](ragg);
  delete [](bagg);
  delete [](iagg);
  delete [](dtype);
  delete [](lines);
  delete(l1afile);
  delete(geoLUTfile);
  //  clo_deleteList(optionList);
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

  // initialize the arrays
  for (int i=0; i<3; i++) {
    a0[i] = 0.0;
    a1[i] = 0.0;
    a2[i] = 0.0;
    a3[i] = 0.0;
  }

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

  // orbit time
  double startOrbTime = torb[0];
  double endOrbTime = torb[n_orb_rec-1];
  const double NAV_FILL = -9999999.0;

  size_t ind = 0;
  for (size_t i = 0; i < sdim; i++) {

    // scantime outside valid orbit time, use fill value
    if (time[i] < startOrbTime || time[i] > endOrbTime) {
      for (int vector=0; vector<3; vector++) {
        posi[i][vector] = NAV_FILL;
        veli[i][vector] = BAD_FLT;
      }
      // skip calulations below
      continue; 
    }
    
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
      if (dt >= 3.5) {
        a2[j] = 3 * p[ind + 1][j] - 3 * p[ind][j] - 2 * v[ind][j] * dt -
          v[ind + 1][j] * dt;
        a3[j] = 2 * p[ind][j] - 2 * p[ind + 1][j] + v[ind][j] * dt +
          v[ind + 1][j] * dt;
      } else {
        a2[j] = (v[ind + 1][j] - v[ind][j]) / 2;
        a3[j] = 0.0;
        // a2(*) = (v(*,ind+1)-v(*,ind))*dt/2
        //         a3(*) = 0.d0
      }
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

  // Attitude
  double startAttTime = tq[0];
  double endAttTime = tq[n_att_rec-1];
  const double NAV_FILL = -32767.0;
  size_t ind = 0;
  for (size_t i = 0; i < sdim; i++) {

    // if scan time is outside the valid att time, use fill value
    if (time[i] < startAttTime || time[i] > endAttTime) {
      for (int element=0; element<4; element++) {
        qi[i][element] = NAV_FILL;

      }
      continue;
    }

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


int tilt_interp(size_t n_tilts, size_t sdim, double *ttilt, float *tiltin,
                double *time, float *tilt) {

  double startTiltTime = ttilt[0];
  double endTiltTime = ttilt[n_tilts-1];
  const double NAV_FILL = -32767.0;

  size_t ind = 0;
  for (size_t i = 0; i < sdim; i++) {

    if (time[i] < startTiltTime || time[i] > endTiltTime) {
      tilt[i] = NAV_FILL;
      continue;
    }

    //  Find input tilt vectors bracketing scan
    for (size_t j = ind; j < n_tilts; j++) {
      if (time[i] > ttilt[j] && time[i] <= ttilt[j + 1]) {
        ind = j;
        break;
      }
    }

    double x = (time[i] - ttilt[ind])/(ttilt[ind+1] - ttilt[ind]);
    tilt[i] = (1-x)*tiltin[ind] + x*tiltin[ind+1];
  }
  
  return 0;
}


int l_sun(size_t sdim, int32_t iyr, int32_t iday, double *sec,
          orb_array *sunr, double *rs) {
  //  Computes unit Sun vector in geocentric rotating coordinates, using
  //  subprograms to compute inertial Sun vector and Greenwich hour angle

  //  Get unit Sun vector in geocentric inertial coordinates
  sun2000(sdim, iyr, iday, sec, sunr, rs);

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
            orb_array *sun, double *rs) {
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
    rs[i]
      = 1.00014 - 0.01671 * cos(gs / RADEG) - 0.00014 * cos(2. * gs / RADEG);

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
    double xlsa = xlsg + dpsi - xk / rs[i];

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

int oci_geonav(const char* dem_file, float pos[3], float vel[3],
               double smat[3][3], double coef[10], float sunr[3],
               float **pview, size_t npix, double *delt, float *xlat,
               float *xlon, short *solz, short *sola,
               short *senz, short *sena, short *height,
               float& lonwest, float& loneast, float& latsouth, float& latnorth)
               {
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
  //  dem_file  char*   I       Digital elevation map file
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
  //  height(*) I*2      O      Pixel terrain height
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
    double o =
      coef[0] * pview[i][0] * pview[i][0] +
      coef[1] * pview[i][1] * pview[i][1] +
      coef[2] * pview[i][2] * pview[i][2] +
      coef[3] * pview[i][0] * pview[i][1] +
      coef[4] * pview[i][0] * pview[i][2] +
      coef[5] * pview[i][1] * pview[i][2];

    double p =
      coef[6] * pview[i][0] + coef[7] * pview[i][1] + coef[8] * pview[i][2];

    double q = coef[9];

    double r = p * p - 4 * q * o;

    xlat[i] = -32767;
    xlon[i] = -32767;

    solz[i] = -32767;
    sola[i] = -32767;
    senz[i] = -32767;
    sena[i] = -32767;
    height[i] = -32767;
    //  Check for scan past edge of Earth
    if (r >= 0) {

      //  Solve for magnitude of sensor-to-pixel vector and compute components
      double d = (-p - sqrt(r)) / (2 * o);
      double x1[3];
      for (size_t j = 0; j < 3; j++) x1[j] = d * pview[i][j];

      //  Convert velocity vector to ground speed
      float re = 6378.137;
      //      double omegae = 7.29211585494e-05;
      double pm = sqrt(pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2]);
      double clatg = sqrt(pos[0]*pos[0] + pos[1]*pos[1]) / pm;
      double rg = re*(1.-f)/sqrt(1.-(2.-f)*f*clatg*clatg);
      double v[3];
      //      v[0] = (vel[0] - pos[1]*omegae) * rg / pm;
      //v[1] = (vel[1] - pos[0]*omegae) * rg / pm;
      v[0] = vel[0] * rg / pm;
      v[1] = vel[1] * rg / pm;
      v[2] = vel[2] * rg / pm;
      
      //  Transform vector from sensor to geocentric frame
      gsl_matrix_view A = gsl_matrix_view_array((double *)smat, 3, 3);
      gsl_vector_view B = gsl_vector_view_array(x1, 3);

      gsl_blas_dgemv(CblasTrans, 1.0, &A.matrix, &B.vector, 0.0, C);

      float rh[3], geovec[3];
      double *ptr_C = gsl_vector_ptr(C, 0);
      for (size_t j = 0; j < 3; j++) {
        rh[j] = ptr_C[j];
        geovec[j] = pos[j] + rh[j] + v[j]*delt[i];
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
      float xlat_ = RADEG * asin(up[2]);
      float xlon_ = RADEG * atan2(up[1], up[0]);
      //check if it is in the range
      if( std::abs(xlat_) > 90.0e0 || std::abs(xlon_) > 180.0e0){
          continue;
        }
      xlat[i] = xlat_;
      xlon[i] = xlon_;
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

      // terrain correct the pixel
      float tmp_sena = sena[i] / 100.0;
      float tmp_senz = senz[i] / 100.0;
      float tmp_height = 0;
      get_nc_height(dem_file, &xlon[i], &xlat[i],
                    &tmp_senz, &tmp_sena, &tmp_height);
      // there is also a bug in terraing correction
      if( std::abs(xlat[i]) > 90.0e0 || std::abs(xlon[i]) > 180.0e0){
        xlat[i] = -32767;
        xlon[i] = -32767;
        solz[i] = -32767;
        sola[i] = -32767;
        senz[i] = -32767;
        sena[i] = -32767;
        continue;
      }

      if (xlon[i] > loneast) loneast = xlon[i];
      if (xlon[i] < lonwest) lonwest = xlon[i];
      if (xlat[i] < latsouth) latsouth = xlat[i];
      if (xlat[i] > latnorth) latnorth = xlat[i];

      sena[i] = (short)(tmp_sena * 100.0);
      senz[i] = (short)(tmp_senz * 100.0);
      height[i] = (short)tmp_height;
    }  // if on-earth
  }    // pixel loop

  gsl_vector_free(C);

  return 0;
}

int mat2rpy(float pos[3], float vel[3], double smat[3][3], float rpy[3]) {
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
  double om[3][3];
  memcpy(&om[0][0], &x, 3 * sizeof(double));
  memcpy(&om[1][0], &y, 3 * sizeof(double));
  memcpy(&om[2][0], &z, 3 * sizeof(double));

  // Compute orbital-to-spacecraft matrix
  double rm[3][3];
  gsl_matrix_view rm_mat = gsl_matrix_view_array(&rm[0][0], 3, 3);

  int s;

  gsl_permutation *perm = gsl_permutation_alloc(3);

  // Compute the LU decomposition of this matrix
  gsl_matrix_view B = gsl_matrix_view_array(&om[0][0], 3, 3);
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
  double cosp = sqrt(rm[2][1] * rm[2][1] + rm[2][2] * rm[2][2]);
  if (rm[2][2] < 0) cosp *= -1;
  rpy[1] = RADEG * atan2(rm[2][0], cosp);
  rpy[2] = RADEG * atan(-rm[1][0] / rm[0][0]);

  return 0;
}


/*----------------------------------------------------------------- */
/* Create a Generic NETCDF4 file                                   */
/* ---------------------------------------------------------------- */
int geoFile::createFile( const char* filename, const char* cdlfile,
                         size_t sdim, int *ncid, int *gid,
                         uint32_t bbb, uint32_t rbb,
                         int16_t pcdim, int16_t psdim, size_t nswband, int32_t *rta_nadir,
                         bool radiance) {

   try {
     geofile = new NcFile( filename, NcFile::replace);
   }
   catch ( NcException& e) {
     e.what();
     cerr << "Failure creating OCI GEO file: " << filename << endl;
     exit(1);
   }

   ifstream output_data_structure;
   string line;
   string dataStructureFile;

   dataStructureFile.assign( cdlfile);
   expandEnvVar( &dataStructureFile);

   output_data_structure.open( dataStructureFile.c_str(), ifstream::in);
   if ( output_data_structure.fail() == true) {
     cout << "\"" << dataStructureFile.c_str() << "\" not found" << endl;
     exit(1);
   }

   // Find "dimensions" section of CDL file
   while(1) {
     getline( output_data_structure, line);
     size_t pos = line.find("dimensions:");
     if ( pos == 0) break;
   }

   // Define dimensions from "dimensions" section of CDL file
   ndims = 0;
   //   int dimid[1000];
   while(1) {
     getline( output_data_structure, line);
     boost::trim(line);
     if (line.substr(0,2) == "//") continue;

     size_t pos = line.find(" = ");
     if ( pos == string::npos) break;

     uint32_t dimSize;
     istringstream iss(line.substr(pos+2, string::npos));
     iss >> dimSize;

     iss.clear(); 
     iss.str( line);
     iss >> skipws >> line;

     if (line.compare("number_of_scans") == 0) {
       dimSize = sdim;
     }

     if (line.compare("blue_bands") == 0) {
       dimSize = bbb;
     }
          
     if (line.compare("red_bands") == 0) {
       dimSize = rbb;
     }

     if (line.compare("ccd_pixels") == 0) {
       dimSize = pcdim;
     }

     if (line.compare("SWIR_pixels") == 0) {
       dimSize = psdim;
     }
     
     if (line.compare("SWIR_bands") == 0) {
     dimSize = nswband;
     }
     try {
       ncDims[ndims++] = geofile->addDim( line, dimSize);
     }
     catch ( NcException& e) {
       e.what();
       cerr << "Failure creating dimension: " << line.c_str() << endl;
       exit(1);
     }

   } // while loop

   // Read global attributes (string attributes only) 
   while(1) {
     getline( output_data_structure, line);
     boost::trim(line);
     size_t pos = line.find("// global attributes");
     if ( pos == 0) break;
   }

   while(1) {
     getline( output_data_structure, line);
     size_t pos = line.find(" = ");
     if ( pos == string::npos) break;

     string attValue = line.substr(pos+3);

     // Remove any leading and trailing quotes
     attValue.erase(attValue.length()-2); // skip final " ;"
     size_t begQuote = attValue.find('"');
     size_t endQuote = attValue.find_last_of('"');
     if (begQuote == string::npos) continue; // Skip non-string global attributes
     attValue = attValue.substr(begQuote+1, endQuote-begQuote-1);

     istringstream iss(line.substr(pos+2));
     iss.clear(); 
     iss.str( line);
     iss >> skipws >> line;

     // Skip commented out attributes
     if (line.compare("//") == 0) continue;

     string attName;
     attName.assign(line.substr(1).c_str());

     try {
       geofile->putAtt(attName, attValue);
     }
     catch ( NcException& e) {
       e.what();
       cerr << "Failure creating attribute: " + attName << endl;
       exit(1);
     }
     
   } // while(1)

   geofile->putAtt("rta_nadir", NC_LONG, 2, rta_nadir);

   ngrps = 0;
   // Loop through groups
   while(1) {
     getline( output_data_structure, line);

     // Check if end of CDL file
     // If so then close CDL file and return
     if (line.substr(0,1).compare("}") == 0) {
       output_data_structure.close();
       return 0;
     }

     // Check for beginning of new group
     size_t pos = line.find("group:");

     // If found then create new group and variables
     if ( pos == 0) {

       // Parse group name
       istringstream iss(line.substr(6, string::npos));
       iss >> skipws >> line;

       // Create NCDF4 group
       ncGrps[ngrps++] = geofile->addGroup(line);

       int numDims=0;
       //       int varDims[NC_MAX_DIMS];
       //size_t dimSize[NC_MAX_DIMS];
       //char dimName[NC_MAX_NAME+1];
       string sname;
       string lname;
       string standard_name;
       string units;
       string description;
       string flag_masks;
       string flag_meanings;
       double valid_min=0.0;
       double valid_max=0.0;
       double scale=1.0;
       double offset=0.0;
       double fill_value=0.0;
       string coordinates;

       vector<NcDim> varVec;

       int ntype=0;
       NcType ncType;

       // Loop through datasets in group
       getline( output_data_structure, line); // skip "variables:"
       while(1) {
         getline( output_data_structure, line);
         boost::trim(line);

         if (line.substr(0,2) == "//") continue;
         if (line.length() == 0) continue;
         if (line.substr(0,1).compare("\r") == 0) continue;
         if (line.substr(0,1).compare("\n") == 0) continue;

         if (radiance) {
           if (line.find("rhot_") != std::string::npos) continue;
         } else {
           if (line.find("Lt_") != std::string::npos) continue;
         }
         
         size_t pos = line.find(":");

         // No ":" found, new dataset or empty line or end-of-group
         if ( pos == string::npos) {

           if ( numDims > 0) {
             // Create previous dataset
             createField( ncGrps[ngrps-1], sname.c_str(), lname.c_str(),
                          standard_name.c_str(), units.c_str(),
                          description.c_str(),
                          (void *) &fill_value, 
                          flag_masks.c_str(), flag_meanings.c_str(),
                          valid_min, valid_max, 
                          scale, offset,
                          ntype, varVec, coordinates);

             flag_masks.clear();
             flag_meanings.clear();
             units.clear();
             description.clear();
             varVec.clear();
             coordinates.clear();
           }

           valid_min=0.0;
           valid_max=0.0;
           scale=1.0;
           offset=0.0;
           fill_value=0.0;

           if (line.substr(0,10).compare("} // group") == 0) break;

           // Parse variable type
           string varType;
           istringstream iss(line);
           iss >> skipws >> varType;

           // Get corresponding NC variable type
           if ( varType.compare("char") == 0) ntype = NC_CHAR;
           else if ( varType.compare("byte") == 0) ntype = NC_BYTE;
           else if ( varType.compare("short") == 0) ntype = NC_SHORT;
           else if ( varType.compare("int") == 0) ntype = NC_INT;
           else if ( varType.compare("long") == 0) ntype = NC_INT;
           else if ( varType.compare("float") == 0) ntype = NC_FLOAT;
           else if ( varType.compare("real") == 0) ntype = NC_FLOAT;
           else if ( varType.compare("double") == 0) ntype = NC_DOUBLE;
           else if ( varType.compare("ubyte") == 0) ntype = NC_UBYTE;
           else if ( varType.compare("ushort") == 0) ntype = NC_USHORT;
           else if ( varType.compare("uint") == 0) ntype = NC_UINT;
           else if ( varType.compare("ulong") == 0) ntype = NC_UINT;
           else if ( varType.compare("int64") == 0) ntype = NC_INT64;
           else if ( varType.compare("uint64") == 0) ntype = NC_UINT64;

           // Parse short name (sname)
           pos = line.find("(");
           size_t posSname = line.substr(0, pos).rfind(" ");
           sname.assign(line.substr(posSname+1, pos-posSname-1));

           // Parse variable dimension info
           this->parseDims( line.substr(pos+1, string::npos), varVec);
           numDims = varVec.size();

         } else {
           // Parse variable attributes
           size_t posEql = line.find("=");
           size_t pos1qte = line.find("\"");
           size_t pos2qte = line.substr(pos1qte+1, string::npos).find("\"");

           string attrName = line.substr(pos+1, posEql-pos-2);

           // Get long_name
           if ( attrName.compare("long_name") == 0) {
             lname.assign(line.substr(pos1qte+1, pos2qte));
           }

           // Get units
           else if ( attrName.compare("units") == 0) {
             units.assign(line.substr(pos1qte+1, pos2qte));
           }

           // Get description
           else if ( attrName.compare("description") == 0) {
             description.assign(line.substr(pos1qte+1, pos2qte));
           }

           // Get _FillValue
           else if ( attrName.compare("_FillValue") == 0) {
             iss.clear(); 
             iss.str( line.substr(posEql+1, string::npos));
             iss >> fill_value;
           }

           // Get flag_masks
           else if ( attrName.compare("flag_masks") == 0) {
             flag_masks.assign(line.substr(pos1qte+1, pos2qte));
           }
           else if ( attrName.compare("flag_masks") == 0) {
             flag_masks.assign(line.substr(pos1qte+1, pos2qte));
           }

           // Get flag_meanings
           else if ( attrName.compare("flag_meanings") == 0) {
             flag_meanings.assign(line.substr(pos1qte+1, pos2qte));
           }

           // Get valid_min
           else if ( attrName.compare("valid_min") == 0) {
             iss.clear(); 
             iss.str( line.substr(posEql+1, string::npos));
             iss >> valid_min;
           }

           // Get valid_max
           else if ( attrName.compare("valid_max") == 0) {
             iss.clear(); 
             iss.str( line.substr(posEql+1, string::npos));
             iss >> valid_max;
           }

           // Get scale
           else if ( attrName.compare("scale_factor") == 0) {
             iss.clear(); 
             iss.str( line.substr(posEql+1, string::npos));
             iss >> scale;
           }

           // Get offset
           else if ( attrName.compare("add_offset") == 0) {
             iss.clear(); 
             iss.str( line.substr(posEql+1, string::npos));
             iss >> offset;
           }

           // Get coordinates
           else if ( attrName.compare("coordinates") == 0) {
             coordinates.assign(line.substr(pos1qte+1, pos2qte));
           }

         } // if ( pos == string::npos)
       } // datasets in group loop
     } // New Group loop
   } // Main Group loop

   return 0;
}


int geoFile::parseDims( string dimString, vector<NcDim>& varDims) {

  size_t curPos=0;
  //  char dimName[NC_MAX_NAME+1];
  string dimName;

  while(1) {
    size_t pos = dimString.find(",", curPos);
    if ( pos == string::npos) 
      pos = dimString.find(")");

    string varDimName;
    istringstream iss(dimString.substr(curPos, pos-curPos));
    iss >> skipws >> varDimName;

    for (int i=0; i<ndims; i++) {
      
      try {      
        dimName = ncDims[i].getName();
      }
      catch ( NcException& e) {
        e.what();
        cerr << "Failure accessing dimension: " + dimName << endl;
        exit(1);
      }
      
      if ( varDimName.compare(dimName) == 0) {
        varDims.push_back(ncDims[i]);
        break;
      }
    }
    if ( dimString.substr(pos, 1).compare(")") == 0) break;

    curPos = pos + 1;
  }

  return 0;
}

int geoFile::close() {

  try { 
    geofile->close();
  }
  catch ( NcException& e) {
    cout << e.what() << endl;
    cerr << "Failure closing: " + fileName << endl;
    exit(1);
  }
  
  return 0;
}

