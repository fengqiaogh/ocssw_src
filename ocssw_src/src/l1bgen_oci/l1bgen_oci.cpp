#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>

#include <sstream>
#include <fstream>
#include <iomanip>
#include <getopt.h>
#include <libgen.h>
#include <dirent.h>
#include <sys/stat.h>

#include "nc4utils.h"
#include "global_attrs.h"
#include "l1bgen_oci.h"

#include <clo.h>
#include <genutils.h>

using namespace std;
using namespace netCDF;
using namespace netCDF::exceptions;

#define VERSION "1.2"

// NOTE: we are comparing to the int val of solz (scale=0.01)
#define MAX_SOLZ (88 * 100) 

//    Modification history:
//  Programmer     Organization   Date     Ver   Description of change
//  ----------     ------------   ----     ---   ---------------------
//  Joel Gales     SAIC           04/29/20 0.01  Original development
//  Joel Gales     SAIC           01/13/21 0.02  Add support for SWIR
//  Joel Gales     SAIC           08/11/21 0.03  Add support for handling
//                                               fill values in science data
//  Joel Gales     SAIC           08/12/21 0.04  Add support for hysteresis
//                                               correction
//  Joel Gales     SAIC           09/23/21 0.045 Initialize uninitialized
//                                               variables
//  Joel Gales     SAIC           05/29/23 0.060 Implemented F.Patt code changes
//                                               (05/08/23) #272
//  Joel Gales     SAIC           06/21/23 0.061 Fix hysteresis bugs
//  Joel Gales     SAIC           06/23/23 0.062 Move get_agg_mat to common.cpp
//                                               Split off get_nib_nbb
//  Gwyn Fireman   SAIC           07/10/23 0.063  Read global metadata from json file
//  Joel Gales     SAIC           07/24/23 0.063 Add support for saturation
//  Joel Gales     SAIC           08/09/23 0.065 Add check_scan_times
//  Joel Gales     SAIC           08/25/23 0.700 Add support for reflectance
//                                               output
//  Joel Gales     SAIC           08/29/23 0.710 Convert thetap,thetas todeg
//  Joel Gales     SAIC           09/29/23 0.800 Call geolocation as function
//                                               rather than run beforehand
//  Joel Gales     SAIC           10/16/23 0.810 Add planarity correction
//  Joel Gales     SAIC           10/19/23 0.820 Fix metadata issues
//  Joel Gales     SAIC           12/06/23 0.823 Add rhot description attribute
//  Joel Gales     SAIC           12/10/23 0.840 Add tilt_home
//  Joel Gales     SAIC           01/02/24 0.860 Fix encoder interpolation bug
//  Joel Gales     SAIC           02/09/24 0.870 Add scan angle & polarization
//                                               coefficients

ofstream tempOut;

int main (int argc, char* argv[])
{

  cout << "l1bgen_oci " << VERSION << " (" 
       <<  __DATE__ << " " << __TIME__ << ")" << endl;

  clo_optionList_t *optionList = clo_createList();
  clo_addOption(optionList, "ifile", CLO_TYPE_IFILE, NULL, "Input L1A file");
  clo_addOption(optionList, "ofile", CLO_TYPE_OFILE, NULL, "Output L1B file");
  clo_addOption(optionList, "cal_lut", CLO_TYPE_IFILE, NULL, "CAL LUT file");
  clo_addOption(optionList, "geo_lut", CLO_TYPE_IFILE, NULL, "GEO LUT file");
  clo_addOption(optionList, "doi", CLO_TYPE_STRING, NULL,
                "Digital Object Identifier (DOI) string");
  clo_addOption(optionList, "pversion", CLO_TYPE_STRING, "Unspecified",
            "processing version string");
  clo_addOption(optionList, "demfile", CLO_TYPE_STRING,
                "$OCDATAROOT/common/gebco_ocssw_v2020.nc",
                "Digital elevation map file");
  clo_addOption(optionList, "radiance", CLO_TYPE_BOOL, "false",
                "Generate radiances");
  clo_addOption(optionList, "ephfile", CLO_TYPE_STRING, nullptr,
                "Definitive ephemeris file name");

  string l1a_filename;
  string l1b_filename;
  string cal_lut_filename;
  string geo_lut_filename;
  string dem_file;
  string doi;
  string pversion;
  string ephFile;

  bool radiance;
  
  if (argc == 1) {
    clo_printUsage(optionList);
    exit(EXIT_FAILURE);
  }
  clo_readArgs(optionList, argc, argv);
  
  // Grab the OPER LUTs
  string CALDIR = getenv("OCVARROOT");
  CALDIR.append("/oci/cal/OPER/");
  vector<string> lut_files;

  DIR *caldir;
  struct dirent *caldirptr;
  if ((caldir  = opendir(CALDIR.c_str())) != NULL) {
    while ((caldirptr = readdir(caldir)) != NULL) {
        lut_files.push_back(string(caldirptr->d_name));
    }
    closedir(caldir);
  }

  if (clo_isSet(optionList, "ifile")) {
    l1a_filename = clo_getString(optionList, "ifile");
  } else {
    clo_printUsage(optionList);
    exit(EXIT_FAILURE);
  }
  if (clo_isSet(optionList, "ofile")) {
    l1b_filename = clo_getString(optionList, "ofile");
  } else {
    clo_printUsage(optionList);
    exit(EXIT_FAILURE);
  }
  struct stat fstat_buffer;
  if (clo_isSet(optionList, "cal_lut")) {
    cal_lut_filename = clo_getString(optionList, "cal_lut");
  } else {
    for(const string& lut : lut_files) {
      if (lut.find("PACE_OCI_L1B_LUT") != std::string::npos) {
        cal_lut_filename = CALDIR;
        cal_lut_filename.append(lut);
        break;
      }
    }
  }
  if (cal_lut_filename.empty() || (stat (cal_lut_filename.c_str(), &fstat_buffer) != 0)) {
    cout << "Error: input CAL LUT file: " << cal_lut_filename.c_str() << " does not exist\n";
    exit(EXIT_FAILURE);
  }
  if (clo_isSet(optionList, "geo_lut")) {
    geo_lut_filename = clo_getString(optionList, "geo_lut");
  } else {
    for(const string& lut : lut_files) {
      if (lut.find("PACE_OCI_GEO_LUT") != std::string::npos) {
        geo_lut_filename = CALDIR;
        geo_lut_filename.append(lut);
        break;
      }
    }
  }
  if (geo_lut_filename.empty() || (stat (geo_lut_filename.c_str(), &fstat_buffer) != 0)) {
    cout << "Error: input GEO LUT file: " << geo_lut_filename.c_str() << " does not exist\n";
    exit(EXIT_FAILURE);
  }
  radiance = clo_getBool(optionList, "radiance");
  if (clo_isSet(optionList, "doi")) {
    doi = clo_getString(optionList, "doi");
    if (doi == "None"){
      doi.clear();
    }
  }

  if (clo_isSet(optionList, "ephfile")) {
    ephFile = clo_getString(optionList, "ephfile");
  }

  if (clo_isSet(optionList, "pvesion")) {
    pversion = clo_getString(optionList, "pversion");
  }

  char tmp_filename[FILENAME_MAX];
  parse_file_name(clo_getString(optionList, "demfile"), tmp_filename);
  dem_file = tmp_filename;
  if ((stat (dem_file.c_str(), &fstat_buffer) != 0)) {
    cout << "Error: DEM file: " << dem_file.c_str() << " does not exist\n";
    exit(EXIT_FAILURE);
  }
  // ********************************* //
  // *** Read calibration LUT file *** //
  // ********************************* //
  nc_set_chunk_cache(CHUNK_CACHE_SIZE, CHUNK_CACHE_NELEMS, CHUNK_CACHE_PREEMPTION);
  NcFile *calLUTfile = new NcFile( cal_lut_filename, NcFile::read);

  NcGroup gidCommon, gidBlue, gidRed, gidSWIR;
  gidCommon = calLUTfile->getGroup( "common");
  gidBlue = calLUTfile->getGroup( "blue");
  gidRed = calLUTfile->getGroup( "red");
  gidSWIR = calLUTfile->getGroup( "SWIR");
    
  float bwave[NBWAVE];
  float rwave[NRWAVE];
  float swave[NIWAVE];// = {939.716, 1038.317, 1250.375, 1248.550, 1378.169,
  // 1619.626, 1618.036, 2130.593, 2258.432};
  float bf0[NBWAVE];
  float rf0[NRWAVE];
  float sf0[NIWAVE];// = {81.959, 67.021, 44.350, 44.490, 35.551, 23.438, 23.476,
  //                       9.122, 7.427}; 
  float spass[NIWAVE] = {45, 80, 30, 30, 15, 75, 75, 50, 75};

  NcDim timeDim = calLUTfile->getDim("number_of_times");
  if(timeDim.isNull()) {
    cout << "Error: could not read number_of_times from " << cal_lut_filename << "\n";
    exit(EXIT_FAILURE);
  }
  size_t numTimes = timeDim.getSize();

  double *K2t = (double*)malloc(numTimes * sizeof(double));
  
  NcVar var;

  var = gidCommon.getVar( "blue_wavelength");
  var.getVar( bwave);

  var = gidCommon.getVar( "blue_F0");
  var.getVar( bf0);

  var = gidCommon.getVar( "red_wavelength");
  var.getVar( rwave);

  var = gidCommon.getVar( "red_F0");
  var.getVar( rf0);

  var = gidCommon.getVar( "SWIR_wavelength");
  var.getVar( swave);

  var = gidCommon.getVar( "SWIR_F0");
  var.getVar( sf0);

  var = gidCommon.getVar( "K2t");
  var.getVar( K2t);

  string tag;

  cal_lut_struct blue_lut;
  cal_lut_struct red_lut;
  cal_lut_struct swir_lut;

  uint32_t bbanddim, rbanddim, sbanddim, nldim, poldim;

  tag.assign("blue");
  read_oci_cal_lut( calLUTfile, tag, gidBlue, bbanddim, 1, nldim, poldim,
                    blue_lut);

  NcDim K3Tdim = calLUTfile->getDim("number_of_temperatures");

  float *K3T = new float[K3Tdim.getSize()];
  var = gidCommon.getVar( "K3T");
  var.getVar( K3T);

  tag.assign("red");
  read_oci_cal_lut( calLUTfile, tag, gidRed, rbanddim, 1, nldim, poldim,
                    red_lut);

  tag.assign("SWIR");
  read_oci_cal_lut( calLUTfile, tag, gidSWIR, sbanddim, 2, nldim, poldim,
                    swir_lut);

  // Read hysterisis parameters
  float hysttime[9][4];
  float hystamp[9][4];
  
  var = gidSWIR.getVar( "hyst_time_const");
  var.getVar( &hysttime[0][0]);
  var = gidSWIR.getVar( "hyst_amplitude");
  var.getVar( &hystamp[0][0]);
  
  calLUTfile->close();
  
  NcGroupAtt att;

  geo_struct geoLUT;
  geolocate_oci( l1a_filename, geo_lut_filename, geoLUT, l1b_filename, dem_file,
                 radiance, doi, ephFile, pversion);
  
  static l1bFile outfile;
  outfile.l1bfile = new NcFile( l1b_filename, NcFile::write);

  outfile.ncGrps[0] = outfile.l1bfile->getGroup( "sensor_band_parameters");
  outfile.ncGrps[1] = outfile.l1bfile->getGroup( "scan_line_attributes");
  outfile.ncGrps[2] = outfile.l1bfile->getGroup( "geolocation_data");
  outfile.ncGrps[3] = outfile.l1bfile->getGroup( "navigation_data");
  outfile.ncGrps[4] = outfile.l1bfile->getGroup( "observation_data");
  
  
  // Append call sequence to existing history
  string history = get_history(outfile.l1bfile);  // if updating existing file
  history.append(call_sequence(argc, argv));

  double distcorr;
  att = outfile.l1bfile->getAtt("earth_sun_distance_correction");
  att.getValues(&distcorr);

  NcDim nscanl1b_dim = outfile.l1bfile->getDim("number_of_scans");
  uint32_t nscanl1b = nscanl1b_dim.getSize();

  NcDim ccd_pixels_dim = outfile.l1bfile->getDim("ccd_pixels");
  uint32_t ccd_pixels = ccd_pixels_dim.getSize();

  double *evtime = new double[nscanl1b];
  var = outfile.ncGrps[1].getVar( "time");
  var.getVar( evtime);

  short *solz;
  float *csolz;
  unsigned char *qualFlag;
  
  // var = outfile.ncGrps[4].getVar( "rhot_blue");
  // if ( !var.isNull()) reflectance=true; else reflectance=false;

  if (!radiance) {
    solz = new short[nscanl1b*ccd_pixels];
    var = outfile.ncGrps[2].getVar( "solar_zenith");
    var.getVar( solz);

    csolz = new float[nscanl1b*ccd_pixels];
    for (size_t i=0; i<nscanl1b*ccd_pixels; i++)
      csolz[i] = cos(solz[i]*M_PI/180/100);
    
    qualFlag = new unsigned char[nscanl1b*ccd_pixels];
    var = outfile.ncGrps[2].getVar("quality_flag");
    var.getVar(qualFlag);
  }
  
  // Open and read data from L1Afile
  NcFile *l1afile = new NcFile( l1a_filename, NcFile::read);
   
  NcGroup ncGrps[4];

  ncGrps[0] = l1afile->getGroup( "scan_line_attributes");
  ncGrps[1] = l1afile->getGroup( "spatial_spectral_modes");
  ncGrps[2] = l1afile->getGroup( "engineering_data");
  ncGrps[3] = l1afile->getGroup( "science_data");

  // Get date (change this when year and day are added to time field)
  string tstart, tend;
  att = l1afile->getAtt("time_coverage_start");
  att.getValues(tstart);
  cout << "time_coverage_start: " << tstart << endl;

  att = l1afile->getAtt("time_coverage_end");
  att.getValues(tend);
  cout << "time_coverage_end:   " << tend << endl;

  uint16_t iyr, imn, idom;
  istringstream iss;

  iss.str(tstart.substr(0,4));
  iss >> iyr; iss.clear();
  iss.str(tstart.substr(5,2));
  iss >> imn; iss.clear();
  iss.str(tstart.substr(8,2));
  iss >> idom;
  int32_t jd = jday(iyr, imn, idom);

  int32_t iyr32, idy32;
  jdate( jd, &iyr32, &idy32);

  // Get numbers of blue and red bands
  NcDim blue_dim = l1afile->getDim("blue_bands");
  uint32_t bbands = blue_dim.getSize();
  NcDim red_dim = l1afile->getDim("red_bands");
  uint32_t rbands = red_dim.getSize();

  // Scan time, spin ID and HAM side
  NcDim nscan_dim = l1afile->getDim("number_of_scans");
  uint32_t nscan = nscan_dim.getSize();

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
  
  uint32_t nscan_good=0;
  for (size_t i=0; i<nscan; i++) {
    if ( spin[i] > 0) {
      sstime[nscan_good] = sstime[i];
      hside[nscan_good] = hside[i];
      nscan_good++;
    }
  }

  // Check for and fill in missing scan times
  short *tfl = new short[nscan_good]();
  check_scan_times( nscan_good, sstime, tfl);
 
  // ******************************************** //
  // *** Get spatial and spectral aggregation *** //
  // ******************************************** //
  NcDim ntaps_dim = l1afile->getDim("number_of_taps");
  uint32_t ntaps = ntaps_dim.getSize();
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

  int16_t *bagg = new int16_t[ntaps];
  var = ncGrps[1].getVar( "blue_spectral_mode");
  var.getVar( bagg);

  int16_t *ragg = new int16_t[ntaps];
  var = ncGrps[1].getVar( "red_spectral_mode");
  var.getVar( ragg);

  
  // ********************************************************************* //
  // *** Get # of EV lines/offset from scan start time to EV mid-time  *** //
  // ********************************************************************* //
  // This will be done by geolocation when integrated into L1B
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
                secpline, board_id, mspin, ot_10us, enc_count,
                &hamenc[0], rtaenc, iret);

  float clines[32400], slines[4050];
  uint16_t pcdim, psdim;
  //  int16_t iret;
  double ev_toff, deltc[32400], delts[4050];
  bool dark = false;
  get_ev( secpline, dtype, lines, iagg, pcdim, psdim, ev_toff, clines, slines,
          deltc, delts, dark, iret);
  if ( iret < 0) {
    cout << "No science collect in file: " << l1a_filename.c_str() << endl;
    l1afile->close();
    return 1;
  }

  size_t ka;
  for (size_t i=0; i<spatzn; i++) {
    if ( dtype[i] != 0 && dtype[i] != 2 && dtype[i] != 10) {
      ka = i;
      break;
    }
  }


  // *********************************************************** //
  // *** Generate matrices for spectral and gain aggregation *** //
  // *********************************************************** //
  size_t *ia = new size_t[ntaps];
  uint32_t ntb[16];
  
  // Blue bands
  uint32_t bib=1, bbb=1;

  if (get_nib_nbb( ntaps, ia, ntb, bagg, bib, bbb) == 2)
    cout << "All blue taps disabled" << endl;
  
  // Note: bgmat & rgmat not necessarily contiguous
  float **bamat = new float*[bib];  
  float **bgmat = new float*[512];

  if (bib != 1)
    get_agg_mat( ia, bagg, ntb, bib, bbb, bamat, bgmat);

  if (bib != bbands) {
    cout << "Number of blue bands in file: " << l1a_filename.c_str() <<
      " not consistent with spectral aggregation" << endl;
    l1afile->close();
    return 1;
  } else if ( bib < 4) cout << "No blue bands in file: " <<
                         l1a_filename.c_str() << endl;
    
    
  // Red bands
  uint32_t rib=1, rbb=1;

  if (get_nib_nbb( ntaps, ia, ntb, ragg, rib, rbb) == 2)
    cout << "All red taps disabled" << endl;
  
  float **ramat = new float*[rib];  
  float **rgmat = new float*[512];

  if (rib != 1)
    get_agg_mat( ia, ragg, ntb, rib, rbb, ramat, rgmat);

  if (rib != rbands) {
    cout << "Number of red bands in file: " << l1a_filename.c_str() <<
      " not consistent with spectral aggregation" << endl;
    l1afile->close();
    return 1;
  } else if ( rib < 4) cout << "No red bands in file: " << l1a_filename.c_str()
                            << endl;

  uint16_t swb = 9;


  // ********************************* //
  // *** Get dark collect location *** //
  // ********************************* //
  int16_t kd=-1;
  for (size_t i=0; i<spatzn; i++) {
    if ( dtype[i] == 2) {
      kd = (int16_t) i;
    }
  }
  if ( kd == -1) {
    cout << "No dark collect in file: " << l1a_filename.c_str() << endl;
    l1afile->close();
    return 1;
  }

  int16_t ldc=0, lds=0;
  for (size_t i=0; i<(size_t) kd; i++) {
    if ( dtype[i] != 0 && dtype [1] != 10) {
      ldc += lines[i] / iagg[i];
      lds += lines[i] / 8;
    }
  }
  int16_t ndc = lines[kd] / iagg[kd];
  int16_t nds = lines[kd] / 8;


  // *********************************************************************** //
  // *** Generate band gain structs from LUTs, date/time & gain matrices *** //
  // *********************************************************************** //
  gains_struct bgains;
  gains_struct rgains;
  gains_struct sgains;

  // Note: bgmat & rgmat not necessarily contiguous
  // That's why we pass the location of the array of pointers
  if ( bib >= 4)
    make_oci_gains( bib, bbanddim, iyr, jd, evtime[0], numTimes, K2t, -1,
                    iagg[ka], bagg, blue_lut, &bgmat[0], bgains);

  if ( rib >= 4)
    make_oci_gains( rib, rbanddim, iyr, jd, evtime[0], numTimes, K2t, -1,
                    iagg[ka], ragg, red_lut, &rgmat[0], rgains);
  
  float **sgmat = new float*[swb];
  for (size_t i=0; i<swb; i++) {
    sgmat[i] = new float[swb];
    for (size_t j=0; j<swb; j++) {
      if (i == j) sgmat[i][j] = 1.0; else sgmat[i][j] = 0.0;
    }
  }
  make_oci_gains( swb, swb, iyr, jd, evtime[0], numTimes, K2t, board_id, -1, NULL,
                  swir_lut, &sgmat[0], sgains);

  // Compute bibf0,b1bf0
  float *bibf0 = new float[bib];
  for (size_t i=0; i<bib; i++) {
    bibf0[i] = 0.0;
    for (size_t j=0; j<512; j++) {
      bibf0[i] += bf0[j]*bgmat[j][i];
    }
  }
    
  float *b1bf0 = new float[bbb];
  for (size_t i=0; i<bbb; i++) {
    b1bf0[i] = 0.0;
    for (size_t j=0; j<bib; j++) {
      b1bf0[i] += bibf0[j]*bamat[j][i];
    }
  }

  // Compute ribf0,r1bf0
  float *ribf0 = new float[rib];
  for (size_t i=0; i<rib; i++) {
    ribf0[i] = 0.0;
    for (size_t j=0; j<512; j++) {
      ribf0[i] += rf0[j]*rgmat[j][i];
    }
  }

  float *r1bf0 = new float[rbb];
  for (size_t i=0; i<rbb; i++) {
    r1bf0[i] = 0.0;
    for (size_t j=0; j<rib; j++) {
      r1bf0[i] += ribf0[j]*ramat[j][i];
    }
  }
  
  // Read selected temperature fields and interpolate to scan times
  uint16_t ntemps = bgains.ldims[1] + sgains.ldims[1];

  float **caltemps = new float *[nscan_good];
  caltemps[0] = new float[ntemps*nscan_good];
  for (size_t i=1; i<nscan_good; i++) caltemps[i] = caltemps[i-1] + ntemps;

  for (size_t i=0; i<nscan_good; i++)
    for (size_t j=1; j<ntemps; j++)
      caltemps[i][j] = 0.0;
  get_oci_cal_temps( l1afile, ncGrps[2], ntemps, nscan_good, evtime, caltemps);
  uint16_t nctemps = bgains.ldims[1];


  // Read dark collects from science data arrays
  vector<size_t> start, count;
  start.push_back(0);
  start.push_back(0);
  start.push_back(ldc);

  size_t dims[3];
  dims[0] = nscan_good; dims[1] = bib; dims[2] = ndc;
  uint16_t ***bdark = make3dT<uint16_t>(dims);
  dims[0] = nscan_good; dims[1] = rib; dims[2] = ndc;
  uint16_t ***rdark = make3dT<uint16_t>(dims);

  uint16_t bfill;
  uint16_t rfill;
  
  if ( bib > 4) {
    count.push_back(nscan_good);
    count.push_back(bib);
    count.push_back(ndc);
    
    var = ncGrps[3].getVar( "sci_blue");
    var.getVar( start, count, &bdark[0][0][0]);

    var.getAtt("_FillValue").getValues(&bfill);
  }

  if ( rib > 4) {
    count.clear();
    count.push_back(nscan_good);
    count.push_back(rib);
    count.push_back(ndc);
    
    var = ncGrps[3].getVar( "sci_red");
    var.getVar( start, count, &rdark[0][0][0]);

    var.getAtt("_FillValue").getValues(&rfill);
  }

  dims[0] = nscan_good; dims[1] = swb; dims[2] = nds;
  uint32_t ***sdark = make3dT<uint32_t>(dims);

  uint32_t sfill;
  
  start.clear();
  start.push_back(0);
  start.push_back(0);
  start.push_back(lds);

  count.clear();
  count.push_back(nscan_good);
  count.push_back(swb);
  count.push_back(nds);

  var = ncGrps[3].getVar( "sci_SWIR");
  var.getVar( start, count, &sdark[0][0][0]);

  var.getAtt("_FillValue").getValues(&sfill);

  uint32_t fill32;
  
  // number of scans of dark data to average; will make this an input parameter
  uint16_t ndsc = 1;
  //number of dark pixels to skip; will make this an input parameter
  uint16_t nskp = 0;

  uint32_t bicount[3] = {1,bib,(uint32_t) pcdim};
  uint32_t ricount[3] = {1,rib,(uint32_t) pcdim};
  uint32_t sicount[3] = {1,swb,(uint32_t) psdim};
  uint32_t bcount[3] = {bbb,1,(uint32_t) pcdim};
  uint32_t rcount[3] = {rbb,1,(uint32_t) pcdim};
  uint32_t scount[3] = {swb,1,(uint32_t) psdim};

  // Calibrated data variables
  float **bdn = new float *[bib];
  bdn[0] = new float[pcdim*bib];
  for (size_t i=1; i<bib; i++) bdn[i] = bdn[i-1] + pcdim;

  float **rdn = new float *[rib];
  rdn[0] = new float[pcdim*rib];
  for (size_t i=1; i<rib; i++) rdn[i] = rdn[i-1] + pcdim;

  float **sdn = new float *[swb];
  sdn[0] = new float[psdim*swb];
  for (size_t i=1; i<swb; i++) sdn[i] = sdn[i-1] + psdim;

  float **bcal = new float *[bib];
  bcal[0] = new float[pcdim*bib];
  for (size_t i=1; i<bib; i++) bcal[i] = bcal[i-1] + pcdim;

  float **rcal = new float *[rib];
  rcal[0] = new float[pcdim*rib];
  for (size_t i=1; i<rib; i++) rcal[i] = rcal[i-1] + pcdim;

  float **scal = new float *[swb];
  scal[0] = new float[psdim*swb];
  for (size_t i=1; i<swb; i++) scal[i] = scal[i-1] + psdim;

  // Saturation arrays
  uint8_t *bsat = new uint8_t[bib];
  uint8_t **bqual = new uint8_t *[bbb];
  bqual[0] = new uint8_t[pcdim*bbb];
  for (size_t i=1; i<bbb; i++) bqual[i] = bqual[i-1] + pcdim;

  uint8_t *rsat = new uint8_t[rib];
  uint8_t **rqual = new uint8_t *[rbb];
  rqual[0] = new uint8_t[pcdim*rbb];
  for (size_t i=1; i<rbb; i++) rqual[i] = rqual[i-1] + pcdim;

  uint8_t *ssat = new uint8_t[swb];
  uint8_t **squal = new uint8_t *[swb];
  squal[0] = new uint8_t[psdim*swb];
  for (size_t i=1; i<swb; i++) squal[i] = squal[i-1] + psdim;

  
  double *thetap = new double[pcdim];
  double *thetas = new double[psdim];

  float **pview = new float *[pcdim];
  pview[0] = new float[3*pcdim];
  for (size_t i=1; i<pcdim; i++) pview[i] = pview[i-1] + 3;

  float **sview = new float *[psdim];
  sview[0] = new float[3*psdim];
  for (size_t i=1; i<psdim; i++) sview[i] = sview[i-1] + 3;

  uint16_t **bsci = new uint16_t *[bib];
  bsci[0] = new uint16_t[pcdim*bib];
  for (size_t i=1; i<bib; i++) bsci[i] = bsci[i-1] + pcdim;

  uint16_t **rsci = new uint16_t *[rib];
  rsci[0] = new uint16_t[pcdim*rib];
  for (size_t i=1; i<rib; i++) rsci[i] = rsci[i-1] + pcdim;

  uint32_t **ssci = new uint32_t *[swb];
  ssci[0] = new uint32_t[psdim*swb];
  for (size_t i=1; i<swb; i++) ssci[i] = ssci[i-1] + psdim;

  float **bcalb = new float *[bbb];
  bcalb[0] = new float[pcdim*bbb];
  for (size_t i=1; i<bbb; i++) bcalb[i] = bcalb[i-1] + pcdim;

  float **rcalb = new float *[rbb];
  rcalb[0] = new float[pcdim*rbb];
  for (size_t i=1; i<rbb; i++) rcalb[i] = rcalb[i-1] + pcdim;

  ///////////////
  // Main loop //
  ///////////////
  
  // Read, calibrate and write science data
  for (size_t iscn=0; iscn<nscan_good; iscn++) {

    if ((iscn % 50) == 0) cout << "Calibrating scan: " << iscn << endl;

    // Check for valid mirror side
    if ( hside[iscn] == 0 || hside[iscn] == 1) {

      // Get scan angle
      get_oci_vecs( nscan, pcdim, geoLUT.as_planarity, geoLUT.at_planarity,
                    geoLUT.rta_nadir, geoLUT.ham_ct_angles,
                    ev_toff, spin[iscn], hside[iscn],
                    clines, deltc, revpsec, ppr_off, board_id, nmcescan, mspin,
                    enc_count, &hamenc[0], &rtaenc[0], pview, thetap, iret);

      for (size_t k=0; k<pcdim; k++) thetap[k] *= 180/M_PI;

      get_oci_vecs( nscan, psdim, geoLUT.as_planarity, geoLUT.at_planarity,
                    geoLUT.rta_nadir, geoLUT.ham_ct_angles,
                    ev_toff, spin[iscn], hside[iscn],
                    slines, delts, revpsec, ppr_off, board_id, nmcescan, mspin,
                    enc_count, &hamenc[0], &rtaenc[0], sview, thetas, iret);

      for (size_t k=0; k<psdim; k++) thetas[k] *= 180/M_PI;

      // Write scan angles
      start.clear();
      start.push_back(iscn);
      start.push_back(0);

      count.clear();
      count.push_back(1);
      count.push_back(pcdim);

      var = outfile.ncGrps[3].getVar( "CCD_scan_angles");
      var.putVar( start, count, thetap);

      count[1] = psdim;
      var = outfile.ncGrps[3].getVar( "SWIR_scan_angles");
      var.putVar( start, count, thetas);

      //  Blue bands
      if (bib >= 4) {

        start.clear();
        start.push_back(iscn);
        start.push_back(0);
        start.push_back(0);

        count.clear();
        count.push_back(bicount[0]);
        count.push_back(bicount[1]);
        count.push_back(bicount[2]);

        
        var = ncGrps[3].getVar( "sci_blue");
        var.getVar( start, count, bsci[0]);

        // Compute dark offset, correct data, and apply absolute and
        // temporal gain and temperature correction
        float *bdc = new float[bib];

        fill32 = bfill;
        int16_t iret;
        get_oci_dark<uint16_t>( iscn, nscan_good, hside, ndsc, nskp,
                                iagg[ka], iagg[kd], ntaps, bagg, fill32,
                                ndc, bdark, bib, bdc, iret);

        float *k3 = new float[bib];
        get_oci_temp_corr( bib, bgains, K3T, caltemps[iscn], nscan_good, k3);
        for (size_t j=0; j<bib; j++) {
          for (size_t k=0; k<pcdim; k++) {

            // Handle fill value
            if (bsci[j][k] == bfill) {
              bdn[j][k] = -32767;
              bcal[j][k] = -32767;
              continue;
            }

            // Need to save dn for linearity correction
            bdn[j][k] = bsci[j][k] - bdc[j];
            bcal[j][k] = k3[j] * bgains.K1K2[j][hside[iscn]] * bdn[j][k];
          }
        }

        delete [] k3;
        delete [] bdc;
        
        // Compute and apply RVS and linearity
        float **k4 = new float *[bib];
        k4[0] = new float[pcdim*bib];
        for (size_t i=1; i<bib; i++) k4[i] = k4[i-1] + pcdim;
        get_oci_rvs_corr( bib, pcdim, hside[iscn], bgains, thetap, k4);

        float **k5 = new float *[bib];
        k5[0] = new float[pcdim*bib];
        for (size_t i=1; i<bib; i++) k5[i] = k5[i-1] + pcdim;
        get_oci_lin_corr( bib, pcdim, nldim, bgains, bdn, k5);

        for (size_t j=0; j<bib; j++) {
          for (size_t k=0; k<pcdim; k++) {
            if(bcal[j][k] != -32767)
              bcal[j][k] *= k4[j][k] * k5[j][k];
          }
        }

        delete [] k4[0];
        delete [] k4;
        delete [] k5[0];
        delete [] k5;

        // Aggregate to L1B bands
        //bcalb = transpose(bamat#transpose(bcal))
        for (size_t j=0; j<bbb; j++) {
          for (size_t k=0; k<pcdim; k++) {
            float sum = 0.0;
            for (size_t l=0; l<bib; l++)
              if (bcal[l][k] != -32767) sum += bamat[l][j]*bcal[l][k];
            bcalb[j][k] = sum;
            if (!radiance) {
              if(solz[iscn*pcdim+k] < MAX_SOLZ) // NOTE: solz is short with scale of 0.01
                bcalb[j][k] *= M_PI*distcorr/(b1bf0[j]*csolz[iscn*pcdim+k]);
              else
                bcalb[j][k] = -32767;
            }
          }
        }

        // Check for saturation
        for (size_t k=0; k<pcdim; k++) {
          for (size_t j=0; j<bib; j++) {
            bsat[j] = 0;
            if ( bdn[j][k] >= bgains.sat_thres[j]) bsat[j] = 1;
          }
          for (size_t j=0; j<bbb; j++) {
            float sum = 0.0;
            bqual[j][k] = 0;
            for (size_t l=0; l<bib; l++) sum += bamat[l][j]*bsat[l];
            if ( sum > 0) bqual[j][k] = 1;
          }
        }
        
        start.clear();
        start.push_back(0);
        start.push_back(iscn);
        start.push_back(0);

        count.clear();
        count.push_back(bcount[0]);
        count.push_back(bcount[1]);
        count.push_back(bcount[2]);
        
        // Output to L1B file
        if (!radiance)
          var = outfile.ncGrps[4].getVar( "rhot_blue");
        else
          var = outfile.ncGrps[4].getVar( "Lt_blue");
        var.putVar( start, count, &bcalb[0][0]);

        var = outfile.ncGrps[4].getVar( "qual_blue");
        var.putVar( start, count, &bqual[0][0]);
      } // End blue


      //  Red bands
      if (rib >= 4) {

        start.clear();
        start.push_back(iscn);
        start.push_back(0);
        start.push_back(0);

        count.clear();
        count.push_back(ricount[0]);
        count.push_back(ricount[1]);
        count.push_back(ricount[2]);

        var = ncGrps[3].getVar( "sci_red");
        var.getVar( start, count, rsci[0]);

        // Compute dark offset, correct data, and apply absolute and
        // temporal gain and temperature correction
        float *rdc = new float[rib];
        fill32 = rfill;
        int16_t iret;
        get_oci_dark<uint16_t>( iscn, nscan_good, hside, ndsc, nskp,
                                iagg[ka], iagg[kd], ntaps, ragg, fill32,
                                ndc, rdark, rib, rdc, iret);

        float *k3 = new float[rib];
        get_oci_temp_corr( rib, rgains, K3T, caltemps[iscn], nscan_good, k3);

        for (size_t j=0; j<rib; j++) {
          for (size_t k=0; k<pcdim; k++) {

            // Handle fill value
            if (rsci[j][k] == rfill) {
              rdn[j][k] = -32767;
              rcal[j][k] = -32767;
              continue;
            }
            
            // Need to save dn for linearity correction
            rdn[j][k] = rsci[j][k] - rdc[j];
            rcal[j][k] = k3[j] * rgains.K1K2[j][hside[iscn]] * rdn[j][k];
          }
        }

        delete [] k3;
        delete [] rdc;
        
        // Compute and apply RVS and linearity
        float **k4 = new float *[rib];
        k4[0] = new float[pcdim*rib];
        for (size_t i=1; i<rib; i++) k4[i] = k4[i-1] + pcdim;
        get_oci_rvs_corr( rib, pcdim, hside[iscn], rgains, thetap, k4);

        float **k5 = new float *[rib];
        k5[0] = new float[pcdim*rib];
        for (size_t i=1; i<rib; i++) k5[i] = k5[i-1] + pcdim;

        get_oci_lin_corr( rib, pcdim, nldim, rgains, rdn, k5);

        for (size_t j=0; j<rib; j++) {
          for (size_t k=0; k<pcdim; k++) {
            if(rcal[j][k] != -32767)
              rcal[j][k] *= k4[j][k] * k5[j][k];
          }
        }
        
        delete [] k4[0];
        delete [] k4;
        delete [] k5[0];
        delete [] k5;

        // Aggregate to L1B bands
        for (size_t j=0; j<rbb; j++) {
          for (size_t k=0; k<pcdim; k++) {
            float sum = 0.0;
            for (size_t l=0; l<rib; l++)
              if (rcal[l][k] != -32767) sum += ramat[l][j]*rcal[l][k];
            rcalb[j][k] = sum;
            if (!radiance) {
              if(solz[iscn*pcdim+k] < MAX_SOLZ) // NOTE: solz is short with scale of 0.01
                rcalb[j][k] *= M_PI*distcorr/(r1bf0[j]*csolz[iscn*pcdim+k]);
              else
                rcalb[j][k] = -32767;
            }
          }
        }

        // Check for saturation
        for (size_t k=0; k<pcdim; k++) {
          for (size_t j=0; j<rib; j++) {
            rsat[j] = 0;
            if ( rdn[j][k] >= rgains.sat_thres[j]) rsat[j] = 1;
          }
          for (size_t j=0; j<rbb; j++) {
            float sum = 0.0;
            rqual[j][k] = 0;
            for (size_t l=0; l<rib; l++) sum += ramat[l][j]*rsat[l];
            if ( sum > 0) rqual[j][k] = 1;
          }
        }

        start.clear();
        start.push_back(0);
        start.push_back(iscn);
        start.push_back(0);

        count.clear();
        count.push_back(rcount[0]);
        count.push_back(rcount[1]);
        count.push_back(rcount[2]);
  
        // Output to L1B file
        if (!radiance)
          var = outfile.ncGrps[4].getVar( "rhot_red");
        else
          var = outfile.ncGrps[4].getVar( "Lt_red");
        var.putVar( start, count, &rcalb[0][0]);
        
        var = outfile.ncGrps[4].getVar( "qual_red");
        var.putVar( start, count, &rqual[0][0]);

      } // End red

      //  SWIR bands
      start.clear();
      start.push_back(iscn);
      start.push_back(0);
      start.push_back(0);

      count.clear();
      count.push_back(sicount[0]);
      count.push_back(sicount[1]);
      count.push_back(sicount[2]);

      var = ncGrps[3].getVar( "sci_SWIR");
      var.getVar( start, count, ssci[0]);

      // Compute dark offset, correct data, and apply absolute and
      // temporal gain and temperature correction
      float *sdc = new float[swb];
      int16_t sagg = 1;
      int16_t iret;

      float *k3 = new float[swb];

      float **k4 = new float *[swb];
      k4[0] = new float[psdim*swb];
      for (size_t i=1; i<swb; i++) k4[i] = k4[i-1] + psdim;

      float **k5 = new float *[swb];
      k5[0] = new float[psdim*swb];
      for (size_t i=1; i<swb; i++) k5[i] = k5[i-1] + psdim;


      // initialize k4, k5 and sdn
      // for (int band=0; band<swb; band++) {
      //   for (int pix=0; pix<psdim; pix++) {
      //     k4[band][pix] = 0.0;
      //     k5[band][pix] = 0.0;
      //     sdn[band][pix] = 0.0;
      //   }
      // }
      for (int pix=0; pix<psdim*swb; pix++) {
          k4[0][pix] = 0.0;
          k5[0][pix] = 0.0;
          sdn[0][pix] = 0.0;
      }

      get_oci_dark<uint32_t>( iscn, nscan_good, hside, ndsc, nskp,
                              1, 1, 1, &sagg, sfill, nds, sdark,
                              swb, sdc, iret);

      if ( iret != -1) {
        get_oci_temp_corr( swb, sgains, &K3T[nctemps], &caltemps[iscn][nctemps],
                           nscan_good, k3);

        // Compute and apply RVS and linearity
        get_oci_rvs_corr( swb, psdim, hside[iscn], sgains, thetas, k4);
        get_oci_lin_corr( swb, psdim, nldim, sgains, sdn, k5);

        
        // SWIR band loop
        for (size_t j=0; j<swb; j++) {

          uint32_t goodcnt=0;
          int32_t *indx = new int32_t[psdim];
          for (size_t k=0; k<psdim; k++) {

            // radiance == true, do not check qualFlag bc it does not retrieve it
            if (radiance && ssci[j][k] != sfill) {
              indx[goodcnt++] = k;
            }
            // radiance == false, check quality flag 
            else if (!radiance && ssci[j][k] != sfill && (int)qualFlag[iscn*pcdim+k] == 0) {
              indx[goodcnt++] = k;
            } 
            else {
              // Handle fill value
              sdn[j][k] = -32767;
              scal[j][k] = -32767;
            }
          }

          // Hysteresis correction
          float hc_prev[4]={0,0,0,0};
          float hc[4]={0,0,0,0};
          float *hyst = new float[psdim];

          for (size_t k=0; k<goodcnt; k++) {
            // Need to save dn for linearity correction
            sdn[j][indx[k]] = ssci[j][indx[k]] - sdc[j];

            if ( k > 0) {
              hyst[k] = 0.0;
              for (size_t l=0; l<3; l++) {
                // Compute exponential decay constants
                float e = exp(-1.0/hysttime[j][l]);

                hc[l] = hc_prev[l]*e + sdn[j][indx[k-1]]*hystamp[j][l];
                hyst[k] += hc[l];

                hc_prev[l] = hc[l];
              } // l-loop
            } else {
              hyst[k] = 0.0;
            }
          } // (reduced) k-loop

          for (size_t k=0; k<goodcnt; k++) {
            scal[j][indx[k]] =
              k3[j] * sgains.K1K2[j][hside[iscn]] * (sdn[j][indx[k]] - hyst[k]);

            scal[j][indx[k]] *= k5[j][indx[k]]*k4[j][indx[k]];

            if (!radiance) {
              if(solz[iscn*pcdim+indx[k]] < MAX_SOLZ) // NOTE: solz is short with scale of 0.01
                scal[j][indx[k]] *= M_PI*distcorr/(sf0[j]*csolz[iscn*psdim+indx[k]]);
              else
                scal[j][indx[k]] = -32767;
            }
          } // (reduced) k-loop

          delete [] indx;
          delete [] hyst;
        } // j-loop (band)
      } // iret != -1

      
      // Check for saturation
      for (size_t k=0; k<psdim; k++) {
        for (size_t j=0; j<swb; j++) {
          ssat[j] = 0;
          if ( ssci[j][k] >= sgains.sat_thres[j]) ssat[j] = 1;
        }
        for (size_t j=0; j<swb; j++) {
          squal[j][k] = ssat[j];
        }
      }

      delete [] k3;
      delete [] sdc;
      
      delete [] k4[0];
      delete [] k4;
      delete [] k5[0];
      delete [] k5;

      start.clear();
      start.push_back(0);
      start.push_back(iscn);
      start.push_back(0);

      count.clear();
      count.push_back(scount[0]);
      count.push_back(scount[1]);
      count.push_back(scount[2]);
  
      // Output to L1B file
      if (!radiance)
        var = outfile.ncGrps[4].getVar( "rhot_SWIR");
      else
        var = outfile.ncGrps[4].getVar( "Lt_SWIR");
      var.putVar( start, count, &scal[0][0]);

      var = outfile.ncGrps[4].getVar( "qual_SWIR");
      var.putVar( start, count, &squal[0][0]);
        
    } else {
      cout << "No mirror side index for scan: " << iscn << endl;
    } // Check for valid mirror side
  } // Scan loop

  // End Main loop
  
  // Write spectral band information
  // Calculate band centers for aggregated hyperspectral bands
  
  if (bib >= 4) {
    // bgmat#bwave
    float *b1bwave = new float[bbb];
    float *sum = new float[bib];
    for (size_t i=0; i<bib; i++) {
      sum[i] = 0.0;
      for (size_t j=0; j<512; j++) {
        sum[i] += bwave[j]*bgmat[j][i];
      }
    }

    // bamat#sum
    for (size_t i=0; i<bbb; i++) {
      b1bwave[i] = 0.0;
      for (size_t j=0; j<bib; j++) {
        b1bwave[i] += bamat[j][i]*sum[j];
      }
    }

    start.clear();
    start.push_back(0);

    count.clear();
    count.push_back(bbb);
    var = outfile.ncGrps[0].getVar( "blue_wavelength");
    var.putVar( start, count, b1bwave);

    var = outfile.ncGrps[0].getVar( "blue_solar_irradiance");
    var.putVar( start, count, b1bf0);

    // Add m12_coef
    // b1bm12[ip,im,*] = bamat#bgmat#transpose(blue_lut.m12_coef[ip,im,*])
    // blue_m12_coef(blue_bands, HAM_sides, polarization_coefficients) ;

    dims[0] = bbb; dims[1] = 2; dims[2] = 3;
    float ***b1bm12 = make3dT<float>(dims);
    float ***b1bm13 = make3dT<float>(dims);

    for (size_t l=0; l<3; l++) {
      for (size_t m=0; m<2; m++) {

        for (size_t i=0; i<bib; i++) {
          sum[i] = 0.0;
          for (size_t j=0; j<NBWAVE; j++) {
            sum[i] += blue_lut.m12_coef[j][m][l]*bgmat[j][i];
          }
        }

        // bamat#sum
        for (size_t i=0; i<bbb; i++) {
          b1bm12[i][m][l] = 0.0;
          for (size_t j=0; j<bib; j++) {
            b1bm12[i][m][l] += bamat[j][i]*sum[j];
          }
        }

        for (size_t i=0; i<bib; i++) {
          sum[i] = 0.0;
          for (size_t j=0; j<NBWAVE; j++) {
            sum[i] += blue_lut.m13_coef[j][m][l]*bgmat[j][i];
          }
        }

        // bamat#sum
        for (size_t i=0; i<bbb; i++) {
          b1bm13[i][m][l] = 0.0;
          for (size_t j=0; j<bib; j++) {
            b1bm13[i][m][l] += bamat[j][i]*sum[j];
          }
        }
        
      }
    }

    start.clear();
    start.push_back(0); start.push_back(0); start.push_back(0);

    count.clear();
    count.push_back(bbb); count.push_back(2); count.push_back(3);
    
    var = outfile.ncGrps[0].getVar( "blue_m12_coef");
    var.putVar( start, count, &b1bm12[0][0][0]);
    var = outfile.ncGrps[0].getVar( "blue_m13_coef");
    var.putVar( start, count, &b1bm13[0][0][0]);
    
    delete [] b1bwave;
    delete [] sum;

    delete [] b1bm12[0][0];
    delete [] b1bm13[0][0];
  }

  if (rib >= 4) {
    float *r1bwave = new float[rbb];
    float *sum = new float[rib];
    for (size_t i=0; i<rib; i++) {
      sum[i] = 0.0;
      for (size_t j=0; j<512; j++) {
        sum[i] += rwave[j]*rgmat[j][i];
      }
    }

    for (size_t i=0; i<rbb; i++) {
      r1bwave[i] = 0.0;
      for (size_t j=0; j<rib; j++) {
        r1bwave[i] += ramat[j][i]*sum[j];
      }
    }

    start.clear();
    start.push_back(0);

    count.clear();
    count.push_back(rbb);
    var = outfile.ncGrps[0].getVar( "red_wavelength");
    var.putVar( start, count, r1bwave);

    var = outfile.ncGrps[0].getVar( "red_solar_irradiance");
    var.putVar( start, count, r1bf0);

    dims[0] = rbb; dims[1] = 2; dims[2] = 3;
    float ***r1bm12 = make3dT<float>(dims);
    float ***r1bm13 = make3dT<float>(dims);

    for (size_t l=0; l<3; l++) {
      for (size_t m=0; m<2; m++) {

        for (size_t i=0; i<rib; i++) {
          sum[i] = 0.0;
          for (size_t j=0; j<NRWAVE; j++) {
            sum[i] += red_lut.m12_coef[j][m][l]*rgmat[j][i];
          }
        }

        // bamat#sum
        for (size_t i=0; i<rbb; i++) {
          r1bm12[i][m][l] = 0.0;
          for (size_t j=0; j<rib; j++) {
            r1bm12[i][m][l] += ramat[j][i]*sum[j];
          }
        }

        for (size_t i=0; i<rib; i++) {
          sum[i] = 0.0;
          for (size_t j=0; j<NRWAVE; j++) {
            sum[i] += red_lut.m13_coef[j][m][l]*rgmat[j][i];
          }
        }

        // bamat#sum
        for (size_t i=0; i<rbb; i++) {
          r1bm13[i][m][l] = 0.0;
          for (size_t j=0; j<rib; j++) {
            r1bm13[i][m][l] += ramat[j][i]*sum[j];
          }
        }
        
      }
    }

    start.clear();
    start.push_back(0); start.push_back(0); start.push_back(0);

    count.clear();
    count.push_back(rbb); count.push_back(2); count.push_back(3);
    
    var = outfile.ncGrps[0].getVar( "red_m12_coef");
    var.putVar( start, count, &r1bm12[0][0][0]);
    var = outfile.ncGrps[0].getVar( "red_m13_coef");
    var.putVar( start, count, &r1bm13[0][0][0]);
    
    delete [] r1bwave;
    delete [] sum;

    delete [] r1bm12[0][0];
    delete [] r1bm13[0][0];
  }

  // SWIR wavelengths/bandpass
  start.clear();
  start.push_back(0);
  count.clear();
  count.push_back(NIWAVE);

  var = outfile.ncGrps[0].getVar( "SWIR_wavelength");
  var.putVar( start, count, swave);
  var = outfile.ncGrps[0].getVar( "SWIR_bandpass");
  var.putVar( start, count, spass);

  var = outfile.ncGrps[0].getVar( "SWIR_solar_irradiance"); 
  var.putVar( start, count, sf0);

  start.clear();
  start.push_back(0); start.push_back(0); start.push_back(0);

  count.clear();
  count.push_back(9); count.push_back(2); count.push_back(3);
    
  var = outfile.ncGrps[0].getVar( "SWIR_m12_coef");
  var.putVar( start, count, &swir_lut.m12_coef[0][0][0]);
  var = outfile.ncGrps[0].getVar( "SWIR_m13_coef");
  var.putVar( start, count, &swir_lut.m13_coef[0][0][0]);

  
  //  l1b_filename.assign(argv[4]);
  outfile.write_granule_metadata( tstart, tend, l1b_filename);

  // write global attributes, including history and date_created
  set_global_attrs(outfile.l1bfile, history, doi, pversion);
  outfile.close();
  
  delete [] sstime;
  delete [] spin;
  delete [] hside;
  delete [] tfl;
  delete [] dtype;
  delete [] lines;
  delete [] iagg;
  delete [] bagg;
  delete [] ragg;
  delete [] mspin;
  delete [] ot_10us;
  delete [] enc_count;
  delete [] sgmat;
  delete [] thetap;
  delete [] thetas;
  delete [] ia;
  delete [] evtime;
  
  delete [] hamenc[0];
  delete [] hamenc;
  delete [] rtaenc[0];
  delete [] rtaenc;
  delete [] caltemps[0];
  delete [] caltemps;
  delete [] bdn[0];
  delete [] bdn;
  delete [] rdn[0];
  delete [] rdn;
  delete [] sdn[0];
  delete [] sdn;  
  delete [] bcal[0];
  delete [] bcal;
  delete [] rcal[0];
  delete [] rcal;
  delete [] scal[0];
  delete [] scal;  
  delete [] pview[0];
  delete [] pview;
  delete [] sview[0];
  delete [] sview;  

  delete [] bsci[0];
  delete [] bsci;
  delete [] rsci[0];
  delete [] rsci;
  delete [] ssci[0];
  delete [] ssci;

  delete [] bcalb[0];
  delete [] bcalb;
  delete [] rcalb[0];
  delete [] rcalb;

  delete [] blue_lut.K1[0];
  delete [] blue_lut.K1;
  delete [] blue_lut.K5_coef[0];
  delete [] blue_lut.K5_coef;

  delete [] K3T;
  
  if (bgains.K1K2 != NULL) delete [] bgains.K1K2[0];
  if (bgains.K1K2 != NULL) delete [] bgains.K1K2;
  if (bgains.K5_coef != NULL) delete [] bgains.K5_coef[0];
  if (bgains.K5_coef != NULL) delete [] bgains.K5_coef;

  // delete [] arrT[0][0]; delete [] arrT[0]; delete [] arrT;

  delete [] red_lut.K1[0];
  delete [] red_lut.K1;
  delete [] red_lut.K5_coef[0];
  delete [] red_lut.K5_coef;

  if (rgains.K1K2 != NULL) delete [] rgains.K1K2[0];
  if (rgains.K1K2 != NULL) delete [] rgains.K1K2;
  if (rgains.K5_coef != NULL) delete [] rgains.K5_coef[0];
  if (rgains.K5_coef != NULL) delete [] rgains.K5_coef;
  
  delete [] swir_lut.K1[0];
  delete [] swir_lut.K1;
  delete [] swir_lut.K5_coef[0];
  delete [] swir_lut.K5_coef;

  delete [] sgains.K1K2[0];
  delete [] sgains.K1K2;
  delete [] sgains.K5_coef[0];
  delete [] sgains.K5_coef;

  // Add delete for dark arrays
  delete [] bdark[0][0];

  clo_deleteList(optionList);
  
  return 0;
}


int read_oci_cal_lut( NcFile *calLUTfile, string tag, NcGroup gidLUT,
                      uint32_t& banddim, uint32_t mcedim, uint32_t& nldim,
                      uint32_t& poldim, cal_lut_struct& cal_lut) {

  size_t dims[3],dims4[4];
  NcDim ncDIM;
  uint32_t timedim, tempdim, tcdim, rvsdim, msdim=2;
  string bandname;
  bandname.assign( tag);
  bandname.append( "_bands");
  
  ncDIM = calLUTfile->getDim(bandname.c_str());
  banddim = ncDIM.getSize();
  
  ncDIM = calLUTfile->getDim("number_of_times");
  timedim = ncDIM.getSize();

  if (tag.compare("blue") == 0 || tag.compare("red") == 0)
    ncDIM = calLUTfile->getDim("number_of_CCD_temperatures");
  else
    ncDIM = calLUTfile->getDim("number_of_SWIR_temperatures");
  tempdim = ncDIM.getSize();

  ncDIM = calLUTfile->getDim("number_of_T_coefficients");
  tcdim = ncDIM.getSize();
  ncDIM = calLUTfile->getDim("number_of_RVS_coefficients");
  rvsdim = ncDIM.getSize();
  ncDIM = calLUTfile->getDim("number_of_nonlinearity_coefficients");
  nldim = ncDIM.getSize();
  ncDIM = calLUTfile->getDim("number_of_polarization_coefficients");
  poldim = ncDIM.getSize();
  
  float**K1 = new float *[banddim];
  K1[0] = new float[msdim*banddim];
  for (size_t i=1; i<banddim; i++) K1[i] = K1[i-1] + msdim;
  cal_lut.K1 = K1;

  dims[0] = banddim; dims[1] = msdim; dims[2] = timedim;
  float ***K2 = make3dT<float>(dims);
  cal_lut.K2 = K2;
  dims[0] = banddim; dims[1] = tempdim; dims[2] = tcdim;
  float ***K3_coef = make3dT<float>(dims);
  cal_lut.K3_coef = K3_coef;
  dims4[0] = banddim; dims4[1] = msdim; dims4[2] = mcedim; dims4[3] = rvsdim;
  float ****K4_coef = make4dT<float>(dims4);
  cal_lut.K4_coef = K4_coef;

  double **K5_coef = new double *[banddim];
  K5_coef[0] = new double[nldim*banddim];
  for (size_t i=1; i<banddim; i++) K5_coef[i] = K5_coef[i-1] + nldim;
  cal_lut.K5_coef = K5_coef;

  uint32_t *sat_thres = new uint32_t [banddim];
  cal_lut.sat_thres = sat_thres;

  dims[0] = banddim; dims[1] = msdim; dims[2] = poldim;
  float ***m12_coef = make3dT<float>(dims);
  cal_lut.m12_coef = m12_coef;
  float ***m13_coef = make3dT<float>(dims);
  cal_lut.m13_coef = m13_coef;

  // Assign cal lut dimensions array
  cal_lut.ldims[0] = timedim; cal_lut.ldims[1] = tempdim;
  cal_lut.ldims[2] = tcdim; cal_lut.ldims[3] = rvsdim;
  cal_lut.ldims[4] = nldim; cal_lut.ldims[5] = msdim;
  cal_lut.ldims[6] = poldim;

  NcVar var;
  var = gidLUT.getVar( "K1");
  var.getVar( &cal_lut.K1[0][0]);
  var = gidLUT.getVar( "K2");
  var.getVar( &cal_lut.K2[0][0][0]);
  var = gidLUT.getVar( "K3_coef");
  var.getVar( &cal_lut.K3_coef[0][0][0]);
  var = gidLUT.getVar( "K4_coef");
  var.getVar( &cal_lut.K4_coef[0][0][0][0]);
  var = gidLUT.getVar( "K5_coef");
  var.getVar( &cal_lut.K5_coef[0][0]);
  var = gidLUT.getVar( "sat_thres");
  var.getVar( &cal_lut.sat_thres[0]);

  var = gidLUT.getVar( "m12_coef");
  var.getVar( &cal_lut.m12_coef[0][0][0]);
  var = gidLUT.getVar( "m13_coef");
  var.getVar( &cal_lut.m13_coef[0][0][0]);

  return 0;
}

// float  K1K2[nib][msdim]: absolute and temporal gain
// float  K3_coef[nib][tempdim][tcdim]: temperature correction
// float  K4_coef[nib][msdim][mcedim][rvsdim]: RVS correction
// double K5_coef[nib][nldim]: linearity correction

int make_oci_gains( uint32_t nib, uint32_t banddim, uint16_t iyr, uint32_t jd,
                    double stime, size_t numTimes, double *K2t, int16_t board_id,
                    int16_t iagg, int16_t *jagg, cal_lut_struct& cal_lut,
                    float **gmat, gains_struct& gains) {

  for (size_t i=0; i<6; i++) gains.ldims[i] = cal_lut.ldims[i];

  uint16_t timedim = gains.ldims[0];
  uint16_t tempdim = gains.ldims[1];
  uint16_t tcdim = gains.ldims[2];
  uint16_t rvsdim = gains.ldims[3];
  uint16_t nldim = gains.ldims[4];
  uint16_t msdim = gains.ldims[5];

  bool hyper = false;
  uint16_t bd_id;
  int16_t *iaf = NULL;

  // Hyperspectral bands
  if ( board_id == -1) {
    hyper = true;
    iaf = new int16_t[nib];
    int ib = 0;
    for (size_t i=0; i<16; i++) {
      if ( jagg[i] > 0) {
        uint32_t nb = 32/jagg[i];
        for (size_t j=0; j<nb; j++) {
          if (iagg*jagg[i] < 4)
            iaf[ib+j] = 4/(iagg*jagg[i]);
          else
            iaf[ib+j] = 4/4;
        }
        ib += nb;
      }
    }
  } else bd_id = board_id % 2;
  
  size_t dims[3];
  
  float **K1K2 = new float *[nib];
  K1K2[0] = new float[nib*msdim];
  for (size_t i=1; i<nib; i++) K1K2[i] = K1K2[i-1] + msdim;
  gains.K1K2 = K1K2;
  dims[0] = nib; dims[1] = tempdim; dims[2] = tcdim;
  float ***K3_coef = make3dT<float>(dims);
  gains.K3_coef = K3_coef;
  dims[0] = nib; dims[1] = msdim; dims[2] = rvsdim;
  float ***K4_coef = make3dT<float>(dims);
  gains.K4_coef = K4_coef;

  double **K5_coef = new double *[nib];
  K5_coef[0] = new double[nib*nldim];
  for (size_t i=1; i<nib; i++) K5_coef[i] = K5_coef[i-1] + nldim;
  gains.K5_coef = K5_coef;

  uint32_t *sat_thres = new uint32_t[nib];
  gains.sat_thres = sat_thres;
  
  // Mirror-side dependent gains
  double *K2 = new double[banddim];
  for (size_t ms=0; ms<msdim; ms++) {

    // Get temporal gain and combine with absolute gain
    double d2 = jd - 2451545 + stime/86400.0;

    size_t kd=0;
    for (size_t j=numTimes-1; j>=0; j--) {
      if ( d2 > K2t[j]) {
        kd = j;
        break;
      }
    }
    if ( kd < (size_t) (timedim-1)) {
      double ff = (d2 - K2t[kd]) / (K2t[kd+1] - K2t[kd]);
      for (size_t j=0; j<banddim; j++)
        K2[j] = cal_lut.K2[j][ms][kd]*(1.0-ff) + cal_lut.K2[j][ms][kd+1]*ff;
    } else {
      for (size_t j=0; j<banddim; j++) K2[j] = cal_lut.K2[j][ms][kd];
    }

    int16_t iaf0 = 1;
    for (size_t i=0; i<nib; i++) {
      gains.K1K2[i][ms] = 0;
      if (hyper) iaf0 = iaf[i];
      for (size_t j=0; j<banddim; j++)
        gains.K1K2[i][ms] += gmat[j][i] * cal_lut.K1[j][ms]*K2[j]*iaf0;
    }

    // Generate RVS coefficents
    for (size_t i=0; i<nib; i++) {
      for (size_t k=0; k<rvsdim; k++) {
        gains.K4_coef[i][ms][k] = 0;
        for (size_t j=0; j<banddim; j++) {
          if (hyper)
            gains.K4_coef[i][ms][k] +=
              gmat[j][i] * cal_lut.K4_coef[j][ms][0][k];
          else
            gains.K4_coef[i][ms][k] +=
              gmat[j][i] * cal_lut.K4_coef[j][ms][bd_id][k];
        }
      }
    }
  }
  delete [] K2;
  
  // Generate temperature coefficients
  for (size_t i=0; i<nib; i++) {
    for (size_t k=0; k<tcdim; k++) {
      for (size_t l=0; l<tempdim; l++) {
        gains.K3_coef[i][l][k] = 0;
        for (size_t j=0; j<banddim; j++)
          gains.K3_coef[i][l][k] += gmat[j][i] * cal_lut.K3_coef[j][l][k];
      }
    }
  }

  // Generate linearity coefficients
  for (size_t i=0; i<nib; i++) {
    for (size_t k=0; k<nldim; k++) {
      gains.K5_coef[i][k] = 0;
      for (size_t j=0; j<banddim; j++) {
        if (hyper)
          gains.K5_coef[i][k] +=
            gmat[j][i] * cal_lut.K5_coef[j][k] * powf(iaf[i], (float) i);
        else
          gains.K5_coef[i][k] +=
            gmat[j][i] * cal_lut.K5_coef[j][k];
      }
    }
  }

  // Generate saturation thresholds
  for (size_t i=0; i<nib; i++) {
    gains.sat_thres[i] = 0;
    for (size_t j=0; j<banddim; j++) {
      if (hyper)
        gains.sat_thres[i] += gmat[j][i] * cal_lut.sat_thres[j] / iaf[i];
      else
        gains.sat_thres[i] += gmat[j][i] * cal_lut.sat_thres[j];
    }
  }

  if (hyper) delete[] iaf;
  
  return 0;
}


template <typename T>
int get_oci_dark( size_t iscn, uint32_t nscan, uint8_t *hside, uint16_t ndsc,
                  uint16_t nskp, int16_t iags, int16_t iagd, uint32_t ntaps,
                  int16_t *jagg, uint32_t dfill, int16_t ndc, T ***dark,
                  uint32_t nib, float *dc, int16_t& iret) {

  // Program to generate dark corrections for OCI data by averaging the
  // dark collect data and correcting for bit shift/truncation if necessary

  // Determine number of bands per tap for hyperspectral data
  int16_t *nbndt = new int16_t[ntaps];

  if (ntaps == 16) {
    // hyperspectral bands
    for (size_t i=0; i<ntaps; i++)
      if ( jagg[i] > 0) nbndt[i] = 32 / jagg[i]; else nbndt[i] = 0;
  } else {
    for (size_t i=0; i<ntaps; i++) nbndt[i] = 9;
  }
  int16_t nbnd = 0;
  for (size_t i=0; i<ntaps; i++) nbnd += nbndt[i];

  // Select data for HAM side and determine scan indices
  int32_t *kh = new int32_t[nscan];
  int32_t nkh=0;
  for (size_t i=0; i<nscan; i++) {
    if ( hside[i] == hside[iscn]) {
      kh[i] = (int32_t) i;
      nkh++;
    } else {
      kh[i] = -1;
    }
  }

  int32_t js=0;
  for (size_t i=0; i<nscan; i++) {
    if ( kh[i] == (int32_t) iscn) {
      js = (int32_t) i;
      break;
    }
  }

  // Check for valid dark collect data within specified range
  uint16_t ndscl = ndsc;
  bool valid_dark_found = false;
  int32_t is1=js, is2=js;

  while (!valid_dark_found && ndscl <= nkh) {
    if ( ndsc > 1) {
      is1 = js - ndsc/2;
      is2 = js + ndsc/2;
      // Check for start or end of granule
      if (is1 < 0) {
        is1 = 0;
        is2 = ndsc - 1;
      }
      if (is2 >= nkh) {
        is1 = nkh - ndsc;
        is2 = nkh - 1;
      }
    }

    // If no valid dark data, expand scan range
    for (size_t i=is1; i<=(size_t) is2; i++) {
      for (size_t j=nskp; j<(size_t) ndc; j++) {
        if ( dark[kh[i]][0][j] != dfill) {
          valid_dark_found = true;
          break;
        }
      }
    }
    if ( !valid_dark_found) ndscl += 2;
  }

  if ( !valid_dark_found) {
    iret =-1;
    return 0;
  }
  
  // Loop through taps and compute dark correction
  int16_t ibnd=0;
  for (size_t i=0; i<ntaps; i++) {
    if ( jagg[i] > 0) {
      float ddiv = 1.0;
      float doff = 0.0;
      if ( iags*jagg[i] > 4) {
        ddiv = iagd * jagg[i] / 4.0;
        doff = (ddiv-1) / (2*ddiv);
      }

      for (size_t j=0; j<(size_t) nbndt[i]; j++) {
        float sum = 0.0;
        int nv = 0;
        for (size_t k=kh[is1]; k<=(size_t) kh[is2]; k++) {
          for (size_t l=nskp; l<(size_t) ndc; l++) {
            if (dark[k][ibnd+j][l] != dfill) {
              sum += dark[k][ibnd+j][l];
              nv++;
            }
          }
        }
        dc[ibnd+j] = ((sum/(nv))/ddiv) - doff;
      } // j loop
    } // if ( 
    ibnd += nbndt[i];
  } // i loop

  delete [] kh;
  delete [] nbndt;

  iret = 0;
  if (ndscl > ndsc) iret = 1;
  
  return 0;
}


int get_oci_temp_corr( uint32_t nib, gains_struct gains, float *K3T,
                       float *caltemps, uint32_t nscan, float *k3) {

  uint16_t tempdim = gains.ldims[1];
  uint16_t tcdim = gains.ldims[2];

  for (size_t i=0; i<nib; i++) k3[i] = 1.0;

  for (size_t i=0; i<tempdim; i++) {
    float td = caltemps[i] - K3T[i];
      for (size_t j=0; j<tcdim; j++) {
        for (size_t k=0; k<nib; k++) {
          k3[k] -= gains.K3_coef[k][i][j] * powf(td, j+1);
        }
      }
  }

  return 0;
}

int get_oci_rvs_corr( uint32_t nib, uint16_t pdim, uint8_t hside,
                      gains_struct gains, double *theta, float **k4) {

  // Program to compute RVS correction from coefficients and scan angle

  uint16_t rvsdim = gains.ldims[3];

  for (size_t i=0; i<nib; i++)
    for (size_t j=0; j<pdim; j++)
      k4[i][j] = 1.0;

  for (size_t i=0; i<rvsdim; i++)
    for (size_t j=0; j<nib; j++)
      for (size_t k=0; k<pdim; k++)
        k4[j][k] += gains.K4_coef[j][hside][i] * powf(theta[k], i+1);
  
  return 0;
}

int get_oci_lin_corr( uint32_t nib, uint16_t pdim, uint32_t nldim,
                      gains_struct gains, float **dn, float **k5) {

  for (size_t i=0; i<pdim; i++) {
    // Zeroth-order correction
    for (size_t j=0; j<nib; j++) {
      k5[j][i] = gains.K5_coef[j][0];

      for (size_t k=1; k<nldim; k++) {
        k5[j][i] += gains.K5_coef[j][k]*powf(dn[j][i], k);
      }
    }
  }

  return 0;
}


int get_oci_cal_temps( NcFile *l1afile, NcGroup egid,
                       uint16_t ntemps, uint32_t nscan, double *evtime,
                       float **caltemps) {

  // Program to read temperatures used in calibration from L1A file and
  // interpolate to scan times
  // The order of temperatures is:
  //   Lens housings: blue CCD side, blue grating side,
  //                  red CCD side,  red grating side
  //   CCDs: red right,  red left,
  //         blue right, blue left
  //   SDA detector 1 - 16
  //   AOB 1, 2, 3, 4, 7, 8
  //   MOSB near MLA

  NcDim ncDIM;
  uint32_t tlmdim, daudim, icdudim;
  ncDIM = l1afile->getDim("tlm_packets");
  tlmdim = ncDIM.getSize();
  ncDIM = l1afile->getDim("DAUC_temps");
  daudim = ncDIM.getSize();
  ncDIM = l1afile->getDim("ICDU_therm");
  icdudim = ncDIM.getSize();

  double *dauctime = new double[tlmdim];
  
  float **dauctemp = new float *[tlmdim];
  dauctemp[0] = new float[daudim*tlmdim];
  for (size_t i=1; i<tlmdim; i++) dauctemp[i] = dauctemp[i-1] + daudim;

  double *icdutime = new double[tlmdim];
  
  float **icdutherm = new float *[tlmdim];
  icdutherm[0] = new float[icdudim*tlmdim];
  for (size_t i=1; i<tlmdim; i++) icdutherm[i] = icdutherm[i-1] + icdudim;

  NcVar var;
  var = egid.getVar( "DAUC_temp_time");
  var.getVar( dauctime);
  var = egid.getVar( "DAUC_temperatures");
  var.getVar( &dauctemp[0][0]);
  var = egid.getVar( "TC_tlm_time");
  var.getVar( icdutime);
  var = egid.getVar( "ICDU_thermisters");
  var.getVar( &icdutherm[0][0]);

  /*
DAUCTEMP        FLOAT     = Array[69, 301]
DAUCTIME        DOUBLE    = Array[301]
ICDUTHERM       FLOAT     = Array[74, 301]
ICDUTIME        DOUBLE    = Array[301]

ICDU_therm = 74 ;
  */

  // Indices of required temperatures
  // MLA and lens housings (blue CCD, blue grating, red CCD, red grating)
  uint32_t iicdu[5] = {23,24,25,26,11};
  // Red and blue CCDs, SDA detectors, AOB 1, 2, 3, 4, 7, 8
  uint32_t idauc[26] = {5,6,12,13,
                    14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,
                    30,31,32,33,36,37};

  // Loop through temperatures
  double *dstime = new double[nscan];
  for (size_t i=0; i<nscan; i++) dstime[i] = evtime[i] - evtime[0];

  // ICDU thermistors
  double *ditime = new double[tlmdim];
  double *dummy = new double[tlmdim];
  double c0, c1, cov00, cov01, cov11, sumsq;
  for (size_t i=0; i<=3; i++) {
    uint32_t k=0;
    for (size_t j=0; j<tlmdim; j++)
      if ( icdutime[j] > 0) {
        ditime[k] = icdutime[j] - evtime[0];
        dummy[k++] = icdutherm[j][iicdu[i]];
      }
    
    gsl_fit_linear(ditime, 1, dummy, 1, k,
                   &c0, &c1, &cov00, &cov01, &cov11, &sumsq);

    for (size_t j=0; j<nscan; j++) caltemps[j][i] = c0 + c1*dstime[j];
  }

  size_t k=0;
  for (size_t j=0; j<tlmdim; j++)
    if ( icdutime[j] > 0)
      dummy[k++] = icdutherm[j][iicdu[4]];

  gsl_fit_linear(ditime, 1, dummy, 1, k,
                 &c0, &c1, &cov00, &cov01, &cov11, &sumsq);
  for (size_t j=0; j<nscan; j++) caltemps[j][30] = c0 + c1*dstime[j];

  // DAUC temperatures
  double *ddtime = new double[tlmdim];
  for (size_t i=0; i<=25; i++) {
    size_t k=0;
    for (size_t j=0; j<tlmdim; j++)
      if ( dauctime[j] > 0) {
        ddtime[k] = dauctime[j] - evtime[0];
        dummy[k++] = dauctemp[j][idauc[i]];
      }
    
    gsl_fit_linear(ddtime, 1, dummy, 1, k,
                   &c0, &c1, &cov00, &cov01, &cov11, &sumsq);

    for (size_t j=0; j<nscan; j++) caltemps[j][i+4] = c0 + c1*dstime[j];
  }
    
  delete [] dauctime;
  delete [] dauctemp[0];
  delete [] dauctemp;

  delete [] icdutime;
  delete [] icdutherm[0];
  delete [] icdutherm;

  delete [] ditime;
  delete [] ddtime;
  
  return 0;
}


int l1bFile::write_granule_metadata( std::string tstart, std::string tend,
                                     std::string l1b_name) {

  l1bfile->putAtt("time_coverage_start", tstart.c_str());
  l1bfile->putAtt("time_coverage_end", tend.c_str());

  // Write product file name
  l1bfile->putAtt("product_name", l1b_name.c_str());

  return 0;
}


int l1bFile::close() {

  try { 
    l1bfile->close();
  }
  catch ( NcException& e) {
    cout << e.what() << endl;
    cerr << "Failure closing: " + fileName << endl;
    exit(1);
  }
  
  return 0;
}


template <typename T>
T*** make3dT( size_t dims[3]) {

  T ***arr3d = new T **[dims[0]];

  arr3d[0] = new T*[dims[0]*dims[1]];
  arr3d[0][0] = new T[dims[0]*dims[1]*dims[2]];

  for (size_t i=1; i<dims[0]; i++) arr3d[i] = arr3d[i-1] + dims[1];

  for (size_t i=0; i<dims[0]; i++) {
    if ( i > 0) arr3d[i][0] = arr3d[i-1][0] + dims[1]*dims[2];
    for (size_t j=1; j<dims[1]; j++)
      arr3d[i][j] = arr3d[i][j-1] + dims[2];
  }

  return arr3d;
}

template <typename T>
T**** make4dT( size_t dims[4]) {

  T ****arr4d = new T ***[dims[0]];

  arr4d[0] = new T**[dims[0]*dims[1]];
  arr4d[0][0] = new T*[dims[0]*dims[1]*dims[2]];
  arr4d[0][0][0] = new T[dims[0]*dims[1]*dims[2]*dims[3]];
   
  for (size_t i=0; i<dims[0]; i++) {
    if ( i > 0) {
      arr4d[i] = arr4d[i-1] + dims[1];
      arr4d[i][0] = arr4d[i-1][0] + dims[1]*dims[2];
      arr4d[i][0][0] = arr4d[i-1][0][0] + dims[1]*dims[2]*dims[3];
    }
    for (size_t j=0; j<dims[1]; j++) {
      if ( j > 0) {
        arr4d[i][j] = arr4d[i][j-1] + dims[2];
        arr4d[i][j][0] = arr4d[i][j-1][0] + dims[2]*dims[3];
      }
      for (size_t k=1; k<dims[2]; k++) {
        arr4d[i][j][k] = arr4d[i][j][k-1] + dims[3];
      }
    }
  }

  return arr4d;
}
