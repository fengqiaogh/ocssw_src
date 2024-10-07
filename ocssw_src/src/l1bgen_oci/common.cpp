#include "common.h"

using namespace std;
using namespace netCDF;
using namespace netCDF::exceptions;
#include "nc4utils.h"

int read_mce_tlm( NcFile *l1afile, geo_struct& geo_lut, NcGroup egid,
                  uint32_t nmcescan, uint32_t nenc, int32_t& ppr_off,
                  double& revpsec, double&secpline,
                  int16_t& board_id, int32_t *mspin, int32_t *ot_10us,
                  uint8_t *enc_count, float **hamenc, float **rtaenc,
                  int16_t &iret) {

  iret = 0;
  
  NcDim mce_dim = l1afile->getDim("MCE_block");
  uint32_t mce_blk = mce_dim.getSize();
  NcDim ddc_dim = l1afile->getDim("DDC_tlm");
  uint32_t nddc = ddc_dim.getSize();
  NcDim tlm_dim = l1afile->getDim("tlm_packets");
  uint32_t ntlmpack = tlm_dim.getSize();
  
  uint8_t **mtlm = new uint8_t *[nmcescan];
  mtlm[0] = new uint8_t[mce_blk*nmcescan];
  for (size_t i=1; i<nmcescan; i++) mtlm[i] = mtlm[i-1] + mce_blk;

  int16_t **enc = new int16_t *[nmcescan];
  enc[0] = new int16_t[4*nenc*nmcescan];
  for (size_t i=1; i<nmcescan; i++) enc[i] = enc[i-1] + 4*nenc;
  
  uint8_t **ddctlm = new uint8_t *[ntlmpack];
  ddctlm[0] = new uint8_t[nddc*ntlmpack];
  for (size_t i=1; i<ntlmpack; i++) ddctlm[i] = ddctlm[i-1] + nddc;

  NcVar var;
  vector<size_t> start, count;
  start.push_back(0);
  start.push_back(0);
  count.push_back(1);

  var = egid.getVar( "MCE_telemetry");
  var.getVar( &mtlm[0][0]);
  
  var = egid.getVar( "MCE_encoder_data");
  count.pop_back();
  count.push_back(4*nenc);
  var.getVar( &enc[0][0]);

  var = egid.getVar( "MCE_spin_ID");
  var.getVar( mspin);

  var = egid.getVar( "DDC_telemetry");
  var.getVar( &ddctlm[0][0]);

  // Check for missing MCE telemetry
  bool noMCE = true;
  for (size_t i=0; i<nmcescan; i++) {
    if (mspin[i] >= 0) {
      noMCE = false;
      break;
    }
  }

  if (noMCE == true) {
    cout << "No MCE telemetry in L1A file" << endl;
    exit(1);
  }
  
  int32_t max_enc_cts = 131072; // 2^17
  double clock[2];
  clock[0] = geo_lut.master_clock;
  clock[1] = geo_lut.MCE_clock;

  // Get ref_pulse_divider and compute commanded rotation rate
  uint32_t ui32;
  uint32_t ref_pulse_div[2];
  memcpy( &ui32, &mtlm[0][0], 4);
  ref_pulse_div[0] = SWAP_4( ui32) % 16777216; // 2^24
  memcpy( &ui32, &mtlm[0][4], 4);
  ref_pulse_div[1] = SWAP_4( ui32) % 16777216; // 2^24

  int32_t ref_pulse_sel = mtlm[0][9] / 128;
  
  revpsec = clock[ref_pulse_sel] / 2 / max_enc_cts /
    (ref_pulse_div[ref_pulse_sel]/256.0 + 1);

  // Check for static or reverse scan  
  int32_t avg_step_spd = mtlm[0][428]*256+mtlm[0][429];
  if (abs(avg_step_spd) < 1000) {
    // Static mode
    iret = 2;
    revpsec = 0.0;
    cout << "OCI static mode" << endl;
  } else if (avg_step_spd < 0) {
    // Reverse spin
    revpsec = -revpsec;
    iret = 3;
    cout << "OCI reverse scan" << endl;
  }
                   
  // Get PPR offset and on-time_10us
  memcpy( &ui32, &mtlm[0][8], 4);
  ppr_off = SWAP_4( ui32) % max_enc_cts;
  for (size_t i=0; i<nmcescan; i++) {
    memcpy( &ui32, &mtlm[i][368], 4);
    ot_10us[i] = SWAP_4( ui32) % 4;
  }

  // Get MCE board ID  
  board_id = (int16_t) mtlm[0][322] / 16;

  // Get TDI time and compute time increment per line
  uint16_t ui16;
  memcpy( &ui16, &ddctlm[0][346], 2);
  int32_t tditime = SWAP_2( ui16);
  secpline = (tditime+1) / clock[0]; 

  // Get valid encoder count, HAM and RTA encoder data
  for (size_t i=0; i<nmcescan; i++) enc_count[i] = mtlm[i][475];
  //  float enc_s = 81.0 / 2560;
  for (size_t i=0; i<nmcescan; i++) {
    for (size_t j=0; j<nenc; j++) {
      hamenc[i][j] = enc[i][4*j+0] * geo_lut.ham_enc_scale;
      rtaenc[i][j] = enc[i][4*j+1] * geo_lut.rta_enc_scale;
    }
  }

  delete[] mtlm[0];
  delete[] mtlm;

  delete [] enc[0];
  delete [] enc;

  delete [] ddctlm[0];
  delete [] ddctlm;

  return 0;
}


int get_ev( double secpline, int16_t *dtype, int16_t *lines, int16_t *iagg,
            uint16_t& pcdim, uint16_t& psdim, double& ev_toff,
            float *clines, float *slines, double *deltc, double *delts,
            bool dark, int16_t &iret) {

  // Find end of no-data zone
  int16_t iz=0, line0=0;
  iret = -1;
  while ( dtype[iz] == 0) {
    line0 += lines[iz];
    iz++;
  }
  if (iz == 10) return 0;
  
  // Find number of pixels in Earth views
  pcdim = 0;
  psdim = 0;
  int16_t linen = line0;

  uint8_t ltype[] = {0,1,0,1,1,1,1,1,1,1,0,0,0};
  if (dark) ltype[2] = 1;
  
  for (size_t i=iz; i<9; i++) {
    // Check for not dark or no-data
    if ( ltype[dtype[i]] == 1) { 
      //    if ( dtype[i] != 0 && dtype[i] != 2 && dtype[i] < 10) {

      uint16_t np = lines[i] / iagg[i];
      for (size_t j=0; j<np; j++) {
        clines[pcdim+j] = linen + j*iagg[i] + 0.5*iagg[i] - 64;
      }
      pcdim += np;
      uint16_t ns = lines[i] / 8;
      for (size_t j=0; j<ns; j++) {
        slines[psdim+j] = linen + j*8 + 4 - 64;
      }
      psdim += ns;
      iret = 0;
    }
    linen += lines[i];
  }

  // Calculate times
  for (size_t i=0; i<(size_t) pcdim; i++) deltc[i] = secpline * clines[i];
  ev_toff = 0.5 * (deltc[0] + deltc[pcdim-1]);
  for (size_t i=0; i<(size_t) pcdim; i++) deltc[i] -= ev_toff;
    
  for (size_t i=0; i<(size_t) psdim; i++)
    delts[i] = secpline * slines[i] - ev_toff;

  return 0;
}

int get_oci_vecs( uint32_t nscan, uint16_t pdim, double as_planarity[5],
                  double at_planarity[5], int32_t *rta_nadir,
                  double ham_ct_angles[2], double ev_toff, int32_t spin,
                  uint8_t hside, float *clines,
                  double *delt, double revpsec, int32_t ppr_off,
                  int16_t board_id, uint32_t nmcescan, int32_t *mspin,
                  uint8_t *enc_count, float **hamenc, float **rtaenc,
                  float **pview, double *theta, int16_t& iret) {

  // This program generates the OCI Earth view vectors for one spin.  
  // It uses MCE telemetry and encoder data.  Further refinements will be made
  // as the instrument optics model and test results become available.  
  // Reference: "Use of OCI Telemetry to Determine Pixel Line-of-Sight",
  // F. Patt, 2020-05-18

  int32_t max_enc_cts = 131072; // 2^17
  double dtenc = 0.001;
  constexpr double pi = 3.14159265358979323846;
  int16_t bd_id = board_id % 2;
  
  double rad2asec = (180/pi) * 3600;
  iret = 0;

  // Compute scan angle corresponding to PPR
  //  float pprang = 2 * pi * (ppr_off - geoLUT.rta_nadir[bd_id]) / max_enc_cts;
  float pprang = 2 * pi * (ppr_off - rta_nadir[bd_id]) / max_enc_cts;
  if (pprang > pi) pprang -= 2*pi;
  
  // Compute ideal scan angles for science pixels
  double *toff = new double[pdim];
  for (size_t i=0; i<pdim; i++) toff[i] = delt[i] + ev_toff;

  for (size_t i=0; i<pdim; i++)
    theta[i] = pprang + 2 * pi * revpsec * toff[i];
  // Interpolate encoder data to pixel times and add to scan angles
  // RTA only for now, include HAM when we know how.

  double *thetacor = new double[pdim];

  int isp = -1;
  for (size_t i=0; i<nmcescan; i++) {
    if ( mspin[i] == spin) {
      isp = (int) i;
      break;
    }
  }

  // Interpolate encoder data to pixel times and add to scan angles
  // RTA only for now, include HAM when we know how.
  if ( isp == -1) {
    cout << "No MCE encoder data for spin: " << spin << endl;
    iret = 1;
  } else {
    size_t ip = 0, ke;
    double tenc_ke;
    while ( ip < pdim) {
      // Encoder sample times at 1 KHz
      // ke = where(tenc le toff[ip] AND tenc[1:*] gt toff[ip])
      bool found=false;
      for (size_t j=0; j<enc_count[isp]; j++) {
        double tenc = j*dtenc;
        if ( tenc <= toff[ip] && (tenc+dtenc) > toff[ip]) {
          ke = j;
          tenc_ke = tenc;
          found = true;
          break;
        }
      }
      
      size_t njp=0;
      if (found) {
        for (size_t i=0; i<pdim; i++) {
          if ( toff[i] >= tenc_ke && toff[i] < tenc_ke+dtenc) {
            double ft = (toff[i] - tenc_ke) / dtenc;
            thetacor[i] =
              (1-ft)*rtaenc[isp][ke] + ft*rtaenc[isp][ke+1] -
              ((1-ft)*hamenc[isp][ke] + ft*hamenc[isp][ke+1] +
               ham_ct_angles[hside])*0.236;
               
            njp++;
          }
        }
      } else {
        cout << "Encoder interpolation error" << endl;
        exit(-1);
      }
      ip += njp;
    } // while ( ip < pdim)

    // Calculate planarity deviations and view vectors
    double *ascan  = new double[pdim];
    double *atrack = new double[pdim];
    for (size_t i=0; i<pdim; i++) {
      theta[i] = theta[i] - thetacor[i] / rad2asec;
      
      ascan[i]  = as_planarity[0];
      atrack[i] = at_planarity[0];
      for (size_t k=1; k<5; k++) {
        ascan[i]  += as_planarity[k]*pow(theta[i], k);
        atrack[i] += at_planarity[k]*pow(theta[i], k);
      }

      pview[i][0] = -sin(atrack[i]/rad2asec);
      pview[i][1] = sin(theta[i] - ascan[i]/rad2asec);
      pview[i][2] = cos(theta[i] - ascan[i]/rad2asec);
    }

    delete [] ascan;
    delete [] atrack;
  }
  
  delete [] thetacor;
  delete [] toff;

  return 0;
}

int createField( NcGroup &ncGrp, const char *sname, const char *lname, 
                 const char *standard_name, const char *units,
                 const char *description,
                 void *fill_value, const char *flag_masks,
                 const char *flag_meanings,
                 double low, double high, 
                 double scale, double offset,
                 int nt, vector<NcDim>& varVec, string coordinates) {

  /* Create the NCDF dataset */
  NcVar ncVar;

  try {
    ncVar = ncGrp.addVar(sname, nt, varVec);   
  }
  catch ( NcException& e) {
    cout << e.what() << endl;
    cerr << "Failure creating variable: " << sname << endl;
    exit(1);
  }
  int var_id = ncVar.getId();
  int nc_id = ncVar.getParentGroup().getId();
  std::string grp_name = ncVar.getParentGroup().getName();
  int32_t dimids[5];
  for (size_t idim = 0 ; idim < ncVar.getDims().size(); idim ++){
      dimids[idim] = ncVar.getDims()[idim].getId();
  }
  // Set fill value
  double fill_value_dbl;
  memcpy( &fill_value_dbl, fill_value, sizeof(double));

  int8_t i8;
  uint8_t ui8;
  int16_t i16;
  uint16_t ui16;
  int32_t i32;
  uint32_t ui32;
  float f32;

  if ( low != fill_value_dbl) {
    if ( nt == NC_BYTE) {
      i8 = fill_value_dbl;
      ncVar.setFill(true, (void *) &i8);
    } else if ( nt == NC_UBYTE) {
      ui8 = fill_value_dbl;
      ncVar.setFill(true, (void *) &ui8);
    } else if ( nt == NC_SHORT) {
      i16 = fill_value_dbl;
      ncVar.setFill(true, (void *) &i16);
    } else if ( nt == NC_USHORT) {
      ui16 = fill_value_dbl;
      ncVar.setFill(true, (void *) &ui16);
    } else if ( nt == NC_INT) {
      i32 = fill_value_dbl;
      ncVar.setFill(true, (void *) &i32);
    } else if ( nt == NC_UINT) {
      ui32 = fill_value_dbl;
      ncVar.setFill(true, (void *) &ui32);
    } else if ( nt == NC_FLOAT) {
      f32 = fill_value_dbl;
      ncVar.setFill(true, (void *) &f32);
    } else {
      ncVar.setFill(true, (void *) &fill_value_dbl);
    }
  }

  /* vary chunck size based on dimensions */ 
  int deflate_level = 5;
  std::vector<size_t> chunkVec;
  bool set_compression = (grp_name == "geolocation_data") || (grp_name == "observation_data");
  if (set_compression) {
    if(varVec.size() == 2)
      chunkVec = {512, 1272};  // 256 lines, all pixels(1272)
    if(varVec.size() == 3)
      chunkVec = {40, 16, 1272};  // 40 bands, 16 lines, all pixels(1272)
    nc_init_compress(nc_id, var_id, dimids, varVec.size(),
                         chunkVec.data(), deflate_level);

  }

  /* Add a "long_name" attribute */
  try {
    ncVar.putAtt("long_name", lname);
  }
  catch ( NcException& e) {
    e.what();
    cerr << "Failure creating 'long_name' attribute: " << lname << endl;
    exit(1);
  }

  if ( strcmp( flag_masks, "") != 0) {

    size_t curPos=0;

    string fv;
    fv.assign( flag_masks);
    size_t pos = fv.find("=", curPos);
    fv = fv.substr(pos+1);

    size_t semicln = fv.find(";");
    pos = 0;

    int8_t vec[1024];
    int n = 0;
    while(pos != semicln) {
      pos = fv.find(",", curPos);
      if ( pos == string::npos) 
        pos = semicln;

      string flag_mask;
      istringstream iss(fv.substr(curPos, pos-curPos));
      iss >> skipws >> flag_mask;
      vec[n++] = atoi( flag_mask.c_str());
      curPos = pos + 1;
    }

    try {
      ncVar.putAtt("flag_masks", nt, n, vec);
    }
    catch ( NcException& e) {
      e.what();
      cerr << "Failure creating 'flag_masks' attribute: " << lname << endl;
      exit(1);
    }
  }

  /* Add a "flag_meanings" attribute if specified*/
  if ( strcmp( flag_meanings, "") != 0) {

    try {
      ncVar.putAtt("flag_meanings", flag_meanings);
    }
    catch ( NcException& e) {
      e.what();
      cerr << "Failure creating 'flag_meanings' attribute: "
           << flag_meanings << endl;
      exit(1);
    }
  }

  /* Add "valid_min/max" attributes if specified */
  if (low < high) {
    switch(nt) {              /* Use the appropriate number type */
    case NC_BYTE:
      {
	uint8_t vr[2];
	vr[0] = (uint8_t)low;
	vr[1] = (uint8_t)high;

        try {
          ncVar.putAtt("valid_min", NC_BYTE, 1, &vr[0]);
        }
        catch ( NcException& e) {
          e.what();
          cerr << "Failure creating 'valid_min' attribute: " << vr[0] << endl;
          exit(1);
        }

        try {
          ncVar.putAtt("valid_max", NC_BYTE, 1, &vr[1]);
        }
        catch ( NcException& e) {
          e.what();
          cerr << "Failure creating 'valid_max' attribute: " << vr[1] << endl;
          exit(1);
        }
      }
      break;
    case NC_UBYTE:
      {
	uint8_t vr[2];
	vr[0] = (uint8_t)low;
	vr[1] = (uint8_t)high;

        try {
          ncVar.putAtt("valid_min", NC_UBYTE, 1, &vr[0]);
        }
        catch ( NcException& e) {
          e.what();
          cerr << "Failure creating 'valid_min' attribute: " << vr[0] << endl;
          exit(1);
        }

        try {
          ncVar.putAtt("valid_max", NC_UBYTE, 1, &vr[1]);
        }
        catch ( NcException& e) {
          e.what();
          cerr << "Failure creating 'valid_max' attribute: " << vr[1] << endl;
          exit(1);
        }
      }
      break;
    case NC_SHORT:
      {
	int16_t vr[2];
	vr[0] = (int16_t)low;
	vr[1] = (int16_t)high;

        try {
          ncVar.putAtt("valid_min", NC_SHORT, 1, &vr[0]);
        }
        catch ( NcException& e) {
          e.what();
          cerr << "Failure creating 'valid_min' attribute: " << vr[0] << endl;
          exit(1);
        }

        try {
          ncVar.putAtt("valid_max", NC_SHORT, 1, &vr[1]);
        }
        catch ( NcException& e) {
          e.what();
          cerr << "Failure creating 'valid_max' attribute: " << vr[1] << endl;
          exit(1);
        }
      }
      break;
    case NC_USHORT:
      {
	uint16_t vr[2];
	vr[0] = (uint16_t)low;
	vr[1] = (uint16_t)high;

        try {
          ncVar.putAtt("valid_min", NC_USHORT, 1, &vr[0]);
        }
        catch ( NcException& e) {
          e.what();
          cerr << "Failure creating 'valid_min' attribute: " << vr[0] << endl;
          exit(1);
        }

        try {
          ncVar.putAtt("valid_max", NC_USHORT, 1, &vr[1]);
        }
        catch ( NcException& e) {
          e.what();
          cerr << "Failure creating 'valid_max' attribute: " << vr[1] << endl;
          exit(1);
        }
              }
      break;
    case NC_INT:
      {
	int32_t vr[2];
	vr[0] = (int32_t)low;
	vr[1] = (int32_t)high;

        try {
          ncVar.putAtt("valid_min", NC_INT, 1, &vr[0]);
        }
        catch ( NcException& e) {
          e.what();
          cerr << "Failure creating 'valid_min' attribute: " << vr[0] << endl;
          exit(1);
        }

        try {
          ncVar.putAtt("valid_max", NC_INT, 1, &vr[1]);
        }
        catch ( NcException& e) {
          e.what();
          cerr << "Failure creating 'valid_max' attribute: " << vr[1] << endl;
          exit(1);
        }

      }
      break;
    case NC_UINT:
      {
	uint32_t vr[2];
	vr[0] = (uint32_t)low;
	vr[1] = (uint32_t)high;

        try {
          ncVar.putAtt("valid_min", NC_UINT, 1, &vr[0]);
        }
        catch ( NcException& e) {
          e.what();
          cerr << "Failure creating 'valid_min' attribute: " << vr[0] << endl;
          exit(1);
        }

        try {
          ncVar.putAtt("valid_max", NC_UINT, 1, &vr[1]);
        }
        catch ( NcException& e) {
          e.what();
          cerr << "Failure creating 'valid_max' attribute: " << vr[1] << endl;
          exit(1);
        }
        
      }
      break;
    case NC_FLOAT:
      {
	float vr[2];
	vr[0] = (float)low;
	vr[1] = (float)high;

        try {
          ncVar.putAtt("valid_min", NC_FLOAT, 1, &vr[0]);
        }
        catch ( NcException& e) {
          e.what();
          cerr << "Failure creating 'valid_min' attribute: " << vr[0] << endl;
          exit(1);
        }

        try {
          ncVar.putAtt("valid_max", NC_FLOAT, 1, &vr[1]);
        }
        catch ( NcException& e) {
          e.what();
          cerr << "Failure creating 'valid_max' attribute: " << vr[1] << endl;
          exit(1);
        }        
      }
      break;
    case NC_DOUBLE:
      {
	double vr[2];
	vr[0] = low;
	vr[1] = high;

        try {
          ncVar.putAtt("valid_min", NC_DOUBLE, 1, &vr[0]);
        }
        catch ( NcException& e) {
          e.what();
          cerr << "Failure creating 'valid_min' attribute: " << vr[0] << endl;
          exit(1);
        }

        try {
          ncVar.putAtt("valid_max", NC_DOUBLE, 1, &vr[1]);
        }
        catch ( NcException& e) {
          e.what();
          cerr << "Failure creating 'valid_max' attribute: " << vr[1] << endl;
          exit(1);
        }
      }
      break;
    default:
      fprintf(stderr,"-E- %s line %d: ",__FILE__,__LINE__);
      fprintf(stderr,"Got unsupported number type (%d) ",nt);
      fprintf(stderr,"while trying to create NCDF variable, \"%s\", ",sname);
      return(1);
    }
  }           
    
  /* Add "scale_factor" and "add_offset" attributes if specified */
  if(scale != 1.0 || offset != 0.0) {
    try {
      ncVar.putAtt("scale_factor", NC_DOUBLE, 1, &scale);
    }
    catch ( NcException& e) {
      e.what();
      cerr << "Failure creating 'scale_factor' attribute: " << scale << endl;
      exit(1);
    }

    try {
      ncVar.putAtt("add_offset", NC_DOUBLE, 1, &offset);
    }
    catch ( NcException& e) {
      e.what();
      cerr << "Failure creating 'add_offset' attribute: " << offset << endl;
      exit(1);
    }
  }

  /* Add a "units" attribute if one is specified */
  if(units != NULL && *units != 0) {

    try {
      ncVar.putAtt("units", units);
    }
    catch ( NcException& e) {
      e.what();
      cerr << "Failure creating 'units' attribute: " << units << endl;
      exit(1);
    }
  }

  /* Add a "description" attribute if one is specified */
  if(description != NULL && *description != 0) {

    try {
      ncVar.putAtt("description", description);
    }
    catch ( NcException& e) {
      e.what();
      cerr << "Failure creating 'description' attribute: " << description << endl;
      exit(1);
    }
  }

  /* Add a "standard_name" attribute if one is specified */
  if(standard_name != NULL && *standard_name != 0) {
    try {
      ncVar.putAtt("standard_name", standard_name);
    }
    catch ( NcException& e) {
      e.what();
      cerr << "Failure creating 'standard_name' attribute: "
           << standard_name << endl;
      exit(1);
    }
  }
  
  /* Add a "coordinates" attribute if specified*/
  if ( !coordinates.empty()) {

    try {
      ncVar.putAtt("coordinates", coordinates);
    }
    catch ( NcException& e) {
      e.what();
      cerr << "Failure creating 'coordinates' attribute: "
           << coordinates << endl;
      exit(1);
    }
  }

  return 0;
}

int get_nib_nbb( uint32_t ntaps, size_t *ia, uint32_t ntb[16], int16_t jagg[16],
                 uint32_t& nib, uint32_t& nbb) {

  uint32_t iia;
  
  iia=0;
  for (size_t i=0; i<ntaps; i++) {
    if ( jagg[i] > 0) {
      ia[iia] = i;
      iia++;
    }
  }

  if ( iia == 0) {
    return 2;
  } else {
    for (size_t i=0; i<16; i++) ntb[i] = 0;
    for (size_t i=0; i<iia; i++) ntb[ia[i]] = 32 / jagg[ia[i]];
    nib = 0;
    for (size_t i=0; i<16; i++) nib += ntb[i];
  }

  // Compute number of bands for 8x aggregation with overlapping bands
  nbb = (ntb[ia[0]]*3) / 4 + 1;
  //amat(nbb,nib)
  for (size_t i=1; i<iia; i++) {
    if (jagg[ia[i]] >= jagg[ia[i]-1])
      nbb += ntb[ia[i]];
    else
      nbb += (ntb[ia[i]]*3) / 4 + ntb[ia[i]-1]/4;
  }

  return 0;
}

int get_agg_mat( size_t *ia, int16_t jagg[16], uint32_t ntb[16],
                 uint32_t nib, uint32_t nbb, float **amat, float **gmat) {
  
  for (size_t i=0; i<512; i++) {
    gmat[i] = new float[nib];
    for (size_t j=0; j<nib; j++) gmat[i][j] = 0.0;
  }
  
  // Populate gain aggregation matrix
  int16_t ii = 0;
  int16_t itt[16][2];
  for (size_t i=0; i<16; i++) {
    itt[i][0] = itt[i][1] = 0;
    if ( jagg[i] > 0) {
      for (size_t k=0; k<ntb[i]; k++) {
        size_t ic = 32*i;
        size_t kj = k*jagg[i];
        for (size_t j=0; j<(size_t) jagg[i]; j++)
          gmat[ic+kj+j][ii+k] = (1.0/jagg[i]);
      }
      itt[i][0] = ii;
      ii += ntb[i];
      itt[i][1] = ii - 1;
    }
  }

  for (size_t i=0; i<nib; i++) {
    amat[i] = new float[nbb];
    for (size_t j=0; j<nbb; j++) amat[i][j] = 0;
  }

  // First tap
  int16_t ib;
  for (size_t k=0; k<=(ntb[ia[0]]*3)/4; k++) {
    for (size_t l=0; l<=ntb[ia[0]]/4-1; l++) amat[k+l][k] = jagg[ia[0]]/8.0;
  }
  ib = (ntb[ia[0]]*3)/4 + 1;
  
  // Remaining taps
  uint16_t nr;
  for (size_t i=ia[1]; i<16; i++) {
    if (ntb[i] > 0) {
      if (ntb[i] >= ntb[i-1]) {
        // Transition resolution determined by preceding tap
      	nr = ntb[i-1]/4 - 1;
        // Remaining bands using preceding tap
      	if (nr > 0) {
          for (size_t k=0; k<nr; k++) {
	    uint16_t k1 = nr - k - 1;
            uint16_t k2 = ((k+1)*ntb[i]) / ntb[i-1] - 1;
            for (size_t j=0; j<=k1; j++)
              amat[itt[i-1][1]-k1+j][ib+k] = jagg[i-1] / 8.0;
            for (size_t j=0; j<=k2; j++)
              amat[itt[i][0]+j][ib+k] = jagg[i] / 8.0;
          }
      	  ib += nr;
        }
      } else {
        // Transition resolution determined using current tap
        nr = ntb[i]/4 - 1;
        // Remaining bands using previous tap
        if (nr > 0) {
          for (size_t k=0; k<nr; k++) {
            uint16_t k1 = ((nr-k)*ntb[i-1]) / ntb[i] - 1;
	    uint16_t k2 = k;
            for (size_t j=0; j<=k1; j++)
              amat[itt[i-1][1]-k1+j][ib+k] = jagg[i-1] / 8.0;
            for (size_t j=0; j<=k2; j++)
              amat[itt[i][0]+j][ib+k] = jagg[i] / 8.0;
          }            
   	  ib += nr;
        }
      }
      // Remaining bands using this tap
      for (size_t k=0; k<=(ntb[i]*3)/4; k++) {
        for (size_t j=0; j<ntb[i]/4; j++)
        amat[itt[i][0]+k+j][ib+k] = jagg[i] / 8.0;
      }
      ib += (ntb[i]*3) / 4 + 1;
    }
  }
                            
  return 0;
}

// sfl = array of scan qual flags 2 = missu=ing time
int check_scan_times( uint32_t nscan, double *sstime, short *sfl) {

  vector<uint32_t> kf, kv;
  for (size_t i=0; i<nscan; i++) {
    if (sstime[i] == -999 || sstime[i] == -32767)
      kf.push_back(i);
    else
      kv.push_back(i);
  }
  
  size_t nf = kf.size();
  if ( nf == 0) return 0;

  size_t nv = kv.size();
  
  // Interpolate valid times to fill missing values
  for (size_t i=0; i<nf; i++) {

    // Check for missing time before valid time
    if (kf[i] < kv[0]) { //if there are no good previous times
      sstime[kf[i]] = sstime[kv[0]] -
        (sstime[kv[1]]-sstime[kv[0]])*(kv[0]-kf[i])/(kv[1]-kv[0]);
    } else if (kf[i] > kv[nv-1]) { //if there are no good next times
      sstime[kf[i]] = sstime[kv[nv-1]] +
        (sstime[kv[nv-1]]-sstime[kv[nv-2]])*(kf[i]-kv[nv-1])/
        (kv[nv-1]-kv[nv-2]);
    } else {
      size_t iv=0;
      for (int j=nv-1; j>=0; j--) {
        if (kv[j] < kf[i]) {
          iv = j;
          break;
        }
      }
      sstime[kf[i]] = sstime[kv[iv]] +
        (sstime[kv[iv+1]]-sstime[kv[iv]])*(kf[i]-kv[iv])/(kv[iv+1]-kv[iv]);
    }
    *(sfl+kf[i]) |= 2;  // set missing time flag for all nf(ailed)
  }
  return 0;

}
