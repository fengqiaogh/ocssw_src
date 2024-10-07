#include <stdio.h>
#include <math.h>
#include <libgen.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <netcdf>

#include <cgal_interp.h>
#include <allocate2d.h>
#include <allocate3d.h>
#include <timeutils.h>
#include <clo.h>
#include <genutils.h>

#include <gsl/gsl_matrix_float.h>
#include <gsl/gsl_blas.h>

#define NCACHE 20
const double PI = 3.14159265358979323846;
const double RADEG = 180 / PI;

using namespace std;
using namespace netCDF;
using namespace netCDF::exceptions;

//    Modification history:
//  Programmer     Organization   Date     Ver   Description of change
//  ----------     ------------   ----     ---   ---------------------
//  Joel Gales     SAIC           01/14/24 1.50  Fix type bug in the
//                                               calculation of the water cloud
//                                               average values for the L1C bins
//
//  Joel Gales     SAIC           02/01/24 1.51  Make CLD processing optional
//  Martin Montes    SSAI           07/22/24 1.52  Fix water/ice cloud mask

int accum( uint32_t totlines, uint32_t npixels,
           uint32_t nalong, uint32_t nacross,
           short *brow, short *bcol,
           float **cldprod, int8_t **cld_phase,
           float **binval_ic, short **nobs_ic,
           float **binval_wc, short **nobs_wc) {

  // Clear accumulation arrays
  for (size_t i=0; i<nalong; i++) {
    for (size_t j=0; j<nacross; j++) {
      binval_ic[i][j] = 0.0;
      binval_wc[i][j] = 0.0;
      nobs_ic[i][j] = 0;
      nobs_wc[i][j] = 0;
    }
  }
    
  // Accumulate bin values and counts
  short *brptr=brow;
  short *bcptr=bcol;
  for (size_t i=0; i<totlines; i++) {
    for (size_t j=0; j<npixels; j++) {
      if (*brptr != -1) {

        if (cldprod[i][j] == -32767) {
          brptr++;
          bcptr++;
          continue;
        }
        
        if (cld_phase[i][j] == 1) {
          binval_ic[*brptr][*bcptr] += cldprod[i][j];
          nobs_ic[*brptr][*bcptr]++;
        } else {
          binval_wc[*brptr][*bcptr] += cldprod[i][j];
          nobs_wc[*brptr][*bcptr]++;
        }
      }
      brptr++;
      bcptr++;
    }
  }

  return 0;
}

// accumulate the ice and water cloud fraction
int accum_frac( uint32_t totlines, uint32_t npixels,
           uint32_t nalong, uint32_t nacross,
           short *brow, short *bcol,
           float **cldprod, int8_t **cld_phase,
           float **binval_ic, short **nobs_ic,
           float **binval_wc, short **nobs_wc) {

  // Clear accumulation arrays
  for (size_t i=0; i<nalong; i++) {
    for (size_t j=0; j<nacross; j++) {
      binval_ic[i][j] = 0.0;
      binval_wc[i][j] = 0.0;
      nobs_ic[i][j] = 0;
      nobs_wc[i][j] = 0;
    }
  }
    
  // Accumulate bin values and counts
  short *brptr=brow;
  short *bcptr=bcol;
  for (size_t i=0; i<totlines; i++) {
    for (size_t j=0; j<npixels; j++) {
      if (*brptr != -1) {

        if(cldprod[i][j] != -32767) {

          nobs_ic[*brptr][*bcptr]++;
          nobs_wc[*brptr][*bcptr]++;

          if (cldprod[i][j] == 1) {
            if (cld_phase[i][j] == 1) {
              binval_ic[*brptr][*bcptr]++;
            } else {
              binval_wc[*brptr][*bcptr]++;
            }
          }
        }

      }
      brptr++;
      bcptr++;
    }
  }

  return 0;
}


int accum_wm( uint32_t nwm, short *brow, short *bcol, float *wlval,
           float **binval, short **nobs) {

  // Accumulate bin values and counts
  short *brptr=brow;
  short *bcptr=bcol;
  for (size_t i=0; i<nwm; i++) {
    if (*brptr != -1) {

      if (wlval[i] == -32767) {
        brptr++;
        bcptr++;
        continue;
      }
        
      binval[*brptr][*bcptr] += wlval[i];
      nobs[*brptr][*bcptr]++;
    }
    brptr++;
    bcptr++;
  }

  return 0;
}


int lonlat2rowcol( uint32_t nalong, uint32_t nacross, uint32_t ncm,
                   float gridres, float *lonL1C, float *latL1C,
                   float *lonCM, float *latCM, short *brow, short *bcol);

int copyVarAtts( NcVar *varin, NcVar *varout, string override[],
                 string overridetype[]);

int copyVarAtts( NcVar *varin, NcVar *varout);

template<typename T>
int readCLD( NcGroup ncGrp[], const char *fldname, T *array,
             uint32_t npixels, uint32_t nlines[]) {

  vector<size_t> start, count;
  NcVar var;

  uint32_t nlines0 = 0;
  
  // Read from trailing granule if specified
  if ( !ncGrp[0].isNull()) {
    if(nlines[0] < 250) {
      nlines0 = nlines[0];
    } else {
      nlines0 = 250;
    }

    start.clear();
    start.push_back(nlines[0]-nlines0);
    start.push_back(0);

    count.clear();
    count.push_back(nlines0);
    count.push_back(npixels);

    var = ncGrp[0].getVar( fldname);
    var.getVar( start, count, &array[0]);
    
  }

  var = ncGrp[1].getVar( fldname);
  var.getVar( &array[nlines0*npixels]);

  // Read from following granule if specified
  if ( !ncGrp[2].isNull()) {
    start.clear();
    start.push_back(0);
    start.push_back(0);

    count.clear();
    if(nlines[2] < 250)
      count.push_back(nlines[2]);
    else
      count.push_back(250);
    count.push_back(npixels);

    var = ncGrp[2].getVar( fldname);
    var.getVar( start, count, &array[(nlines0+nlines[1])*npixels]);
  }

  return 0;
}

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


int main (int argc, char* argv[]) {

  NcVar var, varin, varout;

  clo_optionList_t *list;
  clo_option_t *option;
  char *strVal;
  char keyword [50];
  char targetfilename[FILENAME_MAX];
  char chlfilename[FILENAME_MAX];
  char profilename[FILENAME_MAX];
  char metfilename[FILENAME_MAX];
  char aerfilename[FILENAME_MAX];
  char geoscffilename[FILENAME_MAX];
  char cldmask1[FILENAME_MAX];
  char cldmask2[FILENAME_MAX];
  char cldmask3[FILENAME_MAX];
  char cldprod1[FILENAME_MAX];
  char cldprod2[FILENAME_MAX];
  char cldprod3[FILENAME_MAX];
  char albedofilename[FILENAME_MAX];
  char camsch4filename[FILENAME_MAX];
  char camsco2filename[FILENAME_MAX];
  char camsn2ofilename[FILENAME_MAX];
  char gebcofilename[FILENAME_MAX];
  char ancoutfilename[FILENAME_MAX];
  
  list = clo_createList();
  
  option = clo_addOption(list, "targetfile", CLO_TYPE_IFILE, NULL,
                         "Input target file");
  option = clo_addOption(list, "ancfile", CLO_TYPE_IFILE, NULL,
                         "Output ancillary file");
  option = clo_addOption(list, "merraprofile", CLO_TYPE_IFILE, NULL,
                         "Input MERRA2 profile ancillary file");
  option = clo_addOption(list, "merrametfile", CLO_TYPE_IFILE, NULL,
                         "Input MERRA2 MET ancillary file");
  option = clo_addOption(list, "merraaerfile", CLO_TYPE_IFILE, NULL,
                         "Input MERRA2 AER ancillary file (optional)");
  option = clo_addOption(list, "geoscffile", CLO_TYPE_IFILE, NULL,
                         "Input GEOS CF ancillary file");

  option = clo_addOption(list, "cldmask1", CLO_TYPE_IFILE, NULL,
                         "Input trailing cloud mask L2 file");
  option = clo_addOption(list, "cldmask2", CLO_TYPE_IFILE, NULL,
                         "Input current cloud mask L2 file");
  option = clo_addOption(list, "cldmask3", CLO_TYPE_IFILE, NULL,
                         "Input following cloud mask L2 file");

  option = clo_addOption(list, "cldprod1", CLO_TYPE_IFILE, NULL,
                         "Input trailing cloud product L2 file");
  option = clo_addOption(list, "cldprod2", CLO_TYPE_IFILE, NULL,
                         "Input current cloud product L2 file");
  option = clo_addOption(list, "cldprod3", CLO_TYPE_IFILE, NULL,
                         "Input following cloud product L2 file");
    

  option = clo_addOption(list, "albedofile", CLO_TYPE_IFILE, NULL,
                         "Input albedo ancillary file");
  option = clo_addOption(list, "chlfile", CLO_TYPE_IFILE, NULL,
                         "Input chlor_a L3 map file (optional)");
  option = clo_addOption(list, "ch4file", CLO_TYPE_IFILE, NULL,
                         "Input CAMS CH4 ancillary file");
  option = clo_addOption(list, "co2file", CLO_TYPE_IFILE, NULL,
                         "Input CAMS CO2 ancillary file");
  option = clo_addOption(list, "n2ofile", CLO_TYPE_IFILE, NULL,
                         "Input CAMS N2O ancillary file");

  option = clo_addOption(list, "gebcofile", CLO_TYPE_IFILE, NULL,
                         "Input GEBCO landmask file");
  
  clo_setVersion2("ancgen", "1.52");
  
  if (argc == 1) {
        clo_printUsage(list);
        exit(1);
  }

  clo_readArgs(list, argc, argv);
  
  strVal = clo_getString(list, "targetfile");
  strcpy(targetfilename, strVal);

  strVal = clo_getString(list, "ancfile");
  strcpy(ancoutfilename, strVal);
  
  strVal = clo_getString(list, "merraprofile");
  strcpy(profilename, strVal);
    
  strVal = clo_getString(list, "merrametfile");
  strcpy(metfilename, strVal);

  strVal = clo_getString(list, "geoscffile");
  strcpy(geoscffilename, strVal);

  //strVal = clo_getString(list, "cldmask2");
  //strcpy(cldmask2, strVal);

  //strVal = clo_getString(list, "cldprod2");
  //strcpy(cldprod2, strVal);

  strVal = clo_getString(list, "albedofile");
  strcpy(albedofilename, strVal);

  strVal = clo_getString(list, "ch4file");
  strcpy(camsch4filename, strVal);

  strVal = clo_getString(list, "co2file");
  strcpy(camsco2filename, strVal);

  strVal = clo_getString(list, "n2ofile");
  strcpy(camsn2ofilename, strVal);
  
  int numOptions, optionId;

  numOptions = clo_getNumOptions(list);
  for (optionId = 0; optionId < numOptions; optionId++) {
    option = clo_getOption(list, optionId);

    strcpy(keyword, option->key);
    
    if (strcmp(keyword, "chlfile") == 0) {
      if (clo_isOptionSet(option))
        parse_file_name(clo_getOptionString(option), chlfilename);
      else
        strcpy(chlfilename, "");

    } else if (strcmp(keyword, "merraaerfile") == 0) {
      if (clo_isOptionSet(option))
        parse_file_name(clo_getOptionString(option), aerfilename);
      else
        strcpy(aerfilename, "");

    } else if (strcmp(keyword, "cldmask1") == 0) {
      if (clo_isOptionSet(option))
        parse_file_name(clo_getOptionString(option), cldmask1);
      else
        strcpy(cldmask1, "");

    } else if (strcmp(keyword, "cldmask2") == 0) {
      if (clo_isOptionSet(option))
        parse_file_name(clo_getOptionString(option), cldmask2);
      else
        strcpy(cldmask2, "");

    } else if (strcmp(keyword, "cldmask3") == 0) {
      if (clo_isOptionSet(option))
        parse_file_name(clo_getOptionString(option), cldmask3);
      else
        strcpy(cldmask3, "");

    } else if (strcmp(keyword, "cldprod1") == 0) {
      if (clo_isOptionSet(option))
        parse_file_name(clo_getOptionString(option), cldprod1);
      else
        strcpy(cldprod1, "");

    } else if (strcmp(keyword, "cldprod2") == 0) {
      if (clo_isOptionSet(option))
        parse_file_name(clo_getOptionString(option), cldprod2);
      else
        strcpy(cldprod2, "");

    } else if (strcmp(keyword, "cldprod3") == 0) {
      if (clo_isOptionSet(option))
        parse_file_name(clo_getOptionString(option), cldprod3);
      else
        strcpy(cldprod3, "");
    } else if (strcmp(keyword, "gebcofile") == 0) {
      if (clo_isOptionSet(option))
        parse_file_name(clo_getOptionString(option), gebcofilename);
      else
        strcpy(gebcofilename, "$OCDATAROOT/common/gebco_ocssw_v2020.nc");
    }

  }

  ////////////////////////////////
  ///////// Open L1C file ////////
  ////////////////////////////////
  
  // Read lon/lat of l1c
  cout << endl << "Opening L1C file" << endl;
  NcFile *L1Cfile = new NcFile( targetfilename, NcFile::read);
  // get the start time
  NcGroupAtt att = L1Cfile->getAtt("time_coverage_start");
  if(att.isNull()) {
    cout << "Error - Could not find time_coverage_start global attribute.\n";
    exit(1);
  }
  string time_coverage_start;
  att.getValues(time_coverage_start);

  // get the end time
  att = L1Cfile->getAtt("time_coverage_end");
  if(att.isNull()) {
    cout << "Error - Could not find time_coverage_end global attribute.\n";
    exit(1);
  }
  string time_coverage_end;
  att.getValues(time_coverage_end);

  // get instrument
  att = L1Cfile->getAtt("instrument");
  if(att.isNull()) {
    cout << "Error - Could not find instrument global attribute.\n";
    exit(1);
  }
  string instrument;
  att.getValues(instrument);

  NcDim across_dim = L1Cfile->getDim("bins_across_track");
  uint32_t nacross = across_dim.getSize();
  NcDim along_dim = L1Cfile->getDim("bins_along_track");
  uint32_t nalong = along_dim.getSize();
  //  cout << "lines=" << nalong << ", pixels=" << nacross << endl;

  float **out_l1c_2d = allocate2d_float(nalong, nacross);
  float **latL1C = allocate2d_float(nalong, nacross);
  float **lonL1C = allocate2d_float(nalong, nacross);

  int bins_along_track = nalong;
  int bins_across_track = nacross;

  size_t nl1c = bins_along_track * bins_across_track;

  vector<NcDim> dims;

  float gridres=5.2;
  
  ////////////////////////////////
  /////// Create ANC file ////////
  ////////////////////////////////
  
  // Create output ancillary file
  cout << "Creating ancillary file" << endl;
  NcFile* nc_output;
  nc_output = new NcFile( ancoutfilename, NcFile::replace);
  // "PACE_SPEXONE_ANC.20170115T225222.L1C.5km.nc",

  // Global metadata
  string date_created = string(unix2isodate(now(), 'G'));

  nc_output->putAtt( "date_created", date_created);
  nc_output->putAtt( "time_coverage_start", time_coverage_start);
  nc_output->putAtt( "time_coverage_end", time_coverage_end);
  nc_output->putAtt( "title", instrument + " L1C ancillary file");
  nc_output->putAtt( "product_name", basename(ancoutfilename));
  nc_output->putAtt( "creator_name", "NASA/GSFC/OBPG");
  nc_output->putAtt( "creator_url", "http://oceandata.sci.gsfc.nasa.gov");
  nc_output->putAtt( "creator_email", "data@oceancolor.gsfc.nasa.gov");
  nc_output->putAtt( "project",
                     "Ocean Biology Processing Group (NASA/GSFC/OBPG)");
  nc_output->putAtt( "publisher_name", "NASA/GSFC/OBPG");
  nc_output->putAtt( "publisher_url", "http://oceandata.sci.gsfc.nasa.gov");
  nc_output->putAtt( "publisher_email", "data@oceancolor.gsfc.nasa.gov");
  

  // Create netCDF dimensions
  NcDim rows = nc_output->addDim("bins_along_track", bins_along_track);
  NcDim cols = nc_output->addDim("bins_across_track", bins_across_track);

  // Add target lat/lon fields
  dims.push_back(rows);
  dims.push_back(cols);

  // Get L1C lon/lat
  NcGroup gidGEO = L1Cfile->getGroup( "geolocation_data");

  float *fptr;
  
  var = gidGEO.getVar( "latitude");
  var.getVar( &latL1C[0][0]);
  varout = nc_output->addVar("latitude", ncFloat, dims);

  string override_ll[] = {"-32767", "=", "=", "=", "="};
  string overridetype_ll[] = {"F", "=", "=", "=", "="};
  copyVarAtts( &var, &varout, override_ll, overridetype_ll);

  // Convert fill values to -32767
  fptr = &latL1C[0][0];
  for (size_t i=0; i<nl1c; i++) {
    if (*fptr < -900) *fptr = -32767;
    fptr++;
  }
  varout.putVar( &latL1C[0][0]);
  
  var = gidGEO.getVar( "longitude");
  var.getVar( &lonL1C[0][0]);
  varout = nc_output->addVar("longitude", ncFloat, dims);

  copyVarAtts( &var, &varout, override_ll, overridetype_ll);

  // Convert fill values to -32767
  fptr = &lonL1C[0][0];
  for (size_t i=0; i<nl1c; i++) {
    if (*fptr < -900) *fptr = -32767;
    fptr++;
  }
  varout.putVar( &lonL1C[0][0]);

  
  // Get lon/lat min/max for L1C granule
  float lonmax=-180, lonmin=+180;
  fptr = &lonL1C[0][0];
  for (size_t i=0; i< (size_t) nl1c; i++) {
    if (*fptr < lonmin) lonmin = *fptr;
    if (*fptr > lonmax) lonmax = *fptr;
    fptr++;
  }
  
  float latmax=-90, latmin=+90;
  fptr = &latL1C[0][0];
  for (size_t i=0; i< (size_t) nl1c; i++) {
    if (*fptr < latmin) latmin = *fptr;
    if (*fptr > latmax) latmax = *fptr;
    fptr++;
  }
  
  vector<size_t> start, count;
  vector<ptrdiff_t> stride;
 
  ////////////////////////////////
  //////////// ALBEDO ////////////
  ////////////////////////////////
  size_t albedo_npix = 43200;
  size_t albedo_nlines = 21600;
  size_t albedo_nfields = 5;
  float albedo_pixel_size = 0.00833333;

  cout << endl << "Opening ALBEDO file" << endl;
  NcFile *albedofile = new NcFile( albedofilename, NcFile::read);

  uint8_t **albedorow = new uint8_t*[albedo_nlines];

  std::vector<string> ALBEDOfields
    = {"Albedo_Map_0.659", "Albedo_Map_0.858",
       "Albedo_Map_1.24", "Albedo_Map_1.64",
       "Albedo_Map_2.13"};

  int32_t syear;
  int32_t sday;
  int32_t smsec;
  isodate2ydmsec((char*)time_coverage_start.c_str(), &syear, &sday, &smsec);
  int32_t period_8day = sday / 8;
  start.clear();
  start.push_back(period_8day);
  start.push_back(0);
  start.push_back(0);

  count.clear();
  count.push_back(1);
  count.push_back(NCACHE);
  count.push_back(albedo_npix);

  // Loop through albedo fields
  for (size_t k=0; k<albedo_nfields; k++) {
    var = albedofile->getVar( ALBEDOfields[k].c_str());

    for (size_t i=0; i<albedo_nlines; i++) 
      albedorow[i] = NULL;

    int first=-1;
    for (size_t i=0; i<nalong; i++) {
      for (size_t j=0; j<nacross; j++) {
        //        cout << "i,j: " << i << " " << j << endl;
        int latrow = (latL1C[i][j] + 90) / albedo_pixel_size;

        if (albedorow[latrow] == NULL) {

          // align the row to the cache size boundry
          size_t cache_boundary = ((size_t)(latrow / NCACHE)) * NCACHE;
          if (first == -1) first = cache_boundary;
          albedorow[cache_boundary] = new uint8_t[albedo_npix*NCACHE];
          for (size_t k=1; k<NCACHE; k++)
            albedorow[cache_boundary+k]
              = albedorow[cache_boundary] + k*albedo_npix;
          
          start[1] = cache_boundary;
          // cout << "Reading: " << latrow << " " << start[1] << " " << count[1]
          //   << endl;
          var.getVar( start, count, albedorow[cache_boundary]);
        } 
        int loncol = (lonL1C[i][j] + 180) / albedo_pixel_size;
        //  cout << "latrow, loncol: " << latrow << " " << loncol << endl;
        uint8_t albedo = albedorow[latrow][loncol];
        if (albedo == 255)
          out_l1c_2d[i][j] = -32767;
        else
          out_l1c_2d[i][j] = albedo * 0.004;
      }
    }

    cout << "Writing " << ALBEDOfields[k].c_str() << endl;
    varout = nc_output->addVar(ALBEDOfields[k].c_str(), ncFloat, dims);

    string override[] = {"-32767", "", "=", "-32767", "", "=", "="};
    string overridetype[] = {"F", "", "=", "F", "", "=", "="};
    copyVarAtts( &var, &varout, override, overridetype);

    varout.putVar(&out_l1c_2d[0][0]);

    for (size_t i=first; i<albedo_nlines; i+=NCACHE)
      if (albedorow[i] != NULL) delete [] albedorow[i];

  } // albedo field loop

  delete[] albedorow;

  // SW to NE (landmask)

  
  // Determine if L1C granule crosses the dateline
  // Adjust L1C lon to 0-360 if dateline crossed
  bool datelinecross=false;
  if (lonmax - lonmin > 180) {
    cout << endl << "Correcting longitude for dataline crossing" << endl;
    datelinecross = true;
    lonmax=0; lonmin=360;
    fptr = &lonL1C[0][0];
    for (size_t i=0; i< (size_t) nl1c; i++) {
      if (*fptr < 0) *fptr += 360;
      if (*fptr < lonmin) lonmin = *fptr;
      if (*fptr > lonmax) lonmax = *fptr;
      fptr++;
    }
  }

  
  ////////////////////////////////
  //////////// MERRA2 ////////////
  ////////////////////////////////
  
  // Read MERRA2 ancillary file
  int nlon=576;
  int nlat=361;
  
  cout << endl << "Opening MERRA2 file" << endl;
  NcFile *MERRAfile = new NcFile( profilename, NcFile::read);
  // "N202200903_PROFILE_GMAOFP_3h.nc", NcFile::read);
  
  float** lat_merra = allocate2d_float(361, 576); 
  float** lon_merra = allocate2d_float(361, 576);

  for (size_t i=0; i<361; i++) {
    for (size_t j=0; j<576; j++) {
      lat_merra[i][j] = i*0.50 - 90;
      lon_merra[i][j] = j*0.625 - 180;

      if (datelinecross && lon_merra[i][j] < 0) lon_merra[i][j] += 360;
    }
  }
  int nmerra = nlon*nlat;

  // Get nearest neighbor lon/lat bin numbers for MERRA
  short **latbinL1C = allocate2d_short(nalong, nacross);
  short **lonbinL1C = allocate2d_short(nalong, nacross);

  fptr = &latL1C[0][0];
  short *slatptr = &latbinL1C[0][0];
  for (size_t i=0; i< (size_t) nl1c; i++) {
    *slatptr = (short) ((*fptr + 90) / 0.50);
    fptr++;
    slatptr++;
  }

  fptr = &lonL1C[0][0];
  short *slonptr = &lonbinL1C[0][0];
  for (size_t i=0; i< (size_t) nl1c; i++) {
    *slonptr = (short) ((*fptr + 180) / 0.625);
    fptr++;
    slonptr++;
  }

  cout << "Computing MERRA to L1C interpolation mapping" << endl;
  int nnc_tag_merra = cgal_nnc( nmerra, &lat_merra[0][0], &lon_merra[0][0],
                                nl1c, &latL1C[0][0], &lonL1C[0][0]);

  size_t nlev = 42;

  NcDim levs = nc_output->addDim("levels", nlev);

  dims.clear();
  dims.push_back(levs);
  dims.push_back(rows);
  dims.push_back(cols);

  float ***out_l1c_3d = allocate3d_float(nlev, nalong, nacross);
  float ***merra_data = allocate3d_float(42, 361, 576);

  // Geo-potential height profile
  // Temperature vertical profile
  // Relative humidity profile
  // Specific humidity profile
  // Ozone profile

  string MERRAfields[5]={"H", "T", "RH", "QV", "O3"};

  // Product loop
  for (size_t k=0; k<5; k++) {
    
    varin = MERRAfile->getVar( MERRAfields[k].c_str());
    varin.getVar(&merra_data[0][0][0]);

    // Convert MERRA2 fill values to -1e14
    float *fptr = &merra_data[0][0][0];
    for (size_t i=0; i<42*361*576; i++) {
      if (*fptr > 1e14) *fptr = -1e14;
      fptr++;
    }

    // Level loop
    for (size_t i=0; i<nlev; i++) {
      cgal_interp2( nmerra,
                    &lat_merra[0][0], &lon_merra[0][0], &merra_data[i][0][0],
                    nl1c, &out_l1c_3d[i][0][0], nnc_tag_merra);

      // Convert interpolation failures to -32767
      fptr = &out_l1c_3d[0][0][0];
      for (size_t i=0; i<nlev*nalong*nacross; i++) {
        if (*fptr < 0) *fptr = -32767;
        fptr++;
      }

      // Attempt to correct fill values with nearest-neighbor values
      fptr = &out_l1c_3d[i][0][0];
      slatptr = &latbinL1C[0][0];
      slonptr = &lonbinL1C[0][0];

      for (size_t j=0; j< (size_t) nl1c; j++) {
        if (*fptr == -32767)
          *fptr = merra_data[i][*slatptr][*slonptr];

        fptr++;
        slatptr++;
        slonptr++;
      }
    }

    cout << "Writing " << MERRAfields[k] << endl;

    varout = nc_output->addVar(MERRAfields[k].c_str(), ncFloat, dims);

    string override[] = {"-32767", "=", "-32767", "=", "="};
    string overridetype[] = {"F", "=", "F", "=", "="};
    copyVarAtts( &varin, &varout, override, overridetype);
    
    varout.putVar(&out_l1c_3d[0][0][0]);
    
  } // fields loop

  delete MERRAfile;

  
  ////////////////////////////////
  ////////// MERRA2 MET //////////
  ////////////////////////////////

  cout << endl << "Opening MERRA2 MET file" << endl;
  MERRAfile = new NcFile( metfilename, NcFile::read);
  // "GMAO_MERRA2.20220321T000000.MET.nc"

  dims.clear();
  dims.push_back(rows);
  dims.push_back(cols);

  string MERRA_MET_fields[10]={"PS", "QV10M", "SLP", "T10M", "TO3",
                               "TQV", "U10M", "V10M", "FRSNO", "FRSEAICE"};

  for (size_t k=0; k<10; k++) {
    
    varin = MERRAfile->getVar( MERRA_MET_fields[k].c_str());
    varin.getVar(&merra_data[0][0][0]);

    // Convert MERRA2 fill values to -1e14
    float *fptr = &merra_data[0][0][0];
    for (size_t i=0; i<361*576; i++) {
      if (*fptr > 1e14) *fptr = -1e14;
      fptr++;
    }
    
    cgal_interp2( nmerra,
                  &lat_merra[0][0], &lon_merra[0][0], &merra_data[0][0][0],
                  nl1c, &out_l1c_3d[0][0][0], nnc_tag_merra);

    // Convert interpolation failures to -32767
    fptr = &out_l1c_3d[0][0][0];
    for (size_t i=0; i<nalong*nacross; i++) {
    if (*fptr < 0) *fptr = -32767;
    fptr++;
    }
    

    // Attempt to correct fill values with nearest-neighbor values
    fptr = &out_l1c_3d[0][0][0];
    slatptr = &latbinL1C[0][0];
    slonptr = &lonbinL1C[0][0];

    for (size_t j=0; j< (size_t) nl1c; j++) {
      if (*fptr == -32767)
      *fptr = merra_data[0][*slatptr][*slonptr];

      fptr++;
      slatptr++;
      slonptr++;
    }
    
    cout << "Writing " << MERRA_MET_fields[k] << endl;

    varout = nc_output->addVar(MERRA_MET_fields[k].c_str(), ncFloat, dims);

    string override[] = {"-32767", "=", "-32767", "=", "="};
    string overridetype[] = {"F", "=", "F", "=", "="};
    copyVarAtts( &varin, &varout, override, overridetype);
    varout.putVar(&out_l1c_3d[0][0][0]);

  } // fields loop


  ////////////////////////////////
  ////////// MERRA2 AER //////////
  ////////////////////////////////

  if (strcmp(aerfilename, "") != 0) {
    cout << endl << "Opening MERRA2 AER_file" << endl;
    MERRAfile = new NcFile( aerfilename, NcFile::read);
    // "GMAO_MERRA2.20220321T000000.AER.nc"

    dims.clear();
    dims.push_back(rows);
    dims.push_back(cols);

    string MERRA_AER_fields[13]={"BCEXTTAU", "BCSCATAU", "DUEXTTAU",
                                 "DUSCATAU", "SSEXTTAU", "SSSCATAU",
                                 "SUEXTTAU", "SUSCATAU", "OCEXTTAU",
                                 "OCSCATAU", "TOTEXTTAU","TOTSCATAU",
                                 "TOTANGSTR"};

    for (size_t k=0; k<13; k++) {
    
      varin = MERRAfile->getVar( MERRA_AER_fields[k].c_str());
      varin.getVar(&merra_data[0][0][0]);

      // Convert MERRA2 fill values to -1e14
      float *fptr = &merra_data[0][0][0];
      for (size_t i=0; i<361*576; i++) {
        if (*fptr > 1e14) *fptr = -1e14;
        fptr++;
      }
    
      cgal_interp2( nmerra,
                    &lat_merra[0][0], &lon_merra[0][0], &merra_data[0][0][0],
                    nl1c, &out_l1c_3d[0][0][0], nnc_tag_merra);

      // Convert interpolation failures to -32767
      fptr = &out_l1c_3d[0][0][0];
      for (size_t i=0; i<nalong*nacross; i++) {
        if (*fptr < 0) *fptr = -32767;
        fptr++;
      }
      

      // Attempt to correct fill values with nearest-neighbor values
      fptr = &out_l1c_3d[0][0][0];
      slatptr = &latbinL1C[0][0];
      slonptr = &lonbinL1C[0][0];

      for (size_t j=0; j< (size_t) nl1c; j++) {
        if (*fptr == -32767)
        *fptr = merra_data[0][*slatptr][*slonptr];

        fptr++;
        slatptr++;
        slonptr++;
      }
    
      cout << "Writing " << MERRA_AER_fields[k] << endl;

      varout = nc_output->addVar(MERRA_AER_fields[k].c_str(), ncFloat, dims);

      string override[] = {"-32767", "=", "-32767", "=", "="};
      string overridetype[] = {"F", "=", "F", "=", "="};
      copyVarAtts( &varin, &varout, override, overridetype);
      varout.putVar(&out_l1c_3d[0][0][0]);

    } // fields loop

  
    //  cgal_release_tag( nnc_tag_merra);
    //  free2d_float(lon_merra);
    //free2d_float(lat_merra);
    free2d_short(lonbinL1C);
    free2d_short(latbinL1C);
    free3d_float(merra_data);
  }

  /////////////////////////////////
  //////////// CAMS SAT ///////////
  /////////////////////////////////

  NcFile *CAMSfile = new NcFile( camsch4filename, NcFile::read);

  NcDim camslat_dim = CAMSfile->getDim("latitude_bins");
  uint32_t camslat = camslat_dim.getSize();

  NcDim camslon_dim = CAMSfile->getDim("longitude_bins");
  uint32_t camslon = camslon_dim.getSize();

  float **lat_cams = allocate2d_float(camslat, camslon);
  float **lon_cams = allocate2d_float(camslat, camslon); 

  for (size_t i=0; i<camslat; i++) {
    for (size_t j=0; j<camslon; j++) {
      lat_cams[i][j] = i*2 - 89.0;
      lon_cams[i][j] = j*3 - 178.5;

      if (datelinecross && lon_cams[i][j] < 0) lon_cams[i][j] += 360;
    }
  }
  int ncams = camslon*camslat;
    
  cout << endl <<"Computing CAMS SATELLITE to L1C interpolation mapping"
       << endl;
  int nnc_tag_cams = cgal_nnc( ncams, &lat_cams[0][0], &lon_cams[0][0],
                               nl1c, &latL1C[0][0], &lonL1C[0][0]);

  dims.clear();
  dims.push_back(levs);
  dims.push_back(rows);
  dims.push_back(cols);

  float ***cams_data = allocate3d_float(nlev, camslat, camslon);

  int month = atoi(time_coverage_start.substr(5,2).c_str());

  string s = to_string(month);
  unsigned int number_of_zeros = 2 - s.length();
  s.insert(0, number_of_zeros, '0');
  s.insert(0, "CH4_");
      
  varin = CAMSfile->getVar( s.c_str());
  varin.getVar(&cams_data[0][0][0]);

  // Convert CAMS fill values to -32767
  fptr = &cams_data[0][0][0];
  for (size_t i=0; i<nlev*camslat*camslon; i++) {
    if (*fptr == -999) *fptr = -32767;
    fptr++;
  }

  for (size_t l=0; l<42; l++) {
    cgal_interp2( ncams,
                  &lat_cams[0][0], &lon_cams[0][0], &cams_data[l][0][0],
                  nl1c, &out_l1c_3d[l][0][0], nnc_tag_cams);
  }

  // Convert interpolation failures to -32767
  fptr = &out_l1c_3d[0][0][0];
  for (size_t i=0; i<nlev*nalong*nacross; i++) {
    if (*fptr < 0) *fptr = -32767;
    fptr++;
  }
    
  s.assign("CH4");
  cout << "Writing " << s << endl;
    
  varout = nc_output->addVar(s.c_str(), ncFloat, dims);

  string override[] = {"-32767", "=", "=", "=", "=", "=", "=", "="};
  string overridetype[] = {"F", "=", "=", "=", "=", "=", "=", "="};
  copyVarAtts( &varin, &varout, override, overridetype);

  varout.putVar(&out_l1c_3d[0][0][0]);

  delete CAMSfile;
  
  free2d_float(lat_cams);
  free2d_float(lon_cams);
  free3d_float(cams_data);

  cgal_release_tag( nnc_tag_cams);

    
  /////////////////////////////////
  //////////// CAMS INST //////////
  /////////////////////////////////

  CAMSfile = new NcFile( camsco2filename, NcFile::read);

  camslat_dim = CAMSfile->getDim("latitude_bins");
  camslat = camslat_dim.getSize();

  camslon_dim = CAMSfile->getDim("longitude_bins");
  camslon = camslon_dim.getSize();

  lat_cams = allocate2d_float(camslat, camslon);
  lon_cams = allocate2d_float(camslat, camslon);

  for (size_t i=0; i<camslat; i++) {
    for (size_t j=0; j<camslon; j++) {
      lat_cams[i][j] = i*1.8947368421 - 90.0;
      lon_cams[i][j] = j*3.75 - 180.0;

      if (datelinecross && lon_cams[i][j] < 0) lon_cams[i][j] += 360;
    }
  }
  ncams = camslon*camslat;

  cout << endl << "Computing CAMS INST to L1C interpolation mapping"
       << endl;
  nnc_tag_cams = cgal_nnc( ncams, &lat_cams[0][0], &lon_cams[0][0],
                           nl1c, &latL1C[0][0], &lonL1C[0][0]);

  cams_data = allocate3d_float(nlev, camslat, camslon);

  s = to_string(month);
  number_of_zeros = 2 - s.length();
  s.insert(0, number_of_zeros, '0');
  s.insert(0, "CO2_");

  varin = CAMSfile->getVar( s.c_str());
  varin.getVar(&cams_data[0][0][0]);

  // Convert CAMS fill values to -32767
  fptr = &cams_data[0][0][0];
  for (size_t i=0; i<nlev*camslat*camslon; i++) {
    if (*fptr == -999) *fptr = -32767;
    fptr++;
  }
    
  for (size_t l=0; l<42; l++) {
    cgal_interp2( ncams,
                  &lat_cams[0][0], &lon_cams[0][0], &cams_data[l][0][0],
                  nl1c, &out_l1c_3d[l][0][0], nnc_tag_cams);
  }

  // Convert interpolation failures to -32767
  fptr = &out_l1c_3d[0][0][0];
  for (size_t i=0; i<nlev*nalong*nacross; i++) {
    if (*fptr < 0) *fptr = -32767;
    fptr++;
  }
    
  s.assign("CO2");
  cout << "Writing " << s << endl;

  varout = nc_output->addVar(s.c_str(), ncFloat, dims);
  copyVarAtts( &varin, &varout, override, overridetype);
    
  varout.putVar(&out_l1c_3d[0][0][0]);

  delete CAMSfile;

  CAMSfile = new NcFile( camsn2ofilename, NcFile::read);

  s = to_string(month);
  number_of_zeros = 2 - s.length();
  s.insert(0, number_of_zeros, '0');
  s.insert(0, "N2O_");

  varin = CAMSfile->getVar( s.c_str());
  varin.getVar(&cams_data[0][0][0]);

  // Convert CAMS fill values to -32767
  fptr = &cams_data[0][0][0];
  for (size_t i=0; i<nlev*camslat*camslon; i++) {
    if (*fptr == -999) *fptr = -32767;
    fptr++;
  }
    
  for (size_t l=0; l<42; l++) {
    cgal_interp2( ncams,
                  &lat_cams[0][0], &lon_cams[0][0], &cams_data[l][0][0],
                  nl1c, &out_l1c_3d[l][0][0], nnc_tag_cams);
  }

  // Convert interpolation failures to -32767
  fptr = &out_l1c_3d[0][0][0];
  for (size_t i=0; i<nlev*nalong*nacross; i++) {
    if (*fptr < 0) *fptr = -32767;
    fptr++;
  }
    
  s.assign("N2O");
  cout << "Writing " << s << endl;

  varout = nc_output->addVar(s.c_str(), ncFloat, dims);
  copyVarAtts( &varin, &varout, override, overridetype);

  varout.putVar(&out_l1c_3d[0][0][0]);

  delete CAMSfile;
  
  free2d_float(lat_cams);
  free2d_float(lon_cams);
  free3d_float(cams_data);

  cgal_release_tag( nnc_tag_cams);


  /////////////////////////////////
  //////////// GEOS-CF ////////////
  /////////////////////////////////
  
  // Read and process GEOS-CF data

  cout << endl << "Opening GEOS-CF file" << endl;
  NcFile *GEOSfile = new NcFile( geoscffilename, NcFile::read);
  //"GMAO_GEOS-CF.20210428T020000.PROFILE.nc"
  
  nlon=1440;
  nlat=721;
  int ngeos = nlon*nlat;

  float **lat_geos = allocate2d_float(721, 1440);
  float **lon_geos = allocate2d_float(721, 1440);

  // Compute lon/lat geos
  for (size_t i=0; i<721; i++) {
    for (size_t j=0; j<1440; j++) {
      lat_geos[i][j] = i*0.25 - 90;
      lon_geos[i][j] = j*0.25 - 180;

      if (datelinecross && lon_geos[i][j] < 0) lon_geos[i][j] += 360;
    }
  }

  // Write level field
  dims.clear();
  dims.push_back(levs);
  varout = nc_output->addVar("levels", ncFloat, dims);

  float levels[nlev];
  varin = GEOSfile->getVar( "lev");
  varin.getVar(levels);
  copyVarAtts( &varin, &varout);
  varout.putVar(levels);

  float ***geos_data = allocate3d_float(42, 721, 1440);

  cout << "Computing GEOS to L1C interpolation mapping" << endl;
  int nnc_tag_geos = cgal_nnc( ngeos, &lat_geos[0][0], &lon_geos[0][0],
                                  nl1c, &latL1C[0][0], &lonL1C[0][0]);

  string GEOSfields[5]={"NO2", "SO2",
                        "TOTCOL_NO2", "STRATCOL_NO2", "TROPCOL_NO2"};

  for (size_t k=0; k<5; k++) {
    
    varin = GEOSfile->getVar( GEOSfields[k].c_str());
    varin.getVar(&geos_data[0][0][0]);

    if (GEOSfields[k].compare("NO2") == 0 ||
        GEOSfields[k].compare("SO2") == 0)
      nlev = 42; else nlev = 1;

    // Convert GOES fill values to -1e14
    float *fptr = &geos_data[0][0][0];
    for (size_t i=0; i<nlev*721*1440; i++) {
      if (*fptr > 1e14) *fptr = -1e14;
      fptr++;
    }

    for (size_t i=0; i<nlev; i++) {
      cgal_interp2( ngeos,
                       &lat_geos[0][0], &lon_geos[0][0], &geos_data[i][0][0],
                       nl1c, &out_l1c_3d[i][0][0], nnc_tag_geos);
    }

    // Convert interpolation failures to -32767
    fptr = &out_l1c_3d[0][0][0];
    for (size_t i=0; i<nlev*nalong*nacross; i++) {
      if (*fptr < 0) *fptr = -32767;
      fptr++;
    }
    
    dims.clear();
    if (GEOSfields[k].compare("NO2") == 0  ||
        GEOSfields[k].compare("SO2") == 0) {
      dims.push_back(levs);
      dims.push_back(rows);
      dims.push_back(cols);
    } else {
      dims.push_back(rows);
      dims.push_back(cols);
    }

    cout << "Writing " << GEOSfields[k].c_str() << endl;

    varout = nc_output->addVar(GEOSfields[k].c_str(), ncFloat, dims);

    string overrideGEOS[] = {"-32767", "=", "-32767", "=", "="};
    string overridetypeGEOS[] = {"F", "=", "F", "=", "="};
    copyVarAtts( &varin, &varout, overrideGEOS, overridetypeGEOS);
    
    varout.putVar(&out_l1c_3d[0][0][0]);
  } // fields loop
  
  cgal_release_tag( nnc_tag_geos);
  free2d_float(lat_geos);
  free2d_float(lon_geos);
  free3d_float(geos_data);
  free3d_float(out_l1c_3d);

  
  ///////////////////////////////
  //////////// CLOUD ////////////
  ///////////////////////////////

  if ( strcmp(cldmask2, "") != 0 &&
       strcmp(cldprod2, "") != 0) {
  
    cout << endl << "Opening CLDMSK files" << endl;

    NcFile *CLDMSKfile[3]={NULL, NULL, NULL};
  
    uint32_t nlines[3]={0,0,0};
    uint32_t npixels;
    NcDim lines_dim;
    NcDim pixel_dim;

    uint32_t nlines0 = 0;
  
    // Open trailing/following L2 CLDMASK datasets
    if ( strcmp(cldmask1, "") != 0) {
      CLDMSKfile[0] = new NcFile( cldmask1, NcFile::read);
      lines_dim = CLDMSKfile[0]->getDim("number_of_lines");
      nlines[0] = lines_dim.getSize();
    }

    if ( strcmp(cldmask3, "") != 0) {
      CLDMSKfile[2] = new NcFile( cldmask3, NcFile::read);
      lines_dim = CLDMSKfile[2]->getDim("number_of_lines");
      nlines[2] = lines_dim.getSize();
    }

    // Open current L2 LDMASK dataset
    CLDMSKfile[1] = new NcFile( cldmask2, NcFile::read);
    lines_dim = CLDMSKfile[1]->getDim("number_of_lines");
    nlines[1] = lines_dim.getSize();

    pixel_dim = CLDMSKfile[1]->getDim("pixels_per_line");
    npixels = pixel_dim.getSize();

    // Use only 250 (or less) lines from trailing/following datasets
    uint32_t totlines = nlines[1];
    if (nlines[0] != 0) {
      if (nlines[0] < 250)
        totlines += nlines[0];
      else
        totlines += 250;
    }
        
    if (nlines[2] != 0) {
      if (nlines[2] < 250)
        totlines += nlines[2];
      else
        totlines += 250;
    }
  
    float **latCM = allocate2d_float(totlines, npixels);
    float **lonCM = allocate2d_float(totlines, npixels);
    int8_t **adj_mask = allocate2d_schar(totlines, npixels);

    for (size_t i=0; i<totlines; i++) {
      for (size_t j=0; j<npixels; j++) {
        latCM[i][j] = 0.0;
        lonCM[i][j] = 0.0;
        adj_mask[i][j] = 0;
      }
    }
  
    NcGroup ncGroup;
    // Read from trailing granule if specified
    if ( CLDMSKfile[0] != NULL) {
      if(nlines[0] < 250) {
        nlines0 = nlines[0];
        start[0] = 0;
      } else {
        nlines0 = 250;
        start[0] = nlines[0] - 250;
      }
      count[0] = nlines0;
      start[1] = 0;
      count[1] = npixels;

      ncGroup = CLDMSKfile[0]->getGroup("navigation_data");
      if(ncGroup.isNull()) {
        var = CLDMSKfile[0]->getVar( "latitude");
        var.getVar( start, count, &latCM[0][0]);

        var = CLDMSKfile[0]->getVar( "longitude");
        var.getVar( start, count, &lonCM[0][0]);

        var = CLDMSKfile[0]->getVar( "ADJ_MASK");
        var.getVar( start, count, &adj_mask[0][0]);
      } else {
        var = ncGroup.getVar( "latitude");
        var.getVar( start, count, &latCM[0][0]);

        var = ncGroup.getVar( "longitude");
        var.getVar( start, count, &lonCM[0][0]);

        ncGroup = CLDMSKfile[0]->getGroup("geophysical_data");
        var = ncGroup.getVar("cloud_flag");
        var.getVar( start, count, &adj_mask[0][0]);
      }
    }

    ncGroup = CLDMSKfile[1]->getGroup("navigation_data");
    if(ncGroup.isNull()) {
      var = CLDMSKfile[1]->getVar( "latitude");
      var.getVar( &latCM[nlines0][0]);

      var = CLDMSKfile[1]->getVar( "longitude");
      var.getVar( &lonCM[nlines0][0]);

      var = CLDMSKfile[1]->getVar( "ADJ_MASK");
      var.getVar( &adj_mask[nlines0][0]);
    } else {
      var = ncGroup.getVar( "latitude");
      var.getVar( &latCM[nlines0][0]);

      var = ncGroup.getVar( "longitude");
      var.getVar( &lonCM[nlines0][0]);

      ncGroup = CLDMSKfile[1]->getGroup("geophysical_data");
      var = ncGroup.getVar( "cloud_flag");
      var.getVar( &adj_mask[nlines0][0]);
    }    

    // Read from following granule if specified
    if ( CLDMSKfile[2] != NULL) {
      start[0] = 0;
      if(nlines[2] < 250) {
        count[0] = nlines[2];
      } else {
        count[0] = 250;
      }
      start[1] = 0;
      count[1] = npixels;

      ncGroup = CLDMSKfile[2]->getGroup("navigation_data");
      if(ncGroup.isNull()) {
        var = CLDMSKfile[2]->getVar( "latitude");
        var.getVar( start, count, &latCM[nlines0+nlines[1]][0]);

        var = CLDMSKfile[2]->getVar( "longitude");
        var.getVar( start, count, &lonCM[nlines0+nlines[1]][0]);

        var = CLDMSKfile[2]->getVar( "ADJ_MASK");
        var.getVar( start, count, &adj_mask[nlines0+nlines[1]][0]);
      } else {
        var = ncGroup.getVar( "latitude");
        var.getVar( start, count, &latCM[nlines0+nlines[1]][0]);

        var = ncGroup.getVar( "longitude");
        var.getVar( start, count, &lonCM[nlines0+nlines[1]][0]);

        ncGroup = CLDMSKfile[2]->getGroup("geophysical_data");
        var = ncGroup.getVar( "cloud_flag");
        var.getVar( start, count, &adj_mask[nlines0+nlines[1]][0]);
      }
    }

    uint32_t ncm = totlines*npixels;


    // Determine L1C row/col of L2 CLD lon/lat
    // Ported from code developed by F.Patt

    short *brow = new short[ncm];
    short *bcol = new short[ncm];
    lonlat2rowcol( nalong, nacross, ncm, gridres, &lonL1C[0][0], &latL1C[0][0],
                   &lonCM[0][0], &latCM[0][0], brow, bcol);

    ofstream tempOut;

    //  tempOut.open ("rc.bin", ios::out | ios::trunc | ios::binary);
    //tempOut.write((char *) brow, sizeof(short)*ncm);
    //tempOut.write((char *) bcol, sizeof(short)*ncm);
    //tempOut.close();


    // Open CLD products file

    // Note: CLD files uses same lon/lat as CLDMSK

    NcFile *CLDfile[3];
    NcGroup ncGrp[3];

    if ( strcmp(cldprod1, "") != 0) {
      CLDfile[0] = new NcFile( cldprod1, NcFile::read);
      ncGrp[0] = CLDfile[0]->getGroup( "geophysical_data");
    }

    if ( strcmp(cldprod3, "") != 0) {
      CLDfile[2] = new NcFile( cldprod3, NcFile::read);
      ncGrp[2] = CLDfile[2]->getGroup( "geophysical_data");
    }

    CLDfile[1] = new NcFile( cldprod2, NcFile::read);
    ncGrp[1] = CLDfile[1]->getGroup( "geophysical_data");

    int8_t **cld_phase = allocate2d_schar(totlines, npixels);

    // Read cloud phase
    readCLD<int8_t>( ncGrp, "cld_phase_21", &cld_phase[0][0], npixels, nlines);
      
    // Set ice cloud pixels to 1 otherwise 0
    for (size_t i=0; i<totlines; i++)
      for (size_t j=0; j<npixels; j++)
        if (cld_phase[i][j] == 3)
          cld_phase[i][j] = 1;
        else
          cld_phase[i][j] = 0;

    // Allocate bin variables
    float **binval_ic = allocate2d_float(nalong, nacross);
    float **binval_wc = allocate2d_float(nalong, nacross);
    short **nobs_ic = allocate2d_short(nalong, nacross);
    short **nobs_wc = allocate2d_short(nalong, nacross);

    float **cldprod = allocate2d_float(totlines, npixels);

    // Average ADJ_MASK for water/ice cloud pixels

    // Convert BYTE to FLOAT
    for (size_t i=0; i<totlines; i++)
      for (size_t j=0; j<npixels; j++)
        if(adj_mask[i][j] == -128)
          cldprod[i][j] = -32767;
        else
          cldprod[i][j] = (float) adj_mask[i][j];
  
    free2d_schar(adj_mask);

    //  tempOut.open ("adj_mask.bin", ios::out | ios::trunc | ios::binary);
    //tempOut.write((char *) &cldprod[0][0], sizeof(float)*totlines*npixels);
    //tempOut.close();
  
    accum_frac( totlines, npixels, nalong, nacross, brow, bcol, cldprod, cld_phase,
           binval_ic, nobs_ic, binval_wc, nobs_wc);

    // Compute average bin values for ice clouds
    for (size_t i=0; i<nalong; i++) {
      for (size_t j=0; j<nacross; j++) {
        if (nobs_ic[i][j] != 0)
          out_l1c_2d[i][j] = binval_ic[i][j] / nobs_ic[i][j];
        else
          out_l1c_2d[i][j] = -32767;
      }
    }

    dims.clear();
    dims.push_back(rows);
    dims.push_back(cols);

    cout << "Writing ice cloud fraction" << endl;
    varout = nc_output->addVar("ice_cloud_fraction", ncFloat, dims);

    string override_icf[] =
      {"-32767", "", "", "", "", "Ice cloud fraction", "", "1.0", "0.0"};
    string overridetype_icf[] = {"F", "", "", "", "", "=", "", "F", "F"};
    copyVarAtts( &var, &varout, override_icf, overridetype_icf);

    varout.putVar(&out_l1c_2d[0][0]);

    // Compute average bin values for water clouds
    for (size_t i=0; i<nalong; i++) {
      for (size_t j=0; j<nacross; j++) {
        if (nobs_wc[i][j] != 0)
          out_l1c_2d[i][j] = binval_wc[i][j] / nobs_wc[i][j];
        else 
          out_l1c_2d[i][j] = -32767;
      }
    }

    cout << "Writing water cloud fraction" << endl;
    varout = nc_output->addVar("water_cloud_fraction", ncFloat, dims);

    string override_wcf[] =
      {"-32767", "", "", "", "", "Water cloud fraction", "", "1.0", "0.0"}; 
    string overridetype_wcf[] = {"F", "", "", "", "", "=", "", "F", "F"};
    copyVarAtts( &var, &varout, override_wcf, overridetype_wcf);

    varout.putVar(&out_l1c_2d[0][0]);
  

    // Cloud Top Pressure (CTP)
    // Cloud Top Temperature (CTT)
    // Cloud Top Height (CTH)
    // Particle effective radius cer_21 (micron)
    // COT cot_21

    string CLDfields[5]={"ctp", "ctt", "cth", "cer_21", "cot_21"};

    for (size_t k=0; k<5; k++) {
      cout << "Writing " << CLDfields[k].c_str() << endl;
    
      var = ncGrp[1].getVar( CLDfields[k].c_str());

      for (size_t i=0; i<totlines; i++)
        for (size_t j=0; j<npixels; j++) cldprod[i][j] = 0.0;
    
      readCLD<float>( ncGrp, CLDfields[k].c_str(), &cldprod[0][0],
                      npixels, nlines);

      accum( totlines, npixels, nalong, nacross, brow, bcol, cldprod, cld_phase,
             binval_ic, nobs_ic, binval_wc, nobs_wc);

      // Compute average bin values for ice clouds
      for (size_t i=0; i<nalong; i++) {
        for (size_t j=0; j<nacross; j++) {
          if (nobs_ic[i][j] != 0)
            out_l1c_2d[i][j] = binval_ic[i][j] / nobs_ic[i][j];
          else
            out_l1c_2d[i][j] = -32767;
        }
      }

      string outname;

      outname.assign(CLDfields[k]);
      outname.append("_ice_cloud");
    
      varout = nc_output->addVar(outname.c_str(), ncFloat, dims);
      copyVarAtts( &var, &varout);
      varout.putVar(&out_l1c_2d[0][0]);

      // Compute average bin values for water clouds
      for (size_t i=0; i<nalong; i++) {
        for (size_t j=0; j<nacross; j++) {
          if (nobs_wc[i][j] != 0)
            out_l1c_2d[i][j] = binval_wc[i][j] / nobs_wc[i][j];
          else
            out_l1c_2d[i][j] = -32767;
        }
      }
    
      outname.assign(CLDfields[k]);
      outname.append("_water_cloud");
    
      varout = nc_output->addVar(outname.c_str(), ncFloat, dims);
      copyVarAtts( &var, &varout);    
      varout.putVar(&out_l1c_2d[0][0]);
    } // CLD field loop

    delete [] brow;
    delete [] bcol;

    free2d_float(latCM);
    free2d_float(lonCM);
    free2d_float(cldprod);

    free2d_float(binval_ic);
    free2d_float(binval_wc);
    free2d_short(nobs_ic);
    free2d_short(nobs_wc);
  } // if CLD processing

  
  ///////////////////////////////
  ////////// WATERMASK //////////
  ///////////////////////////////

  float wm_pixel_size = 360.0 / 86400;
  int wm_latrow = (latmin +  90) / wm_pixel_size;
  int wm_loncol = (lonmin + 180) / wm_pixel_size;

  //cout << latmin << " " << latmax << " " << lonmin << " " << lonmax << endl;
  //cout << wm_latrow << " " << wm_loncol << endl;
  
  start.clear();
  start.push_back(wm_latrow);
  start.push_back(wm_loncol);

  size_t n_wm_latrow = ((latmax +  90) / wm_pixel_size - wm_latrow + 1)/4;
  size_t n_wm_loncol = ((lonmax + 180) / wm_pixel_size - wm_loncol + 1)/4;
  
  count.clear();
  count.push_back(n_wm_latrow);
  count.push_back(n_wm_loncol);

  stride.clear();
  stride.push_back(4);
  stride.push_back(4);

  float *lon_wm = new float[n_wm_latrow*n_wm_loncol];
  float *lat_wm = new float[n_wm_latrow*n_wm_loncol];

  //cout << n_wm_latrow << " " << n_wm_loncol << endl;
  
  int k=0;
  for (size_t i=0; i<(size_t) n_wm_latrow; i++) {
    for (size_t j=0; j<(size_t) n_wm_loncol; j++) {
      lat_wm[k] = (4*i+wm_latrow)*wm_pixel_size - 90;
      lon_wm[k] = (4*j+wm_loncol)*wm_pixel_size - 180;
      
      if (datelinecross && lon_wm[k] < 0) lon_wm[k] += 360;
      k++;
    }
  }
  
  uint32_t nwm = n_wm_latrow * n_wm_loncol;

  // Reading landmask
  string gebcofilestring = gebcofilename;
  expandEnvVar(&gebcofilestring);
  
  NcFile *WMfile = new NcFile( gebcofilestring.c_str(), NcFile::read);

  float *wm = new float[n_wm_latrow*n_wm_loncol];
  uint8_t *wmbyte = new uint8_t[n_wm_latrow*n_wm_loncol];
  var = WMfile->getVar( "watermask");
  if (datelinecross) {
      uint8_t *wm0 = new uint8_t[n_wm_latrow*n_wm_loncol];
      count[1] = (86400 - wm_loncol) / 4;
      var.getVar( start, count, stride, &wm0[0]);
      for (size_t i=0; i<(size_t) n_wm_latrow; i++)
        memcpy(&wmbyte[i*n_wm_loncol], &wm0[i*count[1]],
               count[1]*sizeof(uint8_t));

      start[1] = 0;
      count[1] = n_wm_loncol - count[1];
      //cout << count[1] << endl;
      var.getVar( start, count, stride, &wm0[0]);
      for (size_t i=0; i<(size_t) n_wm_latrow; i++)
        memcpy(&wmbyte[i*n_wm_loncol+(86400-wm_loncol)/4],
               &wm0[i*count[1]], count[1]*sizeof(uint8_t));
      delete [] wm0;
  } else {
    var.getVar( start, count, stride, &wmbyte[0]);
  }

  // Convert byte to float
  k=0;
  for (size_t i=0; i<(size_t) n_wm_latrow; i++) {
    for (size_t j=0; j<(size_t) n_wm_loncol; j++) {
      wm[k] = (float) wmbyte[k];
      k++;
    }
  }


  // Determine L1C row/col of L2 CLD lon/lat

  short *brow = new short[nwm];
  short *bcol = new short[nwm];

  float **binval_wm = allocate2d_float(nalong, nacross);
  short **nobs_wm = allocate2d_short(nalong, nacross);

  // Clear accumulation arrays
  for (size_t i=0; i<nalong; i++) {
    for (size_t j=0; j<nacross; j++) {
      binval_wm[i][j] = 0.0;
      nobs_wm[i][j] = 0;
    }
  }

  size_t N = 2147483648 / (n_wm_loncol*nalong);
  size_t M = n_wm_latrow/N;
  if (n_wm_latrow % N != 0) M++;
  for (size_t i=0; i<M; i++) {
    // cout << i*N << " " << n_wm_latrow << endl;
    size_t L=N;
    if (n_wm_latrow - i*N < N) L = n_wm_latrow - i*N;
    //    cout << "L: " << L << endl;
    lonlat2rowcol( nalong, nacross, L*n_wm_loncol, gridres,
                   &lonL1C[0][0], &latL1C[0][0],
                   &lon_wm[i*N*n_wm_loncol], &lat_wm[i*N*n_wm_loncol],
                   brow, bcol);

    accum_wm( L*n_wm_loncol, brow, bcol, &wm[i*N*n_wm_loncol],
              binval_wm, nobs_wm);
  }

  for (size_t i=0; i<nalong; i++) {
    for (size_t j=0; j<nacross; j++) {

      if (nobs_wm[i][j] != 0)
        out_l1c_2d[i][j] = binval_wm[i][j] / nobs_wm[i][j];
      else
        out_l1c_2d[i][j] = -32767;
    }
  }
    
  cout << "Writing waterfraction" << endl;
  varout = nc_output->addVar("waterfraction", ncFloat, dims);

  string override_wf[] =
    {"-32767", "=", "Water fraction", "waterfraction", "1.0", "0.0"};
  string overridetype_wf[] = {"F", "", "", "=", "F", "F"};
  copyVarAtts( &var, &varout, override_wf, overridetype_wf);
  varout.putVar(&out_l1c_2d[0][0]);
  
  delete [] lat_wm;
  delete [] lon_wm;
  delete [] wm;

  delete [] brow;
  delete [] bcol;
  
  free2d_float(binval_wm);
  free2d_short(nobs_wm);

  
  ///////////////////////////////
  //////////// CHLOR ////////////
  ///////////////////////////////
 
  k=0;
  
  if (strcmp(chlfilename, "") != 0) {
    float chl_pixel_size = 360.0 / 8640;
    int chl_latrow = (90 - latmax) / chl_pixel_size;
    int chl_loncol = (lonmin + 180) / chl_pixel_size;

    cout << latmin << " " << latmax << " " << lonmin << " " << lonmax << endl;
    cout << chl_latrow << " " << chl_loncol << endl;
  
    start.clear();
    start.push_back(chl_latrow);
    start.push_back(chl_loncol);

    int n_chl_latrow = (90 - latmin) / chl_pixel_size - chl_latrow + 1;
    int n_chl_loncol = (lonmax + 180) / chl_pixel_size - chl_loncol + 1;

    count.clear();
    count.push_back(n_chl_latrow);
    count.push_back(n_chl_loncol);

    float *lon_chl = new float[n_chl_latrow*n_chl_loncol];
    float *lat_chl = new float[n_chl_latrow*n_chl_loncol];

    cout << n_chl_latrow << " " << n_chl_loncol << endl;
  
    k=0;
    for (size_t i=0; i<(size_t) n_chl_latrow; i++) {
      for (size_t j=0; j<(size_t) n_chl_loncol; j++) {
        lat_chl[k] = 90 - (i+chl_latrow)*chl_pixel_size;
        lon_chl[k] = (j+chl_loncol)*chl_pixel_size - 180;
      
        if (datelinecross && lon_chl[k] < 0) lon_chl[k] += 360;
        k++;
      }
    }
  
    int nchl = n_chl_latrow * n_chl_loncol;

    cout << "Computing CHLOR_A map to L1C interpolation mapping" << endl;
    int nnc_tag_chl = cgal_nnc( nchl, &lat_chl[0], &lon_chl[0],
                                nl1c, &latL1C[0][0], &lonL1C[0][0]);

    // Reading chlor_a from L3M file
    NcFile *CHLfile = new NcFile( chlfilename, NcFile::read);

    float *chl = new float[n_chl_latrow*n_chl_loncol];
    var = CHLfile->getVar( "chlor_a");

    if (datelinecross) {
      float *chl0 = new float[n_chl_latrow*n_chl_loncol];
      count[1] = 8640 - chl_loncol;
      var.getVar( start, count, &chl0[0]);
      for (size_t i=0; i<(size_t) n_chl_latrow; i++)
        memcpy(&chl[i*n_chl_loncol], &chl0[i*count[1]], count[1]*sizeof(float));

      start[1] = 0;
      count[1] = n_chl_loncol + chl_loncol - 8640;
      cout << count[1] << endl;
      var.getVar( start, count, &chl0[0]);
      for (size_t i=0; i<(size_t) n_chl_latrow; i++)
        memcpy(&chl[i*n_chl_loncol+(8640-chl_loncol)],
               &chl0[i*count[1]], count[1]*sizeof(float));
      delete [] chl0;
    } else {
      var.getVar( start, count, &chl[0]);
    }
  
    cgal_interp2( nchl, &lat_chl[0], &lon_chl[0], &chl[0],
                  nl1c, &out_l1c_2d[0][0], nnc_tag_chl);

    // Convert interpolation failures to -32767
    fptr = &out_l1c_2d[0][0];
    for (size_t i=0; i<nalong*nacross; i++) {
      if (*fptr < 0) *fptr = -32767;
      fptr++;
    }

    cout << "Writing chlor_a" << endl;
    varout = nc_output->addVar("chlor_a", ncFloat, dims);
    copyVarAtts( &var, &varout);
    varout.putVar(&out_l1c_2d[0][0]);
  
    cgal_release_tag(nnc_tag_chl);
    delete [] lat_chl;
    delete [] lon_chl;
    delete [] chl;
  }

  free2d_float(out_l1c_2d);
  free2d_float(latL1C);
  free2d_float(lonL1C);

  return 0;
}

int lonlat2rowcol( uint32_t nalong, uint32_t nacross, uint32_t ncm,
                   float gridres, float *lonL1C, float *latL1C,
                   float *lonCM, float *latCM, short *brow, short *bcol) {

  ofstream tempOut;
    
  // Determine L1C row/col of L2 CLD lon/lat
  // Ported from code developed by F.Patt

  float *fptrlon, *fptrlat;

  // Convert lon/lat to unit vectors (L1C)
  float ***gvec = allocate3d_float(nalong, nacross, 3);
  fptrlon = lonL1C;
  fptrlat = latL1C;
  for (size_t i=0; i<nalong; i++) {
    for (size_t j=0; j<nacross; j++) {
      gvec[i][j][0] = cos(*fptrlon/RADEG) * cos(*fptrlat/RADEG);
      gvec[i][j][1] = sin(*fptrlon/RADEG) * cos(*fptrlat/RADEG);
      gvec[i][j][2] = sin(*fptrlat/RADEG);

      fptrlon++;
      fptrlat++;
    }
  }

  // Convert lon/lat to unit vectors (L2 CLD)
  float **bvec = allocate2d_float(ncm, 3);
  fptrlon = lonCM;
  fptrlat = latCM;
  for (size_t i=0; i<ncm; i++) {
    bvec[i][0] = cos(*fptrlon/RADEG)*cos(*fptrlat/RADEG);
    bvec[i][1] = sin(*fptrlon/RADEG)*cos(*fptrlat/RADEG);
    bvec[i][2] = sin(*fptrlat/RADEG);

    fptrlon++;
    fptrlat++;
  }

  //tempOut.open ("gb.bin", ios::out | ios::trunc | ios::binary);
  //tempOut.write((char *) &gvec[0][0][0], sizeof(float)*nalong*nacross*3);
  //tempOut.write((char *) &bvec[0][0], sizeof(float)*ncm*3);
  //tempOut.close();

  // Compute normal vectors for L1C rows
  // (cross product of first and last vector in each row)
  float **gnvec = allocate2d_float(nalong, 3);
  float *gnvm = new float[nalong];
  float vecm[3];
  for (size_t i=0; i<nalong; i++) {
    gnvec[i][0] =
      gvec[i][nacross-1][1]*gvec[i][0][2] -
      gvec[i][nacross-1][2]*gvec[i][0][1];

    gnvec[i][1] =
      gvec[i][nacross-1][2]*gvec[i][0][0] -
      gvec[i][nacross-1][0]*gvec[i][0][2];

    gnvec[i][2] =
      gvec[i][nacross-1][0]*gvec[i][0][1] -
      gvec[i][nacross-1][1]*gvec[i][0][0];
    
    vecm[0] = gnvec[i][0];
    vecm[1] = gnvec[i][1];
    vecm[2] = gnvec[i][2];
    
    gnvm[i] = sqrt(vecm[0]*vecm[0] + vecm[1]*vecm[1] + vecm[2]*vecm[2]);
  }
  for (size_t i=0; i<nalong; i++) {
    gnvec[i][0] /= gnvm[i];
    gnvec[i][1] /= gnvm[i];
    gnvec[i][2] /= gnvm[i];
  }
  delete [] gnvm;

  //tempOut.open ("gnvec.bin", ios::out | ios::trunc | ios::binary);
  //tempOut.write((char *) &gnvec[0][0], sizeof(float)*nalong*3);
  //tempOut.close();
  
  // Compute dot products of L1B vectors with normal vectors

  gsl_matrix_float_view G =
    gsl_matrix_float_view_array(&gnvec[0][0], nalong, 3);
  gsl_matrix_float_view B =
    gsl_matrix_float_view_array(&bvec[0][0], ncm, 3);

  float **bdotgn = allocate2d_float(ncm, nalong);
  gsl_matrix_float_view bdotgn_mat =
    gsl_matrix_float_view_array(&bdotgn[0][0], ncm, nalong);

  //cout << "matrix multiplication" << endl;
  gsl_blas_sgemm(CblasNoTrans, CblasTrans, 1.0, &B.matrix, &G.matrix,
                 0.0, &bdotgn_mat.matrix);


  //  tempOut.open ("bd.bin", ios::out | ios::trunc | ios::binary);
  //tempOut.write((char *) &bdotgn[0][0], sizeof(float)*ncm*nalong);
  //tempOut.close();
  
  // Determine row and column for each pixel within grid
  // Compute normals to center columns
  uint32_t ic = nacross / 2;
  float **cnvec = allocate2d_float(nalong, 3);
  for (size_t i=0; i<nalong; i++) {
    cnvec[i][0] = gvec[i][ic][1]*gnvec[i][2] - gvec[i][ic][2]*gnvec[i][1];
    cnvec[i][1] = gvec[i][ic][2]*gnvec[i][0] - gvec[i][ic][0]*gnvec[i][2];
    cnvec[i][2] = gvec[i][ic][0]*gnvec[i][1] - gvec[i][ic][1]*gnvec[i][0];
  }
  
  // Compute grid row nadir resolution
  float *dcm = new float[nalong];
  for (size_t i=0; i<nalong; i++) {
    vecm[0] = gvec[i][ic+1][0] - gvec[i][ic][0];
    vecm[1] = gvec[i][ic+1][1] - gvec[i][ic][1];
    vecm[2] = gvec[i][ic+1][2] - gvec[i][ic][2];

    dcm[i] = sqrt(vecm[0]*vecm[0] + vecm[1]*vecm[1] + vecm[2]*vecm[2]);
  }


  float db = gridres/6311/2;

  //  tempOut.open ("ibad.bin", ios::out | ios::trunc | ios::binary);

  for (size_t i=0; i<ncm; i++) {
    if (bdotgn[i][0] > db || bdotgn[i][nalong-1] < -db) {
      brow[i] = -1;
      bcol[i] = -1;
      //      tempOut.write((char *) &i, sizeof(size_t));

    } else {
      brow[i] = -1;
      bcol[i] = -1;
      //if (i % 100000 == 0) cout << i << endl;
      for (size_t j=0; j<nalong; j++) {
        if (bdotgn[i][j] > -db && bdotgn[i][j] <= db) {
          brow[i] = j;

          float bdotcn = 0.0;
          for (size_t k=0; k<3; k++)
            bdotcn += (bvec[i][k] * cnvec[j][k]);

          bcol[i] = (short) round(ic + bdotcn/dcm[j]);

          if (bcol[i] < 0 || bcol[i] >= (short) nacross) {
            bcol[i] = -1;
            brow[i] = -1;
          }
          break;
        }
      }
    }
  }

  //  tempOut.close();

  free3d_float(gvec);
  free2d_float(bvec);
  free2d_float(gnvec);
  free2d_float(cnvec);
  free2d_float(bdotgn);
  
  delete [] dcm;
  
  return 0;
}

