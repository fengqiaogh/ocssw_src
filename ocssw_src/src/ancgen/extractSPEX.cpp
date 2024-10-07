#include <stdio.h>
#include <math.h>
#include <libgen.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <netcdf>

#include <allocate2d.h>
#include <allocate3d.h>
#include <timeutils.h>
#include <clo.h>

using namespace std;
using namespace netCDF;
using namespace netCDF::exceptions;

int copyVarAtts( NcVar *varin, NcVar *varout, string override[],
                 string overridetype[]);

int copyVarAtts( NcVar *varin, NcVar *varout);

int main (int argc, char* argv[]) {

  NcVar var, varin, varout;

  clo_optionList_t *list;
  clo_option_t *option;
  char *strVal;
  char keyword [50];
  char inancfilename[FILENAME_MAX];
  char outancfilename[FILENAME_MAX];
  
  list = clo_createList();
  
  option = clo_addOption(list, "inancfile", CLO_TYPE_IFILE, NULL,
                         "Input ancillary file file");

  option = clo_addOption(list, "outancfile", CLO_TYPE_OFILE, NULL,
                         "Output ancillary file file");

  clo_setVersion("0.1");
  
  if (argc == 1) {
        clo_printUsage(list);
        exit(1);
  }

  clo_readArgs(list, argc, argv);
  
  strVal = clo_getString(list, "inancfile");
  strcpy(inancfilename, strVal);

  strVal = clo_getString(list, "outancfile");
  strcpy(outancfilename, strVal);

  int numOptions, optionId;

  numOptions = clo_getNumOptions(list);
  for (optionId = 0; optionId < numOptions; optionId++) {
    option = clo_getOption(list, optionId);

    strcpy(keyword, option->key);

  }

  ////////////////////////////////
  ////// Open input ANC file /////
  ////////////////////////////////
  
  cout << endl << "Opening ancillary file" << endl;
  NcFile *inANCfile = new NcFile( inancfilename, NcFile::read);

  NcDim across_dim = inANCfile->getDim("bins_across_track");
  uint32_t nacross = across_dim.getSize();
  NcDim along_dim = inANCfile->getDim("bins_along_track");
  uint32_t nalong = along_dim.getSize();

  NcDim levs_dim = inANCfile->getDim("levels");
  uint32_t nlev = levs_dim.getSize();
  //  cout << "lines=" << nalong << ", pixels=" << nacross << endl;

  uint32_t nacross_out=25;
  
  float **out_2d = allocate2d_float(nalong, nacross_out);
  float ***out_3d = allocate3d_float(nlev, nalong, nacross_out);

  int bins_along_track = nalong;

  ////////////////////////////////
  //// Create output ANC file ////
  ////////////////////////////////
  
  // Create output ancillary file
  cout << "Creating ancillary file" << endl;
  NcFile* nc_output;
  nc_output = new NcFile( outancfilename, NcFile::replace);

  // Global metadata
  string date_created = string(unix2isodate(now(), 'G'));

  char buffer[1000];
  multimap<string, NcGroupAtt> ancillaryGAtts = inANCfile->getAtts();
  multimap<string, NcGroupAtt>::iterator itratt;
  for (itratt = ancillaryGAtts.begin(); itratt != ancillaryGAtts.end();
       ++itratt) {
    string attrname = itratt->first;
    NcGroupAtt attin = itratt->second;

    NcType type = attin.getType();
    size_t attlen = attin.getAttLength();
    
    cout << itratt->first << " " << attlen << endl;

    if (attrname.compare("date_created") == 0) {
      nc_output->putAtt( "date_created", date_created);
    } else if (attrname.compare("product_name") == 0) {
      nc_output->putAtt( "product_name", basename(outancfilename));
    } else if (attrname.compare("title") == 0) {
      nc_output->putAtt( "title", "SPEXone L1C ancillary file");
    } else {
      attin.getValues((void *) buffer);
      nc_output->putAtt( attrname, type, attlen, buffer);
    }
  }

  // Create netCDF dimensions
  NcDim rows = nc_output->addDim("bins_along_track", bins_along_track);
  NcDim cols = nc_output->addDim("bins_across_track", nacross_out);
  NcDim levs = nc_output->addDim("levels", nlev);

  vector<NcDim> dims;
  vector<size_t> start, count;

  multimap<string, NcVar> ancillaryVars = inANCfile->getVars();

  multimap<string, NcVar>::iterator itr;
  for (itr = ancillaryVars.begin(); itr != ancillaryVars.end(); ++itr) {
    string fieldname = itr->first;
    varin = itr->second;
    
    cout << itr->first << endl;

    int ndims = varin.getDimCount();

    if (ndims == 2) {
      start.clear();
      start.push_back(0);
      start.push_back((nacross - nacross_out) / 2);

      count.clear();
      count.push_back(nalong);
      count.push_back(nacross_out);
      
      varin.getVar( start, count, &out_2d[0][0]);

      dims.clear();
      dims.push_back(rows);
      dims.push_back(cols);

      varout = nc_output->addVar(fieldname, ncFloat, dims);
      copyVarAtts( &varin, &varout);
      varout.putVar( &out_2d[0][0]);
    } else if (ndims == 3) {
      start.clear();
      start.push_back(0);
      start.push_back(0);
      start.push_back((nacross - nacross_out) / 2);

      count.clear();
      count.push_back(nlev);
      count.push_back(nalong);
      count.push_back(nacross_out);
      
      varin.getVar( start, count, &out_3d[0][0][0]);

      dims.clear();
      dims.push_back(levs);
      dims.push_back(rows);
      dims.push_back(cols);

      varout = nc_output->addVar(fieldname, ncFloat, dims);
      copyVarAtts( &varin, &varout);
      varout.putVar( &out_3d[0][0][0]);
    }
  }
  
  free2d_float(out_2d);
  free3d_float(out_3d);
  
  return 0;
}
