#include <stdio.h>

#include <iostream>
#include <netcdf>
#include <cstdint>

using namespace std;
using namespace netCDF;
using namespace netCDF::exceptions;

int copyVarAtts( NcVar *varin, NcVar *varout, string override[],
                 string overridetype[]) {

  char buffer[1000];
    
  // Copy field attributes
  map <string, NcVarAtt> attributes;
  NcVarAtt vattr;
  attributes = varin->getAtts();

  int k=0;
  for (map<string, NcVarAtt>:: iterator it=attributes.begin();
       it!=attributes.end(); ++it) {

    vattr = (*it).second;
    NcType type = vattr.getType();
    size_t attlen = vattr.getAttLength();

    //    cout << (*it).first << " " << attlen << " "
    //        << type.getSize() << endl;
    //if (override != NULL)       cout << override[k].c_str() << endl;

    string attrname = (*it).first;;

    if (override == NULL || override[k].compare("=") == 0) {
      vattr.getValues((void *) buffer);
    } else if (override[k].compare("") == 0) {
      k++;
      continue;
    } else {
      //cout << override[k].c_str() << " " << overridetype[k].c_str() << endl;
      if ( overridetype[k] == "B") {
        int8_t temp = (uint8_t) stoi(override[k].c_str());
        memcpy( buffer, &temp, 1);
      } else if ( overridetype[k] == "F") {
        float temp = stof(override[k].c_str());
        memcpy( buffer, &temp, 4);
      } else {
        memcpy( buffer, override[k].c_str(), override[k].length());
      }
    }
    
    if (attrname.compare("_FillValue") == 0) {
      varout->setFill(true, (void *) buffer);
    } else {
      if (override == NULL || override[k].compare("=") == 0)
        varout->putAtt( attrname, type, attlen, buffer);
      else if ( overridetype[k] == "F")
        varout->putAtt( attrname, NC_FLOAT, attlen, buffer);
      else
        varout->putAtt( attrname, type, override[k].length(), buffer);
    }

    k++;
  }
  
  return 0;
}

int copyVarAtts( NcVar *varin, NcVar *varout) {

  copyVarAtts( varin, varout, NULL, NULL);

  return 0;
}

