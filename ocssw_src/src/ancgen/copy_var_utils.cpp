#include <stdio.h>

#include <iostream>
#include "copy_var_utils.hpp"

using namespace std;
using namespace netCDF;
using namespace netCDF::exceptions;


// CAMS file have variable attributes as strings so we can just copy it over to the new ANC file
// string CH4:cell_measures = "area: cell_area" ;   (NO) This is what's in CAMS files.
// CH4:cell_measures = "area: cell_area" ;          (YES) 
int copyVarAttsCams(netCDF::NcVar *camsVar, netCDF::NcVar *outVar) {

    // grab all the attributes names
    map<string, NcVarAtt> camsVarAttributes = camsVar->getAtts();

    for (const auto &attribute : camsVarAttributes) {
        string attributeName = attribute.first;
        NcVarAtt attributeData = attribute.second;
        
        // skip global attributes. Dont need to be in variables
        if (attributeName == "time_coverage_start" || attributeName == "time_coverage_end") {
            continue;
        }
        // catch fill value because the fill value was changed from CAMS to ANC
        else if (attributeName == "_FillValue") {
            outVar->setFill(true, (void*)&FLOAT_FILL_VALUE);
        }
        // everything else
        else {

            // pointer to the first char index of the attribute string
            char* attributeStr = nullptr;
            attributeData.getValues(&attributeStr);

            /*
                Need to store it as a NC_CHAR because if you use NC_STRING it will look like:
                    string CH4:cell_methods = "area: mean time: mean"
                We want it to look like:
                    CH4:cell_methods = "area: mean time: mean"
                So it needs to be stored as an NC_CHAR and specifying the length
            */
            outVar->putAtt(attributeName, NC_CHAR, strlen(attributeStr), attributeStr);
            nc_free_string(1, &attributeStr);
        }
    }
    return 0;
}


/**
 * @brief copy attributes from input to output netcdf variable
 * @param varin input pointer to variable
 * @param varout output pointer to variable
 * @param override fill values with number of fields length
 * @param overrideType fill value types with number of fields length
 * @return
 */
int copyVarAtts(NcVar *varin, NcVar *varout, string override[], string overrideType[]) {
    char buffer[1000];

    map<string, NcVarAtt> attributes;
    NcVarAtt vattr;
    attributes = varin->getAtts();

    int k = 0;
    for (map<string, NcVarAtt>::iterator it = attributes.begin(); it != attributes.end(); ++it) {
        vattr = (*it).second;
        NcType type = vattr.getType();
        size_t attlen = vattr.getAttLength();

        string attrName = (*it).first;
        ;

        if (override == NULL || override[k].compare("=") == 0) {
            vattr.getValues((void *)buffer);
        } else if (override[k].compare("") == 0) {
            k++;
            continue;
        } else {
            if (overrideType[k] == "B") {
                int8_t temp = (uint8_t)stoi(override[k].c_str());
                memcpy(buffer, &temp, 1);
            } else if (overrideType[k] == "F") {
                float temp = stof(override[k].c_str());
                memcpy(buffer, &temp, 4);
            } else {
                memcpy(buffer, override[k].c_str(), override[k].length());
            }
        }

        if (attrName.compare("_FillValue") == 0) {
            varout->setFill(true, (void *)buffer);
        } else {
            if (override == NULL || override[k].compare("=") == 0)
                varout->putAtt(attrName, type, attlen, buffer);
            else if (overrideType[k] == "F")
                varout->putAtt(attrName, NC_FLOAT, attlen, buffer);
            else
                varout->putAtt(attrName, type, override[k].length(), buffer);
        }

        k++;
    }

    return 0;
}


// copy attributes from input to output netcdf variable but without ovrride
int copyVarAtts(NcVar *varin, NcVar *varout) {
    copyVarAtts(varin, varout, NULL, NULL);

    return 0;
}
