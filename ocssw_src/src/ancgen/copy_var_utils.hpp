#ifndef _COPY_VAR_UTILS_HPP_
#define _COPY_VAR_UTILS_HPP_
#include <stdint.h>
#include <netcdf>

const float FLOAT_FILL_VALUE = -32767.0;


/**
 * @brief copy CAMS file variable (store as all strings) info into the output var 
 * @param camsVar pointer to variable from the CAMS file
 * @param outVar pointer to variable from the output ANC file
 * @return 
 */
int copyVarAttsCams(netCDF::NcVar *camsVar, netCDF::NcVar *outVar);

/**
 * @brief copy attributes from input to output netcdf variable
 * @param varin input pointer to variable
 * @param varout output pointer to variable
 * @param override fill values with number of fields length
 * @param overrideType fill value types with number of fields length
 * @return
 */
int copyVarAtts(netCDF::NcVar *varin, netCDF::NcVar *varout, std::string override[], std::string overrideType[]);

/**
 * @brief calls copyVarAtts but without override and overrideType
 * @param varin 
 * @param varout 
 * @return 
 */
int copyVarAtts(netCDF::NcVar *varin, netCDF::NcVar *varout);

#endif