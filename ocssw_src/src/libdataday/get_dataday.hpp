#ifndef GET_DATADAY_HPP
#define GET_DATADAY_HPP
#include <utility>
#include <cstdlib>
#include <iostream>
#include <netcdf>
#include "get_dataday.h"
std::pair<int32_t, int32_t> get_datadays(const netCDF::NcFile& nc_input,float deltaeqcross = 0.0f, int night_flag = 0);
int32_t get_plusday();
#endif
