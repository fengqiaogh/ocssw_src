#ifndef _GLOBAL_ATTRS_
#define _GLOBAL_ATTRS_

#include <netcdf>
#include <timeutils.h>

std::string call_sequence(int argc, char* argv[]);

std::string get_history(netCDF::NcFile *ncfile);

void set_global_attrs(std::string filename,
                      std::string call_sequence,
                      std::string doi="",
                      std::string pversion="");

void set_global_attrs(netCDF::NcFile *ncfile,
                      std::string call_sequence,
                      std::string doi="",
                      std::string pversion="");

#endif /* _GLOBAL_ATTRS_ */
