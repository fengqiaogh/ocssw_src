#include "scaled_nc_var.hpp"
#include "l2_reader.hpp"
// Simple program to test l2_reader
// IT is not a part of CTEST. It needs to be developed later.
int main(int argc, char** argv) {
    std::string filename;
    std::string varname;
    if (argc > 2) {
        filename = argv[1];
        varname = argv[2];
    }
    else return 1;
    netCDF::NcFile file_nc(filename,netCDF::NcFile::read );
    auto var_nc = find_nc_variable_cpp(varname,file_nc);
    std::cout << get_full_nc_path(var_nc) << std::endl;
    return 0;
}