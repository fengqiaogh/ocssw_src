#ifndef FIND_VARIABLE_HPP
#define FIND_VARIABLE_HPP
#include <netcdf>
template <typename T>
netCDF::NcVar find_nc_variable_cpp(const std::string& var_name, T& nc_file) {
    netCDF::NcVar var = nc_file.getVar(var_name);
    if (var.isNull()) {
        auto grps = nc_file.getGroups();
        for (const auto& grp : grps) {
            var = find_nc_variable_cpp(var_name, grp.second);
            if (!var.isNull()) {
                return var;
            }
        }
    }
    return var;
};

template <typename T>
netCDF::NcGroupAtt find_nc_grp_attribute_cpp(const std::string& att_name, T& nc_file) {
    netCDF::NcGroupAtt attr = nc_file.getAtt(att_name);
    if (attr.isNull()) {
        auto grps = nc_file.getGroups();
        for (const auto& grp : grps) {
            attr = find_nc_grp_attribute_cpp(att_name, grp.second);
            if (!attr.isNull()) {
                return attr;
            }
        }
    }
    return attr;
};

template <typename T>
std::multimap<std::string, netCDF::NcVar> find_all_variables(T& nc_id) {
    auto grps = nc_id.getGroups();
    std::multimap<std::string, netCDF::NcVar> vars = nc_id.getVars();
    for (const auto& grp : grps) {
        std::multimap<std::string, netCDF::NcVar> temp_var = find_all_variables(grp.second);
        vars.insert(temp_var.begin(), temp_var.end());
    }
    return vars;
}
#endif