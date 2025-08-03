#ifndef FIND_VARIABLE_HPP
#define FIND_VARIABLE_HPP
#include <netcdf>
#include <unordered_set>
#define NC_CHECK(...)                                                                                    \
    {                                                                                                    \
        try {                                                                                            \
            __VA_ARGS__;                                                                                 \
        } catch (netCDF::exceptions::NcException const& e) {                                             \
            std::cerr << "--Error--: " << __FILE__ << ":" << __LINE__ << " : " << e.what() << std::endl; \
            exit(EXIT_FAILURE);                                                                          \
        } catch (std::exception const& e) {                                                              \
            std::cerr << "--Error--:" << __FILE__ << ":" << __LINE__ << " : " << e.what() << std::endl;  \
            exit(EXIT_FAILURE);                                                                          \
        }                                                                                                \
    }
/**
 * @brief Find a netCDF variable by name in a file or group, recursively searching subgroups
 * @tparam T Type of netCDF file/group object
 * @tparam NcVarType Type of netCDF variable (defaults to netCDF::NcVar)
 * @param var_name Name of variable to find
 * @param nc_file netCDF file/group to search in
 * @return Found variable, or null variable if not found
 */
template <typename T, typename NcVarType = netCDF::NcVar>
NcVarType find_nc_variable_cpp(const std::string& var_name, T& nc_file) {
    NcVarType var;
    netCDF::NcVar temp_var = nc_file.getVar(var_name);
    if (temp_var.isNull()) {
        std::multimap<std::string, netCDF::NcGroup> grps = nc_file.getGroups();
        for (const auto& [_,grp] : grps) {
            temp_var = find_nc_variable_cpp(var_name, grp);
            if (!temp_var.isNull()) {
                return temp_var;
            }
        }
    } else {
        var = temp_var;
    }
    return var;
}

/**
 * @brief Find a netCDF group attribute by name, recursively searching subgroups
 * @tparam T Type of netCDF file/group object
 * @param att_name Name of attribute to find
 * @param nc_file netCDF file/group to search in
 * @return Found attribute, or null attribute if not found
 */
template <typename T>
netCDF::NcGroupAtt find_nc_grp_attribute_cpp(const std::string& att_name, T& nc_file) {
    netCDF::NcGroupAtt attr = nc_file.getAtt(att_name);
    if (attr.isNull()) {
        std::multimap<std::string, netCDF::NcGroup>  grps = nc_file.getGroups();
        for (const auto& [_,grp]  : grps) {
            attr = find_nc_grp_attribute_cpp(att_name, grp);
            if (!attr.isNull()) {
                return attr;
            }
        }
    }
    return attr;
};

/**
 * @brief Find all variables in a netCDF file/group and its subgroups
 * @tparam T Type of netCDF file/group object
 * @tparam NcVarType Type of netCDF variable (defaults to netCDF::NcVar)
 * @param nc_id netCDF file/group to search in
 * @return Multimap of variable names to variables
 */
template <typename T, typename NcVarType = netCDF::NcVar>
std::multimap<std::string, NcVarType> find_all_variables(T& nc_id) {
    std::multimap<std::string, netCDF::NcGroup> grps = nc_id.getGroups();
    std::multimap<std::string, NcVarType> vars = {};
    std::multimap<std::string, netCDF::NcVar> temp_var = nc_id.getVars();
    vars.insert(temp_var.begin(), temp_var.end());
    for (const auto& [_,grp] : grps) {
        temp_var = find_all_variables(grp);
        vars.insert(temp_var.begin(), temp_var.end());
    }
    return vars;
}

/**
 * @brief Find all variables with specified names in a netCDF file/group and its subgroups
 * @tparam T Type of netCDF file/group object
 * @tparam NcVarType Type of netCDF variable (defaults to netCDF::NcVar)
 * @param nc_id netCDF file/group to search in
 * @param names Set of variable names to search for
 * @return Multimap of matching variable names to variables
 */
template <typename T, typename NcVarType = netCDF::NcVar>
std::multimap<std::string, NcVarType> find_all_variables(T& nc_id,
                                                         const std::unordered_set<std::string>& names) {
    std::multimap<std::string, netCDF::NcGroup>  grps = nc_id.getGroups();
    std::multimap<std::string, NcVarType> vars{};
    std::multimap<std::string, netCDF::NcVar> temp_var = nc_id.getVars();
    for (const auto& nc_var : temp_var) {
        if (names.count(nc_var.first))
            vars.insert(nc_var);
    }
    for (const auto& [_,grp] : grps) {
        temp_var = find_all_variables(grp);
        for (const auto& nc_var : temp_var) {
            if (names.count(nc_var.first))
                vars.insert(nc_var);
        }
    }
    return vars;
}
#endif