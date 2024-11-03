#include "find_variable.hpp"

const std::unordered_set<std::string> non_data_names = {"l2_flags",
                                                        "latitude",
                                                        "Latitude",
                                                        "lat",
                                                        "Lat",
                                                        "longitude",
                                                        "Longitude",
                                                        "lon",
                                                        "Lon",
                                                        "solzen",
                                                        "solar_zenith_angle",
                                                        "solz",
                                                        "Solar Zenith Angle",
                                                        "csol_z",
                                                        "Center Solar Zenith Angle",
                                                        "center_solar_zenith_angle",
                                                        "satzen",
                                                        "senzen",
                                                        "sensor_zenith_angle",
                                                        "satellite_zenith_angle",
                                                        "satz",
                                                        "senz",
                                                        "Sensor Zenith Angle",
                                                        "Satellite Zenith Angle",
                                                        "senzen"};

std::multimap<std::string, netCDF::NcVar> find_variables_geo_physical(const std::string &nc_path,
                                                                      std::string &prod_list) {
    netCDF::NcFile file_;
    file_.open(nc_path, netCDF::NcFile::read);
    if (file_.isNull())
        throw std::runtime_error("--Error--: Could not open file " + nc_path);
    std::multimap<std::string, netCDF::NcVar> res = find_all_variables(file_, non_data_names);
    for (auto &var : res) {
        std::string grp_name = var.second.getParentGroup().getName();
        if(grp_name!="geophysical_data")
            continue;
        if (prod_list.empty())
            prod_list += var.first;
        else
            prod_list += "," + var.first;
    }
    return res;
}
