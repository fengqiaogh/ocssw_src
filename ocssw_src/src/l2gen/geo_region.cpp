#include "geo_region.h"
#include "l12_proto.h"
#include <netcdf>
#include <string>
#include <iostream>
#include <vector>

#define ERROR_EXIT                                                                         \
    {                                                                                      \
        std::cerr << "\n See " << __LINE__ << " in file " << __FILE__ << ". Exiting...\n"; \
        georegion_file.close();                                                              \
        exit(EXIT_FAILURE);                                                                \
    }
#define CHECK_BOUNDS(max_val, min_val)                                                                   \
    {                                                                                                    \
        if (min_val > max_val) {                                                                         \
            std::cerr << "--Error--: " << #min_val << "=" << min_val << " is larger than " << #max_val << "=" \
                      << max_val;                                                                        \
            ERROR_EXIT;                                                                                  \
        }                                                                                                \
    }
namespace {
netCDF::NcFile georegion_file;
netCDF::NcVar georegion_variable;
netCDF::NcVar lat_variable;
netCDF::NcVar lon_variable;
size_t nlat;
size_t nlon;
std::vector<double> latitude;
std::vector<double> longitude;
double step_lat;
double step_lon;
double min_lon;
double max_lon;
double min_lat;
double max_lat;
}  // namespace

void read_variable(const netCDF::NcFile &file, netCDF::NcVar &var, const std::string &varname) {
    var = file.getVar(varname);
    if (var.isNull()) {
        std::cerr << "--Error--: variable " << varname << " is not found in the georegion file";
        ERROR_EXIT;
    }
}

void close_georegion_file() {
    georegion_file.close();
}

int get_georegion(float lat, float lon) {
    if (georegion_variable.isNull()) {
        std::cout << "Loading the georegion file: " << input->georegionfile << std::endl;
        georegion_file.open(input->georegionfile, netCDF::NcFile::read);
        read_variable(georegion_file, georegion_variable, "georegion");
        read_variable(georegion_file, lat_variable, "lat");
        read_variable(georegion_file, lon_variable, "lon");
        std::vector<netCDF::NcDim> dims = georegion_variable.getDims();
        if (dims.size() != 2) {
            std::cerr << "--Error-- : The georegion variable should be a two dimensional variable.";
            ERROR_EXIT;
        }
        if (dims.at(0).getName() != "lat" || dims.at(1).getName() != "lon") {
            std::cerr << "--Error-- :  The first dimension of the georegion should be lat, the second one "
                         "should be lon. The dimensions found in the geofile are: first \""
                      << dims.at(0).getName() << "\" and second \"" << dims.at(1).getName() << "\"";
            ERROR_EXIT;
        }
        nlat = dims.at(0).getSize();
        nlon = dims.at(1).getSize();
        std::vector<netCDF::NcDim> dim_lat = lat_variable.getDims();
        std::vector<netCDF::NcDim> dim_lon = lon_variable.getDims();
        if (dim_lat.size() != 1) {
            std::cerr << "--Error-- : The latitude variable should be a one dimensional variable.";
            ERROR_EXIT;
        }
        if (dim_lon.size() != 1) {
            std::cerr << "--Error-- : The latitude variable should be a one dimensional variable.";
            ERROR_EXIT;
        }
        if (dim_lat.at(0).getName() != "lat" || dim_lon.at(0).getName() != "lon") {
            std::cerr << "--Error-- : The dimension mismatch between latitude/longitude and the georegion";
            ERROR_EXIT;
        }
        // checking for datatype
        if (netCDF::NcType::nc_BYTE != georegion_variable.getType().getTypeClass()) {
            std::cerr << "--Error-- : The datatype of the georegion should be signed byte";
            ERROR_EXIT;
        };
        if (netCDF::NcType::nc_DOUBLE != lat_variable.getType().getTypeClass()) {
            std::cerr << "--Error-- : The datatype of the latitude should be double";
            ERROR_EXIT;
        };
        if (netCDF::NcType::nc_DOUBLE != lon_variable.getType().getTypeClass()) {
            std::cerr << "--Error-- : The datatype of the longitude should be double";
            ERROR_EXIT;
        };
        latitude = std::vector<double>(nlat);
        longitude = std::vector<double>(nlon);
        lat_variable.getVar(latitude.data());
        lon_variable.getVar(longitude.data());
        step_lat = latitude[1] - latitude[0];
        step_lon = longitude[1] - longitude[0];
        if (step_lat > 0) {
            max_lat = latitude[nlat - 1];
            min_lat = latitude[0];
        } else {
            min_lat = latitude[nlat - 1];
            max_lat = latitude[0];
        }
        if (step_lon > 0) {
            max_lon = longitude[nlon - 1];
            min_lon = longitude[0];
        } else {
            min_lon = longitude[nlon - 1];
            max_lon = longitude[0];
        }
        netCDF::NcGroupAtt geospatial_lat_min = georegion_file.getAtt("geospatial_lat_min");
        netCDF::NcGroupAtt geospatial_lat_max = georegion_file.getAtt("geospatial_lat_max");
        netCDF::NcGroupAtt geospatial_lon_min = georegion_file.getAtt("geospatial_lon_min");
        netCDF::NcGroupAtt geospatial_lon_max = georegion_file.getAtt("geospatial_lon_max");
        if (!geospatial_lat_min.isNull()) {
            geospatial_lat_min.getValues(&min_lat);
        }
        if (!geospatial_lat_max.isNull()) {
            geospatial_lat_max.getValues(&max_lat);
        }
        if (!geospatial_lon_min.isNull()) {
            geospatial_lon_min.getValues(&min_lon);
        }
        if (!geospatial_lon_max.isNull()) {
            geospatial_lon_max.getValues(&max_lon);
        }
        CHECK_BOUNDS(max_lat,min_lat);
        CHECK_BOUNDS(max_lon,min_lon);
    }
    if (lat > max_lat || lat < min_lat)
        return 0;
    if (lon > max_lon || lon < min_lon)
        return 0;
    // should be uniform meshgrid
    size_t ilat = static_cast<size_t>(std::round((lat - latitude[0]) / step_lat));
    size_t ilon = static_cast<size_t>(std::round((lon - longitude[0]) / step_lon));
    // takes care of the edge cases

    std::vector<size_t> start = {ilat, ilon};
    std::vector<size_t> count = {1, 1};
    int8_t flag = 0;
    georegion_variable.getVar(start, count, &flag);
    //    if (flag == 0)
    //    {
    //        std::cout << ilat << " " <<  ilon << " " << lat << " " << lon << " " << latitude[0] << " " <<
    //        longitude[0] << std::endl;
    //    }
    return flag;
}
