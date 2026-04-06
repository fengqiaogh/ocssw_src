#ifndef NCFILE_UTILS_HPP
#define NCFILE_UTILS_HPP

#include <iostream>
#include <iterator>
#include <netcdf>
#include <string>
#include <vector>

#define lineinfo std::cerr << __FILE__ << ":" << __LINE__ << " (" << __FUNCTION__ << ")\n"

// NetCDF general utils

/// Group attributes

bool hasAtt(const netCDF::NcGroup grp, const std::string& name) {
    auto attlist = grp.getAtts();
    auto iter = attlist.find(name);
    return (iter != attlist.end());
}
// read scalar attribute
template <typename T>
void readAtt(const netCDF::NcGroup grp, const std::string& name, T* attval) {
    if (hasAtt(grp, name))
        grp.getAtt(name).getValues(attval);
}
// read attribute array into vector
template <typename T>
void readAtt(const netCDF::NcGroup grp, const std::string& name, std::vector<T>* attvals) {
    if (hasAtt(grp, name)) {
        netCDF::NcGroupAtt att = grp.getAtt(name);
        size_t attlen = att.getAttLength();
        attvals->resize(attlen);
        att.getValues(attvals->data());
    }
}

/// Variable attributes

bool hasAtt(const netCDF::NcVar var, const std::string& name) {
    auto attlist = var.getAtts();
    auto iter = attlist.find(name);
    return (iter != attlist.end());
}
// read scalar attribute
template <typename T>
void readAtt(const netCDF::NcVar var, const std::string& name, T* attval) {
    if (hasAtt(var, name))
        var.getAtt(name).getValues(attval);
}
// read attribute array into vector
template <typename T>
void readAtt(const netCDF::NcVar var, const std::string& name, std::vector<T>* attvals) {
    if (hasAtt(var, name)) {
        netCDF::NcVarAtt att = var.getAtt(name);
        size_t attlen = att.getAttLength();
        attvals->resize(attlen);
        att.getValues(attvals->data());
    }
}

/// get scaling values, if they exist
template <typename T>
struct Scaling {
    bool is_scaled;
    T scale;
    T offset;
};
template <typename T>
Scaling<T> varscaling(const netCDF::NcVar var) {
    Scaling<T> s;
    s.scale = 1;
    s.offset = 0;
    s.is_scaled = hasAtt(var, "scale_factor") || hasAtt(var, "add_offset");
    if (s.is_scaled) {
        readAtt(var, "scale_factor", &s.scale);
        readAtt(var, "add_offset", &s.offset);
    }
    return s;
}

/// read data into scalar
template <typename T>
void readVar(const netCDF::NcVar var, T scalar) {
    try {
        // read data
        var.getVar(&scalar);

        // apply traditional (y=mx+b) scaling as needed
        auto varscales = varscaling<T>(var);
        if (varscales.is_scaled) {
            scalar *= varscales.scale;
            scalar += varscales.offset;
        }
    } catch (netCDF::exceptions::NcException& e) {
        lineinfo;
        std::cerr << e.what() << std::endl;
        throw;
    }
}

/// read data into pointer
template <typename T>
void readVar(const netCDF::NcVar var, T* dataptr) {
    try {
        // read data
        var.getVar(&dataptr[0]);

        // apply traditional (y=mx+b) scaling as needed
        auto varscales = varscaling<T>(var);
        if (varscales.is_scaled) {
            size_t nvals = 1;
            for (auto& dim : var.getDims())
                nvals *= dim.getSize();
            for (size_t i = 0; i != nvals; i++) {
                dataptr[i] *= varscales.scale;
                dataptr[i] += varscales.offset;
            }
        }
    } catch (netCDF::exceptions::NcException& e) {
        lineinfo;
        std::cerr << e.what() << std::endl;
        throw;
    }
}

/// read data into vector
template <typename T>
void readVar(const netCDF::NcVar var, std::vector<T>& datavec) {
    try {
        // allocate memory as needed
        if (datavec.empty()) {
            size_t nvals = 1;
            for (auto& dim : var.getDims())
                nvals *= dim.getSize();
            datavec.resize(nvals);
        }

        // read data
        var.getVar(datavec.data());

        // apply traditional (y=mx+b) scaling as needed
        auto varscales = varscaling<T>(var);
        if (varscales.is_scaled) {
            for (auto& val : datavec) {
                val *= varscales.scale;
                val += varscales.offset;
            }
        }
    } catch (netCDF::exceptions::NcException& e) {
        lineinfo;
        std::cerr << e.what() << std::endl;
        throw;
    }
}

/// find variable, then read data
template <typename T>
void readVar(const netCDF::NcGroup grp, const std::string& name, T& data) {
    try {
        netCDF::NcVar var = grp.getVar(name, netCDF::NcGroup::ChildrenAndCurrent);
        readVar(var, data);
    } catch (netCDF::exceptions::NcException& e) {
        lineinfo;
        std::cerr << e.what() << std::endl;
        throw;
    }
}

#endif
