/**
 *
 * @author Jakob C. Lindo (SSAI)
 * @brief Implementation of ScaledNcFile
 * @date January 22 2025
 *
 */

#include "scaled_nc_file.hpp"
#include <type_traits>

using namespace netCDF;
using namespace std;

ScaledNcFile::ScaledNcFile() : netCDF::NcFile() {
}

ScaledNcFile::ScaledNcFile(const string fileName, NcFile::FileMode openMode) : NcFile(fileName, openMode) {
}
ScaledNcFile::~ScaledNcFile() {
}

multimap<string, ScaledNcVar> ScaledNcFile::getVars(NcGroup::Location location) const {
    multimap<string, ScaledNcVar> scaledVars;
    multimap<string, NcVar> originalVars = NcFile::getVars(location);

    for (const auto &var : originalVars) {
        scaledVars.emplace(var.first, ScaledNcVar(var.second));
    }

    return scaledVars;
}

set<ScaledNcVar> ScaledNcFile::getVars(const string &name, NcGroup::Location location) const {
    set<ScaledNcVar> scaledVars;
    set<NcVar> originalVars = NcFile::getVars(name, location);

    for (const auto &var : originalVars) {
        scaledVars.emplace(ScaledNcVar(var));
    }

    return scaledVars;
}

ScaledNcVar ScaledNcFile::getVar(const string &name, NcGroup::Location location) const {
    return ScaledNcVar(NcFile::getVar(name, location));
}

ScaledNcVar ScaledNcFile::addVar(const string &name, const netCDF::NcType &ncType) const {
    return ScaledNcVar(NcFile::addVar(name, ncType));
}

ScaledNcVar ScaledNcFile::addVar(const string &name, const string &typeName, const string &dimName) const {
    return ScaledNcVar(NcFile::addVar(name, typeName, dimName));
}

ScaledNcVar ScaledNcFile::addVar(const string &name, const netCDF::NcType &ncType,
                                 const netCDF::NcDim &ncDim) const {
    return ScaledNcVar(NcFile::addVar(name, ncType, ncDim));
}

ScaledNcVar ScaledNcFile::addVar(const string &name, const string &typeName,
                                 const vector<string> &dimNames) const {
    return ScaledNcVar(NcFile::addVar(name, typeName, dimNames));
}

ScaledNcVar ScaledNcFile::addVar(const string &name, const netCDF::NcType &ncType,
                                 const vector<netCDF::NcDim> &dims) const {
    return ScaledNcVar(NcFile::addVar(name, ncType, dims));
}

ScaledNcGroup ScaledNcFile::getGroup(const string &name, NcGroup::GroupLocation location) const {
    return ScaledNcGroup(NcFile::getGroup(name, location));
}
multimap<string, ScaledNcGroup> ScaledNcFile::getGroups(NcGroup::GroupLocation location) const {
    multimap<string, ScaledNcGroup> scaledGroups;
    multimap<string, NcGroup> originalGroups = NcFile::getGroups(location);

    for (const auto &group : originalGroups) {
        scaledGroups.emplace(group.first, ScaledNcGroup(group.second));
    }

    return scaledGroups;
}

set<ScaledNcGroup> ScaledNcFile::getGroups(const string &name, NcGroup::GroupLocation location) const {
    set<ScaledNcGroup> scaledGroups;
    set<NcGroup> originalGroups = NcFile::getGroups(name, location);

    for (const auto &group : originalGroups) {
        scaledGroups.emplace(ScaledNcGroup(group));
    }

    return scaledGroups;
}
ScaledNcGroup ScaledNcFile::getParentGroup() const {
    return ScaledNcGroup(NcFile::getParentGroup());
}

ScaledNcGroup ScaledNcFile::addGroup(const string &name) const {
    return ScaledNcGroup(NcFile::addGroup(name));
}
