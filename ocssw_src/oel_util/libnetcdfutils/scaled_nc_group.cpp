/**
 *
 * @author Jakob C. Lindo (SSAI)
 * @brief Implementation of ScaledNcGroup
 * @date January 22 2025
 *
 */

#include "scaled_nc_group.hpp"
#include <type_traits>

using namespace netCDF;
using namespace std;

ScaledNcGroup::ScaledNcGroup() : NcGroup() {
}
ScaledNcGroup::ScaledNcGroup(const NcGroup &copied) : netCDF::NcGroup(copied) {
}
ScaledNcGroup::~ScaledNcGroup() {  // netCDF::NcGroup will handle destruction
}
multimap<string, ScaledNcVar> ScaledNcGroup::getVars(NcGroup::Location location) const {
    multimap<string, ScaledNcVar> scaledVars;
    multimap<string, NcVar> originalVars = NcGroup::getVars(location);

    for (const auto &pair : originalVars) {
        ScaledNcVar scaledVar(pair.second);
        scaledVars.insert({pair.first, scaledVar});
    }

    return scaledVars;
}

set<ScaledNcVar> ScaledNcGroup::getVars(const string &name, NcGroup::Location location) const {
    set<ScaledNcVar> scaledVars;
    set<NcVar> originalVars = netCDF::NcGroup::getVars(name, location);

    for (const NcVar &var : originalVars) {
        scaledVars.insert(ScaledNcVar(var));
    }

    return scaledVars;
}

ScaledNcVar ScaledNcGroup::getVar(const string &name, NcGroup::Location location) const {
    return ScaledNcVar(NcGroup::getVar(name, location));
}

ScaledNcVar ScaledNcGroup::addVar(const string &name, const netCDF::NcType &ncType) const {
    return ScaledNcVar(NcGroup::addVar(name, ncType));
}

ScaledNcVar ScaledNcGroup::addVar(const string &name, const string &typeName, const string &dimName) const {
    return ScaledNcVar(NcGroup::addVar(name, typeName, dimName));
}

ScaledNcVar ScaledNcGroup::addVar(const string &name, const netCDF::NcType &ncType,
                                  const netCDF::NcDim &ncDim) const {
    return ScaledNcVar(NcGroup::addVar(name, ncType, ncDim));
}

ScaledNcVar ScaledNcGroup::addVar(const string &name, const string &typeName,
                                  const vector<string> &dimNames) const {
    return ScaledNcVar(NcGroup::addVar(name, typeName, dimNames));
}

ScaledNcVar ScaledNcGroup::addVar(const string &name, const netCDF::NcType &ncType,
                                  const vector<netCDF::NcDim> &dims) const {
    return ScaledNcVar(NcGroup::addVar(name, ncType, dims));
}

ScaledNcGroup ScaledNcGroup::getGroup(const string &name, NcGroup::GroupLocation location) const {
    NcGroup originalGroup = NcGroup::getGroup(name, location);
    if (originalGroup.isNull()) {
        return ScaledNcGroup();
    }
    return ScaledNcGroup(originalGroup);
}
multimap<string, ScaledNcGroup> ScaledNcGroup::getGroups(NcGroup::GroupLocation location) const {
    multimap<string, ScaledNcGroup> scaledGroups;
    multimap<string, NcGroup> originalGroups = NcGroup::getGroups(location);

    for (const auto &pair : originalGroups) {
        ScaledNcGroup scaledGroup(pair.second);
        scaledGroups.insert({pair.first, scaledGroup});
    }

    return scaledGroups;
}

set<ScaledNcGroup> ScaledNcGroup::getGroups(const string &name, NcGroup::GroupLocation location) const {
    set<ScaledNcGroup> scaledGroups;
    set<NcGroup> originalGroups = NcGroup::getGroups(name, location);

    for (const NcGroup &group : originalGroups) {
        scaledGroups.insert(ScaledNcGroup(group));
    }

    return scaledGroups;
}
ScaledNcGroup ScaledNcGroup::getParentGroup() const {
    NcGroup parentGroup = NcGroup::getParentGroup();
    if (parentGroup.isNull()) {
        return ScaledNcGroup();
    }
    return ScaledNcGroup(parentGroup);
}

ScaledNcGroup ScaledNcGroup::addGroup(const string &name) const {
    NcGroup newGroup = NcGroup::addGroup(name);
    if (newGroup.isNull()) {
        return ScaledNcGroup();
    }
    return ScaledNcGroup(newGroup);
}