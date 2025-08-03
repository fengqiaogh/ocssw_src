/**
 *
 * @author Jakob C. Lindo (SSAI)
 * @date January 2025
 * @brief An extension of netCDF::NcGroup to help integrate some functions for use with ScaledNcVar
 *
 */

#ifndef SCALED_NCGROUP_H
#define SCALED_NCGROUP_H

#include <netcdf>
#include "scaled_nc_var.hpp"

class ScaledNcGroup : public netCDF::NcGroup {
   public:
    /**
     * @brief Constructor. Generates a null object
     */
    ScaledNcGroup();

    /**
     * @brief Create a ScaledNcGroup out of an existing netCDF::NcGroup
     * @param copied The group to be copied
     */
    ScaledNcGroup(const netCDF::NcGroup &copied);

    /**
     * @brief ScaledNcGroup destructor
     */
    ~ScaledNcGroup();

    std::multimap<std::string, ScaledNcVar> getVars(NcGroup::Location location = Current) const;
    std::set<ScaledNcVar> getVars(const std::string &name, NcGroup::Location location = Current) const;

    ScaledNcVar getVar(const std::string &name, NcGroup::Location location = Current) const;
    ScaledNcVar addVar(const std::string &name, const netCDF::NcType &ncType) const;
    ScaledNcVar addVar(const std::string &name, const std::string &typeName,
                       const std::string &dimName) const;
    ScaledNcVar addVar(const std::string &name, const netCDF::NcType &ncType,
                       const netCDF::NcDim &ncDim) const;
    ScaledNcVar addVar(const std::string &name, const std::string &typeName,
                       const std::vector<std::string> &dimNames) const;
    ScaledNcVar addVar(const std::string &name, const netCDF::NcType &ncType,
                       const std::vector<netCDF::NcDim> &dims) const;

    ScaledNcGroup getGroup(const std::string &name, NcGroup::GroupLocation location = ChildrenGrps) const;
    std::multimap<std::string, ScaledNcGroup> getGroups(NcGroup::GroupLocation location = ChildrenGrps) const;
    std::set<ScaledNcGroup> getGroups(const std::string &name,
                                      NcGroup::GroupLocation location = ChildrenGrps) const;
    ScaledNcGroup getParentGroup() const;
    ScaledNcGroup addGroup(const std::string &name) const;
};

#endif
