/**
 *
 * @author Jakob C. Lindo (SSAI)
 * @date January 2025
 * @brief An extension of netCDF::NcFile to help integrate some functions for use with ScaledNcVar
 *
 */

#ifndef SCALED_NCFILE_H
#define SCALED_NCFILE_H

#include <netcdf>
#include "scaled_nc_group.hpp"

class ScaledNcFile : public netCDF::NcFile {
   public:
    /**
     * @brief Constructor. Generates a null object
     */
    ScaledNcFile();

    /**
     * @brief Create a ScaledNcFile out of an existing netCDF::NcFile
     * @param copied The file to be copied
     */
    ScaledNcFile(const netCDF::NcFile &existingFile, netCDF::NcFile::FileMode openMode);

    /**
     * @brief Constructor that opens an existing netCDF file
     * @param fileName The name of the file to open
     * @param openMode The mode in which to open the file
     */
    ScaledNcFile(const std::string fileName, netCDF::NcFile::FileMode openMode);

    /**
     * @brief Destructor. Closes the file if open
     */
    ~ScaledNcFile();

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

    ScaledNcGroup getParentGroup() const;
    std::multimap<std::string, ScaledNcGroup> getGroups(NcGroup::GroupLocation location = ChildrenGrps) const;
    std::set<ScaledNcGroup> getGroups(const std::string &name,
                                      NcGroup::GroupLocation location = ChildrenGrps) const;
    ScaledNcGroup getGroup(const std::string &name, NcGroup::GroupLocation location = ChildrenGrps) const;
    ScaledNcGroup addGroup(const std::string &name) const;
};

#endif
