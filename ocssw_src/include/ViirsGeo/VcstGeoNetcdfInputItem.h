/*******************************************************************************
 *
 * NAME: VcstGeoNetcdfInputItem
 *
 * DESCRIPTION: Subclass of the VcstGeoLutInputItem object class to support
 * reading netCDF LUT files.
 *
 * Created on: June 18, 2015
 *     Author: Sam Anderson, VCST
 *
 *
 *******************************************************************************/

#ifndef VcstGeoNetcdfInputItem_h
#define VcstGeoNetcdfInputItem_h

#include <VcstGeoLutInputItem.h>

#include <netcdf>

using namespace std;
using namespace netCDF;
using namespace netCDF::exceptions;

class VcstGeoNetcdfInputItem : public VcstGeoLutInputItem {
public:

    /**
     * Constructor
     *
     * @param lutName The LUT name that identifies the data buffer.
     */
    VcstGeoNetcdfInputItem(const std::string& groupName, size_t size);

    /**
     * Destructor
     */
    ~VcstGeoNetcdfInputItem();


    /**
     * Read data from LUT or netCDF file depending on availability
     *
     */
    int getData();


protected:

    /**
     * Read netCDF LUT file.
     */

    int readNetCDF();

    int read_GEO_PARAM(NcFile* nc_input);

};



#endif
