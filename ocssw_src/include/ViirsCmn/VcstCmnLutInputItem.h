/*******************************************************************************
 *
 * NAME: VcstCmnLutInputItem
 *
 * DESCRIPTION: Subclass of the VcstLutInputItem object class to support
 * reading netCDF LUT files.
 *
 * Created on: June 18, 2015
 *     Author: Sam Anderson, VCST
 *
 *
 *******************************************************************************/

#ifndef VcstCmnLutInputItem_h
#define VcstCmnLutInputItem_h

#include <netcdf>
#include <VcstLutInputItem.h>

using namespace std;
using namespace netCDF;
using namespace netCDF::exceptions;

class VcstCmnLutInputItem : public VcstLutInputItem {
public:

    /**
     * Class Constructor
     */
    VcstCmnLutInputItem(const std::string& groupName, size_t size);

    /**
     * Class Destructor
     */
    virtual ~VcstCmnLutInputItem();

    /**
     * Read data from LUT or netCDF file depending on availability
     */
    int getData();

protected:

    /**
     * Functions to read from netCDF LUT data used by Common Geolocation.
     */

    int readNetCDF();

    int read_CMNGEO_PLATFORM_LUT(NcFile* nc_input);

    int read_CMNGEO_JPL_EPHEM(NcFile* nc_input);

    int read_CMNGEO_SAA_COEFF(NcFile* nc_input);

    int read_CMNGEO_PARAM(NcFile* nc_input);

    int read_VIIRS_SOLAR_DIFF_VOLT_LUT(NcFile* nc_input);

    int read_VIIRS_SOLAR_DIFF_ROT_MATRIX_LUT(NcFile* nc_input);

    int read_SDR_QA_LUT(NcFile* nc_input);

    /**
     * Functions to convert endianness of LUT files used by Common Geolocation.
     */

    int convertEndianness();

    int convert_CMNGEO_JPL_EPHEM();

    int convert_CMNGEO_SAA_COEFF();

    int convert_CMNGEO_PARAM();

    int convert_VIIRS_SOLAR_DIFF_VOLT_LUT();

    int convert_VIIRS_SOLAR_DIFF_ROT_MATRIX_LUT();

    int convert_SDR_QA_LUT();
};

#endif
