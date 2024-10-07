/*******************************************************************************
 *
 * NAME: VcstCalNetcdfInputItem.h
 *
 * DESCRIPTION: Subclass of the VcstCalLutInputItem object class to support
 * reading a netCDF LUT file.
 *
 * Created on: June 18, 2015
 *     Author: Sam Anderson, VCST
 *
 *  Modified: October 15, 2015, for generation and processing of time-dependent
 *  netCDF4 LUTs
 *
 *******************************************************************************/

#ifndef VcstCalNetcdfInputItem_h
#define VcstCalNetcdfInputItem_h

#include <VcstCalLutInputItem.h>
#include <VcstCalLutStructures.h>

#include <netcdf>

using namespace std;
using namespace netCDF;
using namespace netCDF::exceptions;

class VcstCalNetcdfInputItem : public VcstCalLutInputItem {
public:

    /**
     * Constructor
     *
     */

    VcstCalNetcdfInputItem(const std::string& groupName, size_t size);

    /**
     * Destructor
     *
     */

    virtual ~VcstCalNetcdfInputItem();

    /**
     * Read data from LUT or netCDF file depending on availability
     *
     */

    virtual int getData();

    /**
     * Read data from netCDF time-dependent dynamic LUT file
     *
     */

    int getData(double tai);


protected:

    /**
     * Functions to from netCDF LUT data used by Common Geolocation.
     */

    int readNetCDF();

    int read_DG_ANOMALY_DN_LIMITS_LUT(NcFile* nc_input);
    int read_DNB_STRAY_LIGHT_CORRECTION_LUT(NcFile* nc_input);
    int read_SDR_DNB_DN0_LUT(NcFile* nc_input);
    int read_SDR_DNB_RVS(NcFile* nc_input);
    int read_SDR_DNB_FRAME_TO_ZONE(NcFile* nc_input);
    int read_SDR_DNB_LGS_GAINS_LUT(NcFile* nc_input);
    int read_SDR_DNB_GAIN_RATIOS_LUT(NcFile* nc_input);
    int read_SDR_F_PREDICTED_LUT(NcFile* nc_input);
    int read_SDR_GAIN_LUT(NcFile* nc_input);
    int read_SDR_HAM_ER_TABLE(NcFile* nc_input);
    int read_SDR_RTA_ER_TABLE(NcFile* nc_input);
    int read_SDR_OBC_ER_TABLE(NcFile* nc_input);
    int read_SDR_OBC_RR_TABLE(NcFile* nc_input);
    int read_SDR_EBBT_TABLE(NcFile* nc_input);
    int read_SDR_TELE_COEFFS(NcFile* nc_input);
    int read_SDR_SOLAR_IRAD_LUT(NcFile* nc_input);
    int read_SDR_OBS_TO_PIXELS(NcFile* nc_input);
    int read_SDR_RADIOMETRIC_PARAMETERS(NcFile* nc_input);
    int read_SDR_EMISSIVE_LUT(NcFile* nc_input);
    int read_SDR_REFLECTIVE_LUT(NcFile* nc_input);
    int read_SDR_RVS_LUT(NcFile* nc_input);
    int read_SDR_BB_TEMP_COEFFS(NcFile* nc_input);
    int read_SDR_DELTA_C_LUT(NcFile* nc_input);
    int read_SDR_COEFF_A_LUT(NcFile* nc_input);
    int read_SDR_COEFF_B_LUT(NcFile* nc_input);
    int read_SDR_RELATIVE_SPECTRAL_RESPONSE_LUT(NcFile* nc_input);
    int read_SDR_SOLAR_SPECTRAL_IRAD_LUT(NcFile* nc_input);


    /**
     * expand_StrayLight_LUT()
     *
     * Expand DNB Stray Light LUT through linear interpolation of values.
     *
     * @return FAIL or SUCCESS
     */

    int expand_StrayLight_LUT(shortViirsCalDnbStrayLightLUT* short_lut,
            proSdrViirsCalDnbStrayLightLUT* long_lut);


    /**
     * Functions to from time-dependent reprocess LUT data used by Calibration.
     */

    int readDynamicNetcdf();

    int read_Dynamic_DNB_STRAY_LIGHT_CORRECTION_LUT(NcFile* nc_input);
    int read_Dynamic_SDR_DNB_DN0_LUT(NcFile* nc_input);
    int read_Dynamic_SDR_DNB_LGS_GAINS_LUT(NcFile* nc_input);
    int read_Dynamic_SDR_DNB_GAIN_RATIOS_LUT(NcFile* nc_input);
    int read_Dynamic_SDR_F_PREDICTED_LUT(NcFile* nc_input);
    //    int read_Dynamic_SDR_RELATIVE_SPECTRAL_RESPONSE_LUT( NcFile* nc_input );
    int read_Dynamic_SDR_SOLAR_SPECTRAL_IRAD_LUT(NcFile* nc_input);

    long long iet_;

};



#endif
