/**
 * @name l1b_file.cpp
 *
 * @brief Functions definitions for a Level 1B (L1B or l1b) file.
 *
 * @authors Joel Gales (SAIC), Jakob C Lindo (SSAI)
 * @date Aug 2024
 */

#include "l1b_file.hpp"
#include <genutils.h>
#include <netcdf>
#include <nc4utils.h>
#include <sstream>
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/multi_array.hpp>
#include <allocate3d.h>

using namespace std;
using namespace netCDF;
using namespace netCDF::exceptions;

int DEFLATE_LEVEL = 5;

int Level1bFile::createFile(size_t numScans, size_t numBlueBands, size_t numRedBands, size_t numHyperSciPix,
                            size_t numSwirPix, size_t numSwirBands, int32_t *pprOffset,
                            bool radianceGenerationEnabled, LocatingContext locatingContext,
                            int deflateLevel) {

    if (deflateLevel < 0 || 9 < deflateLevel) { // Should never happen; should be checked by L1bOptions
        cerr << "Invalid deflate level: " << to_string(deflateLevel) << endl;
        exit(EXIT_FAILURE);
    }

    DEFLATE_LEVEL = deflateLevel;

    // Put global dimensions
    map<string, size_t> dims;
    if (numHyperSciPix == numSwirPix) {
        dims = {{"scans", numScans},          {"pixels", numHyperSciPix},       {"red_bands", numRedBands},
                {"blue_bands", numBlueBands}, {"SWIR_bands", numSwirBands},     {"vector_elements", 3},
                {"quaternion_elements", 4},   {"polarization_coefficients", 3}, {"HAM_sides", 2}};
    } else {
        dims = {{"scans", numScans},
                {"ccd_pixels", numHyperSciPix},
                {"SWIR_pixels", numSwirPix},
                {"red_bands", numRedBands},
                {"blue_bands", numBlueBands},
                {"SWIR_bands", numSwirBands},
                {"vector_elements", 3},
                {"quaternion_elements", 4},
                {"polarization_coefficients", 3},
                {"HAM_sides", 2}};
    }

    for (auto dim : dims) {
        try {
            ncDims[numDimensions++] = l1bFile->addDim(dim.first, dim.second);
        } catch (NcException const &e) {
            cerr << "-E- Failed to create NcDim " << dim.first << ". NetCDF error " << e.what() << endl;
            return EXIT_FAILURE;
        }
    }

    // Put global attributes
    try {
        l1bFile->putAtt("title", "PACE OCI Level-1B Data");
        l1bFile->putAtt("instrument", "OCI");
        l1bFile->putAtt("platform", "PACE");
        l1bFile->putAtt("product_name", "PACE_OCI_yyyymmddThhmmss.L1B.nc");
        l1bFile->putAtt("processing_level", "L1B");
        l1bFile->putAtt("cdm_data_type", "swath");
        l1bFile->putAtt("geospatial_lat_units", "degrees_north");
        l1bFile->putAtt("geospatial_lon_units", "degrees_east");
        l1bFile->putAtt("geospatial_lat_min", NC_FLOAT, BAD_FLT);
        l1bFile->putAtt("geospatial_lat_max", NC_FLOAT, BAD_FLT);
        l1bFile->putAtt("geospatial_lon_min", NC_FLOAT, BAD_FLT);
        l1bFile->putAtt("geospatial_lon_max", NC_FLOAT, BAD_FLT);
        l1bFile->putAtt("rta_nadir", NC_LONG, 2, pprOffset);
    } catch (NcException const &e) {
        cerr << "-E- Failed to create global attributes. NetCDF Error " << e.what() << endl;
        return EXIT_FAILURE;
    }

    // Put top-level groups
    try {
        sensorBandParameters = l1bFile->addGroup("sensor_band_parameters");
        scanLineAttributes = l1bFile->addGroup("scan_line_attributes");
        geolocationData = l1bFile->addGroup("geolocation_data");
        navigationData = l1bFile->addGroup("navigation_data");
        observationData = l1bFile->addGroup("observation_data");
    } catch (NcException const &e) {
        cerr << "-E- Failed to create top-level groups. NetCDF Error " << e.what() << endl;
        return EXIT_FAILURE;
    }

    // Put variables in groups
    double constexpr fillDouble = BAD_FLT;

    // sensor_band_parameters
    vector<NcDim> dimVec{l1bFile->getDim("blue_bands")};
    createField(sensorBandParameters, "blue_wavelength", "Band center wavelengths for bands from blue CCD",
                NULL, "nm", "", fillDouble, "", "", 305.f, 610.f, 1.0, 0.0, NC_FLOAT, dimVec, "");
    createField(sensorBandParameters, "blue_solar_irradiance",
                "Mean extraterrestrial solar irradiance at 1 astronomical unit for the wavelengths of "
                "the blue CCD",
                NULL, "W m^-2 um^-1", "", fillDouble, "", "", 0.f, 2500.f, 1.0, 0.0, NC_FLOAT, dimVec, "");
    dimVec = {l1bFile->getDim("red_bands")};
    createField(sensorBandParameters, "red_wavelength", "Band center wavelengths for bands from red CCD",
                NULL, "nm", "", fillDouble, "", "", 595.f, 900.f, 1.0, 0.0, NC_FLOAT, dimVec, "");
    createField(sensorBandParameters, "red_solar_irradiance",
                "Mean extraterrestrial solar irradiance at 1 astronomical unit for the wavelengths of "
                "the red CCD",
                NULL, "W m^-2 um^-1", "", fillDouble, "", "", 0.f, 2500.f, 1.0, 0.0, NC_FLOAT, dimVec, "");
    dimVec = {l1bFile->getDim("SWIR_bands")};
    createField(sensorBandParameters, "SWIR_wavelength", "Band center wavelengths for SWIR bands", NULL, "nm",
                "", fillDouble, "", "", 900.f, 2260.f, 1.0, 0.0, NC_FLOAT, dimVec, "");
    createField(sensorBandParameters, "SWIR_bandpass", "Bandpasses for SWIR bands", NULL, "nm", "",
                fillDouble, "", "", 0.f, 100.f, 1.0, 0.0, NC_FLOAT, dimVec, "");
    createField(sensorBandParameters, "SWIR_solar_irradiance",
                "Mean extraterrestrial solar irradiance at 1 astronomical unit for the SWIR "
                "wavelengths",
                NULL, "W m^-2 um^-1", "", fillDouble, "", "", 0.f, 2500.f, 1.0, 0.0, NC_FLOAT, dimVec, "");
    dimVec = {l1bFile->getDim("blue_bands"), l1bFile->getDim("HAM_sides"),
              l1bFile->getDim("polarization_coefficients")};
    createField(sensorBandParameters, "blue_m12_coef", "Blue band M12/M11 polynomial coefficients", NULL,
                NULL, "", fillDouble, "", "", fillDouble, fillDouble, 1.0, 0.0, NC_FLOAT, dimVec, "");
    createField(sensorBandParameters, "blue_m13_coef", "Blue band M13/M11 polynomial coefficients", NULL,
                NULL, "", fillDouble, "", "", fillDouble, fillDouble, 1.0, 0.0, NC_FLOAT, dimVec, "");
    dimVec = {l1bFile->getDim("red_bands"), l1bFile->getDim("HAM_sides"),
              l1bFile->getDim("polarization_coefficients")};
    createField(sensorBandParameters, "red_m12_coef", "Red band m12/M11 polynomial coefficients", NULL, NULL,
                "", fillDouble, "", "", fillDouble, fillDouble, 1.0, 0.0, NC_FLOAT, dimVec, "");
    createField(sensorBandParameters, "red_m13_coef", "Red band M13/M11 polynomial coefficients", NULL, NULL,
                "", fillDouble, "", "", fillDouble, fillDouble, 1.0, 0.0, NC_FLOAT, dimVec, "");
    dimVec = {l1bFile->getDim("SWIR_bands"), l1bFile->getDim("HAM_sides"),
              l1bFile->getDim("polarization_coefficients")};
    createField(sensorBandParameters, "SWIR_m12_coef", "SWIR band M12/M11 polynomial coefficients", NULL,
                NULL, "", fillDouble, "", "", fillDouble, fillDouble, 1.0, 0.0, NC_FLOAT, dimVec, "");
    createField(sensorBandParameters, "SWIR_m13_coef", "SWIR band M13/M11 polynomial coefficients", NULL,
                NULL, "", fillDouble, "", "", fillDouble, fillDouble, 1.0, 0.0, NC_FLOAT, dimVec, "");

    // scan_line_attributes
    vector<NcDim> scanLineAttrs{l1bFile->getDim("scans")};
    createField(scanLineAttributes, "time", "time", NULL, "seconds since YYYY-MM-DD 00:00:00Z",
                "Earth view mid time in seconds of day", fillDouble, "", "", 0., 172802., 1.0, 0.0, NC_DOUBLE,
                scanLineAttrs, "");
    createField(scanLineAttributes, "HAM_side", "Half-angle mirror side", NULL, NULL, NULL, BAD_UBYTE, "", "",
                0, 1, 1.0, 0.0, NC_UBYTE, scanLineAttrs, "");
    createField(scanLineAttributes, "scan_quality_flags", "Scan quality flags ", NULL, NULL, NULL, BAD_UBYTE,
                "1UB, 2UB, 4UB", "tilt_change missing_time missing_encoder", 0, 0, 1.0, 0.0, NC_UBYTE,
                scanLineAttrs, "");

    vector<NcDim> scansByPixels;
    if (numHyperSciPix == numSwirPix)
        scansByPixels = {l1bFile->getDim("scans"), l1bFile->getDim("pixels")};
    else
        scansByPixels = {l1bFile->getDim("scans"), l1bFile->getDim("ccd_pixels")};

    string latString = "Latitudes of pixel locations";
    string lonString = "Longitudes of pixel locations";
    if (locatingContext == SELENO) {
        latString = "Selenographic latitudes of pixel locations";
        lonString = "Selenographic longitudes of pixel locations";
    }

    createField(geolocationData, "latitude", latString.c_str(), NULL, "degrees_north", NULL, fillDouble, "",
                "", -90.f, 90.f, 1.0, 0.0, NC_FLOAT, scansByPixels, "");
    createField(geolocationData, "longitude", lonString.c_str(), NULL, "degrees_east", NULL, fillDouble, "",
                "", -180.f, 180.f, 1.0, 0.0, NC_FLOAT, scansByPixels, "");
    createField(geolocationData, "height", "Terrain height at pixel locations", NULL, "meters", NULL,
                fillDouble, "", "", static_cast<short>(-10000), static_cast<short>(10000), 1.0, 0.0, NC_SHORT,
                scansByPixels, "longitude latitude");
    createField(geolocationData, "watermask", "Watermask", NULL, "", "0 for land, 1 for water", BAD_UBYTE, "",
                "", 0, 1, 1, 0, NC_UBYTE, scansByPixels, "longitude latitude");
    createField(geolocationData, "sensor_azimuth", "Sensor azimuth angle at pixel locations", NULL, "degrees",
                NULL, fillDouble, "", "", -18000, 18000, 0.01, 0.f, NC_SHORT, scansByPixels,
                "longitude latitude");
    createField(geolocationData, "sensor_zenith", "Sensor zenith angle at pixel locations", NULL, "degrees",
                NULL, fillDouble, "", "", 0, 18000, 0.01, 0.f, NC_SHORT, scansByPixels, "longitude latitude");
    createField(geolocationData, "solar_azimuth", "Solar azimuth angle at pixel locations", NULL, "degrees",
                NULL, fillDouble, "", "", -18000, 18000, 0.01, 0.f, NC_SHORT, scansByPixels,
                "longitude latitude");
    createField(geolocationData, "solar_zenith", "Solar zenith angle at pixel locations", NULL, "degrees",
                NULL, fillDouble, "", "", 0, 18000, 0.01, 0.f, NC_SHORT, scansByPixels, "longitude latitude");

    // Uncomment these and delete above when time comes to use floats in the output instead of shorts
    // createField(geolocationData, "sensor_azimuth", "Sensor azimuth angle at pixel locations", NULL,
    // "degrees",
    //             NULL, fillDouble, "", "", -180, 180, 1.0, 0.f, NC_FLOAT, scansByPixels,
    //             "longitude latitude");
    // createField(geolocationData, "sensor_zenith", "Sensor zenith angle at pixel locations", NULL,
    // "degrees",
    //             NULL, fillDouble, "", "", 0, 180, 1.0, 0.f, NC_FLOAT, scansByPixels, "longitude latitude");
    // createField(geolocationData, "solar_azimuth", "Solar azimuth angle at pixel locations", NULL,
    // "degrees",
    //             NULL, fillDouble, "", "", -180, 180, 1.0, 0.f, NC_FLOAT, scansByPixels,
    //             "longitude latitude");
    // createField(geolocationData, "solar_zenith", "Solar zenith angle at pixel locations", NULL, "degrees",
    //             NULL, fillDouble, "", "", 0, 180, 1.0, 0.f, NC_FLOAT, scansByPixels, "longitude latitude");

    createField(geolocationData, "quality_flag", "Geolocation pixel quality flags", NULL, NULL, NULL,
                BAD_UBYTE, "1UB, 2UB, 4UB", "Off_Earth Solar_eclipse Terrain_bad", 0, 0, 1.0, 0.0, NC_UBYTE,
                scansByPixels, "longitude latitude");

    vector<NcDim> scansVsQuatElems{l1bFile->getDim("scans"), l1bFile->getDim("quaternion_elements")};
    vector<NcDim> scansVsVecElems{l1bFile->getDim("scans"), l1bFile->getDim("vector_elements")};
    createField(navigationData, "att_quat", "Attitude quaternions at EV mid-times", NULL, NULL, NULL,
                fillDouble, "", "", -1.f, 1.f, 1.0, 0.0, NC_FLOAT, scansVsQuatElems, "");
    createField(navigationData, "att_ang", "Attitude angles (roll, pitch, yaw) at EV mid-times", NULL,
                "degrees", NULL, fillDouble, "", "", -180.f, 180.f, 1.0, 0.0, NC_FLOAT, scansVsVecElems, "");

    double orbFill = -9999999.0;
    createField(navigationData, "orb_pos", "Orbit position vectors at EV mid-times (ECR)", NULL, "kilometers",
                NULL, orbFill, "", "", -7100.f, 7100.f, 1.0, 0.0, NC_FLOAT, scansVsVecElems, "");
    createField(navigationData, "orb_vel", "Orbit velocity vectors at EV mid-times (ECR)", NULL,
                "kilometers/second", NULL, fillDouble, "", "", -7.6f, 7.6f, 1.0, 0.0, NC_FLOAT,
                scansVsVecElems, "");
    createField(navigationData, "sun_ref", "Solar unit vectors in J2000 frame", NULL, NULL, NULL, fillDouble,
                "", "", -1.f, 1.f, 1.0, 0.0, NC_FLOAT, scansVsVecElems, "");
    createField(navigationData, "tilt_angle", "OCI tilt angles at EV mid-times", NULL, "degrees", NULL,
                fillDouble, "", "", -22.5f, 22.5f, 1.0, 0.0, NC_FLOAT, scanLineAttrs, "");
    createField(navigationData, "CCD_scan_angles", "Scan angles for blue and red band science pixels", NULL,
                "degrees", NULL, fillDouble, "", "", -110.f, 250.f, 1.0, 0.0, NC_FLOAT, scansByPixels, "");

    vector<NcDim> scansVsSwirPix;
    if (numHyperSciPix == numSwirPix)
        scansVsSwirPix = {l1bFile->getDim("scans"), l1bFile->getDim("pixels")};
    else
        scansVsSwirPix = {l1bFile->getDim("scans"), l1bFile->getDim("SWIR_pixels")};

    createField(navigationData, "SWIR_scan_angles", "Scan angles for SWIR band science pixels", NULL,
                "degrees", NULL, fillDouble, "", "", -110.f, 250.f, 1.0, 0.0, NC_FLOAT, scansVsSwirPix, "");

    vector<NcDim> bluePicture;
    vector<NcDim> redPicture;
    vector<NcDim> swirPicture;
    if (numHyperSciPix == numSwirPix) {
        bluePicture = {l1bFile->getDim("blue_bands"), l1bFile->getDim("scans"), l1bFile->getDim("pixels")};
        redPicture = {l1bFile->getDim("red_bands"), l1bFile->getDim("scans"), l1bFile->getDim("pixels")};
        swirPicture = {l1bFile->getDim("SWIR_bands"), l1bFile->getDim("scans"), l1bFile->getDim("pixels")};
    } else if (numHyperSciPix != numSwirPix) {
        bluePicture = {l1bFile->getDim("blue_bands"), l1bFile->getDim("scans"),
                       l1bFile->getDim("ccd_pixels")};
        redPicture = {l1bFile->getDim("red_bands"), l1bFile->getDim("scans"), l1bFile->getDim("ccd_pixels")};
        swirPicture = {l1bFile->getDim("SWIR_bands"), l1bFile->getDim("scans"),
                       l1bFile->getDim("SWIR_pixels")};
    }

    if (radianceGenerationEnabled) {
        createField(observationData, "Lt_blue", "Top of Atmosphere Blue Band Radiance", NULL, "W m-2 sr-1",
                    "", BAD_FLT, "", "", 0.0f, 800.0f, 1.0, 0.0, NC_FLOAT, bluePicture, "longitude latitude");
        createField(observationData, "Lt_red", "Top of Atmosphere Red Band Radiance", NULL, "W m-2 sr-1", "",
                    BAD_FLT, "", "", 0.0f, 700.0f, 1.0, 0.0, NC_FLOAT, redPicture, "longitude latitude");
        createField(observationData, "Lt_SWIR", "Top of Atmosphere SWIR Band Reflectance", NULL, "W m-2 sr-1",
                    "", BAD_FLT, "", "", 0.0f, 300.0f, 1.0, 0.0, NC_FLOAT, swirPicture, "longitude latitude");

    } else {
        createField(observationData, "rhot_blue", "Top of Atmosphere Blue Band Reflectance", NULL, NULL,
                    "rhot = Lt * Pi * earth_sun_distance_correction/(solar_irradiance * cos(solar_zenith))",
                    BAD_FLT, "", "", 0.0f, 1.3f, 1.0, 0.0, NC_FLOAT, bluePicture, "longitude latitude");
        createField(observationData, "rhot_red", "Top of Atmosphere Red Band Reflectance", NULL, NULL,
                    "rhot = Lt * Pi * earth_sun_distance_correction/(solar_irradiance * cos(solar_zenith))",
                    fillDouble, "", "", 0.0f, 1.3f, 1.0, 0.0, NC_FLOAT, redPicture, "longitude latitude");
        createField(observationData, "rhot_SWIR", "Top of Atmosphere SWIR Band Reflectance", NULL, NULL,
                    "rhot = Lt * Pi * earth_sun_distance_correction/(solar_irradiance * cos(solar_zenith))",
                    BAD_FLT, "", "", 0.0f, 1.3f, 1.0, 0.0, NC_FLOAT, swirPicture, "longitude latitude");
    }
    createField(observationData, "qual_blue", "Blue Band Quality Flag", NULL, NULL, NULL, BAD_UBYTE, "1UB",
                "saturation", 0, 0, 1.0, 0.0, NC_UBYTE, bluePicture, "longitude latitude");
    createField(observationData, "qual_red", "Red Band Quality Flag", NULL, NULL, NULL, BAD_UBYTE, "1UB, 2UB",
                "saturation HAM-B_striping", 0, 0, 1.0, 0.0, NC_UBYTE, redPicture, "longitude latitude");
    createField(observationData, "qual_SWIR", "SWIR Band Quality Flag", NULL, NULL, NULL, BAD_UBYTE, "1UB",
                "saturation", 0, 0, 1.0, 0.0, NC_UBYTE, swirPicture, "longitude latitude");

    return EXIT_SUCCESS;
}

int Level1bFile::parseDims(string dimString, vector<NcDim> &varDims) {
    istringstream dimStream(dimString);
    string varDimName;

    while (getline(dimStream, varDimName, ',')) {
        varDimName.erase(0, varDimName.find_first_not_of(" \t"));
        varDimName.erase(varDimName.find_last_not_of(" \t)") + 1);

        if (varDimName.empty())
            continue;

        bool dimensionFound = false;
        for (const auto &dim : ncDims) {
            try {
                if (varDimName == dim.getName()) {
                    varDims.push_back(dim);
                    dimensionFound = true;
                    break;
                }
            } catch (const NcException &e) {
                cerr << "Error accessing dimension: " << varDimName << endl;
                cerr << "Exception: " << e.what() << endl;
                return EXIT_FAILURE;
            }
        }

        if (!dimensionFound) {
            cerr << "Dimension not found: " << varDimName << endl;
            return EXIT_FAILURE;
        }
    }

    return EXIT_SUCCESS;
}

int Level1bFile::createField(NcGroup &parentGroup, const char *shortName, const char *longName,
                             const char *stdName, const char *units, const char *description,
                             double fillValue, const char *flagMasks, const char *flagMeanings,
                             double validMin, double validMax, double scale, double offset, int ncType,
                             vector<NcDim> &dimVec, string coordinates) {
    /* Create the NCDF dataset */
    NcVar var;

    try {
        var = parentGroup.addVar(shortName, ncType, dimVec);
    } catch (NcException &e) {
        cout << e.what() << endl;
        cerr << "Failure creating variable: " << shortName << endl;
        exit(1);
    }
    int varId = var.getId();
    int parentId = var.getParentGroup().getId();
    std::string parentGroupName = var.getParentGroup().getName();
    int32_t dimensionIds[5];
    for (size_t idim = 0; idim < var.getDims().size(); idim++) {
        dimensionIds[idim] = var.getDims()[idim].getId();
    }

    if (validMin != fillValue) {
        if (ncType == NC_BYTE) {
            int8_t i8 = fillValue;
            var.setFill(true, (void *)&i8);
        } else if (ncType == NC_UBYTE) {
            uint8_t ui8 = fillValue;
            var.setFill(true, (void *)&ui8);
        } else if (ncType == NC_SHORT) {
            int16_t i16 = fillValue;
            var.setFill(true, (void *)&i16);
        } else if (ncType == NC_USHORT) {
            uint16_t ui16 = fillValue;
            var.setFill(true, (void *)&ui16);
        } else if (ncType == NC_INT) {
            int32_t i32 = fillValue;
            var.setFill(true, (void *)&i32);
        } else if (ncType == NC_UINT) {
            uint32_t ui32 = fillValue;
            var.setFill(true, (void *)&ui32);
        } else if (ncType == NC_FLOAT) {
            float f32 = fillValue;
            var.setFill(true, (void *)&f32);
        } else {
            var.setFill(true, (void *)&fillValue);
        }
    }

    /* vary chunck size based on dimensions */
    vector<size_t> chunkVec;
    if (parentGroupName == "geolocation_data" || parentGroupName == "observation_data") {
        if (dimVec.size() == 2)
            chunkVec = {512, 1272};  // 256 lines, all pixels(1272)
        if (dimVec.size() == 3)
            chunkVec = {40, 16, 1272};  // 40 bands, 16 lines, all pixels(1272)
        nc_init_compress(parentId, varId, dimensionIds, dimVec.size(), chunkVec.data(), DEFLATE_LEVEL);
    }

    /* Add a "long_name" attribute */
    try {
        var.putAtt("long_name", longName);
    } catch (NcException &e) {
        e.what();
        cerr << "Failure creating 'long_name' attribute: " << longName << endl;
        exit(1);
    }

    vector<int8_t> flagValues = parseFlagValues(flagMasks);
    if (flagValues.size() > 0) {
        try {
            var.putAtt("flag_masks", ncType, flagValues.size(), flagValues.data());
        } catch (NcException &e) {
            e.what();
            cerr << "Failure creating 'flag_masks' attribute: " << longName << endl;
            exit(1);
        }
    }

    /* Add a "flag_meanings" attribute if specified*/
    if (strcmp(flagMeanings, "") != 0) {
        try {
            var.putAtt("flag_meanings", flagMeanings);
        } catch (NcException &e) {
            e.what();
            cerr << "Failure creating 'flag_meanings' attribute: " << flagMeanings << endl;
            exit(1);
        }
    }

    /* Add "valid_min/max" attributes if specified */
    if (validMin < validMax) {
        try {
            var.putAtt("valid_min", ncType, 1, &validMin);
            var.putAtt("valid_max", ncType, 1, &validMax);
            // assignValidBounds(var, validMin, validMax, ncType, shortName);
        } catch (NcException const &e) {
            cerr << e.what() << "\nError assigning validMin or validMax" << endl;
        }
    }

    /* Add "scale_factor" and "add_offset" attributes if specified */
    if (scale != 1.0 || offset != 0.0) {
        try {
            var.putAtt("scale_factor", NC_DOUBLE, 1, &scale);
        } catch (NcException &e) {
            e.what();
            cerr << "Failure creating 'scale_factor' attribute: " << scale << endl;
            exit(1);
        }

        try {
            var.putAtt("add_offset", NC_DOUBLE, 1, &offset);
        } catch (NcException &e) {
            e.what();
            cerr << "Failure creating 'add_offset' attribute: " << offset << endl;
            exit(1);
        }
    }

    /* Add a "units" attribute if one is specified */
    if (units != NULL && *units != 0) {
        try {
            var.putAtt("units", units);
        } catch (NcException &e) {
            e.what();
            cerr << "Failure creating 'units' attribute: " << units << endl;
            exit(1);
        }
    }

    /* Add a "description" attribute if one is specified */
    if (description != NULL && *description != 0) {
        try {
            var.putAtt("description", description);
        } catch (NcException &e) {
            e.what();
            cerr << "Failure creating 'description' attribute: " << description << endl;
            exit(1);
        }
    }

    /* Add a "standard_name" attribute if one is specified */
    if (stdName != NULL && *stdName != 0) {
        try {
            var.putAtt("standard_name", stdName);
        } catch (NcException &e) {
            e.what();
            cerr << "Failure creating 'standard_name' attribute: " << stdName << endl;
            exit(1);
        }
    }

    /* Add a "coordinates" attribute if specified*/
    if (!coordinates.empty()) {
        try {
            var.putAtt("coordinates", coordinates);
        } catch (NcException &e) {
            e.what();
            cerr << "Failure creating 'coordinates' attribute: " << coordinates << endl;
            exit(1);
        }
    }

    return 0;
}

vector<int8_t> Level1bFile::parseFlagValues(string flagMasks) {
    if (flagMasks.empty()) {
        return {};
    }

    using namespace boost;
    using namespace boost::algorithm;
    typedef tokenizer<char_separator<char>> Tokenizer;
    char_separator<char> space(" ");

    vector<int8_t> integerFlagValues;

    Tokenizer tokenMaster(flagMasks, space);  // Hello I am the token master are you the gatekeeeeeper?

    for (auto &flagValueUb : tokenMaster) {
        string flagValue = trim_left_copy(string(flagValueUb.begin(), flagValueUb.begin() + 1));
        int value = lexical_cast<int>(flagValue);
        integerFlagValues.push_back(static_cast<int8_t>(value));
    }

    return integerFlagValues;
}

int Level1bFile::writeGranuleMetadata(std::string tstart, std::string tend, std::string l1b_name) {
    try {
        l1bFile->putAtt("time_coverage_start", tstart.c_str());
        l1bFile->putAtt("time_coverage_end", tend.c_str());

        // Write product file name
        l1bFile->putAtt("product_name", l1b_name.c_str());
    } catch (NcException const &e) {
        cerr << e.what() << "\n Failure writing metadata for " << fileName << " (AKA " << l1b_name << ")"
             << endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

int Level1bFile::writeBandInfo(Device deviceType, const vector<float> &l1aWavelengths,
                               const vector<double> &irradiances, const size_t numBands,
                               const size_t numInsBands, float **gainAggMat, float **insAggMat,
                               float ***m12Coefs, float ***m13Coefs) {
    vector<float> l1bWavelengths = vector<float>(numBands);
    vector<float> sum = vector<float>(numInsBands);
    string color = determineColor(deviceType);

    for (size_t i = 0; i < numInsBands; i++) {
        sum[i] = 0.0;
        for (size_t j = 0; j < NUM_BLUE_WAVELENGTHS; j++) {
            sum[i] += l1aWavelengths[j] * gainAggMat[j][i];
        }
    }

    // bamat#sum
    for (size_t i = 0; i < numBands; i++) {
        l1bWavelengths[i] = 0.0;
        for (size_t j = 0; j < numInsBands; j++) {
            l1bWavelengths[i] += insAggMat[j][i] * sum[j];
        }
    }

    const size_t HAM_SIDES = 2;
    const size_t POLARIZATION_COEFS = 3;
    float ***l1bM12 = allocate3d_float(numBands, HAM_SIDES, POLARIZATION_COEFS);
    float ***l1bM13 = allocate3d_float(numBands, HAM_SIDES, POLARIZATION_COEFS);

    for (size_t l = 0; l < POLARIZATION_COEFS; l++) {
        for (size_t m = 0; m < HAM_SIDES; m++) {
            for (size_t k = 0; k < 2; k++) {  // k=0 for M12, k=1 for M13
                // Calculate M12 or m13
                for (size_t i = 0; i < numBands; i++) {
                    float &result = (k == 0 ? l1bM12[i][m][l] : l1bM13[i][m][l]);
                    result = 0.0;
                    for (size_t j = 0; j < numInsBands; j++) {
                        // Calculate sum
                        float sum = 0.0;
                        for (size_t n = 0; n < NUM_CCD_WAVELENGTHS; n++) {
                            sum += (k == 0 ? m12Coefs[n][m][l] : m13Coefs[n][m][l]) * gainAggMat[n][j];
                        }
                        result += insAggMat[j][i] * sum;
                    }
                }
            }
        }
    }

    sensorBandParameters.getVar(color + "_solar_irradiance").putVar({0}, {numBands}, irradiances.data());

    vector<size_t> start(3, 0);
    vector<size_t> count(3, 0);
    count[0] = numBands;
    count[1] = 2;
    count[2] = 3;

    if (deviceType == SWIR) {
        sensorBandParameters.getVar(color + "_wavelength").putVar({0}, {numBands}, l1aWavelengths.data());
        sensorBandParameters.getVar("SWIR_bandpass").putVar({0}, {NUM_SWIR_WAVELENGTHS}, SWIR_BANDPASS);
        sensorBandParameters.getVar(color + "_m12_coef").putVar(start, count, &m12Coefs[0][0][0]);
        sensorBandParameters.getVar(color + "_m13_coef").putVar(start, count, &m13Coefs[0][0][0]);
    } else {
        sensorBandParameters.getVar(color + "_wavelength").putVar({0}, {numBands}, l1bWavelengths.data());
        sensorBandParameters.getVar(color + "_m12_coef").putVar(start, count, &l1bM12[0][0][0]);
        sensorBandParameters.getVar(color + "_m13_coef").putVar(start, count, &l1bM13[0][0][0]);
    }

    free3d_float(l1bM12);
    free3d_float(l1bM13);

    return EXIT_SUCCESS;
}

Level1bFile::Level1bFile() {
    cout << "Level1bFile::Level1bFile() - Empty constructor will not work" << endl;
    exit(EXIT_FAILURE);
}

Level1bFile::Level1bFile(string filename) {
    l1bFile = new NcFile(filename, NcFile::replace);
    this->fileName = filename;
}

Level1bFile::~Level1bFile() {
    try {
        l1bFile->close();
    } catch (NcException &e) {
        cout << e.what() << endl;
        cerr << "Failure closing: " + fileName << endl;
        exit(EXIT_FAILURE);
    }
    delete l1bFile;
}
