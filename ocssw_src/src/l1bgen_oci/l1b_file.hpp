/**
 * @name l1b_file.hpp
 *
 * @brief Functions and variables for Level 1B (L1B or l1b) file.
 *
 * @authors Joel Gales (SAIC), Jakob C Lindo (SSAI)
 * @date Aug 2024
 */

#ifndef __L1B_FILE_H__
#define __L1B_FILE_H__

#include <netcdf>
#include "vecxD.hpp"
#include "device.hpp"

#define NUM_BLUE_WAVELENGTHS 512
#define NUM_RED_WAVELENGTHS 512
#define NUM_CCD_WAVELENGTHS 512
#define NUM_SWIR_WAVELENGTHS 9

const float SWIR_BANDPASS[NUM_SWIR_WAVELENGTHS] = {45, 80, 30, 30, 15, 75, 75, 50, 75};

class Level1bFile {
   public:
    int numGroups;
    int numDimensions;
    netCDF::NcDim ncDims[1000];
    std::string fileName;
    Level1bFile();
    Level1bFile(std::string name);
    ~Level1bFile();

    netCDF::NcFile *l1bFile;
    netCDF::NcGroup sensorBandParameters;
    netCDF::NcGroup scanLineAttributes;
    netCDF::NcGroup geolocationData;
    netCDF::NcGroup navigationData;
    netCDF::NcGroup observationData;
    netCDF::NcGroup spatialSpectralModes;

    /**
     * @brief Create a new NetCDF file for geolocation data
     *
     * @param l1aFilename The name of the L1A file
     * @param numScans The number of scans in the L1A file
     * @param numBlueBands The number of blue bands
     * @param numRedBands The number of red bands
     * @param numHyperSciPix The number of hyperspectral science pixels
     * @param numSwirPixels The number of SWIR pixels
     * @param numSwirBands The number of SWIR bands
     * @param pprOffset Pulse per revolution (PPR) offset of the rotating telescope assembly (RTA) in encoder
     * counts
     * @param radianceGenerationEnabled Boolean flag indicating if radiance generation is enabled
     *
     * @return int Status code (0 for success, non-zero for error)
     */
    int createFile(const char *l1aFilename, size_t numScans, size_t numBlueBands, size_t numRedBands,
                   size_t numHyperSciPix, size_t numSwirPixels, size_t numSwirBands, int32_t *pprOffset,
                   bool radianceGenerationEnabled);

    /**
     * @brief Parse dimension string and populate a vector with corresponding NetCDF dimensions
     *
     * @param dimString A string containing dimension names separated by spaces
     * @param varDims A vector to be populated with the corresponding NetCDF dimensions
     * @return int Status code (0 for success, non-zero for error)
     */
    int parseDims(std::string dimString, std::vector<netCDF::NcDim> &varDims);

    /**
     * @brief Create a field in the NetCDF file
     *
     * @param parentGroup The parent group to create the field in
     * @param shortName The short name of the field
     * @param longName The long name of the field
     * @param stdName The standard name of the field
     * @param units The units of the field
     * @param description The description of the field
     * @param fillValue The fill value for the field
     * @param flagMasks The flag masks for the field
     * @param flagMeanings The flag meanings for the field
     * @param validMin The minimum valid value for the field
     * @param validMax The maximum valid value for the field
     * @param scale The scale factor for the field
     * @param offset The offset for the field
     * @param ncType The NetCDF data type for the field
     * @param varVec The vector of dimensions for the field
     * @param coordinates The coordinate variables for the field. E.g., "latitude longitude"
     *
     * @return int Status code (0 for success, non-zero for error)
     */
    int createField(netCDF::NcGroup &parentGroup, const char *shortName, const char *longName,
                    const char *stdName, const char *units, const char *description, double fillValue,
                    const char *flagMasks, const char *flagMeanings, double validMin, double validMax,
                    double scale, double offset, int ncType, std::vector<netCDF::NcDim> &varVec,
                    std::string coordinates);

    /**
     * @brief Write granule metadata to the NetCDF file
     *
     * @param startTime The start time of the granule
     * @param endTime The end time of the granule
     * @param productName The name of the L1B product in metadata
     *
     * @return int Status code (0 for success, non-zero for error)
     */

    int writeGranuleMetadata(std::string startTime, std::string endTime, std::string productName);

    /**
     * @brief Parse flag values from a string of flag masks, whose format is expected to be <number><type
     * specifier like 'UB', etc>
     *
     * @param flagMasks A string containing flag mask values
     * @return std::vector<int8_t> A vector of parsed flag values
     */
    std::vector<int8_t> parseFlagValues(std::string flagMasks);

    /**
     * @brief Write spectral band information.
     * Calculate band centers for aggregated hyperspectral bands and write them out to an NcFile. Takes care
     * of red and blue, not SWIR
     */
    /**
     * @brief Writes band information to the L1B file
     *
     * @param l1bFile Reference to the Level 1B output file
     * @param deviceType Enum indicating the device type (BLUE or RED)
     * @param l1aWavelengths Vector of L1A wavelengths
     * @param irradiances Vector of solar irradiance values
     * @param numL1aBands Number of L1A bands
     * @param numInsBands Number of instrument bands
     * @param gainAggMat Gain aggregation matrix
     * @param insAggMat Instrument aggregation matrix
     * @param calLut Calibration lookup table
     *
     * @return int Status code (EXIT_SUCCESS or EXIT_FAILURE)
     *
     * This function calculates and writes band-specific information to the L1B file, including:
     * - Wavelengths
     * - Solar irradiances
     * - Mueller matrix coefficients (m12 and m13)
     */

    int writeBandInfo(Device deviceType, const std::vector<float> &l1aWavelengths,
                      const std::vector<double> &irradiances, const size_t numL1aBands,
                      const size_t numInsBands, float **gainAggMat, float **insAggMat, float ***m12Coefs,
                      float ***m13Coefs);

    int close();
};

#endif