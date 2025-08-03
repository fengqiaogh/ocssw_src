#ifndef _L1AFILE_OCI_HPP_
#define _L1AFILE_OCI_HPP_

#include <netcdf>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <boost/algorithm/string.hpp>

#include "global_attrs.h"
#include "l0PacketUtil_oci.h"
#include "dimension_shape.hpp"



/**
   Header for the L1aFile class that takes care of generating an L1A file and also writing to it.
 */

// CONSTANTS -- to set cache size, nelems and preemption policy for Nc file
const size_t CHUNK_CACHE_SIZE = 96 * 1024 * 1024;
const size_t CHUNK_CACHE_NELEMS = 251;
const float CHUNK_CACHE_PREEMPTION = 1.0;

// Compression settings for Nc variable
const size_t CHUNK_BANDS = 36;
const size_t CHUNK_PIXELS = 256;
const size_t CHUNK_LINES = 128;

class L1aFile {
   public:
    L1aFile();
    ~L1aFile();


    std::string getFileName();

    // free's the l1afile pointer
    void freeFile();

    /**
     * @brief Write global attributes to the current NetCDF file.
     * @param history command line parameters that was used to run the program
     * @param doi
     * @param pversion
     * @return
     */
    int writeGlobalAttributes(std::string history, std::string doi, std::string pversion);

    /**
     * @brief create the l1aOutputFile and initialize it so it can be written into.
     * @param l1aFileName
     * @param maxScans
     * @param numCcdPixels
     * @param numBlueBands
     * @param numRedBands
     * @param numSwirPixels
     * @param numDarkCcdPixels
     * @return
     */
    int initializeL1aFile(char *l1aFileName, uint16_t maxScans, uint16_t numCcdPixels, uint16_t numBlueBands,
                          uint16_t numRedBands, uint16_t numSwirPixels, uint16_t numDarkCcdPixels);

    /**
     * @brief write 1 line of science data for the "science_data" group with variables:
     *          sci_blue, sci_red, sci_SWIR, frm_type_SWIR
     * @param scanNum
     * @param numBlueBands
     * @param numRedBands
     * @param numSwirBands
     * @param numCcdPixels
     * @param numSwirPixels
     * @param blueScienceData
     * @param redScienceData
     * @param swirScienceData
     * @param swirFrameTypesSciData
     * @return
     */
    int writeScienceData(DimensionShape* dimShape, uint32_t scanNum, uint16_t numBlueBands, uint16_t numRedBands, uint16_t numSwirBands,
                         uint16_t numCcdPixels, uint16_t numSwirPixels, uint16_t **blueScienceData,
                         uint16_t **redScienceData, uint32_t **swirScienceData,
                         int8_t *swirFrameTypesSciData);

    /**
     * @brief write  to "processing_control" group and "input_parameter" group that specifies all L0
     *          and HKT files used to generate the L1A file.
     * @param hktList
     * @param l0List
     * @param time_start
     * @param maxgap
     * @param nametag
     * @param swir_loff_set
     * @param outlist
     * @param outfile
     * @param doi
     * @param pversion
     * @param isSPW
     * @param VERSION
     * @return
     */
    int writeProcessingControl(std::string hktList, std::string l0List, std::string time_start,
                               std::string maxgap, std::string nametag, std::string swir_loff_set,
                               std::string outlist, std::string outfile, std::string doi,
                               std::string pversion, std::string isSPW, std::string VERSION);

    /**
     * @brief   Writes dark calibration data for red, blue and swir pixels to a NetCDF file.
     *          Also writes swir's dark calibration frame type data.
     * @param dimShape - current size of dimensions for the current l1a file
     * @param isc scan number
     * @param numBlueBands
     * @param numRedBands
     * @param numSwirBands
     * @param numDarkCcdPixels
     * @param numDarkSwirPixels
     * @param blueDarkCalibrationData
     * @param redDarkCalibrationData
     * @param swirDarkCalibrationData
     * @param swirDarkCalFrameTypeData
     * @return 0 on success
     */
    int writeCalibrationData(DimensionShape *dimShape, uint32_t isc, uint16_t numBlueBands, uint16_t numRedBands, uint16_t numSwirBands,
                             uint16_t numDarkCcdPixels, uint16_t numDarkSwirPixels,
                             uint16_t *blueDarkCalibrationData, uint16_t *redDarkCalibrationData,
                             uint32_t *swirDarkCalibrationData, int8_t *swirDarkCalFrameTypeData);

    /**
     * @brief Write scan line attributes to the netcdf file. Grabbing ccsds scan times and start times.
     * @param dimShape - current size of dimensions for the current l1a file
     * @param scanNum
     * @param ancillaryData
     * @param sciPacketSequenceError
     * @param ccdLineError
     * @param spinID
     * @param starttime struct that contains year, day and seconds
     * @return
     */
    int writeScanMetaData(DimensionShape* dimShape, uint32_t scanNum, uint8_t *ancillaryData, uint8_t *sciPacketSequenceError,
                          int8_t *ccdLineError, int32_t *spinID, AncillaryPktTimeStamp &starttime);

    /**
     * @brief Write the files start, end times, file name and basic information about the l1a file
     * @param starttime
     * @param endtime
     * @param l1aFileName
     * @param startDirectionStr
     * @param endDirectionStr
     * @param dataType
     * @param swirModeIndex
     * @param cdsModeIndex
     * @param outlistFile
     * @return
     */
    int writeGlobalMetaData(AncillaryPktTimeStamp &starttime, AncillaryPktTimeStamp &endtime,
                            std::string l1aFileName, std::string startDirectionStr,
                            std::string endDirectionStr, short dataType, uint16_t swirModeIndex,
                            uint16_t cdsModeIndex, std::ofstream &outlistFile);

    /**
     * @brief Extract spatial and spectral data and write them to the output NetCDF file
     * @param dimShape - current size of dimensions for the current l1a file
     * @param scanNum
     * @param ancillaryData
     * @return
     */
    int writeAncillaryData(DimensionShape* dimShape, uint32_t scanNum, uint8_t *ancillaryData);

    /**
     * @brief For each telemetry packet, extract the data and group them based on apid. Then write them
     *          to the NetCDF file
     * @param dimShape - current size of dimensions for the current l1a file
     * @param spatialAggList
     * @param numTelemetryPackets
     * @param telemetryData
     * @param spinID
     * @param cdsMode
     * @param scanNum
     * @param starttime
     * @return
     */
    int writeTelemetryData(DimensionShape* dimShape, spatialAggTable *spatialAggList, uint32_t numTelemetryPackets,
                           uint8_t (*telemetryData)[TLMSIZE], int32_t *spinID, uint16_t &cdsMode,
                           uint32_t scanNum, const AncillaryPktTimeStamp &starttime);

    /**
     * @brief Write spacecrafts navigation data like orbit position, speed, etc.
     * @param dimShape - current size of dimensions for the current l1a file
     * @param hktList
     * @param startTime
     * @param endTime
     * @return
     */
    int writeNavigationData(DimensionShape* dimShape, std::string hktList, AncillaryPktTimeStamp &startTime,
                            AncillaryPktTimeStamp &endTime);

    /**
     * @brief Close the l1afile referenced by this class
     * @return 0 on success
     */
    int close();

   private:
    netCDF::NcFile *l1afile = nullptr;
    std::string fileName;

    /**
     * @brief Apply chunking to number of lines, bands and pixels for science data. Compression default to
     *          level 5.
     * @param sciVar the variable to apply the chunking and compression
     * @param currDimSizes current dim list for the variable. Will check to see if it is smaller than the
     *                      default set
     */
    void applyChunkingAndCompressionToSciData(netCDF::NcVar &sciVar,
                                              std::vector<netCDF::NcDim> &currDimSizes);

    /**
     * @brief Create a NcVariable
     * @param groupRef group where the variable belongs to
     * @param varName
     * @param varType Nc types are enums from 0-12 that defines the type in the file (ie. NC_INT == 4)
     * @param varLongName longer name of varName if it has one.
     * @param varDims dimensions of the variable
     * @param fillValue given as int, the API will automatically convert it based on varType
     *                   (ie. (int)255 -> 255UB if the varType was ubyte)
     * @param validMin also given as int and API will convert
     * @param validMax
     * @param units
     */
    void createVariable(netCDF::NcGroup &groupRef, std::string varName, int varType, std::string varLongName,
                        std::vector<netCDF::NcDim> &varDims, double fillValue, double validMin,
                        double validMax, std::string units, std::string reference);

    /**
     * @brief overloaded function to make a flag variable that is a byte.
     * @param groupRef group where the variable belongs to
     * @param varName
     * @param varType
     * @param varLongName
     * @param varDims
     * @param hasFillValue - indicate if fillValue needs to be set
     * @param fillValue if hasFillValue, set this number
     * @param flagValues vector of char where each element defines a possible flag meaning
     * @param flagMeaning string that lists definition of what each flag values mean, separated by space
     */
    void createFlagVariable(netCDF::NcGroup &groupRef, std::string varName, int varType,
                            std::string varLongName, std::vector<netCDF::NcDim> &varDims, bool hasFillValue,
                            double fillValue, std::vector<double> &flagValues, std::string flagMeaning);

    /**
     * @brief Sets the variables fill value. Called from the createVariable or createFlagVariable functions
     * @param varRef variable to set the fill value
     * @param varType NcType reference
     * @param fillValue fill value to use and referenced
     */
    void setVariableFillValue(netCDF::NcVar &varRef, int &varType, double fillValue);

    /**
     * @brief find the start and end index of the current navigation data
     * @param startTimeStr
     * @param startTimeInUnixSecs
     * @param endTimeInUnixSecs
     * @param navigationTimeData
     * @param navigationDataSize
     * @param startIndex
     * @param endIndex
     * @return
     */
    int findNavigationIndex(std::string startTimeStr, double startTimeInUnixSecs, double endTimeInUnixSecs,
                            double *navigationTimeData, size_t navigationDataSize, int &startIndex,
                            int &endIndex);

    /**
     * @brief Reverse the byte order of floats
     * @param inFloat num to change big to little or little to big endian
     * @return
     */
    float reverseFloat(const float inFloat);

    /**
     * @brief add the difference between the hkt time and ancillary epoch time to the navigation arr
     * @param navigationTimeArr time for tilt, orbit or attitude
     * @param arrSize number of records for tilt, orbit and attitude. Or just the size of time arr
     * @param hktFileEpochTime
     * @param ancillaryEpochTime
     */
    void synchronizeEpochTime(std::vector<double> &navigationTimeArr, size_t arrSize,
                              double &hktFileEpochTime, double &ancillaryEpochTime);
};

#endif