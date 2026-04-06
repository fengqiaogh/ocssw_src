#ifndef L3MAPGEN_H /* avoid re-inclusion */
#define L3MAPGEN_H

#include <clo.h>
#include <genutils.h>
#include <nc_gridutils.h>
#include "version.h"
#include <unordered_map>
#include <productInfo.h>
#include "OutFile.h"

extern clo_optionList_t* optionList;

/**
 * @brief Adds all of the accepted command line options to the list.
 * @param list The option list to add the options to.
 * @param softwareVersion The software version string.
 * @return int Returns 0 on success.
 */
int l3mapgenInitOptions(clo_optionList_t* list, const char* softwareVersion);

/**
 * @brief Reads the command line options and all of the default parameter files.
 * 
 * Order for loading the options: (1)loads the main program defaults file, 
 * (2) load the command line (including specified par files)
 * (3) re-load the command line disabling file descending so command
 * line arguments will over ride
 * 
 * @param list The option list to read options into.
 * @param argc The number of command line arguments.
 * @param argv The command line arguments.
 * @return int Returns 0 on success, -1 if the OCDATAROOT environment variable is not defined.
 */
int l3mapgenReadOptions(clo_optionList_t* list, int argc, char* argv[]);

/**
 * @brief Retrieves the mapping from 2D names to 3D expansions.
 * @return A constant reference to the unordered map containing 2D to 3D expansions.
 */
const std::unordered_map<std::string, std::vector<float>>& getWv3dName2dTo3dExpansion();

/**
 * @brief Retrieves the mapping from 3D names to 2D names.
 * @return A constant reference to the unordered map containing 3D to 2D name mappings.
 */
const std::unordered_map<std::string, std::string>& getWv3d3dNameTo2D();


/**
 * @brief Gets the length of the WV3D data.
 * @return The length of the WV3D data.
 */
std::size_t getLenWv3d();

void setupWavelengthProduct(const std::string& productName, productInfo_t* productInfo, int sensorId);

#define EXIT_LOG(...)                                                                     \
    {                                                                                     \
        __VA_ARGS__;                                                                      \
        cerr << "Exiting. See " << __FILE__ << ":" << __LINE__ << endl; \
        exit(EXIT_FAILURE);                                                               \
    }

#endif /* L3MAPGEN_H */
