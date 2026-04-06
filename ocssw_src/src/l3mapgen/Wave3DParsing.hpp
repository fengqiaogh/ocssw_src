#ifndef WAVE3DPARSING_HPP
#define WAVE3DPARSING_HPP

#include <unordered_map>
#include <string>
#include <vector>
#include <sensorInfo.h>
#include <productInfo.h>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <unordered_set>
#include <boost/algorithm/string.hpp>
#include <L3File.h>
#include <clo.h>

/**
 * @brief Gets the 3D expansion mapping for 2D names.
 * @return A reference to the unordered map of 3D expansion mapping.
 */
const std::unordered_map<std::string, std::vector<float>>& getWv3dName2dTo3dExpansion();

/**
 * @brief get the map of 2D product's floating point wavelengths
 * @return A reference to unordered_map<2dProductName, wavelength>
 */
const std::unordered_map<std::string, float>& getWavelength2dMap();

/**
 * @brief Gets the mapping from 3D names to 2D names.
 * @return A reference to the unordered map of 3D to 2D name mapping.
 */
const std::unordered_map<std::string, std::string>& getWv3d3dNameTo2D();

/**
 * @brief Gets the length of the 3D wavelength.
 * @return The size of the 3D wavelength.
 */
size_t getLenWv3d();

/**
 * @brief Get the Product Names List
 * 
 * @param productName an input string of all products
 * @param productNameList the output vector of products, expanded over 3D wavelenght, if applicable.
 * @param l3File l3 file handle
 * @param optionList clo option list
 */
void getProductNames(const std::string &productName ,std::vector<std::string>& productNameList, l3::L3File* l3File, clo_optionList_t* optionList);

#endif
