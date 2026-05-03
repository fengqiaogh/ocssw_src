#include "Wave3DParsing.hpp"
#include "l3mapgen.h"
#include <L3FileSMI.h>

using namespace std;

/** Namespace to hold 3d expansion/mapping between 3D and 2D */
namespace wv3d {
/** List of wavelengths separated */
vector<string> waveLength3DListSeparated;
/** Mapping from 2D names to 3D expansion */
unordered_map<string, vector<float>> wv3dName2dTo3dExpansion;
/** Mapping from 3D names to 2D names */
unordered_map<string, string> wv3d3dNameTo2D;
/** List of output products with 3D expansion */
vector<string> outputProductsWith3d;
/** Size of wavelength 3D */
size_t wavelength3dSize;
/** map of 2d products and its floating point wavelength */
unordered_map<string, float> wavelength2dMap;

bool numericalOrder(const std::string& a, const std::string& b) {
    return std::stoi(a) < std::stoi(b);
}
}  // namespace wv3d

/**
 * @brief Gets the 3D expansion mapping for 2D names.
 * @return A reference to the unordered map of 3D expansion mapping.
 */
const unordered_map<string, vector<float>>& getWv3dName2dTo3dExpansion() {
    return wv3d::wv3dName2dTo3dExpansion;
}

/**
 * @brief Gets the mapping from 3D names to 2D names.
 * @return A reference to the unordered map of 3D to 2D name mapping.
 */
const unordered_map<string, string>& getWv3d3dNameTo2D() {
    return wv3d::wv3d3dNameTo2D;
}

/**
 * @brief Gets the length of the 3D wavelength.
 * @return The size of the 3D wavelength.
 */
size_t getLenWv3d() {
    return wv3d::wavelength3dSize;
}

/**
 * @brief get the map of 2D product's floating point wavelengths
 * @return A reference to unordered_map<2dProductName, wavelength>
 */
const unordered_map<string, float>& getWavelength2dMap() {
    return wv3d::wavelength2dMap;
}

void getProductNames(const string &productName ,vector<string>& productNameList, l3::L3File* l3File, clo_optionList_t* optionList) {
    // make temp unordered set for quick lookup;
    unordered_set<string> table_product_available_list;
    const int number_of_products = l3File->getNumProducts();
    for (int i = 0; i < number_of_products; i++) {
        const string& name = l3File->getProductName(i);
        table_product_available_list.insert(name);
    }

    boost::split(productNameList, productName, boost::is_any_of(","));
    // all the wavelength for the sensor
    int32_t* wave_array_of_the_sensor;
    float* fwave_array_of_the_sensor;
    // quick look up for wv
    unordered_set<int32_t> look_up_for_wv;
    int sensorID = (*l3File->getMetaData()).sensorID;
    const int total_num_bands = rdsensorinfo(sensorID, 0, "Lambda", (void**)&wave_array_of_the_sensor);
    if (total_num_bands != rdsensorinfo(sensorID, 0, "fwave", (void**)&fwave_array_of_the_sensor)) {
        fprintf(stderr, "-E-: Couldn't retrive floating point wavelength\n. See %s:%d", __FILE__, __LINE__);
        exit(1);
    }
    // to lookup a rank
    productInfo_t* productInfo;
    productInfo = allocateProductInfo();
    for (int i = 0; i < total_num_bands; i++) {
        look_up_for_wv.insert(wave_array_of_the_sensor[i]);
    }
    if (clo_isSet(optionList, "wavelength_3d")) {
        const string wavelengths_3d_list = clo_getRawString(optionList, "wavelength_3d");
        boost::split(wv3d::waveLength3DListSeparated, wavelengths_3d_list, boost::is_any_of(", "),
                     boost::algorithm::token_compress_on);

        vector<string> colon_expanded_list;
        for (const auto& wv_par : wv3d::waveLength3DListSeparated) {
            if (boost::contains(wv_par, ":")) {
                vector<string> pars;
                boost::split(pars, wv_par, boost::is_any_of(":"));
                if (pars.size() != 2) {
                    EXIT_LOG(cerr << "--Error-: Wrong range specifier: " << wv_par << endl;)
                }
                try {
                    int wav_st = boost::lexical_cast<int32_t>(pars.at(0));
                    int wav_end = boost::lexical_cast<int32_t>(pars.at(1));
                    if (look_up_for_wv.count(wav_st) == 0) {
                        EXIT_LOG(cerr << "--Error--: The start wavelength " << wav_st
                                      << " is not found in sensor wv list.\n Check "
                                         "the range "
                                      << wv_par << endl);
                    }
                    if (look_up_for_wv.count(wav_end) == 0) {
                        EXIT_LOG(cerr << "--Error--: The end wavelength " << wav_end
                                      << " is not found in sensor wv list.\n Check "
                                         "the range "
                                      << wv_par << endl);
                    }
                    for (int32_t i = wav_st; i <= wav_end; i++) {
                        if (look_up_for_wv.count(i) == 0)
                            continue;
                        colon_expanded_list.push_back(boost::lexical_cast<string>(i));
                    }
                } catch (const boost::bad_lexical_cast& e) {
                    EXIT_LOG(cerr << e.what() << '\n'; cerr << "--Error--: Provided wavelength are not valid "
                                                               "numbers. "
                                                            << wv_par << endl;)
                }
            } else {
                colon_expanded_list.push_back(wv_par);
            }
        }
        wv3d::waveLength3DListSeparated = std::move(colon_expanded_list);
    }

    vector<string> temp_prod_name_list;
    unordered_set<string> already_set_wv;
    // check that there are no duplicates in the wv list
    for (const auto& wv : wv3d::waveLength3DListSeparated)
        if (already_set_wv.count(wv) == 0) {
            already_set_wv.insert(wv);
        } else {
            EXIT_LOG(cerr << "--Error--: A duplicate found in the wavelength_3d list  " << wv << endl);
        }

    for (size_t i = 0; i < productNameList.size(); i++) {
        vector<string> names;
        boost::split(names, productNameList.at(i), boost::is_any_of(":"));
        string cleanName = names.at(0);
        int result = findProductInfo(cleanName.c_str(), sensorID, productInfo);
        if (result != 1) {
            EXIT_LOG(cerr << "--Error--: Could not find the product: " << cleanName << endl);
        }
        int prod_rank = productInfo->rank;
        const string suffix = productInfo->suffix;
        string prefix = productInfo->prefix;

        // fill up wavelength map for 2D products
        if(prod_rank == 2 && table_product_available_list.count(cleanName) != 0 ) {
            float wavelength_floating;
            if (l3File->getProductAttribute(cleanName, "wavelength", &wavelength_floating,sizeof(float))) {
                wv3d::wavelength2dMap.insert({cleanName,wavelength_floating});
            }
            temp_prod_name_list.push_back(productNameList.at(i));

        } else if(prod_rank == 3) {
           {
                if (wv3d::waveLength3DListSeparated.empty()) {
                    for (const auto& product_in_l3in : table_product_available_list) {
                        int res = findProductInfo(product_in_l3in.c_str(), sensorID, productInfo);

                        if (res != 1) {
                            EXIT_LOG(cerr << "--Error--: Could not read the product info " << product_in_l3in
                                    << endl);
                        }
                        const string local_suffix = productInfo->suffix;
                        const string localName = productInfo->productName;
                        if (boost::contains(cleanName, localName)) {
                            if (!local_suffix.empty()) {
                                if (!boost::contains(cleanName, local_suffix))
                                    continue;
                            } else {
                                if (cleanName != localName)
                                    continue;
                            }

                            const string wave_length = boost::lexical_cast<string>(productInfo->prod_ix);
                            if (wave_length.empty()) {
                                EXIT_LOG(cerr << "--Error--: Not valid 2D slice of " << cleanName
                                        << "in the l3bin file " << endl);
                            }
                            wv3d::waveLength3DListSeparated.push_back(wave_length);
                        }
                    }
                    sort(wv3d::waveLength3DListSeparated.begin(), wv3d::waveLength3DListSeparated.end(),
                            wv3d::numericalOrder);
                }
                bool prod_3d_expand_found = false;
                for (size_t i = 0; i < wv3d::waveLength3DListSeparated.size(); i++) {
                    string wv = wv3d::waveLength3DListSeparated.at(i);
                    string prod_3d_name = prefix + "_" + wv + suffix;
                    if (table_product_available_list.count(prod_3d_name) == 0) {
                        EXIT_LOG(cerr << "--Error--: Neither product " << cleanName
                                << " or its wavelength 3d slice " << prod_3d_name
                                << " are found. \nExiting ... " << endl);
                    } else {
                        prod_3d_expand_found = true;
                        names.at(0) = prod_3d_name;
                        int32_t wavelength;
                        try {
                            wavelength = boost::lexical_cast<int32_t>(wv);
                        } catch (const boost::bad_lexical_cast& e) {
                            EXIT_LOG(cerr << e.what() << '\n'; cerr
                                    << "--Error--: Provided wavelength are not valid "
                                    "numbers. \nExiting...");
                        }
                        float wavelength_floating = static_cast<float>(wavelength);
                        if (l3File->getProductAttribute(prod_3d_name, "wavelength", &wavelength_floating,sizeof(float))) {
                            if (wavelength != std::lroundf(wavelength_floating)) {
                                fprintf(stderr, "-E-: The wavelength %d from the product name %s. \n", wavelength,
                                        prod_3d_name.c_str());
                                fprintf(stderr, "The wavelength %f from the product attribute. See %s:%d \n",
                                        wavelength_floating, __FILE__, __LINE__);
                                exit(1);
                            }
                        } else {
                            auto it = std::find(wave_array_of_the_sensor,
                                    wave_array_of_the_sensor + total_num_bands, wavelength);
                                    // check if it is valid and wavelength found
                            if (it == wave_array_of_the_sensor + total_num_bands) {
                                EXIT_LOG(cerr << "--Error--: Wavelength " << wavelength
                                        << " not found for the sensor " << sensorID << endl);
                            }
                            int index = std::distance(wave_array_of_the_sensor, it);
                            wavelength_floating = fwave_array_of_the_sensor[index];
                        }
                        wv3d::wv3dName2dTo3dExpansion[cleanName].push_back(wavelength_floating);
                        wv3d::wv3d3dNameTo2D[prod_3d_name] = cleanName;
                        string temp_prod_3d_name;
                        for (const auto& name : names) {
                            if (!temp_prod_3d_name.empty())
                                temp_prod_3d_name += ":";
                            temp_prod_3d_name += name;
                        }
                        temp_prod_name_list.push_back(temp_prod_3d_name);
                    }
                }
                if (!prod_3d_expand_found) {
                    EXIT_LOG(cerr << "--Error--: Product not found : " << cleanName << endl);
                }

            }
            wv3d::outputProductsWith3d.push_back(cleanName);
        }
    }
    freeProductInfo(productInfo);

    productNameList = std::move(temp_prod_name_list);
    wv3d::wavelength3dSize = wv3d::waveLength3DListSeparated.size();
    if (wv3d::wavelength3dSize >= 1) {
        if (wv3d::wv3dName2dTo3dExpansion.size() > 0) {
            if (clo_isSet(optionList, "ofile2")) {
                const string oformatStr2 = getFileFormatName(clo_getString(optionList, "oformat2"));
                if (oformatStr2.compare("netCDF4") != 0) {
                    EXIT_LOG(cerr << "The user supplied a 3D product "
                                  << " and the output format for ofile2 is not netCDF4.\n"
                                  << "Exiting ... " << endl);
                }
            }
        }
    }
}
