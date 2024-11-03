#include "expand3D.hpp"
#include "sensorInfo.h"
#include <limits>
#include "find_variable.hpp"

namespace {
std::unordered_map<std::string, std::string> l3d_to_l2_prod_map;
std::unordered_map<std::string, float> l3d_wave_attr;
std::unordered_map<std::string, std::pair<int, int>> l3d_indexes;
std::unordered_map<std::string, std::vector<int>> three_dims_prod_wavelength;
std::unordered_map<std::string, std::vector<float>> three_d_original_prod_wavelength;
std::unordered_set<std::string> original_l2_products;
std::unordered_map<std::string, std::pair<float, float>> min_max_values;
std::unordered_set<int> wavelength_all;
std::string platform, instrument;
int sensorID;
bool ini_3d{false};
bool use_l2_flags{false};

template <class T, class F>
std::vector<F> vector_rounder(const std::vector<T> &inp) {
    std::vector<F> result(inp.size());
    std::transform(inp.begin(), inp.end(), result.begin(),
                   [](const auto &c) { return static_cast<F>(std ::round(c)); });
    return result;
}

}  // namespace

void parse_colon_range(const std::string &wv_range, const std::unordered_map<int, int> &wv3d_list,
                       std::vector<int> &wv_requested_vals, std::vector<int> &wv_requested_indexes,
                       std::unordered_set<int> &wv_requested_vals_lookup) {
    std::vector<std::string> result_colon;
    boost::split(result_colon, wv_range, boost::is_any_of(":"), boost::algorithm::token_compress_on);
    if (result_colon.size() != 2) {
        EXIT_LOG(std::cerr << "--Error-- : The range must be supplied in the form nnn:nnn. The input "
                           << wv_range << " does not match the format" << std::endl)
    }
    try {
        int min_wv = boost::lexical_cast<int>(result_colon.at(0));
        int max_wv = boost::lexical_cast<int>(result_colon.at(1));
        if (wv3d_list.count(min_wv) == 0) {
            EXIT_LOG(std::cerr << "--Error-- : Start range value not found : " << min_wv << std::endl)
        }
        if (wv3d_list.count(max_wv) == 0) {
            EXIT_LOG(std::cerr << "--Error-- : End range value not found : " << max_wv << std::endl)
        }
        for (int wv = min_wv; wv <= max_wv; wv++) {
            if (wv3d_list.count(wv) > 0) {
                if (wv_requested_vals_lookup.count(wv) == 0) {
                    wv_requested_vals.push_back(wv);
                    wv_requested_indexes.push_back(wv3d_list.at(wv));
                    wv_requested_vals_lookup.insert(wv);
                }
            }
        }
    } catch (const boost::bad_lexical_cast &error_boost) {
        EXIT_LOG(std::cerr << "--Error-- : Conversion failed for WV3D : " << error_boost.what() << std::endl)
    }
};
/**
 * @brief
 *
 * @param wv_list user supplied wv3d list
 * @param wv3d_list  wavelength and their indexes;
 * @param wv_requested_vals user requested wavelength list
 * @param wv_requested_indexes slices indexes of the 3D variable
 * @param exit_on_not_found hard exit if wavelength not found
 */
void parse_wv_list(const std::string &wv_list, const std::unordered_map<int, int> &wv3d_list,
                   std::vector<int> &wv_requested_vals, std::vector<int> &wv_requested_indexes,
                   bool exit_on_not_found) {
    std::vector<std::string> result_comma;
    std::unordered_set<int> wv_requested_vals_lookup;
    boost::split(result_comma, wv_list, boost::is_any_of(", "), boost::algorithm::token_compress_on);
    try {
        std::for_each(result_comma.begin(), result_comma.end(), [&](const std::string &wv_str) {
            if (boost::contains(wv_str, ":"))
                parse_colon_range(wv_str, wv3d_list, wv_requested_vals, wv_requested_indexes,
                                  wv_requested_vals_lookup);
            else {
                int wv = boost::lexical_cast<int>(wv_str);
                if (wv3d_list.count(wv) > 0) {
                    if (wv_requested_vals_lookup.count(wv) == 0) {
                        wv_requested_vals.push_back(wv);
                        wv_requested_indexes.push_back(wv3d_list.at(wv));
                        wv_requested_vals_lookup.insert(wv);
                    }
                } else if (exit_on_not_found) {
                    EXIT_LOG(std::cerr << "--Error-- : Wavelength " << wv << " nm "
                                       << "not found");
                }
            }
        });
    } catch (const boost::bad_lexical_cast &error_boost) {
        EXIT_LOG(std::cerr << "--Error-- : Conversion failed for WV3D : " << error_boost.what() << std::endl)
    }
    std::cout << "\n";
}

std::vector<std::string> create_min_max_values(const std::vector<std::string> &products_requested_separated,
                                               productInfo_t *info) {
    std::vector<std::string> output_expanded_products;

    for (const auto &prod_with_min_max : products_requested_separated) {
        std::vector<std::string> prod_data;
        float max_val_preset, min_val_preset;
        boost::algorithm::split(prod_data, prod_with_min_max, boost::is_any_of("="),
                                boost::token_compress_on);  // products_requested}
        {
            if (findProductInfo(prod_data.at(0).c_str(), sensorID, info)) {
                min_val_preset = info->validMin;
                max_val_preset = info->validMax;
            } else {
                min_val_preset = 0.0;
                max_val_preset = std::numeric_limits<float>::max();
            }
        }
        if (prod_data.size() > 1) {
            output_expanded_products.push_back(prod_data.at(0));
            std::vector<std::string> min_max_sep;
            boost::algorithm::split(min_max_sep, prod_data.at(1), boost::is_any_of(":"),
                                    boost::token_compress_on);  // products_requested}
            if (min_max_sep.size() == 1) {
                min_max_values[prod_data.at(0)] = {std::stof(min_max_sep.at(0)), max_val_preset};
            } else if (min_max_sep.size() == 2) {
                float max_val =
                    min_max_sep.at(1).length() > 0 ? std::stof(min_max_sep.at(1)) : max_val_preset;
                float min_val =
                    min_max_sep.at(0).length() > 0 ? std::stof(min_max_sep.at(0)) : min_val_preset;
                min_max_values[prod_data.at(0)] = {min_val, max_val};
            } else {
                EXIT_LOG(std::cerr << "--Error-- : Can't parse the product " << prod_with_min_max
                                   << std::endl)
            }
        } else {
            output_expanded_products.push_back(prod_with_min_max);
            min_max_values[prod_with_min_max] = {min_val_preset, max_val_preset};
        }
    }

    return output_expanded_products;
}

std::string get_3d_product_name(const std::string &prod_3d, const std::string &wv,
                                productInfo_t *info = nullptr) {
    if (findProductInfo(prod_3d.c_str(), sensorID, info)) {
        std::string prefix = info->prefix;
        std::string suffix = info->suffix;
        return prefix + "_" + wv + suffix;
    }
    return prod_3d + "_" + wv;
}

std::string get_3d_product_name(const std::string &prod_3d, int wv, productInfo_t *info = nullptr) {
    return get_3d_product_name(prod_3d, std::to_string(wv), info);
}

std::string get_3d_product_name(const std::string &prod_3d, int *wv, productInfo_t *info = nullptr) {
    std::vector<std::string> parsed_product;
    boost::algorithm::split(parsed_product, prod_3d, boost::is_any_of("_"));
    std::string return_string;
    for (size_t i = 0; i < parsed_product.size(); i++) {
        try {
            int wv_3d = std::stoi(parsed_product.at(i));
            if (wavelength_all.find(wv_3d) != wavelength_all.end()) {
                // check if the name does exist
                std::string name_l2 = parsed_product.at(0);
                for (size_t k = 1; k < i; k++) {
                    name_l2 += "_" + parsed_product.at(k);
                }
                for (size_t k = i + 1; k < parsed_product.size(); k++) {
                    name_l2 += "_" + parsed_product.at(k);
                }
                return_string = name_l2;
                if (findProductInfo(return_string.c_str(), sensorID, info)) {
                    std::string prefix = info->prefix;
                    std::string suffix = info->suffix;
                    return_string = prefix + suffix;
                }
                if (original_l2_products.count(name_l2) || original_l2_products.count(return_string)) {
                    *wv = wv_3d;
                    break;
                }
            }
        } catch (...) {
        }
    }
    return return_string;
}

/**
 * @brief
 * The user can supply a 3D product such as  "Lt" or "Rrs" and they will be expanded: Rrs_449 etc
 * @param products_requested_separated - list of all products.
 */
void expand_l2_requested(std::vector<std::string> &products_requested_separated, productInfo_t *info) {
    std::unordered_set<std::string> remove;
    std::vector<std::string> add_3d_expanded;
    for (auto &prod : products_requested_separated) {
        if (three_dims_prod_wavelength.count(prod) > 0) {
            remove.insert(prod);
            std::vector<int> &wavelengths_3d = three_dims_prod_wavelength[prod];
            for (const auto &wv : wavelengths_3d) {
                std::string exp3d_name =
                    get_3d_product_name(prod, wv, info);  // prod + "_" + std::to_string(wv);
                add_3d_expanded.push_back(exp3d_name);
                min_max_values[exp3d_name] = min_max_values.at(prod);
            }
        }
    }

    if (!remove.empty()) {
        for (const auto &rem : remove) {
            // std::cout << "All wavelength will be binned and sent to output for " << rem << std::endl;
            auto it =
                std::find(products_requested_separated.begin(), products_requested_separated.end(), rem);
            products_requested_separated.erase(it);
        }
        for (const auto &prod_exp_3d : add_3d_expanded) {
            products_requested_separated.push_back(prod_exp_3d);
        }
    }
}

/**
 * @brief
 *  readL2, openL2 uses the maps from this module. If the module is not initialized, the  function read and
 * process an unexpanded L2 products
 * @param input_file_name path to nc file
 * @param products_requested products requested by the user, with delimiters
 * @param products_requested_separated list of the separated products
 * @param requested_wavelengths list of requested wavelength
 */
void set_mapper(const std::string &input_file_name, std::string &products_requested,
                std::vector<std::string> &products_requested_separated,
                const std::string &requested_wavelengths) {
    ini_3d = true;
    netCDF::NcFile l2_file(input_file_name, netCDF::NcFile::read);
    l2_file.getAtt("instrument").getValues(instrument);
    l2_file.getAtt("platform").getValues(platform);
    sensorID = instrumentPlatform2SensorId(instrument.c_str(), platform.c_str());
    // set 3d and 2d products
    {
        std::multimap<std::string, netCDF::NcVar> vars = find_all_variables(l2_file);
        for (const auto &var : vars) {
            const auto dims = var.second.getDims();
            if (dims.size() <= 1)
                continue;
            if (dims.size() == 3) {
                std::string dim_name = dims[2].getName();
                int dim_id = dims[2].getId();
                size_t dim_size = dims[2].getSize();
                // get variable name using dimension
                const std::string &variable_3d = dim_name;
                netCDF::NcVar nc_var = find_nc_variable_cpp(variable_3d, l2_file);
                int dim_3d_id = -1;
                size_t dim_3d_size = 0;
                if (!nc_var.isNull()) {
                    dim_3d_id = nc_var.getDims().at(0).getId();
                    dim_3d_size = nc_var.getDims().at(0).getSize();
                }
                if (nc_var.isNull() || dim_3d_id != dim_id || dim_3d_size != dim_size) {
                    netCDF::NcGroup parent_grp = var.second.getParentGroup();
                    nc_var = find_nc_variable_cpp("wavelength", parent_grp);
                }
                // third dimension is not wavelength
                if (nc_var.isNull()) {
                    std::cout << "Variable " << var.second.getName()
                              << " is 3D dimensional but no suitable wavelength set was found" << std::endl;
                    continue;
                }
                std::vector<int> waves;
                std::vector<float> waves_original;
                if (nc_var.getType() == netCDF::ncFloat) {
                    std::vector<float> temp_(dim_size);
                    nc_var.getVar(temp_.data());
                    waves = vector_rounder<float, int>(temp_);
                    std::copy(temp_.begin(), temp_.end(), std::back_inserter(waves_original));
                } else if (nc_var.getType() == netCDF::ncDouble) {
                    std::vector<double> temp_(dim_size);
                    nc_var.getVar(temp_.data());
                    waves = vector_rounder<double, int>(temp_);
                    std::copy(temp_.begin(), temp_.end(), std::back_inserter(waves_original));
                } else if (nc_var.getType() == netCDF::ncInt) {
                    waves.resize(dim_size);
                    nc_var.getVar(waves.data());
                    std::copy(waves.begin(), waves.end(), std::back_inserter(waves_original));
                }
                three_dims_prod_wavelength[var.first] = waves;
                three_d_original_prod_wavelength[var.first] = waves_original;
                for (auto &wv_ : waves)
                    wavelength_all.insert(wv_);
            }
            original_l2_products.insert(var.first);
        }
    }
    // read dimensions and set wavelength needed.
    for (auto &var_dim : three_dims_prod_wavelength) {
        if (requested_wavelengths != "ALL") {
            std::unordered_map<int, int> wv_3d_val_index;
            std::vector<int> &wavelengths_3d = var_dim.second;
            int count = 0;
            std::for_each(wavelengths_3d.begin(), wavelengths_3d.end(), [&](const int &val) {
                wv_3d_val_index.insert({val, count});
                count++;
            });
            std::vector<int> wv_requested_indexes;
            wavelengths_3d.clear();
            parse_wv_list(requested_wavelengths, wv_3d_val_index, wavelengths_3d, wv_requested_indexes,
                          false);
        }
    }
    {
        boost::algorithm::split(products_requested_separated, products_requested, boost::is_any_of(", "),
                                boost::token_compress_on);  // products_requested
        productInfo_t *info = allocateProductInfo();
        products_requested_separated = create_min_max_values(products_requested_separated, info);
        // l2->l2 3d wavelength expansion
        expand_l2_requested(products_requested_separated, info);
        {
            for (const auto &prod : products_requested_separated) {
                if (original_l2_products.count(prod)) {
                    l3d_to_l2_prod_map[prod] = prod;
                    l3d_indexes[prod] = {0, 1};
                    continue;
                }
                {
                    int wv_num = -1;
                    std::string original_2_prod_name = get_3d_product_name(prod, &wv_num, info);
                    // original_2_prod_name should be a 3D poduct
                    if (original_l2_products.count(original_2_prod_name)) {
                        // const auto wv_str = parsed_product.at(len_of_underscore_sep - 1);
                        try {
                            // const int wv_num = std::stoi(wv_str);
                            int index = -1;
                            if (three_dims_prod_wavelength.count(original_2_prod_name)) {
                                std::vector<int> &wavelengths_3d =
                                    three_dims_prod_wavelength[original_2_prod_name];
                                auto it = std::find(wavelengths_3d.begin(), wavelengths_3d.end(), wv_num);
                                if (it == wavelengths_3d.end()) {
                                    EXIT_LOG(std::cerr << "--Error-- : Wavelength not found; product = "
                                                       << prod << ", wavelength =  " << wv_num
                                                       << ", L2 file = " << input_file_name << std::endl)
                                }
                                index = it - wavelengths_3d.begin();
                                l3d_indexes[prod] = {index, index + 1};
                                l3d_to_l2_prod_map[prod] = original_2_prod_name;
                                std::vector<float> &wavelengths_original =
                                    three_d_original_prod_wavelength[original_2_prod_name];
                                l3d_wave_attr[prod] = wavelengths_original[index];
                            } else {
                                EXIT_LOG(std::cerr << "--Error-- : product = " << prod << " must be 3d"
                                                   << std::endl)
                            }

                        } catch (const std::exception &e) {
                            EXIT_LOG(std::cerr << e.what() << '\n')
                        }
                    } else {
                        EXIT_LOG(std::cerr << "--Error-- : Product not found = " << prod
                                           << ", L2 file = " << input_file_name << std::endl)
                    }
                }
            }
        }
        freeProductInfo(info);
    }

    l2_file.close();
}

int32_t get_l2prod_index(const l2_prod &l2, const char *prodname) {
    int32_t index;
    for (index = 0; index < l2.nprod; index++)
        if (strcmp(prodname, l2.prodname[index]) == 0)
            break;
    if (index == l2.nprod)
        index = -1;
    return index;
}

void set_prodname_3d_to_l2(const std::vector<std::string> &prodparam, l2_prod &l2_str,
                           std::vector<std::string> &l2_prodname, std::vector<std::string> &l3_prodname,
                           std::vector<int32_t> &thirdDimId, std::vector<float> &min_value,
                           std::vector<float> &max_value) {
    for (const auto &prod_name : prodparam) {
        int32_t l2_iprod = get_l2prod_index(l2_str, prod_name.c_str());
        if (l2_iprod == -1) {
            printf("--Error-- Product: %s was not found in the L2 file: %s. See line %d in %s\n",
                   prod_name.c_str(), l2_str.filename, __LINE__, __FILE__);
            exit(EXIT_FAILURE);
        }
        int32_t tmpThirdDim = l2_str.thirdDim[l2_iprod];
        if (tmpThirdDim == 1) {
            l2_prodname.push_back(prod_name);
            l3_prodname.push_back(prod_name);
            thirdDimId.push_back(0);
            min_value.push_back(min_max_values.at(prod_name).first);
            max_value.push_back(min_max_values.at(prod_name).second);
        } else {
            EXIT_LOG(std::cerr << "--Error--: All products must be 2D. Error for " << prod_name << std::endl)
        }
    } /* iprod loop */
}

bool set_l2_flags_use(const std::string &flagsuse) {
    if (!flagsuse.empty())
        use_l2_flags = true;
    return use_l2_flags;
}

void l3_l2_conversion(char **inp3, char **out2) {
    if (ini_3d) {
        *out2 = &l3d_to_l2_prod_map.at(*inp3)[0];
    } else {
        *out2 = *inp3;
    }
}

float l3_attr(const std::string &inp3) {
    if (l3d_wave_attr.count(inp3))
        return l3d_wave_attr.at(inp3);
    else
        return BAD_FLT;
}

void l3_l2_index(char **inp3, int *start, int *count) {
    if (ini_3d) {
        *start = l3d_indexes.at(*inp3).first;
        *count = l3d_indexes.at(*inp3).second - l3d_indexes.at(*inp3).first;
    }
}
int get_set_flag() {
    return ini_3d ? 1 : 0;
}
int get_l2_flag_use() {
    int ans = use_l2_flags ? 1 : 0;
    if (!ini_3d)
        return true;
    return ans;
}
