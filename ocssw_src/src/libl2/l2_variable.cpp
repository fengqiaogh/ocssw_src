#include "l2_variable.hpp"
#include <sensorInfo.h>
#include "find_variable.hpp"
#include "l2_utils.hpp"
#include "expand3D.hpp"

L2Variable::L2Variable(netCDF::NcFile &l2_file) : L2Variable() {
    read(l2_file);
}

L2Variable::L2Variable(netCDF::NcFile &l2_file, std::string &products_requested,
                       std::vector<std::string> &products_requested_separated,
                       const std::string &requested_wavelengths, const std::string &file_name) : L2Variable() {
    read(l2_file, products_requested, products_requested_separated, requested_wavelengths, file_name);
}

void L2Variable::read(netCDF::NcFile &l2_file) {
    // get the instrument and platform from the file
    std::string instrument, platform;
    try {
        l2_file.getAtt("instrument").getValues(instrument);
        l2_file.getAtt("platform").getValues(platform);
    } catch (netCDF::exceptions::NcException const &e) {
        printf("Warning: failed getting instrument and platform from file: %s\n", e.what());
    }
    // get the sensor ID
    sensorID = instrumentPlatform2SensorId(instrument.c_str(), platform.c_str());
    // set 3d and 2d products
    std::multimap<std::string, netCDF::NcVar> vars = find_all_variables(l2_file);
    for (const auto &[var_name, var_nc] : vars) {
        const auto dims = var_nc.getDims();
        if (dims.size() <= 1)
            continue;
        if (dims.size() == 3) {
            std::string dim_name = dims[2].getName();
            int dim_id = dims[2].getId();
            size_t dim_size = dims[2].getSize();
            // get variable name using dimension
            const std::string &variable_3d = dim_name;
            netCDF::NcGroup parent_grp = var_nc.getParentGroup();
            netCDF::NcVar nc_var = find_nc_variable_cpp(variable_3d, parent_grp);
            int dim_3d_id = -1;
            size_t dim_3d_size = 0;
            if (!nc_var.isNull()) {
                dim_3d_id = nc_var.getDims().at(0).getId();
                dim_3d_size = nc_var.getDims().at(0).getSize();
            }
            if (nc_var.isNull() || dim_3d_id != dim_id || dim_3d_size != dim_size) {
                nc_var = find_nc_variable_cpp(variable_3d, l2_file);
            }
            if (!nc_var.isNull()) {
                dim_3d_id = nc_var.getDims().at(0).getId();
                dim_3d_size = nc_var.getDims().at(0).getSize();
            }
            // third dimension is not wavelength
            if (nc_var.isNull() || dim_3d_id != dim_id || dim_3d_size != dim_size) {
#ifdef DEBUG_LOG
                std::cout << "WARNING: Variable " << var_nc.getName()
                          << " is 3D dimensional but no suitable wavelength set was found" << std::endl;
#endif
                continue;
            } else {
                prod_wavenames[get_full_nc_path(nc_var)].insert(var_name);
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
            three_dims_prod_wavelength[var_name] = waves;
            three_d_original_prod_wavelength[var_name] = waves_original;
            for (auto &wv_ : waves)
                wavelength_all.insert(wv_);
        }
        std::map<std::string, netCDF::NcVarAtt> attributes = var_nc.getAtts();
        // check if valid_min, valid_max are defined
        products_units[var_name] = "undefined";  // set undefined
        if (sensorID != -1) {
            if (findProductInfo(var_name.c_str(), sensorID, info)) {
                products_units[var_name] = info->units;
            }
        }
        if (attributes.count("units")) {
            attributes.at("units").getValues(products_units.at(var_name));
        }
        if (attributes.count("wavelength")) {
            float wave_var;
            attributes.at("wavelength").getValues(&wave_var);
            l3d_wave_attr[var_name] = wave_var;
        }
        if (attributes.count("valid_min") > 0 &&
            attributes.count("valid_max") > 0) {  // need to think about it
            auto validMinAttr = attributes.at("valid_min");
            auto validMaxAttr = attributes.at("valid_max");
            float scale_factor = 1.0;
            float offset = 0.0;
            if (attributes.count("scale_factor")) {
                scale_factor = get_attr<float>(attributes.at("scale_factor"));
                offset = get_attr<float>(attributes.at("add_offset"));
            }
            if (validMaxAttr.getType() == netCDF::ncDouble) {
                validMin_validMax_nc_file[var_name].first =
                    get_attr<double>(validMinAttr) * scale_factor + offset;
                validMin_validMax_nc_file[var_name].second =
                    get_attr<double>(validMaxAttr) * scale_factor + offset;
            }
            if (validMaxAttr.getType() == netCDF::ncFloat) {
                validMin_validMax_nc_file[var_name].first =
                    get_attr<float>(validMinAttr) * scale_factor + offset;
                validMin_validMax_nc_file[var_name].second =
                    get_attr<float>(validMaxAttr) * scale_factor + offset;
            }
            if (validMaxAttr.getType() == netCDF::ncInt) {
                validMin_validMax_nc_file[var_name].first =
                    get_attr<int>(validMinAttr) * scale_factor + offset;
                validMin_validMax_nc_file[var_name].second =
                    get_attr<int>(validMaxAttr) * scale_factor + offset;
            }
            if (validMaxAttr.getType() == netCDF::ncShort) {
                validMin_validMax_nc_file[var_name].first =
                    get_attr<short>(validMinAttr) * scale_factor + offset;
                validMin_validMax_nc_file[var_name].second =
                    get_attr<short>(validMaxAttr) * scale_factor + offset;
            }
        }
        original_l2_products.insert(var_name);
    }
}

void L2Variable::read(netCDF::NcFile &l2_file, std::string &products_requested,
                      std::vector<std::string> &products_requested_separated,
                      const std::string &requested_wavelengths, const std::string &file_name) {
    read(l2_file);
    // get the product name from the file
    try {
        l2_file.getAtt("product_name").getValues(input_file_name);
    } catch (netCDF::exceptions::NcException const &e) {
        input_file_name = file_name;
    }

    // copy the integerized wavelengths set for further usage
    std::unordered_map<std::string, std::vector<int>> three_dims_prod_wavelength_all = three_dims_prod_wavelength;
    // set the requested subsets for each product
    if (requested_wavelengths != "ALL") {
        std::unordered_map<std::string, std::string> prod_subset_waves =
            parse_wave_subsets(requested_wavelengths);
        for (auto &var_dim : three_dims_prod_wavelength) {
            if (prod_subset_waves[var_dim.first] != "ALL") {
                std::string requestd_subset = prod_subset_waves[var_dim.first];
                std::unordered_map<int, int> wv_3d_val_index;
                std::vector<int> &wavelengths_3d = var_dim.second;
                int count = 0;
                std::for_each(wavelengths_3d.begin(), wavelengths_3d.end(), [&](const int &val) {
                    wv_3d_val_index.insert({val, count});
                    count++;
                });
                std::vector<int> wv_requested_indexes;
                wavelengths_3d.clear();
                parse_wv_list(requestd_subset, wv_3d_val_index, wavelengths_3d, wv_requested_indexes,
                              false);
            }
        }
    }

    boost::algorithm::split(products_requested_separated, products_requested,
                            boost::is_any_of(product_delimiter),
                            boost::token_compress_on);  // products_requested
    products_requested_separated = create_min_max_values(products_requested_separated);

    expand_l2_requested(products_requested_separated);

    for (const auto &prod : products_requested_separated) {
        if (original_l2_products.count(prod)) {
            l3d_to_l2_prod_map[prod] = prod;
            l3d_indexes[prod] = {0, 1};
            continue;
        }
        {
            int wv_num = -1;
            std::string original_2_prod_name = get_3d_product_name(prod, &wv_num);
            if (wv_num == -1) {
                EXIT_LOG(std::cerr << "--Error-- : Wavelength not found or product = " << prod
                                   << ",  does not exist in L2 file = " << input_file_name << std::endl)
            }
            // original_2_prod_name should be a 3D product
            if (original_l2_products.count(original_2_prod_name)) {
                try {
                    int index = -1;
                    if (three_dims_prod_wavelength.count(original_2_prod_name)) {
                        std::vector<int> &wavelengths_3d =
                            three_dims_prod_wavelength_all[original_2_prod_name];
                        auto it = std::find(wavelengths_3d.begin(), wavelengths_3d.end(), wv_num);
                        if (it == wavelengths_3d.end()) {
                            EXIT_LOG(std::cerr << "--Error-- : Wavelength not found; product = " << prod
                                               << ", wavelength =  " << wv_num
                                               << ", L2 file = " << input_file_name << std::endl)
                        }
                        index = it - wavelengths_3d.begin();
                        l3d_indexes[prod] = {index, index + 1};
                        l3d_to_l2_prod_map[prod] = original_2_prod_name;
                        std::vector<float> &wavelengths_original =
                            three_d_original_prod_wavelength[original_2_prod_name];
                        l3d_wave_attr[prod] = wavelengths_original[index];
                    } else {
                        EXIT_LOG(std::cerr << "--Error-- : product = " << prod << " must be 3d" << std::endl)
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

std::vector<std::string> L2Variable::create_min_max_values(
    const std::vector<std::string> &products_requested_separated) {
    std::vector<std::string> output_expanded_products;

    for (const auto &prod_with_min_max : products_requested_separated) {
        std::vector<std::string> prod_data;
        float max_val_preset, min_val_preset;
        boost::algorithm::split(prod_data, prod_with_min_max, boost::is_any_of(min_max_token),
                                boost::token_compress_on);  // products_requested}

        if (findProductInfo(prod_data.at(0).c_str(), sensorID, info)) {
            min_val_preset = info->validMin;
            max_val_preset = info->validMax;
        } else {
            min_val_preset = 0;
            max_val_preset = std::numeric_limits<float>::max();
        }
        std::string prodname = prod_data.at(0);
        if (validMin_validMax_nc_file.count(prodname)) {
            min_val_preset = validMin_validMax_nc_file.at(prodname).first;
            max_val_preset = validMin_validMax_nc_file.at(prodname).second;
        }
        if (prod_data.size() > 1) {
            output_expanded_products.push_back(prod_data.at(0));
            std::vector<std::string> min_max_sep;
            boost::algorithm::split(min_max_sep, prod_data.at(1), boost::is_any_of(min_max_range_delimiter),
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

void L2Variable::expand_l2_requested(std::vector<std::string> &products_requested_separated) {
    std::unordered_set<std::string> remove;
    std::vector<std::string> add_3d_expanded;
    for (auto &prod : products_requested_separated) {
        if (three_dims_prod_wavelength.count(prod) > 0) {
            remove.insert(prod);
            std::vector<int> &wavelengths_3d = three_dims_prod_wavelength[prod];
            for (const auto &wv : wavelengths_3d) {
                std::string exp3d_name = get_3d_product_name(prod, wv);  // prod + "_" + std::to_string(wv);
                add_3d_expanded.push_back(exp3d_name);
                min_max_values[exp3d_name] = min_max_values.at(prod);
                products_units[exp3d_name] = products_units.at(prod);
                if (products_units[exp3d_name] == "undefined") {
                    if (findProductInfo(exp3d_name.c_str(),sensorID,info))
                        products_units[exp3d_name] = info->units;
                }
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

std::string L2Variable::get_3d_product_name(const std::string &prod_3d, const std::string &wv) {
    if (findProductInfo(prod_3d.c_str(), sensorID, info)) {
        std::string prefix = info->prefix;
        std::string suffix = info->suffix;
        return prefix + "_" + wv + suffix;
    }
    return prod_3d + "_" + wv;
}

std::string L2Variable::get_3d_product_name(const std::string &prod_3d, int wv) {
    return get_3d_product_name(prod_3d, std::to_string(wv));
}

std::string L2Variable::get_3d_product_name(const std::string &prod_3d, int *wv) {
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
std::string L2Variable::l3_l2_conversion(const std::string &l3_name) {
    return l3d_to_l2_prod_map.at(l3_name);
}

std::unordered_map<std::string, std::string> L2Variable::get_units_list() {
    return products_units;
}
std::unordered_map<std::string, std::pair<float, float>> L2Variable::get_min_max_values() {
    return min_max_values;
}

size_t L2Variable::slice_3d_index(std::string l3_name) {
    return l3d_indexes.at(l3_name).first;
}

std::map<std::string, float> L2Variable::l3_attrs(const std::string &prod_name) {
    std::map<std::string, float> out;
    if (l3d_wave_attr.count(prod_name)) {
        out["wavelength"] = l3d_wave_attr.at(prod_name);
    }
    if (min_max_values.count(prod_name)) {
        out["valid_min"] = min_max_values.at(prod_name).first;
        out["valid_max"] = min_max_values.at(prod_name).second;
    }
    return out;
}

std::unordered_map<std::string,std::string> L2Variable:: parse_wave_subsets(
    const std::string &subsets_wave_requested) {
    // Map to store product->wavelength subset mappings
    std::unordered_map<std::string, std::string> prod_subset_waves;
    
    // Split input string on semicolons to get individual subset specifications
    std::vector<std::string> subsets_wave_requested_list;
    boost::split(subsets_wave_requested_list, subsets_wave_requested, boost::is_any_of(";"),
                 boost::algorithm::token_compress_on);

    // Process each subset specification
    for (auto &subset_wave_requested : subsets_wave_requested_list) {
        // Split on equals to separate wavelength name from subset value
        std::vector<std::string> subset_wave_requested_list;
        boost::split(subset_wave_requested_list, subset_wave_requested, boost::is_any_of("="),
                     boost::algorithm::token_compress_on);

        // Validate format
        if (subset_wave_requested_list.size() != 2 && subsets_wave_requested_list.size() != 1) {
            EXIT_LOG(std::cerr << "-E- : For " << subset_wave_requested
                               << " subsetting format should contain one equals sign =";);
        } else if (subset_wave_requested_list.size() == 1) {
            // If no equals, apply subset to all products
            for (auto &[prod,_]: three_dims_prod_wavelength) {
                prod_subset_waves[prod] = subset_wave_requested;
            }
        } else {
            // Process wavelength name and subset value
            std::string wave_name = subset_wave_requested_list[0];
            std::string subset = subset_wave_requested_list[1];
            
            // Validate wavelength exists
            std::string wave_path = find_key_nc_path(prod_wavenames,wave_name);
            if (wave_path.empty()) {
                EXIT_LOG(std::cout << "-E- : Wavelength or parent group " << wave_name << " not found in the file "
                                   << input_file_name);
            }
            // Apply subset to all products with this wavelength
            for (auto &prod : prod_wavenames[wave_path]) {
                prod_subset_waves[prod] = subset;
            }
        }
    }

    // Set default "ALL" for any products without explicit subset
    for (auto &[prod,_] : three_dims_prod_wavelength) {
        if (prod_subset_waves.count(prod) == 0) {
            prod_subset_waves[prod] = "ALL";
        }
    }
    return prod_subset_waves;
}