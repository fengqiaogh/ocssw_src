#include "l2_reader.hpp"
#include "find_variable.hpp"
#include "timeutils.h"
#include <boost/algorithm/string.hpp>

L2_Reader::~L2_Reader() {
}

void L2_Reader::iniGeolocation() {
    std::string error_message = "In file: " + file_path;
    find_lat_lon<netCDF::NcFile, ScaledNcVar>(file_nc, lat_nc, lon_nc, error_message);
    NC_CHECK(first_dimension = lat_nc.getDims()[0].getSize();
             second_dimension = lat_nc.getDims()[1].getSize();)
    NC_CHECK(if (first_dimension != lon_nc.getDims()[0].getSize() ||
                 second_dimension != lon_nc.getDims()[1].getSize()) {EXIT_LOG("--Error--: Latitude/Longitude dimensions mismatch in file " +  file_path) })
}

int L2_Reader::setVariables(std::string& product_list, std::string& wave_list,
                            std::vector<std::string>& products_requested) {
    if (first_dimension == 0 || second_dimension == 0) {
       iniGeolocation();
    }
    if (product_list.empty() || product_list == "ALL") {
        product_list = "";
        // look for geophysical_data group first
        netCDF::NcGroup geophysical_data = file_nc.getGroup("geophysical_data");
        std::multimap<std::string, netCDF::NcVar> geo_vars;
        if (!geophysical_data.isNull()) {
            geo_vars = find_all_variables(geophysical_data);
        } else {
            geo_vars = find_all_variables(file_nc);
        }

        // check dimensions:
        for (const auto& [var_name, var_nc] : geo_vars) {
            std::vector<netCDF::NcDim> dims = var_nc.getDims();
            if (dims[0].getSize() == first_dimension && dims[1].getSize() == second_dimension &&
                dims.size() >= 2 && dims.size() <= 3 && var_name != "l2_flags" &&
                lat_possible_names.count(var_name) == 0 && lon_possible_names.count(var_name) == 0) {
                if (product_list.empty()) {
                    product_list = var_name;
                } else {
                    product_list += "," + var_name;
                }
            }
        }
    }
    if (product_list.empty()) {
        std::cerr << "Warning: No L2 products are found.\n";
        return 1;
    }
#ifdef DEBUG_LOG
    std::cout << "Products to read: " << product_list << std::endl;
#endif
    L2Variable l2_variable(file_nc, product_list, products_requested, wave_list, file_path);
    std::unordered_map<std::string, std::string> product_units_list = l2_variable.get_units_list();
    std::unordered_map<std::string, std::pair<float, float>> min_max_product_list =
        l2_variable.get_min_max_values();
    for (auto& var_l3_name : products_requested) {
        std::string l2_name_original = l2_variable.l3_l2_conversion(var_l3_name);
        index.push_back(l2_variable.slice_3d_index(var_l3_name));
        l2_products.push_back(find_nc_variable_cpp<netCDF::NcFile, ScaledNcVar>(l2_name_original, file_nc));
        products_list.push_back(l2_name_original);
        min_values.push_back(min_max_product_list.at(var_l3_name).first);
        max_values.push_back(min_max_product_list.at(var_l3_name).second);
        units.push_back(product_units_list.at(l2_name_original));
        product_attributes[var_l3_name] = l2_variable.l3_attrs(var_l3_name);
    }
    geo_bounds = Geospatialbounds(file_nc);
    if (!geo_bounds.get_time_coverage_start().empty() && !geo_bounds.get_time_coverage_end().empty()) {
        start_time_unix = isodate2unix(convert_string(geo_bounds.get_time_coverage_start()));
        end_time_unix = isodate2unix(convert_string(geo_bounds.get_time_coverage_end()));
    }
    // set number of l3 products
    number_of_l3_products = products_requested.size();
    // resize cached products
    start_cached_products.resize(number_of_l3_products);
    end_cached_products.resize(number_of_l3_products);
    data_products_cached.resize(number_of_l3_products);
    return 0;
}

int L2_Reader::ini_l2_flags() {
    l2_flags_nc = find_nc_variable_cpp("l2_flags", file_nc);
    if (l2_flags_nc.isNull())
        return 1;
    NC_CHECK(l2_meaning_attr = l2_flags_nc.getAtt("flag_meanings");
             l2_mask_attr = l2_flags_nc.getAtt("flag_masks"); l2_meaning_attr.getValues(l2_meaning);
             l2_mask_bits.resize(l2_mask_attr.getAttLength()); l2_mask_attr.getValues(l2_mask_bits.data()))
    flag_l2_set = true;
    return 0;
}

int L2_Reader::ini_quality_flags(const std::string& qual_product) {
    quality_flags_nc = find_nc_variable_cpp(qual_product, file_nc);
    if (quality_flags_nc.isNull())
        return 1;
    quality_flag_set = true;
    qual_product_name = qual_product;
    return 0;
}

void L2_Reader::readL2dataScan(std::vector<float*>& data, size_t scan, const std::vector<uint8_t>& mask) {
    std::vector<size_t> start = {scan, 0};
    std::vector<size_t> count = {1, second_dimension};
    readL2data(data, start, count, mask, false);
}

void L2_Reader::getDimensions(size_t& lines, size_t& pixels) const {
    lines = first_dimension;
    pixels = second_dimension;
}

std::string L2_Reader::get_platform() const {
    return geo_bounds.get_platform();
};
std::string L2_Reader::get_instrument() const {
    return geo_bounds.get_instrument();
};

void L2_Reader::get_escan_bscan_row(std::vector<int>& bscan_row, std::vector<int>& escan_row, float dlat) {
    std::vector<float> slat(first_dimension);
    std::copy(geo_bounds.get_slat(), geo_bounds.get_slat() + first_dimension, slat.begin());
    std::vector<float> elat(first_dimension);
    std::copy(geo_bounds.get_elat(), geo_bounds.get_elat() + first_dimension, elat.begin());
    escan_row.resize(first_dimension);
    bscan_row.resize(first_dimension);
    for (size_t jsrow = 0; jsrow < first_dimension; jsrow++) {
        escan_row[jsrow] = (int32_t)((90 + elat[jsrow]) / dlat);
        bscan_row[jsrow] = (int32_t)((90 + slat[jsrow]) / dlat);

        if (escan_row[jsrow] > bscan_row[jsrow]) {
            int k = escan_row[jsrow];
            escan_row[jsrow] = bscan_row[jsrow];
            bscan_row[jsrow] = k;
        }
        escan_row[jsrow] -= padding;
        bscan_row[jsrow] += padding;
    }
}

Geospatialbounds& L2_Reader::get_geospatial() {
    return geo_bounds;
}

double L2_Reader::get_start_time() const {
    return start_time_unix;
}
double L2_Reader::get_end_time() const {
    return end_time_unix;
}

void L2_Reader::reset_cache() {
    std::for_each(start_cached_products.begin(), start_cached_products.end(),
                  [](cache_bounds& el) { el.reset(); });
    std::for_each(end_cached_products.begin(), end_cached_products.end(),
                  [](cache_bounds& el) { el.reset(); });
    start_cached_latitude.reset();
    end_cached_latitude.reset();
    start_cached_longitude.reset();
    end_cached_longitude.reset();
    start_cached_l2flags.reset();
    end_cached_l2flags.reset();
    start_cached_quality.reset();
    end_cached_quality.reset();
    std::for_each(data_products_cached.begin(), data_products_cached.end(),
                  [](std::vector<float>& el) { free_vector(el); });
    free_vector(l2flags_cached);
    free_vector(qualityflags_cached);
    free_vector(latitude_cached);
    free_vector(longitude_cached);
    if (!file_nc.isNull()) {
        file_nc.close();
    }
}

L2_Reader::Cache_data L2_Reader::set_cached_data(int prodId, std::vector<size_t>& start,
                                                 std::vector<size_t>& count, cache_bounds& start_cached,
                                                 cache_bounds& end_cached, const std::vector<uint8_t>& mask) {
    Cache_data cached_data(start, count);
    bool& read_l2_file = cached_data.read_l2_file;
    size_t& bscan = cached_data.bscan;
    size_t& bpixel = cached_data.bpixel;
    size_t& escan = cached_data.escan;
    size_t& epixel = cached_data.epixel;
    size_t& bscan_requested = cached_data.bscan_requested;
    size_t& escan_requested = cached_data.escan_requested;
    size_t& bpixel_requested = cached_data.bpixel_requested;
    size_t& epixel_requested = cached_data.epixel_requested;
    std::vector<size_t>& start_nc = cached_data.start_nc;
    std::vector<size_t>& count_nc = cached_data.count_nc;
    if (prodId >= 0) {
        start_nc.push_back(index[prodId]);
        count_nc.push_back(1);
    }
    if (start_cached.has_value() && end_cached.has_value()) {
        bscan = start_cached.value()[0];
        bpixel = start_cached.value()[1];
        escan = end_cached.value()[0];
        epixel = end_cached.value()[1];

        if (bscan > bscan_requested || escan_requested > escan || bpixel_requested < bpixel ||
            epixel_requested > epixel) {
            if (!mask.empty()) {
                auto itb = std::find(mask.begin(), mask.end(), scan_within_the_group);
                if (itb != mask.end()) {
                    bscan = itb - mask.begin();
                    auto ite = std::find(mask.rbegin(), mask.rend(), scan_within_the_group);
                    escan = mask.rend() - ite;
                }
            }
            bscan = std::min(bscan_requested, bscan);
            bpixel = std::min(bpixel_requested, bpixel);
            escan = std::max(escan_requested, escan);
            epixel = std::max(epixel_requested, epixel);

            start_cached = std::array<size_t, 2>{bscan, bpixel};
            end_cached = std::array<size_t, 2>{escan, epixel};
            start_nc[0] = bscan;
            start_nc[1] = bpixel;
            count_nc[0] = escan - bscan;
            count_nc[1] = epixel - bpixel;
            read_l2_file = true;
        }
    } else {
        escan = escan_requested;
        bscan = bscan_requested;
        if (!mask.empty()) {
            auto itb = std::find(mask.begin(), mask.end(), scan_within_the_group);
            if (itb != mask.end()) {
                bscan = itb - mask.begin();
                auto ite = std::find(mask.rbegin(), mask.rend(), scan_within_the_group);
                escan = mask.rend() - ite;
            }
            bscan = std::min(bscan_requested, bscan);
            escan = std::max(escan_requested, escan);
        }
        bpixel = bpixel_requested;
        epixel = epixel_requested;
        start_nc[0] = bscan;
        count_nc[0] = escan - bscan;
        start_cached = std::array<size_t, 2>{bscan, bpixel};
        end_cached = std::array<size_t, 2>{escan, epixel};
        read_l2_file = true;
    }
    return cached_data;
}

void L2_Reader::readLatitudeScan(float** data, size_t scan, const std::vector<uint8_t>& mask) {
    std::vector<size_t> start = {scan, 0};
    std::vector<size_t> count = {1, second_dimension};
    readLatitude(data, start, count, mask, false);
}
void L2_Reader::readLongitudeScan(float** data, size_t scan, const std::vector<uint8_t>& mask) {
    std::vector<size_t> start = {scan, 0};
    std::vector<size_t> count = {1, second_dimension};
    readLongitude(data, start, count, mask, false);
}

void L2_Reader::readL2FlagsScan(int** data, size_t scan, const std::vector<uint8_t>& mask) {
    std::vector<size_t> start = {scan, 0};
    std::vector<size_t> count = {1, second_dimension};
    readL2Flags(data, start, count, mask, false);
}

void L2_Reader::readQualityFlagsScan(int32_t** data, size_t scan, const std::vector<uint8_t>& mask) {
    std::vector<size_t> start = {scan, 0};
    std::vector<size_t> count = {1, second_dimension};
    readQualityFlags(data, start, count, mask, false);
}

std::unordered_map<std::string, int> L2_Reader::get_l2_meaning_bit_dict() const {
    std::vector<std::string> separated_meanings;
    std::unordered_map<std::string, int> out;
    boost::split(separated_meanings, l2_meaning, boost::is_any_of(" "), boost::algorithm::token_compress_on);
    if (separated_meanings.size() != l2_mask_bits.size()) {
        EXIT_LOG(std::cerr << "-Error-: the number of L2 flag meanings doesn't equal to the number of bits")
    }
    for (size_t i = 0; i < separated_meanings.size(); ++i) {
        out[separated_meanings[i]] = l2_mask_bits[i];
    }
    return out;
}

std::vector<float> L2_Reader::get_min_value_product() const {
    return min_values;
};
std::vector<float> L2_Reader::get_max_value_product() const {
    return max_values;
};
std::vector<std::string> L2_Reader::get_units() const {
    return units;
};

void L2_Reader::reopenL2() {
    if (!file_nc.isNull()) {
        file_nc.close();
    }
    NC_CHECK(file_nc.open(file_path, netCDF::NcFile::read));
    for (size_t i = 0; i < l2_products.size(); i++) {
        l2_products[i] = find_nc_variable_cpp<netCDF::NcFile, ScaledNcVar>(products_list[i], file_nc);
    }
    NC_CHECK(find_lat_lon<netCDF::NcFile, ScaledNcVar>(file_nc, lat_nc, lon_nc));
    if (flag_l2_set) {
        NC_CHECK(l2_flags_nc = find_nc_variable_cpp("l2_flags", file_nc));
    }
    if (quality_flag_set) {
        NC_CHECK(quality_flags_nc = find_nc_variable_cpp(qual_product_name, file_nc));
    }
}

int32_t L2_Reader::find_product_index(const std::string& product_name, size_t& index) {
    auto it = std::find(products_list.begin(), products_list.end(), product_name);
    if (it == products_list.end())
        return 1;
    else {
        index = std::distance(products_list.begin(), it);
        return 0;
    }
}

std::map<std::string, float> L2_Reader::get_product_attributes(const std::string& product_name) const {
    return product_attributes.at(product_name);
}

netCDF::NcVar L2_Reader::get_variable(const std::string& name) {
    return find_nc_variable_cpp(name, file_nc);
}

void L2_Reader::set_cache_flag(bool cache_usage) {
    use_cache = cache_usage;
}
