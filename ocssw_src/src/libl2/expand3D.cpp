#include "expand3D.hpp"
#include "l2_utils.hpp"

/**
 * @brief Parses a colon-separated range of wavelength values and adds them to the requested lists
 * 
 * @param wv_range String containing wavelength range in format "nnn:nnn"
 * @param wv3d_list Map of wavelength values to their indices
 * @param wv_requested_vals Vector to store requested wavelength values
 * @param wv_requested_indexes Vector to store indices of requested wavelengths
 * @param wv_requested_vals_lookup Set to track which wavelengths have been requested
 * @throws Exits if range format is invalid or wavelengths not found
 */
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
}