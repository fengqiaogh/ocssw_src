#include "l2extract_3d_wv.h"

#include <expand3D.hpp>


/**
 * @brief Get the wv3 indexes object
 *
 * @param inp_file - l2 file to process
 * @param wv_list - wv list supplied by the user.
 * @param wv_vals - wv values to output
 * @param wv_indexes_out - wv indexes to output
 * @param prod_list - product list (reserved for the future use)
 */
void get_wv3_indexes(const char *inp_file, const char *wv_list,
                     int **wv_vals_out, int **wv_indexes_out,
                     int *wv_num_to_pass, const char *prod_list) {
    try {
        netCDF::NcFile l2_file(inp_file, netCDF::NcFile::read);
        auto grp = l2_file.getGroup("sensor_band_parameters");
        {
            if (grp.getVar("wavelength_3d").isNull()) {
                if (std::string(wv_list).length() > 0) {
                    std::cerr << "No wavelength 3D are found in the L2 file "
                              << inp_file << " but the user supplied it";
                    exit(EXIT_FAILURE);
                }
                return;
            }
            if (std::string(wv_list).length() == 0) {
                std::cout << "No wavelength_3d list is supplied by the user. "
                             "All wavelength 3D will be in output "
                          << std::endl;
                return;
            }
        }
        auto dims_wv_3d = grp.getVar("wavelength_3d").getDims();
        size_t wv_3d_len = 1;
        std::for_each(
            dims_wv_3d.begin(), dims_wv_3d.end(),
            [&](const netCDF::NcDim &dim) { wv_3d_len *= dim.getSize(); });
        std::vector<int> wv_3d_vals(wv_3d_len, 0);
        grp.getVar("wavelength_3d").getVar(wv_3d_vals.data());
        std::unordered_map<int, int> wv_3d_vals_lut;
        int count = 0;
        std::for_each(wv_3d_vals.begin(), wv_3d_vals.end(),
                      [&](const int &val) {
                          wv_3d_vals_lut.insert({val, count});
                          count++;
                      });
        static std::vector<int> wv_requested_vals,wv_requested_indexes;
        parse_wv_list(wv_list, wv_3d_vals_lut,wv_requested_vals,wv_requested_indexes);
        l2_file.close();
        *wv_vals_out = wv_requested_vals.data();
        *wv_indexes_out = wv_requested_indexes.data();
        *wv_num_to_pass = (int)wv_requested_vals.size();
    } catch (netCDF::exceptions::NcException const &e) {
        e.what();
        exit(EXIT_FAILURE);
    } catch (std::exception const &e) {
        std::cerr << "Error opening file for reading: " << e.what()
                  << std::endl;
        exit(EXIT_FAILURE);
    }
}
