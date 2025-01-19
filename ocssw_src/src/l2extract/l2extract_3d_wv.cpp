#include "l2extract_3d_wv.h"

#include <expand3D.hpp>


/**
 * @brief Get the wv3 indexes object
 *
 * @param inFile - l2 file to process
 * @param waveList - wv list supplied by the user.
 * @param wv_vals - wv values to output
 * @param waveIndexesOut - wv indexes to output
 * @param prodList - product list (reserved for the future use)
 */
void getWave3Indexes(const char *inFile, const char *waveList, int **waveOut, int **waveIndexesOut,
                     int *waveNumPass, const char *prodList) {
    try {
        netCDF::NcFile l2_file(inFile, netCDF::NcFile::read);
        auto grp = l2_file.getGroup("sensor_band_parameters");
        {
            if (grp.getVar("wavelength_3d").isNull()) {
                if (std::string(waveList).length() > 0) {
                    std::cerr << "No wavelength 3D are found in the L2 file " << inFile
                              << " but the user supplied it";
                    exit(EXIT_FAILURE);
                }
                return;
            }
            if (std::string(waveList).length() == 0) {
                std::cout << "No wavelength_3d list is supplied by the user. "
                             "All wavelength 3D will be in output "
                          << std::endl;
                return;
            }
        }
        auto dimsWave3d = grp.getVar("wavelength_3d").getDims();
        size_t wave3dLen = 1;
        std::for_each(dimsWave3d.begin(), dimsWave3d.end(),
                      [&](const netCDF::NcDim &dim) { wave3dLen *= dim.getSize(); });
        std::vector<int> wave3dVals(wave3dLen, 0);
        grp.getVar("wavelength_3d").getVar(wave3dVals.data());
        std::unordered_map<int, int> wave3dValsLut;
        int count = 0;
        std::for_each(wave3dVals.begin(), wave3dVals.end(), [&](const int &val) {
            wave3dValsLut.insert({val, count});
            count++;
        });
        static std::vector<int> waveRequestedVals, waveRequestedIndexes;
        parse_wv_list(waveList, wave3dValsLut, waveRequestedVals, waveRequestedIndexes);
        l2_file.close();
        *waveOut = waveRequestedVals.data();
        *waveIndexesOut = waveRequestedIndexes.data();
        *waveNumPass = (int)waveRequestedVals.size();
    } catch (netCDF::exceptions::NcException const &e) {
        e.what();
        exit(EXIT_FAILURE);
    } catch (std::exception const &e) {
        std::cerr << "Error opening file for reading: " << e.what() << std::endl;
        exit(EXIT_FAILURE);
    }
}
