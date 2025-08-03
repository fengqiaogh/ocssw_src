#ifndef SST_CLOUD_MASK_H
#define SST_CLOUD_MASK_H

#include "l2_struc.h"
#include "l1q_struc.h"
#include <vector>
#include <iostream>
#include "l2_struc.h"
#include "l1q_struc.h"
#include <boost/variant.hpp>
#include <string>
#include <utility>
#include "flags_sst.h"
#include <assert.h>
#include <algorithm>
#include <fstream>
#include <cstdlib>
#include <rapidjson/document.h>
#include <rapidjson/istreamwrapper.h>
#include <rapidjson/writer.h>
#include <rapidjson/stringbuffer.h>
#include <rapidjson/ostreamwrapper.h>
#include <unordered_set>
#include <chrono>
#include <unordered_map>
#include <memory>
#include <netcdf>
#include "sst_dsdi.h"
#include <boost/algorithm/string.hpp>
#include <boost/unordered_map.hpp>
// parameters
namespace envset {
    /**
     * @brief Returns OCDATAROOT
     */
    std::string get_ocdata_root();

    /**
     * @brief returns the sst coefs paths to netCDF files depending on requested product
     */
    const std::unordered_map<std::string, std::string (*)(const instr *)> get_sst_coeffs_path = {{std::string("SST"),
                                                                                                         [](const instr *input) {
                                                                                                             return std::string(
                                                                                                                     input->sstcoeffile);
                                                                                                         }},
                                                                                                 {std::string("SST4"),
                                                                                                         [](const instr *input) {
                                                                                                             return std::string(
                                                                                                                     input->sst4coeffile);
                                                                                                         }},
                                                                                                 {std::string("SST3"),
                                                                                                         [](const instr *input) {
                                                                                                             return std::string(
                                                                                                                     input->sst3coeffile);
                                                                                                         }}};

    /**
     * @brief returns the SSES LUTs paths to netCDF files depending on requested product
     */
    const std::unordered_map<std::string, std::string (*)(const instr *)> get_sses_coeffs_path = {{std::string("SST"),
                                                                                                          [](const instr *input) {
                                                                                                              return std::string(
                                                                                                                      input->sstssesfile);
                                                                                                          }},
                                                                                                  {std::string("SST4"),
                                                                                                          [](const instr *input) {
                                                                                                              return std::string(
                                                                                                                      input->sst4ssesfile);
                                                                                                          }},
                                                                                                  {std::string("SST3"),
                                                                                                          [](const instr *input) {
                                                                                                              return std::string(
                                                                                                                      input->sst3ssesfile);
                                                                                                          }}};

}
/**
 * @brief Contains supportive functions and constants.
 *
 */
namespace cldmsk {
    const float cldthresh = -1.0;
    const std::unordered_map<std::string, float> cldthresh_list = {{"modis", 0.01},
                                                                   {"viirs", 0.04}}; // threshold for cold SST test
    const std::unordered_map<std::string, const std::unordered_map<std::string, int>> bands_set =     // bands WV for VIIRS/MODIS needed to produce SST/Cloud mask
            {{"modis",
                     {{"ibred", 678},
                             {"ib07", 748},
                             {"ib16", 1640},
                             {"ib37", 3750},
                             {"ib39", 3959},
                             {"ib40", 4050},
                             {"ib67", 6715},
                             {"ib73", 7325},
                             {"ib85", 8550},
                             {"ib11", 11000},
                             {"ib12", 12000}}},
             {"viirs",
                     {{"ibred", 672},
                             {"ib07", 748},
                             {"ib16", 1601},
                             {"ib37", 3750},
                             {"ib40", 4050},
                             {"ib85", 8550},
                             {"ib11", 11000},
                             {"ib12", 12000}}}};
    const std::unordered_map<std::string, size_t> bt_box_sizes = {{"modis", 3},
                                                                  {"viirs", 5}}; // bt box size
    const float hisenz = 55.0;                                                                 // senz thresholds
    const float vhisenz = 75.0;
    const float vhisenzv2 = 65.0;
    const float solznight = 85.0;
    const float Btmin = -4.0; // BT ranges for BT range cloud mask tests
    const float Btmax = 37.0;
    const float Btmax40 = 35.0;
    const float SSTmin = -1.8; // SST ranges for SST cloud mask tests
    const float SSTmax = 40.0;
    const float SSTmaxn = 37.0;
    const float glintmax = 0.005;
    const float dBtmin = 0.0; // dBT ranges for the cloud mask tests
    const float dBtmax = 3.6;
    const float dBt4min = 0.0;
    const float dBt4max = 8.0;
    const float SSTdiff = 3.0;
    const float SSTvdiff = 5.0; // dSST ranges (thresholds) for the cloud mask tests
    const float SST4diff1 = -0.8;
    const float SST4diff2 = -1.0;
    const float SST3diff1 = 0.8;
    const float SST3diff2 = 1.0;
    const float Bt11unif1 = 0.7; // BT thresholds for BT uniformity tests
    const float Bt12unif1 = 0.7;
    const float Bt11unif2 = 1.9;
    const float Bt12unif2 = 1.9;
    const float Bt37unif1 = 0.7;
    const float Bt37unif2 = 1.9;
    const float Bt39unif1 = 0.7;
    const float Bt40unif1 = 0.7;
    const float Bt39unif2 = 1.9;
    const float Bt40unif2 = 1.9;
    const float dBtrefmin = -1.1;
    const float dBtrefmax = 10.0;
    const float equatorialNorth = 30.0; // the ranges where SST diff tests are applied
    const float equatorialSouth = -10.0;
    const float equatorialWest = -105.0;
    const float equatorialEast = 105.0;
    const float invalid_val = (BAD_FLT + 0.1); // bad vals ranges
    const float max_bt = BT_HI - 0.1;
    const float min_bt = BT_LO + 0.1;
    const float min_lt = BAD_FLT + 1.0;
    const float max_lt = BT_HI - 1.0;
    const std::unordered_set<int> modis_sensors = {MODIST,
                                                   MODISA}; // sensors/platforms INT values and corresponding string values
    const std::unordered_set<int> viirs_sensors = {VIIRSJ1, VIIRSJ2, VIIRSN};
    const std::unordered_map<std::string, std::unordered_set<int>> all_sensors = {{"modis", modis_sensors},
                                                                                  {"viirs", viirs_sensors}};
    const std::unordered_map<int, std::string> platforms = {{MODIST,  "terra"},
                                                            {MODISA,  "aqua"},
                                                            {VIIRSN,  "npp"},
                                                            {VIIRSJ1, "j1"},
                                                            {VIIRSJ2, "j2"}};
    const double scan_time_modis_t_day1 = 972777600; // 29 Oct 2000, MODIS TERRA correction dates
    const double scan_time_modis_t_day2 = 993945600; //  1 Jul 2001 
    const float el_corr_modis_t_1 = 0.4452; // MODIS TERRA correction values
    const float el_corr_modis_t_2 = 0.2593;

    /**
     * @brief month variable for the desicion tree traversal
     * @return
     */
    std::vector<float> & month_data();

    /**
     * @brief Get the sensor string
     * @param key - integer sensor ID
     * @return std::string
     */
    std::string get_sensor(int key);

    /**
     * @brief
     * prints a vector
     * @tparam T
     * @param stream - out stream
     * @param data - vector
     * @return std::ostream&
     */
    template<class T>
    std::ostream &operator<<(std::ostream &stream, const std::vector<T> &data) {
        for (const auto &el: data) {
            stream << el << " ";
        }
        stream << "\n";
        return stream;
    }

    /**
     * @brief Legacy BT diff test with interpolation adjutment.
     *
     * @param ip - pixel
     * @param BT39 - Bt39
     * @param BT40 - Bt40
     * @param l1rec - l1 record
     * @param fullscanpix - full scan size
     * @return float
     */
    float btrefdiffv6(int32_t ip, float BT39, float BT40, const l1str *l1rec, int fullscanpix);

    enum product_types {
        SST,
        SST3,
        SST4
    };

    // Functions ptr below get the values of BT/RHOs. Compute a value of BT/RHO or its mask at a given pixel
    typedef bool (*get_valid)(const l1str &, int, int, int); //
    typedef float (*get_value)(const l1str &, int, int, int);

    const get_valid cirrus_mask = [](const l1str &l1str_, int pixel, int nbands, int ib) {
        return l1str_.rho_cirrus[pixel] > invalid_val;
    };
    const get_value cirrus_value = [](const l1str &l1str_, int pixel, int nbands,
                                int ib) { return l1str_.rho_cirrus[pixel]; };
    // bts
    const get_valid bt_mask = [](const l1str &l1str_, int pixel, int nbands, int ib) {
        const float bt = l1str_.Bt[pixel * nbands + ib];
        return bt > min_bt && bt < max_bt;
    };
    const get_value bt_value = [](const l1str &l1str_, int pixel, int nbands, int ib) {
        return l1str_.Bt[pixel * nbands + ib];
    };
    // rho (RSMAS)
    const get_valid rho_mask = [](const l1str &l1str_, int pixel, int nbands, int ib) {
        const float lt = l1str_.Lt[pixel * nbands + ib];
        return lt > min_lt && lt < max_lt;
    };
    const get_value rho_value = [](const l1str &l1str_, int pixel, int nbands, int ib) {
        const float lt = l1str_.Lt[pixel * nbands + ib];
        const float fo = l1str_.Fo[ib];
        const float csolz = l1str_.csolz[pixel];
        const float rho = OEL_PI * lt / fo / csolz;
        return rho;
    };

    const get_value cldrh_value = [](const l1str &l1str_, int pixel, int nbands, int ib) {
        const float lt = l1str_.Lt[pixel * nbands + ib];
        const float fo = l1str_.Fo[ib];
        // const float csolz = l1str_.csolz[pixel];
        const float tg_sol = l1str_.tg_sol[pixel * nbands + ib];
        const float tg_sen = l1str_.tg_sen[pixel * nbands + ib];
        const float t_sen = l1str_.t_sen[pixel * nbands + ib];
        const float t_sol = l1str_.t_sol[pixel * nbands + ib];
        const float cldrh = OEL_PI * lt / fo / tg_sol / tg_sen / t_sol / t_sen; // /csolz
        return cldrh;
    };
    const get_valid cldrh_mask = [](const l1str &l1str_, int pixel, int nbands, int ib) {
        const float lt = l1str_.Lt[pixel * nbands + ib];
        return lt > 0.0 && lt < max_lt;
    };
    // other vars (not bandd)
    const std::unordered_map<std::string, float *(*)(const l1str &)> get_non_BT_vars =
            {
                    {"solz",      [](const l1str &l1rec) { return l1rec.solz; }},
                    {"senz",      [](const l1str &l1rec) { return l1rec.senz; }},
                    {"wv",        [](const l1str &l1rec) { return l1rec.wv; }},
                    {"glintcoef", [](const l1str &l1rec) { return l1rec.glint_coef; }},
                    {"lat",       [](const l1str &l1rec) { return l1rec.lat; }},
                    {"lon",       [](const l1str &l1rec) { return l1rec.lon; }},
                    {"month",     [](const l1str &l1rec) { return month_data().data(); }}};

    /**
     * @brief Get the mask (line) of a BT
     *
     * @param mask - mask
     * @param l1str - l1 record
     * @param npix - number of pixel in a line
     * @param nbands - number of bands
     * @param ib - band number
     * @param get_mask - mask function ptr
     */
    void get_var_mask(int *mask, const l1str &l1str, size_t npix, int nbands, int ib, get_valid get_mask);

    /**
     * @brief the  same as get_var_mask but gets values
     */
    void get_var_vals(float *BT, const l1str &l1str, size_t npix, int nbands, int ib, get_value get_val);

    /**
     * @brief Compute max values within a 1D window in one line in a queue
     *
     * @param maxs_1d - output queue
     * @param inp1d - inputd queu
     * @param mask1d - input mask
     * @param npix - number of pixels
     * @param radius - radius of the sliding window
     */
    void get_window_1D_max(float *maxs_1d, const float *inp1d, const int *mask1d, size_t npix, size_t radius);

    /**
     * @brief Get the max value of center line with 2D sliding window using a 1D max queue computed with get_window_1D_max
     *
     * @param max_global - max values, out
     * @param inp2d - input queue, obtained from get_window_1D_max
     * @param npix - number of pixels
     * @param nscan - queue size/ nubmer of scans
     * @param center - center of the line
     * @param radius - window size
     */
    void
    get_total_2d_max(float *max_global, const float *inp2d, size_t npix, size_t nscan, size_t center, size_t radius);

    /**
     * @brief the same as get_window_1D_max but for min
     */
    void get_window_1D_min(float *mins_1d, const float *inp1d, const int *mask1d, size_t npix, size_t radius);

    /**
     * @brief the same as get_total_2d_max but for min
     */
    void
    get_total_2d_min(float *min_global, const float *inp2d, size_t npix, size_t nscan, size_t center, size_t radius);

    /**
     * @brief Compute standard deviation in a window
     *
     * @param inp2d - input array
     * @param mask2d - mask
     * @param npix - pixel in a line (x dim)
     * @param nscan - queue size ( y dim)
     * @param radius_x - radius along npix
     * @param radius_y - radiis along nscan
     * @param out - STD
     * @param center - Center of the queu
     * @param l1qrec - queue struct
     */
    void get_std_box(const float *inp2d, const int *mask2d, size_t npix, size_t nscan, size_t radius_x, size_t radius_y,
                     float *out, size_t center, l1qstr *l1qrec);

    /**
     * @brief
     * SSES data structure. Contains the LUTs and a look up routine. Reads NetCDF file LUTs file
     *
     */
    struct SSESData {
        size_t nsst;
        size_t nday;
        size_t nquar;
        size_t nsenz;
        size_t ndiff;
        size_t nlat;
        size_t nqual;
        size_t total_size;
        size_t npix;
        std::vector<float> sst;
        std::vector<float> senz;
        std::vector<float> diff;
        std::vector<float> lat;
        std::vector<float> qual;
        // LUTS
        std::vector<float> bias;
        std::vector<float> stdv;
        std::vector<float> bias_mean;
        std::vector<int16_t> counts;

        std::vector<float> bias_out;
        std::vector<float> stdv_out;
        std::vector<float> bias_mean_out;
        std::vector<int16_t> counts_out;

        /**
         * @brief Construct a new SSESData object
         *
         * @param nc_path_LUTs - path to netcdf file
         * @param npix - number of pixel in the line
         */
        SSESData(const std::string &nc_path_LUTs, size_t npix) : npix(npix) {
            try {
                netCDF::NcFile sst_sses_lut_bias(nc_path_LUTs, netCDF::NcFile::read);
                {
                    std::cout << "Reading SSES LUT netCDF file : " << nc_path_LUTs << std::endl;
                }
                const auto vars = sst_sses_lut_bias.getVars();
                for (const auto &var: vars) {
                    netCDF::NcVar var_sses = sst_sses_lut_bias.getVar(var.first);
                    auto dims = var_sses.getDims();
                    for (const auto &dim: dims) {
                        if (dim.getName() == "nlat") {
                            nlat = dim.getSize();
                        }
                        if (dim.getName() == "nsatz") {
                            nsenz = dim.getSize();
                        }
                        if (dim.getName() == "ndiff") {
                            ndiff = dim.getSize();
                        }
                        if (dim.getName() == "nquar") {
                            nquar = dim.getSize();
                        }
                        if (dim.getName() == "nqual") {
                            nqual = dim.getSize();
                        }
                        if (dim.getName() == "nsst") {
                            nsst = dim.getSize();
                        }
                        if (dim.getName() == "nday") {
                            nday = dim.getSize();
                        }
                    }
                }
                {
                    total_size = nday * nsst * nquar * nqual * nlat * ndiff * nsenz;
                }
                bias.resize(total_size);
                stdv.resize(total_size);

                sst.resize(nsst);
                senz.resize(nsenz);
                diff.resize(ndiff);
                lat.resize(nlat);
                qual.resize(nqual);
                {
                    netCDF::NcVar var_sses = sst_sses_lut_bias.getVar("bias");
                    var_sses.getVar(bias.data());
                }
                {
                    netCDF::NcVar var_sses = sst_sses_lut_bias.getVar("bias_mean");
                    if (!var_sses.isNull()) {
                        bias_mean.resize(total_size);
                        var_sses.getVar(bias_mean.data());
                    }
                }
                {
                    netCDF::NcVar var_sses = sst_sses_lut_bias.getVar("stdv");
                    var_sses.getVar(stdv.data());
                }
                {
                    netCDF::NcVar var_sses = sst_sses_lut_bias.getVar("counts");
                    if (!var_sses.isNull()) {
                        counts.resize(total_size);
                        var_sses.getVar(counts.data());
                    }
                }
                {
                    netCDF::NcVar var_sses = sst_sses_lut_bias.getVar("sst");
                    var_sses.getVar(sst.data());
                }
                {
                    netCDF::NcVar var_sses = sst_sses_lut_bias.getVar("senz");
                    var_sses.getVar(senz.data());
                }
                {
                    netCDF::NcVar var_sses = sst_sses_lut_bias.getVar("lat");
                    var_sses.getVar(lat.data());
                }
                {
                    netCDF::NcVar var_sses = sst_sses_lut_bias.getVar("qual");
                    var_sses.getVar(qual.data());
                }
                {
                    netCDF::NcVar var_sses = sst_sses_lut_bias.getVar("BTdiff");
                    var_sses.getVar(diff.data());
                }
                {
                    bias_out = std::vector<float>(npix, BAD_FLT);
                    stdv_out = std::vector<float>(npix, BAD_FLT);

                    if (!counts.empty())
                        counts_out = std::vector<int16_t>(npix, 0);
                    if (!bias_mean.empty())
                        bias_mean_out = std::vector<float>(npix, BAD_FLT);
                }
            }
            catch (...) {
                std::cerr << "Error opening and reading SSES LUT file : " << nc_path_LUTs << std::endl;
                std::cerr << "Exiting ... " << nc_path_LUTs << std::endl;
                exit(EXIT_FAILURE);
            }
        };

        size_t
        operator()(size_t iqual, size_t ilat, size_t idiff, size_t isenz, size_t iquar, size_t iday, size_t isst) {
            size_t index = isst + iday * nsst + iquar * nday * nsst + isenz * nquar * nsst * nday +
                           idiff * nday * nsst * nquar * nsenz + ilat * nday * nsst * nquar * nsenz * ndiff +
                           iqual * nday * nsst * nquar * nsenz * nlat * ndiff;
            return index;
        }

        /**
         * @brief Calculates the SSES bias vars
         *
         * @param diff - input BT diff
         * @param sst - input SST
         * @param solz - input solar zenith angle
         * @param senz - input sattelite zenith angle
         * @param day - input night or day
         * @param lat - input latitude
         * @param qual - qual level, input
         * @param flags - cloud mask flag, input
         * @param glint - input, glint
         */
        void
        calculate_sses_bias_stats(const float *diff, const float *sst, const float *solz, const float *senz, size_t day,
                                  const float *lat, const int8_t *qual, const int16_t *flags, const float *glint) {
            if (!bias_mean_out.empty())
                std::fill(bias_mean_out.begin(), bias_mean_out.end(), BAD_FLT);
            std::fill(stdv_out.begin(), stdv_out.end(), BAD_FLT);
            std::fill(bias_out.begin(), bias_out.end(), BAD_FLT);
            if (!counts_out.empty())
                std::fill(counts_out.begin(), counts_out.end(), 0);
            size_t iquar, iday, ilat, isenz, isst, idiff;
            size_t i;
            for (size_t i_p = 0; i_p < npix; i_p++) {
                size_t iqual = *(qual + i_p);
                if (iqual >= nqual)
                    continue;
                if (iqual == 1 && (glint[i_p] > glintmax) && ((flags[i_p] & SSTF_BTNONUNIF) == 0) &&
                    ((flags[i_p] & SSTF_HISENZ) == 0)) {
                    iqual = 0;
                }
                iquar = day / 91;
                iday = 0;
                if (solz[i_p] < solznight) {
                    iday = 1;
                }
                if (iday > nday - 1)
                    continue;
                for (i = 1; i < nlat; i++)
                    if (*(lat + i_p) < this->lat.at(i))
                        break;
                ilat = i - 1;
                for (i = 1; i < nsenz; i++)
                    if (*(senz + i_p) < this->senz.at(i))
                        break;
                isenz = i - 1;
                for (i = 1; i < nsst; i++)
                    if (*(sst + i_p) < this->sst.at(i))
                        break;
                isst = i - 1;
                for (i = 1; i < ndiff; i++)
                    if (*(diff + i_p) <= this->diff.at(i))
                        break;
                idiff = i - 1;
                const size_t index = this->operator()(iqual, ilat, idiff, isenz, iquar, iday, isst);
                bias_out.at(i_p) = bias.at(index);
                if (!bias_mean_out.empty())
                    bias_mean_out.at(i_p) = bias_mean.at(index);
                stdv_out.at(i_p) = stdv.at(index);
                if (!counts_out.empty())
                    counts_out.at(i_p) = counts.at(index);
            }
        }
    };

    void read_sst_bands(const std::string &sat, std::unordered_map<std::string, int> &bands);
}

/**
 * @brief Contains the data/stats for variables (BTs and SSTs) needed for cloud masks tests / ADT
 *
 */
namespace bstats {
    /**
     * @brief A base class for a pair of vars.
     * Holds the variables difference and ratio.
     */
    struct PairOfVars {
        std::vector<float> difference;
        std::vector<float> ratio;
        size_t npix;
        int current_scan = 0;

        PairOfVars() {
        }

        PairOfVars(size_t npix) : npix(npix) {
        }

        float *operator()() {
            return difference.data();
        }

        float *get_Tdeflong() {
            return ratio.data();
        }
    };

    /**
     * @brief A derived class for a pair of BTs
     */
    struct PairOfVarsBT : public PairOfVars {
        int ib1, ib2;
        float offset = 0.00001; // in case if BT_2 == 0
        PairOfVarsBT() {

        };

        PairOfVarsBT(size_t npix, int ib1, int ib2, l1str &l1str, bool has_difference = false, bool has_ratio = false)
                : PairOfVars(npix), ib1(ib1), ib2(ib2) {
            if (has_difference) {
                difference = std::vector<float>(npix, BAD_FLT);
                get_difference(l1str);
            }
            if (has_ratio) {
                ratio = std::vector<float>(npix, BAD_FLT);
                get_ratio(l1str);
            }
        }

        void get_difference(l1str &l1str) {
            for (size_t i_p = 0; i_p < npix; i_p++) {
                difference.at(i_p) =
                        cldmsk::bt_value(l1str, i_p, NBANDSIR, ib1) - cldmsk::bt_value(l1str, i_p, NBANDSIR, ib2);
            }
            if (!ratio.empty())
                get_ratio(l1str);
        }

        void get_ratio(l1str &l1str) {
            if (ratio.empty())
                ratio = std::vector<float>(npix, BAD_FLT);
            for (size_t i_p = 0; i_p < npix; i_p++) {
                if (cldmsk::bt_value(l1str, i_p, NBANDSIR, ib1) != 0)
                    ratio.at(i_p) = 1.0 - cldmsk::bt_value(l1str, i_p, NBANDSIR, ib2) /
                                          cldmsk::bt_value(l1str, i_p, NBANDSIR, ib1);
                else if (cldmsk::bt_value(l1str, i_p, NBANDSIR, ib2) == 0) {
                    ratio.at(i_p) = 0.0;
                } else {
                    ratio.at(i_p) = 1.0 - cldmsk::bt_value(l1str, i_p, NBANDSIR, ib2) / offset;
                }
            }
        }
    };

    /**
     * @brief A derived class for a pair of SSTs
     */
    struct PairOfSST : public PairOfVars {
        float *sst1;
        float *sst2;

        PairOfSST() {
        }

        PairOfSST(size_t npix, float *sst1, float *sst2, bool has_difference) : PairOfVars(npix), sst1(sst1),
                                                                                sst2(sst2) {
            if (has_difference) {
                difference = std::vector<float>(npix, BAD_FLT);
                get_difference();
            }
        }

        void get_difference() {
            std::transform(sst1, sst1 + npix, sst2, difference.begin(), std::minus<float>());
        }
    };

    /**
     * @brief A base class for a variable.
     * Holds values for STD, MAX, MIN within a sliding windows
     */
    struct Stats {
        int current_scan = 0;
        size_t npix;
        size_t nscan;
        size_t rad_x, rad_y, center;
        size_t i_s, i_e;
        l1qstr *l1qrec;
        std::vector<float> var_box;
        std::vector<int> mask_box;
        std::vector<float> var_max, var_min, var_std, var_minmax;
        std::vector<float> var_max_box, var_min_box;
        std::unordered_map<std::string, float *> return_vals;

        /**
         * @brief Init and compute STD
         */
        void set_std() {
            if (var_std.empty()) {
                var_std = std::vector<float>(npix, BAD_FLT);
                cldmsk::get_std_box(var_box.data(), mask_box.data(), npix, nscan, rad_x, rad_y, var_std.data(), center,
                                    l1qrec);
                return_vals["STD"] = var_std.data();
            }
        }

        /**
         * @brief Init and compute MIN
         */
        void set_min() {
            if (var_min.empty()) {
                const float max_possible_val = std::abs(BAD_FLT);
                var_min = std::vector<float>(npix, max_possible_val);
                var_min_box = std::vector<float>(npix * nscan, max_possible_val);
                return_vals["MIN"] = var_min.data();
                for (size_t i_q = i_s; i_q <= i_e; i_q++) {
                    const size_t index = npix * i_q;
                    cldmsk::get_window_1D_min(var_min_box.data() + index, var_box.data() + index,
                                              mask_box.data() + index, npix, rad_x);
                }
                cldmsk::get_total_2d_min(var_min.data(), var_min_box.data(), npix, nscan, center, rad_y);
            }
        }

        /**
         * @brief Init and compute MAX
         */
        void set_max() {
            if (var_max.empty()) {
                var_max = std::vector<float>(npix, BAD_FLT);
                var_max_box = std::vector<float>(npix * nscan, BAD_FLT);
                return_vals["MAX"] = var_max.data();
                for (size_t i_q = i_s; i_q <= i_e; i_q++) {
                    const size_t index = npix * i_q;
                    cldmsk::get_window_1D_max(var_max_box.data() + index, var_box.data() + index,
                                              mask_box.data() + index, npix, rad_x);
                }
                cldmsk::get_total_2d_max(var_max.data(), var_max_box.data(), npix, nscan, center, rad_y);
            }
        }

        /**
         * @brief Init and compute difference between MAX and MIN
         */
        void set_maxmin() {
            if (var_minmax.empty()) {
                var_minmax = std::vector<float>(npix, BAD_FLT);
                return_vals["MAXMIN"] = var_minmax.data();
                set_max();
                set_min();
                std::transform(var_max.begin(), var_max.end(), var_min.begin(), var_minmax.begin(),
                               std::minus<float>());
            }
        }

        /**
         * @brief when a l1 queue is updated, we need to updated the supporting arrays/queus as well
         *
         */
        void rearrange_stats() {
            if (!var_std.empty()) {
                cldmsk::get_std_box(var_box.data(), mask_box.data(), npix, nscan, rad_x, rad_y, var_std.data(), center,
                                    l1qrec);
            };
            if (!var_max.empty()) {
                for (size_t i_q = i_s; i_q < i_e; i_q++) {
                    const size_t index = (i_q + 1) * npix;
                    std::copy(var_max_box.begin() + index, var_max_box.begin() + index + npix,
                              var_max_box.begin() + index - npix);
                }
                const size_t index = npix * i_e;
                cldmsk::get_window_1D_max(var_max_box.data() + index, var_box.data() + index, mask_box.data() + index,
                                          npix, rad_x);
                cldmsk::get_total_2d_max(var_max.data(), var_max_box.data(), npix, nscan, center, rad_y);
            };
            if (!var_min.empty()) {
                for (size_t i_q = i_s; i_q < i_e; i_q++) {
                    const size_t index = (i_q + 1) * npix;
                    std::copy(var_min_box.begin() + index, var_min_box.begin() + index + npix,
                              var_min_box.begin() + index - npix);
                }
                const size_t index = npix * i_e;
                cldmsk::get_window_1D_min(var_min_box.data() + index, var_box.data() + index, mask_box.data() + index,
                                          npix, rad_x);
                cldmsk::get_total_2d_min(var_min.data(), var_min_box.data(), npix, nscan, center, rad_y);
            };
            if (!var_minmax.empty()) {
                std::transform(var_max.begin(), var_max.end(), var_min.begin(), var_minmax.begin(),
                               std::minus<float>());
            };
        };

        Stats() {};

        Stats(size_t npix, size_t nscan, size_t rad_x, size_t rad_y, size_t center, l1qstr *l1qrec) : npix(npix),
                                                                                                      nscan(nscan),
                                                                                                      rad_x(rad_x),
                                                                                                      rad_y(rad_y),
                                                                                                      center(center),
                                                                                                      l1qrec(l1qrec) {
            i_s = std::min(std::max(0, (int) center - (int) rad_y), (int) nscan - 1);
            i_e = std::max(std::min(nscan - 1, center + rad_y), 0ul);
        }

        float *operator()(const std::string &key) {
            return return_vals.at(key);
        }

        int *get_valid_mask() {
            return mask_box.data() + npix * center;
        }
    };

    /**
     * @brief A derived class to compute stats for an SST
     */
    struct StatsSST : public Stats {
        float *sst;

        StatsSST() {};

        StatsSST(size_t npix, size_t nscan, size_t rad_x, size_t rad_y, size_t center, l1qstr *l1qrec, float *sst,
                 bool std_requested = false, bool max_req = false, bool min_req = false, bool min_max_req = false)
                : Stats(npix, nscan, rad_x, rad_y, center, l1qrec), sst(sst) {
            return_vals["VAR"] = sst;
            if (std_requested || max_req || min_req || min_max_req) {
                var_box = std::vector<float>(npix * nscan, BAD_FLT);
                mask_box = std::vector<int>(npix * nscan, 0);
                for (size_t index = 0; index < npix * nscan; index++) {
                    var_box.at(index) = *(sst - npix * center + index);
                    mask_box.at(index) = (var_box.at(index) >= cldmsk::SSTmin) && (var_box.at(index) <= cldmsk::SSTmax);
                }
                if (std_requested)
                    set_std();
                if (min_max_req)
                    set_maxmin();
                else if (min_req)
                    set_min();
                else if (max_req)
                    set_max();
            }
        }

        void rearrange() {
            if (!var_box.empty()) {
                for (size_t index = 0; index < npix * nscan; index++) {
                    var_box.at(index) = *(sst - npix * center + index);
                    mask_box.at(index) = var_box.at(index) >= cldmsk::SSTmin && var_box.at(index) <= cldmsk::SSTmax;
                }
            }
            rearrange_stats();
        }
    };

    /**
     * @brief A derived class to compute stats for a BT
     */
    struct StatsVarBand : public Stats {
        int ib;
        size_t nbands;
        cldmsk::get_valid get_mask;
        cldmsk::get_value get_val;

        StatsVarBand() {};

        StatsVarBand(size_t npix, size_t nscan, size_t rad_x, size_t rad_y, size_t center, l1qstr *l1qrec, int ib,
                     size_t nbands, cldmsk::get_valid get_mask, cldmsk::get_value get_val, bool std_requested = false,
                     bool max_requested = false, bool min_requested = false, bool min_max_requested = false) : Stats(
                npix, nscan, rad_x, rad_y, center, l1qrec), ib(ib), nbands(nbands), get_mask(get_mask),
                                                                                                               get_val(get_val) {
            var_box = std::vector<float>(npix * nscan, BAD_FLT);
            mask_box = std::vector<int>(npix * nscan, 0);
            for (size_t i_q = i_s; i_q <= i_e; i_q++) {
                const l1str &q_ref = l1qrec->r[i_q];
                cldmsk::get_var_vals(var_box.data() + npix * i_q, q_ref, npix, nbands, ib, get_val);
                cldmsk::get_var_mask(mask_box.data() + npix * i_q, q_ref, npix, nbands, ib, get_mask);
            }
            return_vals["VAR"] = var_box.data() + npix * center;
            if (std_requested) {
                set_std();
            }
            if (min_max_requested) {
                max_requested = true;
                min_requested = true;
                set_maxmin();
            }
            if (max_requested) {
                set_max();
            }
            if (min_requested) {
                set_min();
            }
        }

        void rearrange() {
            for (size_t i_q = i_s; i_q < i_e; i_q++) {
                const size_t index = (i_q + 1) * npix;
                std::copy(var_box.begin() + index, var_box.begin() + index + npix, var_box.begin() + index - npix);
                std::copy(mask_box.begin() + index, mask_box.begin() + index + npix, mask_box.begin() + index - npix);
            }
            {
                const l1str &q_ref = l1qrec->r[i_e];
                cldmsk::get_var_vals(var_box.data() + npix * i_e, q_ref, npix, nbands, ib, get_val);
                cldmsk::get_var_mask(mask_box.data() + npix * i_e, q_ref, npix, nbands, ib, get_mask);
            }
            rearrange_stats();
        }
    };
}
#endif
