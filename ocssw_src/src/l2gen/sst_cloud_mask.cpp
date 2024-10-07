#include "sst_adt.hpp"
#include "sst_cloud_mask_utils.hpp"
#include "sst.h"
#include "timeutils.h"
/**
 * @brief
 * anonym namespace to define some parameters/local fucntions
 */
namespace {
    const double CtoK = 273.15e0;  // Celcius to Kelvin
/**
 * @brief
 * Prints a vector into a stream
 * @tparam T
 * @param stream stream
 * @param data vector
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
 * @brief
 * parsing the cloud mask test name variable; determines which BT/RHO/SST needs
 * to be used in order to compute the variable used for both ADT and quality
 * flags settings
 * @param test_name - the test name. For example, d_SST_SST4 will be interpreted
 * as SST - SST4 and RHOred_MAXMIN as rho_red_max - rho_red_min
 * @return std::vector<std::string> - strings of keywords
 */
    std::vector<std::string> parse_test_name(const std::string &test_name) {
        std::vector<std::string> out;
        boost::split(out, test_name, boost::is_any_of("_"));
        return out;
    }

/**
 * @brief
 * find the interval where a given point belongs to; assumed that data is sorted
 * @tparam T
 * @param data - the interval
 * @param val - the point
 * @param start -  the min bound search
 * @param end - the max bound searcg
 */
    template<class T>
    void binary_search(const std::vector<T> &data, const T &val, size_t &start,
                       size_t &end) {
        while (end - start > 1) {
            size_t middle = (end + start) / 2;
            if (data.at(middle) >= val) {
                end = middle;
                binary_search(data, val, start, end);
            } else {
                start = middle;
                binary_search(data, val, start, end);
            }
        }
    }

// bands
    std::unordered_map<std::string, int> ibands;
    size_t n_bands;
// variables needed for quality flags setting
    const boost::unordered_map<std::pair<std::string, std::string>,
            std::set<std::string>>
            needed_stats_for_cloud_mask = {
            {{"viirs", "SST3"},
                    {"BT37_MAXMIN", "BT12_MAXMIN", "BT11_MAXMIN", "BT40", "d_BT37_BT12",
                            "d_BT11_BT12"}},
            {{"viirs", "SST"},
                    {"BT12_MAXMIN", "BT11_MAXMIN", "BT40",        "CLDRHred_MAXMIN"}},
            {{"modis", "SST"},
                    {"BT12_MAXMIN", "BT11_MAXMIN", "BT40_MAXMIN", "BT39_MAXMIN",
                                                                          "CLDRHred_MAXMIN"}},
            {{"modis", "SST4"},
                    {"BT12_MAXMIN", "BT11_MAXMIN", "BT40_MAXMIN", "BT39_MAXMIN",
                                                                          "d_BT39_BT40"}}};

}  // namespace

/**
 * @brief *
 * Class to compute the SST products (the sst, cloud mask etc)
 * @param this_product - product type. Can be SST, SST3 or SST4
 * @param sst - contatins the SST queue
 * @param valid_mask - ocean/vlaid bt mask
 * @param flags_sst - bit flags SST
 * @param qual_sst - quality flags (0-5)
 * @param treesum - the treesum array
 * @param scan_time - scan time
 * @param npix - number of pixels in line
 * @param n_q_size - queu size
 * @param fullscanpix - number of pixels for full scan
 * @param btbox - size of the sliding windwow
 * @param latwin - the width of latitude window used in interpolation for the
 * NLSST
 * @param processed_scan_sst - processed scan sst @param  processed_scan_flags -
 * processed scan for flags, @param  processed_scan_sses - processed scan for
 * sses
 * @param regression_model_name - regression model used @param nlatbands - size
 * of the latitude arrays used in NLSST @param ncoeffs - number of regression
 * coefficients
 * @param bounds_search - bounds for latitude used in NLSST
 * @param sses_bias_data - the object which contains SSES data
 * @param desicion_tree - ADT
 * @param test_data - contains pointers to the data used for in ADT/Cloud
 * masking/SSES Bias calculations
 * @param test_masks - contains pointers to the data(masks) used for in
 * ADT/Cloud masking
 * @param StatsBTs - contains the data/stats needed for ADT/Cloudmasking/SSES
 * Bias calculations. Each map element Stores the data for a Single BT
 * @param PairOfVarsBTs - containts the ratio/ difference of two BT
 * @param StatsSSTs - contains the data/stats for single SST product
 * @param PairOfSSTs - containts the difference of two SST
 * @param modis_ref - the modis reference. if product SST4, than it uses the
 * CMC. if the product SST, at night it uses SST4 as a reference
 * @param dust_correction_data - dust correction for MODIST SST
 * @param external_products - Calculations of SST product (the sst and the cloud
 * mask) may require some variables from the SST3 and SST4 (and vice versa).
 * This variable stores the references to the products
 */
class SeaSurfaceTemperatureCalcuations {
private:
    std::string this_product;
    std::vector<float> sst;
    std::vector<int> valid_mask;
    std::vector<int16_t> flags_sst;
    std::vector<int8_t> qual_sst;
    std::vector<float> treesum;
    double scan_time;
    int16_t day;
    size_t npix;
    size_t n_q_size;
    size_t fullscanpix = 1354;
    const float latwin = 2.5f;
    size_t bt_box = 5;  // default is 5
    size_t i_center, i_s, i_e;
    int32_t processed_scan_sst = -1;
    int32_t processed_scan_flags = -1;
    int32_t processed_scan_sses = -1;
    std::string sensor;
    std::string platform;
    std::vector<float> sst_regression_coefficients;
    std::vector<float> bounds_lat, bounds_search;
    std::unique_ptr<cldmsk::SSESData> sses_bias_data;
    // regression model name
    std::string regression_model;
    size_t nlatbands, ncoeffs;
    int ib11, ib12, ib37, ib40, ib39, ib85;
    std::unordered_map<std::string, SeaSurfaceTemperatureCalcuations *>
            external_products;
    adt::Treenode *desicion_tree = nullptr;
    std::set<std::string> test_names;
    std::unordered_map<std::string, std::vector<std::string>>
            test_classifications;
    std::unordered_map<std::string, float *> test_data;
    std::unordered_map<std::string, int *> test_masks;
    // maps
    boost::unordered_map<std::string, bstats::StatsVarBand> StatsBTs;
    boost::unordered_map<std::pair<std::string, std::string>,
            bstats::PairOfVarsBT>
            PairOfVarsBTs;
    boost::unordered_map<std::string, bstats::StatsSST> StatsSSTs;
    boost::unordered_map<std::pair<std::string, std::string>, bstats::PairOfSST>
            PairOfSSTs;
    std::string path_to_cloud_mask_config;
    // modis reference
    std::vector<float> modis_ref;
    // modis dust correction
    std::vector<float> dust_correction_data;

    // modis A BT40 data for detector #0
    std::vector<float> modis_a_bt40;
    l1qstr *l1qstr_modis_a = nullptr;
    // modis T electrionic correction
    float elecor_terra = 0.0f;
    bool is_aqua = false;
    bool is_terra = false;

    /**
     * @brief Set the correction to the modis 4 micron BT for detector #0
     * @param l1rec - l1 record
     * @param l1qrec - l1 queue
     * @param input - input struct
     * @param product_type - product type
     */
    void set_modis_a_bt_40_data(l1str *l1rec, l1qstr *l1qrec, instr *input,
                                const std::string &product_type) {
        if (!is_aqua) return;
        std::vector<float> &var_box = StatsBTs.at("40").var_box;
        const size_t center = StatsBTs.at("40").center;
        if (l1rec->detnum == 0 && l1qrec->r[i_center].detnum == 0 && i_center > 0) {
            modis_a_bt40 = std::vector<float>(npix, BAD_FLT);
            const size_t index_prev = (i_center - 1) * npix;
            const size_t index_next = (i_center + 1) * npix;
            for (size_t i_p = 0; i_p < npix; i_p++) {
                float sum = 0.0f;
                int count = 0;
                const float prev_val = var_box.at(index_prev + i_p);
                const float next_val = var_box.at(index_next + i_p);
                const int mask_prev = cldmsk::bt_mask(l1qstr_modis_a->r[i_center - 1], i_p, NBANDSIR, ib40);
                const int mask_next = cldmsk::bt_mask(l1qstr_modis_a->r[i_center + 1], i_p, NBANDSIR, ib40);
                if (mask_prev) {
                    count++;
                    sum += prev_val;
                }
                if (mask_next) {
                    count++;
                    sum += next_val;
                }
                if (count > 0) {
                    modis_a_bt40.at(i_p) = sum / count;
                } else {
                    const int mask_cur = cldmsk::bt_mask(l1qstr_modis_a->r[i_center], i_p, NBANDSIR, ib40);
                    if (mask_cur) {
                        modis_a_bt40.at(i_p) =
                            cldmsk::bt_mask(l1qstr_modis_a->r[i_center], i_p, NBANDSIR, ib40);
                    }
                }
            }
            StatsBTs.at("40").return_vals.at("VAR") = modis_a_bt40.data();
            for (auto &pair_bt : PairOfVarsBTs) {
                if (pair_bt.first.first == "40") {
                    const auto ib_pair = ibands.at("ib" + pair_bt.first.second);
                    for (size_t i_p = 0; i_p < npix; i_p++) {
                        pair_bt.second.difference.at(i_p) =
                            modis_a_bt40.at(i_p) - cldmsk::bt_value(*l1rec, i_p, NBANDSIR, ib_pair);
                    }
                }
                if (pair_bt.first.second == "40") {
                    const auto ib_pair = ibands.at("ib" + pair_bt.first.first);
                    for (size_t i_p = 0; i_p < npix; i_p++) {
                        pair_bt.second.difference.at(i_p) =
                            cldmsk::bt_value(*l1rec, i_p, NBANDSIR, ib_pair) - modis_a_bt40.at(i_p);
                    }
                }
            }
        } else {
            StatsBTs.at("40").return_vals.at("VAR") =
                    var_box.data() + npix * center;
        }
    };

    void modis_dust_correction(const l1str &q_ref, size_t i_q) {
        // get dust correction
        if (!dust_correction_data.empty() && this_product != "SST4") {
            const int sensor_id = q_ref.l1file->sensorID;
            for (size_t i_p = 0; i_p < npix; i_p++)
                if (q_ref.solz[i_p] >= cldmsk::solznight) {
                    const float dustExtinction =
                            q_ref.anc_aerosol->dust_ext[i_p];
                    float csenz = q_ref.csenz[i_p];
                    const float Bt39 =
                            cldmsk::bt_value(q_ref, i_p, NBANDSIR, ib39);
                    const float Bt85 =
                            cldmsk::bt_value(q_ref, i_p, NBANDSIR, ib85);
                    const float Bt11 =
                            cldmsk::bt_value(q_ref, i_p, NBANDSIR, ib11);
                    const float Bt12 =
                            cldmsk::bt_value(q_ref, i_p, NBANDSIR, ib12);
                    dust_correction_data.at(i_p) =
                            dust_correction(dustExtinction, csenz, Bt39, Bt85, Bt11,
                                            Bt12, sensor_id);
                    sst.at(i_p + i_q * npix) += dust_correction_data.at(i_p);
                }
        }
    }

    /**
     * @brief  Get the SST reference for the NLSST
     * @param l1rec - l1 record
     * @param l1qrec - l1 queue
     * @param input - input struct
     * @param product_type - product type
     * @param is - index in the queue
     */
    float *get_reference(l1str *l1rec, l1qstr *l1qrec, instr *input,
                         const std::string &product_type, size_t is)  //
    {
        if (sensor != "modis" || this_product == "SST4")
            return l1qrec->r[is].sstref;
        else {
            if (modis_ref.empty())
                modis_ref = std::vector<float>(npix, BAD_FLT);
            for (size_t i_p = 0; i_p < npix; i_p++) {
                const float solz = l1qrec->r[is].solz[i_p];
                if (solz >= cldmsk::solznight) {
                    modis_ref.at(i_p) = external_products.at("SST4")->get_sst(
                            l1rec, l1qrec, input,
                            "SST4")[i_p + (is - i_center) * npix];
                } else {
                    modis_ref.at(i_p) = l1qrec->r[is].sstref[i_p];
                }
            }
            return modis_ref.data();
        }
    }

    /**
     * @brief  Get the SST valid mask
     * @param q_ref - l1 record
     * @param mask - mask
     * @param sst_ref - sst_ref
     */
    void get_valid_mask_sst(const l1str &q_ref, int *mask,
                            const float *sst_ref) {
        for (size_t i_p = 0; i_p < npix; i_p++) {
            mask[i_p] = cldmsk::bt_mask(q_ref, i_p, NBANDSIR, ib11) &&
                        cldmsk::bt_mask(q_ref, i_p, NBANDSIR, ib12) &&
                        ((q_ref.flags[i_p] & LAND) == 0) &&
                        ((q_ref.flags[i_p] & NAVFAIL) == 0) &&
                        (sst_ref[i_p] >= BAD_FLT + 1.0f);
        }
    }

    /**
     * @brief Same as get_valid_mask_sst but for SST3
     */
    void get_valid_mask_sst3(const l1str &q_ref, int *mask,
                             const float *sst_ref) {
        for (size_t i_p = 0; i_p < npix; i_p++) {
            mask[i_p] = cldmsk::bt_mask(q_ref, i_p, NBANDSIR, ib37) &&
                        cldmsk::bt_mask(q_ref, i_p, NBANDSIR, ib11) &&
                        cldmsk::bt_mask(q_ref, i_p, NBANDSIR, ib12) &&
                        ((q_ref.flags[i_p] & LAND) == 0) &&
                        ((q_ref.flags[i_p] & NAVFAIL) == 0) &&
                        (sst_ref[i_p] >= BAD_FLT + 1.0f);
        }
    }

    /**
     * @brief Same as get_valid_mask_sst but for SST4
     */
    void get_valid_mask_sst4(const l1str &q_ref, int *mask,
                             const float *sst_ref) {
        for (size_t i_p = 0; i_p < npix; i_p++) {
            mask[i_p] = cldmsk::bt_mask(q_ref, i_p, NBANDSIR, ib39) &&
                        cldmsk::bt_mask(q_ref, i_p, NBANDSIR, ib40) &&
                        ((q_ref.flags[i_p] & LAND) == 0) &&
                        ((q_ref.flags[i_p] & NAVFAIL) == 0) &&
                        (sst_ref[i_p] >= BAD_FLT + 1.0f);
        }
    }

    /**
     * @brief
     * Returns the iterval index in bounds_search where a given lat value falls
     * @param q_ref - l1 rec
     * @param i_p - pixel #
     * @return size_t
     */
    size_t find_lat_band_interpolation_point(const l1str &q_ref, size_t i_p) {
        size_t ic;
        {
            const float lat = q_ref.lat[i_p];
            size_t start = 0;
            size_t end = bounds_search.size() - 1;
            if (lat - latwin > bounds_search[0])
                binary_search(bounds_search, lat - latwin, start, end);
            ic = start;
        }
        return ic;
    }

    /**
     * @brief Get the SST using 6 regression coefficients for SST product
     */
    void get_nlsst_v6_cmc(const l1str &q_ref, float *sst_ptr, const int *mask,
                          const float *ref) {
        for (size_t i_p = 0; i_p < npix; i_p++) {
            sst_ptr[i_p] = BAD_FLT;
            if (mask[i_p]) {
                const size_t interp_pnt =
                        find_lat_band_interpolation_point(q_ref, i_p);
                const float lat = q_ref.lat[i_p];
                const float mu = q_ref.csenz[i_p];
                const float satzen = q_ref.senz[i_p];
                const float Bt11 =
                        cldmsk::bt_value(q_ref, i_p, NBANDSIR, ib11) + CtoK;
                const float dBT =
                        Bt11 - cldmsk::bt_value(q_ref, i_p, NBANDSIR, ib12) - CtoK;
                const float cmc = std::max(ref[i_p], 0.0f);
                const float satzdir =
                        (q_ref.pixnum[i_p] < (int32_t)fullscanpix / 2) ? -1.0 : 1.0;
                float *coeffptr =
                        sst_regression_coefficients.data() + ncoeffs * interp_pnt;
                const float lsst = coeffptr[0] + coeffptr[1] * Bt11 +
                                   coeffptr[2] * dBT * cmc +
                                   coeffptr[3] * dBT * ((1.0 / mu) - 1.0) +
                                   coeffptr[4] * satzen * satzdir +
                                   coeffptr[5] * satzen * satzen;
                if (lat < bounds_search[interp_pnt + 1] - latwin ||
                    interp_pnt == (nlatbands - 1)) {
                    sst_ptr[i_p] = lsst - CtoK;
                } else {
                    coeffptr += ncoeffs;  // ncoeffs must be 6 for viirs
                    const float dBtlo = bounds_search[interp_pnt + 1] - latwin;
                    const float hsst = coeffptr[0] + coeffptr[1] * Bt11 +
                                       coeffptr[2] * dBT * cmc +
                                       coeffptr[3] * dBT * ((1.0 / mu) - 1.0) +
                                       coeffptr[4] * satzen * satzdir +
                                       coeffptr[5] * satzen * satzen;
                    sst_ptr[i_p] =
                            lsst + ((lat - dBtlo) / latwin / 2.0) * (hsst - lsst) -
                            CtoK;
                }
            }
        }
    }

    /**
     * @brief Get the SST using 7 regression coefficients for SST product
     */
    void get_nlsst_v7_cmc(const l1str &q_ref, float *sst_ptr, const int *mask,
                          const float *ref) {
        const int mside = q_ref.mside;
        for (size_t i_p = 0; i_p < npix; i_p++) {
            sst_ptr[i_p] = BAD_FLT;
            if (mask[i_p]) {
                const size_t interp_pnt =
                        find_lat_band_interpolation_point(q_ref, i_p);
                const float lat = q_ref.lat[i_p];
                const float mu = q_ref.csenz[i_p];
                const float satzen = q_ref.senz[i_p];
                const float Bt11 = cldmsk::bt_value(q_ref, i_p, NBANDSIR, ib11);
                const float dBT =
                        Bt11 - cldmsk::bt_value(q_ref, i_p, NBANDSIR, ib12);
                //const float solz = q_ref.solz[i_p];
                const float cmc = std::max(ref[i_p], 0.0f);
                const float satzdir =
                        (q_ref.pixnum[i_p] < (int32_t)fullscanpix / 2) ? -1.0 : 1.0;
                float *coeffptr =
                        sst_regression_coefficients.data() + ncoeffs * interp_pnt;
                const float lsst =
                        coeffptr[0] + coeffptr[1] * Bt11 + coeffptr[2] * dBT * cmc +
                        coeffptr[3] * dBT * ((1.0 / mu) - 1.0) +
                        coeffptr[4] * mside + coeffptr[5] * satzen * satzdir +
                        coeffptr[6] * satzen * satzen;
                if (lat < bounds_search[interp_pnt + 1] - latwin ||
                    interp_pnt == (nlatbands - 1)) {
                    sst_ptr[i_p] = lsst;
                } else {
                    coeffptr += ncoeffs;  // ncoeffs must be 7 for modis
                    const float dBtlo = bounds_search[interp_pnt + 1] - latwin;
                    const float hsst = coeffptr[0] + coeffptr[1] * Bt11 +
                                       coeffptr[2] * dBT * cmc +
                                       coeffptr[3] * dBT * ((1.0 / mu) - 1.0) +
                                       coeffptr[4] * mside +
                                       coeffptr[5] * satzen * satzdir +
                                       coeffptr[6] * satzen * satzen;
                    sst_ptr[i_p] =
                            lsst + ((lat - dBtlo) / latwin / 2.0) * (hsst - lsst);
                }
            }
        }
    }
    /**
     * @brief Get the SST b3 using 7 regression coefficients for SST product
     */
    void get_mcsst3_v7(const l1str &q_ref, float *sst_ptr, const int *mask, const float *ref) {
        const int mside = q_ref.mside;
        for (size_t i_p = 0; i_p < npix; i_p++) {
            sst_ptr[i_p] = BAD_FLT;
            if (mask[i_p]) {
                const size_t interp_pnt = find_lat_band_interpolation_point(q_ref, i_p);
                const float lat = q_ref.lat[i_p];
                const float mu = q_ref.csenz[i_p];
                const float satzen = q_ref.senz[i_p];
                const float Bt37 = cldmsk::bt_value(q_ref, i_p, NBANDSIR, ib37);
                const float Bt11 = cldmsk::bt_value(q_ref, i_p, NBANDSIR, ib11);
                const float dBT = Bt37 - cldmsk::bt_value(q_ref, i_p, NBANDSIR, ib12);
                float *coeffptr = sst_regression_coefficients.data() + ncoeffs * interp_pnt;
                const float satzdir =
                        (q_ref.pixnum[i_p] < (int32_t)fullscanpix / 2) ? -1.0 : 1.0;
                const float lsst = coeffptr[0] + coeffptr[1] * Bt11 + coeffptr[2] * dBT +
                                   coeffptr[3] * ((1.0 / mu) - 1.0) + coeffptr[5] * satzen * satzdir +
                                   coeffptr[6] * satzen * satzen + coeffptr[4] * mside;
                if (lat < bounds_search[interp_pnt + 1] - latwin || interp_pnt == (nlatbands - 1)) {
                    sst_ptr[i_p] = lsst;
                } else {
                    coeffptr += ncoeffs;  // ncoeffs must be 7 for modis
                    const float dBtlo = bounds_search[interp_pnt + 1] - latwin;
                    const float hsst = coeffptr[0] + coeffptr[1] * Bt11 + coeffptr[2] * dBT +
                                       coeffptr[3] * ((1.0 / mu) - 1.0) + coeffptr[5] * satzen * satzdir+
                                       coeffptr[6] * satzen * satzen + coeffptr[4] * mside;
                    sst_ptr[i_p] = lsst + ((lat - dBtlo) / latwin / 2.0) * (hsst - lsst);
                }
            }
        }
    }
    /**
     * @brief Get the SST4 using 7 regression coefficients for SST product
     */
    void get_nlsst_modis4(const l1str &q_ref, float *sst_ptr, const int *mask,
                          const float *ref) {
        const int mside = q_ref.mside;
        const int iscan = q_ref.iscan;
        int center = i_center;
        if (iscan > 0) center = i_e;
        bool modis_a_correction = false;
        if (is_aqua)
            modis_a_correction =
                    (l1qstr_modis_a->r[center].detnum == 0) && q_ref.detnum == 0 && center > 0;
        for (size_t i_p = 0; i_p < npix; i_p++) {
            sst_ptr[i_p] = BAD_FLT;
            if (mask[i_p]) {
                const size_t interp_pnt =
                        find_lat_band_interpolation_point(q_ref, i_p);
                const float lat = q_ref.lat[i_p];
                const float mu = q_ref.csenz[i_p];
                const float satzen = q_ref.senz[i_p];
                const float Bt39 = cldmsk::bt_value(q_ref, i_p, NBANDSIR, ib39);
                float Bt40 = cldmsk::bt_value(q_ref, i_p, NBANDSIR, ib40);
                if (modis_a_correction) {
                    float sum = 0.0f;
                    int count = 0;
                    const float Bt40_prev = cldmsk::bt_value(
                            l1qstr_modis_a->r[center - 1], i_p, NBANDSIR, ib40);
                    const float Bt40_next = cldmsk::bt_value(
                            l1qstr_modis_a->r[center + 1], i_p, NBANDSIR, ib40);
                    const int mask_prev = cldmsk::bt_mask(
                            l1qstr_modis_a->r[center - 1], i_p, NBANDSIR, ib40);
                    const int mask_next = cldmsk::bt_mask(
                            l1qstr_modis_a->r[center + 1], i_p, NBANDSIR, ib40);
                    if (mask_prev) {
                        sum += Bt40_prev;
                        count++;
                    }
                    if (mask_next) {
                        sum += Bt40_next;
                        count++;
                    }
                    if (count > 0) {
                        Bt40 = sum / count;
                    } else
                        continue;
                }

                const float dBT = Bt39 - Bt40;
                const float satzdir =
                        (q_ref.pixnum[i_p] < (int32_t)fullscanpix / 2) ? -1.0 : 1.0;
                float *coeffptr =
                        sst_regression_coefficients.data() + ncoeffs * interp_pnt;
                const float lsst =
                        coeffptr[0] + coeffptr[1] * Bt39 + coeffptr[2] * dBT +
                        coeffptr[3] * ((1.0 / mu) - 1.0) + coeffptr[4] * mside +
                        coeffptr[5] * satzen * satzdir +
                        coeffptr[6] * satzen * satzen;
                if (lat < bounds_search[interp_pnt + 1] - latwin ||
                    interp_pnt == (nlatbands - 1)) {
                    sst_ptr[i_p] = lsst;
                } else {
                    coeffptr += ncoeffs;  // ncoeffs must be 7 for modis
                    const float dBtlo = bounds_search[interp_pnt + 1] - latwin;
                    const float hsst =
                            coeffptr[0] + coeffptr[1] * Bt39 + coeffptr[2] * dBT +
                            coeffptr[3] * ((1.0 / mu) - 1.0) + coeffptr[4] * mside +
                            coeffptr[5] * satzen * satzdir +
                            coeffptr[6] * satzen * satzen;
                    sst_ptr[i_p] =
                            lsst + ((lat - dBtlo) / latwin / 2.0) * (hsst - lsst);
                }
            }
            if (is_terra) {
                sst_ptr[i_p] += elecor_terra;
            }
        }
    }

    /**
     * @brief Get the SST3 using 6 regression coefficients for SST product
     */
    void get_nlsst_viirs3(const l1str &q_ref, float *sst_ptr, const int *mask,
                          const float *ref) {
        for (size_t i_p = 0; i_p < npix; i_p++) {
            sst_ptr[i_p] = BAD_FLT;
            if (mask[i_p]) {
                const size_t interp_pnt =
                        find_lat_band_interpolation_point(q_ref, i_p);
                const float lat = q_ref.lat[i_p];
                const float mu = q_ref.csenz[i_p];
                const float satzen = q_ref.senz[i_p];
                const float Bt11 =
                        cldmsk::bt_value(q_ref, i_p, NBANDSIR, ib11) + CtoK;
                const float dBT = cldmsk::bt_value(q_ref, i_p, NBANDSIR, ib37) -
                                  cldmsk::bt_value(q_ref, i_p, NBANDSIR, ib12);
                const float cmc = std::max(ref[i_p], 0.0f);
                const float satzdir =
                        (q_ref.pixnum[i_p] < (int32_t)fullscanpix / 2) ? -1.0 : 1.0;
                float *coeffptr =
                        sst_regression_coefficients.data() + ncoeffs * interp_pnt;
                const float lsst = coeffptr[0] + coeffptr[1] * Bt11 +
                                   coeffptr[2] * dBT * cmc +
                                   coeffptr[3] * ((1.0 / mu) - 1.0) +
                                   coeffptr[4] * satzen * satzdir +
                                   coeffptr[5] * satzen * satzen;
                if (lat < bounds_search[interp_pnt + 1] - latwin ||
                    interp_pnt == (nlatbands - 1)) {
                    sst_ptr[i_p] = lsst - CtoK;
                } else {
                    coeffptr += ncoeffs;  // ncoeffs must be 6 for viirs
                    const float dBtlo = bounds_search[interp_pnt + 1] - latwin;
                    const float hsst = coeffptr[0] + coeffptr[1] * Bt11 +
                                       coeffptr[2] * dBT * cmc +
                                       coeffptr[3] * ((1.0 / mu) - 1.0) +
                                       coeffptr[4] * satzen * satzdir +
                                       coeffptr[5] * satzen * satzen;
                    sst_ptr[i_p] =
                            lsst + ((lat - dBtlo) / latwin / 2.0) * (hsst - lsst) -
                            CtoK;
                }
            }
        }
    }

    /**
     * @brief Process the cloud mask - set flags    *
     * @param q_ref - l1 record
     * @param input - input
     */
    void cloud_mask(const l1str &q_ref, instr *input) {
        cldmsk::product_types case_switcher;
        if (this_product == "SST")
            case_switcher = cldmsk::SST;
        else if (this_product == "SST3")
            case_switcher = cldmsk::SST3;
        else
            case_switcher = cldmsk::SST4;
        // only VIIRS and ONLY SST
        const size_t day2 = day;
        const float xdoy = day2 + 284.0;
        const float xrad = (360.0 / 365.0) * xdoy / RADEG;
        const float subsolar = 23.45 * std::sin(xrad);
        const bool is_viirs =
                sensor ==
                "viirs";  // technically cthe compiler should optimize it out
        const bool is_modis = sensor == "modis";
        for (size_t i_p = 0; i_p < npix; i_p++) {
            if (((q_ref.flags[i_p] & LAND) > 0) ||
                ((q_ref.flags[i_p] & NAVFAIL) > 0)) {
                flags_sst[i_p] |= SSTF_ISMASKED;
                qual_sst[i_p] = 4;
                continue;
            }
            // bt bad
            {
                bool good_bt = true;
                switch (case_switcher) {
                    case cldmsk::SST3:
                        good_bt &= StatsBTs.at("37").get_valid_mask()[i_p];
                    case cldmsk::SST:
                        good_bt &= (StatsBTs.at("11").get_valid_mask()[i_p] &&
                                    StatsBTs.at("12").get_valid_mask()[i_p]);
                        break;
                    case cldmsk::SST4:
                        good_bt &= (StatsBTs.at("39").get_valid_mask()[i_p] &&
                                    StatsBTs.at("40").get_valid_mask()[i_p]);
                        break;
                    default:
                        break;
                }
                if (!good_bt) {
                    flags_sst[i_p] |= SSTF_BTBAD;
                    qual_sst[i_p] = 4;
                    continue;
                }
            }
            bool night = q_ref.solz[i_p] >= cldmsk::solznight;
            bool glint = q_ref.glint_coef[i_p] <= cldmsk::glintmax && !night;
            {
                if (!night) {
                    switch (case_switcher) {
                        case cldmsk::SST4:
                        case cldmsk::SST3:
                            qual_sst[i_p] = 3;
                            break;
                        case cldmsk::SST:
                            if (q_ref.glint_coef[i_p] > cldmsk::glintmax)
                                qual_sst[i_p] = 1;
                        default:
                            break;
                    }
                }
            }
            // bt range
            {
                bool bt_range = true;
                bool bt_range_max_sst4 = true;
                bool bt_range_min_sst4 = true;
                bool bt_range_max_sst3 = true;
                bool bt_range_min_sst3 = true;
                bool bt_range_sst = true;
                if (case_switcher != cldmsk::SST4) {
                    float bt11 = StatsBTs.at("11")("VAR")[i_p];
                    float bt12 = StatsBTs.at("12")("VAR")[i_p];
                    bt_range_sst &=
                            (bt11 <= cldmsk::Btmax) && (bt12 <= cldmsk::Btmax) &&
                            (bt11 >= cldmsk::Btmin) && (bt12 >= cldmsk::Btmin);
                }
                bt_range &= bt_range_sst;
                if (is_modis) {
                    float bt39 = StatsBTs.at("39")("VAR")[i_p];
                    float bt40 = StatsBTs.at("40")("VAR")[i_p];
                    bt_range_min_sst4 &=
                            (bt39 >= cldmsk::Btmin) && (bt40 >= cldmsk::Btmin);
                    bt_range_max_sst4 &=
                            (bt39 <= cldmsk::Btmax) && (bt40 <= cldmsk::Btmax40);
                }

                if (is_viirs) {
                    float bt37 = StatsBTs.at("37")("VAR")[i_p];
                    float bt40 = StatsBTs.at("40")("VAR")[i_p];
                    bt_range_min_sst3 &=
                            (bt37 >= cldmsk::Btmin) && (bt40 >= cldmsk::Btmin);
                    bt_range_max_sst3 &=
                            (bt37 <= cldmsk::Btmax) && (bt40 <= cldmsk::Btmax40);
                }
                switch (case_switcher) {
                    case cldmsk::SST4:
                        bt_range &= bt_range_min_sst4 && bt_range_max_sst4;
                        break;
                    case cldmsk::SST3:
                        bt_range &= bt_range_min_sst3 && bt_range_max_sst3;
                        break;
                    case cldmsk::SST:
                        if (night) {
                            bt_range &= bt_range_min_sst3 && bt_range_max_sst3;
                        } else if (glint) {
                            bt_range &= bt_range_min_sst3;
                        }
                        if (is_modis) {
                            if (night) {
                                bt_range &=
                                        bt_range_min_sst4 && bt_range_max_sst4;
                            } else if (glint) {
                                bt_range &= bt_range_min_sst4;
                            }
                        }
                        break;
                    default:
                        break;
                }
                if (!bt_range) {
                    flags_sst[i_p] |= SSTF_BTRANGE;
                    qual_sst[i_p] = std::max(qual_sst[i_p], (int8_t) 3);
                }
            }
            // BT diff test for SST/SST3
            {
                bool bt_diff = true;
                switch (case_switcher) {
                    case cldmsk::SST3:
                    case cldmsk::SST: {
                        float diffBT = PairOfVarsBTs.at({"11", "12"})()[i_p];
                        if (diffBT < cldmsk::dBtmin || diffBT > cldmsk::dBtmax)
                            bt_diff = false;
                    }
                        break;
                    case cldmsk::SST4: {
                        float diffBT = PairOfVarsBTs.at({"39", "40"})()[i_p];
                        if (diffBT < cldmsk::dBt4min ||
                            diffBT > cldmsk::dBt4max)
                            bt_diff = false;
                    }
                        break;
                    default:
                        break;
                }
                if (!bt_diff) {
                    flags_sst[i_p] |= SSTF_BTDIFF;
                }
            }
            // BT diff for SST4/SST
            {
                switch (case_switcher) {
                    case cldmsk::SST3:
                        break;
                    case cldmsk::SST:
                    case cldmsk::SST4:
                        if (is_modis) {
                            const float BT39 = StatsBTs.at("39")("VAR")[i_p];
                            const float BT40 = StatsBTs.at("40")("VAR")[i_p];
                            const float d3940ref = cldmsk::btrefdiffv6(
                                    i_p, BT39, BT40, &q_ref, fullscanpix);
                            if ((d3940ref < cldmsk::dBtrefmin ||
                                 d3940ref > cldmsk::dBtrefmax) &&
                                night) {
                                flags_sst[i_p] |= SSTF_BT4REFDIFF;
                                qual_sst[i_p] =
                                        std::max(qual_sst[i_p], (int8_t) 3);
                            }
                        }
                        break;
                    default:
                        break;
                }
            }
            // sst range
            {
                bool sst_range = true;
                float sst_val = sst.at(i_center * npix + i_p);
                bool sst_min = sst_val >= cldmsk::SSTmin;
                bool sst_max_night = sst_val <= cldmsk::SSTmaxn;
                bool sst_max_day = sst_val <= cldmsk::SSTmax;

                switch (case_switcher) {
                    case cldmsk::SST4:
                    case cldmsk::SST3:
                        sst_range &= sst_min && sst_max_night;
                        break;
                    case cldmsk::SST:
                        sst_range &= sst_min && ((sst_max_night && night) ||
                                                 (sst_max_day && !night));
                        break;
                    default:
                        break;
                }
                if (!sst_range) {
                    flags_sst[i_p] |= SSTF_SSTRANGE;
                    qual_sst[i_p] = std::max(qual_sst[i_p], (int8_t) 3);
                }
            }
            // sst diff continue
            {
                const float dt_analysis =
                        sst.at(i_center * npix + i_p) - q_ref.sstref[i_p];
                if ((l1_input->evalmask & SSTMODS) == 0) {
                    float dsst = cldmsk::SSTdiff;
                    dsst = input->sstrefdif;
                    if ((dt_analysis < -dsst) &&
                        q_ref.lat[i_p] >= cldmsk::equatorialSouth &&
                        q_ref.lat[i_p] <= cldmsk::equatorialNorth &&
                        q_ref.lon[i_p] >= cldmsk::equatorialWest &&
                        q_ref.lon[i_p] <= cldmsk::equatorialEast) {
                        flags_sst[i_p] |= SSTF_SSTREFDIFF;
                        qual_sst[i_p] = std::max(qual_sst[i_p], (int8_t) 2);
                    }
                };
                if (night) {
                    if (std::abs(dt_analysis) > cldmsk::SSTdiff) {
                        flags_sst[i_p] |= SSTF_SSTREFDIFF;
                        qual_sst[i_p] = std::max(qual_sst[i_p], (int8_t) 2);
                    }
                    if (std::abs(dt_analysis) > cldmsk::SSTvdiff) {
                        flags_sst[i_p] |= SSTF_SSTREFVDIFF;
                        qual_sst[i_p] = std::max(qual_sst[i_p], (int8_t) 3);
                    }
                } else {
                    if (dt_analysis < -cldmsk::SSTdiff ||
                        q_ref.sstref[i_p] < cldmsk::invalid_val) {
                        flags_sst[i_p] |= SSTF_SSTREFDIFF;
                        qual_sst[i_p] = std::max(qual_sst[i_p], (int8_t) 2);
                    }
                    if (dt_analysis < -cldmsk::SSTvdiff ||
                        q_ref.sstref[i_p] < cldmsk::invalid_val) {
                        flags_sst[i_p] |= SSTF_SSTREFVDIFF;
                        qual_sst[i_p] = std::max(qual_sst[i_p], (int8_t) 3);
                    }
                }
            }
            // difference between sst
            {
                float dsst, lim1, lim2;
                bool flag1(false), flag2(false);
                if (is_viirs) {
                    dsst = PairOfSSTs.at({"SST", "SST3"})()[i_p];
                    lim1 = cldmsk::SST3diff1;
                    lim2 = cldmsk::SST3diff2;
                    flag1 = std::abs(dsst) > lim1;
                    flag2 = std::abs(dsst) > lim2;
                }
                if (is_modis) {
                    dsst = PairOfSSTs.at({"SST", "SST4"})()[i_p];
                    // this is actually a bug
                    // dsst = std::abs(dsst);
                    lim1 = cldmsk::SST4diff1;
                    lim2 = cldmsk::SST4diff2;
                    flag1 = dsst < lim1;
                    flag2 = dsst < lim2;
                }
                switch (case_switcher) {
                    case cldmsk::SST4:
                        if (flag1) flags_sst[i_p] |= SSTF_SST4DIFF;
                        if (flag2) flags_sst[i_p] |= SSTF_SST4VDIFF;
                        break;
                    case cldmsk::SST3:
                        if (flag1) flags_sst[i_p] |= SSTF_SST3DIFF;
                        if (flag2) flags_sst[i_p] |= SSTF_SST3VDIFF;
                        break;
                    case cldmsk::SST:
                        if (night) {
                            if (is_modis) {
                                if (flag1) flags_sst[i_p] |= SSTF_SST4DIFF;
                                if (flag2) flags_sst[i_p] |= SSTF_SST4VDIFF;
                            }
                            if (is_viirs) {
                                if (flag1) flags_sst[i_p] |= SSTF_SST3DIFF;
                                if (flag2) flags_sst[i_p] |= SSTF_SST3VDIFF;
                            }
                        }
                        break;
                    default:
                        break;
                }
                if (night) {
                    if (flag1) {
                        qual_sst[i_p] = std::max(qual_sst[i_p], (int8_t) 2);
                    }
                    if (flag2) {
                        qual_sst[i_p] = std::max(qual_sst[i_p], (int8_t) 1);
                    }
                }
            }
            // some vza tests
            {
                if (q_ref.senz[i_p] > cldmsk::hisenz) {
                    flags_sst[i_p] |= SSTF_HISENZ;
                    qual_sst[i_p] = std::max(qual_sst[i_p], (int8_t) 1);
                }
                if (q_ref.pixnum[i_p] < 2 ||
                    q_ref.pixnum[i_p] > ((int32_t)fullscanpix - 3))
                    flags_sst[i_p] |= SSTF_VHISENZ;
                float vhisenz = cldmsk::vhisenz;
                if (is_viirs) {
                    vhisenz = cldmsk::vhisenzv2;
                }
                if (is_terra && is_modis) {
                    if (q_ref.pixnum[i_p] > 1349)
                        flags_sst[i_p] |= SSTF_VHISENZ;
                }
                switch (case_switcher) {
                    case cldmsk::SST4:
                    case cldmsk::SST3:
                        if (q_ref.senz[i_p] > cldmsk::vhisenz) {
                            flags_sst[i_p] |= SSTF_VHISENZ;
                        }
                        break;
                    case cldmsk::SST:
                        if (q_ref.senz[i_p] > vhisenz) {
                            flags_sst[i_p] |= SSTF_VHISENZ;
                        }
                        break;
                    default:
                        break;
                }
                {
                    if ((flags_sst[i_p] & SSTF_VHISENZ) > 0) {
                        qual_sst[i_p] = std::max(qual_sst[i_p], (int8_t) 3);
                    }
                }
            }
            {
                // normalize rho test
                if (case_switcher == cldmsk::SST) {
                    {
                        const float cldthreshold =
                                cldmsk::cldthresh_list.at(sensor);
                        const float val_to_threshold =
                                StatsBTs.at("RHred")("MAXMIN")[i_p];
                        if (glint) {
                            const float dt_analysis =
                                    sst.at(i_center * npix + i_p) -
                                    q_ref.sstref[i_p];
                            const float cos_solz = q_ref.csolz[i_p];
                            if (dt_analysis < cldmsk::cldthresh &&
                                (val_to_threshold == 0.0 ||
                                 val_to_threshold / cos_solz > cldthreshold ||
                                 val_to_threshold < cldmsk::invalid_val ||
                                 val_to_threshold > cldmsk::max_bt)) {
                                flags_sst[i_p] |= SSTF_REDNONUNIF;
                                qual_sst[i_p] =
                                        std::max(qual_sst[i_p], (int8_t) 2);
                            }
                        }
                    }
                    if (is_viirs) {
                        const float rhotRED = StatsBTs.at("red")("VAR")[i_p];
                        const float rhot16 = StatsBTs.at("16")("VAR")[i_p];
                        if (rhotRED > 0.3 && rhot16 >= 0.006 && rhot16 < 0.1 &&
                            ((q_ref.lat[i_p] > (subsolar + 30.0)) ||
                             (q_ref.lat[i_p] < (subsolar - 30.0))) &&
                            glint) {
                            flags_sst[i_p] |= SSTF_CLOUD;
                            qual_sst[i_p] = std::max(qual_sst[i_p], (int8_t) 3);
                        }
                    }
                }
            }
            // uniformity tests
            {
                switch (case_switcher) {
                    case cldmsk::SST3:
                        if (StatsBTs.at("37")("MAXMIN")[i_p] >
                            cldmsk::Bt37unif1)
                            flags_sst[i_p] |= SSTF_BTNONUNIF;
                        if (StatsBTs.at("37")("MAXMIN")[i_p] >
                            cldmsk::Bt37unif2)
                            flags_sst[i_p] |= SSTF_BTVNONUNIF;
                    case cldmsk::SST:
                        if (StatsBTs.at("11")("MAXMIN")[i_p] >
                            cldmsk::Bt11unif1)
                            flags_sst[i_p] |= SSTF_BTNONUNIF;
                        if (StatsBTs.at("11")("MAXMIN")[i_p] >
                            cldmsk::Bt11unif2)
                            flags_sst[i_p] |= SSTF_BTVNONUNIF;
                        if (StatsBTs.at("12")("MAXMIN")[i_p] >
                            cldmsk::Bt12unif1)
                            flags_sst[i_p] |= SSTF_BTNONUNIF;
                        if (StatsBTs.at("12")("MAXMIN")[i_p] >
                            cldmsk::Bt12unif2)
                            flags_sst[i_p] |= SSTF_BTVNONUNIF;
                        break;
                    case cldmsk::SST4:
                        if (StatsBTs.at("39")("MAXMIN")[i_p] >
                            cldmsk::Bt39unif1)
                            flags_sst[i_p] |= SSTF_BTNONUNIF;
                        if (StatsBTs.at("39")("MAXMIN")[i_p] >
                            cldmsk::Bt39unif2)
                            flags_sst[i_p] |= SSTF_BTVNONUNIF;
                        if (StatsBTs.at("40")("MAXMIN")[i_p] >
                            cldmsk::Bt40unif1)
                            flags_sst[i_p] |= SSTF_BTNONUNIF;
                        if (StatsBTs.at("40")("MAXMIN")[i_p] >
                            cldmsk::Bt40unif2)
                            flags_sst[i_p] |= SSTF_BTVNONUNIF;
                        break;
                    default:
                        break;
                }
                {
                    if ((flags_sst[i_p] & SSTF_BTVNONUNIF) > 0) {
                        qual_sst[i_p] = std::max(qual_sst[i_p], (int8_t) 3);
                    }
                    if ((flags_sst[i_p] & SSTF_BTNONUNIF) > 0) {
                        qual_sst[i_p] = std::max(qual_sst[i_p], (int8_t) 1);
                    }
                }
            }
            {
                if (treesum.at(i_p) <= 0.0) {
                    flags_sst[i_p] |= SSTF_CLOUD;
                    qual_sst[i_p] = std::max(qual_sst[i_p], (int8_t) 3);
                }
            }
            {
                if (case_switcher == cldmsk::SST) {
                    if (night && is_modis) {
                        bool decrease_night_qual_non_uniform =
                                (StatsBTs.at("39")("MAXMIN")[i_p] >
                                 cldmsk::Bt39unif1) ||
                                (StatsBTs.at("40")("MAXMIN")[i_p] >
                                 cldmsk::Bt40unif1);
                        if (decrease_night_qual_non_uniform) {
                            qual_sst[i_p] =
                                    std::min((int8_t) 3, ++qual_sst[i_p]);
                        }
                    }
                }
            }
        }
    }

    /**
     * @brief Set the treesum     *
     * @param iscan - scan #
     */
    void set_treesum(int iscan) {
        adt::VarsAtPixel vars;
        {
            for (const auto &test_data_: test_data) {
                vars.vars[test_data_.first] = BAD_FLT;
                vars.masks[test_data_.first] = true;
            }
        }
        std::fill(treesum.begin(), treesum.end(), 0.0);
        std::fill(flags_sst.begin(), flags_sst.end(), 0);
        std::fill(qual_sst.begin(), qual_sst.end(), 0);
        for (size_t i_p = 0; i_p < npix; i_p++) {
            if (valid_mask.at(
                    i_p +
                    i_center * npix))  //(valid_mask.at(i_p + i_center * npix))
            {
                for (const auto &test_var: test_data) {
                    vars.vars.at(test_var.first) = test_var.second[i_p];
                }
                for (const auto &test_mask: test_masks) {
                    vars.masks.at(test_mask.first) = test_mask.second[i_p];
                }
                vars.i_p = i_p;
                vars.i_scan = iscan;
                adt::tree_traversal(treesum.data() + i_p,
                                    flags_sst.data() + i_p, vars,
                                    desicion_tree);
                adt::tree_traversal(desicion_tree);
            }
        }
    };

    /**
     * @brief Ini the ADT
     * @param product_type - product
     */
    void desicion_tree_ini(const std::string &product_type) {
        treesum = std::vector<float>(npix, 0.0f);
        if (desicion_tree == nullptr) {
            desicion_tree = new adt::Treenode();
            const std::string ocdata_root = envset::get_ocdata_root();
            if (product_type == "SST3")
                path_to_cloud_mask_config = ocdata_root + "/" + sensor + "/" +
                                            platform +
                                            "/cal/cloud_mask_sst3.json";
            if (product_type == "SST")
                path_to_cloud_mask_config = ocdata_root + "/" + sensor + "/" +
                                            platform +
                                            "/cal/cloud_mask_sst.json";
            if (product_type == "SST4")
                path_to_cloud_mask_config = ocdata_root + "/" + sensor + "/" +
                                            platform +
                                            "/cal/cloud_mask_sst4.json";
            adt::build_tree(path_to_cloud_mask_config, desicion_tree,
                            test_names);
        }
        for (const auto &test_name: test_names) {
            test_classifications[test_name] = parse_test_name(test_name);
        }
        tree_traversal(desicion_tree);
    };

    /**
     * @brief Ini the data struct later used for cloud masking
     * @param l1rec - l1 rec
     * @param l1qrec - l1 queue
     * @param input - input struct
     * @param product_type - product type
     * @param tests_classifiers - map that contains test  classifiers; each key
     * - test name, each value - keywords corresponding to the test name
     * @param set_adt_data - true or false. Used the data for ADT or only for
     * flag setting.
     */
    void set_data_cldmsk(
            l1str *l1rec, l1qstr *l1qrec, instr *input,
            const std::string &product_type,
            const std::unordered_map<std::string, std::vector<std::string>>
            &tests_classifiers,
            bool set_adt_data) {
        for (const auto &test_pair: tests_classifiers) {
            const std::string &name = test_pair.first;
            const std::vector<std::string> &classifier = test_pair.second;
            const std::vector<std::string> key_vars = {"BT", "RHO", "CLD"};
            std::string pos;
            size_t pos_c = 0;
            bool key_ided = false;
            for (const auto &key_: key_vars) {
                if (boost::algorithm::contains(classifier.at(0), key_)) {
                    pos = key_;
                    key_ided = true;
                    pos_c = classifier.at(0).find(pos);
                    { auto bt = classifier.at(0).substr(pos_c + pos.size()); }
                    break;
                }
            }
            const size_t size_of_classifier = classifier.size();
            if (size_of_classifier == 3)  // difference
            {
                if (boost::algorithm::contains(classifier.at(1), "SST") &&
                    boost::algorithm::contains(classifier.at(2), "SST")) {
                    auto product_name_1 = classifier.at(1);
                    auto sst1 =
                            external_products.at(product_name_1)
                                    ->get_sst(l1rec, l1qrec, input, product_name_1);
                    auto product_name_2 = classifier.at(2);
                    auto sst2 =
                            external_products.at(product_name_2)
                                    ->get_sst(l1rec, l1qrec, input, product_name_2);
                    if (PairOfSSTs.count({product_name_1, product_name_2}) ==
                        0) {
                        PairOfSSTs[{product_name_1, product_name_2}] =
                                bstats::PairOfSST(npix, sst1, sst2, true);
                        PairOfSSTs.at({product_name_1, product_name_2})
                                .get_difference();
                    }
                    if (set_adt_data)
                        test_data[name] =
                                PairOfSSTs.at({product_name_1, product_name_2})();
                } else if (boost::algorithm::contains(classifier.at(1), "BT") &&
                           boost::algorithm::contains(classifier.at(2), "BT")) {
                    assert((classifier.at(1).size() == 4) == 1);
                    assert((classifier.at(2).size() == 4) == 1);
                    auto bt1 = classifier.at(1).substr(2, 2);
                    auto bt2 = classifier.at(2).substr(2, 2);
                    auto ib1 = ibands.at("ib" + bt1);
                    auto ib2 = ibands.at("ib" + bt2);
                    if (PairOfVarsBTs.count({bt1, bt2}) == 0) {
                        PairOfVarsBTs[{bt1, bt2}] =
                                bstats::PairOfVarsBT(npix, ib1, ib2, *l1rec, true);
                        PairOfVarsBTs.at({bt1, bt2}).get_difference(*l1rec);
                    }
                    if (set_adt_data)
                        test_data[name] = PairOfVarsBTs.at({bt1, bt2})();
                } else {
                    adt::print_error_message_for_adt(
                            path_to_cloud_mask_config, name);
                }
            } else if (size_of_classifier == 2)  // stats
            {
                if (boost::algorithm::contains(classifier.at(0), "SST")) {
                    auto product_name = classifier.at(0);
                    auto sst =
                            external_products.at(product_name)
                                    ->get_sst(l1rec, l1qrec, input, product_name);
                    bool std = classifier.at(1) == "STD";
                    bool maxmin = classifier.at(1) == "MAXMIN";
                    bool max_ = classifier.at(1) == "MAX";
                    bool min_ = classifier.at(1) == "MIN";
                    if (StatsSSTs.count(product_name) == 0) {
                        StatsSSTs[product_name] = bstats::StatsSST(
                                npix, n_q_size, bt_box / 2, bt_box / 2, i_center,
                                l1qrec, sst, std, max_, min_, maxmin);
                    }
                    test_data[name] =
                            StatsSSTs.at(product_name)(classifier.at(1));
                } else if (key_ided) {
                    auto bt = classifier.at(0).substr(pos_c + pos.size());
                    int ib = -1;
                    if (!boost::algorithm::contains(classifier.at(0),
                                                    "RHOCIRRUS") &&
                        !boost::algorithm::contains(classifier.at(0), "CLD"))
                        ib = ibands.at("ib" + bt);
                    if (boost::algorithm::contains(classifier.at(0), "CLD"))
                        ib = ibands.at("ibred");
                    bool std = classifier.at(1) == "STD";
                    bool maxmin = classifier.at(1) == "MAXMIN";
                    bool max_ = classifier.at(1) == "MAX";
                    bool min_ = classifier.at(1) == "MIN";
                    if (StatsBTs.count(bt) == 0) {
                        if (boost::algorithm::contains(classifier.at(0), "BT"))
                            StatsBTs[bt] = bstats::StatsVarBand(
                                    npix, n_q_size, bt_box / 2, bt_box / 2,
                                    i_center, l1qrec, ib, NBANDSIR, cldmsk::bt_mask,
                                    cldmsk::bt_value, std, max_, min_, maxmin);
                        if (boost::algorithm::contains(classifier.at(0), "RHO"))
                            StatsBTs[bt] = bstats::StatsVarBand(
                                    npix, n_q_size, bt_box / 2, bt_box / 2,
                                    i_center, l1qrec, ib, n_bands, cldmsk::rho_mask,
                                    cldmsk::rho_value, std, max_, min_, maxmin);
                        if (boost::algorithm::contains(classifier.at(0),
                                                       "RHOCIRRUS"))
                            StatsBTs[bt] = bstats::StatsVarBand(
                                    npix, n_q_size, bt_box / 2, bt_box / 2,
                                    i_center, l1qrec, ib, n_bands,
                                    cldmsk::cirrus_mask, cldmsk::cirrus_value, std,
                                    max_, min_, maxmin);
                        if (boost::algorithm::contains(classifier.at(0), "CLD"))
                            StatsBTs[bt] = bstats::StatsVarBand(
                                    npix, n_q_size, bt_box / 2, bt_box / 2,
                                    i_center, l1qrec, ib, n_bands,
                                    cldmsk::cldrh_mask, cldmsk::cldrh_value, std,
                                    max_, min_, maxmin);
                    } else {
                        if (std) StatsBTs.at(bt).set_std();
                        if (maxmin) StatsBTs.at(bt).set_maxmin();
                        if (min_) StatsBTs.at(bt).set_min();
                        if (max_) StatsBTs.at(bt).set_max();
                    }
                    if (set_adt_data) {
                        test_data[name] = StatsBTs.at(bt)(classifier.at(1));
                        test_masks[name] = StatsBTs.at(bt).get_valid_mask();
                    }
                } else {
                    adt::print_error_message_for_adt(
                            path_to_cloud_mask_config, name);
                }
            } else if (size_of_classifier == 1)  // variable
            {
                if (boost::algorithm::contains(classifier.at(0), "SST")) {
                    auto product_name = classifier.at(0);
                    auto sst =
                            external_products.at(product_name)
                                    ->get_sst(l1rec, l1qrec, input, product_name);
                    if (StatsSSTs.count(product_name) == 0) {
                        StatsSSTs[product_name] =
                                bstats::StatsSST(npix, n_q_size, bt_box / 2,
                                                 bt_box / 2, i_center, l1qrec, sst);
                    }
                    test_data[name] = StatsSSTs.at(product_name)("VAR");
                    continue;
                } else if (key_ided) {
                    auto bt = classifier.at(0).substr(pos_c + pos.size());
                    int ib = -1;
                    if (!boost::algorithm::contains(classifier.at(0),
                                                    "RHOCIRRUS"))
                        ib = ibands.at("ib" + bt);
                    if (StatsBTs.count(bt) == 0) {
                        if (boost::algorithm::contains(classifier.at(0), "BT"))
                            StatsBTs[bt] = bstats::StatsVarBand(
                                    npix, n_q_size, bt_box / 2, bt_box / 2,
                                    i_center, l1qrec, ib, NBANDSIR, cldmsk::bt_mask,
                                    cldmsk::bt_value);
                        if (boost::algorithm::contains(classifier.at(0), "RHO"))
                            StatsBTs[bt] = bstats::StatsVarBand(
                                    npix, n_q_size, bt_box / 2, bt_box / 2,
                                    i_center, l1qrec, ib, n_bands, cldmsk::rho_mask,
                                    cldmsk::rho_value);
                        if (boost::algorithm::contains(classifier.at(0),
                                                       "RHOCIRRUS"))
                            StatsBTs[bt] = bstats::StatsVarBand(
                                    npix, n_q_size, bt_box / 2, bt_box / 2,
                                    i_center, l1qrec, ib, n_bands,
                                    cldmsk::cirrus_mask, cldmsk::cirrus_value);
                    }
                    if (set_adt_data) {
                        test_data[name] = StatsBTs.at(bt)("VAR");
                        test_masks[name] = StatsBTs.at(bt).get_valid_mask();
                    }
                    continue;
                } else if (boost::algorithm::contains(classifier.at(0),
                                                      "Tdeflong")) {
                    std::string bt1 = "11";
                    std::string bt2 = "12";
                    auto ib1 = ibands.at("ib" + bt1);
                    auto ib2 = ibands.at("ib" + bt2);
                    if (PairOfVarsBTs.count({bt1, bt2}) == 0) {
                        PairOfVarsBTs[{bt1, bt2}] = bstats::PairOfVarsBT(
                                npix, ib1, ib2, *l1rec, true, true);
                        PairOfVarsBTs.at({bt1, bt2}).get_difference(*l1rec);
                        PairOfVarsBTs.at({bt1, bt2}).get_ratio(*l1rec);
                    } else {
                        PairOfVarsBTs.at({bt1, bt2}).get_ratio(*l1rec);
                    }
                    if (set_adt_data)
                        test_data[name] =
                                PairOfVarsBTs.at({bt1, bt2}).get_Tdeflong();
                    continue;
                } else {
                    if (cldmsk::get_non_BT_vars.count(name) == 0)
                        adt::print_error_message_for_adt(
                                path_to_cloud_mask_config, name);
                }
                if (set_adt_data)
                    test_data[name] = cldmsk::get_non_BT_vars.at(name)(*l1rec);
            }
        }
    }

    typedef void (SeaSurfaceTemperatureCalcuations::*nlsst_ptr)(const l1str &,
                                                                float *,
                                                                const int *,
                                                                const float *);

    typedef void (SeaSurfaceTemperatureCalcuations::*valid_mask_ptr)(
            const l1str &, int *, const float *);

    typedef float *(
            SeaSurfaceTemperatureCalcuations::*get_desicion_tree_value_ptr)(
            l1str *, l1qstr *, instr *, const std::string &, size_t);

    std::unordered_map<std::string, nlsst_ptr> sst_calls = {
        {"SST_v6_cmc", &SeaSurfaceTemperatureCalcuations::get_nlsst_v6_cmc},
        {"SST3", &SeaSurfaceTemperatureCalcuations::get_nlsst_viirs3},
        {"SST_v7_cmc", &SeaSurfaceTemperatureCalcuations::get_nlsst_v7_cmc},
        {"SST3_v7", &SeaSurfaceTemperatureCalcuations::get_mcsst3_v7},
        {"SST4",
         &SeaSurfaceTemperatureCalcuations::get_nlsst_modis4}};  // ulitmately need to put a model name
    std::unordered_map<std::string, valid_mask_ptr> valid_mask_calls = {
        {"SST", &SeaSurfaceTemperatureCalcuations::get_valid_mask_sst},
        {"SST3", &SeaSurfaceTemperatureCalcuations::get_valid_mask_sst3},
        {"SST4", &SeaSurfaceTemperatureCalcuations::get_valid_mask_sst4}};

   public:
    SeaSurfaceTemperatureCalcuations(l1str *l1rec, l1qstr *l1qrec, instr *input,
                                     const char *product_type) {
        this_product = product_type;
        { std::cout << "Product type to ini " << this_product << std::endl; }
        npix = l1rec->npix;
        const int sensor_id_int = l1rec->l1file->sensorID;
        auto it_s = cldmsk::platforms.find(sensor_id_int);
        if (it_s != cldmsk::platforms.end())
            platform = it_s->second;
        else {
            std::cerr << "Platform is not found. Exiting ... " << sensor_id_int
                      << std::endl;
            exit(EXIT_FAILURE);
        }
        if (platform == "aqua") {
            l1qstr_modis_a = l1qrec;
            is_aqua = true;
            {
                std::cout << "MODIS AQUA correction to BT 4 micron for "
                             "detector #0 will be applied\n";
            }
        } else if (platform == "terra") {
            is_terra = true;
            std::cout << "MODIS TERRA specific corrections will be applied\n";
        }
        {
            // set sensor
            if (cldmsk::modis_sensors.count(sensor_id_int) > 0) {
                sensor = "modis";
                bt_box = cldmsk::bt_box_sizes.at(sensor);
                if (l1_input->resolution == 250) {
                    fullscanpix = 5416;
                } else if (l1_input->resolution == 500) {
                    fullscanpix = 2708;
                }
                if (this_product == "SST") {
                    regression_model = "SST_v7_cmc";
                } else {
                    regression_model = "SST4";
                }
                // dust correction
                {
                    if (l1rec->anc_aerosol && strlen(input->dsdicoeffile)) {
                        dust_correction_data = std::vector<float>(npix, 0.0);
                        {
                            std::cout << "Applying dust correction "
                                      << std::endl;
                        }
                    }
                }
            } else if (cldmsk::viirs_sensors.count(sensor_id_int) > 0) {
                sensor = "viirs";
                bt_box = cldmsk::bt_box_sizes.at(sensor);
                fullscanpix = 3200;
                if (this_product == "SST") {
                    if (platform == "npp") {
                        regression_model = "SST_v6_cmc";
                    } else {
                        regression_model = "SST_v7_cmc";
                    }
                } else {
                    regression_model = "SST3";
                }
            } else {
                std::cerr << "Sensor is not found " << sensor_id_int
                          << std::endl;
                exit(EXIT_FAILURE);
            }
        }
        if (ibands.empty()) {
            cldmsk::read_sst_bands(sensor, ibands);
            n_bands = l1rec->l1file->nbands;
            const std::set<std::string> not_adj_bands = {"ib16", "ibred",
                                                         "ib07", "ib08"};
            for (auto &pair: ibands) {
                if (not_adj_bands.find(pair.first) == not_adj_bands.end())
                    pair.second -= n_bands;
            }
        }
        {
            if (ibands.count("ib11") > 0) {
                ib11 = ibands.at("ib11");
            }
            if (ibands.count("ib12") > 0) {
                ib12 = ibands.at("ib12");
            }
            if (ibands.count("ib37") > 0) {
                ib37 = ibands.at("ib37");
            }
            if (ibands.count("ib39") > 0) {
                ib39 = ibands.at("ib39");
            }
            if (ibands.count("ib40") > 0) {
                ib40 = ibands.at("ib40");
            }
            if (ibands.count("ib85") > 0) {
                ib85 = ibands.at("ib85");
            }
        }
        { desicion_tree_ini(product_type); }
        {
            const std::string sst_coef_file =
                    envset::get_sst_coeffs_path.at(product_type)(input);
            {
                std::cout << "SST Coefficients file for product "
                          << product_type << " " << sst_coef_file << std::endl;
            }
            {
                netCDF::NcFile sst_reg_coefficients_file(sst_coef_file,
                                                         netCDF::NcFile::read);
                netCDF::NcGroupAtt reg_model =
                        sst_reg_coefficients_file.getAtt("regression_model");
                if (!reg_model.isNull()) {
                    size_t attr_len = reg_model.getAttLength();
                    std::string model_name;
                    model_name.reserve(attr_len);
                    reg_model.getValues(model_name);
                    regression_model = model_name;
                    std::cout
                            << "The regression model from the SST coefficient file "
                            << regression_model << std::endl;
                }
                std::vector<float> coeffs;
                netCDF::NcVar coeffs_nc =
                        sst_reg_coefficients_file.getVar("coefficients");
                // bounds
                netCDF::NcVar bounds_nc =
                        sst_reg_coefficients_file.getVar("latbands");
                auto dims = coeffs_nc.getDims();
                auto lbdims = bounds_nc.getDims();
                size_t reserved_size = 1;
                size_t reserved_size_ld = 1;
                std::vector<size_t> dims_all;
                std::vector<size_t> lddims_all;
                for (const auto &dim: dims) {
                    reserved_size *= dim.getSize();
                    dims_all.push_back(dim.getSize());
                }
                for (const auto &dim: lbdims) {
                    reserved_size_ld *= dim.getSize();
                    lddims_all.push_back(dim.getSize());
                }

                coeffs.resize(reserved_size);
                bounds_lat.resize(reserved_size_ld);
                coeffs_nc.getVar(coeffs.data());
                bounds_nc.getVar(bounds_lat.data());
                double scantime = l1rec->scantime;
                int16_t year, month;
                double secs;
                unix2ymds(scantime, &year, & month, & day,& secs);
                   {
                    std::cout << "Scantime - unix time: " << scantime << "; ";
                    std::cout << "scan day : "
                              << day << "; ";
                    std::cout
                            << "scan month : "
                            << month
                            << "; ";
                    cldmsk::month_data() = std::vector<float>(npix, -1.0f);
                    for (size_t i_p = 0; i_p < npix; i_p++) {
                        cldmsk::month_data().at(i_p) =
                                month;
                    }
                    std::cout << "scan year : "
                              << year << "\n";
                }
                const auto month_counter =
                        month - 1;

                if (is_terra) {
                    if (scan_time < cldmsk::scan_time_modis_t_day1) {
                        std::cout << "Applying electronic correction for MODIS "
                                     "TERRA, 2000 October 29"
                                  << std::endl;
                        elecor_terra = cldmsk::el_corr_modis_t_1;
                    } else if (scan_time < cldmsk::scan_time_modis_t_day2) {
                        std::cout << "Applying electronic correction for MODIS "
                                     "TERRA, 2001 July 1"
                                  << std::endl;
                        elecor_terra = cldmsk::el_corr_modis_t_2;
                    }
                }
                {
                    std::cout << "Targeted month " << month_counter + 1
                              << std::endl;
                    sst_regression_coefficients = std::vector<float>(
                            coeffs.data() +
                            dims_all.at(1) * dims_all.at(2) * month_counter,
                            coeffs.data() + dims_all.at(1) * dims_all.at(2) *
                                            (month_counter + 1));
                }
                // sst regression coefficients bounds
                ncoeffs = dims_all.at(2);
                { std::cout << "The regression coefficients used: \n"; }
                for (size_t j = 0; j < dims_all.at(1); j++) {
                    {
                        const size_t i_lat_st = j * lddims_all.at(1);
                        const size_t i_lat_e = j * lddims_all.at(1) + 1;
                        const float lat_st = bounds_lat.at(i_lat_st);
                        const float lat_e = bounds_lat.at(i_lat_e);
                        std::cout << "Lat_start : " << lat_st
                                  << "; lat_end : " << lat_e
                                  << "; Coefficients : ";
                    }
                    for (size_t k = 0; k < dims_all.at(2); k++) {
                        const size_t index = j * dims_all.at(2) + k;
                        std::cout << "C" << k << " = "
                                  << sst_regression_coefficients.at(index)
                                  << "; ";
                    }
                    std::cout << "\n";
                }
                bounds_search.push_back(bounds_lat.at(0));
                nlatbands = lddims_all.at(0);
                for (size_t j = 0; j < lddims_all.at(0); j++) {
                    for (size_t k = 0; k < lddims_all.at(1); k++) {
                        const size_t index = j * lddims_all.at(1) + k;
                        if (k == 1) {
                            bounds_search.push_back(bounds_lat.at(index));
                        }
                    }
                    std::cout << "\n";
                }
                // reading the SSES bias
                {
                    const std::string sses_luts_path =
                            envset::get_sses_coeffs_path.at(product_type)(input);
                    sses_bias_data = std::make_unique<cldmsk::SSESData>(
                            sses_luts_path, npix);
                }
            }
        }
    };

    float *get_sst(l1str *l1rec, l1qstr *l1qrec, instr *input,
                   const std::string &product_type) {
        const int iscan = l1rec->iscan;
        if (sst.empty()) {
            {
                n_q_size = l1qrec->nq;
                i_center = n_q_size / 2;
                i_s = std::min(std::max(0, (int) i_center - (int) bt_box / 2),
                               (int) n_q_size - 1);
                i_e = std::max(std::min((int) n_q_size - 1, (int)i_center + (int) bt_box / 2),
                               (int)0);
            }
            sst = std::vector<float>(n_q_size * npix, 0.0f);
            valid_mask = std::vector<int>(n_q_size * npix, 0);
            for (size_t i_q = i_s; i_q <= i_e; i_q++) {
                l1str &q_ref = l1qrec->r[i_q];
                float *sst_ref =
                        get_reference(l1rec, l1qrec, input, product_type, i_q);
                (this->*valid_mask_calls.at(product_type))(
                        q_ref, valid_mask.data() + i_q * npix, sst_ref);
                (this->*sst_calls.at(regression_model))(
                        q_ref, sst.data() + i_q * npix,
                        valid_mask.data() + i_q * npix, sst_ref);
                modis_dust_correction(q_ref, i_q);
            }
            processed_scan_sst = iscan;
        } else {
            if (processed_scan_sst == iscan ||
                processed_scan_sst != iscan - 1) {
                return sst.data() + i_center * npix;
            }
            for (size_t i_q = i_s; i_q < i_e; i_q++) {
                const size_t index = (i_q + 1) * npix;
                std::copy(valid_mask.begin() + index,
                          valid_mask.begin() + index + npix,
                          valid_mask.begin() + index - npix);
            }
            for (size_t i_q = i_s; i_q < i_e; i_q++) {
                const size_t index = (i_q + 1) * npix;
                std::copy(sst.begin() + index, sst.begin() + index + npix,
                          sst.begin() + index - npix);
            }
            {
                size_t i_q = i_e;
                l1str &q_ref = l1qrec->r[i_q];
                float *sst_ref =
                        get_reference(l1rec, l1qrec, input, product_type, i_q);
                (this->*valid_mask_calls.at(product_type))(
                        q_ref, valid_mask.data() + i_q * npix, sst_ref);
                (this->*sst_calls.at(regression_model))(
                        q_ref, sst.data() + i_q * npix,
                        valid_mask.data() + i_q * npix, sst_ref);
                modis_dust_correction(q_ref, i_q);
            }
            processed_scan_sst = iscan;
        }
        return sst.data() + i_center * npix;

    }

    float *get_std(l1str *l1rec, l1qstr *l1qrec, instr *input,
                   const std::string &product_type) {
        get_bias(l1rec, l1qrec, input, product_type);
        return sses_bias_data->stdv_out.data();
    }

    float *get_bias(l1str *l1rec, l1qstr *l1qrec, instr *input,
                    const std::string &product_type) {
        const int iscan = l1rec->iscan;
        if (processed_scan_sses == iscan || processed_scan_sses != iscan - 1) {
            return sses_bias_data->bias_out.data();
        } else {
            float *diff;
            const int16_t *flags_ptr = get_flags_sst(l1rec, l1qrec, input, product_type);
            if (this_product == "SST")
                diff = PairOfVarsBTs.at({"11", "12"})();
            else if (this_product == "SST3")
                diff = PairOfVarsBTs.at({"37", "12"})();
            else if (this_product == "SST4")
                diff = PairOfVarsBTs.at({"39", "40"})();
            sses_bias_data->calculate_sses_bias_stats(
                    diff, sst.data() + npix * i_center, l1rec->solz, l1rec->senz,
                    day, l1rec->lat, qual_sst.data(), flags_ptr,
                    l1rec->glint_coef);
        }
        processed_scan_sses = iscan;
        return sses_bias_data->bias_out.data();
    }

    float *get_bias_mean(l1str *l1rec, l1qstr *l1qrec, instr *input,
                         const std::string &product_type) {
        get_bias(l1rec, l1qrec, input, product_type);
        return sses_bias_data->bias_mean_out.data();
    }

    int16_t *get_flags_sst(l1str *l1rec, l1qstr *l1qrec, instr *input,
                           const std::string &product_type) {
        const int iscan = l1rec->iscan;
        if (flags_sst.empty()) {
            flags_sst = std::vector<int16_t>(npix, 0);
            qual_sst = std::vector<int8_t>(npix, 0);
            external_products[this_product] =
                    this;  //  cloud mask may need other's SSTs
            {
                std::cout << "Setting up flags_sst for " << product_type
                          << std::endl;
            }
            {
                set_data_cldmsk(l1rec, l1qrec, input, product_type,
                                test_classifications, true);
                {
                    const auto &vars =
                            needed_stats_for_cloud_mask.at({sensor, product_type});
                    std::unordered_map<std::string, std::vector<std::string>>
                            non_tree_test_classifications;
                    for (const auto &test_name: vars) {
                        non_tree_test_classifications[test_name] =
                                parse_test_name(test_name);
                    }
                    set_data_cldmsk(l1rec, l1qrec, input, product_type,
                                    non_tree_test_classifications, false);
                    {
                        std::cout << "The cloud mask has been set for "
                                  << product_type << std::endl;
                    }
                }
                {
                    set_modis_a_bt_40_data(l1rec, l1qrec, input, product_type);
                    set_treesum(iscan);
                    cloud_mask(*l1rec, input);
                }
            }
            processed_scan_flags = iscan;
        } else {
            if (processed_scan_flags == iscan ||
                iscan - 1 != processed_scan_flags) {
                return flags_sst.data();
            }
            for (auto &pairSST: PairOfSSTs) {
                pairSST.second.sst1 =
                        external_products.at(pairSST.first.first)
                                ->get_sst(l1rec, l1qrec, input, pairSST.first.first);
                pairSST.second.sst2 =
                        external_products.at(pairSST.first.second)
                                ->get_sst(l1rec, l1qrec, input, pairSST.first.second);
                pairSST.second.get_difference();
            }
            for (auto &pairBT: PairOfVarsBTs) {
                pairBT.second.get_difference(*l1rec);
            }
            for (auto &statSST: StatsSSTs) {
                statSST.second.sst =
                        external_products.at(statSST.first)
                                ->get_sst(l1rec, l1qrec, input, statSST.first);
                statSST.second.rearrange();
            }
            for (auto &statBT: StatsBTs) {
                statBT.second.rearrange();
            }
            {
                set_modis_a_bt_40_data(l1rec, l1qrec, input, product_type);
                set_treesum(iscan);
                cloud_mask(*l1rec, input);
            }
            processed_scan_flags = iscan;
        }
        return flags_sst.data();
    }

    void set_ext_products(
            std::unique_ptr<SeaSurfaceTemperatureCalcuations> &ext_product,
            const std::string &ext_product_name) {
        external_products[ext_product_name] = ext_product.get();
        {
            std::cout << "External product reference added : "
                      << ext_product_name << std::endl;
        }
    };

    size_t get_number_of_external_products() {
        return external_products.size();
    }

    int8_t *get_quality_flags(l1str *l1rec, l1qstr *l1qrec, instr *input,
                              const std::string &product_type) {
        get_flags_sst(l1rec, l1qrec, input, product_type);
        return qual_sst.data();
    }

    int16_t *get_counts(l1str *l1rec, l1qstr *l1qrec, instr *input,
                        const std::string &product_type) {
        get_bias(l1rec, l1qrec, input, product_type);
        return sses_bias_data->counts_out.data();
    }

    float *get_dust_correction(l1str *l1rec, l1qstr *l1qrec, instr *input,
                               const std::string &product_type) {
        if (dust_correction_data.empty()) {
            std::cerr << "Dust correction is not availible. Exiting ... "
                      << std::endl;
            exit(EXIT_FAILURE);
        }
        return dust_correction_data.data();
    }

    ~SeaSurfaceTemperatureCalcuations() { deallocate_tree(desicion_tree); }
};

/**
 * @brief
 * The driver class to interface with the pure legacy C code.
 * @param sst - the products stored
 * @param ini_product - initialize the product
 */
class SST {
public:
    std::unordered_map<std::string,
            std::unique_ptr<SeaSurfaceTemperatureCalcuations>>
            sst;
    std::set<std::string> initilized_products;
    std::unordered_map<
            std::string, std::unordered_map<std::string, std::vector<std::string>>>
            products_needed = {{"SST",  {{"viirs", {"SST3"}}, {"modis", {"SST4"}}}},
                               {"SST3", {{"viirs", {"SST"}}}},
                               {"SST4", {{"modis", {"SST"}}}}};

    void ini_product(l1str *l1rec, l1qstr *l1qrec, instr *input,
                     const char *product_type) {
        if (initilized_products.count(product_type) > 0) return;
        std::vector<std::string> needed_products;
        // setting up needed products
        {
            const auto sensorId = l1rec->l1file->sensorID;
            auto it = products_needed.find(product_type);
            if (it != products_needed.end()) {
                const auto &entry = it->second;
                auto it_e = entry.find(cldmsk::get_sensor(sensorId));
                if (it_e != entry.end()) {
                    needed_products = it_e->second;
                }
            }
        }
        if (sst.count(product_type) == 0) {
            sst[product_type] =
                    std::make_unique<SeaSurfaceTemperatureCalcuations>(
                            l1rec, l1qrec, input, product_type);
            for (const auto &needed_product: needed_products) {
                ini_product(l1rec, l1qrec, input, needed_product.c_str());
            }
        }

        if (sst.at(product_type)->get_number_of_external_products() <
            needed_products.size()) {
            for (const auto &needed_product: needed_products) {
                sst.at(product_type)
                        ->set_ext_products(sst.at(needed_product), needed_product);
            }
        }
        initilized_products.insert(product_type);
    }

    float *get_sst(l1str *l1rec, l1qstr *l1qrec, instr *input,
                   const char *product_type) {
        ini_product(l1rec, l1qrec, input, product_type);
        return sst.at(product_type)
                ->get_sst(l1rec, l1qrec, input, std::string(product_type));
    }

    float *get_sst_stdev(l1str *l1rec, l1qstr *l1qrec, instr *input,
                         const char *product_type) {
        ini_product(l1rec, l1qrec, input, product_type);
        return sst.at(product_type)
                ->get_std(l1rec, l1qrec, input, std::string(product_type));
    }

    float *get_sst_bias(l1str *l1rec, l1qstr *l1qrec, instr *input,
                        const char *product_type) {
        ini_product(l1rec, l1qrec, input, product_type);
        return sst.at(product_type)
                ->get_bias(l1rec, l1qrec, input, std::string(product_type));
    }

    float *get_sst_mean_bias(l1str *l1rec, l1qstr *l1qrec, instr *input,
                             const char *product_type) {
        ini_product(l1rec, l1qrec, input, product_type);
        return sst.at(product_type)
                ->get_bias_mean(l1rec, l1qrec, input, std::string(product_type));
    }

    float *get_dust_correction(l1str *l1rec, l1qstr *l1qrec, instr *input,
                               const char *product_type) {
        ini_product(l1rec, l1qrec, input, product_type);
        return sst.at(product_type)
                ->get_dust_correction(l1rec, l1qrec, input,
                                      std::string(product_type));
    }

    int16_t *get_flags_sst(l1str *l1rec, l1qstr *l1qrec, instr *input,
                           const char *product_type) {
        ini_product(l1rec, l1qrec, input, product_type);
        return sst.at(product_type)
                ->get_flags_sst(l1rec, l1qrec, input, std::string(product_type));
    }

    int8_t *get_qual_flags(l1str *l1rec, l1qstr *l1qrec, instr *input,
                           const char *product_type) {
        ini_product(l1rec, l1qrec, input, product_type);
        return sst.at(product_type)
                ->get_quality_flags(l1rec, l1qrec, input,
                                    std::string(product_type));
    }

    int16_t *get_counts(l1str *l1rec, l1qstr *l1qrec, instr *input,
                        const char *product_type) {
        ini_product(l1rec, l1qrec, input, product_type);
        return sst.at(product_type)
                ->get_counts(l1rec, l1qrec, input, std::string(product_type));
    }
};
namespace cld_mask_sst {
    SST sst_obj;
}

void call_sst(l1str *l1rec, l1qstr *l1qrec, instr *input, const char *product_type, float **sst_out_2) {
    *sst_out_2 = cld_mask_sst::sst_obj.get_sst(l1rec, l1qrec, input, product_type);
}
void call_sses_bias(l1str *l1rec, l1qstr *l1qrec, instr *input, const char *product_type,
                    float **sses_bias_out) {
    *sses_bias_out = cld_mask_sst::sst_obj.get_sst_bias(l1rec, l1qrec, input, product_type);
}
void call_sses_std(l1str *l1rec, l1qstr *l1qrec, instr *input, const char *product_type,
                   float **sses_std_out) {
    *sses_std_out = cld_mask_sst::sst_obj.get_sst_stdev(l1rec, l1qrec, input, product_type);
}
void call_sst_flags(l1str *l1rec, l1qstr *l1qrec, instr *input, const char *product_type,
                    int16_t **sst_flags) {
    *sst_flags = cld_mask_sst::sst_obj.get_flags_sst(l1rec, l1qrec, input, product_type);
}
void call_qual_flags(l1str *l1rec, l1qstr *l1qrec, instr *input, const char *product_type,
                     int8_t **qual_flags) {
    *qual_flags = cld_mask_sst::sst_obj.get_qual_flags(l1rec, l1qrec, input, product_type);
}
void call_sses_bias_mean(l1str *l1rec, l1qstr *l1qrec, instr *input, const char *product_type,
                         float **sses_bias_mean_out) {
    *sses_bias_mean_out = cld_mask_sst::sst_obj.get_sst_mean_bias(l1rec, l1qrec, input, product_type);
}
void call_dust_correction(l1str *l1rec, l1qstr *l1qrec, instr *input, const char *product_type,
                          float **dust_correction) {
    *dust_correction = cld_mask_sst::sst_obj.get_sst_mean_bias(l1rec, l1qrec, input, product_type);
}
void call_sses_counts(l1str *l1rec, l1qstr *l1qrec, instr *input, const char *product_type,
                      int16_t **sses_counts) {
    *sses_counts = cld_mask_sst::sst_obj.get_counts(l1rec, l1qrec, input, product_type);
}
