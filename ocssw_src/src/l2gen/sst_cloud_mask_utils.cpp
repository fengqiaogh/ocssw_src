#include "sst_cloud_mask_utils.hpp"
#include "d3940tref.h"

namespace envset {
std::string get_ocdata_root() {
    if (const char* env_p = std::getenv("OCDATAROOT")) {
        return env_p;
    }
    std::cerr << "--Error-: OCDATAROOT env variable NOT found; exiting ..." << std::endl;
    exit(EXIT_FAILURE);
    return "";
}
std::string get_sst_coeff_file(const std::string& sst_type, const instr* input) {
    if (sst_type == "SST")
        return input->sstcoeffile;
    else if (sst_type == "SST4")
        return input->sst4coeffile;
    else if (sst_type == "SST3")
        return input->sst3coeffile;
    else
        fprintf(stderr, "-E-: not supported SST type %s\n", sst_type.c_str());
    exit(EXIT_FAILURE);
}

std::string get_sst_sses_file(const std::string& sst_type, const instr* input) {
    if (sst_type == "SST")
        return input->sstssesfile;
    else if (sst_type == "SST4")
        return input->sst4ssesfile;
    else if (sst_type == "SST3")
        return input->sst3ssesfile;
    else
        fprintf(stderr, "-E-: not supported SST type %s\n", sst_type.c_str());
    exit(EXIT_FAILURE);
}

}  // namespace envset

namespace cldmsk {
    std::vector<float> month_data_;  // month_data;
    std::vector<float>& month_data() {
        return month_data_;
    }

    std::unordered_map<std::string, float* (*)(const l1str&)> get_non_BT_vars_{};
    std::unordered_map<std::string, std::unordered_set<int>> all_sensors_{};
    std::unordered_set<int> modis_sensors_{};
    std::unordered_set<int> viirs_sensors_{};
    std::unordered_map<int, std::string> platforms_{};
    std::unordered_map<std::string, float> cldthresh_list_{};    
    std::unordered_map<std::string, std::unordered_map<std::string, int>> bands_set_{};
    std::unordered_map<std::string, size_t> bt_box_sizes_{};

    std::string get_sensor(int key) {
        const std::unordered_map<std::string, std::unordered_set<int>> sensors = all_sensors();
        for (const auto& sens_data : sensors) {
            if (sens_data.second.count(key) > 0)
                return sens_data.first;
        }
        std::cerr
            << "--Error-: Sensor is  Not Found; Currently, the code process only VIIRS and MODIS; exiting ...  "
            << std::endl;
        exit(EXIT_FAILURE);
        return "";
    }

    bool cirrus_mask(const l1str& l1str_, int pixel, int nbands, int ib) {
        return l1str_.rho_cirrus[pixel] > invalid_val;
    };
    float cirrus_value(const l1str& l1str_, int pixel, int nbands, int ib) {
        return l1str_.rho_cirrus[pixel];
    };
    // bts
    bool bt_mask(const l1str& l1str_, int pixel, int nbands, int ib) {
        const float bt = l1str_.Bt[pixel * nbands + ib];
        return bt > min_bt && bt < max_bt;
    };
    float bt_value(const l1str& l1str_, int pixel, int nbands, int ib) {
        return l1str_.Bt[pixel * nbands + ib];
    };
    // rho (RSMAS)
    bool rho_mask(const l1str& l1str_, int pixel, int nbands, int ib) {
        const float lt = l1str_.Lt[pixel * nbands + ib];
        return lt > min_lt && lt < max_lt;
    };
    float rho_value(const l1str& l1str_, int pixel, int nbands, int ib) {
        const float lt = l1str_.Lt[pixel * nbands + ib];
        const float fo = l1str_.Fo[ib];
        const float csolz = l1str_.csolz[pixel];
        const float rho = OEL_PI * lt / fo / csolz;
        return rho;
    };

    float cldrh_value(const l1str& l1str_, int pixel, int nbands, int ib) {
        const float lt = l1str_.Lt[pixel * nbands + ib];
        const float fo = l1str_.Fo[ib];
        // const float csolz = l1str_.csolz[pixel];
        const float tg_sol = l1str_.tg_sol[pixel * nbands + ib];
        const float tg_sen = l1str_.tg_sen[pixel * nbands + ib];
        const float t_sen = l1str_.t_sen[pixel * nbands + ib];
        const float t_sol = l1str_.t_sol[pixel * nbands + ib];
        const float cldrh = OEL_PI * lt / fo / tg_sol / tg_sen / t_sol / t_sen;  // /csolz
        return cldrh;
    };
    bool cldrh_mask(const l1str& l1str_, int pixel, int nbands, int ib) {
        const float lt = l1str_.Lt[pixel * nbands + ib];
        return lt > 0.0 && lt < max_lt;
    };

    float btrefdiffv6(int32_t ip, float BT39, float BT40, const l1str *l1rec, int fullscanpix) {
        float diff;
        float satzdir;
        float senz;
        float tref;
        satzdir = (l1rec->pixnum[ip] < fullscanpix / 2) ? -1.0 : 1.0;
        senz = l1rec->senz[ip] * satzdir;
        tref = linterp(sza, trefv6, NSENZ, senz);
        diff = BT39 - BT40 - tref;
        return (diff);
    }

    void get_var_mask(int *mask, const l1str &l1str, size_t npix, int nbands, int ib, get_valid get_mask) {
        for (size_t i = 0; i < npix; i++) {
            mask[i] = get_mask(l1str, i, nbands, ib);
        }
    }

    void get_var_vals(float *BT, const l1str &l1str, size_t npix, int nbands, int ib, get_value get_val) {
        for (size_t i = 0; i < npix; i++) {
            BT[i] = get_val(l1str, i, nbands, ib);
        }
    }

    void get_window_1D_max(float *maxs_1d, const float *inp1d, const int *mask1d, size_t npix, size_t radius) {
        for (size_t i_p = 0; i_p < npix; i_p++) {
            maxs_1d[i_p] = BAD_FLT;
            const size_t i_p_s = std::max(0, (int) i_p - (int) radius);
            const size_t i_p_e = std::min(npix - 1, i_p + radius);
            for (size_t i = i_p_s; i <= i_p_e; i++) {
                if (mask1d[i]) {
                    const float val = inp1d[i];
                    if (val > maxs_1d[i_p]) {
                        maxs_1d[i_p] = val;
                    }
                }
            }
        }
    }

    void
    get_total_2d_max(float *max_global, const float *inp2d, size_t npix, size_t nscan, size_t center, size_t radius) {
        const size_t i_q_s = std::max(0, (int) center - (int) radius);
        const size_t i_q_e = std::min(nscan - 1, center + radius);
        for (size_t i_p = 0; i_p < npix; i_p++) {
            max_global[i_p] = BAD_FLT;
            for (size_t i_q = i_q_s; i_q <= i_q_e; i_q++) {
                const size_t ind = i_p + i_q * npix;
                const float val = inp2d[ind];
                if (val > max_global[i_p]) {
                    max_global[i_p] = val;
                }
            }
        }
    }

    void get_window_1D_min(float *mins_1d, const float *inp1d, const int *mask1d, size_t npix, size_t radius) {
        const float max_possible_val = std::abs(BAD_FLT);
        for (size_t i_p = 0; i_p < npix; i_p++) {
            mins_1d[i_p] = max_possible_val;
            const size_t i_p_s = std::max(0, (int) i_p - (int) radius);
            const size_t i_p_e = std::min(npix - 1, i_p + radius);
            for (size_t i = i_p_s; i <= i_p_e; i++) {
                if (mask1d[i]) {
                    const float val = inp1d[i];
                    if (val < mins_1d[i_p]) {
                        mins_1d[i_p] = val;
                    }
                }
            }
        }
    }

    void
    get_total_2d_min(float *min_global, const float *inp2d, size_t npix, size_t nscan, size_t center, size_t radius) {
        const size_t i_q_s = std::max(0, (int) center - (int) radius);
        const size_t i_q_e = std::min(nscan - 1, center + radius);
        const float max_possible_val = std::abs(BAD_FLT);
        for (size_t i_p = 0; i_p < npix; i_p++) {
            min_global[i_p] = max_possible_val;
            for (size_t i_q = i_q_s; i_q <= i_q_e; i_q++) {
                const size_t ind = i_p + i_q * npix;
                const float val = inp2d[ind];
                if (val < min_global[i_p]) {
                    min_global[i_p] = val;
                }
            }
        }
    }

    void get_std_box(const float *inp2d, const int *mask2d, size_t npix, size_t nscan, size_t radius_x, size_t radius_y,
                     float *out, size_t center, l1qstr *l1qrec) {
        std::vector<double> sums(npix * nscan, 0.0);
        std::vector<double> sums2(npix * nscan, 0.0);
        std::vector<size_t> counts(npix * nscan, 0);
        for (size_t i_q = 0; i_q < nscan; i_q++) {
            for (size_t i_p = 0; i_p < npix; i_p++) {
                const size_t i_p_s = std::max(0, (int) i_p - (int) radius_x);
                const size_t i_p_e = std::min(npix - 1, i_p + radius_x);
                const size_t ind = i_p + i_q * npix;
                for (size_t i = i_p_s; i <= i_p_e; i++) {
                    const size_t index = i + i_q * npix;
                    if (mask2d[index]) {
                        sums.at(ind) += (inp2d[index]);
                        sums2.at(ind) += (inp2d[index]) * (inp2d[index]);
                        counts.at(ind) += 1;
                    }
                }
            }
        }
        const size_t i_q_s = std::max(0, (int) center - (int) radius_y);
        const size_t i_q_e = std::min(nscan - 1, center + radius_y);
        for (size_t i_p = 0; i_p < npix; i_p++) {
            double sum_local = 0.0f;
            double sum_local2 = 0.0f;
            size_t count_local = 0;
            for (size_t i_q = i_q_s; i_q <= i_q_e; i_q++) {
                const size_t ind = i_p + i_q * npix;
                count_local += counts.at(ind);
                sum_local += sums.at(ind);
                sum_local2 += sums2.at(ind);
            }
            if (count_local > 2) {
                double diff = (sum_local2 - sum_local * sum_local / (int) count_local) / ((int) count_local - 1);
                if (diff <= 0) {
                    out[i_p] = BAD_FLT;
                } else {
                    out[i_p] = std::sqrt(diff);
                }
            } else {
                out[i_p] = BAD_FLT;
            }
        }
    }

    void read_sst_bands(const std::string& sat, std::unordered_map<std::string, int>& bands) {
        std::unordered_map<std::string, std::unordered_map<std::string, int>> bands_sets = bands_set();
        if (bands_sets.count(sat) == 0) {
            std::cerr << "--Error-: Satellite " << sat << " is not supported; exiting ...  " << std::endl;
            exit(EXIT_FAILURE);
        }
        const std::unordered_map<std::string, int> bands_sat = bands_sets.at(sat);
        for (const auto& var_pair : bands_sat) {
            const std::string key = var_pair.first;
            const int value = var_pair.second;
            bands[key] = bindex_get(value);
        }
    }

    const std::unordered_map<std::string, float> & cldthresh_list() {
        return cldthresh_list_;
    }  // threshold for cold SST test
   const std::unordered_map<std::string, std::unordered_map<std::string, int>> &  bands_set() {
        return bands_set_;  // bands WV for VIIRS/MODIS needed to produce SST/Cloud mask
    }
   const std::unordered_map<std::string, size_t>  &bt_box_sizes() {
        return bt_box_sizes_;
    }  // bt box size

   const std::unordered_set<int>& modis_sensors() {
        return modis_sensors_;
    }  // sensors/platforms INT values and corresponding string values
   const std::unordered_set<int>& viirs_sensors() {
        return viirs_sensors_;
    }
    const std::unordered_map<std::string, std::unordered_set<int>>& all_sensors() {
        return all_sensors_;
    }
    const std::unordered_map<int, std::string> & platforms() {
        return platforms_;
    }
    const std::unordered_map<std::string, float* (*)(const l1str&)>& get_non_BT_vars() {
        return get_non_BT_vars_;
    }
    void init_parameters() {
        get_non_BT_vars_ = {{{"solz", [](const l1str& l1rec) { return l1rec.solz; }},
                            {"senz", [](const l1str& l1rec) { return l1rec.senz; }},
                            {"wv", [](const l1str& l1rec) { return l1rec.wv; }},
                            {"glintcoef", [](const l1str& l1rec) { return l1rec.glint_coef; }},
                            {"lat", [](const l1str& l1rec) { return l1rec.lat; }},
                            {"lon", [](const l1str& l1rec) { return l1rec.lon; }},
                            {"month", [](const l1str& l1rec) { return month_data().data(); }}}};

        modis_sensors_ = {MODIST, MODISA};
        viirs_sensors_ = {VIIRSJ1, VIIRSJ2, VIIRSN};
        all_sensors_ = {{"modis", modis_sensors_}, {"viirs", viirs_sensors_}};
        platforms_ = {{MODIST, "terra"}, {MODISA, "aqua"}, {VIIRSN, "npp"}, {VIIRSJ1, "j1"}, {VIIRSJ2, "j2"}};
        cldthresh_list_ ={{"modis", 0.01}, {"viirs", 0.04}};
        bands_set_ = {{"modis",
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
        bt_box_sizes_ = {{"modis", 3}, {"viirs", 5}};
    }
}  // namespace cldmsk