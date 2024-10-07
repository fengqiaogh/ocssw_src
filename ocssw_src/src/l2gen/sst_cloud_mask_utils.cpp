#include "sst_cloud_mask_utils.hpp"
#include "d3940tref.h"

namespace envset {
    std::string get_ocdata_root() {
        if (const char *env_p = std::getenv("OCDATAROOT")) {
            return env_p;
        }
        std::cerr << "--Error-: OCDATAROOT env variable NOT found; exiting ..." << std::endl;
        exit(EXIT_FAILURE);
        return "";
    }

}

namespace cldmsk {
    std::vector<float> month_data_; // month_data;
    std::vector<float> &month_data() { return month_data_; }

    std::string get_sensor(int key) {
        for (const auto &sens_data: all_sensors) {
            if (sens_data.second.count(key) > 0)
                return sens_data.first;
        }
        std::cerr << "--Error-: Sensor is  Not Found; Currently, the code process only VIIRS and MODIS; exiting ...  "
                  << std::endl;
        exit(EXIT_FAILURE);
        return "";
    }

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

    void read_sst_bands(const std::string &sat, std::unordered_map<std::string, int> &bands) {
        for (const auto &var_pair: bands_set.at(sat)) {
            const std::string key = var_pair.first;
            const int value = var_pair.second;
            bands[key] = bindex_get(value);
        }
    }
}