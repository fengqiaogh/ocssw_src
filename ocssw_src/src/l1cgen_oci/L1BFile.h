#ifndef L1BFILE_H
#define L1BFILE_H
#include "bin_l1b.h"
#include "filehandle.h"
#include "l1.h"
#include "utils.h"
#include "netcdf.h"
#include <cassert>
#include "l1_oci_private.h"
class L1BFile {
    l1str l1rec;
    filehandle l1file;
    std::vector<int> bin_indexes;
    std::vector<std::vector<std::pair<int, double>>> area_weights{};
    std::vector<float> height_data_corrected{}, senz_corrected{}, sena_corrected{};
    std::string ifile;
    int qual_red_id, qual_blue_id, qual_swir_id, observation_grp_id;
    std::vector<u_char> qual_red, qual_blue, qual_SWIR;
    size_t n_blue_bands, n_red_bands, n_swir_bands;
    size_t n_pixels, n_lines, n_bands;
    file_format format;  // FT_OCIL1B
    int ncid;
    int32_t max_blue_index, max_red_index;
    std::vector<u_char> quality_flags;
    std::vector<u_char> quality_flags_blue, quality_flags_red, quality_flags_swir;

   public:
    [[nodiscard]] const l1str& get_l1rec() const {
        return l1rec;
    }
    [[nodiscard]] const filehandle& get_l1file() const {
        return l1file;
    }
    [[nodiscard]] const file_format get_format() const {
        return format;
    }
    friend int bin_l1b_file(L1BFile& l1bfile, L1CFile& l1cfile);
    L1BFile() = delete;
    explicit L1BFile(const std::string& ifile, std::vector<int>&& bin_indexes,
                     std::vector<std::vector<std::pair<int, double>>>&& area_weights,
                     std::vector<float>&& height_data_corrected, std::vector<float>&& senz_corrected,
                     std::vector<float>&& sena_corrected);
    explicit L1BFile(const std::string& ifile, std::vector<int>&& bin_indexes);
    int read_line(int line);
    ~L1BFile() {
        closel1(&l1file);
        free_l1(&l1rec);
    }
};
#endif  // L1CFILE_H
