#ifndef L1C_INPUT_STR_H
#define L1C_INPUT_STR_H
#include "clo.h"
#include <string>
#include "genutils.h"
#include "timeutils.h"
#include <netcdf>
#include <algorithm>
#include "L1BFile.h"
#include "L1CFile.h"
#include "search_l1c.h"
#include <scaled_nc_var.hpp>
#include <memory>
#include "utils.h"
#include <optional>
#include <boost/algorithm/string.hpp>

class L1CProcessor {
    std::vector<std::string> input_l1bs{};
    std::vector<std::string> input_l2_anc{};
    std::string l1_grid_file{};
    std::string ofile{};
    std::string dem_file{};
    std::vector<std::tuple<std::string, double, double>> l1b_unix_time_start_end{};
    std::vector<std::tuple<std::string, double, double>> l2_unix_time_start_end{};
    std::unique_ptr<L1CFile> l1c_ofile;
    std::optional<float> cloud_top_height;
    bool verbose = false;
    int area_weighting = 1;
    bool cloud_correct = false;
    InputAttributes input_attributes;
    bool parse_clo(clo_optionList_t* list);

   public:
    L1CProcessor() = default;
    explicit L1CProcessor(const std::string& program_name, const std::string& version,
                          const std::vector<std::string>& argv);
    ~L1CProcessor() = default;
    explicit L1CProcessor(clo_optionList_t* list);

    bool load_input(clo_optionList_t* list);
    std::vector<short> read_dem(std::vector<float>& lat_data, std::vector<float>& lon_data, size_t n_lines, size_t n_pixels);
    std::tuple<std::vector<float>, std::vector<int8_t>,std::vector<float>,std::vector<float>> get_cth(size_t l2_file_count,size_t & n_lines_l2, size_t & n_pixels_l2);
    bool validate_l1c_geolocation(const std::vector<float> & lat, const std::vector<float> & lon, size_t npixels);
    void binL1Bgranules();
};

bool l1c_input(int argc, char** argv, L1CProcessor& input, const std::string& program_name,
               const std::string& version);
#endif  // L1C_INPUT_STR_H
