#include "l1c_processor.h"
#include "area_weighting.h"
#include "surface_correction.h"
#include <boost/algorithm/string/trim.hpp>
#include "utils.h"
/*
 * @brief Constructor for L1CProcessor
 * @param program_name - name of the program
 * @param version - version of the program
 * @param argv - command line arguments
 */
L1CProcessor::L1CProcessor(const std::string& program_name, const std::string& version,
                           const std::vector<std::string>& argv) {
    std::string history{};
    for (const auto& arg : argv) {
        history += arg + " ";
    }
    boost::trim(history);
    input_attributes["history"] = history;
    input_attributes["software_name"] = program_name;
    input_attributes["software_version"] = version;
};

L1CProcessor::L1CProcessor(clo_optionList_t* list) {
    load_input(list);
}

/**
 * @brief Parse command line options and populate input attributes.
 * @param list - command line options
 * @return true if parsing is successful, false otherwise
 */
bool L1CProcessor::parse_clo(clo_optionList_t* list) {
    const int numOptions = clo_getNumOptions(list);
    for (int optionId = 0; optionId < numOptions; optionId++) {
        clo_option_t* option = clo_getOption(list, optionId);
        // ignore options of type CLO_TYPE_HELP
        if (option->dataType == CLO_TYPE_HELP)
            continue;
        std::string keyword = option->key;
        input_attributes[keyword] = "";
        if (keyword == "ifile") {
            if (clo_isOptionSet(option)) {
                std::string parm_str = clo_getOptionString(option);
                std::vector<std::string> input_files = readFileList(parm_str);
                // test first file for format
                StdoutRedirector redirector;
                redirector.redirect();
                file_format format = getFormat(const_cast<char*>(input_files[0].c_str()));
                redirector.restore();
                std::string input_str = "";
                if (format.type == FT_OCIL1B) {
                    input_l1bs = input_files;
                    input_str = boost::algorithm::join(input_l1bs, ",");
                } else {
                    fprintf(stderr, "-E-: %s:%d input file format not recognized\n", __FILE__, __LINE__);
                    return false;
                }
                input_attributes["ifile"] = input_str;
            }
        }
        if (keyword == "ofile") {
            if (clo_isOptionSet(option)) {
                ofile = clo_getOptionString(option);
                input_attributes["ofile"] = ofile;
            }
        }
        if (keyword == "l1c_grid") {
            if (clo_isOptionSet(option)) {
                l1_grid_file = clo_getOptionString(option);
                input_attributes["l1c_grid"] = l1_grid_file;
            }
        }
        if (keyword == "par") {
            if (clo_isOptionSet(option)) {
                input_attributes["par"] = clo_getOptionString(option);
            }
        }
        if (keyword == "doi") {
            if (clo_isOptionSet(option)) {
                input_attributes["doi"] = clo_getOptionString(option);
            }
        }
        if (keyword == "pversion") {
            if (clo_isOptionSet(option)) {
                input_attributes["pversion"] = clo_getOptionString(option);
            }
        }
        if (keyword == "cloud_height") {
            if (clo_isOptionSet(option)) {
                cloud_top_height = clo_getOptionFloat(option);
                input_attributes["cloud_height"] = cloud_top_height.value();
                cloud_correct = true;
            }
        }
        if (keyword == "cloud_anc_files") {
            if (clo_isOptionSet(option)) {
                std::string parm_str = clo_getOptionString(option);
                input_l2_anc = readFileList(parm_str);
                // input string of files
                std::string cloud_anc_files_str = boost::algorithm::join(input_l2_anc, ",");
                input_attributes["cloud_anc_files"] = cloud_anc_files_str;
                cloud_correct = true;
            }
        }
        if (keyword == "verbose") {
            if (clo_isOptionSet(option)) {
                verbose = clo_getOptionBool(option);
                input_attributes["verbose"] = verbose;
            }
        }
        if (keyword == "area_weighting") {
            area_weighting = clo_getOptionInt(option);
            input_attributes["area_weighting"] = area_weighting;
        }
        if (keyword == "demfile") {
            dem_file = expand_env_variable_path(clo_getOptionString(option));
            input_attributes["demfile"] = dem_file;
        }
    }
    if (cloud_correct && input_l2_anc.empty()) {
        fprintf(stderr, "-E-: %s:%d Cloud correction requested but no L2 ancillary files provided\n",
                __FILE__, __LINE__);
        return false;
    }
    if (ofile.empty()) {
        fprintf(stderr, "-E-: %s:%d No output file specified\n", __FILE__, __LINE__);
        return false;
    }
    if (!input_l1bs.empty() && l1_grid_file.empty()) {
        fprintf(stderr, "-E-: %s:%d No L1C grid file specified\n", __FILE__, __LINE__);
        return false;
    }
    if (input_l1bs.empty()) {
        fprintf(stderr, "-E-: %s:%d No input L1B or L1C files specified\n", __FILE__, __LINE__);
        return false;
    }
    // using std::filessystem to check if input files exist
    for (const auto& file : input_l1bs) {
        if (!std::filesystem::exists(file)) {
            fprintf(stderr, "-E-: %s:%d Input file '%s' does not exist\n", __FILE__, __LINE__, file.c_str());
            return false;
        }
    }
    if (!l1_grid_file.empty() && !std::filesystem::exists(l1_grid_file)) {
        fprintf(stderr, "-E-: %s:%d L1C grid file '%s' does not exist\n", __FILE__, __LINE__,
                l1_grid_file.c_str());
        return false;
    }
    for (const auto& file : input_l2_anc) {
        if (!std::filesystem::exists(file)) {
            fprintf(stderr, "-E-: %s:%d L2 Ancillary file '%s' does not exist\n", __FILE__, __LINE__,
                    file.c_str());
            return false;
        }
    }
    return true;
}
/**
 * @brief Load input files and check consistency of time coverage between L1B and L2 ancillary files.
 * @param list - command line options
 * @return true if input is valid, false otherwise
 */
bool L1CProcessor::load_input(clo_optionList_t* list) {
    // parse command line options
    if (!parse_clo(list))
        return false;
    // read time coverage from L1B files
    if (!input_l1bs.empty()) {
        l1b_unix_time_start_end = get_time_coverage(input_l1bs);
    }
    if (!input_l2_anc.empty()) {
        l2_unix_time_start_end = get_time_coverage(input_l2_anc);
    }
    if (!input_l1bs.empty() && !input_l2_anc.empty()) {
        size_t number_of_files = input_l1bs.size();
        if (input_l2_anc.size() != number_of_files) {
            fprintf(
                stderr,
                "-E-: %s:%d Number of L1B files (%lu) does not match number of L2 Ancillary files (%lu)\n",
                __FILE__, __LINE__, input_l1bs.size(), input_l2_anc.size());
            return false;
        }
        // check exact match of time coverage within a tolerance of 0.5 second
        const double time_tolerance = 0.5;  // seconds
        for (size_t i = 0; i < number_of_files; ++i) {
            const auto& [l1b_file, l1b_start, l1b_end] = l1b_unix_time_start_end[i];
            const auto& [l2_file, l2_start, l2_end] = l2_unix_time_start_end[i];
            if (std::abs(l1b_start - l2_start) > time_tolerance ||
                std::abs(l1b_end - l2_end) > time_tolerance) {
                fprintf(stderr,
                        "-E-: %s:%d Time coverage of L1B file '%s' (%.3f - %.3f) does not match time "
                        "coverage of L2 Ancillary file '%s' (%.3f - %.3f) within tolerance of %.3f seconds\n",
                        __FILE__, __LINE__, l1b_file.c_str(), l1b_start, l1b_end, l2_file.c_str(), l2_start,
                        l2_end, time_tolerance);
                return false;
            }
        }
    }

    return true;
}
/**
 * @brief Bin input L1B granules to L1C grid.
 */
void L1CProcessor::binL1Bgranules() {

    netCDF::NcFile nc_l1c(l1_grid_file, netCDF::NcFile::read);
    if (nc_l1c.isNull()) {
        fprintf(stderr, "-E-: %s:%d Failed to open L1C grid file %s\n", __FILE__, __LINE__,
                l1_grid_file.c_str());
        exit(EXIT_FAILURE);
    }
    ScaledNcVar lat_var = nc_l1c.getGroup("geolocation_data").getVar("latitude");
    ScaledNcVar lon_var = nc_l1c.getGroup("geolocation_data").getVar("longitude");
    ScaledNcVar nadir_view_time_var = nc_l1c.getGroup("bin_attributes").getVar("nadir_view_time");
    // check if variables are null
    if (lat_var.isNull() || lon_var.isNull() || nadir_view_time_var.isNull()) {
        fprintf(stderr, "-E-: %s:%d the file %s doesn't contain required variables", __FILE__, __LINE__,
                l1_grid_file.c_str());
        exit(EXIT_FAILURE);
    }
    auto lat_dims = lat_var.getDims();
    size_t n_lines = lat_dims[0].getSize();
    size_t n_pixels = lat_dims[1].getSize();
    std::vector<float> lat_data(n_lines * n_pixels);
    std::vector<float> lon_data(n_lines * n_pixels);
    std::vector<short> height_data(n_lines * n_pixels);
    std::vector<double> nadir_view_time_data(n_lines);
    NCPP_ERROR(nadir_view_time_var.getVar(nadir_view_time_data.data()));
    NCPP_ERROR(lat_var.getVar(lat_data.data()));
    NCPP_ERROR(lon_var.getVar(lon_data.data()));

    if (!validate_l1c_geolocation(lat_data, lon_data, n_pixels)) {
        fprintf(stderr, "-E-:%s:%d Invalid geolocation found in l1c grid input file %s ", __FILE__, __LINE__,
                l1_grid_file.c_str());
        exit(EXIT_FAILURE);
    }

    if (dem_file.empty()) {
        std::cout << "Reading l1c grid file " << l1_grid_file << " for height data..." << std::endl;
        ScaledNcVar height_var = nc_l1c.getGroup("geolocation_data").getVar("height");
        NCPP_ERROR(height_var.getVar(height_data.data()));
    } else {
        std::cout << "Reading DEM file " << dem_file << " for height data..." << std::endl;
        height_data = read_dem(lat_data, lon_data, n_lines, n_pixels);
    }

    SearchL1C searcher(lat_data, lon_data, n_lines, n_pixels);
    size_t file_count = 0;
    for (const auto& [file_l1b, time_st, time_end] : l1b_unix_time_start_end) {
        netCDF::NcFile nc_l1b(file_l1b, netCDF::NcFile::read);
        if (nc_l1b.isNull()) {
            fprintf(stderr, "-E-: %s:%d Failed to open L1B file %s\n", __FILE__, __LINE__, file_l1b.c_str());
            exit(EXIT_FAILURE);
        }
        ScaledNcVar lat_var_b = nc_l1b.getGroup("geolocation_data").getVar("latitude");
        ScaledNcVar lon_var_b = nc_l1b.getGroup("geolocation_data").getVar("longitude");
        if (lat_var_b.isNull() || lon_var_b.isNull()) {
            fprintf(stderr, "-E-: %s:%d the file %s doesn't contain latitude/longitude", __FILE__, __LINE__,
                    file_l1b.c_str());
            exit(EXIT_FAILURE);
        }
        auto lat_dims_b = lat_var_b.getDims();
        size_t n_lines_b = lat_dims_b[0].getSize();
        size_t n_pixels_b = lat_dims_b[1].getSize();
        std::vector<float> lat_data_b(n_lines_b * n_pixels_b);
        std::vector<float> lon_data_b(n_lines_b * n_pixels_b);
        NCPP_ERROR(lat_var_b.getVar(lat_data_b.data()));
        NCPP_ERROR(lon_var_b.getVar(lon_data_b.data()));
        std::vector<float> height_data_b{};
        std::vector<float> sensor_zenith_angle_data_b{};
        std::vector<float> sensor_azimuth_angle_data_b{};
        // here we need do cloud top height correction if requested.
        if (cloud_correct) {
            size_t n_lines_l2, n_pixels_l2;
            const auto [cth_data, cloud_mask, lat_data_l2, lon_data_l2] =
                get_cth(file_count, n_lines_l2, n_pixels_l2);
            ScaledNcVar height_b = nc_l1b.getGroup("geolocation_data").getVar("height");
            if (height_b.isNull()) {
                fprintf(stderr, "-E-: %s:%d the file %s doesn't contain height", __FILE__, __LINE__,
                        file_l1b.c_str());
                exit(EXIT_FAILURE);
            } else {
                height_data_b.resize(n_lines_b * n_pixels_b);
                std::vector<short> height_meteres(n_lines_b * n_pixels_b);
                NCPP_ERROR(height_b.getVar(height_meteres.data()));
                for (size_t ip = 0; ip < n_lines_b * n_pixels_b; ip++) {
                    height_data_b[ip] = static_cast<float>(height_meteres[ip]) / 1000.0;  // to km
                }
            }
            ScaledNcVar sensor_zenith_angle_b = nc_l1b.getGroup("geolocation_data").getVar("sensor_zenith");
            sensor_zenith_angle_data_b.resize(n_lines_b * n_pixels_b);
            ScaledNcVar sensor_azimuth_angle_b = nc_l1b.getGroup("geolocation_data").getVar("sensor_azimuth");
            sensor_azimuth_angle_data_b.resize(n_lines_b * n_pixels_b);
            if (sensor_zenith_angle_b.isNull() || sensor_azimuth_angle_b.isNull()) {
                fprintf(stderr,
                        "-E-: %s:%d the file %s doesn't contain sensor_azimuth_angle/sensor_zenith_angle",
                        __FILE__, __LINE__, file_l1b.c_str());
                exit(EXIT_FAILURE);
            } else {
                NCPP_ERROR(sensor_azimuth_angle_b.getVar(sensor_azimuth_angle_data_b.data()));
                NCPP_ERROR(sensor_zenith_angle_b.getVar(sensor_zenith_angle_data_b.data()));
            }
            std::cout << "Performing cloud top height correction for L1B file " << file_l1b << " ..."
                      << std::endl;
            srf::surface_correction(lat_data_b, lon_data_b, height_data_b, sensor_zenith_angle_data_b,
                                    sensor_azimuth_angle_data_b, lat_data_l2, lon_data_l2, cth_data,
                                    cloud_mask, n_lines_b, n_pixels_b, n_lines_l2, n_pixels_l2);
            file_count++;
        };
        nc_l1b.close();
        // get bin indexes for the l1c grid
        std::vector<int> bin_indexes = searcher.query_lat_lon(lat_data_b, lon_data_b);
        std::vector<std::vector<std::pair<int, double>>> area_weights{};
        if (area_weighting) {
            std::cout << "Calculating area weights for " << file_l1b << std::endl;
            area_weights = get_area_weights(lat_data, lon_data, lat_data_b, lon_data_b, n_lines, n_pixels,
                                            n_lines_b, n_pixels_b, bin_indexes, verbose);
        }
        // debug test
        if (verbose && area_weighting) {
            size_t number_of_valid = 0;
            for (const auto& aw : area_weights) {
                number_of_valid += aw.size() > 0;
            }
            std::cout << "Number of valid indices from area_weights: " << number_of_valid << std::endl;
            size_t count_data = 0;
            for (size_t ib = 0; ib < area_weights.size(); ib++) {
                const auto& area_weight = area_weights[ib];
                if (area_weight.empty())
                    continue;
                count_data++;
                std::cout << "OCI L1B pixel " << ib << " : area = " << get_l1b_area(ib) << std::endl;
                for (const auto& [index_l1c, area_l1c] : area_weight) {
                    std::cout << "L1C " << index_l1c << " ; Area =  " << area_l1c << std::endl;
                }
                if (count_data > 20)
                    break;
            }
        }
        if (verbose) {
            size_t number_of_valid = count_valid(bin_indexes);
            std::cout << "Number of valid indices from bin_indexes: " << number_of_valid << std::endl;
            // for first  10 indices, which are not -1, print the index, and haversine distance
            size_t printed = 0;
            for (size_t i = 0; i < bin_indexes.size() && printed < 10; ++i) {
                if (bin_indexes[i] != -1) {
                    int index_l1c = bin_indexes[i];
                    size_t line = index_l1c / n_pixels;
                    size_t pixel = index_l1c % n_pixels;
                    double dist = haversineDistance(lat_data_b[i], lon_data_b[i], lat_data[index_l1c],
                                                    lon_data[index_l1c]);
                    std::cout << "Index: " << index_l1c << ", Line: " << line << ", Pixel: " << pixel
                              << ", Haversine Distance: " << dist << " km" << std::endl;
                    printed++;
                }
            }
        }
        L1BFile l1bfile(file_l1b, std::move(bin_indexes), std::move(area_weights), std::move(height_data_b),
                        std::move(sensor_zenith_angle_data_b), std::move(sensor_azimuth_angle_data_b));
        if (!l1c_ofile) {
            std::multimap<std::string, netCDF::NcGroupAtt> global_attributes = nc_l1c.getAtts();
            l1c_ofile =
                std::make_unique<L1CFile>(ofile, global_attributes, l1bfile, n_lines, n_pixels, lat_data,
                                          lon_data, height_data, nadir_view_time_data, input_attributes);
            nc_l1c.close();
        }
        std::cout << "Binning L1B file " << file_l1b << " ..." << std::endl;
        int status = bin_l1b_file(l1bfile, *l1c_ofile);
        if (status != 0) {
            std::cerr << "Error binning L1B file " << file_l1b << std::endl;
            exit(EXIT_FAILURE);
        }
    }
    l1c_ofile->close();
}
/*
 * @brief Get time coverage information from a list of NetCDF files.
 * @param files - vector of NetCDF file paths
 * @return vector of tuples containing file name, start time, and end time
 */
std::vector<std::tuple<std::string, double, double>> get_time_coverage(
    const std::vector<std::string>& files) {
    std::vector<std::tuple<std::string, double, double>> time_coverage;
    for (const auto& file : files) {
        netCDF::NcFile nc_file(file, netCDF::NcFile::read);
        if (nc_file.isNull()) {
            fprintf(stderr, "-E-: %s:%d opening file '%s'\n", __FILE__, __LINE__, file.c_str());
            exit(EXIT_FAILURE);
        }
        std::string time_coverage_start, time_coverage_end;
        netCDF::NcGroupAtt time_coverage_start_att = nc_file.getAtt("time_coverage_start");
        netCDF::NcGroupAtt time_coverage_end_att = nc_file.getAtt("time_coverage_start");
        if (time_coverage_end_att.isNull() || time_coverage_start_att.isNull()) {
            fprintf(stderr, "-E-: %s:%d opening file '%s'\n", __FILE__, __LINE__, file.c_str());
            exit(EXIT_FAILURE);
        }
        try {
            time_coverage_start_att.getValues(time_coverage_start);
            time_coverage_end_att.getValues(time_coverage_end);
        } catch (const netCDF::exceptions::NcException& e) {
            fprintf(stderr, "-E-: %s:%d getting attributes from  file '%s'. %s\n", __FILE__, __LINE__,
                    file.c_str(), e.what());
            exit(EXIT_FAILURE);
        }
        double start_time, end_time;
        start_time = isodate2unix(time_coverage_start.c_str());
        end_time = isodate2unix(time_coverage_end.c_str());
        time_coverage.emplace_back(std::make_tuple(file, start_time, end_time));
    }
    // sort by start time
    std::sort(
        time_coverage.begin(), time_coverage.end(),
        [](const std::tuple<std::string, double, double>& a,
           const std::tuple<std::string, double, double>& b) { return std::get<1>(a) < std::get<1>(b); });
    return time_coverage;
}

std::tuple<std::vector<float>, std::vector<int8_t>, std::vector<float>, std::vector<float>>
L1CProcessor::get_cth(size_t l2_file_count, size_t& n_lines_l2, size_t& n_pixels_l2) {
    const auto& [l2_file, l2_start, l2_end] = l2_unix_time_start_end[l2_file_count];
    std::cout << "Reading cloud top height data from L2 ancillary file " << l2_file << std::endl;

    netCDF::NcFile nc_l2(l2_file, netCDF::NcFile::read);
    if (nc_l2.isNull()) {
        fprintf(stderr, "-E-: %s:%d Failed to open L2 ancillary cloud file %s\n", __FILE__, __LINE__,
                l2_file.c_str());
        exit(EXIT_FAILURE);
    }
    ScaledNcVar lat_var_l2 = nc_l2.getGroup("navigation_data").getVar("latitude");
    ScaledNcVar lon_var_l2 = nc_l2.getGroup("navigation_data").getVar("longitude");
    ScaledNcVar l2_flags = nc_l2.getGroup("geophysical_data").getVar("l2_flags");
    // check if variables are null
    if (lat_var_l2.isNull()) {
        fprintf(stderr, "-E-: %s:%d Latitude variable not found in L2 ancillary file '%s'\n", __FILE__,
                __LINE__, l2_file.c_str());
        exit(EXIT_FAILURE);
    }
    if (lon_var_l2.isNull()) {
        fprintf(stderr, "-E-: %s:%d Longitude variable not found in L2 ancillary file '%s'\n", __FILE__,
                __LINE__, l2_file.c_str());
        exit(EXIT_FAILURE);
    }
    if (l2_flags.isNull()) {
        fprintf(stderr, "-E-: %s:%d l2_flags variable not found in L2 ancillary file '%s'\n", __FILE__,
                __LINE__, l2_file.c_str());
        exit(EXIT_FAILURE);
    }
    auto lat_dims_l2 = lat_var_l2.getDims();
    n_lines_l2 = lat_dims_l2[0].getSize();
    n_pixels_l2 = lat_dims_l2[1].getSize();
    if (verbose)
        std::cout << "L2 ancillary file " << l2_file << " has dimensions: " << n_lines_l2 << " x "
                  << n_pixels_l2 << std::endl;
    std::vector<float> lat_data_l2(n_lines_l2 * n_pixels_l2);
    std::vector<float> lon_data_l2(n_lines_l2 * n_pixels_l2);
    std::vector<int> l2_flags_data(n_lines_l2 * n_pixels_l2);
    NCPP_ERROR(lat_var_l2.getVar(lat_data_l2.data()));
    NCPP_ERROR(lon_var_l2.getVar(lon_data_l2.data()));
    NCPP_ERROR(l2_flags.getVar(l2_flags_data.data()));

    // calculate cloud mask using CLOUD flag bit and cloud top height threshold
    std::vector<int8_t> cloud_mask(n_lines_l2 * n_pixels_l2, 0);
    int cloud_flag_bit = 0;
    // read flag_masks and flag_meanings variable attributes to get the correct cloud_flag_bit
    netCDF::NcVarAtt flag_masks_att = l2_flags.getAtt("flag_masks");
    netCDF::NcVarAtt flag_meanings_att = l2_flags.getAtt("flag_meanings");
    // find the index of the CLOUD flag in flag_meanings
    // use boost to split string flag_meanings by space and find the index of CLOUD
    std::string flag_meanings;
    std::vector<std::string> flag_meanings_vec;
    NCPP_ERROR(flag_meanings_att.getValues(flag_meanings));
    std::vector<std::string> flags = boost::split(flag_meanings_vec, flag_meanings, boost::is_any_of(" "));
    std::vector<int> flag_masks_values(flags.size());
    NCPP_ERROR(flag_masks_att.getValues(flag_masks_values.data()));
    for (size_t i = 0; i < flags.size(); ++i) {
        if (flags[i] == "CLOUD" || flags[i] == "CLDICE") {
            cloud_flag_bit = flag_masks_values[i];
            break;
        }
    }
    if (verbose)
        std::cout << "Cloud flag bit: " << cloud_flag_bit << std::endl;
    if (cloud_flag_bit == 0) {
        fprintf(stderr, "-E-: %s:%d CLOUD/CLDICE flag not found in l2_flags of file '%s'\n", __FILE__,
                __LINE__, l2_file.c_str());
        exit(EXIT_FAILURE);
    }
    size_t cloud_pixel_count = 0;
    for (size_t ip = 0; ip < n_lines_l2 * n_pixels_l2; ip++) {
        if ((l2_flags_data[ip] & cloud_flag_bit) != 0) {
            cloud_mask[ip] = 1;
            cloud_pixel_count++;
        }
    }
    if (verbose)
        std::cout << "Number of cloudy pixels in L2 ancillary file " << l2_file << " : " << cloud_pixel_count
                  << std::endl;
    std::vector<float> cth_data(n_lines_l2 * n_pixels_l2);
    if (!cloud_top_height.has_value()) {
        ScaledNcVar cth = nc_l2.getGroup("geophysical_data").getVar("cth");
        if (cth.isNull()) {
            fprintf(stderr, "-E-: %s:%d Cloud top height variable 'cth' not found in file '%s'\n", __FILE__,
                    __LINE__, l2_file.c_str());
            exit(EXIT_FAILURE);
        }
        NCPP_ERROR(cth.getVar(cth_data.data()));
        for (size_t ip = 0; ip < n_lines_l2 * n_pixels_l2; ip++) {
            if (cth_data[ip] == BAD_FLT)
                cloud_mask[ip] = 0;
        }
    } else {
        std::fill(cth_data.begin(), cth_data.end(), cloud_top_height.value());
    }
    nc_l2.close();
    if (verbose)
        std::cout << "Finished reading cloud top height data from L2 ancillary file " << l2_file << std::endl;
    return {cth_data, cloud_mask, lat_data_l2, lon_data_l2};
}

/*
 * @brief Read DEM GEBCO data for given lat/lon coordinates and return surface height.
 * @param lat_data - vector of latitude values
 * @param lon_data - vector of longitude values
 * @param n_lines - number of lines in the L1C grid
 * @param n_pixels - number of pixels in the L1C grid
 * @return vector of surface height values corresponding to the input lat/lon coordinates, using nearest neighbor interpolation from the DEM data
 */
std::vector<short> L1CProcessor::read_dem(std::vector<float>& lat_data, std::vector<float>& lon_data,
                                          size_t n_lines, size_t n_pixels) {
    netCDF::NcFile nc_dem(dem_file, netCDF::NcFile::read);
    if (nc_dem.isNull()) {
        fprintf(stderr, "-E-: %s:%d Failed to open DEM file %s\n", __FILE__, __LINE__, dem_file.c_str());
        exit(EXIT_FAILURE);
    }
    // read lat / lon double, 1D variables
    ScaledNcVar lat_var = nc_dem.getVar("lat");
    ScaledNcVar lon_var = nc_dem.getVar("lon");
    if (lat_var.isNull() || lon_var.isNull()) {
        fprintf(stderr, "-E-: %s:%d the DEM file %s doesn't contain lat/lon variables\n", __FILE__, __LINE__,
                dem_file.c_str());
        exit(EXIT_FAILURE);
    }
    auto lat_dims = lat_var.getDims();
    auto lon_dims = lon_var.getDims();
    if (lat_dims.size() != 1 || lon_dims.size() != 1) {
        fprintf(stderr, "-E-: %s:%d the DEM file %s lat/lon variables are not 1D\n", __FILE__, __LINE__,
                dem_file.c_str());
        exit(EXIT_FAILURE);
    }

    size_t n_lat = lat_dims[0].getSize();
    size_t n_lon = lon_dims[0].getSize();
    std::vector<double> lat_data_dem(n_lat);
    std::vector<double> lon_data_dem(n_lon);
    NCPP_ERROR(lat_var.getVar(lat_data_dem.data()));
    NCPP_ERROR(lon_var.getVar(lon_data_dem.data()));
    double step_lat = lat_data_dem[1] - lat_data_dem[0];
    double step_lon = lon_data_dem[1] - lon_data_dem[0];
    ScaledNcVar water_surface_height_var = nc_dem.getVar("water_surface_height");
    ScaledNcVar height_var = nc_dem.getVar("height");
    ScaledNcVar watermask_var = nc_dem.getVar("watermask");
    if (water_surface_height_var.isNull()) {
        fprintf(stderr, "-E-: %s:%d the DEM file %s doesn't contain water_surface_height variable\n",
                __FILE__, __LINE__, dem_file.c_str());
        exit(EXIT_FAILURE);
    }
    if (height_var.isNull()) {
        fprintf(stderr, "-E-: %s:%d the DEM file %s doesn't contain height variable\n", __FILE__, __LINE__,
                dem_file.c_str());
        exit(EXIT_FAILURE);
    }
    if (watermask_var.isNull()) {
        fprintf(stderr, "-E-: %s:%d the DEM file %s doesn't contain watermask variable\n", __FILE__, __LINE__,
                dem_file.c_str());
        exit(EXIT_FAILURE);
    }
    std::vector<short> surface_height_data(n_lines * n_pixels);
    for (size_t i = 0; i < n_lines * n_pixels; i++) {
        float lat_l1b = lat_data[i];
        float lon_l1b = lon_data[i];
        size_t ilat = static_cast<size_t>(std::round((lat_l1b - lat_data_dem[0]) / step_lat));
        size_t ilon = static_cast<size_t>(std::round((lon_l1b - lon_data_dem[0]) / step_lon));
        int8_t water_flag;
        // handle DEM file wrap around the anitmeridian or poles
        if (ilat >= n_lat)
            ilat = 0;
        if (ilon >= n_lon)
            ilon = 0;
        watermask_var.getVar({ilat, ilon}, {1, 1}, &water_flag);
        if (water_flag)
            water_surface_height_var.getVar({ilat, ilon}, {1, 1}, &surface_height_data[i]);
        else
            height_var.getVar({ilat, ilon}, {1, 1}, &surface_height_data[i]);
    }
    nc_dem.close();
    return surface_height_data;
}
/*
 * @brief Validate geolocation coordinates for L1C products.
 * @param lat - vector of latitude values
 * @param lon - vector of longitude values
 * @param npixels - number of pixels in the product
 * @return true if all coordinates are valid, false otherwise
 */
bool L1CProcessor::validate_l1c_geolocation(const std::vector<float>& lat, const std::vector<float>& lon,
                                            size_t npixels) {
    size_t nelems = lat.size();
    for (size_t ip = 0; ip < nelems; ip++) {
        if (!valid_coordinates(lat[ip], lon[ip])) {
            printf("Invalid geolocation found at line %ld, pixel %ld\n", ip / npixels, ip % npixels);
            return false;
        }
    }
    return true;
}