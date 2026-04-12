#include "L1CFile.h"
#include "time.h"
/**
 * @brief Constructs an L1CF file and initializes a netCDF output file for binned L1C data.
 *
 *
 * @param[in] output_file_name Path to the output netCDF file to be created.
 * @param[in] l1c_global_attributes Map of global attributes from the user input
 * @param[in] l1bfile Reference to the input L1BFile file containing source data and metadata.
 * @param[in] n_lines Number of bins along track dimension. (bins_along_track)
 * @param[in] n_pixels Number of bins across track dimension. (bins_across_track)
 * @param[in] latitude Vector of latitude values for bin grid (bins_along_track x bins_across_track).
 * @param[in] longitude Vector of longitude values for bin grid (bins_along_track x bins_across_track).
 * @param[in] height Vector of surface height values for bin grid (bins_along_track x bins_across_track).
 * @param[in] nadir_view_time Vector of nadir view times for each line (bins_along_track).
 * @param[in] input_attributes Map of processing parameters including area weighting mode,
 *                             cloud correction settings, product name, version, and DOI.
 */
L1CFile::L1CFile(const std::string& output_file_name,
                 const std::multimap<std::string, netCDF::NcGroupAtt>& l1c_global_attributes,
                 const L1BFile& l1bfile, size_t n_lines, size_t n_pixels, const std::vector<float>& latitude,
                 const std::vector<float>& longitude, const std::vector<short>& height,
                 const std::vector<double>& nadir_view_time, const InputAttributes& input_attributes)
    : n_lines(n_lines),
      n_pixels(n_pixels),
      latitude_data(latitude),
      longitude_data(longitude),
      height_data(height),
      nadir_view_time_data(nadir_view_time),
      o_file(output_file_name) {
    if (std::filesystem::exists(output_file_name))
        std::remove(output_file_name.c_str());
    NCPP_ERROR(file.open(output_file_name, netCDF::NcFile::newFile))
    n_bands = l1bfile.get_l1file().nbands;
    std::vector<size_t> chunking = {32, n_pixels, 1, n_bands};
    netCDF::NcDim bins_along_track_dim = file.addDim("bins_along_track", n_lines);
    netCDF::NcDim bins_across_track_dim = file.addDim("bins_across_track", n_pixels);
    netCDF::NcDim intensity_bands_per_view_dim = file.addDim("intensity_bands_per_view", n_bands);
    netCDF::NcDim number_of_views_dim = file.addDim("number_of_views", 1);
    netCDF::NcGroup geolocation_data_group = file.addGroup("geolocation_data");
    netCDF::NcGroup observation_data_group = file.addGroup("observation_data");
    netCDF::NcGroup bin_attributes_group = file.addGroup("bin_attributes");
    netCDF::NcGroup sensor_views_bands_group = file.addGroup("sensor_views_bands");
    netCDF::NcGroup processing_control_group = file.addGroup("processing_control");
    netCDF::NcGroup input_parameters_group = processing_control_group.addGroup("input_parameters");
    i_data = std::vector<float>(n_bands * n_lines * n_pixels, BAD_FLT);
    i_stdev_data = std::vector<float>(n_bands * n_lines * n_pixels, BAD_FLT);
    number_of_observations = std::vector<short>(n_lines * n_pixels, 0);
    number_of_observations_band = std::vector<short>(n_bands * n_lines * n_pixels, 0);
    // allocated angles
    sensor_zenith_angles_data = std::vector<float>(n_lines * n_pixels, BAD_FLT);
    solar_zenith_angles_data = std::vector<float>(n_lines * n_pixels, BAD_FLT);
    sensor_azimuth_angles_data = std::vector<float>(n_lines * n_pixels, BAD_FLT);
    solar_azimuth_angles_data = std::vector<float>(n_lines * n_pixels, BAD_FLT);
    scattering_angles_data = std::vector<float>(n_lines * n_pixels, BAD_FLT);
    height_stdev = std::vector<float>(n_lines * n_pixels, BAD_FLT);
    // time data:
    view_time_offsets_data = std::vector<double>(n_lines * n_pixels, BAD_FLT);

    int area_weighting = std::get<int>(input_attributes.at("area_weighting"));
    if (area_weighting) {
        area_total = std::vector<float>(n_lines * n_pixels, 0);
        area_total_per_band = std::vector<float>(n_bands * n_lines * n_pixels, 0);
        area_second_momentum = std::vector<float>(n_bands * n_lines * n_pixels, 0);
    }
    void* attr_storage = malloc(1024);
    std::string time_coverage_start;
    for (const auto& [name, att] : l1c_global_attributes) {
        if (name == "history" || name == "product_name" || name == "date_created") {
            continue;
        }
        NCPP_ERROR(att.getValues(attr_storage));
        netCDF::NcType type = att.getType();
        size_t len = att.getAttLength();
        NCPP_ERROR((void)file.putAtt(name, type, len, attr_storage));
        if (name == "time_coverage_start")
            att.getValues(time_coverage_start);
    }
    // get date_created

    char time_string[512];                     //
    time_t now = time(NULL);                  //
    struct tm* local_time = localtime(&now);  // Convert to local time structure

    // Format the time into the string buffer (e.g., "YYYY-MM-DD HH:MM:SS")
    strftime(time_string, sizeof(time_string), "%Y-%m-%dT%H:%M:%SZ", local_time);
    // put date_created attribute
    NCPP_ERROR((void)file.putAtt("date_created", time_string))
    // set history and product_name, set processing control
    Visitor visitor;
    std::string history = std::visit(visitor, input_attributes.at("history"));
    std::string product_name = std::visit(visitor, input_attributes.at("ofile"));
    std::string doi = std::visit(visitor, input_attributes.at("doi"));
    std::string processing_version = std::visit(visitor, input_attributes.at("pversion"));
    NCPP_ERROR((void)file.putAtt("history", history));
    NCPP_ERROR((void)file.putAtt("product_name", product_name));
    NCPP_ERROR((void)file.putAtt("doi", doi));
    NCPP_ERROR((void)file.putAtt("processing_version", processing_version));
    std::string software_name = std::visit(visitor, input_attributes.at("software_name"));
    std::string software_version = std::visit(visitor, input_attributes.at("software_version"));

    NCPP_ERROR((void)processing_control_group.putAtt("software_name", software_name));
    NCPP_ERROR((void)processing_control_group.putAtt("software_version", software_version));

    for (const auto& [name, att] : input_attributes) {
        if (name == "history" || name == "software_name" || name == "software_version")
            continue;
        std::string att_val = std::visit(visitor, att);
        NCPP_ERROR((void)input_parameters_group.putAtt(name, att_val));
    }

    free(attr_storage);
    int16_t syear, smon, sday;
    double secs;
    unix2ymds(isodate2unix(time_coverage_start.c_str()), &syear, &smon, &sday, &secs);
    offset_time = ymds2unix(syear, smon, sday, 0.0);
    sprintf(time_string, "seconds since %04d-%02d-%02d", syear, smon, sday);
    units_since = time_string;
    // set geolocation
    NCPP_ERROR(lat_var = geolocation_data_group.addVar("latitude", netCDF::ncFloat,
                                                       {bins_along_track_dim, bins_across_track_dim}));
    NCPP_ERROR((void)lat_var.putAtt("_FillValue", netCDF::ncFloat, 1, &fillValue));
    NCPP_ERROR((void)lat_var.putAtt("units", "degrees_north"));
    NCPP_ERROR((void)lat_var.putAtt("long_name", "Latitudes of bin locations"));
    constexpr float max_lat = 90.f;
    constexpr float min_lat = -90.f;
    NCPP_ERROR((void)lat_var.putAtt("valid_min", netCDF::ncFloat, 1, &min_lat));
    NCPP_ERROR((void)lat_var.putAtt("valid_max", netCDF::ncFloat, 1, &max_lat));
    NCPP_ERROR(lat_var.setCompression(true, true, 5))
    NCPP_ERROR(lat_var.setChunking(netCDF::NcVar::ChunkMode::nc_CHUNKED, chunking))
    NCPP_ERROR(lat_var.putVar(latitude_data.data()))
    NCPP_ERROR(lon_var = geolocation_data_group.addVar("longitude", netCDF::ncFloat,
                                                       {bins_along_track_dim, bins_across_track_dim}));
    NCPP_ERROR((void)lon_var.putAtt("_FillValue", netCDF::ncFloat, 1, &fillValue))
    NCPP_ERROR((void)lon_var.putAtt("units", "degrees_east"))
    NCPP_ERROR((void)lat_var.putAtt("long_name", "Longitude of bin locations"))
    constexpr float max_lon = 180.f;
    constexpr float min_lon = -180.f;
    NCPP_ERROR((void)lon_var.putAtt("valid_min", netCDF::ncFloat, 1, &min_lon));
    NCPP_ERROR((void)lon_var.putAtt("valid_max", netCDF::ncFloat, 1, &max_lon));
    NCPP_ERROR(lon_var.setCompression(true, true, 5))
    NCPP_ERROR(lon_var.setChunking(netCDF::NcVar::ChunkMode::nc_CHUNKED, chunking))
    NCPP_ERROR(lon_var.putVar(longitude_data.data()))

    NCPP_ERROR(height_var = geolocation_data_group.addVar("height", netCDF::ncShort,
                                                          {bins_along_track_dim, bins_across_track_dim}));
    NCPP_ERROR((void)height_var.putAtt("_FillValue", netCDF::ncShort, 1, &fillValueInt))
    NCPP_ERROR((void)height_var.putAtt("units", "meters"))
    NCPP_ERROR((void)height_var.putAtt("long_name", "Altitude at bin locations"))
    NCPP_ERROR((void)height_var.putAtt("valid_min", netCDF::ncShort, 1, &min_height));
    NCPP_ERROR((void)height_var.putAtt("valid_max", netCDF::ncShort, 1, &max_height));
    NCPP_ERROR(height_var.setCompression(true, true, 5))
    NCPP_ERROR(height_var.setChunking(netCDF::NcVar::ChunkMode::nc_CHUNKED, chunking))
    NCPP_ERROR(height_var.putVar(height_data.data()))

    NCPP_ERROR(height_stdev_var = geolocation_data_group.addVar(
                   "height_stdev", netCDF::ncShort, {bins_along_track_dim, bins_across_track_dim}));
    NCPP_ERROR((void)height_stdev_var.putAtt("_FillValue", netCDF::ncShort, 1, &fillValueInt))
    NCPP_ERROR((void)height_stdev_var.putAtt("units", "meters"))
    NCPP_ERROR(
        (void)height_stdev_var.putAtt("long_name", "Standard deviation of terrain altitude within bin"))
    NCPP_ERROR((void)height_stdev_var.putAtt("valid_min", netCDF::ncShort, 1, &min_height_stdev));
    NCPP_ERROR((void)height_stdev_var.putAtt("valid_max", netCDF::ncShort, 1, &max_height));
    NCPP_ERROR(height_stdev_var.setCompression(true, true, 5))
    NCPP_ERROR(height_stdev_var.setChunking(netCDF::NcVar::ChunkMode::nc_CHUNKED, chunking))
    // create variables in observation data
    NCPP_ERROR(i_var = observation_data_group.addVar("i", netCDF::ncFloat,
                                                     {bins_along_track_dim, bins_across_track_dim,
                                                      number_of_views_dim, intensity_bands_per_view_dim}))
    NCPP_ERROR(
        i_stdev_var = observation_data_group.addVar(
            "i_stdev", netCDF::ncFloat,
            {bins_along_track_dim, bins_across_track_dim, number_of_views_dim, intensity_bands_per_view_dim}))
    NCPP_ERROR(number_of_observations_var = observation_data_group.addVar(
                   "number_of_observations", netCDF::ncShort,
                   {bins_along_track_dim, bins_across_track_dim, number_of_views_dim}))
    NCPP_ERROR(sensor_azimuth_angles_var = geolocation_data_group.addVar(
                   "sensor_azimuth_angle", netCDF::ncShort,
                   {bins_along_track_dim, bins_across_track_dim, number_of_views_dim}))
    NCPP_ERROR(sensor_zenith_angles_var = geolocation_data_group.addVar(
                   "sensor_zenith_angle", netCDF::ncShort,
                   {bins_along_track_dim, bins_across_track_dim, number_of_views_dim}))
    NCPP_ERROR(solar_azimuth_angles_var = geolocation_data_group.addVar(
                   "solar_azimuth_angle", netCDF::ncShort,
                   {bins_along_track_dim, bins_across_track_dim, number_of_views_dim}))
    NCPP_ERROR(solar_zenith_angles_var = geolocation_data_group.addVar(
                   "solar_zenith_angle", netCDF::ncShort,
                   {bins_along_track_dim, bins_across_track_dim, number_of_views_dim}))
    NCPP_ERROR(scattering_angles_var = geolocation_data_group.addVar(
                   "scattering_angle", netCDF::ncShort,
                   {bins_along_track_dim, bins_across_track_dim, number_of_views_dim}))
    int deflateLevel = 5;
    NCPP_ERROR((void)i_var.putAtt("_FillValue", netCDF::ncFloat, 1, &fillValue))
    NCPP_ERROR((void)i_stdev_var.putAtt("_FillValue", netCDF::ncFloat, 1, &fillValue))
    NCPP_ERROR((void)number_of_observations_var.putAtt("_FillValue", netCDF::ncShort, 1, &fillValueInt))
    NCPP_ERROR((void)sensor_azimuth_angles_var.putAtt("_FillValue", netCDF::ncShort, 1, &fillValueInt))
    NCPP_ERROR((void)sensor_zenith_angles_var.putAtt("_FillValue", netCDF::ncShort, 1, &fillValueInt))
    NCPP_ERROR((void)solar_azimuth_angles_var.putAtt("_FillValue", netCDF::ncShort, 1, &fillValueInt))
    NCPP_ERROR((void)solar_zenith_angles_var.putAtt("_FillValue", netCDF::ncShort, 1, &fillValueInt))
    NCPP_ERROR((void)scattering_angles_var.putAtt("_FillValue", netCDF::ncShort, 1, &fillValueInt))
    NCPP_ERROR((void)sensor_azimuth_angles_var.putAtt("scale_factor", netCDF::ncFloat, 1, &scale_factor))
    NCPP_ERROR((void)sensor_azimuth_angles_var.putAtt("add_offset", netCDF::ncFloat, 1, &add_offset))
    NCPP_ERROR((void)solar_azimuth_angles_var.putAtt("scale_factor", netCDF::ncFloat, 1, &scale_factor))
    NCPP_ERROR((void)solar_azimuth_angles_var.putAtt("add_offset", netCDF::ncFloat, 1, &add_offset))
    NCPP_ERROR((void)sensor_zenith_angles_var.putAtt("scale_factor", netCDF::ncFloat, 1, &scale_factor))
    NCPP_ERROR((void)sensor_zenith_angles_var.putAtt("add_offset", netCDF::ncFloat, 1, &add_offset))
    NCPP_ERROR((void)solar_zenith_angles_var.putAtt("scale_factor", netCDF::ncFloat, 1, &scale_factor))
    NCPP_ERROR((void)solar_zenith_angles_var.putAtt("add_offset", netCDF::ncFloat, 1, &add_offset))
    NCPP_ERROR((void)scattering_angles_var.putAtt("scale_factor", netCDF::ncFloat, 1, &scale_factor))
    NCPP_ERROR((void)scattering_angles_var.putAtt("add_offset", netCDF::ncFloat, 1, &add_offset))

    NCPP_ERROR(i_var.setCompression(true, true, deflateLevel))
    NCPP_ERROR(i_var.setChunking(netCDF::NcVar::ChunkMode::nc_CHUNKED, chunking))
    NCPP_ERROR(i_stdev_var.setCompression(true, true, deflateLevel))
    NCPP_ERROR(i_stdev_var.setChunking(netCDF::NcVar::ChunkMode::nc_CHUNKED, chunking))
    NCPP_ERROR(sensor_azimuth_angles_var.setCompression(true, true, deflateLevel))
    NCPP_ERROR(sensor_zenith_angles_var.setCompression(true, true, deflateLevel))
    NCPP_ERROR(solar_azimuth_angles_var.setCompression(true, true, deflateLevel))
    NCPP_ERROR(solar_zenith_angles_var.setCompression(true, true, deflateLevel))
    NCPP_ERROR(scattering_angles_var.setCompression(true, true, deflateLevel))
    NCPP_ERROR(number_of_observations_var.setCompression(true, true, deflateLevel))
    NCPP_ERROR(sensor_azimuth_angles_var.setChunking(netCDF::NcVar::ChunkMode::nc_CHUNKED, chunking))
    NCPP_ERROR(sensor_zenith_angles_var.setChunking(netCDF::NcVar::ChunkMode::nc_CHUNKED, chunking))
    NCPP_ERROR(solar_azimuth_angles_var.setChunking(netCDF::NcVar::ChunkMode::nc_CHUNKED, chunking))
    NCPP_ERROR(solar_zenith_angles_var.setChunking(netCDF::NcVar::ChunkMode::nc_CHUNKED, chunking))
    NCPP_ERROR(scattering_angles_var.setChunking(netCDF::NcVar::ChunkMode::nc_CHUNKED, chunking));
    NCPP_ERROR(number_of_observations_var.setChunking(netCDF::NcVar::ChunkMode::nc_CHUNKED, chunking));
    // set valid min, valid max
    NCPP_ERROR((void)i_var.putAtt("valid_min", netCDF::ncFloat, 1, &min_i))
    NCPP_ERROR((void)i_var.putAtt("valid_max", netCDF::ncFloat, 1, &max_i))
    NCPP_ERROR((void)i_stdev_var.putAtt("valid_min", netCDF::ncFloat, 1, &min_i_stdev))
    NCPP_ERROR((void)i_stdev_var.putAtt("valid_max", netCDF::ncFloat, 1, &max_i_stdev))
    NCPP_ERROR((void)sensor_azimuth_angles_var.putAtt("valid_min", netCDF::ncShort, 1, &min_azimuth))
    NCPP_ERROR((void)sensor_azimuth_angles_var.putAtt("valid_max", netCDF::ncShort, 1, &max_azimuth))
    NCPP_ERROR((void)sensor_zenith_angles_var.putAtt("valid_min", netCDF::ncShort, 1, &min_zenith))
    NCPP_ERROR((void)sensor_zenith_angles_var.putAtt("valid_max", netCDF::ncShort, 1, &max_zenith))
    NCPP_ERROR((void)solar_azimuth_angles_var.putAtt("valid_min", netCDF::ncShort, 1, &min_azimuth))
    NCPP_ERROR((void)solar_azimuth_angles_var.putAtt("valid_max", netCDF::ncShort, 1, &max_azimuth))
    NCPP_ERROR((void)solar_zenith_angles_var.putAtt("valid_min", netCDF::ncShort, 1, &min_zenith))
    NCPP_ERROR((void)solar_zenith_angles_var.putAtt("valid_max", netCDF::ncShort, 1, &max_zenith))
    NCPP_ERROR((void)scattering_angles_var.putAtt("valid_min", netCDF::ncShort, 1, &min_scattering))
    NCPP_ERROR((void)scattering_angles_var.putAtt("valid_max", netCDF::ncShort, 1, &max_scattering))
    NCPP_ERROR(
        (void)number_of_observations_var.putAtt("valid_min", netCDF::ncShort, 1, &min_number_of_observations))
    NCPP_ERROR(
        (void)number_of_observations_var.putAtt("valid_max", netCDF::ncShort, 1, &max_number_of_observations))
    // set CF convection attributes
    NCPP_ERROR((void)i_var.putAtt("long_name", "I Stokes vector component"))
    NCPP_ERROR((void)i_var.putAtt("coordinates", "geolocation_data/longitude geolocation_data/latitude"))
    NCPP_ERROR((void)i_var.putAtt("units", "W m^-2 sr^-1 um^-1"))
    NCPP_ERROR((void)i_stdev_var.putAtt("long_name", "Random stdev of i in bin"))
    NCPP_ERROR(
        (void)i_stdev_var.putAtt("coordinates", "geolocation_data/longitude geolocation_data/latitude"))
    NCPP_ERROR((void)i_stdev_var.putAtt("units", "W m^-2 sr^-1 um^-1"))
    // "Observations contributing to bin from each view"
    NCPP_ERROR((void)number_of_observations_var.putAtt(
        "long_name", "Number of observations contributing to bin from each view"))
    NCPP_ERROR((void)number_of_observations_var.putAtt(
        "coordinates", "geolocation_data/longitude geolocation_data/latitude"));
    // angles
    NCPP_ERROR((void)sensor_azimuth_angles_var.putAtt("long_name", "Sensor azimuth angle at bin locations"));
    NCPP_ERROR((void)sensor_zenith_angles_var.putAtt("long_name", "Sensor zenith angle at bin locations"));
    NCPP_ERROR((void)solar_azimuth_angles_var.putAtt("long_name", "Solar azimuth angle at bin locations"));
    NCPP_ERROR((void)solar_zenith_angles_var.putAtt("long_name", "Solar zenith angle at bin locations"));
    NCPP_ERROR((void)scattering_angles_var.putAtt("long_name", "Scattering angle at bin locations"));
    // put units, degrees all
    NCPP_ERROR((void)sensor_azimuth_angles_var.putAtt("units", "degrees"));
    NCPP_ERROR((void)sensor_zenith_angles_var.putAtt("units", "degrees"));
    NCPP_ERROR((void)solar_azimuth_angles_var.putAtt("units", "degrees"));
    NCPP_ERROR((void)solar_zenith_angles_var.putAtt("units", "degrees"));
    NCPP_ERROR((void)scattering_angles_var.putAtt("units", "degrees"));

    // bin_attributes
    NCPP_ERROR(nadir_view_time_var =
                   bin_attributes_group.addVar("nadir_view_time", netCDF::ncDouble, {bins_along_track_dim}));
    NCPP_ERROR(view_time_offsets_var = bin_attributes_group.addVar(
                   "view_time_offsets", netCDF::ncDouble,
                   {bins_along_track_dim, bins_across_track_dim, number_of_views_dim}));
    NCPP_ERROR((void)nadir_view_time_var.putAtt("_FillValue", netCDF::ncDouble, 1, &fillValue));
    NCPP_ERROR((void)nadir_view_time_var.putAtt("long_name", "Time bin was viewed at nadir view"));
    NCPP_ERROR((void)nadir_view_time_var.putAtt("units", units_since));
    NCPP_ERROR((void)nadir_view_time_var.putAtt("valid_min", netCDF::ncDouble, 1, &nadir_time_min));
    NCPP_ERROR((void)nadir_view_time_var.putAtt("valid_max", netCDF::ncDouble, 1, &nadir_time_max));
    NCPP_ERROR((void)view_time_offsets_var.putAtt("_FillValue", netCDF::ncDouble, 1, &fillValue));
    NCPP_ERROR((void)view_time_offsets_var.putAtt("long_name", "Time offsets of views from nadir view"));
    NCPP_ERROR((void)view_time_offsets_var.putAtt("units", "seconds"));
    NCPP_ERROR((void)view_time_offsets_var.putAtt("valid_min", netCDF::ncDouble, 1, &offset_time_min));
    NCPP_ERROR((void)view_time_offsets_var.putAtt("valid_max", netCDF::ncDouble, 1, &offset_time_max));
    NCPP_ERROR(view_time_offsets_var.setCompression(true, true, deflateLevel));
    NCPP_ERROR(view_time_offsets_var.setChunking(netCDF::NcVar::ChunkMode::nc_CHUNKED, chunking));
    netCDF::NcVar intensity_wavelength_var, intensity_f0_var, intensity_bandpass_var;
    NCPP_ERROR(intensity_wavelength_var =
                   sensor_views_bands_group.addVar("intensity_wavelength", netCDF::ncFloat,
                                                   {number_of_views_dim, intensity_bands_per_view_dim}));

    NCPP_ERROR((void)intensity_wavelength_var.putAtt("long_name",
                                                     "Intensity field center wavelengths at each view"));
    NCPP_ERROR((void)intensity_wavelength_var.putAtt("units", "nm"));
    NCPP_ERROR((void)intensity_wavelength_var.putAtt("_FillValue", netCDF::ncFloat, fillValue));
    NCPP_ERROR((void)intensity_wavelength_var.putAtt("valid_min", netCDF::ncFloat, min_wavelength));
    NCPP_ERROR((void)intensity_wavelength_var.putAtt("valid_max", netCDF::ncFloat, max_wavelength));
    NCPP_ERROR(intensity_wavelength_var.putVar(l1bfile.get_l1file().fwave));

    NCPP_ERROR(intensity_f0_var = sensor_views_bands_group.addVar(
                   "intensity_f0", netCDF::ncFloat, {number_of_views_dim, intensity_bands_per_view_dim}));

    NCPP_ERROR((void)intensity_f0_var.putAtt("long_name", "Intensity band solar irradiance"));
    NCPP_ERROR((void)intensity_f0_var.putAtt("units", "W m^-2 um^-1"));
    NCPP_ERROR((void)intensity_f0_var.putAtt("_FillValue", netCDF::ncFloat, fillValue));
    NCPP_ERROR((void)intensity_f0_var.putAtt("valid_min", netCDF::ncFloat, min_f0));
    NCPP_ERROR((void)intensity_f0_var.putAtt("valid_max", netCDF::ncFloat, max_f0));
    std::vector<float> f0(l1bfile.get_l1file().Fobar, l1bfile.get_l1file().Fobar + n_bands);
    // multiply all elements of f0 by 10
    std::transform(f0.begin(), f0.end(), f0.begin(), [](float x) { return x * 10; });
    NCPP_ERROR(intensity_f0_var.putVar(f0.data()));

    NCPP_ERROR(intensity_bandpass_var =
                   sensor_views_bands_group.addVar("intensity_bandpass", netCDF::ncFloat,
                                                   {number_of_views_dim, intensity_bands_per_view_dim}));

    NCPP_ERROR((void)intensity_bandpass_var.putAtt("long_name", "Intensity field bandpass at each view"));
    NCPP_ERROR((void)intensity_bandpass_var.putAtt("units", "nm"));
    NCPP_ERROR((void)intensity_bandpass_var.putAtt("_FillValue", netCDF::ncFloat, fillValue));
    NCPP_ERROR((void)intensity_bandpass_var.putAtt("valid_min", netCDF::ncFloat, min_width));
    NCPP_ERROR((void)intensity_bandpass_var.putAtt("valid_max", netCDF::ncFloat, max_width));
    if (l1bfile.get_format().type == FT_OCIL1B) {
        std::vector<float> fwhm(n_bands, oci_bandpass);
        NCPP_ERROR(intensity_bandpass_var.putVar(fwhm.data()))
    };

    NCPP_ERROR(sensor_view_angle_var = sensor_views_bands_group.addVar("sensor_view_angle", netCDF::ncFloat,
                                                                       {bins_along_track_dim}));
    NCPP_ERROR((void)sensor_view_angle_var.putAtt("long_name", "Tilt angle along track"));
    NCPP_ERROR((void)sensor_view_angle_var.putAtt("units", "degrees"));
    NCPP_ERROR((void)sensor_view_angle_var.putAtt("_FillValue", netCDF::ncFloat, fillValue));
    NCPP_ERROR((void)sensor_view_angle_var.putAtt("valid_min", netCDF::ncFloat, min_tilt));
    NCPP_ERROR((void)sensor_view_angle_var.putAtt("valid_max", netCDF::ncFloat, max_tilt));
    tilt = std::vector<float>(n_lines, fillValue);
    number_of_observations_line = std::vector<size_t>(n_lines, 0);
    scan_quality_var =
        bin_attributes_group.addVar("scan_quality_flags", netCDF::ncUbyte, {bins_along_track_dim});
    NCPP_ERROR((void)scan_quality_var.putAtt("long_name", "Scan quality flags"));
    NCPP_ERROR((void)scan_quality_var.putAtt("_FillValue", netCDF::ncUbyte, fill_value_char));
    NCPP_ERROR((void)scan_quality_var.putAtt("valid_min", netCDF::ncUbyte, min_uchar));
    NCPP_ERROR((void)scan_quality_var.putAtt("valid_max", netCDF::ncUbyte, max_uchar));
    scan_quality_flags = std::vector<u_char>(n_lines, fill_value_char);
    if (l1bfile.get_format().type == FT_OCIL1B) {
        NCPP_ERROR(qc = observation_data_group.addVar("qc", netCDF::ncUbyte,
                                                      {bins_along_track_dim, bins_across_track_dim,
                                                       number_of_views_dim, intensity_bands_per_view_dim}));
        quality_flags = std::vector<u_char>(n_bands * n_pixels * n_lines);
        NCPP_ERROR(qc.setCompression(true, true, deflateLevel));
        NCPP_ERROR(qc.setChunking(netCDF::NcVar::ChunkMode::nc_CHUNKED, chunking));
        NCPP_ERROR((void)qc.putAtt("_FillValue", netCDF::ncUbyte, 1, &fill_value_char));
        NCPP_ERROR((void)qc.putAtt("long_name", "quality indicator"));
        NCPP_ERROR((void)qc.putAtt("valid_min", netCDF::ncUbyte, 1, &min_uchar));
        NCPP_ERROR((void)qc.putAtt("valid_max", netCDF::ncUbyte, 1, &max_uchar));
    }
    cloud_correction = !((std::get<std::string>(input_attributes.at("cloud_anc_files"))).empty());
    if (cloud_correction) {
        NCPP_ERROR(aggregate_height_var = geolocation_data_group.addVar(
                       "aggregated_height", netCDF::ncShort, {bins_along_track_dim, bins_across_track_dim}));
        NCPP_ERROR((void)aggregate_height_var.putAtt("_FillValue", netCDF::ncShort, 1, &fillValueInt))
        NCPP_ERROR((void)aggregate_height_var.putAtt("units", "meters"))
        NCPP_ERROR((void)aggregate_height_var.putAtt("long_name", "Aggregated height at bin locations"))
        constexpr short max_aggregated_height = 30000;
        NCPP_ERROR((void)aggregate_height_var.putAtt("valid_min", netCDF::ncShort, 1, &min_height));
        NCPP_ERROR(
            (void)aggregate_height_var.putAtt("valid_max", netCDF::ncShort, 1, &max_aggregated_height));
        NCPP_ERROR(aggregate_height_var.setCompression(true, true, 5))
        NCPP_ERROR(aggregate_height_var.setChunking(netCDF::NcVar::ChunkMode::nc_CHUNKED, chunking))

        NCPP_ERROR(aggregate_height_stdev_var =
                       geolocation_data_group.addVar("aggregated_height_stdev", netCDF::ncShort,
                                                     {bins_along_track_dim, bins_across_track_dim}));
        NCPP_ERROR((void)aggregate_height_stdev_var.putAtt("_FillValue", netCDF::ncShort, 1, &fillValueInt))
        NCPP_ERROR((void)aggregate_height_stdev_var.putAtt("units", "meters"))
        NCPP_ERROR((void)aggregate_height_stdev_var.putAtt(
            "long_name", "Aggregated height standard deviation at bin locations"))
        constexpr short max_aggregated_height_stdev = 30000;
        NCPP_ERROR(
            (void)aggregate_height_stdev_var.putAtt("valid_min", netCDF::ncShort, 1, &min_height_stdev));
        NCPP_ERROR((void)aggregate_height_stdev_var.putAtt("valid_max", netCDF::ncShort, 1,
                                                           &max_aggregated_height_stdev));
        NCPP_ERROR(aggregate_height_stdev_var.setCompression(true, true, 5))
        NCPP_ERROR(aggregate_height_stdev_var.setChunking(netCDF::NcVar::ChunkMode::nc_CHUNKED, chunking))
        aggregated_height.resize(n_lines * n_pixels, BAD_FLT);
        aggregated_height_stdev.resize(n_lines * n_pixels, BAD_FLT);
    }
}

/**
 * @brief Finalizes and closes the L1C netCDF output file after binning is complete.
 *
 */
void L1CFile::close() {
    if (file.isNull())
        return;
    std::cout << "Writing to file: " << o_file << std::endl;
    NCPP_ERROR(i_var.putVar(i_data.data()))
    // take square root
    for (size_t i = 0; i < n_pixels * n_lines * n_bands; i++) {
        if (i_stdev_data[i] > 0 && number_of_observations_band[i] > 1)
            i_stdev_data[i] = sqrt(i_stdev_data[i]);
        else
            i_stdev_data[i] = BAD_FLT;
    }
    NCPP_ERROR(i_stdev_var.putVar(i_stdev_data.data()))
    // put angles, but rescale them
    std::vector<short> buffer(n_pixels * n_lines);
    for (size_t i = 0; i < n_pixels * n_lines; i++) {
        if (sensor_azimuth_angles_data[i] != BAD_FLT)
            buffer[i] = static_cast<short>((sensor_azimuth_angles_data[i] - add_offset) / scale_factor);
        else
            buffer[i] = fillValueInt;
    }
    NCPP_ERROR(sensor_azimuth_angles_var.putVar(buffer.data()))
    for (size_t i = 0; i < n_pixels * n_lines; i++) {
        if (sensor_zenith_angles_data[i] != BAD_FLT)
            buffer[i] = static_cast<short>((sensor_zenith_angles_data[i] - add_offset) / scale_factor);
        else
            buffer[i] = fillValueInt;
    }
    NCPP_ERROR(sensor_zenith_angles_var.putVar(buffer.data()))
    for (size_t i = 0; i < n_pixels * n_lines; i++) {
        if (solar_azimuth_angles_data[i] != BAD_FLT)
            buffer[i] = static_cast<short>((solar_azimuth_angles_data[i] - add_offset) / scale_factor);
        else
            buffer[i] = fillValueInt;
    }
    NCPP_ERROR(solar_azimuth_angles_var.putVar(buffer.data()))
    for (size_t i = 0; i < n_pixels * n_lines; i++) {
        if (solar_zenith_angles_data[i] != BAD_FLT)
            buffer[i] = static_cast<short>((solar_zenith_angles_data[i] - add_offset) / scale_factor);
        else
            buffer[i] = fillValueInt;
    }
    NCPP_ERROR(solar_zenith_angles_var.putVar(buffer.data()))
    for (size_t i = 0; i < n_pixels * n_lines; i++) {
        if (scattering_angles_data[i] != BAD_FLT)
            buffer[i] = static_cast<short>((scattering_angles_data[i] - add_offset) / scale_factor);
        else
            buffer[i] = fillValueInt;
    }
    NCPP_ERROR(scattering_angles_var.putVar(buffer.data()))
    number_of_observations_var.putVar(number_of_observations.data());
    for (size_t i_line = 0; i_line < n_lines; i_line++) {
        for (size_t i_pixel = 0; i_pixel < n_pixels; i_pixel++) {
            size_t index = i_line * n_pixels + i_pixel;
            if (number_of_observations[index] > 0)
                view_time_offsets_data[index] = nadir_view_time_data[i_line] - view_time_offsets_data[index];
        }
    }
    NCPP_ERROR(nadir_view_time_var.putVar(nadir_view_time_data.data()))
    NCPP_ERROR(view_time_offsets_var.putVar(view_time_offsets_data.data()))
    if (!qc.isNull()) {
        NCPP_ERROR(qc.putVar(quality_flags.data()));
    }
    NCPP_ERROR(sensor_view_angle_var.putVar(tilt.data()));
    NCPP_ERROR(scan_quality_var.putVar(scan_quality_flags.data()));
    // put height stdev
    for (size_t i = 0; i < n_pixels * n_lines; i++) {
        height_stdev[i] -= height_data[i] * height_data[i];
        if (height_stdev[i] > 0 && number_of_observations_band[i] > 1)
            buffer[i] = static_cast<short>(std::sqrt(height_stdev[i]));
        else
            buffer[i] = fillValueInt;
    }
    NCPP_ERROR(height_stdev_var.putVar(buffer.data()));
    if (cloud_correction) {
        for (size_t i = 0; i < n_pixels * n_lines; i++) {
            buffer[i] = static_cast<short>(aggregated_height[i]);
        }
        NCPP_ERROR(aggregate_height_var.putVar(buffer.data()));
        for (size_t i = 0; i < n_pixels * n_lines; i++) {
            if (aggregated_height_stdev[i] > 0 && number_of_observations_band[i] > 1)
                buffer[i] = static_cast<short>(std::sqrt(aggregated_height_stdev[i]));
            else
                buffer[i] = fillValueInt;
        }
        NCPP_ERROR(aggregate_height_stdev_var.putVar(buffer.data()));
    }
    file.close();
    std::cout << "Done." << std::endl;
}