#include "L1BFile.h"

/**
 * @brief Constructs an L1B object file and opens an L1 sensor data file for reading.
 *
 * @param[in] ifile Path to the L1 input file to open.
 * @param[in] bin_indexes Rvalue reference to vector of bin indices mapping L1 pixels to L1C bins.
 *                         Size must equal n_pixels * n_lines. Use -1 for invalid pixels.
 *
 */
L1BFile::L1BFile(const std::string& ifile, std::vector<int>&& bin_indexes) : ifile(ifile) {
    StdoutRedirector redirector;
    redirector.redirect();
    this->bin_indexes = std::move(bin_indexes);
    filehandle_init(&l1file);
    clo_optionList_t* optionList = clo_createList();
    l1_add_options(optionList);

    // ignore extra options in msl12_defaults.par files
    clo_setEnableExtraOptions(1);
    l1_read_default_files(optionList, &l1file, ifile.c_str());
    clo_setEnableExtraOptions(0);
    l1_load_options(optionList, &l1file);
    if (std::getenv("OCDATAROOT") == nullptr) {
        fprintf(stderr, "-E-: %s:%d  Env Var OCDATAROOT is not defined  \n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    strcpy(l1file.name, ifile.c_str());
    format = getFormat(l1file.name);
    redirector.restore();
    if (format.type == FT_INVALID) {
        fprintf(stderr, "-E- %s:%d: Could not find type for file %s.\n", __FILE__, __LINE__, ifile.c_str());
        exit(EXIT_FAILURE);
    }
    l1file.format = format.type;
    l1file.sensorID = format.sensor_id;
    l1file.subsensorID = format.subsensor_id;
    redirector.redirect();
    l1file.nbands = rdsensorinfo(l1file.sensorID, l1_input->evalmask, "Nbands", NULL);
    int status = openl1(&l1file);
    if (status != 0) {
        fprintf(stderr, "-E- %s:%d: Unable to open L1 record.\n", __FILE__, __LINE__);
        exit(1);
    }
    if (alloc_l1(&l1file, &l1rec) == 0) {
        fprintf(stderr, "-E- %s:%d: Unable to allocate L1 record.\n", __FILE__, __LINE__);
        exit(1);
    }
    redirector.restore();
    // set pixel boundaries:
    l1file.spix = 0;
    l1file.epix = l1file.npix;
    ncid = l1file.sd_id;
    n_bands = l1file.nbands;
    n_pixels = l1file.npix;
    n_lines = l1file.nscan;

    if (format.type == FT_OCIL1B) {
        NC_ERROR(nc_inq_ncid(ncid, "observation_data", &observation_grp_id));
        NC_ERROR(nc_inq_varid(observation_grp_id, "qual_red", &qual_red_id));
        NC_ERROR(nc_inq_varid(observation_grp_id, "qual_blue", &qual_blue_id));
        NC_ERROR(nc_inq_varid(observation_grp_id, "qual_SWIR", &qual_swir_id));
        int dimids[NC_MAX_VAR_DIMS];
        int ndims;
        NC_ERROR(nc_inq_var(observation_grp_id, qual_red_id, NULL, NULL, &ndims, dimids, NULL));
        if (ndims != 3) {
            fprintf(stderr, "-E-: Expected 3 dimensions, got %d\n. %s:%d", ndims, __FILE__, __LINE__);
            exit(EXIT_FAILURE);
        }
        NC_ERROR(nc_inq_dimlen(observation_grp_id, dimids[0], &n_red_bands));
        NC_ERROR(nc_inq_dimlen(observation_grp_id, dimids[1], &n_lines));
        NC_ERROR(nc_inq_dimlen(observation_grp_id, dimids[2], &n_pixels));
        assert((int)n_pixels == l1file.npix);
        assert((int)n_lines == l1file.nscan);

        NC_ERROR(nc_inq_var(observation_grp_id, qual_blue_id, NULL, NULL, &ndims, dimids, NULL));
        if (ndims != 3) {
            fprintf(stderr, "-E-: Expected 3 dimensions, got %d\n. %s:%d", ndims, __FILE__, __LINE__);
            exit(EXIT_FAILURE);
        }
        NC_ERROR(nc_inq_dimlen(observation_grp_id, dimids[0], &n_blue_bands));
        NC_ERROR(nc_inq_dimlen(observation_grp_id, dimids[1], &n_lines));
        NC_ERROR(nc_inq_dimlen(observation_grp_id, dimids[2], &n_pixels));
        assert((int)n_pixels == l1file.npix);
        assert((int)n_lines == l1file.nscan);

        NC_ERROR(nc_inq_var(observation_grp_id, qual_swir_id, NULL, NULL, &ndims, dimids, NULL));
        if (ndims != 3) {
            fprintf(stderr, "-E-: Expected 3 dimensions, got %d\n. %s:%d", ndims, __FILE__, __LINE__);
            exit(EXIT_FAILURE);
        }
        NC_ERROR(nc_inq_dimlen(observation_grp_id, dimids[0], &n_swir_bands));
        NC_ERROR(nc_inq_dimlen(observation_grp_id, dimids[1], &n_lines));
        NC_ERROR(nc_inq_dimlen(observation_grp_id, dimids[2], &n_pixels));
        assert((int)n_pixels == l1file.npix);
        assert((int)n_lines == l1file.nscan);
        quality_flags = std::vector<u_char>(n_pixels * n_bands);
        quality_flags_blue = std::vector<u_char>(n_pixels * n_blue_bands);
        quality_flags_red = std::vector<u_char>(n_pixels * n_red_bands);
        qual_SWIR = std::vector<u_char>(n_pixels * n_swir_bands);
        max_blue_index = expected_num_blue_bands - skipped_blue_bands;
        max_red_index = max_blue_index + expected_num_red_bands;
    }
    // find qual red, blue and SWIR:
}

int L1BFile::read_line(int line) {
    if (format.type == FT_OCIL1B) {
        size_t start[3] = {0, (size_t)line, 0}; /* band, line, pixel */
        std::vector<size_t> count = {n_blue_bands, 1, n_pixels};
        nc_get_vara_ubyte(observation_grp_id, qual_blue_id, start, count.data(), quality_flags_blue.data());
        int band = 0;
        for (int sb = 0; sb < max_blue_index; sb++) {
            for (int pixel = 0; pixel < (int)n_pixels; pixel++) {
                int isp = n_pixels * sb + pixel;
                int ibp = n_bands * pixel + band;
                quality_flags[ibp] = quality_flags_blue[isp];
            }
            band++;
        }
        count = {n_red_bands, 1, n_pixels};
        nc_get_vara_ubyte(observation_grp_id, qual_red_id, start, count.data(), quality_flags_red.data());
        for (int sb = 0; sb < (int)n_red_bands; sb++) {
            for (int pixel = 0; pixel < (int)n_pixels; pixel++) {
                int isp = n_pixels * sb + pixel;
                int ibp = n_bands * pixel + band;
                quality_flags[ibp] = quality_flags_red[isp];
            }
            band++;
        }
        count = {n_swir_bands, 1, n_pixels};
        nc_get_vara_ubyte(observation_grp_id, qual_red_id, start, count.data(), qual_SWIR.data());
        for (int sb = 0; sb < (int)n_swir_bands; sb++) {
            for (int pixel = 0; pixel < (int)n_pixels; pixel++) {
                int isp = n_pixels * sb + pixel;
                int ibp = n_bands * pixel + band;
                quality_flags[ibp] = quality_flags_red[isp];
            }
            if (sb == 3 || sb == 6)
                continue;
            band++;
        }
        if (band != (int)n_bands) {
            fprintf(stderr, "-E- : mismatch %d, %d\n", band, (int)n_bands);
            exit(EXIT_FAILURE);
        }
    }

    return readl1(&l1file, line, &l1rec);
}

/**
 * @brief Constructs an L1B File with area weighting and cloud-corrected height data.
 *
 * @param[in] ifile Path to the L1 input file to open.
 * @param[in] bin_indexes vector of bin indices mapping L1 pixels to L1C bins.
 *                         Size must equal n_pixels * n_lines. Use -1 for invalid pixels.
 * @param[in] area_weights 2D vector where each element contains pairs of
 *                          (L1C bin index, area fraction) for sub-pixel mapping.
 *                          Size must match bin_indexes if non-empty.
 *                          Element i corresponds to L1 pixel i.
 * @param[in] height_data_corrected vector of cloud-corrected height values.
 *                                   Size must match bin_indexes if non-empty.
 *                                   Values are in kilometers
 */
L1BFile::L1BFile(const std::string& ifile, std::vector<int>&& bin_indexes,
                 std::vector<std::vector<std::pair<int, double>>>&& area_weights,
                 std::vector<float>&& height_data_corrected, std::vector<float>&& senz_corrected,
                 std::vector<float>&& sena_corrected)
    : L1BFile(ifile, std::move(bin_indexes)) {
    this->area_weights = std::move(area_weights);
    this->height_data_corrected = std::move(height_data_corrected);
    this->senz_corrected = std::move(senz_corrected);
    this->sena_corrected = std::move(sena_corrected);
    if (!this->area_weights.empty()) {
        if (this->area_weights.size() != this->bin_indexes.size()) {
            fprintf(stderr, "-E-: %s:%d area_weights and bin_indexes mismatch, %ld vs %ld  \n", __FILE__,
                    __LINE__, this->area_weights.size(), this->bin_indexes.size());
            exit(EXIT_FAILURE);
        }
    }
    if (!this->height_data_corrected.empty()) {
        if (this->height_data_corrected.size() != this->bin_indexes.size()) {
            fprintf(stderr, "-E-: %s:%d height_data_corrected and bin_indexes mismatch, %ld vs %ld  \n",
                    __FILE__, __LINE__, this->height_data_corrected.size(), this->bin_indexes.size());
            exit(EXIT_FAILURE);
        }
    }
    if (!this->senz_corrected.empty()) {
        if (this->senz_corrected.size() != this->bin_indexes.size()) {
            fprintf(stderr, "-E-: %s:%d senz_corrected and bin_indexes mismatch, %ld vs %ld  \n", __FILE__,
                    __LINE__, this->senz_corrected.size(), this->bin_indexes.size());
            exit(EXIT_FAILURE);
        }
    }
    if (!this->sena_corrected.empty()) {
        if (this->sena_corrected.size() != this->bin_indexes.size()) {
            fprintf(stderr, "-E-: %s:%d sena_corrected and bin_indexes mismatch, %ld vs %ld  \n", __FILE__,
                    __LINE__, this->sena_corrected.size(), this->bin_indexes.size());
            exit(EXIT_FAILURE);
        }
    }
}