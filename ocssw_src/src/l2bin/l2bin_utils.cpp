#include "l2bin_utils.h"
#include "constants.h"
#include <hdfi.h>

#ifdef MAC_MEM
#include <mach/mach.h>
#endif

namespace {
float lonRanges[4] = {0, -180, 180, 0};
float northmost = -90.0, southmost = 90.0, eastmost = -180.0, westmost = 180.0;

}  // namespace

void append_composite_product(std::string &products_requested, const std::string &composite_product_name,
                              const std::string &composite_scheme, int &min_max_scheme) {
    if (products_requested.find(composite_product_name) == std::string::npos) {
        products_requested.append(",");
        products_requested.append(composite_product_name);
    }
    if (composite_scheme == "max") {
        min_max_scheme = 1;
    } else if (composite_scheme == "min") {
        min_max_scheme = 0;
    } else {
        std::cerr << "-E-: unidentified composite scheme: " << composite_scheme << ". See " << __FILE__ << ":"
                  << __LINE__ << std::endl;
        exit(EXIT_FAILURE);
    }
}

void printMemoryInfo(const std::string &message) {
#ifndef MAC_MEM
    std::ifstream file("/proc/self/status");
    std::string line;
    std::cout << message;
    while (std::getline(file, line)) {
        if (line.find("Vm") != std::string::npos) {
            std::cout << line << std::endl;
        }
    }
#else
    task_vm_info_data_t vm_info;
    mach_msg_type_number_t count = TASK_VM_INFO_COUNT;
    kern_return_t result = task_info(mach_task_self(), TASK_VM_INFO, (task_info_t)&vm_info, &count);
    std::cout << message << std::endl;
    if (result == KERN_SUCCESS) {
        std::cout << "Virtual memory size " << vm_info.virtual_size / 1024.0 << " KB" << std::endl;
        std::cout << "Physical footprint memory size " << vm_info.phys_footprint / 1024.0 << " KB"
                  << std::endl;
        std::cout << "RRS memory size " << vm_info.resident_size / 1024.0 << " KB" << std::endl;
        std::cout << "PAGE memory size " << vm_info.page_size / 1024.0 << " KB" << std::endl;
    } else {
        std::cout << "--Error--: couldn't get VM info" << std::endl;
    }

#endif
}
bool binIntersectsPixel(l3::L3ShapeIsine *shape, int32_t row, int32_t col, Box_t &pixelBox,
                        double &areaFrac) {
    // Initialize result and get bin boundaries
    bool result = false;
    double n, s, e, w;
    shape->rowcol2bounds(row, col, n, s, e, w);

    // Create box representing the bin
    Box_t box(Point_t(w, s), Point_t(e, n));
    areaFrac = 0;

    // Handle pixel box crossing 180/-180 meridian
    std::array<Box_t, 2> rotated_boxes{};
    size_t number_of_boxes = rotate_polygon(pixelBox, rotated_boxes);

    // Check intersection with each rotated box
    for (size_t ib = 0; ib < number_of_boxes; ib++) {
        if (Box_t chopped_box = rotated_boxes[ib]; !bg::disjoint(box, chopped_box)) {
            // Calculate intersection if boxes overlap
            if (Box_t output; bg::intersection(box, chopped_box, output)) {
                // If intersection has positive area, calculate area fraction
                if (double intersectArea = bg::area(output); intersectArea > 0) {
                    result = true;
                    double binArea = (n - s) * (e - w);
                    areaFrac = intersectArea / binArea;
                    break;
                }
            }
        }
    }
    return result;
}

bool binIntersectsPixel(l3::L3ShapeIsine *shape, int32_t row, int32_t col, Polygon_t &pixelPoly,
                        double &areaFrac) {
    // Initialize result and get bin boundaries
    bool result = false;
    double n, s, e, w;
    shape->rowcol2bounds(row, col, n, s, e, w);

    // Create box representing the bin
    Box_t box(Point_t(w, s), Point_t(e, n));
    areaFrac = 0;

    // Handle polygon crossing 180/-180 meridian
    std::array<Polygon_t, 2> rotated_polys{};
    size_t number_of_polys = rotate_polygon(pixelPoly, rotated_polys);

    // Check intersection with each rotated polygon
    for (size_t ip = 0; ip < number_of_polys; ip++) {
        if (!bg::disjoint(box, rotated_polys[ip])) {
            // Calculate intersection if shapes overlap
            std::deque<Polygon_t> output;
            if (bg::intersection(box, rotated_polys[ip], output)) {
                double binArea = (n - s) * (e - w);
                // Process each intersection polygon
                BOOST_FOREACH (Polygon_t const &p, output) {
                    double intersectArea = bg::area(p);
                    if (intersectArea > 0.0) {
                        result = true;
                        areaFrac += intersectArea / binArea;
                    }
                }
            }
        }
    }
    return result;
}

void getBins(l3::L3ShapeIsine *shape, AreaWeighting &areaWeight, size_t ipixl,
             std::map<uint64_t, double> &areas) {
    areas.clear();

    double lat0 = areaWeight.get_latitude(ipixl);
    double lon0 = areaWeight.get_longitude(ipixl);
    int32_t row0 = shape->lat2row(lat0);
    lat0 = shape->row2lat(row0);
    double lat = lat0;
    double deltaLat = 180.0 / shape->getNumRows();
    if (areaWeight.get_mode() == 1) {
        Box_t pixelBox;
        // use the pixel box
        pixelBox.min_corner().set<0>(areaWeight.get_longitude(ipixl) - areaWeight.get_delta_lon(ipixl));
        pixelBox.min_corner().set<1>(areaWeight.get_latitude(ipixl) - areaWeight.get_delta_lat(ipixl));
        pixelBox.max_corner().set<0>(areaWeight.get_longitude(ipixl) + areaWeight.get_delta_lon(ipixl));
        pixelBox.max_corner().set<1>(areaWeight.get_latitude(ipixl) + areaWeight.get_delta_lat(ipixl));
        get_areas(lat, lat0, lon0, deltaLat, shape, pixelBox, areas);
    }

    else if (areaWeight.get_mode() == 2) {
        Box_t pixelBox;

        // use the pixel bounding box
        float latMax =
            std::max({areaWeight.get_upperCornerLat(ipixl), areaWeight.get_upperCornerLat(ipixl + 1),
                      areaWeight.get_lowerCornerLat(ipixl), areaWeight.get_lowerCornerLat(ipixl + 1)});

        float latMin =
            std::min({areaWeight.get_upperCornerLat(ipixl), areaWeight.get_upperCornerLat(ipixl + 1),
                      areaWeight.get_lowerCornerLat(ipixl), areaWeight.get_lowerCornerLat(ipixl + 1)});
        float lon1 = areaWeight.get_upperCornerLon(ipixl);
        float lon2 = areaWeight.get_upperCornerLon(ipixl + 1);
        float lon3 = areaWeight.get_lowerCornerLon(ipixl);
        float lon4 = areaWeight.get_lowerCornerLon(ipixl + 1);

        // Check if longitude values need to be rotated
        // This occurs when pixels span the international dateline (180°/-180°)
        // and have longitude values with different signs (+/-).
        // Note: A single pixel's longitude span should never exceed 180°
        if (std::abs(lon1 - lon2) > 180.0 || std::abs(lon1 - lon3) > 180.0 || std::abs(lon1 - lon4) > 180.0) {
            // Some longitude values need to be rotated to handle the international dateline crossing
            // We'll rotate negative longitudes to their positive equivalents (e.g. -170° -> 190°)
            // The direction of rotation doesn't matter since we'll handle box rotation later if needed
            if (lon1 < 0)
                lon1 += 360.0;
            if (lon2 < 0)
                lon2 += 360.0;
            if (lon3 < 0)
                lon3 += 360.0;
            if (lon4 < 0)
                lon4 += 360.0;
        }

        float lonMax = std::max({lon1, lon2, lon3, lon4});
        float lonMin = std::min({lon1, lon2, lon3, lon4});

        pixelBox.min_corner().set<0>(lonMin);
        pixelBox.min_corner().set<1>(latMin);
        pixelBox.max_corner().set<0>(lonMax);
        pixelBox.max_corner().set<1>(latMax);
        get_areas(lat, lat0, lon0, deltaLat, shape, pixelBox, areas);
    } else {
        // use the exact pixel polygon
        Polygon_t pixelPoly;
        bg::append(pixelPoly.outer(),
                   Point_t(areaWeight.get_upperCornerLon(ipixl), areaWeight.get_upperCornerLat(ipixl)));
        bg::append(pixelPoly.outer(),
                   Point_t(areaWeight.get_lowerCornerLon(ipixl), areaWeight.get_lowerCornerLat(ipixl)));
        bg::append(pixelPoly.outer(), Point_t(areaWeight.get_lowerCornerLon(ipixl + 1),
                                              areaWeight.get_lowerCornerLat(ipixl + 1)));
        bg::append(pixelPoly.outer(), Point_t(areaWeight.get_upperCornerLon(ipixl + 1),
                                              areaWeight.get_upperCornerLat(ipixl + 1)));
        bg::append(pixelPoly.outer(),
                   Point_t(areaWeight.get_upperCornerLon(ipixl), areaWeight.get_upperCornerLat(ipixl)));

        // make sure the polygon is defined properly
        bg::correct(pixelPoly);

        get_areas(lat, lat0, lon0, deltaLat, shape, pixelPoly, areas);
    }
}

void addPixelToBin(std::vector<float *> &l2data, int32_t *qual_flags, L2BinStruct &l2binStruct, size_t ifile,
                   size_t krow, size_t n_rows_in_group, size_t n_bins_in_group, size_t n_allocperbin,
                   size_t l3b_nprod, size_t ipixl, uint64_t bin, double areaFrac) {
    auto &basebin = *l2binStruct.basebin;
    auto lastfile = l2binStruct.lastfile;
    auto nscenes = l2binStruct.nscenes;
    auto data_areas = l2binStruct.data_areas;
    auto data_values = l2binStruct.data_values;
    auto file_index = l2binStruct.file_index;
    auto data_quality = l2binStruct.data_quality;
    auto nobs = l2binStruct.nobs;
    auto allocated_space = l2binStruct.allocated_space;
    int32_t ibin = bin - basebin[krow];
    /* if bin not within bin row group then continue */
    /* --------------------------------------------- */
    if (ibin < 0 || ibin >= static_cast<int32_t>(n_bins_in_group) || bin >= basebin[krow + n_rows_in_group] ||
        bin < basebin[krow])
        return;

    if (static_cast<int16_t>(ifile) != lastfile[ibin]) {
        nscenes[ibin]++;
        lastfile[ibin] = static_cast<int16_t>(ifile);
    }
    l2binStruct.allocate(n_allocperbin, l3b_nprod, ibin);
    file_index[ibin][nobs[ibin]] = static_cast<int16_t>(ifile);
    if (qual_flags) {
        data_quality[ibin][nobs[ibin]] = qual_flags[ipixl];
        // need composite products as well
    }
    /* Get data area for pixel */
    /* ----------------------------------- */
    data_areas[ibin][nobs[ibin]] = areaFrac;
    for (size_t l3_iprod = 0; l3_iprod < l3b_nprod; l3_iprod++) {
        float f32 = l2data[l3_iprod][ipixl];
        data_values[ibin][l3b_nprod * nobs[ibin] + l3_iprod] = f32;
    }
    /* iprod loop */

    /* Increment number of observations in bin */
    /* --------------------------------------- */
    nobs[ibin]++;
    /* Reallocate if necessary */
    /* ----------------------- */
    if (nobs[ibin] == allocated_space[ibin]) {
        l2binStruct.reallocate(n_allocperbin, l3b_nprod, ibin);
    }
}

bool updateBinCompositeScheme(float *composite_data, L2BinStruct &l2binStruct, size_t prod_index,
                              int composite_scheme, size_t krow, size_t l3b_nprod, size_t ipixl,
                              uint64_t bin) {
    auto data_values = l2binStruct.data_values;
    auto nobs = l2binStruct.nobs;
    auto &basebin = *l2binStruct.basebin;
    int32_t ibin = bin - basebin[krow];
    if (ibin < 0 || ibin >= static_cast<int32_t>(l2binStruct.number_of_bins))
        return false;
    l2binStruct.allocate(min_per_bin_allocation, l3b_nprod, ibin);
    if (nobs[ibin] != 0) {
        float f32 = composite_data[ipixl];
        if (composite_scheme) {
            if (f32 < data_values[ibin][prod_index])
                return false;
        } else {
            if (f32 > data_values[ibin][prod_index])
                return false;
        }
        nobs[ibin] = 0;
    }
    return true;
}

int64_t getbinnum(l3::L3ShapeIsine *shape, float *latitude, float *longitude, size_t ipixl) {
    return shape->latlon2bin(latitude[ipixl], longitude[ipixl]);
}

void fill_data(netCDF::NcGroup &binned_data, L2BinStruct &l2binStruct, size_t n_filled_bins,
               size_t n_bins_in_group, size_t n_products) {
    auto nobs = l2binStruct.nobs;
    auto data_values = l2binStruct.data_values;
    auto data_areas = l2binStruct.data_areas;
    auto file_index = l2binStruct.file_index;
    /* Allocate sum & sum-squared arrays */
    /* --------------------------------- */
    std::vector<double> sum_bin(n_bins_in_group), sum2_bin(n_bins_in_group);
    for (size_t iprod = 0; iprod < n_products; iprod++) {
        std::memset(sum_bin.data(), 0, sum_bin.size() * sizeof(double));
        std::memset(sum2_bin.data(), 0, sum2_bin.size() * sizeof(double));
        /* Process bins */
        /* ------------ */
        for (size_t ibin = 0; ibin < n_bins_in_group; ibin++) {
            if (nobs[ibin] == 0)
                continue;

            /* Process data file by file */
            /* ------------------------- */
            int32_t npix_file;  // for weighting spatial average
            double sum_file;    // accumulators for each file
            double sum2_file;
            double area_file;
            int16_t prev_file;  // index of previous file considered

            // initialize sum for this bin with first file
            int16_t j = 0;
            float pixval = data_values[ibin][j * n_products + iprod];
            double pixarea = data_areas[ibin][j];
            npix_file = 1;
            (void)npix_file;
            area_file = pixarea;
            sum_file = pixval * pixarea;
            sum2_file = pixval * pixarea * pixval * pixarea;
            prev_file = file_index[ibin][j];

            // add weighted data for rest of observations (files)
            for (j = 1; j < nobs[ibin]; j++) {
                pixval = data_values[ibin][j * n_products + iprod];
                pixarea = data_areas[ibin][j];

                if (file_index[ibin][j] == prev_file) {  // same file
                    npix_file += 1;
                    area_file += pixarea;
                    sum_file += pixval * pixarea;
                    sum2_file += pixval * pixarea * pixval * pixarea;
                } else {  // new file

                    // finalize weighted sum for previous file
                    sum_bin[ibin] += (sum_file / sqrt(area_file));
                    sum2_bin[ibin] += (sum2_file / sqrt(area_file));

                    // initialize sum for current file
                    npix_file = 1;
                    area_file = pixarea;
                    sum_file = pixval * pixarea;
                    sum2_file = pixval * pixarea * pixval * pixarea;
                    prev_file = file_index[ibin][j];
                }
            } /* observation loop */

            // add data from final file
            sum_bin[ibin] += (sum_file / std::sqrt(area_file));
            sum2_bin[ibin] += (sum2_file / std::sqrt(area_file));

        } /* ibin loop */

        /* Write Product Vdatas */
        /* -------------------- */
        std::vector<float> productData(2 * n_filled_bins);
        /* Fill bin data array */
        /* ------------------- */
        size_t i = 0;
        for (size_t ibin = 0; ibin < n_bins_in_group; ibin++) {
            if (nobs[ibin] != 0) {
                productData[i * 2] = sum_bin[ibin];
                productData[i * 2 + 1] = sum2_bin[ibin];
                i++;
            }
        }

        writeBinData_nc(binned_data.getId(), i, iprod, productData.data());
    }
}

std::vector<uint8_t> fill_bins(bool read_quality, L2BinStruct &l2binStruct, l3::L3ShapeIsine *shape,
                               std::vector<binListStruct64_nc> &binList64nc,
                               std::vector<binListStruct_nc> &binList32nc, size_t &n_filled_bins,
                               size_t n_bins_in_group, size_t n_products, size_t krow, size_t nfiles,
                               size_t qual_max_allowed, bool is64bit, float *bounds) {
    /* Fill "Bin List" vdata array */
    /* --------------------------- */
    auto nscenes = l2binStruct.nscenes;
    auto data_areas = l2binStruct.data_areas;
    auto data_values = l2binStruct.data_values;
    auto file_index = l2binStruct.file_index;
    auto data_quality = l2binStruct.data_quality;
    auto nobs = l2binStruct.nobs;
    auto &basebin = *l2binStruct.basebin;
    std::vector<uint8_t> best_qual = std::vector<uint8_t>(n_bins_in_group, 255);
    if (is64bit)
        binList64nc = std::vector<binListStruct64_nc>(n_filled_bins);
    else
        binList32nc = std::vector<binListStruct_nc>(n_filled_bins);

    size_t filled_bin_index = 0;
    for (size_t ibin = 0; ibin < n_bins_in_group; ibin++) {
        /* Adjust for bins with "bad" quality values */
        /* ----------------------------------------- */
        if (read_quality && nobs[ibin] > 0) {
            best_qual[ibin] = 255;
            for (int16_t j = 0; j < nobs[ibin]; j++)
                if (data_quality[ibin][j] < best_qual[ibin])
                    best_qual[ibin] = data_quality[ibin][j];

            int16_t valid_observation_count = 0;
            for (int16_t j = 0; j < nobs[ibin]; j++) {
                if ((data_quality[ibin][j] <= best_qual[ibin]) &&
                    (data_quality[ibin][j] <= qual_max_allowed)) {
                    if (valid_observation_count < j) {
                        data_areas[ibin][valid_observation_count] = data_areas[ibin][j];
                        for (size_t iprod = 0; iprod < n_products; iprod++) {
                            data_values[ibin][valid_observation_count * n_products + iprod] =
                                data_values[ibin][j * n_products + iprod];
                        }
                    }
                    valid_observation_count++;
                }
            }
            nobs[ibin] = valid_observation_count;

            if (nobs[ibin] == 0)
                n_filled_bins--;
        }

        if (nobs[ibin] != 0) {
            uint64_t bin = (uint64_t)ibin + basebin[krow];

            if (is64bit) {
                binList64nc[filled_bin_index].binnum = bin;
                binList64nc[filled_bin_index].nobs = nobs[ibin];
                binList64nc[filled_bin_index].nscenes = nscenes[ibin];
            } else {
                binList32nc[filled_bin_index].binnum = bin;
                binList32nc[filled_bin_index].nobs = nobs[ibin];
                binList32nc[filled_bin_index].nscenes = nscenes[ibin];
            }

            /* weights {=sqrt(# of L2 files in given bin)} */
            /* ------------------------------------------- */
            double wgt = 0.0;
            for (size_t ifile = 0; ifile <= nfiles; ifile++) {
                double area = 0.0;
                for (int16_t j = 0; j < nobs[ibin]; j++) {
                    if (file_index[ibin][j] == static_cast<int16_t>(ifile))
                        area += data_areas[ibin][j];
                }
                wgt += std::sqrt(area);
            }

            if (is64bit) {
                binList64nc[filled_bin_index].weights = wgt;
                binList64nc[filled_bin_index].time_rec = 0;
            } else {
                binList32nc[filled_bin_index].weights = wgt;
                binList32nc[filled_bin_index].time_rec = 0;
            }

            filled_bin_index++;

            /* Update Max/Min Lon/Lat */
            /* ---------------------- */
            double latbin, lonbin;
            shape->bin2latlon(bin, latbin, lonbin);

            if (latbin > northmost)
                northmost = latbin;
            if (latbin < southmost)
                southmost = latbin;

            float minNeg = lonRanges[0];  // min, max in Western hemisphere
            float maxNeg = lonRanges[1];
            float minPos = lonRanges[2];  // min, max in Eastern hemisphere
            float maxPos = lonRanges[3];
            if (lonbin < 0) {  // Western hemisphere
                minNeg = fmin(minNeg, lonbin);
                maxNeg = fmax(maxNeg, lonbin);
            } else {  // Eastern hemisphere
                minPos = fmin(minPos, lonbin);
                maxPos = fmax(maxPos, lonbin);
            }
            // Adjust east and west granule bounding coordinates
            float max_angle = 90.0;
            int8_t lon000 = minPos - maxNeg < max_angle;
            int8_t lon180 = maxPos - minNeg > 360.0 - max_angle;
            int8_t polar = (lon000 && lon180) || (90.0 - fmax(northmost, -1 * southmost) <= FLT_EPSILON);

            if (minNeg >= maxNeg) {  // Entirely in Eastern hemisphere
                eastmost = maxPos;
                westmost = minPos;
            } else if (minPos >= maxPos) {  // Entirely in Western hemisphere
                eastmost = maxNeg;
                westmost = minNeg;
            } else if (polar) {  // Polar granule
                eastmost = 180;
                westmost = -180;
                if (northmost > 0)
                    northmost = 90;  // North pole
                else
                    southmost = -90;  // South pole
            } else if (lon000) {      // Granule crosses Longitude 0
                eastmost = maxPos;
                westmost = minNeg;
            } else if (lon180) {  // Granule crosses Longitude 180
                eastmost = maxNeg;
                westmost = minPos;
            }
            lonRanges[0] = minNeg;
            lonRanges[1] = maxNeg;
            lonRanges[2] = minPos;
            lonRanges[3] = maxPos;
        } /* nobs[ibin] != 0 */
    } /* ibin loop */
    bounds[0] = northmost;
    bounds[1] = southmost;
    bounds[2] = eastmost;
    bounds[3] = westmost;
    return best_qual;
}

void fill_quality_data(netCDF::NcGroup &binned_data, L2BinStruct &l2binStruct, size_t n_filled_bins,
                       size_t n_bins_in_group, std::vector<uint8_t> &best_qual) {
    auto nobs = l2binStruct.nobs;
    std::vector<uint8_t> qualityData(n_filled_bins, 1);

    size_t i = 0;
    for (size_t ibin = 0; ibin < n_bins_in_group; ibin++) {
        if (nobs[ibin] != 0) {
            memcpy(&qualityData[i], &best_qual[ibin], 1);  // is memcpy needed?
            i++;
        }
    }
    writeQuality_nc(binned_data.getId(), n_filled_bins, (VOIDP)qualityData.data());
}

void update_bin_index(L2BinStruct &l2binStruct, l3::L3ShapeIsine *shape, bool is64bit, size_t n_filled_bins,
                      size_t n_bins_in_group, size_t nrows, size_t krow, size_t n_rows_in_group,
                      int64_t &total_filled_bins, std::vector<uint64_t> &beg, std::vector<uint32_t> &ext,
                      std::vector<binIndexStruct64_nc> &binIndex64nc,
                      std::vector<binIndexStruct_nc> &binIndex32nc) {
    /* Compute "begin" & "extent" vdata entries */
    /* ---------------------------------------- */
    auto nobs = l2binStruct.nobs;
    auto &basebin = *l2binStruct.basebin;

    std::vector<int64_t> binnum_data(n_filled_bins);

    size_t i = 0;
    for (size_t ibin = 0; ibin < n_bins_in_group; ibin++) {
        if (nobs[ibin] != 0) {
            binnum_data[i] = (int64_t)ibin + basebin[krow];

            if (i < 0 || i >= n_filled_bins) {
                printf("Error: %zu %zu %zu %zu\n", i, ibin, n_filled_bins, n_bins_in_group);
            }
            i++;
        }
    }

    get_beg_ext(n_filled_bins, binnum_data.data(), (int64_t *)basebin.data(), nrows, (int64_t *)beg.data(),
                (int32_t *)ext.data());

    total_filled_bins += n_filled_bins;
    // free(best_qual);

    /* Fill BinIndex Vdata */
    /* ------------------- */
    for (i = 0; i < n_rows_in_group; i++) {
        size_t i32 = i + krow;
        if (i32 < 0 || i32 >= nrows) {
            printf("Error: %zu %zu\n", i, krow);
            exit(-1);
        }
        set_bin_index(i32, is64bit, shape, basebin, beg, ext, binIndex64nc, binIndex32nc);
    }
}

void write_l3_metadata(std::vector<uint64_t> &basebin, netCDF::NcFile &ofile, meta_l3bType &meta_l3b,
                       instr &input, size_t nfiles, size_t nrows, size_t total_filled_bins,
                       std::vector<L2_Reader> &l2_files, std::vector<std::string> &product_list, int argc,
                       char **argv, std::vector<int> &brk_scan, bool is64bit,
                       float *geospatial_bounds) { /*
                                                    * determine list of files actually used in the bin output
                                                    */
    std::vector<int> fileused(nfiles, 0);
    if (input.fileuse[0] != 0) {
        auto fp2 = fopen(input.fileuse, "w");
        for (size_t ifile = 0; ifile < nfiles; ifile++) {
            if (brk_scan[ifile] != brake_scan_fill_value) {
                fileused[ifile] = 1;
                fprintf(fp2, "%s\n", l2_files[ifile].get_filename().c_str());
            }
        }
        fclose(fp2);
    } else {
        for (size_t ifile = 0; ifile < nfiles; ifile++) {
            if (brk_scan[ifile] != brake_scan_fill_value) {
                fileused[ifile] = 1;
            }
        }
    }
    // writing metadata

    /* Read and write global attributes */
    /* -------------------------------- */
    if (input.verbose)
        printf("Writing Global Attributes\n");

    strncpy(meta_l3b.product_name, input.ofile, sizeof(meta_l3b.product_name) - 1);
    std::set<std::string> platforms;
    std::set<std::string> sensors;
    // reading missions/sensors
    for (size_t ifile = 0; ifile < nfiles; ifile++) {
        platforms.insert(l2_files[ifile].get_platform());
        sensors.insert(l2_files[ifile].get_instrument());
    }
    // setting platform
    if (platforms.empty()) {
        std::cerr << "No platforms found" << std::endl;
    } else {
        if (platforms.size() > 1) {
            std::cerr << "Warning : You are binning data from different platforms. " << __FILE__ << ":"
                      << __LINE__ << std::endl;
        }
        std::string platform = *platforms.begin();  // take first sensor
        strncpy(meta_l3b.mission, platform.c_str(), sizeof(meta_l3b.mission) - 1);
    }

    // setting instrument
    if (sensors.empty()) {
        std::cerr << "Warning : No sensors found. " << __FILE__ << ":" << __LINE__ << std::endl;
    } else {
        if (sensors.size() > 1) {
            std::cerr << "Warning : You are binning data from different sensors. " << __FILE__ << ":"
                      << __LINE__ << std::endl;
        }
        std::string sensor_name = *sensors.begin();  // take first sensor
        strncpy(meta_l3b.sensor_name, sensor_name.c_str(), sizeof(meta_l3b.sensor_name) - 1);
        int sensorID = sensorName2SensorId(sensor_name.c_str());
        if (sensorID < 0) {
            strncpy(meta_l3b.sensor, sensor_name.c_str(), sizeof(meta_l3b.sensor) - 1);
        } else {
            meta_l3b.sensorID = sensorID;
            std::string sensor = sensorId2InstrumentName(sensorID);
            strncpy(meta_l3b.sensor, sensor.c_str(), sizeof(meta_l3b.sensor) - 1);
            strncpy(meta_l3b.mission, sensorId2PlatformName(sensorID), sizeof(meta_l3b.mission) - 1);
        }

        strncpy(meta_l3b.title, meta_l3b.sensor_name, sizeof(meta_l3b.title) - 1);
        strncat(meta_l3b.title, " Level-3 Binned Data", sizeof(meta_l3b.title) - strlen(meta_l3b.title) - 1);
    }

    if (input.suite[0]) {
        strncat(meta_l3b.title, " ", sizeof(meta_l3b.title) - strlen(meta_l3b.title) - 1);
        strncat(meta_l3b.title, input.suite,
                sizeof(meta_l3b.title) - strlen(meta_l3b.title) - strlen(input.suite));
    }

    strncpy(meta_l3b.prod_type, input.prodtype, sizeof(meta_l3b.prod_type) - 1);

    strncpy(meta_l3b.pversion, input.pversion, sizeof(meta_l3b.pversion) - 1);
    strncpy(meta_l3b.soft_name, "l2bin", sizeof(meta_l3b.soft_name) - 1);
    strncpy(meta_l3b.soft_ver, VERSION, sizeof(meta_l3b.soft_ver) - 1);
    // ignore orbit numbers
    meta_l3b.end_orb = 0;
    meta_l3b.start_orb = 0;
    double startTime = yds2unix(2100, 1, 0);
    double endTime = yds2unix(1900, 1, 0);
    double tmpTime;
    strncpy(meta_l3b.infiles, l2_files[0].get_filename().c_str(), sizeof(meta_l3b.infiles) - 1);
    for (size_t i = 0; i < nfiles; i++) {
        if (fileused[i] == 1) {
            // update start/end times
            tmpTime = l2_files[i].get_start_time();
            if (tmpTime < startTime)
                startTime = tmpTime;
            tmpTime = l2_files[i].get_end_time();
            if (tmpTime > endTime)
                endTime = tmpTime;
        }
        // append to the file list - all files input, not just those used in the binning.
        if (i > 0) {
            strncat(meta_l3b.infiles, ",", sizeof(meta_l3b.infiles) - strlen(meta_l3b.infiles) - 1);
            strncat(meta_l3b.infiles, l2_files[i].get_filename().c_str(),
                    sizeof(meta_l3b.infiles) - strlen(meta_l3b.infiles) - 1);
        }
    }

    meta_l3b.startTime = startTime;
    meta_l3b.endTime = endTime;

    strncpy(meta_l3b.binning_scheme, "Integerized Sinusoidal Grid", sizeof(meta_l3b.binning_scheme) - 1);

    meta_l3b.north = geospatial_bounds[0];
    meta_l3b.south = geospatial_bounds[1];
    meta_l3b.east = geospatial_bounds[2];
    meta_l3b.west = geospatial_bounds[3];
    char buf[65535];
    strncpy(buf, argv[0], sizeof(buf));
    for (int i = 1; i < argc; i++) {
        strncat(buf, " ", sizeof(buf) - strlen(buf) - 1);
        strncat(buf, argv[i], sizeof(buf) - strlen(buf) - 1);
    }
    strncpy(meta_l3b.proc_con, buf, sizeof(meta_l3b.proc_con) - 1);
    strncpy(meta_l3b.input_parms, input.parms, sizeof(meta_l3b.input_parms) - 1);
    strncpy(meta_l3b.flag_names, input.flaguse, sizeof(meta_l3b.flag_names) - 1);

    buf[0] = 0;
    std::vector<std::string> units = l2_files[0].get_units();
    for (size_t iprod = 0; iprod < product_list.size(); iprod++) {
        char *prodName = strdup(product_list[iprod].c_str());
        strncat(buf, prodName, sizeof(buf) - strlen(buf) - 1);
        free(prodName);
        strncat(buf, ":", sizeof(buf) - strlen(buf) - 1);
        strncat(buf, units[iprod].c_str(), sizeof(buf) - strlen(buf) - 1);
        strncat(buf, ",", sizeof(buf) - strlen(buf) - 1);
        if (strlen(buf) >= MD_ATTRSZ) {
            printf("-E- units metadata is too long.  Remove some products from product list.\n");
            exit(EXIT_FAILURE);
        }
    }
    buf[strlen(buf) - 1] = 0;
    strcpy(meta_l3b.units, buf);

    meta_l3b.data_bins = total_filled_bins;
    int64_t totalNumBins;
    totalNumBins = basebin[nrows] - 1;
    meta_l3b.pct_databins = (float)(total_filled_bins) / totalNumBins * 100.0;

    strncpy(meta_l3b.doi, input.doi, sizeof(meta_l3b.doi) - 1);

    const char *keywordStr = getGCMDKeywords(input.suite);
    if (keywordStr)
        strncpy(meta_l3b.keywords, keywordStr, sizeof(meta_l3b.keywords) - 1);
    else
        strncpy(meta_l3b.keywords, "", sizeof(meta_l3b.keywords) - 1);
    time_t tnow;
    time(&tnow);
    strncpy(meta_l3b.ptime, unix2isodate(tnow, 'G'), sizeof(meta_l3b.ptime) - 1);
    idDS ds_id;
    ds_id.fid = ofile.getId();
    ds_id.sid = -1;
    ds_id.fftype = DS_NCDF;  // FMT_L2NCDF
    write_l3b_meta_netcdf4(ds_id, &meta_l3b, is64bit);
}

bool exists_test(const std::string &name) {
    std::ifstream f(name.c_str());
    return f.good();
}

bool skip_DL(float lon, int side, int night_flag, time_t end_day, time_t beg_day) {
    if (night_flag == 1) {
        if ((side == -1) && (beg_day == -1) && (lon < 0))
            return true;
        if ((side == +1) && (end_day == 0) && (lon > 0))
            return true;
    } else {
        if ((side == -1) && (beg_day <= 0) && (lon < 0))
            return true;
        if ((side == +1) && (end_day >= 0) && (lon > 0))
            return true;
    }
    return false;
}

std::map<std::string, std::string> get_oprodname(const std::string &oprodname) {
    std::map<std::string, std::string> output_l3_filenames;
    if (!oprodname.empty()) {
        std::vector<std::string> output_pairs_l2_l3_names;
        boost::split(output_pairs_l2_l3_names, std::string(oprodname), boost::is_any_of(","));
        for (const auto &pair_l2_l3 : output_pairs_l2_l3_names) {
            std::vector<std::string> l2_l3_name;
            boost::split(l2_l3_name, pair_l2_l3, boost::is_any_of(":"));
            boost::trim(l2_l3_name[0]);
            boost::trim(l2_l3_name[1]);
            output_l3_filenames[l2_l3_name[0]] = l2_l3_name[1];
        }
    }
    return output_l3_filenames;
}

void set_breakscan(int &brk_scan, L2_Reader &l2file, float deltaeqcross, int night, int startdate,
                   int enddate) {
    auto [dataday0, dataday1] = get_datadays(l2file.get_geospatial(), deltaeqcross, night);
    if (dataday1 == dataday0) {
        if (dataday0 < startdate || dataday1 > enddate)
            brk_scan = brake_scan_fill_value;
        else
            brk_scan = 0;
    } else {
        if (dataday1 < startdate)  // startdate is dataday conversion of input sday
            brk_scan = brake_scan_fill_value;
        else if (dataday0 > enddate) {  // enddate is dataday conversion of input eday
            brk_scan = brake_scan_fill_value;
        } else {
            if (dataday1 > enddate)
                brk_scan = 1;
            else
                brk_scan = -1;
        }
    }
}

std::vector<std::string> set_nc_prodnames(int deflate, std::vector<std::string> &product_list,
                                          std::map<std::string, std::string> &output_l3_filenames,
                                          netCDF::NcGroup &binned_data) {
    char **prodnames = (char **)malloc(product_list.size() * sizeof(char *));
    std::vector<std::string> output_prod_names;
    for (size_t size = 0; size < product_list.size(); size++) {
        // look for user requested product names
        if (output_l3_filenames.find(product_list[size]) != output_l3_filenames.end())
            prodnames[size] = strdup(output_l3_filenames[product_list[size]].c_str());
        else
            prodnames[size] = strdup(product_list[size].c_str());
        output_prod_names.push_back(prodnames[size]);
    }
    int status = defineBinData_nc(deflate, binned_data.getId(), product_list.size(), prodnames);
    if (status) {
        exit(EXIT_FAILURE);
    }
    // free memory from strdup
    for (size_t size = 0; size < product_list.size(); size++) {
        free(prodnames[size]);
    }
    // free memory from malloc
    free(prodnames);
    return output_prod_names;
}

void set_bin_index(int i, bool is64bit, l3::L3ShapeIsine *shape, std::vector<uint64_t> &basebin,
                   std::vector<uint64_t> &beg, std::vector<uint32_t> &ext,
                   std::vector<binIndexStruct64_nc> &binIndex64nc,
                   std::vector<binIndexStruct_nc> &binIndex32nc) {
    if (is64bit) {
        binIndex64nc[i].start_num = basebin[i];
        binIndex64nc[i].begin = beg[i];
        binIndex64nc[i].extent = ext[i];
        binIndex64nc[i].max = shape->getNumCols(i);
    } else {
        binIndex32nc[i].start_num = basebin[i];
        binIndex32nc[i].begin = beg[i];
        binIndex32nc[i].extent = ext[i];
        binIndex32nc[i].max = shape->getNumCols(i);
    }
}

void ini_bin_index_arrays(instr &input, int nrows, bool is64bit, l3::L3ShapeIsine *shape,
                          std::vector<uint64_t> &basebin, std::vector<uint64_t> &beg,
                          std::vector<uint32_t> &ext, std::vector<binIndexStruct64_nc> &binIndex64nc,
                          std::vector<binIndexStruct_nc> &binIndex32nc, int &ngroup, int &n_rows_in_group) {
    if (is64bit) {
        binIndex64nc.resize(nrows);
    } else {
        binIndex32nc.resize(nrows);
    }
    /* Initialize bin_indx array */
    /* ------------------------- */
    for (int i = 0; i < nrows; i++) {
        set_bin_index(i, is64bit, shape, basebin, beg, ext, binIndex64nc, binIndex32nc);
    }
    // Row Group
    n_rows_in_group = input.rowgroup;
    if (n_rows_in_group <= 0) {
        printf("row_group not defined, using 270.\n");
        n_rows_in_group = 270;
    }
    if (input.verbose)
        printf("%d %d %d\n", input.sday, input.eday, n_rows_in_group);

    /* Find row_group that divides nrows */
    /* --------------------------------- */
    for (int i = nrows; i > 0; i--) {
        if ((nrows % i) == 0) {
            if (i <= n_rows_in_group) {
                n_rows_in_group = i;
                break;
            }
        }
    }
    if (input.rowgroup != n_rows_in_group) {
        printf("Input row_group: %d   Actual row_group: %d\n", input.rowgroup, n_rows_in_group);
    }
    ngroup = nrows / n_rows_in_group;
}

size_t rotate_polygon(Box_t &pixelBox, std::array<Box_t, 2> &boxes) {
    // Get min/max longitude of box
    double min_lon = pixelBox.min_corner().get<0>();
    double max_lon = pixelBox.max_corner().get<0>();

    // Check if box crosses 180° meridian
    if (((max_lon - 180) * (min_lon - 180) < 0)) {
        boxes[0] = pixelBox;
        boxes[1] = pixelBox;
        // Create second box shifted west by 360°
        boxes[1].min_corner().set<0>(min_lon - 360);
        boxes[1].max_corner().set<0>(max_lon - 360);
        return 2;
    }

    // Check if box crosses -180° meridian
    if (((max_lon + 180) * (min_lon + 180) < 0)) {
        boxes[0] = pixelBox;
        boxes[1] = pixelBox;
        // Create second box shifted east by 360°
        boxes[1].min_corner().set<0>(min_lon + 360);
        boxes[1].max_corner().set<0>(max_lon + 360);
        return 2;
    }

    // If box is entirely beyond 180°, shift west by 360°
    if (max_lon > 180 && min_lon > 180) {
        boxes[0] = pixelBox;
        boxes[0].min_corner().set<0>(min_lon - 360);
        boxes[0].max_corner().set<0>(max_lon - 360);
        return 1;
    }

    // If box is entirely beyond -180°, shift east by 360°
    if (max_lon < -180 && min_lon < -180) {
        boxes[0] = pixelBox;
        boxes[0].min_corner().set<0>(min_lon + 360);
        boxes[0].max_corner().set<0>(max_lon + 360);
        return 1;
    }

    // Box is within normal longitude range, no rotation needed
    boxes[0] = pixelBox;
    return 1;
}

size_t rotate_polygon(Polygon_t &pixelPoly, std::array<Polygon_t, 2> &polys) {
    // Get longitude of first vertex as reference
    double lon_0 = pixelPoly.outer()[0].get<0>();
    size_t number_of_vertices = pixelPoly.outer().size();
    bool rotate_lon = false;

    // Check if polygon spans more than 180 degrees by comparing each vertex to first
    // This indicates we need to rotate longitudes
    for (size_t iv = 1; iv < number_of_vertices; iv++) {
        double lon = pixelPoly.outer()[iv].get<0>();
        if (std::abs(lon - lon_0) > 180.0) {
            rotate_lon = true;
            break;
        }
    }

    // If rotation needed, shift negative longitudes to positive by adding 360
    if (rotate_lon) {
        for (size_t iv = 0; iv < number_of_vertices; iv++) {
            double lon = pixelPoly.outer()[iv].get<0>();
            if (lon < 0)
                pixelPoly.outer()[iv].set<0>(lon + 360);
        }
    }

    // Find min and max longitude of polygon vertices
    double min_lon = lon_0;
    double max_lon = lon_0;
    for (size_t iv = 0; iv < number_of_vertices; iv++) {
        double lon = pixelPoly.outer()[iv].get<0>();
        if (min_lon > lon)
            min_lon = lon;
        if (max_lon < lon)
            max_lon = lon;
    }

    // Check if polygon crosses 180 degree meridian (dateline) with lon > 180
    if ((max_lon - 180) * (min_lon - 180) < 0) {
        // Create two polygons - one original and one shifted west by 360
        polys[0] = pixelPoly;
        polys[1] = pixelPoly;
        for (size_t iv = 0; iv < number_of_vertices; iv++) {
            double lon = polys[1].outer()[iv].get<0>();
            polys[1].outer()[iv].set<0>(lon - 360);
        }
        bg::correct(polys[1]);
        return 2;
    }

    // Check if polygon crosses -180 degree meridian with lon < -180
    if ((max_lon + 180) * (min_lon + 180) < 0) {
        // Create two polygons - one original and one shifted east by 360
        polys[0] = pixelPoly;
        polys[1] = pixelPoly;
        for (size_t iv = 0; iv < number_of_vertices; iv++) {
            double lon = polys[1].outer()[iv].get<0>();
            polys[1].outer()[iv].set<0>(lon + 360);
        }
        bg::correct(polys[1]);
        return 2;
    }

    // If polygon entirely beyond 180, shift it west by 360
    if (max_lon > 180 && min_lon > 180) {
        polys[0] = pixelPoly;
        for (size_t iv = 0; iv < number_of_vertices; iv++) {
            double lon = polys[0].outer()[iv].get<0>();
            polys[0].outer()[iv].set<0>(lon - 360);
        }
        bg::correct(polys[0]);
        return 1;
    }

    // If polygon entirely beyond -180, shift it east by 360
    if (max_lon < -180 && min_lon < -180) {
        polys[0] = pixelPoly;
        for (size_t iv = 0; iv < number_of_vertices; iv++) {
            double lon = polys[0].outer()[iv].get<0>();
            polys[0].outer()[iv].set<0>(lon + 360);
        }
        bg::correct(polys[0]);
        return 1;
    }

    // No rotation needed, return original polygon
    polys[0] = pixelPoly;
    return 1;
}