#include "map"
#include "l2_reader.hpp"
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <array>
#include "hdf.h"
#include "setupflags.h"
#include "l2bin_utils.h"
#include "constants.h"
#include <cstdio>

// max number of files
constexpr size_t maxnfiles = 512;

int main(int argc, char **argv) {
    // Input parameters struct
    instr input;
    printf("%s %s (%s %s)\n", PROGRAM, VERSION, __DATE__, __TIME__);
    
    // Parse command line arguments
    l2bin_input(argc, argv, &input, PROGRAM, VERSION);
    
    // Check number of input files doesn't exceed max
    size_t nfiles = input.files.size();
    if (nfiles > maxnfiles) {
        printf("--Error--: number of input files %zu exceeds MAXIMUM number of files %ld\n", nfiles,
               maxnfiles);
        std::cerr << "Exiting...";
        exit(EXIT_FAILURE);
    }
    printf("%zu input files\n", nfiles);

    // Processing flags
    bool read_l2_flags = false;
    bool read_quality = false;
    bool use_composite = false;
    int composite_scheme = -1; // 0 - min, 1 - max, -1 - undefined

    // Input parameters
    std::string products_requested = input.l3bprod;                               // List of L3 products requested for binning
    std::string wavelength_requested = input.output_wavelengths;                  // Wavelengths requested for output products
    std::vector<L2_Reader> l2_files(nfiles);                                     // Vector of L2 file readers, one per input file
    std::vector<AreaWeighting> area_weightings(nfiles);                          // Vector of area weighting objects for each file
    std::vector<size_t> composite_product_index;                                 // Indices for composite product, for each file
    std::vector<std::string> product_list;                                       // List of products found in input files
    std::string flaguse = input.flaguse;                                         // L2 flags to use for filtering
    uint8_t qual_max_allowed = input.qual_max;                                   // Maximum allowed quality level
    std::vector<uint32_t> required(nfiles, 0), flagusemask(nfiles, 0);          // Required flags and flag masks per file
    std::map<std::string, std::string> output_l3_filenames = get_oprodname(input.output_product_names); // Map of output L3 filenames

    // Set processing flags based on input
    if (!flaguse.empty())
        read_l2_flags = true;
    if (!std::string(input.qual_prod).empty())
        read_quality = true;
    if (!std::string(input.composite_prod).empty()) {
        append_composite_product(products_requested, input.composite_prod, input.composite_scheme,
                                 composite_scheme);
        use_composite = true;
        composite_product_index.resize(nfiles);
    }

    // Create a l2 reader object for  each input file
    // product_list is an output and should be the same for all the files (given that all the files contain all needed products)
    for (size_t ifile = 0; ifile < nfiles; ifile++) {
        // Initialize L2 reader and area weighting
        l2_files[ifile] = L2_Reader(input.files[ifile]);
        area_weightings[ifile] = AreaWeighting(input.area_weighting);
        area_weightings[ifile].set_l2_reader(l2_files[ifile]);
        
        // Set up geolocation and variables
        l2_files[ifile].iniGeolocation();
        l2_files[ifile].setVariables(products_requested, wavelength_requested, product_list);

        // Handle L2 flags if needed
        if (read_l2_flags) {
            int status = l2_files[ifile].ini_l2_flags();

            if (status) {
                printf("-E- %s:%d Could not read L2 flags in the input file %s\n", __FILE__, __LINE__,
                       l2_files[ifile].get_filename().c_str());
                exit(EXIT_FAILURE);
            }
            status = setupflags(flaguse, l2_files[ifile].get_l2_meaning_bit_dict(), flagusemask[ifile],
                                required[ifile]);
            if (status) {
                printf("-E- %s:%d Could not find needed flags in the input file %s\n", __FILE__, __LINE__,
                       l2_files[ifile].get_filename().c_str());
                exit(EXIT_FAILURE);
            }
        }

        // Handle quality flags if needed
        if (read_quality) {
            int status = l2_files[ifile].ini_quality_flags(input.qual_prod);
            if (status) {
                printf("-E- %s:%d Could not read quality flags in the input file %s\n", __FILE__, __LINE__,
                       l2_files[ifile].get_filename().c_str());
                exit(EXIT_FAILURE);
            }
        }

        // Handle composite products if needed
        if (use_composite) {
            int status =
                l2_files[ifile].find_product_index(input.composite_prod, composite_product_index[ifile]);
            if (status) {
                printf("-E- %s:%d Could not find composite product in the input file %s\n", __FILE__,
                       __LINE__, l2_files[ifile].get_filename().c_str());
                exit(EXIT_FAILURE);
            }
        }
        l2_files[ifile].reset_cache();
    }

    // Initialize metadata struct
    meta_l3bType meta_l3b;
    int64_t total_filled_bins = 0;
    float geospatial_bounds[4] = {-90.0, 90.0, -180.0, 180.0};

    // Calculate date parameters
    int syear = (int32_t)input.sday / 1000.;
    int sday = input.sday - 1000 * syear;
    int startdate = (int32_t)(yds2unix(syear, sday, 0) / 86400);
    int eyear = (int32_t)input.eday / 1000.;
    int eday = input.eday - 1000 * eyear;
    int enddate = (int32_t)(yds2unix(eyear, eday, 0) / 86400);
    int nrows{-1};
    auto diffday_beg = (time_t)startdate;
    auto diffday_end = (time_t)enddate;

    // Set up grid resolution
    resolve2binRows(input.resolve, nrows, meta_l3b.resolution);
    if (nrows == -1) {
        printf("-E- %s:%d Grid resolution not defined\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    // dlat is used calculate which scans of the file fit in krow of the sinusoidal grid
    float dlat = 180. / nrows;
    // check if 64 bit is needed to avoid overflow
    bool is64bit = nrows > 2160 * 16;
    printf("Resolution: %s\n", input.resolve);
    printf("Max Qual Allowed: %d\n", input.qual_max);

    // Arrays for scan row tracking
    std::array<std::vector<int>, maxnfiles> bscan_row, escan_row;
    std::array<std::vector<u_char>, maxnfiles> scan_in_rowgroup;
    std::vector<int> brk_scan(nfiles);
    // if regional is set and valid data will be binned regardless of date
    bool regional = strcasecmp(input.prodtype, "regional") == 0;

    // Process scan rows for each file
    size_t n_active_files{nfiles};
    for (size_t ifile = 0; ifile < nfiles; ifile++) {
        l2_files[ifile].get_escan_bscan_row(bscan_row[ifile], escan_row[ifile], dlat);
        brk_scan[ifile] = 0;
        if (regional)
            continue;
        set_breakscan(brk_scan[ifile], l2_files[ifile], input.deltaeqcross, input.night, startdate, enddate);
        if (brk_scan[ifile] == brake_scan_fill_value)
            n_active_files--; // determine how many files will contribute to the output. If it's set regional, dates will be ignored.
    }

    // Set up bins
    std::vector<uint64_t> basebin(nrows + 1);
    l3::L3ShapeIsine *shape = new l3::L3ShapeIsine(nrows);
    for (int i = 0; i <= nrows; i++) {
        basebin[i] = shape->getBaseBin(i);
    }
    basebin[nrows] = shape->getNumBins() + 1;

    printf("total number of bins: %ld\n", static_cast<long int>(basebin[nrows]) - 1);
    // Create output file
    std::string ofile = input.ofile;
    if (exists_test(ofile))
        std::remove(ofile.c_str());
    if (ofile.empty())
        EXIT_LOG(std::cerr << "-E-: empty output filename\n")
    netCDF::NcFile l3_file(ofile, netCDF::NcFile::newFile);
    // set the group binned data
    netCDF::NcGroup binned_data = l3_file.addGroup("level-3_binned_data");

    // Define bin list variables
    int status{0};
    if (is64bit)
        status = defineBinList64_nc(input.deflate, binned_data.getId());
    else
        status = defineBinList_nc(input.deflate, binned_data.getId());
    if (status) {
        printf("-E- %s:%d Could not define binList variable in output file\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }

    // Set up output product names and attributes
    std::vector<std::string> output_prod_names =
        set_nc_prodnames(input.deflate, product_list, output_l3_filenames, binned_data);
    size_t n_products = product_list.size();
    for (size_t size = 0; size < n_products; size++) {
        auto attrs = l2_files[0].get_product_attributes(product_list[size]);
        netCDF::NcVar var_nc = binned_data.getVar(output_prod_names[size]);
        for (auto &[attr_name, attr_val] : attrs) {
            var_nc.putAtt(attr_name, netCDF::ncFloat, attr_val);
        }
    }

    // Define bin index variables
    if (is64bit)
        status = defineBinIndex64_nc(input.deflate, binned_data.getId());
    else
        status = defineBinIndex_nc(input.deflate, binned_data.getId());
    if (status) {
        printf("-E- %s:%d Could not define binIndex variable in output file\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }

    // Define quality variables if needed
    if (read_quality) {
        status = defineQuality_nc(input.deflate, binned_data.getId());
        if (status) {
            printf("-E- %s:%d Could not define quality variable in output file\n", __FILE__, __LINE__);
            exit(EXIT_FAILURE);
        }
    }

    // Initialize bin index arrays, for beggning and extent of the bins
    std::vector<uint64_t> beg(nrows);
    std::vector<uint32_t> ext(nrows);
    std::vector<binIndexStruct64_nc> binIndex64nc;
    std::vector<binIndexStruct_nc> binIndex32nc;
    int ngroup{-1}, n_rows_in_group{-1};
    ini_bin_index_arrays(input, nrows, is64bit, shape, basebin, beg, ext, binIndex64nc, binIndex32nc, ngroup,
                         n_rows_in_group);

    // Tracking time and skip bad float variables
    time_t tnow;
    tm *tmnow;
    size_t total_skip_bad_float = 0;

    // Main binning loop over sinusoidal grid rows
    // krow is index of the row.
    for (int krow = 0, igroup = 0; igroup < ngroup; igroup++) {
        // Skip rows outside lat bounds
        if (shape->row2lat(krow + n_rows_in_group) < input.latsouth) {
            krow += n_rows_in_group;
            continue;
        }
        if (shape->row2lat(krow) > input.latnorth) {
            krow += n_rows_in_group;
            continue;
        }

        // Print progress info
        {
            time(&tnow);
            tmnow = localtime(&tnow);
            printf("krow:%6d out of %6d  (%6.2f to %6.2f) ", krow, nrows,
                   static_cast<float>(krow) / nrows * 180 - 90,
                   static_cast<float>(krow + n_rows_in_group) / nrows * 180 - 90);
            printf("%s", asctime(tmnow));
        }

        bool within_flag = false;

        // Set what rows of the sinusoidal grid a give file covers
        for (size_t ifile = 0; ifile < nfiles; ifile++) {
            size_t nrec, nsamp;
            l2_files[ifile].getDimensions(nrec, nsamp);
            scan_in_rowgroup[ifile] = std::vector<u_char>(nrec + 1, 0);

            // Mark scans within row group
            for (size_t iscan = 0; iscan < nrec; iscan++) {
                scan_in_rowgroup[ifile][iscan] = scan_within_the_group;
                if (bscan_row[ifile][iscan] < krow ||
                    escan_row[ifile][iscan] >= (krow + n_rows_in_group - 1)) {
                    scan_in_rowgroup[ifile][iscan] = scan_not_within_the_group;
                }
            }

            // Check if any scans are within group
            if (!within_flag) {
                for (size_t iscan = 0; iscan < nrec; iscan++) {
                    if (scan_in_rowgroup[ifile][iscan] == scan_within_the_group) {
                        within_flag = true;
                        break;
                    }
                }
            }
        }

        // Skip if no scans in group
        if (!within_flag) {
            krow += n_rows_in_group;
            continue;
        }

        // Allocate bin arrays
        uint64_t n_bins_in_group = basebin[krow + n_rows_in_group] - basebin[krow];
        size_t n_filled_bins = 0;
        size_t n_allocperbin{0};
        // memory info debug print
        if (input.verbose)
            printMemoryInfo("Before bin Allocation\n");
        L2BinStruct l2binStruct{n_bins_in_group, basebin};
        if (input.verbose)
            printMemoryInfo("After bin Allocation\n");

        // Calculate allocation per bin
        for (size_t ifile = 0; ifile < nfiles; ifile++) {
            size_t nrec, nsamp;
            l2_files[ifile].getDimensions(nrec, nsamp);
            n_allocperbin =
                std::max(n_active_files * nrec * nsamp / bin_allocator_normalization, n_allocperbin);
            if (n_allocperbin < min_per_bin_allocation)
                n_allocperbin = min_per_bin_allocation;
            if (n_allocperbin > max_per_bin_allocation)
                n_allocperbin = max_per_bin_allocation;
        }
        if (input.verbose)
            printf("%-20s:%8zu\n\n", "# allocated per bin", n_allocperbin);

        // Bin each input file
        for (size_t ifile = 0; ifile < nfiles; ifile++) {
            if (brk_scan[ifile] == brake_scan_fill_value)
                continue;

            // Skip if no scans in row group
            const std::vector<u_char> &vec = scan_in_rowgroup[ifile];
            if (std::find(vec.begin(), vec.end(), scan_within_the_group) == vec.end())
                continue;

            size_t nrec, nsamp;
            l2_files[ifile].getDimensions(nrec, nsamp);
            std::vector<float> min_values = l2_files[ifile].get_min_value_product();
            std::vector<float> max_values = l2_files[ifile].get_max_value_product();

            // Calculate date differences
            if (!regional) {
                auto date = (time_t)l2_files[ifile].get_start_time() / 86400;
                diffday_beg = date - startdate;
                diffday_end = date - enddate;
            }

            // Track I/O
            size_t io_reads_count = 0;
            l2_files[ifile].reopenL2();

            // Process each scan
            for (size_t iscan = 0; iscan < nrec; iscan++) {
                if (scan_in_rowgroup[ifile][iscan] != scan_within_the_group)
                    continue;

                // Data arrays
                std::vector<float *> l2_data(n_products);
                float *latitude{nullptr};
                float *longitude{nullptr};
                int32_t *l2_flags{nullptr};
                int32_t *qual_flags{nullptr};

                // Set up area weighting
                area_weightings[ifile].set_corners(iscan, scan_in_rowgroup[ifile]);

                // Read data
                l2_files[ifile].readL2dataScan(l2_data, iscan, scan_in_rowgroup[ifile]);
                l2_files[ifile].readLatitudeScan(&latitude, iscan, scan_in_rowgroup[ifile]);
                l2_files[ifile].readLongitudeScan(&longitude, iscan, scan_in_rowgroup[ifile]);
                if (read_l2_flags) {
                    l2_files[ifile].readL2FlagsScan(&l2_flags, iscan, scan_in_rowgroup[ifile]);
                }
                if (read_quality) {
                    l2_files[ifile].readQualityFlagsScan(&qual_flags, iscan, scan_in_rowgroup[ifile]);
                }

                // Print I/O stats
                size_t total_reads = l2_files[ifile].get_total_io_reads();
                if (io_reads_count != total_reads && input.verbose) {
                    std::cout << "Total reads from the file " << total_reads << " for krow " << krow << " at iscan=" << iscan
                              << std::endl;
                    io_reads_count = total_reads;
                }

                // Check dateline crossing
                int scan_crosses_dateline = 0;
                float scan_lon = BAD_FLT;
                if (brk_scan[ifile] != 0) {
                    int npixls = nsamp - 1;
                    float slon = BAD_FLT;
                    float elon = BAD_FLT;

                    for (int i_p = 0; i_p <= npixls; i_p++) {
                        if (std::abs(slon) <= 180.0)
                            break;
                        slon = longitude[i_p];
                    }
                    for (int i_p = npixls; i_p >= 0; i_p--) {
                        if (std::abs(elon) <= 180.0)
                            break;
                        elon = longitude[i_p];
                    }
                    if (std::abs(elon) > 180.0 || std::abs(slon) > 180.0)
                        continue;
                    scan_lon = elon;
                    if (slon * elon < 0)
                        scan_crosses_dateline = 1;
                }

                // Skip based on dateline
                if (scan_crosses_dateline == 0 && !regional) {
                    if (skip_DL(scan_lon, brk_scan[ifile], input.night, diffday_end, diffday_beg))
                        continue;
                }

                // Print progress
                if ((iscan % 100) == 0 && input.verbose) {
                    printf("ifile:%4zu  iscan:%6zu  nsamp:%8zu\n", ifile, iscan, nsamp);
                }

                // Process each pixel
                for (size_t ipixl = 0; ipixl < nsamp; ipixl++) {
                    // Check L2 flags
                    if (read_l2_flags) {
                        int32_t flag_check = l2_flags[ipixl];
                        if ((flag_check & flagusemask[ifile]) != 0)
                            continue;
                        if ((flag_check & required[ifile]) != required[ifile])
                            continue;
                    }

                    // Check dateline crossing
                    if (scan_crosses_dateline == 1 && !regional) {
                        if (skip_DL(longitude[ipixl], brk_scan[ifile], input.night, diffday_end, diffday_beg))
                            continue;
                    }

                    // Check for bad values
                    bool bad_value = false;
                    for (size_t iprod = 0; iprod < n_products; iprod++) {
                        float f32 = l2_data[iprod][ipixl];
                        if (f32 == BAD_FLT || min_values[iprod] > f32 || max_values[iprod] < f32) {
                            bad_value = true;
                            total_skip_bad_float++;
                            break;
                        }
                    }
                    if (bad_value)
                        continue;

                    // Check longitude bounds
                    if (input.lonwest != 0.0 || input.loneast != 0.0) {
                        if (longitude[ipixl] < input.lonwest)
                            continue;
                        if (longitude[ipixl] > input.loneast)
                            continue;
                    }

                    // Check latitude bounds
                    if (latitude[ipixl] < input.latsouth)
                        continue;
                    if (latitude[ipixl] > input.latnorth)
                        continue;
                    if (!valid_lat(latitude[ipixl]))
                        continue;
                    if (!valid_lon(longitude[ipixl]))
                        continue;
                    if (!area_weightings[ifile].valid_geolocation(ipixl))
                        continue;

                    // Process pixel based on area weighting mode
                    if (area_weightings[ifile].get_mode()) {
                        // Get bins and areas that intersect pixel
                        std::map<uint64_t, double> areas;
                        getBins(shape, area_weightings[ifile], ipixl, areas);

                        // Add pixel data to each intersecting bin
                        for (auto const &[bin,area] : areas) {
                            addPixelToBin(l2_data, qual_flags, l2binStruct, ifile, krow, n_rows_in_group,
                                          n_bins_in_group, n_allocperbin, n_products, ipixl, bin,
                                          area);
                        }
                    } else {
                        // Get single bin for pixel
                        uint64_t bin = getbinnum(shape, latitude, longitude, ipixl);

                        // Check if composite is set and if the bin should be replaced/filled
                        bool addPixelStatus = true;
                        if (use_composite) {
                            size_t comp_index = composite_product_index[ifile];
                            addPixelStatus =
                                updateBinCompositeScheme(l2_data[comp_index], l2binStruct, comp_index,
                                                        composite_scheme, krow, n_products, ipixl, bin);
                        }

                        // Add pixel to bin
                        if (addPixelStatus)
                            addPixelToBin(l2_data, qual_flags, l2binStruct, ifile, krow, n_rows_in_group,
                                          n_bins_in_group, n_allocperbin, n_products, ipixl, bin);
                    }
                }
            }
            //  Reset file cache
            l2_files[ifile].reset_cache();
            // Memory info debug print
            if (input.verbose) {
                std::cout << "#file " << ifile << " has been processed for krow= " << krow << std::endl;
                printMemoryInfo("Memory footprint after the cache reset \n");
            }
        }

        // Print progress
        if (input.verbose) {
            time(&tnow);
            tmnow = localtime(&tnow);
            printf("krow:%5d After data_value fill: %s\n", krow, asctime(tmnow));
        }

        // Count filled bins
        for (size_t ibin = 0; ibin < n_bins_in_group; ibin++) {
            if (l2binStruct.nobs[ibin] > 0 && l2binStruct.nobs[ibin] < input.minobs)
                l2binStruct.nobs[ibin] = 0;

            if (l2binStruct.nobs[ibin] != 0)
                n_filled_bins++;
        }

        // Process filled bins
        std::vector<binListStruct64_nc> binList64nc;
        std::vector<binListStruct_nc> binList32nc;
        if (n_filled_bins > 0) {
            // Fill bins and get quality flags
            std::vector<uint8_t> best_qual = fill_bins(
                read_quality, l2binStruct, shape, binList64nc, binList32nc, n_filled_bins, n_bins_in_group,
                n_products, krow, nfiles, qual_max_allowed, is64bit, geospatial_bounds);

            if (n_filled_bins > 0) {
                printf("%-20s:%8llu\n", "# bins in row group", (unsigned long long)n_bins_in_group);
                printf("%-20s:%8zu\n", "# filled bins", n_filled_bins);

                // Write bin data
                if (is64bit) {
                    writeBinList_nc(binned_data.getId(), n_filled_bins, (VOIDP)binList64nc.data());
                } else {
                    writeBinList_nc(binned_data.getId(), n_filled_bins, (VOIDP)binList32nc.data());
                }

                fill_data(binned_data, l2binStruct, n_filled_bins, n_bins_in_group, n_products);

                // Write quality data
                if (read_quality) {
                    fill_quality_data(binned_data, l2binStruct, n_filled_bins, n_bins_in_group, best_qual);
                }

                // Update bin indexes
                update_bin_index(l2binStruct, shape, is64bit, n_filled_bins, n_bins_in_group, nrows, krow,
                                 n_rows_in_group, total_filled_bins, beg, ext, binIndex64nc, binIndex32nc);

                if (input.verbose) {
                    time(&tnow);
                    tmnow = localtime(&tnow);
                    printf("krow:%5d After bin processing:  %s", krow, asctime(tmnow));
                }
            }
        }
        krow += n_rows_in_group;

        // Print diagnostics
        if (input.verbose)
            std::cout << "Total bad points skipped " << total_skip_bad_float <<" for krow=" << krow <<"\n";
    }

    // Print final progress
    if (input.verbose) {
        time(&tnow);
        tmnow = localtime(&tnow);
        printf("krow:%5d %s", n_rows_in_group * ngroup, asctime(tmnow));
    }
    printf("total_filled_bins: %ld\n", (long int)total_filled_bins);
    // release memory for shape
    delete shape;
    // Exit if no bins filled
    if (total_filled_bins == 0) {
        l3_file.close();
        std::remove(input.ofile);
        int ret_status = 110;
        std::cout << "rm -rf " << ofile << std::endl;
        exit(ret_status);
    }

    // Write final bin index data
    if (is64bit)
        writeBinIndex_nc(binned_data.getId(), nrows, binIndex64nc.data());
    else
        writeBinIndex_nc(binned_data.getId(), nrows, binIndex32nc.data());

    // Write metadata and finish
    write_l3_metadata(basebin, l3_file, meta_l3b, input, nfiles, nrows, total_filled_bins, l2_files,
                      product_list, argc, argv, brk_scan, is64bit, geospatial_bounds);
    return EXIT_SUCCESS;
}