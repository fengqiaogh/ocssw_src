#include "CacheSize.h"
#include <cmath>
#include "OutFileNetcdf4.h"
namespace {
/// Maximum number of elements in netCDF cache
constexpr int n_elems = 2003;
/// Cache preemption ratio (fraction of least recently used data to flush)
constexpr float preemption = .75;
/// The chunk size along the first dimension, number of lints
size_t chunk_along_lines = 16;
/// Maximum allowed cache size in bytes (1 GB)
constexpr size_t max_allowed_cache_size = 1024 * 1024 * 1024;
/// Minimum allowed cache size in bytes (16 MB)
constexpr size_t min_allowed_cache_size = 16 * 1024 * 1024;
/// Cache size for L3B files in bytes (1 MB)
constexpr size_t l3b_cache_size = 1024 * 1024;
}  // namespace

/**
 * @brief Sets the netCDF cache size and chunk size for writing L3 data files.
 *
 * Calculates the appropriate cache size based on the number of products and
 * file dimensions. If the calculated cache size exceeds the maximum allowed size,
 * adjusts the chunk_along_lines parameter to fit within the limit.
 *
 * @param[in] l3file Pointer to the L3File object of the input file
 * @param[in] file Pointer to the OutFile object of the primary output file
 * @param[in] file2 Pointer to optional secondary OutFile object (may be null)
 * @param[in] verbose Flag to enable informational console output
 *
 * @return true if the cache configuration was successfully set, false otherwise
 */
bool set_cache_size_chunk_size_write(l3::L3File* l3file, OutFile* file, OutFile* file2, bool verbose) {
    // check if file can be dynamic casted to a derived class OutFileNetcdf4 which is derived from OutFile
    OutFileNetcdf4* ncfile = dynamic_cast<OutFileNetcdf4*>(file);
    // if file 2 is not null, check if it can be dynamic casted to a derived class OutFileNetcdf4
    if (file2) {
        OutFileNetcdf4* ncfile2 = dynamic_cast<OutFileNetcdf4*>(file2);
        // if both files are not netCDF, return true
        if (!ncfile && !ncfile2) {
            return true;
        }
    } else {
        // if only file1 is not netCDF, return true
        if (!ncfile) {
            return true;
        }
    }
    // Retrieve the number of active products (l3b products, wavelength_3d expanded)
    size_t number_of_products = l3file->getNumActiveProducts();

    // Get the width from the primary file and use the maximum if secondary file exists
    size_t width = file->getWidth();
    if (file2)
        width = std::max(width, (size_t)file2->getWidth());

    // Calculate initial cache size: products * width * lines_per_chunk * 4 bytes (float size)
    size_t cache_size = number_of_products * width * chunk_along_lines * 4;

    // If calculated cache size exceeds the maximum allowed, reduce chunk size and use max cache
    if (cache_size > max_allowed_cache_size) {
        // Adjust chunk_along_lines to fit within the maximum cache size constraint
        chunk_along_lines =
            std::max(1ul, max_allowed_cache_size / (number_of_products * width * 4));
        if (verbose) {
            std::cout << "-Info-: Setting cache size for writing to " << max_allowed_cache_size / 1024 / 1024 << "Mb"
                      << std::endl;
            cache_size = number_of_products * width * chunk_along_lines * 4;
            if (cache_size > max_allowed_cache_size) {
                fprintf(stdout,
                        "-Warning-: The cache size needed for efficient netCDF/HDF5 caching = ( %lu bytes) exceeds the maximum "
                        "cache size\n",
                        cache_size);
            }
        }
        return nc_set_chunk_cache(max_allowed_cache_size, n_elems, preemption) == NC_NOERR;
    }

    // Round cache size up to the next power of 2
    cache_size = std::pow(2, int(std::log2(cache_size)) + 1);

    // Log the cache size if verbose mode is enabled (in megabytes)
    if (verbose) {
        std::cout << "-Info-: Setting cache size for writing to " << cache_size / 1024 / 1024 << "Mb"
                  << std::endl;
    }

    // Apply the calculated cache configuration to netCDF
    return nc_set_chunk_cache(cache_size, n_elems, preemption) == NC_NOERR;
}

/**
 * @brief Sets the netCDF cache size for reading L3 data files based on file structure.
 *
 * Analyzes the L3 netCDF input file structure to determine the number of products and
 * maximum chunk size, then calculates an appropriate cache size for efficient reading.
 * For L3B files without lat/lon coordinates, uses a predefined L3B cache size.
 *
 * @param[in] fileName Path to the netCDF file to analyze
 * @param[in] verbose Flag to enable informational console output
 *
 * @return true if the cache configuration was successfully set, false otherwise
 */
bool set_cache_size_chunk_size_read(const std::string& fileName, bool verbose) {
    // Open the netCDF file in read-only mode
    netCDF::NcFile file(fileName, netCDF::NcFile::read);

    // Retrieve all variables from the file
    auto vars = file.getVars();

    // Check if file is L3M SMI format (must have lat/lon coordinates)
    // If not, use predefined L3B cache size and return
    if (vars.count("lat") == 0 || vars.count("lon") == 0) {
        file.close();
        return nc_set_chunk_cache(l3b_cache_size, n_elems, preemption) == NC_NOERR;
    }

    // Get longitude variable to determine grid width
    auto lon_var = file.getVar("lon");
    // Extract width from the first dimension of lon
    size_t width = lon_var.getDims().at(0).getSize();

    // Initialize variables for calculating cache size
    size_t max_chunk = 0;           // Maximum chunk size across all variables
    size_t number_of_products = 0;  // Total number of products (accounting for 3D layers)

    // Iterate through all variables to gather dimension and chunking information
    for (const auto& pair : vars) {
        const auto& var = pair.second;
        auto dims = var.getDims();

        // Skip 1D variables (coordinate variables)
        if (dims.size() == 1)
            continue;

        // Count 2D variables as single products
        if (dims.size() == 2)
            number_of_products++;

        // Count 3D variables: each layer in the third dimension is a separate product
        if (dims.size() == 3)
            number_of_products += dims.at(2).getSize();

        // Retrieve chunking information for this variable
        std::vector<size_t> chunkSizes;
        netCDF::NcVar::ChunkMode chunkMode;
        var.getChunkingParameters(chunkMode, chunkSizes);

        // Track the maximum chunk size along the first dimension (lines)
        if (chunkMode == netCDF::NcVar::nc_CHUNKED) {
            max_chunk = std::max(max_chunk, chunkSizes.at(0));
        }
    }

    // Calculate cache size: products * width * max_chunk * 4 bytes (float size)
    size_t cache_size = number_of_products * width * max_chunk * 4;

    // Close the file after analysis
    file.close();

    // Round cache size up to the next power of 2
    cache_size = std::pow(2, int(std::log2(cache_size)) + 1);

    // Ensure cache size meets minimum requirement
    cache_size = std::max(cache_size, min_allowed_cache_size);

    // Log the cache size if verbose mode is enabled (in megabytes)
    if (verbose) {
        std::cout << "-Info-:  Setting cache size for reading to " << cache_size / 1024 / 1024 << "Mb"
                  << std::endl;
    }

    // Apply the calculated cache configuration to netCDF
    return nc_set_chunk_cache(cache_size, n_elems, preemption) == NC_NOERR;
}

/**
 * @brief Retrieves the current chunk size along the line/first dimension.
 *
 * Returns the dynamically adjusted chunk_along_lines value, which may have been
 * modified during cache size calculations to fit within cache constraints.
 *
 * @return The current number of lines per chunk
 */
size_t get_chunk_along_lines() {
    // Return the current chunk size (may have been adjusted by write function)
    return chunk_along_lines;
}
