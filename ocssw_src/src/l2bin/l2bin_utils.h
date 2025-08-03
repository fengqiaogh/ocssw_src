#ifndef L2BIN_UTILS
#define L2BIN_UTILS
#include "l2_utils.hpp"
#include "area_weighting.h"
#include "L3Shape.h"
#include "L3ShapeIsine.h"
#include <ncdfbin_utils.h>
#include "l2_reader.hpp"
#include <seaproto.h>
#include "sensorInfo.h"
#include <meta_l3b.h>
#include "l2bin_input.h"
#include <timeutils.h>
#include <genutils.h>
#include <stdlib.h>
#include <fstream>
#include "get_dataday.hpp"
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/foreach.hpp>
namespace bg = boost::geometry;
typedef bg::model::point<double, 2, bg::cs::geographic<bg::degree>> Point_t;
typedef bg::model::polygon<Point_t> Polygon_t;
typedef bg::model::box<Point_t> Box_t;

/**
 * @brief Structure for managing binned L2 data arrays using raw pointers
 * @details Uses raw pointers instead of vectors due to heap fragmentation issues observed with gcc 9-12.
 * Contains arrays for bin quality, scene counts, file indices, observation counts, data values and areas.
 */
struct L2BinStruct {
    size_t number_of_bins;
    int16_t *qual{nullptr};                   ///< Quality flags for each bin
    int16_t *nscenes{nullptr};                ///< Number of scenes contributing to each bin
    int16_t *lastfile{nullptr};               ///< Index of last file contributing to each bin
    int16_t *nobs{nullptr};                   ///< Number of observations in each bin
    int16_t *allocated_space{nullptr};        ///< Allocated space for each bin
    float **data_values{nullptr};             ///< 2D array of data values [bin][product]
    double **data_areas{nullptr};             ///< Area weights for each observation
    int16_t **file_index{nullptr};            ///< File indices for each observation
    uint8_t **data_quality{nullptr};          ///< Quality flags for each observation
    std::vector<uint64_t> *basebin{nullptr};  ///< External base bin numbers

    /** @brief Default constructor */
    L2BinStruct() = default;

    /**
     * @brief Allocates memory for a bin if not already allocated
     * @param n_allocperbin Number of observations to allocate per bin
     * @param l3b_nprod Number of L3 products
     * @param ibin Bin index to allocate
     */
    void allocate(size_t n_allocperbin, size_t l3b_nprod, size_t ibin) {
        if (file_index[ibin] == nullptr) {
            calloc(&file_index[ibin], n_allocperbin);
            calloc(&data_values[ibin], n_allocperbin * l3b_nprod);
            calloc(&data_areas[ibin], n_allocperbin);
            calloc(&data_quality[ibin], n_allocperbin);
            allocated_space[ibin] = n_allocperbin;
        }
    }

    /**
     * @brief Reallocates memory for a bin to accommodate more observations
     * @param n_allocperbin Additional observations to allocate
     * @param l3b_nprod Number of L3 products
     * @param ibin Bin index to reallocate
     */
    void reallocate(size_t n_allocperbin, size_t l3b_nprod, size_t ibin) {
        allocated_space[ibin] += n_allocperbin;
        file_index[ibin] = realloc(&file_index[ibin], nobs[ibin] + n_allocperbin);
        data_values[ibin] = realloc(&data_values[ibin], (nobs[ibin] + n_allocperbin) * l3b_nprod);
        data_quality[ibin] = realloc(&data_quality[ibin], nobs[ibin] + n_allocperbin);
        data_areas[ibin] = realloc(&data_areas[ibin], nobs[ibin] + n_allocperbin);
    }

    /**
     * @brief Constructs and initializes arrays for given number of bins
     * @param number_of_bins Number of bins to allocate
     * @param basebin Reference to vector of base bin numbers
     */
    L2BinStruct(size_t number_of_bins, std::vector<uint64_t> &basebin)
        : number_of_bins(number_of_bins), basebin(&basebin) {
        calloc(&qual, number_of_bins);
        calloc(&nscenes, number_of_bins);
        calloc(&lastfile, number_of_bins);
        calloc(&nobs, number_of_bins);
        calloc(&allocated_space, number_of_bins);
        calloc(&data_values, number_of_bins);
        calloc(&data_areas, number_of_bins);
        calloc(&file_index, number_of_bins);
        calloc(&data_quality, number_of_bins);
        for (size_t i = 0; i < number_of_bins; i++) {
            lastfile[i] = -1;
            qual[i] = 3;
        }
    }

    /**
     * @brief Destructor that frees all allocated memory
     */
    ~L2BinStruct() {
        if (qual != nullptr) {
            free(qual);
        }
        if (nscenes != nullptr) {
            free(nscenes);
        }
        if (lastfile != nullptr) {
            free(lastfile);
        }
        if (nobs != nullptr) {
            free(nobs);
        }
        if (allocated_space != nullptr) {
            free(allocated_space);
        }
        if (data_values != nullptr) {
            for (size_t i = 0; i < number_of_bins; i++) {
                if (data_values[i] != nullptr) {
                    free(data_values[i]);
                }
                if (data_areas[i] != nullptr) {
                    free(data_areas[i]);
                }
                if (file_index[i] != nullptr) {
                    free(file_index[i]);
                }
                if (data_quality[i] != nullptr) {
                    free(data_quality[i]);
                }
            }
            free(data_values);
            free(data_areas);
            free(file_index);
            free(data_quality);
        }
    }
};

/**
 * @brief Checks if a file exists on the filesystem
 *
 * @param name Path to the file to check
 * @return true if file exists
 * @return false if file does not exist
 */
bool exists_test(const std::string &name);

/**
 * @brief Determines if a longitude point should be skipped based on dateline crossing criteria
 *
 * @param lon Longitude value to check
 * @param side Which side of the dateline to include (east or west)
 * @param night_flag Flag indicating if this is night data
 * @param end_day Time difference between scan start time and end day
 * @param beg_day Time difference between scan start time and beginning day
 * @return true if point should be skipped
 * @return false if point should be included
 */
bool skip_DL(float lon, int side, int night_flag, time_t end_day, time_t beg_day);

/**
 * @brief Parses and returns output product name mappings
 *
 * @param oprodname String containing output product name mappings
 * @return std::map<std::string, std::string> Map of input product names to output product names
 */
std::map<std::string, std::string> get_oprodname(const std::string &oprodname);

/**
 * @brief Sets the scan break points based on equator crossing criteria
 *
 * @param brk_scan Output parameter to store break scan indices
 * @param l2file L2 file reader object
 * @param deltaeqcross Maximum allowed equator crossing longitude difference
 * @param night Flag indicating if processing night data
 * @param startdate Start date to consider
 * @param enddate End date to consider
 */
void set_breakscan(int &brk_scan, L2_Reader &l2file, float deltaeqcross, int night, int startdate,
                   int enddate);

/**
 * @brief Sets up product names in the netCDF output file
 *
 * @param deflate Compression level for netCDF variables
 * @param product_list List of product names to process
 * @param output_l3_filenames Map of output L3 filenames
 * @param binned_data NetCDF group to store binned data
 * @return std::vector<std::string> List of configured product names
 */
std::vector<std::string> set_nc_prodnames(int deflate, std::vector<std::string> &product_list,
                                          std::map<std::string, std::string> &output_l3_filenames,
                                          netCDF::NcGroup &binned_data);

/**
 * @brief Sets up bin indexing structures for a row
 *
 * @param i Row index
 * @param is64bit Flag indicating if using 64-bit bin indices
 * @param shape L3 binning shape object
 * @param basebin Vector of base bin numbers
 * @param beg Vector of beginning indices
 * @param ext Vector of bin extents
 * @param binIndex64nc Vector of 64-bit bin index structures
 * @param binIndex32nc Vector of 32-bit bin index structures
 */
void set_bin_index(int i, bool is64bit, l3::L3ShapeIsine *shape, std::vector<uint64_t> &basebin,
                   std::vector<uint64_t> &beg, std::vector<uint32_t> &ext,
                   std::vector<binIndexStruct64_nc> &binIndex64nc,
                   std::vector<binIndexStruct_nc> &binIndex32nc);

/**
 * @brief Initializes bin index arrays for the entire grid
 *
 * @param input Input parameters structure
 * @param nrows Number of rows in grid
 * @param is64bit Flag indicating if using 64-bit bin indices
 * @param shape L3 binning shape object
 * @param basebin Vector of base bin numbers
 * @param beg Vector of beginning indices
 * @param ext Vector of bin extents
 * @param binIndex64nc Vector of 64-bit bin index structures
 * @param binIndex32nc Vector of 32-bit bin index structures
 * @param ngroup Output parameter for number of groups
 * @param n_rows_in_group Output parameter for number of rows per group
 */
void ini_bin_index_arrays(instr &input, int nrows, bool is64bit, l3::L3ShapeIsine *shape,
                          std::vector<uint64_t> &basebin, std::vector<uint64_t> &beg,
                          std::vector<uint32_t> &ext, std::vector<binIndexStruct64_nc> &binIndex64nc,
                          std::vector<binIndexStruct_nc> &binIndex32nc, int &ngroup, int &n_rows_in_group);

/**
 * @brief Prints current memory usage statistics
 *
 * @param message Optional message to print with memory info
 */
void printMemoryInfo(const std::string &message = "");

/**
 * @brief Gets bin numbers and areas for a pixel
 *
 * @param shape L3 binning shape object
 * @param areaWeight Area weighting calculator
 * @param ipixl Pixel index
 * @param areas Output map of bin numbers to areas
 */
void getBins(l3::L3ShapeIsine *shape, AreaWeighting &areaWeight, size_t ipixl,
             std::map<uint64_t, double> &areas);

/**
 * @brief Adds pixel data to a bin
 *
 * @param l2data Vector of L2 data arrays
 * @param qual_flags Quality flags array
 * @param l2binStruct Binning structure
 * @param ifile File index
 * @param krow Row index
 * @param n_rows_in_group Number of rows per group
 * @param n_bins_in_group Number of bins per group
 * @param n_allocperbin Number of allocations per bin
 * @param l3b_nprod Number of L3 products
 * @param ipixl Pixel index
 * @param bin Bin number
 * @param areaFrac Area fraction for this pixel (default 1.0)
 */
void addPixelToBin(std::vector<float *> &l2data, int32_t *qual_flags, L2BinStruct &l2binStruct, size_t ifile,
                   size_t krow, size_t n_rows_in_group, size_t n_bins_in_group, size_t n_allocperbin,
                   size_t l3b_nprod, size_t ipixl, uint64_t bin, double areaFrac = 1.0);

/**
 * @brief Gets bin number for a lat/lon point
 *
 * @param shape L3 binning shape object
 * @param latitude Pointer to latitude value
 * @param longitude Pointer to longitude value
 * @param ipixl Pixel index
 * @return int64_t Bin number
 */
int64_t getbinnum(l3::L3ShapeIsine *shape, float *latitude, float *longitude, size_t ipixl);

/**
 * @brief Fills bins with data for current group
 *
 * @param read_quality Flag indicating if quality data should be read
 * @param l2binStruct Binning structure
 * @param shape L3 binning shape object
 * @param binList64nc Vector of 64-bit bin list structures
 * @param binList32nc Vector of 32-bit bin list structures
 * @param n_filled_bins Output parameter for number of filled bins
 * @param n_bins_in_group Number of bins in group
 * @param n_products Number of products
 * @param krow Row index
 * @param nfiles Number of input files
 * @param qual_max_allowed Maximum allowed quality value
 * @param is64bit Flag indicating if using 64-bit bins
 * @param bounds Geospatial bounds array
 * @return std::vector<uint8_t> Vector of quality values
 */
std::vector<uint8_t> fill_bins(bool read_quality, L2BinStruct &l2binStruct, l3::L3ShapeIsine *shape,
                               std::vector<binListStruct64_nc> &binList64nc,
                               std::vector<binListStruct_nc> &binList32nc, size_t &n_filled_bins,
                               size_t n_bins_in_group, size_t n_products, size_t krow, size_t nfiles,
                               size_t qual_max_allowed, bool is64bit, float *bounds);

/**
 * @brief Fills binned data into netCDF output
 *
 * @param binned_data NetCDF group for binned data
 * @param l2binStruct Binning structure
 * @param n_filled_bins Number of filled bins
 * @param n_bins_in_group Number of bins in group
 * @param n_products Number of products
 */
void fill_data(netCDF::NcGroup &binned_data, L2BinStruct &l2binStruct, size_t n_filled_bins,
               size_t n_bins_in_group, size_t n_products);

/**
 * @brief Fills quality data into netCDF output
 *
 * @param binned_data NetCDF group for binned data
 * @param l2binStruct Binning structure
 * @param n_filled_bins Number of filled bins
 * @param n_bins_in_group Number of bins in group
 * @param best_qual Vector of best quality values
 */
void fill_quality_data(netCDF::NcGroup &binned_data, L2BinStruct &l2binStruct, size_t n_filled_bins,
                       size_t n_bins_in_group, std::vector<uint8_t> &best_qual);

/**
 * @brief Updates bin index arrays after filling bins
 *
 * @param l2binStruct Binning structure
 * @param shape L3 binning shape object
 * @param is64bit Flag indicating if using 64-bit bins
 * @param n_filled_bins Number of filled bins
 * @param n_bins_in_group Number of bins in group
 * @param nrows Number of rows
 * @param krow Current row index
 * @param n_rows_in_group Number of rows per group
 * @param total_filled_bins Running total of filled bins
 * @param beg Vector of beginning indices
 * @param ext Vector of bin extents
 * @param binIndex64nc Vector of 64-bit bin index structures
 * @param binIndex32nc Vector of 32-bit bin index structures
 */
void update_bin_index(L2BinStruct &l2binStruct, l3::L3ShapeIsine *shape, bool is64bit, size_t n_filled_bins,
                      size_t n_bins_in_group, size_t nrows, size_t krow, size_t n_rows_in_group,
                      int64_t &total_filled_bins, std::vector<uint64_t> &beg, std::vector<uint32_t> &ext,
                      std::vector<binIndexStruct64_nc> &binIndex64nc,
                      std::vector<binIndexStruct_nc> &binIndex32nc);

/**
 * @brief Writes L3 metadata to output file
 *
 * @param basebin Vector of base bin numbers
 * @param ofile Output netCDF file
 * @param meta_l3 L3 metadata structure
 * @param input Input parameters
 * @param nfiles Number of input files
 * @param nrows Number of rows
 * @param total_filled_bins Total number of filled bins
 * @param l2_files Vector of L2 file readers
 * @param product_list List of products
 * @param argc Command line argument count
 * @param argv Command line arguments
 * @param brk_scan Vector of scan break points
 * @param is64bit Flag indicating if using 64-bit bins
 * @param geospatial_bounds Array of geospatial bounds
 */
void write_l3_metadata(std::vector<uint64_t> &basebin, netCDF::NcFile &ofile, meta_l3bType &meta_l3,
                       instr &input, size_t nfiles, size_t nrows, size_t total_filled_bins,
                       std::vector<L2_Reader> &l2_files, std::vector<std::string> &product_list, int argc,
                       char **argv, std::vector<int> &brk_scan, bool is64bit, float *geospatial_bounds);

/**
 * @brief Appends a composite product specification to the products list
 *
 * @param products_requested String of requested products
 * @param composite_product_name Name of composite product
 * @param composite_scheme Compositing scheme to use
 * @param min_max_scheme Output parameter for min/max scheme
 */
void append_composite_product(std::string &products_requested, const std::string &composite_product_name,
                              const std::string &composite_scheme, int &min_max_scheme);

/**
 * @brief Updates bin data according to composite scheme
 *
 * @param composite_data Composite data array
 * @param l2binStruct Binning structure
 * @param prod_index Product index
 * @param composite_scheme Compositing scheme to use
 * @param krow Row index
 * @param l3b_nprod Number of L3 products
 * @param ipixl Pixel index
 * @param bin Bin number
 * @return true if bin was updated
 * @return false if bin was not updated
 */
bool updateBinCompositeScheme(float *composite_data, L2BinStruct &l2binStruct, size_t prod_index,
                              int composite_scheme, size_t krow, size_t l3b_nprod, size_t ipixl,
                              uint64_t bin);

/**
 * @brief Checks if a bin intersects with a pixel box
 * 
 * @param shape L3 binning shape object
 * @param row Row index of bin
 * @param col Column index of bin
 * @param pixelBox Box representing the pixel
 * @param areaFrac Output parameter for area fraction of intersection
 * @return true if bin intersects pixel
 * @return false if no intersection
 */
bool binIntersectsPixel(l3::L3ShapeIsine *shape, int32_t row, int32_t col, Box_t &pixelBox, double &areaFrac);

/**
 * @brief Checks if a bin intersects with a pixel polygon
 * 
 * @param shape L3 binning shape object
 * @param row Row index of bin
 * @param col Column index of bin
 * @param pixelPoly Polygon representing the pixel
 * @param areaFrac Output parameter for area fraction of intersection
 * @return true if bin intersects pixel
 * @return false if no intersection
 */
bool binIntersectsPixel(l3::L3ShapeIsine *shape, int32_t row, int32_t col, Polygon_t &pixelPoly,
                        double &areaFrac);

/**
 * @brief Gets all bins that intersect a pixel in a given row
 * 
 * @tparam T Type of pixel geometry (Box_t or Polygon_t)
 * @param shape L3 binning shape object
 * @param lat Latitude of pixel center
 * @param lon Longitude of pixel center
 * @param pixelPoly Pixel geometry
 * @param areas Output map of bin numbers to intersection areas
 * @return true if any bins intersect
 * @return false if no bins intersect
 */
template <class T>
bool getBinsFromRow(l3::L3ShapeIsine *shape, double lat, double lon, T &pixelPoly,
                    std::map<uint64_t, double> &areas) {
    // Initialize variables for row/column tracking
    int32_t row0, col0;
    int32_t col;
    bool result = false;
    double areaFrac;
    uint64_t bin;

    // Convert lat/lon to row/col indices
    shape->latlon2rowcol(lat, lon, row0, col0);

    // Look right from center column
    col = col0;
    while (binIntersectsPixel(shape, row0, col, pixelPoly, areaFrac)) {
        // Found an intersecting bin
        result = true;
        // Convert row/col to bin number
        bin = shape->rowcol2bin(row0, col);
        // Store bin and its intersection area
        areas.emplace(bin, areaFrac);
        // Move to next column right
        col++;
        // Wrap around if needed
        shape->constrainRowCol(row0, col);
        // Stop if we've wrapped all the way around
        if (col == col0)
            break;
    }

    // Look left from center column
    col = col0 - 1;
    while (binIntersectsPixel(shape, row0, col, pixelPoly, areaFrac)) {
        // Found an intersecting bin
        result = true;
        // Convert row/col to bin number
        bin = shape->rowcol2bin(row0, col);
        // Store bin and its intersection area
        areas.emplace(bin, areaFrac);
        // Move to next column left
        col--;
        if(col < 0)
            break;
        // Wrap around if needed
        shape->constrainRowCol(row0, col);
        // Stop if we've wrapped all the way around
        if (col == col0 - 1)
            break;
    }

    return result;
}

/**
 * @brief Gets areas of all bins that intersect with a pixel
 * 
 * @tparam T Type of pixel geometry (Box_t or Polygon_t)
 * @param lat Current latitude being checked
 * @param lat0 Initial latitude of pixel center
 * @param lon0 Longitude of pixel center
 * @param deltaLat Latitude increment for searching
 * @param shape L3 binning shape object
 * @param polygon Pixel geometry
 * @param areas Output map of bin numbers to intersection areas
 */
template <typename T>
void get_areas(double lat, double lat0, double lon0, double deltaLat, l3::L3ShapeIsine *shape, T &polygon,
                  std::map<uint64_t, double> &areas) {
    // Start at initial latitude and scan upward toward North pole
    lat = lat0;
    while (getBinsFromRow(shape, lat, lon0, polygon, areas)) {
        // Increment latitude by bin spacing
        lat += deltaLat;
        // Stop if we've reached the North pole
        if (lat > 90.0)
            break;
    }

    // Start at one bin below initial latitude and scan toward South pole
    lat = lat0 - deltaLat;
    while (getBinsFromRow(shape, lat, lon0, polygon, areas)) {
        // Decrement latitude by bin spacing
        lat -= deltaLat;
        // Stop if we've reached the South pole
        if (lat < -90.0)
            break;
    }
}


/**
 * @brief Handles box rotation across the 180/-180 degree meridian
 * 
 * @param pixelBox Input box that may cross the meridian
 * @param boxes Output array to store rotated boxes - will contain 1 box if no meridian crossing,
 *              or 2 boxes if the input crosses the meridian
 * @return size_t Number of boxes in output array (1 or 2)
 */
size_t rotate_polygon(Box_t &pixelBox, std::array<Box_t,2> & boxes);

/**
 * @brief Handles polygon rotation across the 180/-180 degree meridian
 * 
 * @param pixelPoly Input polygon that may cross the meridian
 * @param polys Output array to store rotated polygons - will contain 1 polygon if no meridian crossing,
 *              or 2 polygons if the input crosses the meridian
 * @return size_t Number of polygons in output array (1 or 2)
 */
size_t rotate_polygon(Polygon_t &pixelPoly, std::array<Polygon_t,2> & polys);
#endif  // L2BIN_UTILS