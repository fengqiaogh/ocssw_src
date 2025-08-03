#ifndef L2_READER_HPP
#define L2_READER_HPP
#include <string>
#include <netcdf>
#include "l2_utils.hpp"
#include <vector>
#include <unordered_map>
#include <map>
#include "l2_variable.hpp"
#include "l2_metadata.hpp"
#include "get_geospatial.hpp"
#include <optional>
#include <type_traits>
#include <memory>
#include <stdint.h>
#include <array>
typedef std::optional<std::array<size_t, 2>> cache_bounds;

class L2_Reader {
    std::string file_path;
    netCDF::NcFile file_nc;
    ScaledNcVar lat_nc, lon_nc;
    netCDF::NcVar l2_flags_nc, quality_flags_nc;
    size_t first_dimension;
    size_t second_dimension;
    size_t padding{10};
    bool use_cache{false};
    // attributes
    netCDF::NcVarAtt l2_mask_attr, l2_meaning_attr;
    std::vector<std::vector<float>> data_products_cached{};
    std::vector<std::string> products_list{};
    std::vector<float> latitude_cached{};
    std::vector<float> longitude_cached{};
    std::vector<int32_t> l2flags_cached{};
    std::vector<int32_t> qualityflags_cached{};
    std::vector<cache_bounds> start_cached_products;
    std::vector<cache_bounds> end_cached_products;
    cache_bounds start_cached_latitude;
    cache_bounds end_cached_latitude;
    cache_bounds start_cached_longitude;
    cache_bounds end_cached_longitude;
    cache_bounds start_cached_quality;
    cache_bounds end_cached_quality;
    cache_bounds start_cached_l2flags;
    cache_bounds end_cached_l2flags;
    std::vector<size_t> index{};
    std::vector<ScaledNcVar> l2_products{};
    std::vector<float> min_values{};
    std::vector<float> max_values{};
    std::vector<std::string> units{};
    std::string platform{}, instrument{};
    size_t number_of_l3_products;
    double start_time_unix{0};
    double end_time_unix{0};
    bool flag_l2_set = false;
    bool quality_flag_set = false;
    std::string qual_product_name;
    Geospatialbounds geo_bounds;
    std::vector<int32_t> l2_mask_bits;
    std::string l2_meaning;
    std::map<std::string, std::map<std::string, float>> product_attributes;
    // diagnostic variable
    size_t total_reads{0};
    std::unique_ptr<MetaL2> l2_metadata;
    /**
     * @brief Class to manage caching parameters for L2 file data reading
     * @details Handles the storage and calculation of cache boundaries and requested data ranges
     *          for both scan lines and pixels when reading L2 satellite data
     */
    class Cache_data {
       public:
        /** @brief Flag indicating if L2 file needs to be read */
        bool read_l2_file{false};
        /** @brief Beginning scan line of cached region */
        size_t bscan{};
        /** @brief Beginning pixel of cached region */
        size_t bpixel{};
        /** @brief Ending scan line of cached region */
        size_t escan{};
        /** @brief Ending pixel of cached region */
        size_t epixel{};
        /** @brief Beginning scan line requested by user */
        size_t bscan_requested{};
        /** @brief Ending scan line requested by user */
        size_t escan_requested{};
        /** @brief Beginning pixel requested by user */
        size_t bpixel_requested{};
        /** @brief Ending pixel requested by user */
        size_t epixel_requested{};
        /** @brief Starting indices for netCDF read */
        std::vector<size_t> start_nc;
        /** @brief Count of elements to read in each dimension */
        std::vector<size_t> count_nc;

        /**
         * @brief Constructs cache data object with specified boundaries
         * @param start Vector containing starting indices [scan, pixel]
         * @param count Vector containing counts for each dimension
         */
        Cache_data(std::vector<size_t>& start, std::vector<size_t>& count) {
            bscan_requested = start[0];
            escan_requested = start[0] + count[0];
            bpixel_requested = start[1];
            epixel_requested = start[1] + count[1];
            start_nc = start;
            count_nc = count;
        }

        /** @brief Default constructor */
        Cache_data() = default;
    };
    /**
     * @brief Sets up cached data parameters for reading a product
     * @param prodId Product ID to read
     * @param start Starting indices for reading
     * @param count Number of elements to read in each dimension
     * @param start_cached Start bounds of the cache
     * @param end_cached End bounds of the cache
     * @param mask Optional mask for selective reading
     * @return Cache_data object containing read parameters
     */
    Cache_data set_cached_data(int prodId, std::vector<size_t>& start, std::vector<size_t>& count,
                               cache_bounds& start_cached, cache_bounds& end_cached,
                               const std::vector<uint8_t>& mask = {});

    /**
     * @brief Copies data from cache to output array
     * @tparam T Type of output data
     * @param cache_data Cache parameters
     * @param cached_data Source data in cache
     * @param data Destination array
     */
    template <typename T>
    void copy_data(Cache_data& cache_data, T* cached_data, T* data) {
        size_t bpixel_requested = cache_data.bpixel_requested;
        size_t epixel_requested = cache_data.epixel_requested;
        size_t bscan_requested = cache_data.bscan_requested;
        size_t escan_requested = cache_data.escan_requested;
        size_t bpixel = cache_data.bpixel;
        size_t epixel = cache_data.epixel;
        size_t bscan = cache_data.bscan;
        for (size_t iscan = bscan_requested - bscan; iscan < escan_requested - bscan; iscan++) {
            std::copy(cached_data + (epixel - bpixel) * iscan + bpixel_requested - bpixel,
                      cached_data + (epixel - bpixel) * iscan + epixel_requested - bpixel,
                      data + (epixel_requested - bpixel_requested) * (iscan + bscan - bscan_requested));
        }
    };

    /**
     * @brief Copies data from cache to output pointer
     * @tparam T Type of output data
     * @param cache_data Cache parameters
     * @param cached_data Source data in cache
     * @param data Destination pointer
     */
    template <typename T>
    void copy_data(Cache_data& cache_data, T* cached_data, T** data) {
        size_t bscan_requested = cache_data.bscan_requested;
        size_t bscan = cache_data.bscan;
        *data = cached_data + second_dimension * (bscan_requested - bscan);
    };

    /**
     * @brief Copies data from cache to output pointer with copy flag
     * @tparam T Type of output pointer
     * @tparam U Type of cached data
     * @param cache_data Cache parameters
     * @param cached_data Source data in cache
     * @param data_ptr Destination pointer
     * @param copy Whether to copy data or just set pointer
     */
    template <typename T, typename U, std::enable_if_t<std::is_pointer<T>::value>* = nullptr>
    void copy_data(Cache_data& cache_data, U* cached_data, T* data_ptr, bool copy) {
        if (copy)
            copy_data(cache_data, cached_data, *data_ptr);  // Dereference the pointer
        else
            copy_data(cache_data, cached_data, data_ptr);  // Pass as-is
    }

    /**
     * @brief Copies data from cache to output pointer
     * @tparam T Type of output data
     * @tparam U Type of cached data
     * @param cache_data Cache parameters
     * @param cached_data Source data in cache
     * @param data_ptr Destination pointer
     * @param copy Whether to copy data or just set pointer
     */
    template <typename T, typename U, std::enable_if_t<!std::is_pointer<T>::value>* = nullptr>
    void copy_data(Cache_data& cache_data, U* cached_data, T* data_ptr, bool copy) {
        copy_data(cache_data, cached_data, data_ptr);  // Dereference the pointer
    }

    /**
     * @brief Reads a 2D variable from netCDF file into cache and copies to output
     * If cache is not used, the data is always read from the file when the function is called
     * @tparam Var Type of netCDF variable
     * @tparam T Type of output data
     * @tparam U Type of cached data
     * @param data_ptr Destination data pointer
     * @param nc_var netCDF variable to read
     * @param data_cache Cache vector for data
     * @param start Starting indices for reading
     * @param count Number of elements to read in each dimension
     * @param start_cache Start bounds of the cache
     * @param end_cache End bounds of the cache
     * @param copy Whether to copy data or just set pointer
     * @param prodId Product ID being read
     * @param mask Optional mask for selective reading
     */
    template <typename Var, typename T, typename U>
    void read2DVariable(T* data_ptr, Var& nc_var, std::vector<U>& data_cache, std::vector<size_t>& start,
                        std::vector<size_t>& count, cache_bounds& start_cache, cache_bounds& end_cache,
                        bool copy, int prodId, const std::vector<uint8_t>& mask = {}) {
        // if cache is not used and mask is empty just read the data from the file
        if (!use_cache && mask.empty()) {
            if constexpr (std::is_pointer_v<T>)  // passed by pointer, don't copy the data, the data is held
                                                 // in the cache
            {
                resize_nc_buffer(data_cache, count);
                nc_var.getVar(start, count, data_cache.data());
                *data_ptr = data_cache.data();
            } else  // the user allocated memory on their side, read it directly
            {
                nc_var.getVar(start, count, data_ptr);
            }
            // diagnostic data
            total_reads++;
            return;
        }
        // if the cache is used, initialize it
        Cache_data cache_data = set_cached_data(prodId, start, count, start_cache, end_cache, mask);
        // set the cache
        bool read_l2_file = cache_data.read_l2_file;
        std::vector<size_t> start_nc = cache_data.start_nc;
        std::vector<size_t> count_nc = cache_data.count_nc;

        if (read_l2_file) {
            resize_nc_buffer(data_cache, count_nc);
            nc_var.getVar(start_nc, count_nc, data_cache.data());
            // diagnostic data
            total_reads++;
        }
        // get the pointer to cached data
        U* cached_data = data_cache.data();
        copy_data(cache_data, cached_data, data_ptr, copy);
    }

   public:
    L2_Reader(const std::string& file_path)
        : file_path(file_path) {
              NC_CHECK(file_nc.open(file_path, netCDF::NcFile::read))
          };
    L2_Reader(const L2_Reader& other) {
        this->file_path = other.file_path;
        this->file_nc.open(this->file_path, netCDF::NcFile::read);
    }
    L2_Reader operator=(const L2_Reader& other) {
        this->file_path = other.file_path;
        this->file_nc.open(this->file_path, netCDF::NcFile::read);
        return *this;
    }
    L2_Reader() = default;
    /**
     * @brief Initialize L2 quality flags from the input file
     * @details Reads the l2_flags variable and its associated flag meanings and masks
     * @return 0 on success, 1 if l2_flags variable not found
     */
    int ini_l2_flags();

    /**
     * @brief Initialize quality flags for a specific product
     * @param qual_product Name of the quality flag product to initialize
     * @return 0 on success, 1 if quality flag variable not found
     */
    int ini_quality_flags(const std::string& qual_product);

    /**
     * @brief Get the geospatial bounds information
     * @return Reference to the Geospatialbounds object containing spatial metadata
     */
    Geospatialbounds& get_geospatial();
    ~L2_Reader();
    /**
     * @brief Initialize geolocation data by finding latitude/longitude variables and checking dimensions
     * If latitude and/or longitude don't exist or have different dimensions it throws an error
     */
    void iniGeolocation();

    /**
     * @brief Get the platform name
     * @return Platform name string
     */
    std::string get_platform() const;

    /**
     * @brief Get the instrument name
     * @return Instrument name string
     */
    std::string get_instrument() const;

    /**
     * @brief Get begin and end scan rows based on latitude difference
     * @param bscan Vector to store begin scan rows
     * @param escan Vector to store end scan rows
     * @param dlat Latitude difference threshold
     */
    void get_escan_bscan_row(std::vector<int>& bscan, std::vector<int>& escan, float dlat);

    /**
     * @brief Set up variables for reading L2 data
     * @param product_list Comma-separated list of products to read
     * @param wave_list List of wavelengths
     * @param products_requested_l3b Vector to store expanded L3B product names
     * @return 0 on success, 1 if no L2 product found
     */
    int setVariables(std::string& product_list, std::string& wave_list,
                     std::vector<std::string>& products_requested_l3b);

    /**
     * @brief Get the dimensions of the L2 data
     * @param lines Output parameter for number of scan lines
     * @param pixels Output parameter for number of pixels per line
     */
    void getDimensions(size_t& lines, size_t& pixels) const;

    /**
     * @brief Read 2D L2 data for a specific product
     * @tparam T Data type of output array
     * @param data Pointer to output data array
     * @param start Starting indices for reading
     * @param count Number of elements to read in each dimension
     * @param prodId Product index
     * @param mask Optional mask array
     * @param copy Whether to copy data (default true)
     */
    template <typename T>
    void readL2data(T* data, std::vector<size_t>& start, std::vector<size_t>& count, size_t prodId,
                    const std::vector<uint8_t>& mask = {}, bool copy = true) {
        read2DVariable(data, l2_products[prodId], data_products_cached[prodId], start, count,
                       start_cached_products[prodId], end_cached_products[prodId], copy, prodId, mask);
    }

    /**
     * @brief Read 2D L2 data for all products
     * @tparam T Data type of output arrays
     * @param data Vector of output data arrays
     * @param start Starting indices for reading
     * @param count Number of elements to read in each dimension
     * @param mask Optional mask array
     * @param copy Whether to copy data
     */
    template <typename T>
    void readL2data(std::vector<T>& data, std::vector<size_t>& start, std::vector<size_t>& count,
                    const std::vector<uint8_t>& mask, bool copy) {
        for (size_t prodId = 0; prodId < l2_products.size(); prodId++) {
            readL2data(&data[prodId], start, count, prodId, mask, copy);
        }
    }

    /**
     * @brief Read latitude data
     * @tparam T Data type of output array
     * @param data Pointer to output data array
     * @param start Starting indices for reading
     * @param count Number of elements to read in each dimension
     * @param mask Optional mask array
     * @param copy Whether to copy data (default true)
     */
    template <typename T>
    void readLatitude(T* data, std::vector<size_t>& start, std::vector<size_t>& count,
                      const std::vector<uint8_t>& mask = {}, bool copy = true) {
        read2DVariable(data, lat_nc, latitude_cached, start, count, start_cached_latitude,
                       end_cached_latitude, copy, -1, mask);
    }

    /**
     * @brief Read longitude data
     * @tparam T Data type of output array
     * @param data Pointer to output data array
     * @param start Starting indices for reading
     * @param count Number of elements to read in each dimension
     * @param mask Optional mask array
     * @param copy Whether to copy data (default true)
     */
    template <typename T>
    void readLongitude(T* data, std::vector<size_t>& start, std::vector<size_t>& count,
                       const std::vector<uint8_t>& mask = {}, bool copy = true) {
        read2DVariable(data, lon_nc, longitude_cached, start, count, start_cached_longitude,
                       end_cached_longitude, copy, -1, mask);
    }

    /**
     * @brief Read L2 flags data
     * @tparam T Data type of output array
     * @param data Pointer to output data array
     * @param start Starting indices for reading
     * @param count Number of elements to read in each dimension
     * @param mask Optional mask array
     * @param copy Whether to copy data (default true)
     */
    template <typename T>
    void readL2Flags(T* data, std::vector<size_t>& start, std::vector<size_t>& count,
                     const std::vector<uint8_t>& mask = {}, bool copy = true) {
        read2DVariable(data, l2_flags_nc, l2flags_cached, start, count, start_cached_l2flags,
                       end_cached_l2flags, copy, -1, mask);
    }

    /**
     * @brief Read quality flags data
     * @tparam T Data type of output array
     * @param data Pointer to output data array
     * @param start Starting indices for reading
     * @param count Number of elements to read in each dimension
     * @param mask Optional mask array
     * @param copy Whether to copy data (default true)
     */
    template <typename T>
    void readQualityFlags(T* data, std::vector<size_t>& start, std::vector<size_t>& count,
                          const std::vector<uint8_t>& mask = {}, bool copy = true) {
        read2DVariable(data, quality_flags_nc, qualityflags_cached, start, count, start_cached_quality,
                       end_cached_quality, copy, -1, mask);
    }
    /**
     * @brief Reads l2 products data (without copy) for a given scan
     * @param data vector of float pointers
     * @param scan Scan
     * @param mask (optional) Reallocate cache if needed. Supplies begin and end scan
     */
    void readL2dataScan(std::vector<float*>& data, size_t scan, const std::vector<uint8_t>& mask = {});

    /**
     * @brief Reads latitude (without copy) for a given scan
     *
     * @param data Pointer to latitude
     * @param scan Scan
     * @param mask (optional) Reallocate cache if needed. Supplies begin and end scan
     */
    void readLatitudeScan(float** data, size_t scan, const std::vector<uint8_t>& mask = {});

    /**
     * @brief Reads longitude (without copy) for a given scan
     *
     * @param data Pointer to longitude
     * @param scan Scan
     * @param mask (optional) Reallocate cache if needed. Supplies begin and end scan
     */
    void readLongitudeScan(float** data, size_t scan, const std::vector<uint8_t>& mask = {});

    /**
     * @brief Reads l2_flags (without copy) for a given scan
     *
     * @param data Pointer to l2_flags
     * @param scan Scan
     * @param mask (optional) Reallocate cache if needed. Supplies begin and end scan
     */
    void readL2FlagsScan(int** data, size_t scan, const std::vector<uint8_t>& mask = {});

    /**
     * @brief Reads quality_flags (without copy) for a given scan
     *
     * @param data Pointer to quality_flags
     * @param scan Scan
     * @param mask (optional) Reallocate cache if needed. Supplies begin and end scan
     */
    void readQualityFlagsScan(int32_t** data, size_t scan, const std::vector<uint8_t>& mask = {});

    /**
     * @brief get granule start time
     * @return
     */
    double get_start_time() const;

    /**
     * @brief get granule end time
     * @return
     */
    double get_end_time() const;

    /**
     * @brief Get the minimum valid values for all products
     * @return Vector of minimum values, one for each product
     */
    std::vector<float> get_min_value_product() const;

    /**
     * @brief Get the maximum valid values for all products
     * @return Vector of maximum values, one for each product
     */
    std::vector<float> get_max_value_product() const;

    /**
     * @brief Get the units for all products
     * @return Vector of unit strings, one for each product
     */
    std::vector<std::string> get_units() const;

    /**
     * @brief Resets the internal data cache
     *
     * Clears any cached data and resets the cache state to empty
     */
    void reset_cache();

    /**
     * @brief Gets the full path of the L2 file being read
     *
     * @return const std::string& The file path as a constant reference
     */
    const std::string& get_filename() const {
        return file_path;
    }

    /**
     * @brief Gets the mapping between L2 flag meanings and their bit positions
     *
     * @return std::unordered_map<std::string, int> Map of flag names to bit positions
     */
    std::unordered_map<std::string, int> get_l2_meaning_bit_dict() const;

    /**
     * @brief Gets the total number of I/O read operations performed
     *
     * Used for diagnostic purposes to track I/O performance
     *
     * @return size_t The total count of read operations
     */
    size_t get_total_io_reads() const {
        return total_reads;
    }

    /**
     * @brief Reopens the L2 file
     */
    void reopenL2();

    /**
     * @brief Finds the index of a product in the internal product list
     * @param product_name Name of the product to find
     * @param index Output parameter to store the found index
     * @return 0 if product found, non-zero if not found
     */
    int32_t find_product_index(const std::string& product_name, size_t& index);

    /**
     * @brief Gets the attributes associated with a product
     * @param product_name Name of the product
     * @return Map of attribute names to values for the product
     */
    std::map<std::string, float> get_product_attributes(const std::string& product_name) const;

    /**
     * @brief Gets the L2 metadata object
     * @return Reference to the L2 metadata
     */
    MetaL2& get_l2_metadata() {
        return *l2_metadata;
    }

    /**
     * @brief Sets up the L2 metadata object from the netCDF file
     *     
     * */
    void set_l2_metadata() {
        l2_metadata = std::make_unique<MetaL2>(file_nc);
    }

    /**
     * @brief return netCDF variable
     * @param name variable name to find
     * @return netCDF var
     */
    netCDF::NcVar get_variable(const std::string& name);

    /**
     * @brief set cache usage
     * @param cache_usage if true, cache is used
     */
    void set_cache_flag(bool cache_usage);
};

#endif
