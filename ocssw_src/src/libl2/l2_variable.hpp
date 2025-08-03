#ifndef L2_VARIABLE_HPP
#define L2_VARIABLE_HPP

#include <netcdf>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <productInfo.h>
class L2Variable {
   private:
   // input netcdf file name
    std::string input_file_name;
    // Maps 3D product names to corresponding L2 product names
    std::unordered_map<std::string, std::string> l3d_to_l2_prod_map;
    // Stores wavelength attributes for 3D products
    std::unordered_map<std::string, float> l3d_wave_attr;
    // Stores start/end indexes for 3D product dimensions
    std::unordered_map<std::string, std::pair<int, int>> l3d_indexes;
    // Maps product names to their wavelength values (as integers)
    std::unordered_map<std::string, std::vector<int>> three_dims_prod_wavelength;
    // Maps product names to their original wavelength values (as floats)
    std::unordered_map<std::string, std::vector<float>> three_d_original_prod_wavelength;
    // Maps wavelength dimension to products
    std::unordered_map<std::string, std::set<std::string>> prod_wavenames;
    // Set of original L2 product names before any processing
    std::unordered_set<std::string> original_l2_products;
    // Maps product names to their min/max value constraints
    std::unordered_map<std::string, std::pair<float, float>> min_max_values;
    // Maps product names to valid min/max ranges from netCDF file
    std::unordered_map<std::string, std::pair<float, float>> validMin_validMax_nc_file;
    // Maps product names to their units of measurement
    std::unordered_map<std::string, std::string> products_units;
    // Set of all wavelength values used across products
    std::unordered_set<int> wavelength_all;
    // Platform and instrument identifiers
    std::string platform, instrument;
    // Numeric sensor identifier
    int sensorID;
    // Pointer to product information structure
    productInfo_t *info = nullptr;

    /**
     * @brief Parses wavelength subsetting options from a string format
     * 
     * Accepts input in format: "Wavelength_Used_all=354:2200;Wavelength_Used_DTDB=870,1640"
     * Each subset is separated by semicolons and contains a wavelength name and value separated by equals
     * If no equals sign is present, applies the subset to all wavelengths
     * 
     * @param subsets_wave_requested String containing the wavelength subset specifications
     * @return Map of product names to their wavelength subset strings
     */
    std::unordered_map<std::string, std::string> parse_wave_subsets(
        const std::string &subsets_wave_requested);

   public:
    L2Variable() {
        info = allocateProductInfo();
    };
    ~L2Variable() {
        if (info) {
            freeProductInfo(info);
        }
    }
    L2Variable(netCDF::NcFile &file);
    L2Variable(netCDF::NcFile &l2_file, std::string &products_requested,
               std::vector<std::string> &products_requested_separated,
               const std::string &requested_wavelengths, const std::string &file_name);
    void read(netCDF::NcFile &l2_file);
    /**
     * @brief
     *  readL2, openL2 uses the maps from this module. If the module is not initialized, the  function read
     * and process an unexpanded L2 products
     * @param file_name  input L2 file name
     * @param products_requested products requested by the user, separated by comma and/or space, input,
     * string
     * @param products_requested_separated list of the separated products, output, vector
     * @param requested_wavelengths list of requested wavelength, input, string
     */
    void read(netCDF::NcFile &l2_file, std::string &products_requested,
              std::vector<std::string> &products_requested_separated,
              const std::string &requested_wavelengths, const std::string &file_name);
    /**
     * @brief
     * Sets arrays for min/max values and returns only product names without min/max specifier
     * @param products_requested_separated - requested products with min and max values
     * @return std::vector<std::string> - expanded products without min and max values
     */
    std::vector<std::string> create_min_max_values(
        const std::vector<std::string> &products_requested_separated);

    /**
     * @brief
     * The user can supply a 3D product such as  "Lt" or "Rrs" and they will be expanded: Rrs_{WAVE} where WAVE is a wavelength
     * @param products_requested_separated - list of all products.
     */
    void expand_l2_requested(std::vector<std::string> &products_requested_separated);
    /**
     * @brief Get the 3D product name using a wavelength string
     * @param prod_3d The 3D product identifier
     * @param wv The wavelength as a string
     * @return The expanded 2D product name
     */
    std::string get_3d_product_name(const std::string &prod_3d, const std::string &wv);

    /**
     * @brief Get the 3D product name using an integer wavelength
     * @param prod_3d The 3D product identifier
     * @param wv The wavelength as an integer
     * @return The expanded 2D product name
     */
    std::string get_3d_product_name(const std::string &prod_3d, int wv);

    /**
     * @brief Get the 3D product name and extract wavelength into pointer
     * @param prod_3d The 3D product identifier
     * @param wv Pointer to store the extracted wavelength
     * @return The full 3D product name
     */
    std::string get_3d_product_name(const std::string &prod_3d, int *wv);

    /**
     * @brief Convert L3 product name to L2 product name
     * @param l3_name The L3 product name
     * @return The corresponding L2 product name
     */
    std::string l3_l2_conversion(const std::string &l3_name);

    /**
     * @brief Get the map of product units
     * @return Map of product names to their units
     */
    std::unordered_map<std::string, std::string> get_units_list();

    /**
     * @brief Get the map of product min/max values
     * @return Map of product names to their min/max value pairs
     */
    std::unordered_map<std::string, std::pair<float, float>> get_min_max_values();

    /**
     * @brief Get the slice index for a 3D product
     * @param l3_name The L3 product name
     * @return The slice index
     */
    size_t slice_3d_index(std::string l3_name);

    /**
     * @brief Get the attributes for an L3 product
     * @param prod_name The product name
     * @return Map of attribute names to values
     */
    std::map<std::string, float> l3_attrs(const std::string &prod_name);
};

#endif
