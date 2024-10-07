/**
 *
 * @author Jakob C. Lindo (SSAI)
 * @date March 2024
 * @brief An extension of NetCDF::NcVar to perform scale and offset when getting a variable
 *
 */

#ifndef SCALED_NCVAR_H
#define SCALED_NCVAR_H

#define DEFAULT_SENSOR 30  // OCI

#include <netcdf>
#include <productInfo.h>
#include <genutils.h>

class ScaledNcVar : public netCDF::NcVar {
   public:

    // Copy constructor
    ScaledNcVar(const NcVar &copied);

    // destructor
    ~ScaledNcVar();

    /**
     * @brief Takes in a pointer whose memory will be modified. Assumes values found there are integers and
     * performs a scale and offset on them, based on this variable's `scale_factor` and `add_offset`
     * @param dataValues An array of values that will be mutated
     */
    void getVar(float *dataValues);

    /**
     * @brief Takes in a pointer whose memory will be modified. Assumes values found there are integers and
     * performs a scale and offset on them, based on this variable's `scale_factor` and `add_offset`
     * @param dataValues An array of values that will be mutated
     */
    void getVar(double *dataValues);

    /**
     * @brief Get a variable's values between specified indices
     * @param dataValues The values to be stored
     */
    void getVar(std::vector<size_t> start, std::vector<size_t> count, float *dataValues);

    /**
     * @brief Get a variable's values between specified indices
     * @param dataValues The values to be stored
     */
    void getVar(std::vector<size_t> start, std::vector<size_t> count, double *dataValues);

    /**
     * @brief put an entire variable
     * @param dataValues The values to be stored
     */
    void putVar(const float *dataValues);

    /**
     * @brief put an entire variable
     * @param dataValues The values to be stored
     */
    void putVar(const double *dataValues);

    /**
     * @brief Put a variable's values between specified indices
     * @param dataValues The values to be stored
     */
    void putVar(std::vector<size_t> start, std::vector<size_t> count, const float *dataValues);

    /**
     * @brief Put a variable's values between specified indices
     * @param dataValues The values to be stored
     */
    void putVar(std::vector<size_t> start, std::vector<size_t> count, const double *dataValues);

    /**
     * @brief Get the fill value of this NcVar
     * @param fillValue A pointer to a double that will be mutated to contain this NcVar's _FillValue
     */
    void getFillValue(double *fillValue);

    /**
     * @brief Get the fill value of this NcVar
     * @param fillValue A pointer to a float that will be mutated to contain this NcVar's _FillValue
     */
    void getFillValue(float *fillValue);

    productInfo_t* getProductInfo() { return prodInfo; }

    /**
     * @brief Populate scale factor and add offset from product.xml. Assumes default sensor (30 == OCI)
     * @return Whether the product was found in product.xml
     */
    bool populateScaleFactors(int sensorID = DEFAULT_SENSOR);

    /**
     * @brief Populate scale factor and add offset manually. Throws an `invalid_argument`
     */
    void setScaleFactors(double scale, double offset, double fillValue = BAD_FLT);

    /**
     * @brief Reassign the value of badValue
     * @param newValue The desired value of badValue
     */
    void assignBadValue(double newValue);

    /**
     * @brief Reassign the fill value
     * @param newValue The desired fill value
     */
    void assignFillValue(double newValue);

    void setProdInfo(productInfo_t* prodInfo) { this->prodInfo = prodInfo; }

    // netCDF data type for the variable in the file
    netCDF::NcType::ncType thisVarType = getType().getTypeClass();

   private:
    // Indicates whether the scale factors of this NcVar have already been obtained
    bool scaleFactorsSet = false;
    // A scalar value that will be added each value stored in this NcVar. 1 by default
    double scaleFactor = 1.0;
    // A scalar value that will be added each value stored in this NcVar. 0 by default
    double addOffset = 0.0;
    // fill value the netCDF file uses
    double fillValue = BAD_FLT;
    // sentinal value for bad data in teh internal memory (usually the fillValue)
    double badValue = BAD_FLT;
    // pointer to the product structure from the product.xml file
    productInfo_t *prodInfo = nullptr;

    /**
     * @brief Indicates whether this variable is a float or a double
     * @return True when this is an nc_FLOAT or nc_DOUBLE
     */
    bool floatingPoint();

    /**
     * @brief For each value of an array, subtract `addOffset` and then divide by `scaleFactor`. Assumes that
     * the vector which will store the compressed version of the data is of the proper size
     * @param toCompress The array to be compressed
     * @param compressed The vector into which the values will be stored.
     * @param count The number of elements to compress
     */
    void compress(const double *toCompress, std::vector<int32_t> &compressed, size_t count);

    /**
     * @brief For each value of an array, subtract `addOffset` and then divide by `scaleFactor`. Assumes that
     * the vector which will store the compressed version of the data is of the proper size
     * @param toCompress The array to be compressed
     * @param compressed The vector into which the values will be stored.
     * @param count The number of elements to compress
     */
    void compress(const float *toCompress, std::vector<int32_t> &compressed, size_t count);

    /**
     * @brief For each value of an array, subtract `addOffset` and then divide by `scaleFactor`. Assumes that
     * the vector which will store the compressed version of the data is of the proper size. This is an
     * in-place operation.
     * @param toUncompress The array to be compressed
     * @param uncompressed The vector into which the values will be stored.
     */
    void uncompress(double *toUncompress, size_t count);

    /**
     * @brief For each value of an array, subtract `addOffset` and then divide by `scaleFactor`. Assumes that
     * the vector which will store the compressed version of the data is of the proper size. This is an
     * in-place operation.
     * @param toUncompress The array to be compressed
     * @param uncompressed The vector into which the values will be stored.
     */
    void uncompress(float *toUncompress, size_t count);

    /**
     * @brief Runs through the data given and replaces every instance of fill value with bad value
     * @param data The values to be checked for fill value
     * @param count The number of values to check
     */
    template <typename T>
    void fillToBad(T *data, size_t count);

    /**
     * @brief Runs through the data given and replaces every instance of bad value with fill value
     * @param data The values to be checked for fill value
     * @param count The number of values to check
     */
    template <typename T, typename E>
    void badToFill(const T *data, const size_t &count, std::vector<E> &compressed);

    /**
     * @brief Obtain fill value, scale factor, and offset from underlying NcVar
     */
    void getScaleFactorsFromFile();

    /**
     * @brief Initialize fill value based upon the type of this variable
     * @return A decent guess at what the fill value should be
     */
    double initFillValue();

    /**
     * @brief Return the lower bound and upper bound of this variable
     * @return A std::pair containing the lower bound and upper bound of this variable
     */
    std::pair<double, double> range();
};

/**
 * @brief Create a ScaledNcVar using the product.xml for definition
 * @param group ncGroup to put the variable in
 * @param name the name of the variable
 * @param dims dimentions of the varaible
 * @param sensorID optional sensor ID
 * @return a newly created ScaledNcVar class
 */
ScaledNcVar newScaledNcVar(const netCDF::NcGroup &group, const std::string &name,
                            const std::vector<netCDF::NcDim> &dims, int sensorID = DEFAULT_SENSOR);

#endif
