/**
 *
 * @author Jakob C. Lindo (SSAI)
 * @date March 2024
 * @brief An extension of NetCDF::NcVar to apply scale & offset when getting a variable as well as other
 * functionality
 *
 */

#ifndef SCALED_NCVAR_H
#define SCALED_NCVAR_H

#define DEFAULT_SENSOR 30  // OCI

#include <netcdf>
#include <productInfo.h>
#include <genutils.h>
#include <type_traits>  // std::enable_if_t and std::is_arithmetic
#include <cmath>        // For std::isnan()

template<typename T>
using enableIfArithmetic = std::enable_if_t<std::is_arithmetic<T>::value>;

/**
 * @brief Adds new functionality to netCDF::NcVar
 *
 * @details Provides new functionality to netCDF:NcVar by applying scale & offset, representing fill values in
 * the file as any custom "bad" value in code, and converting any custom "bad" value to the variable's fill
 * value upon writing
 */
class ScaledNcVar : public netCDF::NcVar {
   public:
    // Default constructor
    ScaledNcVar();

    // NcVar constructor
    ScaledNcVar(const NcVar &copied);

    // destructor
    ~ScaledNcVar();

    /**
     * @brief Dumps all data stored in this ScaledNcVar into the given pointer
     * @tparam T An arithmetic type
     * @param data The storage pointer
     */
    template <typename T, typename = enableIfArithmetic<T>>
    void getVar(T *data) {
        NcVar::getVar(data);

        if (!scaleFactorsSet) {
            getScaleFactorsFromFile();
        }

        if (isFloatingPoint()) { // Floating point data doesn't get compressed
            fillToBad(data, getDimsSize());
        } else {
            uncompress(data, getDimsSize());
        }
    }

    /**
     * @brief Dumps data found along the given dimensions into a given storage pointer
     * @tparam T An arithmetic type
     * @param start A list of starting points for each dimension
     * @param count A list of ending points for each dimension
     * @param data The storage pointer
     */
    template <typename T, typename = enableIfArithmetic<T>>
    void getVar(std::vector<size_t> start, std::vector<size_t> count, T *data) {
        NcVar::getVar(start, count, data);

        if (!scaleFactorsSet) {
            getScaleFactorsFromFile();
        }

        if (isFloatingPoint()) { // Floating point data doesn't get compressed
            fillToBad(data, getDimsSize(count));
        } else {
            uncompress(data, getDimsSize(count));
        }
    }

    /**
     * @brief Dumps all the values stored in this variable to given pointer
     * @param dataValues An array of values that will be mutated
     *
     * @throws `netCDF::exceptions::NcException` if there's an error when reading variable data
     * @throws `std::runtime_error` if there's a general runtime error
     */
    void getVar(float *dataValues);

    /**
     * @brief Dumps all values stored in this variable to the given pointer
     * @param dataValues An array of values that will be mutated
     *
     * @throws `netCDF::exceptions::NcException` if there's an error when reading variable data
     * @throws `std::runtime_error` if there's a general runtime error
     */
    void getVar(double *dataValues);

    /**
     * @brief Dumps a specified contiguous chunk of data to the given pointer
     * @param start The starting indices on each dimension
     * @param count How many values to get on each dimension
     * @param dataValues The values to be stored
     *
     * @throws `netCDF::exceptions::NcException` if there's an error when reading variable data
     * @throws `std::runtime_error` if there's a general runtime error
     */
    void getVar(std::vector<size_t> start, std::vector<size_t> count, float *dataValues);

    /**
     * @brief Dumps a specified contiguous chunk of data to the given pointer
     * @param start The starting indices on each dimension
     * @param count How many values to get on each dimension
     * @param dataValues The values to be stored
     *
     * @throws `netCDF::exceptions::NcException` if there's an error when reading variable data
     * @throws `std::runtime_error` if there's a general runtime error
     */
    void getVar(std::vector<size_t> start, std::vector<size_t> count, double *dataValues);

    /**
     * @brief Writes data to the file
     * @note The size of the memory allocated to the given pointer must equal or exceed the dimensions of this
     * variable
     * @param dataValues The values to be written
     *
     * @throws `netCDF::exceptions::NcException` if there's an error when writing variable data
     * @throws `std::runtime_error` if there's a general runtime error
     */
    void putVar(const float *dataValues);

    /**
     * @brief Writes data to the file
     * @note The size of the memory allocated to the given pointer must equal or exceed the dimensions of this
     * variable
     * @param dataValues The values to be written
     *
     * @throws `netCDF::exceptions::NcException` if there's an error when writing variable data
     * @throws `std::runtime_error` if there's a general runtime error
     */
    void putVar(const double *dataValues);

    template<typename T>
    void putVar(const T *dataValues) {

        size_t sizeToWrite = getDimsSize();

        if (isFloatingPoint() || !scaleFactorsSet) { // Floating point data doesn't get compressed
            if (fillValue != badValue) {
                std::vector<double> buf(sizeToWrite);
                badToFill(dataValues, sizeToWrite, buf);
                NcVar::putVar(buf.data());
            } else {
                NcVar::putVar(dataValues);
            }
            return;
        }

        // scale factors are set and this is not a float
        std::vector<int32_t> buf(dataValues, dataValues + sizeToWrite);
        compress(dataValues, buf, sizeToWrite);
        NcVar::putVar(buf.data());
    }

    /**
     * @brief Write a selected and contiguous range of data to the file at this variable
     * @tparam T An arithmetic type
     * @param start A vector of the indices indicating the beginning of the selected range along each
     * dimension
     * @param count A vector of the indices indicating the ending of the selected range along each dimension
     * @param dataValues The data to be written
     * 
     * @note It is assumed that the memory at the given pointer is valid
     */
    template <typename T, typename = enableIfArithmetic<T>>
    void putVar(std::vector<size_t> start, std::vector<size_t> count, const T *dataValues) {
        size_t sizeToWrite = getDimsSize(count);

        if (isFloatingPoint() || !scaleFactorsSet) {
            if (fillValue != badValue) { // Floating point data doesn't get compressed
                std::vector<double> buf(sizeToWrite);
                badToFill(dataValues, sizeToWrite, buf);
                NcVar::putVar(start, count, buf.data());
            } else {
                NcVar::putVar(start, count, dataValues);
            }
            return;
        }

        // scale factors are set and this is not a float
        std::vector<int32_t> buf(dataValues, dataValues + sizeToWrite);
        compress(dataValues, buf, sizeToWrite);
        NcVar::putVar(start, count, buf.data());
    }

    /**
     * @brief Writes memory from the given pointer to the file along specified dimensions
     *
     * @note The size of the memory allocated to the given pointer must equal or exceed the dimensions of this
     * variable
     *
     * @param dataValues The values to be written
     *
     * @throws `netCDF::exceptions::NcException` if there's an error when writing variable data
     * @throws `std::runtime_error` if there's a general runtime error
     */
    void putVar(std::vector<size_t> start, std::vector<size_t> count, const float *dataValues);

    /**
     * @brief Writes memory from the given pointer to the file along specified dimensions
     *
     * @note The size of the memory allocated to the given pointer must equal or exceed the dimensions of this
     * variable
     *
     * @param dataValues The values to be written
     *
     * @throws `netCDF::exceptions::NcException` if there's an error when writing variable data
     * @throws `std::runtime_error` if there's a general runtime error
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

    productInfo_t *getProductInfo() {
        return prodInfo;
    }

    /**
     * @brief Populate scale factor and add offset from product.xml. Assumes default sensor (30 == OCI)
     *
     * @throws `netCDF::exceptions::NcException` if there's a netCDF error
     *
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
     *
     * @throws `std::out_of_range` When the given fill value is invalid
     * @throws `netCDF::exceptions::NcException` if there's a netCDF error
     * @throws `std::runtime_error` if there's a general runtime error
     */
    void assignFillValue(double newValue);

    void setProdInfo(productInfo_t *prodInfo) {
        this->prodInfo = prodInfo;
    }

    size_t getDimsSize();
    size_t getDimsSize(std::vector<size_t> count);

    // data type for this var.  Default to something we will probably never use in this class
    netCDF::NcType::ncType thisVarType = netCDF::NcType::ncType::nc_COMPOUND;

   private:
    // Indicates whether the scale factors of this NcVar have already been obtained
    bool scaleFactorsSet = false;
    // A scalar value that will multiply each value stored in this NcVar. 1 by default
    double scaleFactor = 1.0;
    // A scalar value that will be added each value stored in this NcVar. 0 by default
    double addOffset = 0.0;
    // fill value the netCDF file uses
    double fillValue = BAD_FLT;
    // sentinal value for bad data in the internal memory (usually the fillValue)
    double badValue = BAD_FLT;
    // pointer to the product structure from the product.xml file
    productInfo_t *prodInfo = nullptr;

    /**
     * @brief Indicates whether this variable is a float or a double
     * @return True when this is an nc_FLOAT or nc_DOUBLE
     */
    bool isFloatingPoint();

    /**
     * @brief For each value of an array, subtract `addOffset` and then divide by `scaleFactor`. Assumes that
     * the vector which will store the compressed version of the data is of the proper size
     * @tparam An arithmetic type
     * @param toCompress The array to be compressed
     * @param compressed The vector into which the values will be stored.
     * @param count The number of elements to compress
     */
    template <typename T, typename = enableIfArithmetic<T>>
    void compress(const T *toCompress, std::vector<int32_t> &compressed, size_t count) {
        compressed.resize(count);
        std::pair<double, double> thisVarRange = range();  // [lowerBound, upperBound]

        for (size_t i = 0; i < count; i++) {
            if (toCompress[i] == badValue || std::isnan(toCompress[i])) {
                compressed[i] = fillValue;
            } else {
                double value = (toCompress[i] - this->addOffset) / this->scaleFactor;

                if (value < thisVarRange.first || thisVarRange.second < value) {
                    compressed[i] = fillValue;
                } else {
                    compressed[i] = static_cast<int32_t>(value);
                }
            }
        }
    }

    /**
     * @brief For each value of an array, multiply by scaleFactor and add addOffset. Assumes that
     * the vector which will store the compressed version of the data is of the proper size. This is an
     * in-place operation.
     * @tparam An arithmetic type
     * @param toUncompress The array to be compressed
     * @param uncompressed The vector into which the values will be stored.
     */
    template <typename T, typename = enableIfArithmetic<T>>
    void uncompress(T *toUncompress, size_t count) {
        for (size_t i = 0; i < count; i++) {
            if (std::isnan(toUncompress[i]) || toUncompress[i] == fillValue)
                toUncompress[i] = badValue;
            else if (scaleFactor != 1.0 || addOffset != 0.0)
                toUncompress[i] = static_cast<T>((toUncompress[i] * scaleFactor) + addOffset);
        }
    }

    /**
     * @brief Modifies `data` in place by replacing every instance of fill value with bad value.
     * Also checks for NaNs, converting those to badValue too
     * @param data The values to be checked for fill value
     * @param count The number of values to check
     */
    template <typename T>
    void fillToBad(T *data, size_t count) {
        if (fillValue != badValue) {
            for (size_t i = 0; i < count; i++) {
                if (data[i] == fillValue || std::isnan(data[i]))
                    data[i] = static_cast<T>(badValue);
            }
        }
    }

    /**
     * @brief Copies values from `data` into `out`, replacing instance of bad value with fill value.
     * Also checks for NaNs, converting those to fill value too.
     *
     * @note `out` will be resized to `count`
     * @param data The values to be checked for bad value
     * @param count The number of values to check
     * @param out The result
     */
    template <typename T, typename E>
    void badToFill(const T *data, const size_t &count, std::vector<E> &out) {
        out.resize(count);
        for (size_t i = 0; i < count; i++) {
            if (data[i] == badValue || std::isnan(data[i]))
                out[i] = static_cast<E>(fillValue);
            else
                out[i] = static_cast<E>(data[i]);
        }
    }

    /**
     * @brief Obtain fill value, scale factor, and offset from underlying NcVar
     *
     * @throws `netCDF::exceptions::NcException`
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
