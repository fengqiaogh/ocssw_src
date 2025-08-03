/**
 *
 *
 * @author Jakob C. Lindo (SSAI)
 * @date March 2024
 * @brief An extension of NetCDF::NcVar to perform scale and offset automatically
 *
 */

#include "scaled_nc_var.hpp"
#include <ncVarAtt.h>
#include <math.h>
#include <cmath>  // For std::isnan()
#include <limits>

using namespace netCDF;
using namespace netCDF::exceptions;
using namespace std;

ScaledNcVar::ScaledNcVar() : netCDF::NcVar() {
}

ScaledNcVar::ScaledNcVar(const NcVar &copied) : netCDF::NcVar(copied) {
    if (!isNull()) {
        thisVarType = getType().getTypeClass();
    }
}

ScaledNcVar::~ScaledNcVar() {
    if (prodInfo) {
        freeProductInfo(prodInfo);
    }
}

bool ScaledNcVar::isFloatingPoint() {
    return thisVarType == NcType::nc_DOUBLE || thisVarType == NcType::nc_FLOAT;
}

double ScaledNcVar::initFillValue() {
    switch (thisVarType) {  // nc_type
        case (NC_BYTE):
            return PRODUCT_DEFAULT_fillValue_byte;
        case (NC_UBYTE):
            return PRODUCT_DEFAULT_fillValue_ubyte;
        case (NC_USHORT):
            return PRODUCT_DEFAULT_fillValue_ushort;
        case (NC_UINT):
            return PRODUCT_DEFAULT_fillValue_uint;
        default:
            return PRODUCT_DEFAULT_fillValue;
    }
}

void ScaledNcVar::assignBadValue(double newValue) {
    this->badValue = newValue;
}

void ScaledNcVar::assignFillValue(double newValue) {
    switch (thisVarType) {  // nc_type
        case (NC_BYTE): {
            int8_t fill;
            if (newValue >= numeric_limits<int8_t>::max() || newValue <= numeric_limits<int8_t>::lowest()) {
                throw out_of_range("Given fill value is outside the range of this variable");
            } else {
                fill = static_cast<int8_t>(newValue);
            }
            setFill(true, fill);
            break;
        }

        case (NC_UBYTE): {
            uint8_t fill;
            if (newValue >= numeric_limits<uint8_t>::max() || newValue <= numeric_limits<uint8_t>::lowest()) {
                throw out_of_range("Given fill value is outside the range of this variable");
            } else {
                fill = static_cast<uint8_t>(newValue);
            }
            setFill(true, fill);
            break;
        }

        case (NC_SHORT): {
            short fill;
            if (newValue >= numeric_limits<short>::max() || newValue <= numeric_limits<short>::lowest()) {
                throw out_of_range("Given fill value is outside the range of this variable");
            } else {
                fill = static_cast<short>(newValue);
            }
            setFill(true, fill);
            break;
        }

        case (NC_USHORT): {
            uint16_t fill;
            if (newValue >= numeric_limits<uint16_t>::max() ||
                newValue <= numeric_limits<uint16_t>::lowest()) {
                throw out_of_range("Given fill value is outside the range of this variable");
            } else {
                fill = static_cast<uint16_t>(newValue);
            }
            setFill(true, fill);
            break;
        }

        case (NC_INT): {
            int fill;
            if (newValue >= numeric_limits<int>::max() || newValue <= numeric_limits<int>::lowest()) {
                throw out_of_range("Given fill value is outside the range of this variable");
            } else {
                fill = static_cast<int>(newValue);
            }
            setFill(true, fill);
            break;
        }

        case (NC_UINT): {
            uint32_t fill;
            if (newValue >= numeric_limits<uint32_t>::max() ||
                newValue <= numeric_limits<uint32_t>::lowest()) {
                throw out_of_range("Given fill value is outside the range of this variable");
            } else {
                fill = static_cast<uint32_t>(newValue);
            }
            setFill(true, fill);
            break;
        }

        case (NC_FLOAT): {
            float fill;
            if (newValue >= numeric_limits<float>::max() || newValue <= numeric_limits<float>::lowest()) {
                throw out_of_range("Given fill value is outside the range of this variable");
            } else {
                fill = static_cast<float>(newValue);
            }
            setFill(true, fill);
            break;
        }

        default:
            setFill(true, newValue);
    }
    this->fillValue = newValue;
}

size_t ScaledNcVar::getDimsSize() {
    const vector<NcDim> dims = this->getDims();
    size_t size = 1;
    for (const auto &dim : dims)
        size *= dim.getSize();
    return size;
}

void ScaledNcVar::getScaleFactorsFromFile() {
    try {
        getAtt("scale_factor").getValues(&scaleFactor);
    } catch (NcException const &e) {
        scaleFactor = 1.0;
    }
    try {
        getAtt("add_offset").getValues(&addOffset);
    } catch (NcException const &e) {
        addOffset = 0.0;
    }
    try {
        getAtt("_FillValue").getValues(&fillValue);
    } catch (NcException const &e) {
        fillValue = BAD_FLT;
    }
    scaleFactorsSet = true;
}

void ScaledNcVar::getVar(float *dataValues) {
    NcVar::getVar(dataValues);

    if (!scaleFactorsSet)
        getScaleFactorsFromFile();

    if (!isFloatingPoint())
        ScaledNcVar::uncompress(dataValues, getDimsSize());
    else
        fillToBad(dataValues, getDimsSize());  // since uncompress already does this
}

void ScaledNcVar::getVar(double *dataValues) {
    NcVar::getVar(dataValues);

    if (!scaleFactorsSet)
        getScaleFactorsFromFile();

    if (!isFloatingPoint())
        ScaledNcVar::uncompress(dataValues, getDimsSize());
    else
        fillToBad(dataValues, getDimsSize());  // since uncompress already does this
}

void ScaledNcVar::getVar(vector<size_t> start, vector<size_t> count, float *dataValues) {
    size_t sizeToGet = 1;
    for (size_t stop : count)
        sizeToGet *= stop;

    NcVar::getVar(start, count, dataValues);

    if (!scaleFactorsSet)
        getScaleFactorsFromFile();

    if (!isFloatingPoint())
        ScaledNcVar::uncompress(dataValues, sizeToGet);
    else
        fillToBad(dataValues, sizeToGet);  // since uncompress already does this
}

void ScaledNcVar::getVar(vector<size_t> start, vector<size_t> count, double *dataValues) {
    size_t sizeToGet = 1;
    for (size_t stop : count)
        sizeToGet *= stop;

    NcVar::getVar(start, count, dataValues);

    if (!scaleFactorsSet)
        getScaleFactorsFromFile();

    if (!isFloatingPoint())
        ScaledNcVar::uncompress(dataValues, sizeToGet);
    else
        fillToBad(dataValues, sizeToGet);  // since uncompress already does this
}

void ScaledNcVar::putVar(const float *dataValues) {
    size_t sizeToWrite = getDimsSize();

    if (isFloatingPoint() || !scaleFactorsSet) {
        if (fillValue != badValue) {
            vector<double> buf(sizeToWrite);
            badToFill(dataValues, sizeToWrite, buf);
            NcVar::putVar(buf.data());
        } else {
            NcVar::putVar(dataValues);
        }
        return;
    }

    // scale factors are set and this is not a float
    vector<int32_t> buf(dataValues, dataValues + sizeToWrite);
    compress(dataValues, buf, sizeToWrite);
    NcVar::putVar(buf.data());
}

void ScaledNcVar::putVar(const double *dataValues) {
    size_t sizeToWrite = getDimsSize();

    if (isFloatingPoint() || !scaleFactorsSet) {
        if (fillValue != badValue) {
            vector<double> buf(sizeToWrite);
            badToFill(dataValues, sizeToWrite, buf);
            NcVar::putVar(buf.data());
        } else {
            NcVar::putVar(dataValues);
        }
        return;
    }

    // scale factors are set and this is not a float
    vector<int32_t> buf(dataValues, dataValues + sizeToWrite);
    compress(dataValues, buf, sizeToWrite);
    NcVar::putVar(buf.data());
}

void ScaledNcVar::putVar(vector<size_t> start, vector<size_t> count, const float *dataValues) {
    size_t sizeToWrite = 1;
    for (size_t stop : count) {
        sizeToWrite *= stop;
    }

    if (isFloatingPoint() || !scaleFactorsSet) {
        if (fillValue != badValue) {
            vector<double> buf(sizeToWrite);
            badToFill(dataValues, sizeToWrite, buf);
            NcVar::putVar(start, count, buf.data());
        } else {
            NcVar::putVar(start, count, dataValues);
        }
        return;
    }

    // scale factors are set and this is not a float
    vector<int32_t> buf(dataValues, dataValues + sizeToWrite);
    compress(dataValues, buf, sizeToWrite);
    NcVar::putVar(start, count, buf.data());
}

void ScaledNcVar::putVar(vector<size_t> start, vector<size_t> count, const double *dataValues) {
    size_t sizeToWrite = 1;
    for (size_t stop : count) {
        sizeToWrite *= stop;
    }

    if (isFloatingPoint() || !scaleFactorsSet) {
        if (fillValue != badValue) {
            vector<double> buf(sizeToWrite);
            badToFill(dataValues, sizeToWrite, buf);
            NcVar::putVar(start, count, buf.data());
        } else {
            NcVar::putVar(start, count, dataValues);
        }
        return;
    }

    // scale factors are set and this is not a float
    vector<int32_t> buf(dataValues, dataValues + sizeToWrite);
    compress(dataValues, buf, sizeToWrite);
    NcVar::putVar(start, count, buf.data());
}

void ScaledNcVar::getFillValue(double *fillValue) {
    *fillValue = this->fillValue;
}
void ScaledNcVar::getFillValue(float *fillValue) {
    *fillValue = this->fillValue;
}

void ScaledNcVar::compress(const double *toCompress, std::vector<int32_t> &compressed, size_t count) {
    compressed.resize(count);
    pair<double, double> thisVarRange = range();  // [lowerBound, upperBound]

    for (size_t i = 0; i < count; i++) {
        if (toCompress[i] == badValue || isnan(toCompress[i])) {
            compressed[i] = fillValue;
        } else {
            double value = (toCompress[i] - this->addOffset) / this->scaleFactor;

            if (value < thisVarRange.first || thisVarRange.second < value) {
                compressed[i] = fillValue;
            } else {
                compressed[i] = value;
            }
        }
    }
}

void ScaledNcVar::compress(const float *toCompress, std::vector<int32_t> &compressed, size_t count) {
    compressed.resize(count);
    pair<double, double> thisVarRange = range();  // [lowerBound, upperBound]

    for (size_t i = 0; i < count; i++) {
        if (toCompress[i] == badValue || isnan(toCompress[i])) {
            compressed[i] = fillValue;
        } else {
            double value = (toCompress[i] - this->addOffset) / this->scaleFactor;

            if (value < thisVarRange.first || thisVarRange.second < value) {
                compressed[i] = fillValue;
            } else {
                compressed[i] = value;
            }
        }
    }
}

void ScaledNcVar::uncompress(double *toUncompress, size_t count) {
    for (size_t i = 0; i < count; i++) {
        if (toUncompress[i] == fillValue || isnan(toUncompress[i]))
            toUncompress[i] = badValue;
        else if (scaleFactor != 1.0 || addOffset != 0.0)
            toUncompress[i] = (toUncompress[i] * scaleFactor) + addOffset;
    }
}

void ScaledNcVar::uncompress(float *toUncompress, size_t count) {
    for (size_t i = 0; i < count; i++) {
        if (toUncompress[i] == fillValue || isnan(toUncompress[i]))
            toUncompress[i] = badValue;
        else if (scaleFactor != 1.0 || addOffset != 0.0)
            toUncompress[i] = (toUncompress[i] * scaleFactor) + addOffset;
    }
}

template <typename T>
void ScaledNcVar::fillToBad(T *data, size_t count) {
    if (fillValue != badValue) {
        for (size_t i = 0; i < count; i++) {
            if (data[i] == fillValue || isnan(data[i]))
                data[i] = badValue;
        }
    }
}

template <typename T, typename E>
void ScaledNcVar::badToFill(const T *data, const size_t &count, vector<E> &out) {
    out.resize(count);
    for (size_t i = 0; i < count; i++) {
        if (data[i] == badValue || isnan(data[i]))
            out[i] = fillValue;
        else
            out[i] = data[i];
    }
}

bool ScaledNcVar::populateScaleFactors(int sensorID) {
    if (!prodInfo) {
        prodInfo = allocateProductInfo();

        if (findProductInfo(this->NcVar::getName().c_str(), sensorID, prodInfo) != 1) {
            freeProductInfo(prodInfo);
            prodInfo = nullptr;
            return false;
        }
    }

    if (prodInfo->description && strcmp(prodInfo->description, PRODUCT_DEFAULT_description))
        putAtt("long_name", prodInfo->description);

    if (prodInfo->units && strcmp(prodInfo->units, PRODUCT_DEFAULT_units) &&
        strcmp(prodInfo->units, "dimensionless"))
        putAtt("units", prodInfo->units);

    if (prodInfo->standardName)
        putAtt("standard_name", prodInfo->standardName);

    assignFillValue(prodInfo->fillValue);

    if (prodInfo->validMin != PRODUCT_DEFAULT_validMin)
        putAtt("valid_min", thisVarType, prodInfo->validMin);

    if (prodInfo->validMax != PRODUCT_DEFAULT_validMax)
        putAtt("valid_max", thisVarType, prodInfo->validMax);

    // If either scaleFactor or addOffset are set, they both will be.
    if ((prodInfo->scaleFactor != PRODUCT_DEFAULT_scaleFactor) ||
        (prodInfo->addOffset != PRODUCT_DEFAULT_addOffset)) {
        putAtt("scale_factor", NC_DOUBLE, prodInfo->scaleFactor);
        putAtt("add_offset", NC_DOUBLE, prodInfo->addOffset);
        scaleFactor = prodInfo->scaleFactor;
        addOffset = prodInfo->addOffset;
    }

    if (prodInfo->reference != PRODUCT_DEFAULT_reference)
        putAtt("reference", prodInfo->reference);

    if (prodInfo->comment != PRODUCT_DEFAULT_comment)
        putAtt("comment", prodInfo->comment);

    scaleFactorsSet = true;

    return true;
}

void ScaledNcVar::setScaleFactors(double scale, double offset, double fillValue) {
    scaleFactorsSet = true;
    assignFillValue(fillValue);

    // need to set scale factors
    if (scale != 1.0 || offset != 0.0) {
        // not allowed to set scale and offset if not unity for floats
        if (thisVarType == NC_DOUBLE || thisVarType == NC_FLOAT) {
            throw invalid_argument("Setting scale/offset for a float/double is not allowed\n");
        }

        this->scaleFactor = scale;
        this->addOffset = offset;

        putAtt("scale_factor", NC_DOUBLE, scale);
        putAtt("add_offset", NC_DOUBLE, offset);
    }
}

pair<double, double> ScaledNcVar::range() {
    switch (thisVarType) {
        case (NC_FLOAT):
            return pair<double, double>(NC_MIN_FLOAT, NC_MAX_FLOAT);
        case (NC_DOUBLE):
            return pair<double, double>(NC_MIN_DOUBLE, NC_MAX_DOUBLE);
        case (NC_BYTE):
            return pair<double, double>(NC_MIN_BYTE, NC_MAX_BYTE);
        case (NC_UBYTE):
            return pair<double, double>(0, NC_MAX_UBYTE);
        case (NC_SHORT):
            return pair<double, double>(NC_MIN_SHORT, NC_MAX_SHORT);
        case (NC_USHORT):
            return pair<double, double>(0, NC_MAX_USHORT);
        case (NC_INT):
            return pair<double, double>(NC_MIN_INT, NC_MAX_INT);
        case (NC_UINT):
            return pair<double, double>(0, NC_MAX_UINT);
        default:
            return pair<double, double>(1, -1);  // Dumb value
    }
}

netCDF::NcType::ncType productInfoType2ncType(string typeStr) {
    if (typeStr == "byte")
        return NcType::nc_BYTE;
    else if (typeStr == "ubyte")
        return NcType::nc_UBYTE;
    else if (typeStr == "short")
        return NcType::nc_SHORT;
    else if (typeStr == "ushort")
        return NcType::nc_USHORT;
    else if (typeStr == "int")
        return NcType::nc_INT;
    else if (typeStr == "uint")
        return NcType::nc_UINT;
    else if (typeStr == "float")
        return NcType::nc_FLOAT;
    else if (typeStr == "double")
        return NcType::nc_DOUBLE;
    else
        throw runtime_error("productInfoType2ncType could not lookup typeStr = " + typeStr);
    return NcType::nc_FLOAT;
}

ScaledNcVar newScaledNcVar(const netCDF::NcGroup &group, const std::string &name,
                           const std::vector<netCDF::NcDim> &dims, int sensorID) {
    productInfo_t *prodInfo = allocateProductInfo();

    if (!findProductInfo(name.c_str(), sensorID, prodInfo)) {
        freeProductInfo(prodInfo);
        throw runtime_error("newScaledNcVar could not find product = " + name);
    }

    ScaledNcVar var = group.addVar(name, productInfoType2ncType(prodInfo->dataType), dims);
    var.setProdInfo(prodInfo);
    var.populateScaleFactors();
    return var;
}
