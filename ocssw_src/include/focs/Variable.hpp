#ifndef FOCS_VARIABLE
#define FOCS_VARIABLE

#include "Product.hpp"
#include "focs/NVector.hpp"
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>
#include "genutils.h"
namespace focs {

//TODO: change attributes to share with product attributes, or at least copy the theory

enum class OcType {
    Char, String, Float, Double, LongDouble,
    Int8, Int16, Int32, Int64,
    Uint8, Uint16, Uint32, Uint64,
    Unknown
};
// Geospatial = the data moves throughout reading tiles
enum class OcGeospatialType {
    Yes, Line, Pixel, No
};

template<typename T>
class VariableScaler {
    public:
        VariableScaler(){}
        virtual ~VariableScaler(){}
        virtual T apply(const T& data) {return data;}
};
template<typename T>
class SlopeInterceptScaler : public VariableScaler<T> {
    public:
        SlopeInterceptScaler(double slope, double intercept) : slope_{slope}, intercept_{intercept} {}
        T apply(const T& data) override {
            return slope_ * data + intercept_;
        }
    private:
        double slope_;
        double intercept_;
};
template<typename T>
class VariableFill {
    public:
        VariableFill(){}
        virtual ~VariableFill(){}
        virtual bool is_fill(const T&) {return false;}
};
template<typename T>
class VariableFillValue : public VariableFill<T> {
    public:
        VariableFillValue(T value) : value_{value} {}
        bool is_fill(const T& data) override {
            return data == value_;
        }
    private:
        T value_;
};
template<typename T>
class VariableFillRange : public VariableFill<T> {
    public:
        VariableFillRange(T min, T max) : min_{min}, max_{max} {}
        bool is_fill(const T& data) override {
            return !(data < min_ || data > max_);
        }
    private:
        T min_;
        T max_;
};
class BaseVariable {
    public:
        // typedef void base_type;
        // typedef void type;

        BaseVariable() = default;
        BaseVariable(const Product& provides) : provides_{provides} {}
        virtual ~BaseVariable(){}

        std::unordered_map<std::string, std::string>::const_iterator attributes_begin() const {return attributes_.begin();}
        std::unordered_map<std::string, std::string>::const_iterator attributes_end() const {return attributes_.end();}

        auto& attributes(){return attributes_;}
        auto& dimensions(){return dimensions_;}

        virtual void* _data(){return nullptr;}
        virtual const void* _data() const { return nullptr; }
        bool has_data() const {return _data() != nullptr;}

        const Product& provides(){ return provides_; }
        Product& get_product(){ return provides_; }
        void provides(const Product& provides) { provides_ = provides; }
        const Product& provides() const { return provides_; }

        const auto& name(){
            return provides().name();
        }

        bool matches(const Product& needs) const { return provides_.matches(needs); }
        bool matches(const BaseVariable& needs) const { return provides_.matches(needs.provides()); }

        virtual OcType get_type() const { return OcType::Unknown; }
        virtual size_t dimension_count() const { return 0; }

        OcGeospatialType is_geospatial() const { return geo_type_; }
        void geospatial(OcGeospatialType geo_type) { geo_type_ = geo_type; }

        virtual void rotate(const size_t number_of_lines){ (void) number_of_lines; }

        template<typename T>
        OcType get_type(T) const {
            if (std::is_same<T, char>::value){
                return OcType::Char;
            } else if (std::is_same<T, std::string>::value){
                return OcType::String;
            } else if (std::is_same<T, float>::value){
                return OcType::Float;
            } else if (std::is_same<T, double>::value){
                return OcType::Double;
            } else if (std::is_same<T, long double>::value){
                return OcType::LongDouble;
            } else if (std::is_same<T, int8_t>::value){
                return OcType::Int8;
            } else if (std::is_same<T, int16_t>::value){
                return OcType::Int16;
            } else if (std::is_same<T, int32_t>::value){
                return OcType::Int32;
            } else if (std::is_same<T, int64_t>::value){
                return OcType::Int64;
            } else if (std::is_same<T, uint8_t>::value){
                return OcType::Uint8;
            } else if (std::is_same<T, uint16_t>::value){
                return OcType::Uint16;
            } else if (std::is_same<T, uint32_t>::value){
                return OcType::Uint32;
            } else if (std::is_same<T, uint64_t>::value){
                return OcType::Uint64;
            }
            return OcType::Unknown;
        }

        netCDF::NcType type_to_nc(){
            switch (get_type()){
            case OcType::Char:
                return netCDF::NcType::nc_CHAR;
            case OcType::String:
                return netCDF::NcType::nc_STRING;
            case OcType::Float:
                return netCDF::NcType::nc_FLOAT;
            case OcType::Double:
                return netCDF::NcType::nc_DOUBLE;
            case OcType::LongDouble:
                return netCDF::NcType::nc_DOUBLE;
            case OcType::Int8:
                return netCDF::NcType::nc_BYTE;
            case OcType::Int16:
                return netCDF::NcType::nc_SHORT;
            case OcType::Int32:
                return netCDF::NcType::nc_INT;
            case OcType::Int64:
                return netCDF::NcType::nc_INT64;
            case OcType::Uint8:
                return netCDF::NcType::nc_UBYTE;
            case OcType::Uint16:
                return netCDF::NcType::nc_USHORT;
            case OcType::Uint32:
                return netCDF::NcType::nc_UINT;
            case OcType::Uint64:
                return netCDF::NcType::nc_UINT64;
            case OcType::Unknown:
                return netCDF::NcType::nc_OPAQUE;
            }
            return {};
        }

    protected:
        Product provides_{};
        std::unordered_map<std::string, std::string> attributes_{};
        std::vector<std::pair<std::string, size_t>> dimensions_{};
        OcGeospatialType geo_type_{OcGeospatialType::Yes}; // only FileReaders will be making data not based on a DataRecord, which makes them geospatial by nature
};

template<typename T, size_t Dims>
class Variable : public BaseVariable {
    public:
        typedef T base_type;
        typedef typename NVector<T, Dims>::type type;

        Variable(){}
        Variable(const Product& provides) : BaseVariable(provides) {}
        Variable(const type&& data) : data_{data} {}
        const void* _data() const override {
            return static_cast<const void*>(&data_);
        }
        void* _data() override {
            return static_cast<void*>(&data_);
        }
        const type& data() const {
            return data_;
        }
        type& data(){
            return data_;
        }
        auto& operator[](const size_t i){
            return data_[i];
        }

        OcType get_type() const override { return BaseVariable::get_type(T{}); }
        size_t dimension_count() const override { return Dims; }

        VariableFill<T>* fill_value() const { return fill_value_.get(); }
        void fill_value(std::unique_ptr<VariableFill<T>>&& fill_value) { fill_value_ = std::move(fill_value); }
        void fill_value(const T& fill_value) { fill_value_ = std::make_unique<VariableFillValue<T>>(fill_value); }

        void scaler(std::unique_ptr<VariableScaler<T>>&& scaler){scaler_ = std::move(scaler);}
        VariableScaler<T>* scaler(){return scaler_.get();}

        bool is_fill_value(const T& value){
            if (fill_value_){
                return fill_value_->is_fill(value);
            }
            return false;
        }
        template<typename TT>
        bool is_fill_value(const TT& value){
            if (fill_value_){
                return fill_value_->is_fill(value);
            }
            return false;
        }
        template<typename TT, typename std::enable_if<std::is_same<TT, double>::value>::type* = nullptr>
        void apply_scaling(TT* v, const size_t count){
            if (scaler_){
                for (size_t i=0;i<count;i++){
                    if (is_fill_value(v[i])){
                        v[i] = BAD_FLT;
                    } else {
                        v[i] = scaler_->apply(v[i]);
                    }
                }
            }
        }
        template<typename TT, typename std::enable_if<std::is_same<TT, long double>::value>::type* = nullptr>
        void apply_scaling(TT* v, const size_t count){
            if (scaler_){
                for (size_t i=0;i<count;i++){
                    if (is_fill_value(v[i])){
                        v[i] = BAD_FLT;
                    } else {
                        v[i] = scaler_->apply(v[i]);
                    }
                }
            }
        }
        template<typename TT, typename std::enable_if<std::is_same<TT, float>::value>::type* = nullptr>
        void apply_scaling(TT* v, const size_t count){
            if (scaler_){
                for (size_t i=0;i<count;i++){
                    if (is_fill_value(v[i])){
                        v[i] = BAD_FLT;
                    } else {
                        v[i] = scaler_->apply(v[i]);
                    }
                }
            }
        }
        template<typename TT, typename std::enable_if<!std::is_floating_point<TT>::value>::type* = nullptr>
        void apply_scaling(TT* v, const size_t count){
            if (scaler_){
                for (size_t i=0;i<count;i++){
                    if (!is_fill_value(v[i])){
                        v[i] = scaler_->apply(v[i]);
                    }
                }
            }
        }
        virtual void rotate(const size_t number_of_lines) override {
            if (Dims == 2){
                switch (is_geospatial()){
                    case OcGeospatialType::Yes:
                    case OcGeospatialType::Line:
                        {
                            using std::swap;
                            const size_t swap_start = data_.size() - number_of_lines;
                            // std::cout << "Rotation, starting at " << swap_start << "\n";
                            for (size_t i=0;i<number_of_lines;i++){
                                swap(data_[i], data_[swap_start + i]);
                            }
                        }
                        break;
                    default:
                        break;
                }
            }
        }
    protected:
        type data_{};
        // FillValueT fill_value_{std::numeric_limits<FillValueT>::max};
        std::unique_ptr<VariableScaler<T>> scaler_{};
        std::unique_ptr<VariableFill<T>> fill_value_{};
};
template<typename T>
std::ostream& operator<<(std::ostream& os, const Variable<T, 2>& v){
    os << "focs::Variable{" << v.provides() << ", ";
    const size_t max_i = v.data().size() - 1;
    const size_t max_j = v.data()[0].size() - 1;
    // os << "(" << max_i << "x" << max_j << "), ";
    os << "[\n";
    for (size_t i=0;i<=max_i;i++){
        os << "\t[";
        for (size_t j=0;j<=max_j;j++){
            os << v.data()[i][j];
            if (j != max_j){
                os << ", ";
            }
        }
        os << "]";
        if (i != max_i){
            os << ",";
        }
        os << "\n";
    }
    os << "]}";
    return os;
}

} // namespace focs

#endif // FOCS_VARIABLE
