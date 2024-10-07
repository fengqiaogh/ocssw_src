
#ifndef FOCS_PRODUCT
#define FOCS_PRODUCT


#include "StringUtils.hpp"

#include <iostream>
#include <rapidjson/document.h>
#include <rapidjson/istreamwrapper.h>
#include <rapidjson/writer.h>
#include <rapidjson/stringbuffer.h>
#include <rapidjson/ostreamwrapper.h>
#include "rapidjson/error/en.h"
#include <algorithm>
#include <complex>
#include <cstdint>
#include <functional>
#include <set>
#include <string>
#include <typeinfo>
#include <vector>
#include <memory>
#include <type_traits>
#include <typeindex>
#include <typeinfo>
// #include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <unordered_set>
#include <utility>
#include <netcdf>
#define ROUND_ERROR 0.000001
#ifndef ANY_NOEXP
    #include <experimental/any>
    #define ANY_CAST std::experimental::any_cast
    #define ANY std::experimental::any
#else
    #include <any>
    #define ANY_CAST std::any_cast
    #define ANY std::any
#endif
namespace focs {
        template <typename T, size_t N>
        struct LinkedList {
            const T* val;
            LinkedList<T, N - 1> next;
            LinkedList() = default;
            LinkedList(T* arr) {
                val = arr;
                next = LinkedList<T, N - 1>(arr + 1);
            }
            template <typename... Ts>
            LinkedList(const T& _val, Ts&... args) {
                val = &_val;
                next = LinkedList<T, N - 1>(args...);
            }

            T get() const { return *val; }
            const T* data() const { return val; }
        };

        template <typename T>
        struct LinkedList<T, 0> {
            const T* val;
            LinkedList() = default;
            LinkedList(T* arr) { val = arr; }
            template <typename... Ts>
            LinkedList(const T& _val, [[maybe_unused]] Ts&... args) {
                val = &_val;
            }
            T get() const { return *val; }
            const T* data() const { return val; }
        };   
        template<typename T>
        void insert(std::unordered_set<std::type_index> & found_types)
        {
            if(found_types.count(std::type_index(typeid(T)))>0)
            {
                std::cout << "Error : types are already defined " <<  typeid(T).name() << std::endl;
                exit(EXIT_FAILURE);
            }
            else
            {
                found_types.insert(std::type_index(typeid(T)));
            }
        }


        template<typename ... Ts>
        struct VisitorAny
        {
            std::unordered_set<std::type_index> found_types;
            VisitorAny()
            {
                ((insert<Ts>(found_types)),...);
            }
            using first_type = std::tuple_element_t<0, std::tuple<Ts...>>;

            template <typename T>
            constexpr auto convert(T) const {
                return first_type{};
            }

            template <typename Callable, size_t S, typename... Ps>
            void visit(Callable&& callable, LinkedList<ANY, S>& list,
                    Ps... params) const {
                if constexpr (S > 0) {
                    const ANY* any = list.data();
                    (void)(((any->type() == typeid(Ts))
                        ? ((visit(callable, list.next, ANY_CAST<Ts>(any),
                                    params...)),
                            false)
                        : true) &&
                    ...);
                } else {
                    std::forward<Callable>(callable)(*params...);
                }
            }

            template <typename Callable, typename R, size_t S, typename... Ps>
            void visit(Callable&& callable, R& ret, LinkedList<ANY, S>& list,
                    Ps... params) const {
                if constexpr (S > 0) {
                    const ANY* any = list.data();
                    (void)(((any->type() == typeid(Ts))
                        ? ((visit(callable, ret, list.next, ANY_CAST<Ts>(any),
                                    params...)),
                            false)
                        : true) &&
                    ...);
                } else {
                    ret = std::forward<Callable>(callable)(*params...);
                }
            }

            template <typename Callable, typename... Tanys>
            auto apply_visitor(Callable&& callable, Tanys&&... anys) const {
                using R = decltype(std::declval<Callable>()((convert(anys))...));
                constexpr size_t N = sizeof...(Tanys);
                LinkedList<ANY, N> link = LinkedList<ANY, N>(anys...);
                if constexpr (!std::is_void_v<R>) {
                    R ret{};
                    visit(callable, ret, link);
                    return ret;
                } else {
                    visit(callable, link);
                }
            }

        };
        template <typename, typename = void>
        constexpr bool is_iterable{};

        template <typename T>
        constexpr bool is_iterable<T, std::void_t<decltype(std::declval<T>().begin()),
            decltype(std::declval<T>().end())>> =
            true;
        template <typename T>
        std::ostream& operator<<(std::ostream& os, const std::vector<T>& data) noexcept {
            for (auto it = data.begin(); it != data.end(); it++) {
                os << *it << " ";
            }

        return os;
    }

        template<class T, class F>
        inline std::pair<const std::type_index, std::function<void(ANY const&)>>
        to_any_visitor(F const &f)
        {
            return
            {
                std::type_index(typeid(T)),
                [g = f](ANY const &a)
                {
                    if constexpr (std::is_void_v<T>)
                        g();
                    else
                        g(ANY_CAST<T const&>(a));
                }
            };
        }
    // TODO: change attributes to use std::multimap, which will take care of the vector crap
    // (assuming I don't need more than one dimension for an attribute)?
    // For now, the NetCdf reader just makes a comma-separated string out of non-scalar attributes.
    // Stuff downstream can use the StringUtils::stov functions to undo it.
        // std::vector<std::string>,  std::vector<long double>,
        // std::vector<int8_t>, std::vector<int16_t>, std::vector<int32_t>, std::vector<int64_t>,
        // std::vector<uint8_t>, std::vector<uint16_t>, std::vector<uint32_t>, std::vector<uint64_t>>;
    static   VisitorAny<char, std::string, float, double, long double, int8_t, int16_t, int32_t, int64_t, uint8_t, uint16_t, uint32_t, uint64_t,std::vector<int8_t>,
            std::vector<int16_t>, std::vector<int32_t>, std::vector<int64_t>, std::vector<float>,  std::vector<double>,std::vector<uint8_t>, std::vector<uint16_t>, std::vector<uint32_t>, std::vector<uint64_t>> _visitorAny{};
    // class that needs to be changed
    class Attribute {
        public:
            Attribute() = default;
            virtual ~Attribute() = default;
            Attribute(const std::string &name, const ANY & value, size_t attr_len) : name_{name}, a_value_{value}, attr_len_(attr_len) {}
            const std::string& name() const {return name_;}
            const ANY & a_value() const {return a_value_;}
            size_t attr_len() const {return attr_len_;}
            static Attribute parse_attribute(const netCDF::NcAtt& att);
            static Attribute parse_attribute(const rapidjson::Value::ConstMemberIterator  &data);
            double as_double() const {
               return _visitorAny.apply_visitor(AsDouble{},a_value_); //return boost::apply_visitor(AsDouble{}, value_);
            }

            std::string as_string() const {
                return _visitorAny.apply_visitor(AsString{}, a_value_);
            }

            template<typename T>
            std::vector<T> as_vector() const {
                return _visitorAny.apply_visitor(AsVector<T>{}, a_value_);
            }

            friend std::ostream& operator<<(std::ostream& out, const Attribute& attribute){
                out << attribute.name_ << "=";
                _visitorAny.apply_visitor(Print{out}, attribute.a_value_);
                return out;
            }

            friend bool operator<(const Attribute& left, const Attribute& right){
                // if (left.name_ != right.name_){
                if (!boost::iequals(left.name_, right.name_)){
                    return left.name_ < right.name_;
                }
                return _visitorAny.apply_visitor(LessThan{}, left.a_value(), right.a_value());
            }

            friend bool operator==(const Attribute& left, const std::string& name){
                // return (left.name_ == name);
                return boost::iequals(left.name_, name);
            }
            friend bool operator==(const Attribute& left, const Attribute& right){
                // if (left.name_ != right.name_){
                if (!boost::iequals(left.name_, right.name_)){
                    return false;
                }
                return _visitorAny.apply_visitor(Equals{}, left.a_value(), right.a_value());
            }
            static Attribute parse(const std::string& input){
                auto equals = input.find_first_of('=');
                size_t name_start = 0;
                while (input[name_start] == ' '){
                    name_start++;
                }
                size_t name_end = equals;
                while (input[name_end - 1] == ' '){
                    name_end--;
                }
                std::string name(input, name_start, name_end - name_start);

                while (input[++equals] == ' ');

                switch (input[equals]){
                    case '"':
                        return Attribute{name, input.substr(equals+1, input.size() - equals - 2),1};
                    case '\'':
                        return Attribute{name, input[equals+1],1};
                    default:
                        // TODO
                        // This should probably utilize value suffixes, like C (e.g., to force float, 5f)
                        // I also should check for bad parses
                        if (input.substr(equals) == "true"){
                            return Attribute{name, true,1};
                        } else if (input.substr(equals) == "false"){
                            return Attribute{name, false,1};
                        } else if (input.find_first_of('.') == std::string::npos){
                            return Attribute{name, std::stoi(input.substr(equals)),1};
                        } else {
                            return Attribute{name, std::stof(input.substr(equals)),1};
                        }
                }
            }


        protected:
            // Attribute name
            std::string name_;
            // AttributeType value_;
            ANY a_value_;
            // Attribute length
            size_t attr_len_;
        private: 
            struct Print  {
                public:
                    Print(std::ostream& out_) : out{out_} {}

                    void operator()(std::string s) const {
                        out << '"' << s << '"';
                    }
                    void operator()(uint8_t i) const {
                        out << static_cast<int>(i);
                    }
                    void operator()(int8_t i) const {
                        out << static_cast<int>(i);
                    }
                    void operator()(char c) const {
                        out << '\'' << c << '\'';
                    }
                    template <typename T>
                    void operator()(T s) const {
                        out << s;
                    }
                private:
                    std::ostream& out{std::cout};
            };
            struct AsDouble  {
            public:
                AsDouble() {}
                template <typename T, std::enable_if_t<!is_iterable<T> && !std::is_same_v<T, std::string>, bool> = true>
                double operator()(T s) const {
                    return s;
                }
                template <typename T, std::enable_if_t<is_iterable<T> && !std::is_same_v<T, std::string>, bool> = true>
                double operator()(T s) const {
                    return *s.begin();
                }
                double operator()(std::string s) const {
                    return std::stod(s);
                }
            };
            struct AsString  {
            public:
                AsString() {}
                template <typename T, std::enable_if_t<!is_iterable<T> && !std::is_same_v<T, std::string>, bool> = true>
                std::string operator()(T s) const {
                    return std::to_string(s);
                }
                template <typename T, std::enable_if_t<is_iterable<T> && !std::is_same_v<T, std::string>, bool> = true>
                std::string operator()(T s) const {
                    std::string out = "";
                    for (const auto& v : s)
                    {
                        out += std::to_string(v);
                    }
                    return out;
                }
                std::string operator()(std::string s) const {
                    return s;
                }
            };
            template<typename T>
            struct AsVector  {
                public:
                    AsVector() {}
                    template <typename I>
                    std::vector<T> operator()(I s) const {
                        return std::vector<T>(s);
                    }
                    std::vector<T> operator()(std::string s) const {
                        return focs::StringUtils::stov<T>(s);
                    }
            };

            struct Equals {
                template <typename T, typename U, std::enable_if_t<(is_iterable<T> && !is_iterable <U>) && (!std::is_same_v<T, std::string> || !std::is_same_v<U, std::string>), bool> = true>
                bool operator() (T a, U b) const  { return operator()( a.at(0), b); }
                template <typename T, typename U, std::enable_if_t<(!is_iterable<T> && is_iterable <U>) && (!std::is_same_v<T, std::string> || !std::is_same_v<U, std::string>), bool> = true>
                bool operator() (T a, U b) const  { return operator()(a , b.at(0)); }
                template <typename T, typename U, std::enable_if_t<(is_iterable<T> && is_iterable <U>) && (!std::is_same_v<T, std::string> || !std::is_same_v<U, std::string>), bool> = true>
                bool operator() (T a, U b) const  { return operator()(a.at(0) , b.at(0)); }
                template <typename T, std::enable_if_t<!is_iterable<T> && !std::is_same_v<T, std::string>, bool> = true>
                bool operator() (T a, T b) const { return a == b; }
                bool operator()(const std::string& a, const std::string& b) const { return a == b; }

                template <typename T, std::enable_if_t<!is_iterable<T> && !std::is_same_v<T, std::string>, bool> = true>
                bool operator()(std::string, T) const { return false; /* throw std::invalid_argument("Can't compare std::string to non-std::string"); */ }
                template <typename T, std::enable_if_t<!is_iterable<T> && !std::is_same_v<T, std::string>, bool> = true>
                bool operator()(T, std::string) const { return false; /* throw std::invalid_argument("Can't compare std::string to non-std::string"); */ }

                // Integer comparison (extra work to avoid signed/unsigned comparison)
                template <typename T, typename U, typename std::enable_if<std::is_integral<T>::value && std::is_integral<U>::value && std::is_signed<T>::value == std::is_signed<U>::value>::type* = nullptr>
                bool operator() (T a, U b) const { return a == b; }
                template <typename T, typename U, typename std::enable_if<std::is_integral<T>::value && std::is_integral<U>::value && !std::is_signed<T>::value && std::is_signed<U>::value>::type* = nullptr>
                bool operator() (T a, U b) const {
                    if (b < 0){
                        return false;
                    }
                    return static_cast<U>(a) == b;
                }
                template <typename T, typename U, typename std::enable_if<std::is_integral<T>::value && std::is_integral<U>::value && std::is_signed<T>::value && !std::is_signed<U>::value>::type* = nullptr>
                bool operator() (T a, U b) const {
                    if (a < 0){
                        return false;
                    }
                    return a == static_cast<T>(b);
                }

                // Non-integer comparison
                template <typename T, typename U, typename std::enable_if<std::is_floating_point<T>::value && std::is_floating_point<U>::value>::type* = nullptr>
                bool operator() (T a, U b) const { return std::abs(a - b) < ROUND_ERROR; }
                template <typename T, typename U, typename std::enable_if<std::is_integral<T>::value && std::is_floating_point<U>::value>::type* = nullptr>
                bool operator() (T a, U b) const { return std::abs(static_cast<double>(a) - b) < ROUND_ERROR; }
                template <typename T, typename U, typename std::enable_if<std::is_integral<U>::value && std::is_floating_point<T>::value>::type* = nullptr>
                bool operator() (T a, U b) const { return std::abs(a - static_cast<double>(b)) < ROUND_ERROR; }
            };
            struct LessThan {
                template <typename T, typename U, std::enable_if_t<(is_iterable<T> && !is_iterable <U>) && (!std::is_same_v<T, std::string> || !std::is_same_v<U, std::string>), bool> = true>
                bool operator() (T a, U b) const  { return operator()( a.at(0), b); }
                template <typename T, typename U, std::enable_if_t<(!is_iterable<T> && is_iterable <U>) && (!std::is_same_v<T, std::string> || !std::is_same_v<U, std::string>), bool> = true>
                bool operator() (T a, U b) const  { return operator()(a , b.at(0)); }
                template <typename T, typename U, std::enable_if_t<(is_iterable<T> && is_iterable <U>) && (!std::is_same_v<T, std::string> || !std::is_same_v<U, std::string>), bool> = true>
                bool operator() (T a, U b) const  { return operator()(a.at(0) , b.at(0)); }

                template <typename T, std::enable_if_t<!is_iterable<T> && !std::is_same_v<T, std::string>, bool> = true>
                bool operator() (T a, T b) const { return a < b; }
                bool operator()(const std::string& a, const std::string& b) const { return a < b; }

                template <typename T, std::enable_if_t<!is_iterable<T> && !std::is_same_v<T, std::string>, bool> = true>
                bool operator()(std::string, T) const { return false; /* throw std::invalid_argument("Can't compare std::string to non-std::string"); */ }
                template <typename T, std::enable_if_t<!is_iterable<T> && !std::is_same_v<T, std::string>, bool> = true>
                bool operator()(T, std::string) const { return false; /* throw std::invalid_argument("Can't compare std::string to non-std::string"); */ }

                // Integer comparison (extra work to avoid signed/unsigned comparison)
                template <typename T, typename U, typename std::enable_if<std::is_integral<T>::value && std::is_integral<U>::value && std::is_signed<T>::value == std::is_signed<U>::value>::type* = nullptr>
                bool operator() (T a, U b) const { return a < b; }
                template <typename T, typename U, typename std::enable_if<std::is_integral<T>::value && std::is_integral<U>::value && !std::is_signed<T>::value && std::is_signed<U>::value>::type* = nullptr>
                bool operator() (T a, U b) const {
                    if (b < 0){
                        return false;
                    }
                    return static_cast<U>(a) < b;
                }
                template <typename T, typename U, typename std::enable_if<std::is_integral<T>::value && std::is_integral<U>::value && std::is_signed<T>::value && !std::is_signed<U>::value>::type* = nullptr>
                bool operator() (T a, U b) const {
                    if (a < 0){
                        return true;
                    }
                    return a < static_cast<T>(b);
                }

                // Non-integer comparison
                template <typename T, typename U, typename std::enable_if<std::is_floating_point<T>::value && std::is_floating_point<U>::value>::type* = nullptr>
                bool operator() (T a, U b) const { return a < b; }
                template <typename T, typename U, typename std::enable_if<std::is_integral<T>::value && std::is_floating_point<U>::value>::type* = nullptr>
                bool operator() (T a, U b) const { return static_cast<double>(a) < b; }
                template <typename T, typename U, typename std::enable_if<std::is_integral<U>::value && std::is_floating_point<T>::value>::type* = nullptr>
                bool operator() (T a, U b) const { return a < static_cast<double>(b); }
            };
    };
    class AttributeCondition {
        public:
            virtual ~AttributeCondition(){}
            virtual bool matches(const Attribute& attribute) const {(void)attribute; return false;}
            friend std::ostream& operator<<(std::ostream& out, const AttributeCondition&){return out << "unknown condition";}
    };
    class AttributeWild : public AttributeCondition {
        public:
            AttributeWild(const std::string& name) : name_{name} {}
            virtual ~AttributeWild() override {}
            bool matches(const Attribute& attribute) const override {
                return attribute.name() == name_;
            }
            friend std::ostream& operator<<(std::ostream& out, const AttributeWild& condition){
                return out << condition.name_ << "=*";
            }
        private:
            std::string name_;
    };
    template<typename T=double>
    class AttributeRange : public AttributeCondition {
        public:
            // AttributeRange(const std::string& name, T min, T max) : name_{name}, min_{min}, max_{max} {}
            AttributeRange(const std::string& name, T min, T max) : name_{name}, between_{Between<T>(min, max)} {}
            virtual ~AttributeRange() override {}
            bool matches(const Attribute& attribute) const override {
                return attribute.name() == name_ && _visitorAny.apply_visitor(between_, attribute.a_value());
            }
            // friend std::ostream& operator<<(std::ostream& out, const AttributeRange& condition){
            //     // return out name_ << ": " << condition.between_.min_ << " - " << condition.between_.max_;
            //     return out name_ << " between stuff";
            // }
        private:
            template<typename B>
            struct Between  {
                Between(B min, B max) : min_{min}, max_{max} {}
                B min_;
                B max_;

                bool operator()(const std::string&) const { return false; }
                template <typename V, std::enable_if_t<!is_iterable<V> && !std::is_same_v<V, std::string>, bool> = true>
                bool operator() (V a) const { return !(a < min_ || a > max_); }
                template <typename V, std::enable_if_t<is_iterable<V> && !std::is_same_v<V, std::string>, bool> = true>
                bool operator() (V a) const { return !(a.at(0) < min_ || a.at(0) > max_); }
            };

            std::string name_;
            Between<T> between_;
    };
    class BaseVariable;


    template <typename T>
    static Attribute parse_attribute_any(const rapidjson::Value::ConstMemberIterator &itr, T (*get_attr_rj)(const rapidjson::Value &))
    {
        const std::string attr_name = itr->name.GetString();
        const rapidjson::Value &data = itr->value;
        const rapidjson::Value  &value = data["value"];
        if (!value.IsArray())
        {
            return  Attribute{attr_name, get_attr_rj(value) , 1};
        }
        else
        {
            auto arr = value.GetArray();
            size_t attr_len = arr.Size();
            std::vector<T> vec;
            std::transform(arr.begin(),arr.end(),std::back_inserter(vec),[&](auto &v){return get_attr_rj(v);});
              return  Attribute{attr_name,vec , attr_len};
        }

    }
    template <typename T>
    static Attribute parse_attribute_any(const netCDF::NcAtt& att)
    {
        if (att.getAttLength() == 1){
            T v;
            att.getValues(&v);
            const ANY val = v;
            return Attribute{att.getName(), val , 1};
        }
        else
        {
            std::vector<T> value{};
            value.resize(att.getAttLength());
            att.getValues(value.data());
            ANY val;
            size_t att_len_to_pass = att.getAttLength();
            if constexpr(std::is_same_v<T,char>)
            {
                val = std::string(value.cbegin(), value.cend());
                att_len_to_pass = 1;
            }
            else
            {
                val = std::move(value);
            }
            
            return Attribute{att.getName(), val , att_len_to_pass};
        }

    }
    class Product {
        public:
            Product(){}
            Product(std::string name) : name_{name} {}

            Product(std::string name, const std::initializer_list<Attribute>&& attributes) : name_{name}, attributes_{attributes} {}
            Product(std::string name, const std::set<Attribute>& attributes) : name_{name}, attributes_{attributes} {}

            Product(std::string name, const std::set<Attribute>&& attributes, std::vector<std::shared_ptr<AttributeCondition>>&& conditions) : name_{name}, attributes_{attributes}, conditions_{std::move(conditions)} {}
            Product(std::string name, const std::set<Attribute>& attributes, const std::vector<std::shared_ptr<AttributeCondition>>&& conditions) : name_{name}, attributes_{attributes}, conditions_{conditions} {}

            Product(std::string name, const std::vector<std::shared_ptr<AttributeCondition>>& conditions) : name_{name}, conditions_{conditions} {}
            Product(std::string name, std::vector<std::shared_ptr<AttributeCondition>>&& conditions) : name_{name}, conditions_{std::move(conditions)} {}


            // void name(const std::string& name) {name_.assign(name);}
            void name(std::string name) {name_ = name;}
            const std::string& name() const {return name_;}

            auto& attributes() {return attributes_;}
            auto& conditions() {return conditions_;}

            auto& attributes() const {return attributes_;}
            auto& conditions() const {return conditions_;}

            // void add_attribute(std::unique_ptr<Attribute> attr){attributes_.insert(attributes_.end(), std::move(attr));}
            void add_attribute(const Attribute& attr){attributes_.insert(attributes_.end(), attr);}

            friend bool operator==(const Product& me, const Product& other){
                // return me.name() == other.name() &&
                return boost::iequals(me.name(), other.name()) &&
                    std::is_permutation(other.attributes().cbegin(), other.attributes().cend(), me.attributes().cbegin(), me.attributes().cend()) &&
                    std::is_permutation(other.conditions().cbegin(), other.conditions().cend(), me.conditions().cbegin(), me.conditions().cend());
            }
            bool matches(const std::vector<Product>& other) const {
                for (const auto& o : other){
                    if (matches(o)){
                        return true;
                    }
                }
                return false;
            }
            bool matches(const Product& other) const {
                // std::cout << "product.matches\n";
                // std::cout << "other name: " << other.name() << "\n";
                // std::cout << "this name length: " << name().length() << "\n";
                // std::cout << "this name substr: " << name().substr(0, 1) << "\n";
                // std::cout << "this name: " << name() << "\n";
                if (boost::iequals(name(), other.name())){
                // if (name() == other.name()){
                    // std::cout << "names match\n";
                    const auto& match_atts = other.attributes();
                    
                    bool all_matches = std::all_of(match_atts.cbegin(), match_atts.cend(), [this](const auto& o){return attributes_.find(o) != attributes_.end();});
                    // bool all_matches = true;
                    // for (const auto& a : match_atts){
                    //     if (attributes_.find(a) == attributes_.end()){
                    //         std::cout << "Missing " << a << "\n";
                    //         return false;
                    //     }
                    // }
                    if (all_matches){
                        // std::cout << "all attributes matched, checking conditions\n";
                        all_matches = std::all_of(conditions_.cbegin(), conditions_.cend(), [&match_atts](const auto& c){
                            auto i = std::find_if(match_atts.cbegin(), match_atts.cend(), [&c](const auto& o){return c->matches(o);});
                            return i != match_atts.end();
                        });
                    }
                    return all_matches;
                }
                // std::cout << "no match\n";
                return false;
            }
            friend std::ostream& operator<<(std::ostream& out, const Product& product){
                out << product.name_;
                if (product.attributes().size() != 0 || product.conditions().size()){
                    out << "[";
                    bool printed_first = false;
                    for (auto it = product.attributes().cbegin(); it != product.attributes().cend(); ++it){
                        if (printed_first){
                            out << ',';
                        }
                        out << *it;
                        printed_first = true;
                    }
                    for (auto it = product.conditions().cbegin(); it != product.conditions().cend(); ++it){
                        if (printed_first){
                            out << ',';
                        }
                        out << **it;
                        printed_first = true;
                    }
                    out << "]";
                }
                return out;
            }
            void variable(BaseVariable* variable) {variable_ = variable;}
            BaseVariable*  variable() const {return variable_;}

            static size_t skip_passed_char(const std::string& input, size_t pos, char c, bool in_quotes=false){
                while (pos < input.size()){
                    if (input[pos] == c){
                        return pos + 1;
                    } else if (input[pos] == '\\'){
                        pos += 2;
                    } else if (in_quotes){
                        pos++;
                    } else {
                        switch (input[pos]){
                            case '\'':
                                pos = skip_passed_char(input, pos + 1, '\'', true);
                                break;
                            case '"':
                                pos = skip_passed_char(input, pos + 1, '"', true);
                                break;
                            default:
                                pos++;
                        }
                    }
                }
                return std::string::npos;
            }
            static std::vector<Product> parse_list(const std::string& input){
                std::vector<Product> ret{};

                size_t product_start = 0;
                auto current_pos = product_start;
                while (current_pos != std::string::npos && current_pos < input.size()){
                    if (input[current_pos] == ','){
                        ret.push_back(parse(input.substr(product_start, current_pos - product_start)));
                        current_pos++;
                        product_start = current_pos;
                    } else if (input[current_pos] == '\''){
                        current_pos = skip_passed_char(input, current_pos + 1, '\'', true);
                    } else if (input[current_pos] == '"'){
                        current_pos = skip_passed_char(input, current_pos + 1, '"', true);
                    } else if (input[current_pos] == '['){
                        current_pos = skip_passed_char(input, current_pos + 1, ']');
                    } else {
                        current_pos++;
                    }
                }
                if (current_pos != product_start){
                    ret.push_back(parse(input.substr(product_start, current_pos - product_start)));
                }
                {
                    std::unordered_set<std::string> prodnames;
                    for(const auto &pr : ret)
                    {
                        prodnames.insert(pr.name());
                    }
                    for (const auto & pr : default_products)
                    {
                        if(prodnames.count(pr) == 0)
                        {
                            ret.push_back(Product(pr));
                        }
                    }
                }
                return ret;
            }
            static Product parse(const std::string& input){
                auto attribute_bracket = input.find_first_of('[');

                size_t name_start = 0;
                while (input[name_start] == ' '){
                    name_start++;
                }
                size_t name_end = attribute_bracket;
                if (attribute_bracket == std::string::npos){
                    name_end = input.size();
                }
                while (input[name_end - 1] == ' '){
                    name_end--;
                }

                if (attribute_bracket == std::string::npos){
                    if (name_start == 0 && name_end == attribute_bracket){
                        return Product{input};
                    }
                    return Product{input.substr(name_start, name_end - name_start)};
                }
                std::string name(input, name_start, name_end - name_start);

                std::set<Attribute> attributes{};

                auto param_start = attribute_bracket + 1;
                auto current_pos = param_start;
                while (current_pos != std::string::npos && current_pos < input.size()){
                    if (input[current_pos] == ']'){
                        if (current_pos != param_start){
                            attributes.insert(Attribute::parse(input.substr(param_start, current_pos - param_start)));
                        }
                        break;
                    }
                    if (input[current_pos] == ','){
                        attributes.insert(Attribute::parse(input.substr(param_start, current_pos - param_start)));
                        current_pos++;
                        param_start = current_pos;
                    } else if (input[current_pos] == '\''){
                        current_pos = skip_passed_char(input, current_pos + 1, '\'', true);
                    } else if (input[current_pos] == '"'){
                        current_pos = skip_passed_char(input, current_pos + 1, '"', true);
                    } else {
                        current_pos++;
                    }
                }
                return Product{name, attributes};
            }
        private:
            std::string name_{"unspecified"};
            std::set<Attribute> attributes_{};
            std::vector<std::shared_ptr<AttributeCondition>> conditions_{}; // pointer for hierarchy, shared for copy-able
            std::set<Attribute> configuration_{};
            BaseVariable* variable_{nullptr};
            const static std::unordered_set<std::string> default_products;

            // static bool attributes_equal( const std::unique_ptr<Attribute>& left, const std::unique_ptr<Attribute>& right){return *left == *right;}
            // static bool attributes_equal( const Attribute& left, const Attribute& right){return left == right;}
            // static bool attributes_equal( const std::unique_ptr<Attribute>& left, const std::unique_ptr<Attribute>& right){return *left == *right;}
            // static bool attributes_equal( const std::unique_ptr<Attribute>& left, const std::unique_ptr<Attribute>& right){return *left == *right;}


    };
} // namespace focs

#endif // FOCS_PRODUCT

