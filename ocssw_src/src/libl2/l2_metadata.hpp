#ifndef L2_METADATA_HPP
#define L2_METADATA_HPP
#include <netcdf>
#include <variant>
#include <string>
#include <map>
#include "l2_utils.hpp"
#include "timeutils.h"
#include "optional"
class MetaL2 {
   public:
    using NCAttribute =
        std::variant<int, float, std::string, double, short, char, long long, std::vector<int>,
                     std::vector<float>, std::vector<double>, std::vector<short>, std::vector<long long>>;
    std::map<std::string, NCAttribute> attributes;
    int16_t syear{0};
    int16_t sday{0};
    double smsec{0};
    int16_t eyear{0};
    int16_t eday{0};
    double emsec{0};
    int16_t nrec{0};
    int16_t nsamp{0};
    std::string flag_names;
    std::string source;
    std::string sensor_name;
    std::string mission;
    std::string title;
    float westlon{-180};
    float eastlon{180};
    float southlat{-90};
    float northlat{90};
    std::vector<int> bits;
    std::vector<double> scan_time;
    std::vector<int> year;
    std::vector<int> day;
    std::vector<int> msec;
    MetaL2() = default;
    std::string time_coverage_start, time_coverage_end;

    /**
     * @brief Reads metadata attributes from a netCDF file into the MetaL2 object
     * @param file Reference to the netCDF file to read from
     * @details Reads global attributes and group-specific attributes/variables from the file.
     *          Handles different netCDF data types and stores them in the attributes map.
     *          Also reads scan time information from scan_line_attributes group if present.
     */
    void read(netCDF::NcFile& file);
    
    /**
     * @brief Constructs a MetaL2 object from a netCDF file
     * @param file Reference to the netCDF file to read metadata from
     * @details Reads metadata from file and initializes object members.
     *          Processes time coverage information, geospatial bounds,
     *          instrument details, and various other metadata attributes.
     */
    MetaL2(netCDF::NcFile& file);

    template <typename T>
    std::optional<T> get_attribute(const std::string& name) {
        auto it = attributes.find(name);
        if (it == attributes.end()) {
            return {};
        }
        return std::get<T>(it->second);
    }

    template <typename T>
    void set_attribute(const std::string& name, T value) {
        attributes[name] = value;
    }
};
#endif