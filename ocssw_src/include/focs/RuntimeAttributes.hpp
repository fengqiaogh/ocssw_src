#ifndef RUNTIMEATTRIBUTES
#define RUNTIMEATTRIBUTES
#include <focs/DataProvider.hpp>
#include <focs/Product.hpp>
#include <functional>
#include <optional>
#include <stack>
namespace focs {

class RuntimeAttribute : public focs::Attribute {
   protected:
    std::unordered_map<std::string, focs::BaseVariable *> available_products;
    std::string _run_time_name{};
    std::string _group_name{};
    std::string _product_name{};
    std::vector<std::string> needed_products{};
    bool are_needed_products_found{false};
    void check_needs();

   public:
    bool is_valid_attribute() const { return are_needed_products_found; };
    RuntimeAttribute() = default;
    RuntimeAttribute(const std::string &name,
                     DataProviderConfiguration &configuration);
    const std::string &get_name() const { return _run_time_name; };
    const std::string &get_grp_name() const { return _group_name; };
    std::unique_ptr<RuntimeAttribute> init_attribute(
        const std::string &name, DataProviderConfiguration &configuration);
    virtual ~RuntimeAttribute() = default;
};

class GroundControlPoints : public focs::RuntimeAttribute {
   private:
    std::vector<float> gringpointlatitude{};
    std::vector<float> gringpointlongitude{};
    std::vector<int32_t> gringpointsequence{};
    std::vector<std::vector<float>> geobox{};
    std::string geospatial_bounds; // wkt format
    std::optional<float> last_lat;
    size_t geobox_cnt{0};
    float geo_box_lat_lim{20.0f};
    void set_gring();

   public:
    GroundControlPoints() = default;
    GroundControlPoints(DataProviderConfiguration &configuration);
    void geo_box(TileParameters& tile);
    void set_gringpointlatitude();
    void set_gringpointlongitude();
    void set_gringpointsequence();
    void set_geobounds();
};

class DayNightAttribute : public focs::RuntimeAttribute {
   private:
    // enum daytime { DAY = 'Day', NIGHT = 'Night', MIXED = 'Mixed', UNKNOWN =
    // 'Unknown' };
    std::string state{"Unknown"};
    float solzen_night = 90;
    std::string get_daytime(float solzen) {
        if (solzen > solzen_night)
            return "Night";
        else
            return "Day";
    }

   public:
    DayNightAttribute() = default;
    DayNightAttribute(DataProviderConfiguration &configuration);
    void scan_solzen();
    void set_flag() { a_value_ = state; }

    // {"Day", "Night", "Mixed", "Unknown"};
};

class GeoSpatialBounds : public focs::RuntimeAttribute {
   private:
    float min_lat{90.0f}, max_lat{-90.0f}, min_lon{180.f}, max_lon{-180.0f};
    std::unordered_map<std::string, float> geospatial_bounds{};
    template <typename T>
    std::pair<T, T> find_min_max(const std::vector<T> &vec) {
        T max_ = *std::max_element(vec.begin(), vec.end());
        T min_ = *std::min_element(vec.begin(), vec.end());
        return {max_, min_};
    }

   public:
    GeoSpatialBounds() = default;
    GeoSpatialBounds(DataProviderConfiguration &configuration);
    void scan_lat_lon();
    void set_file_attrs(std::function<void(focs::Attribute)> set_attrs);
};
}  // namespace focs

#endif