#ifndef SEARCH_L1C_H
#define SEARCH_L1C_H
#include "genutils.h"
#include <string>
#include "utils.h"
#include <boost/geometry.hpp>
#include <boost/geometry/index/rtree.hpp>
namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;
using point_3d = bg::model::point<double, 3, bg::cs::cartesian>;
using value_t = std::pair<point_3d, std::size_t>;
point_3d get_norm_vector(double lat, double lon);
double haversineDistance(double lat1, double lon1, double lat2, double lon2, double radius = EarthRadius);
class SearchL1C {
   private:
    size_t n_lines{};
    size_t n_pixels{};
    std::vector<double> lats_c{}, lons_c{};
    bgi::rtree<value_t, bgi::quadratic<16>> rtree_setl1c;
    double max_dist{std::numeric_limits<double>::min()};

   public:
    SearchL1C() = default;
    explicit SearchL1C(const std::vector<float>& lat_data_c, const std::vector<float>& lon_data_c,
                       size_t n_lines, size_t n_pixels);

    [[nodiscard]] std::vector<int> query_lat_lon(const std::vector<float>& lat,
                                                 const std::vector<float>& lon) const;
    ~SearchL1C() = default;
};
#endif  // SEARCH_L1C_H
