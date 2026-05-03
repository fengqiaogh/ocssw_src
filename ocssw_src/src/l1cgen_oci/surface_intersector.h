#ifndef OCSSW_SURFACE_INTERSECTOR_H
#define OCSSW_SURFACE_INTERSECTOR_H
#include <CGAL/version.h>
#if CGAL_VERSION_MAJOR >= 6
#include <CGAL/AABB_traits_3.h>
#include <CGAL/AABB_triangle_primitive_3.h>
#define AABB_traits AABB_traits_3
#define AABB_triangle_primitive AABB_triangle_primitive_3
#else
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>
#define AABB_traits AABB_traits
#define AABB_triangle_primitive AABB_triangle_primitive
#endif

#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>


#include <vector>
#include <iostream>
#include <optional>
#include <array>
#include "genutils.h"
#include <tuple>
#include "utils.h"
namespace srf {
typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3 Point;
typedef K::Triangle_3 Triangle;
typedef K::Ray_3 Ray;
typedef std::vector<Triangle>::iterator Iterator;
typedef CGAL::AABB_triangle_primitive<K, Iterator> Primitive;
typedef CGAL::AABB_traits<K, Primitive> AABB_triangle_traits;
typedef CGAL::AABB_tree<AABB_triangle_traits> Tree;
constexpr float pace_altitude = 676.5; // doesn't have to be exact, just far enough to compute the intersection between the surface and the ray
Point elevated_position(float latitude, float longitude, float height);
Point sensor_position(float sensor_zenith_angle, float sensor_azimuth_angle, float latitude, float longitude, float height);
std::tuple<std::array<double, 3>, std::array<double, 3>, std::array<double, 3>>
get_local_vectors(float latitude, float longitude);
std::tuple<double, double, double> get_lat_lon_from_xyz(double x, double y, double z);
std::array<std::complex<double>, 4> solve_depressed_quartic(double a, double b, double c);
double solve_cubic(double a, double b, double c);
bool is_bad(const Point& p);
class SurfaceIntersector {
    std::vector<Triangle> triangles;
    Tree tree;
    // [[maybe_unused]] size_t lines, pixels;
    size_t lines, pixels;

   public:
    size_t get_tree_size() const { return triangles.size(); }
    explicit SurfaceIntersector(const std::vector<float>& x, const std::vector<float>& y,
                                const std::vector<float>& z, size_t lines, size_t pixels,
                                const std::vector<int8_t>& mask = {});
    explicit SurfaceIntersector(const std::vector<Point>& points, size_t lines, size_t pixels,
                                const std::vector<int8_t>& mask = {});
    std::optional<Point> findIntersection(const Point& start, const Point& end) const;
    std::vector<Point> findAllIntersections(const Point& start, const Point& end) const;
    bool doesIntersect(const Point& start, const Point& end) const;
    std::optional<std::pair<Point, std::size_t>> findIntersectionIndex(const Point& start, const Point& end);
    const Triangle& getTriangle(std::size_t index) const {
        return triangles[index];
    }
};
}  // namespace srf

#endif  // OCSSW_SURFACE_INTERSECTOR_H
