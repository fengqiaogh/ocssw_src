#include "surface_intersector.h"

#if CGAL_VERSION_MAJOR >= 6
#define get_cgal std::get_if
#else
#define get_cgal boost::get
#endif


namespace srf {
/**
 * @brief Construct a 3D triangel mesh from input vertices with the coordinates x,y and z
 *
 *
 * @param x       X-coordinates of vertices in ECEF frame (kilometers), stl vector
 * @param y       Y-coordinates of vertices in ECEF frame (kilometers), stl vector
 * @param z       Z-coordinates of vertices in ECEF frame (kilometers), stl vector
 *                All three vectors must have size = lines × pixels.
 * @param lines   Number of rows in the grid.
 * @param pixels  Number of columns per row in the grid.
 * @param mask    Optional binary validity mask. If provided, must have size = lines × pixels.
 *                mask[i] = 0 marks the vertex as invalid and skipped;

 */
SurfaceIntersector::SurfaceIntersector(const std::vector<float>& x, const std::vector<float>& y,
                                       const std::vector<float>& z, size_t lines, size_t pixels,
                                       const std::vector<int8_t>& mask)
    : lines(lines), pixels(pixels) {
    triangles.reserve(lines * pixels * 2);  // 2 triangles per cell ( 4 vertices )
    const int8_t* ptr = nullptr;
    assert(x.size() == lines * pixels);
    assert(y.size() == lines * pixels);
    assert(z.size() == lines * pixels);
    if (!mask.empty()) {
        assert(mask.size() == lines * pixels);
        ptr = mask.data();
    }
    for (size_t iline = 0; iline < lines; iline++) {
        for (size_t ipixel = 0; ipixel < pixels; ipixel++) {
            size_t index = iline * pixels + ipixel;
            int iline_up = std::min(iline + 1, lines - 1);
            int iline_down = std::max((int)iline - 1, 0);
            int ipixel_left = std::max((int)ipixel - 1, 0);
            int ipixel_right = std::min(ipixel + 1, pixels - 1);
            int index_left_down = iline_down * pixels + ipixel_left;
            int index_left_up = iline_up * pixels + ipixel_left;
            int index_right_up = iline_up * pixels + ipixel_right;
            int index_right_down = iline_down * pixels + ipixel_right;
            // check if all 4 corners are valid coordinates
            if (x[index] == BAD_FLT || y[index] == BAD_FLT || z[index] == BAD_FLT)
                continue;
            if (x[index_left_down] == BAD_FLT || y[index_left_down] == BAD_FLT ||
                z[index_left_down] == BAD_FLT)
                continue;
            if (x[index_left_up] == BAD_FLT || y[index_left_up] == BAD_FLT || z[index_left_up] == BAD_FLT)
                continue;
            if (x[index_right_up] == BAD_FLT || y[index_right_up] == BAD_FLT || z[index_right_up] == BAD_FLT)
                continue;
            if (x[index_right_down] == BAD_FLT || y[index_right_down] == BAD_FLT ||
                z[index_right_down] == BAD_FLT)
                continue;
            if (ptr) {
                if (mask[index] == 0)
                    continue;
            }
            float x_average = (x[index] + x[index_left_down]) / 2;
            float y_average = (y[index] + y[index_left_down]) / 2;
            float z_average = (z[index] + z[index_left_down]) / 2;
            Point p00(x_average, y_average, z_average);
            x_average = (x[index] + x[index_left_up]) / 2;
            y_average = (y[index] + y[index_left_up]) / 2;
            z_average = (z[index] + z[index_left_up]) / 2;
            Point p10(x_average, y_average, z_average);
            x_average = (x[index] + x[index_right_up]) / 2;
            y_average = (y[index] + y[index_right_up]) / 2;
            z_average = (z[index] + z[index_right_up]) / 2;
            Point p11(x_average, y_average, z_average);
            x_average = (x[index] + x[index_right_down]) / 2;
            y_average = (y[index] + y[index_right_down]) / 2;
            z_average = (z[index] + z[index_right_down]) / 2;
            Point p01(x_average, y_average, z_average);

            triangles.emplace_back(Triangle(p00, p10, p11));
            triangles.emplace_back(Triangle(p00, p11, p01));
        }
    }
    tree.insert(triangles.begin(), triangles.end());
    tree.build();
}
/**
 * @brief Construct a 3D triangel mesh from input vertices with the coordinates x,y and z
 *
 *
 * @param points  coordinates of vertices in ECEF frame (kilometers), vector of GDAL points

 * @param lines   Number of rows in the grid.
 * @param pixels  Number of columns per row in the grid.
 * @param mask    Optional binary validity mask. If provided, must have size = lines × pixels.
 *                mask[i] = 0 marks the vertex as invalid and skipped;

 */
SurfaceIntersector::SurfaceIntersector(const std::vector<Point>& points, size_t lines, size_t pixels,
                                       const std::vector<int8_t>& mask)
    : lines(lines), pixels(pixels) {
    triangles.reserve((lines - 1) * (pixels - 1) * 2);  // 2 triangles per cell ( 4 vertices )
    const int8_t* ptr = nullptr;
    assert(points.size() == lines * pixels);
    if (!mask.empty()) {
        assert(mask.size() == lines * pixels);
        ptr = mask.data();
    }
    for (size_t iline = 0; iline < lines; iline++) {
        for (size_t ipixel = 0; ipixel < pixels; ipixel++) {
            size_t index = iline * pixels + ipixel;
            int iline_up = std::min(iline + 1, lines - 1);
            int iline_down = std::max((int)iline - 1, 0);
            int ipixel_left = std::max((int)ipixel - 1, 0);
            int ipixel_right = std::min(ipixel + 1, pixels - 1);
            int index_left_down = iline_down * pixels + ipixel_left;
            int index_left_up = iline_up * pixels + ipixel_left;
            int index_right_up = iline_up * pixels + ipixel_right;
            int index_right_down = iline_down * pixels + ipixel_right;

            if (points[index].x() == BAD_FLT || points[index].y() == BAD_FLT || points[index].z() == BAD_FLT)
                continue;

            if (points[index_left_down].x() == BAD_FLT || points[index_left_down].y() == BAD_FLT ||
                points[index_left_down].z() == BAD_FLT)
                continue;

            if (points[index_left_up].x() == BAD_FLT || points[index_left_up].y() == BAD_FLT ||
                points[index_left_up].z() == BAD_FLT)
                continue;

            if (points[index_right_up].x() == BAD_FLT || points[index_right_up].y() == BAD_FLT ||
                points[index_right_up].z() == BAD_FLT)
                continue;

            if (points[index_right_down].x() == BAD_FLT || points[index_right_down].y() == BAD_FLT ||
                points[index_right_down].z() == BAD_FLT)
                continue;

            if (ptr) {
                if (mask[index] == 0)
                    continue;
            }
            float x_average = (points[index].x() + points[index_left_down].x()) / 2;
            float y_average = (points[index].y() + points[index_left_down].y()) / 2;
            float z_average = (points[index].z() + points[index_left_down].z()) / 2;
            Point p00(x_average, y_average, z_average);
            x_average = (points[index].x() + points[index_left_up].x()) / 2;
            y_average = (points[index].y() + points[index_left_up].y()) / 2;
            z_average = (points[index].z() + points[index_left_up].z()) / 2;
            Point p10(x_average, y_average, z_average);
            x_average = (points[index].x() + points[index_right_up].x()) / 2;
            y_average = (points[index].y() + points[index_right_up].y()) / 2;
            z_average = (points[index].z() + points[index_right_up].z()) / 2;
            Point p11(x_average, y_average, z_average);
            x_average = (points[index].x() + points[index_right_down].x()) / 2;
            y_average = (points[index].y() + points[index_right_down].y()) / 2;
            z_average = (points[index].z() + points[index_right_down].z()) / 2;
            Point p01(x_average, y_average, z_average);
            triangles.emplace_back(Triangle(p00, p10, p11));
            triangles.emplace_back(Triangle(p00, p11, p01));
        }
    }
    tree.insert(triangles.begin(), triangles.end());
    tree.build();
}

/**
 * @brief Find the first ray-surface intersection point.
 *
 *
 * @param start The starting GDAL point of the ray in 3D space (ECEF coordinates, kilometers).
 * @param end   The ending GDAL point of the ray in 3D space (ECEF coordinates, kilometers).
 *
 * @return An optional Point containing:
 *         - The intersection GDAL point in ECEF coordinates if an intersection is found
 *         - std::nullopt if no intersection exists
 */
std::optional<Point> SurfaceIntersector::findIntersection(const Point& start, const Point& end) const {
    Ray ray(start, end);
    if (tree.do_intersect(ray)) {
        if (auto intersection = tree.any_intersection(ray)) {
            auto p = get_cgal<Point>(&(intersection->first));
            if (p) {
                return *p;
            }
        }
    }

    return std::nullopt;
}

/**
 * @brief Find all ray-surface intersection points along a ray segment.
 *
 * @param start The starting GDAL point of the ray in 3D space (ECEF coordinates, kilometers).
 * @param end   The ending GDAL point of the ray in 3D space (ECEF coordinates, kilometers).
 *
 * @return A vector of intersection GDAL points in ECEF coordinates (kilometers), sorted by
 *         proximity along the ray. Returns an empty vector if no intersections are found.
 */
std::vector<Point> SurfaceIntersector::findAllIntersections(const Point& start, const Point& end) const {
    std::vector<Point> intersections;
    Ray ray(start, end);

    std::vector<Tree::Intersection_and_primitive_id<Ray>::Type> all_intersections;
    tree.all_intersections(ray, std::back_inserter(all_intersections));

    for (const auto& intersection : all_intersections) {
        auto p = get_cgal<Point>(&(intersection.first));
        if (p) {
            intersections.push_back(*p);
        }
    }

    return intersections;
}

/**
 * @brief Check if a ray intersects the surface mesh.
 *
 * @param start The starting GDAL point of the ray in 3D space (ECEF coordinates, kilometers).
 * @param end   The ending GDAL point of the ray in 3D space (ECEF coordinates, kilometers).
 *
 * @return true if the ray segment intersects the surface, false otherwise.
 */
bool SurfaceIntersector::doesIntersect(const Point& start, const Point& end) const {
    Ray ray(start, end);
    return tree.do_intersect(ray);
}

/**
 * @brief Find the first ray-surface intersection point and its associated triangle index.
 *
 * @param start The starting GDAL point of the ray in 3D space (ECEF coordinates, kilometers).
 * @param end   The ending GDAL point of the ray in 3D space (ECEF coordinates, kilometers).
 *
 * @return An optional pair containing:
 *         - The intersection GDAL point in ECEF coordinates and its triangle index if an intersection is found
 *         - std::nullopt if no intersection exists
 *
 * @note This function cannot be const because accessing intersection->second requires non-const iterators.
 */
std::optional<std::pair<Point, std::size_t>> SurfaceIntersector::findIntersectionIndex(
    const Point& start, const Point& end) {  // can't define as const because of intersection->second is
    // non-const and triangles.begin() must be non-const

    Ray ray(start, end);

    if (tree.do_intersect(ray)) {
        if (auto intersection = tree.any_intersection(ray)) {
            auto p = get_cgal<Point>(&(intersection->first));
            if (p) {
                // Cast away const to get non-const iterator
                std::size_t primitive_id = std::distance(triangles.begin(), intersection->second);
                return std::make_pair(*p, primitive_id);
            }
        }
    }

    return std::nullopt;
}

/**
 * @brief Solve a general cubic equation of the form: x³ + a*x² + b*x + c = 0
 *
 * @param a Coefficient of x² term from original cubic (used to convert depressed to original)
 * @param b Coefficient of x term from original cubic
 * @param c Constant term from original cubic
 * @return One real root of the cubic equation (other roots may be complex)
 *
 */
double solve_cubic(double a, double b, double c) {
    double p = b - (a * a) / 3;
    double q = (2 * a * a * a) / 27 - (a * b) / 3 + c;
    double discriminant = -((q * q) / 4 + (p * p * p) / 27);

    if (discriminant < 0) {
        double r = sqrt(-discriminant);
        double u = cbrt(-q / 2 + r);
        double v = cbrt(-q / 2 - r);
        return u + v - a / 3;
    } else if (discriminant == 0) {
        double u = cbrt(-q / 2);
        return 2 * u - a / 3;
    } else {
        double phi = acos(3 * q / (2 * p) * sqrt(-3 / p));
        return 2 * sqrt(-p / 3) * cos(phi / 3) - a / 3;
    }
}

/**
 * @brief Solve a depressed quartic equation of the form: x⁴ + a*x² + b*x + c = 0
 *
 * @param a Coefficient of x² term
 * @param b Coefficient of x term
 * @param c Constant term
 * @return Array of 4 complex roots (may be real or complex pairs)
 *
 */
std::array<std::complex<double>, 4> solve_depressed_quartic(double a, double b, double c) {
    std::array<std::complex<double>, 4> roots;

    double A = -a / 2;
    double B = -2 * c / 2;
    double C = (a * c - (b * b) / 4) / 2;
    double y = solve_cubic(A, B, C);
    std::complex<double> R = sqrt(2 * y - a);

    std::complex<double> D = sqrt(-(2 * y + a + (2 * b) / R));
    std::complex<double> E = sqrt(-(2 * y + a - (2 * b) / R));
    roots[0] = (R + D) / 2.0;
    roots[1] = (R - D) / 2.0;
    roots[2] = (-R + E) / 2.0;
    roots[3] = (-R - E) / 2.0;
    return roots;
}

/**
 * @brief Extract geodetic coordinates (latitude, longitude, height) from Cartesian (x, y, z) position.
 * We solve the inverse geodetic problem by minimizing the distance from the input point (x,y,z)
 * to the ellipsoid surface defined by: (X/a)² + (Y/a)² + (Z/b)² = 1
 *
 * Using Lagrange multipliers to minimize F(X,Y,Z,λ) = (X - x)² + (Y - y)² + (Z - z)² + λ·((X/a)² + (Y/a)² +
 * (Z/b)² - 1)
 *
 * ∂F/∂X = ∂F/∂Y = ∂F/∂Z = ∂F/∂λ = 0 yields X = x / (1 + λ/a²), Y = y / (1 + λ/a²), Z = z / (1 + λ/b²) and
 * (X/a)² + (Y/a)² + (Z/b)² = 1
 *
 * We end up with a quartic equation u⁴ + A·u² + B·u + C = 0
 * @param x X-coordinate in ECEF frame (kilometers)
 * @param y Y-coordinate in ECEF frame (kilometers)
 * @param z Z-coordinate in ECEF frame (kilometers)
 * @return Tuple of (latitude in degrees, longitude in degrees, height in kilometers)
 */
std::tuple<double, double, double> get_lat_lon_from_xyz(double x, double y, double z) {
    double a2 = semi_axis_a * semi_axis_a;
    double b2 = semi_axis_b * semi_axis_b;

    // Compute the substitution parameters for converting to depressed quartic form
    double v = (b2 - a2) / 2.0;  // half the difference of axes squared

    double A = -2 * v * v - (x * x + y * y) * a2 - z * z * b2;
    double B = -2 * (x * x + y * y) * a2 * v + 2 * z * z * b2 * v;
    double C = -(x * x + y * y) * a2 * v * v - z * z * b2 * v * v + v * v * v * v;

    auto roots = solve_depressed_quartic(A, B, C);

    double min_distance = std::numeric_limits<double>::max();
    double best_lambda = 0.0;
    for (const auto& root : roots) {
        // Check if this root is real (imaginary part negligibly small)
        if (std::imag(root) == 0) {
            // Convert from depressed quartic variable u back to Lagrange multiplier λ
            double lambda = std::real(root) - (a2 + b2) / 2.0;

            // Compute the surface point (X, Y, Z) from the Lagrange multiplier solution
            double X = x / (1 + lambda / a2);
            double Y = y / (1 + lambda / a2);
            double Z = z / (1 + lambda / b2);

            // Calculate distance from input point (x,y,z) to this surface point (X,Y,Z)
            double distance = std::sqrt((X - x) * (X - x) + (Y - y) * (Y - y) + (Z - z) * (Z - z));

            // Track the root that yields the minimum distance (closest point on ellipsoid)
            if (distance < min_distance) {
                min_distance = distance;
                best_lambda = lambda;
            }
        }
    }

    double X = x / (1 + best_lambda / a2);
    double Y = y / (1 + best_lambda / a2);
    double Z = z / (1 + best_lambda / b2);

    //   Convert Cartesian ECEF coordinates (X, Y, Z) to geodetic coordinates:
    //   latitude = arctan2(Z, sqrt(X² + Y²) (1-e²))
    //   longitude = arctan2(Y, X)
    //   height = distance from input point to surface point
    double lat = std::atan2(Z, std::sqrt(X * X + Y * Y) * (1 - e2)) / OEL_DEGRAD;
    double lon = std::atan2(Y, X) / OEL_DEGRAD;
    double height = std::sqrt((X - x) * (X - x) + (Y - y) * (Y - y) + (Z - z) * (Z - z));

    return {lat, lon, height};
}

/**
 * @brief Compute a local orthonormal coordinate system (East, North, Up) at a geodetic location.
 *
 * @param latitude  Geodetic latitude in degrees (-90 to 90)
 * @param longitude Geodetic longitude in degrees (-180 to 180)
 * @return Tuple containing (vector_zenith, vector_azimuth, vector_east)
 *         All vectors are unit length and mutually orthogonal (form a right-handed frame)
 *         Position is in kilometers in ECEF frame
 */
std::tuple<std::array<double, 3>, std::array<double, 3>, std::array<double, 3>> get_local_vectors(
    float latitude, float longitude) {
    // Precomputing sine and cosine of latitude and longitude (C compilers might not optimize this out)
    double sin_lat = sin(latitude * OEL_DEGRAD);
    double cos_lat = cos(latitude * OEL_DEGRAD);
    double sin_lon = sin(longitude * OEL_DEGRAD);
    double cos_lon = cos(longitude * OEL_DEGRAD);
    // compute the zenith unit vector
    std::array<double, 3> vector_zenith{};
    vector_zenith[0] = cos_lat * cos_lon;
    vector_zenith[1] = cos_lat * sin_lon;
    vector_zenith[2] = sin_lat;
    // Calculate the North vector pointing toward the North Pole in local topocentric frame
    // This is tangent to the meridian and perpendicular to the zenith vector
    // In local coordinates, this points  northward (with some  positive vertical component)
    std::array<double, 3> vector_azimuth{};
    vector_azimuth[0] = -sin_lat * cos_lon;
    vector_azimuth[1] = -sin_lat * sin_lon;
    vector_azimuth[2] = cos_lat;
    // Calculate the East vector perpendicular to both azimuth and zenith (right-hand rule)
    // Points eastward in the local topocentric frame (tangent to the parallel)
    // This is the cross product: vector_azimuth × vector_zenith
    std::array<double, 3> vector_east{};
    vector_east[0] = -sin_lon;
    vector_east[1] = cos_lon;
    vector_east[2] = 0;

    // Return the complete local orthonormal coordinate frame and the surface position
    return {vector_zenith, vector_azimuth, vector_east};
}

/**
 * @brief Calculate the 3D ECEF position of a sensor given its viewing geometry and ground location.
 *
 * @param sensor_zenith_angle  Zenith angle in degrees (0° = straight up, 90° = horizon)
 * @param sensor_azimuth_angle Azimuth angle in degrees (0° = North, 90° = East, 180° = South, 270° = West)
 * @param latitude  Geodetic latitude in degrees of the ground point being observed
 * @param longitude Geodetic longitude in degrees of the ground point being observed
 * @return Point object with ECEF coordinates (x, y, z) in kilometers from Earth's center
 *
 */
Point sensor_position(float sensor_zenith_angle, float sensor_azimuth_angle, float latitude, float longitude,
                      float height) {
    if (std::abs(latitude) > 90.0 || std::abs(longitude) > 180.0 || height == BAD_FLT ||
        sensor_zenith_angle == BAD_FLT || sensor_azimuth_angle == BAD_FLT)
        return {BAD_FLT, BAD_FLT, BAD_FLT};
    // Get the local orthonormal coordinate system (Up/Zenith, North, East) and surface position
    const auto [vector_zenith, vector_azimuth, vector_east] = get_local_vectors(latitude, longitude);
    std::array<double, 3> position = geodetic_to_ECEF<std::array<double, 3>>(latitude, longitude, height);

    double cos_az = cos(sensor_azimuth_angle * OEL_DEGRAD);
    double sin_az = sin(sensor_azimuth_angle * OEL_DEGRAD);
    double cos_zen = cos(sensor_zenith_angle * OEL_DEGRAD);
    double sin_zen = sin(sensor_zenith_angle * OEL_DEGRAD);


    std::array<double, 3> vector_ray{};

    // - cos(zenith_angle) projects toward zenith (upward)
    // - sin(zenith_angle)*cos(azimuth_angle) projects toward north
    // - sin(zenith_angle)*sin(azimuth_angle) projects toward east
    for (size_t i = 0; i < 3; i++) {
        vector_ray[i] = vector_zenith[i] * cos_zen + vector_azimuth[i] * sin_zen * cos_az +
                        vector_east[i] * sin_zen * sin_az;
        vector_ray[i] *= pace_altitude / cos_zen;  // Scale to satellite altitude in kilometers
    }

 
    return {vector_ray[0] + position[0], vector_ray[1] + position[1], vector_ray[2] + position[2]};
}

/**
 * @brief Calculate the 3D ECEF position at a given height above the ellipsoid surface.
 *
 * @param latitude  Geodetic latitude in degrees
 * @param longitude Geodetic longitude in degrees
 * @param height    Height above the ellipsoid surface in kilometers
 * @return a GDAL point representing the 3D ECEF position in kilometers relative to Earth's center
 */
Point elevated_position(float latitude, float longitude, float height) {
    if (std::abs(latitude) > 90.0 || std::abs(longitude) > 180.0 || height == BAD_FLT)
        return {BAD_FLT, BAD_FLT, BAD_FLT};
    return geodetic_to_ECEF<Point>(latitude, longitude, height);
}

bool is_bad(const Point& p) {
    if (p.x() == BAD_FLT)
        return true;
    if (p.y() == BAD_FLT)
        return true;
    if (p.z() == BAD_FLT)
        return true;

    return false;
}

}  // namespace srf