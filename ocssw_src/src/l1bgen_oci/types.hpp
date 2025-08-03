

#include <boost/multi_array.hpp>

template <typename T>
using matrix2D = boost::multi_array<T, 2>;

template <typename T>
using vec3D = std::vector<std::vector<std::vector<T>>>;
template <typename T>
using vec2D = std::vector<std::vector<T>>;