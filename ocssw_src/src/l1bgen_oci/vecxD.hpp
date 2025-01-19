#ifndef __VECXD__
#define __VECXD__

#define CPLUSPLUS11 201103L

#if __cplusplus >= CPLUSPLUS11
template <typename T>
using vec3D = std::vector<std::vector<std::vector<T>>>;
template <typename T>
using vec2D = std::vector<std::vector<T>>;
#else
template <typename T>
typedef std::vector<std::vector<std::vector<T>>> vec3D;
template <typename T>
typedef std::vector<std::vector<T>> vec2D;
#endif

#endif