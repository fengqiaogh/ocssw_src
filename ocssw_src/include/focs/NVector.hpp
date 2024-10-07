#ifndef FOCS_NVECTOR_HPP
#define FOCS_NVECTOR_HPP
namespace focs {
template <typename T, size_t Dims>
struct NVector {
    typedef std::vector<typename NVector<T, Dims - 1>::type> type;
    int dimensions = Dims;
};
template <typename T>
struct NVector<T, 0> {
    typedef T type;
};
}  // namespace focs
// usage: NVector<2, double>::type == 2D vector of doubles
#endif  // FOCS_NVECTOR_HPP
