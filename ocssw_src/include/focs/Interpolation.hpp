#ifndef INTERP_HPP
#define INTERP_HPP

#include <boost/date_time/posix_time/posix_time.hpp>
#include <cmath>
#include <string>
#include <netcdf>
#include "focs/Log.hpp"
#include "focs/Common.hpp"

namespace interp {

typedef double (*interpolation_method2D)(double *x, double *y, size_t n, size_t ix, size_t iy, double *lut,
                                         double x0, double y0);

// template<typename T>
// using TouchCallBack = void (*)(float*, const double&, std::vector<T >);
// needs to be replaced with a general template function for N-dimensions
constexpr double bad_float = BAD_FLT;
double bilinear2D(double *x, double *y, size_t n, size_t ix, size_t iy, double *lut, double x0, double y0);

const std::map<std::string, interpolation_method2D> methods = {{"bilinear", &bilinear2D}};

class Interpolator {
   public:
    Interpolator() = default;
};

template <class T>
class TimeInterpolator {
   private:
    std::vector<T> *interpolators;
    std::vector<double> times;
    std::vector<double> weights;
    size_t first{}, second{};

   public:
    TimeInterpolator() = default;

    TimeInterpolator(std::vector<T> &interpolators, std::vector<double> &times, double current_time)
        : times(times) {
        this->interpolators = &interpolators;
        if (times.size() != interpolators.size()) {
            EXIT_LOG(std::cerr << "--Error--: size of the time vector and interpolatros vector are different. ")
        }
        for (double &t : times) {
            weights.push_back(std::abs(1.0 / (t - current_time)));
        }
        first = 0;
        second = times.size() - 1;
    }
    template <typename... Args>
    double operator()(Args... args) {
        return (interpolators->operator[](first)(args...) * weights[first] +
                interpolators->operator[](second)(args...) * weights[second]) /
               (weights[first] + weights[second]);
    }
};

// needs to be replaced with a general template class for N-dimensions
class RegularGridInterpolator2D : public Interpolator {
   private:
    double delta_x;
    double delta_y;
    double start_x;
    double start_y;
    std::vector<double> x_grid;
    std::vector<double> y_grid;
    double end_x;
    double end_y;
    size_t nx;
    size_t ny;
    std::vector<double> lut_data;
    netCDF::NcVar *nc_lut;
    interpolation_method2D method;

   public:
    RegularGridInterpolator2D();

    template <typename T, typename U>
    RegularGridInterpolator2D(T *lut, U *x, U *y, size_t nx, size_t ny,
                              const std::string &methodstr = "bilinear")
        : nx(nx), ny(ny) {
        if constexpr (std::is_arithmetic_v<T> || std::is_arithmetic_v<std::decay<T>>) {
            lut_data.resize(nx * ny);
            std::transform(lut, lut + lut_data.size(), lut_data.data(),
                           [](T var) { return static_cast<double>(var); });
        } else if constexpr (std::is_same_v<std::decay<T>, netCDF::NcVar> ||
                             std::is_same_v<T, netCDF::NcVar>) {
            nc_lut = lut;  // to be fully implemented, reading directly from netcdf for binning?
        } else {
            EXIT_LOG(std::cerr << "--Error--: passed type for the luts is not supported. Exiting.")
        }
        x_grid.resize(nx);
        y_grid.resize(ny);
        std::transform(x, x + x_grid.size(), x_grid.begin(), [](U var) { return static_cast<double>(var); });
        std::transform(y, y + y_grid.size(), y_grid.begin(), [](U var) { return static_cast<double>(var); });
        delta_x = static_cast<double>((x[nx - 1] - x[0]) / (nx - 1));
        delta_y = static_cast<double>((y[ny - 1] - y[0]) / (ny - 1));
        start_x = x[0];
        start_y = y[0];
        end_x = x[nx - 1];
        end_y = y[ny - 1];
        method = methods.at(methodstr);
    }

    double operator()(double x, double y);
};
template <typename T>
size_t search(T *arr, size_t s, size_t e, T val, size_t *i_s, size_t *i_e) {
    const bool acsending = arr[s] < arr[e];
    if (acsending) {
        if (val >= arr[e]) {
            *i_s = e;
            *i_e = e;
            return e;
        }
        if (val <= arr[s]) {
            *i_s = s;
            *i_e = s;
            return s;
        }
    } else {
        if (val <= arr[e]) {
            *i_s = e;
            *i_e = e;
            return e;
        }
        if (val >= arr[s]) {
            *i_s = s;
            *i_e = s;
            return s;
        }
    }
    while (e - s > 1) {
        size_t m = (s + e) / 2;  // compute median
        if (acsending) {
            if (arr[m] <= val) {
                s = m;
            } else {
                e = m;
            }
        } else {
            if (arr[m] <= val) {
                e = m;
            } else {
                s = m;
            }
        }
    }
    {
        *i_s = s;
        *i_e = s + 1;
    }
    return s;
}

template <typename T, typename U, size_t N>
class Interp {
   private:
    std::array<size_t, N> dimensions;
    std::array<size_t, N> point;
    std::array<size_t, N> point_s;
    std::array<size_t, N> point_e;
    std::array<size_t, N> pointer_shifts;
    T *lut;
   std::vector<U*>  grid;

   public:
    Interp() = default;

    template <typename... Args>
    Interp(T *lut,std::vector<U*> &  grid, Args &&...dims) : lut(lut), grid(grid) {
        if constexpr (sizeof...(Args) != N) {
            exit(EXIT_FAILURE);
        }
        set_dimensions(dims...);
        set_shifts<N>();
    }

    template <typename Arg, typename... Args>
    void set_dimensions(Arg &&dim, Args &&...dims) {
        constexpr size_t index = N - sizeof...(Args) - 1;
        dimensions[index] = static_cast<size_t>(dim);
        if constexpr (index != N - 1)
            set_dimensions(dims...);
    }

    template <typename Arg, typename... Args>
    void get_point(Arg &&val, Args &&...vals) {
        constexpr size_t index = N - sizeof...(Args) - 1;
        point[index] = search(grid[index], 0, dimensions[index] - 1, val, &point_s[index], &point_e[index]);
        if constexpr (index != N - 1)
            get_point(vals...);
    }

    template <typename Arg, typename... Args>
    T get_value(T *lut_ptr, Arg &&val, Args &&...vals) {
        constexpr size_t index = N - sizeof...(Args) - 1;
        U w1 = val - grid[index][point_s[index]];
        U w2 = w1;
        if (point_s[index] != point_e[index])
            w2 = grid[index][point_e[index]] - val;
        T u1, u2;
        if constexpr (sizeof...(Args) == 0) {
            u1 = lut_ptr[point_s[index]];
            if (point_s[index] != point_e[index])
                u2 = lut_ptr[point_e[index]];
            else
                u2 = u1;

        } else {
            u1 = get_value(lut_ptr + pointer_shifts[index] * point_s[index], vals...);
            if (point_s[index] != point_e[index])
                u2 = get_value(lut_ptr + pointer_shifts[index] * point_e[index], vals...);
            else
                u2 = u1;
        }
        U inw1, inw2;
        inw1 = 1.0 / w1;
        inw2 = 1.0 / w2;
        if (w1 == 0) {
            inw1 = 1.0;
            inw2 = 0;
        } else if (w2 == 0) {
            inw1 = 0;
            inw2 = 1.0;
        }

        return (inw1 * u1 + inw2 * u2) / (inw1 + inw2);
    }

    template <typename... Args>
    T operator()(Args &&...vals) {
        if constexpr (sizeof...(Args) != N) {
            exit(EXIT_FAILURE);
        }
        get_point(vals...);
        return get_value(lut, vals...);
    }

    template <size_t K>
    void set_shifts() {
        if constexpr (N - K == 0) {
            pointer_shifts[K - 1] = 1;
        } else {
            pointer_shifts[K - 1] = pointer_shifts[K] * dimensions[K];
        }
        if constexpr (K > 0)
            set_shifts<K - 1>();
    }

    template <typename... Args>
    friend Interp interp(T *lut,std::vector<U*> &  grid, Args &&...dims);
};

template <typename T, typename U, typename... Args, size_t K = sizeof...(Args)>
Interp<T, U, K> interp(T *lut,std::vector<U*> &  grid, Args &&...dims) {
    return Interp<T, U, K>(lut, grid, dims...);
}

}  // namespace interp

#endif