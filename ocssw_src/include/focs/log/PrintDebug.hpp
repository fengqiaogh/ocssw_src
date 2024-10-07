#ifndef PRINT_DEBUG_HPP
#define PRINT_DEBUG_HPP
#include <iostream>
template <typename, typename = void>
constexpr bool is_iterable{};

template <typename T>
constexpr bool is_iterable<T, std::void_t<decltype(std::declval<T>().begin()),
                                          decltype(std::declval<T>().end())>> =
    true;

template <typename T, typename U>
std::ostream& operator<<(std::ostream& os,
                         const std::pair<T, U>& data) noexcept {
    os << data.first << " : " << data.second;
    return os;
}
template <typename T, std::enable_if_t<is_iterable<T>, bool> = true>
std::ostream& operator<<(std::ostream& os, const T& data) noexcept {
    for (auto it = data.begin(); it != data.end(); it++) {
        os << *it << " ";
    }

    return os;
}

// An iterator trait those value_type is the value_type of the iterated
// container, supports even back_insert_iterator (where value_type is void)

template <typename T, typename = void>
struct iterator_trait : std::iterator_traits<T> {};

template <typename T>
struct iterator_trait<T, std::void_t<typename T::container_type>>
    : std::iterator_traits<typename T::container_type::iterator> {};


inline void print(std::ostream& stream, const char* format) { stream << format; }

template <typename T, typename... Targs>
void print(std::ostream& stream, const char* format, T&& value,
           Targs&&... fargs) {
    for (; *format != '\0'; format++) {
        if (*format == '%') {
            stream << value;
            print(stream, format + 1, fargs...);
            return;
        }
        stream << *format;
    }
}

inline void print(const char* format) { std::cout << format; }

template <typename T, typename... Targs>
void print(const char* format, T&& value, Targs&&... fargs) {
    print(std::cout, format, value, fargs...);
}

#endif