#ifndef FOCS_STOPWATCH
#define FOCS_STOPWATCH

#include <chrono>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <type_traits>

#include <iostream>
namespace focs {

/*!
    \class focs::Stopwatch

    \brief Simple timer for profiling sections of code.

    This class uses `std::chrono::steady_clock`, so marks are inherently
    arbitrary.  To get the elapsed time, it simply uses `steady_clock::now() - mark`.

    All `time_point` types are `steady_clock::time_point`.

    \section stopwatch-ex1 Example
    \snippet tests/stopwatch.cpp Stopwatch example

    \section stopwatch-ex2 Reconstructing Stopwatches
    \snippet tests/stopwatch.cpp Stopwatch reconstructing stopwatches
*/
class Stopwatch {
using clock = std::chrono::steady_clock;
using time_point = clock::time_point;

public:
    /*! \brief Reset stopwatch to current (arbitrary) time. */
    void reset() noexcept {
        mark_ = clock::now();
    }

    /*! \brief Reset stopwatch to specified mark.
         This function is mostly for testing.

         Due to implementation, the time set is actually `steady_clock::now() -
         t` so that, e.g., setting it to a `time_point` constructed from a
         `duration` of 10 seconds will actually cause `seconds()` to return 10.

        \param[in] t New interval start time.
     */
    void reset(const time_point t) noexcept {
        mark_ = time_point{clock::now() - t};
    }

    /*! \brief Reset stopwatch to specified mark.
         This function is mostly for testing.

         Due to implementation, the time set is actually `steady_clock::now() -
         t` so that, e.g., setting it to a `duration` of 10 seconds will actually
         cause `seconds()` to return 10.

        \param[in] d New interval start time.  Must be a valid argument for
         `chrono::seconds` constructor, e.g., another `chrono::duration` or a
         number.
     */
    template<typename T>
    void reset(const T d) noexcept {
        reset(time_point{std::chrono::seconds{d}});
    }

    /*! \brief Return elapsed time in the duration specified via the template.
         Use `count()` to get the actual value of the returned duration).

        \tparam T Duration type to cast the elapsed time to.
     */
    template<class T = clock::duration>
    auto now() const noexcept {
        return std::chrono::duration_cast<T>(clock::now() - mark_);
    }
    /*! \brief Return elapsed time in nanoseconds. */
    auto ns() const noexcept {
        return std::chrono::duration_cast<std::chrono::nanoseconds>(clock::now() - mark_).count();
    }
    /*! \brief Return elapsed time in microseconds. */
    auto us() const noexcept {
        return std::chrono::duration_cast<std::chrono::microseconds>(clock::now() - mark_).count();
    }
    /*! \brief Return elapsed time in milliseconds. */
    auto ms() const noexcept {
        return std::chrono::duration_cast<std::chrono::milliseconds>(clock::now() - mark_).count();
    }
    /*! \brief Return elapsed time in seconds. */
    auto seconds() const noexcept {
        return std::chrono::duration_cast<std::chrono::seconds>(clock::now() - mark_).count();
    }
    /*! \brief Return elapsed time in seconds (alias for `seconds()`). */
    auto sec() const noexcept {
        return seconds();
    }
    /*! \brief Return elapsed time in seconds (alias for `seconds()`). */
    auto s() const noexcept {
        return seconds();
    }
    /*! \brief Return elapsed time in minutes. */
    auto m() const noexcept {
        return std::chrono::duration_cast<std::chrono::minutes>(clock::now() - mark_).count();
    }
    /*! \brief Return elapsed time in hours. */
    auto h() const noexcept {
        return std::chrono::duration_cast<std::chrono::hours>(clock::now() - mark_).count();
    }
    /*! \brief Return elapsed time in days. */
    auto d() const noexcept {
        return h() / 24.;
    }
    /*! \brief Return the mark in its native format.

        This functionality should not be relied upon, as the API is no longer within the scope of this class.
     */
    time_point mark(){
        return mark_;
    }

    /*! \brief Pretty-print the elapsed time.

        \param[in] decimals Optional, how many decimal places to print seconds to (default = 6).
     */
    std::string pretty(int decimals=6) const noexcept {
        return pretty_print(ns(), decimals);
    }

    /*!
        \brief Pretty print time elapsed to stream
    
        \param out Target output stream
        \param in Stopwatch to print
    
        \return Given output stream, for chaining
    */
    friend std::ostream& operator<<(std::ostream& out, const Stopwatch& in) {
        return out << pretty_print(in.ns());
    }

    /*! \brief Pretty-print the given duration as \<days>:\<hours>:\<minutes>:\<seconds>.\<decimals>.

        Days, hours, etc, will be omitted for times too small.  This
        functionality should be placed elsewhere, probably.  Some kind of
        TimeUtils library or into focs::TimePoint itself, maybe.

        \param[in] duration The duration to print.
        \param[in] decimals Optional, how many decimal places to print seconds to (default = 6).
     */
    template<typename T>
    static std::string pretty_print(std::chrono::duration<T> duration, int decimals=6) noexcept {
        return pretty_print( std::chrono::duration_cast<std::chrono::nanoseconds>(duration), decimals);
    }
    /*! \brief Pretty-print the given nanoseconds as \<days>:\<hours>:\<minutes>:\<seconds>.\<decimals>.

        Days, hours, etc, will be omitted for times too small.  This
        functionality should be placed elsewhere, probably.  Some kind of
        TimeUtils library or into focs::TimePoint itself, maybe.

        \param[in] nanoseconds Nanoseconds to print.
        \param[in] decimals Optional, how many decimal places to print seconds to (default = 6).
     */
    static std::string pretty_print(long nanoseconds, int decimals=6) noexcept {
        std::ostringstream os;
        if (nanoseconds < 0){
            nanoseconds *= -1;
            os << '-';
        }
        const auto fraction = nanoseconds % 1'000'000'000;
        nanoseconds = (nanoseconds - fraction) / 1'000'000'000;
        const auto seconds = nanoseconds % 60;
        nanoseconds = (nanoseconds - seconds) / 60;
        const auto minutes = nanoseconds % 60;
        nanoseconds = (nanoseconds - minutes) / 60;
        const auto hours = nanoseconds % 24;
        nanoseconds = (nanoseconds - hours) / 24;

        if (nanoseconds > 0){ // days
            os << nanoseconds << ':';
        }
        if (os.tellp() || hours > 0){ // hours
            if (os.tellp() > 0){
                os << std::setfill('0') << std::setw(2);
            }
            os << hours << ':';
        }
        if (os.tellp() || minutes > 0){ // minutes
            if (os.tellp() > 0){
                os << std::setfill('0') << std::setw(2);
            }
            os << minutes << ':';
        }

        if (os.tellp()){
            os << std::setfill('0') << std::setw(2) << seconds;
        } else {
            os << seconds;
        }

        // .fraction
        if (decimals > 0){
            os << '.' << std::setfill('0') << std::setw(decimals);
            if (decimals == 9){
                os << fraction;
            } else {
                os << static_cast<unsigned long>(fraction * std::pow(10, decimals - 9));
            }
        }

        return os.str();
    }

private:
    time_point mark_{clock::now()}; /*!< Start of time interval. */

};

} // namespace focs

#endif // FOCS_STOPWATCH

