#ifndef FOCS_TIMEPOINT
#define FOCS_TIMEPOINT

#include <cmath>
#include <string>

#include "focs/LeapSecondDatabase.hpp"

#include <boost/date_time/gregorian_calendar.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
namespace tmgrn {
const boost::posix_time::ptime epoch(boost::gregorian::date(1970, 1, 1));

std::string &granule_start_time();

std::string &granule_end_time();

namespace bt = boost::posix_time;
const std::locale formats[] = {
    std::locale(std::locale::classic(), new bt::time_input_facet("%Y-%m-%dT%H:%M:%S")),
    std::locale(std::locale::classic(), new bt::time_input_facet("%Y/%m/%d %H:%M:%S")),
    std::locale(std::locale::classic(), new bt::time_input_facet("%d.%m.%Y %H:%M:%S")),
    std::locale(std::locale::classic(), new bt::time_input_facet("%Y-%m-%dT%H:%M:%SZ")),
    std::locale(std::locale::classic(), new bt::time_input_facet("minutes since %Y-%m-%d %H:%M:%S"))
    };
std::time_t pt_to_time_t(const bt::ptime &pt);

double seconds_from_epoch(const std::string &s);
}  // namespace tmgrn

namespace focs {
    /*!
        \class focs::TimePoint

        \brief Representation of a single point in time, with clock conversions.

        The internal representation is a boost::posix_time::ptime, meaning that
        leap seconds are compensated for manually within this class.  ptime's are 
        stored as a Gregorian date and the time duration that has passed within 
        that date, wrapping to a new date when hitting 23:59:60; this also means 
        that leap seconds can't be represented in UTC (time will "freeze").

        Conversions to and from local time are not supported, largely due to DST.

        \section Example
        \snippet tests/time_point.cpp focs::TimePoint example
    */
    class TimePoint {
        public:
            /*! \enum TimePoint::clock_type

                \brief Type of clock being represented
            */
            enum class ClockType {
                local, ///< Local UTC clock with time zone offset
                utc,   ///< UTC Unix time
                tai,   ///< TAI clock 1958
                tai93, ///< TAI clock 1993
                gps,   ///< GPS clock 1980, Jan 6
                tt     ///< Terrestrial Time clock 1977
            };

            /*!
                \brief Construct a TimePoint with the current time.
            */
            explicit TimePoint() : ptime_{boost::posix_time::microsec_clock::universal_time()} {};

            /*!
                \brief Construct a clock with offset from the clock epic in seconds.

                This is the epic for each clock:

                Clock | Epoch
                ----- | -----
                local | Jan 1, 1970, 00:00:00 (standard Unix epoch)
                utc   | Jan 1, 1970, 00:00:00 (standard Unix epoch)
                tai   | Jan 1, 1958, 00:00:00
                tai93 | Jan 1, 1993, 00:00:00
                gps   | Jan 6, 1980, 00:00:00
                tt    | Jan 1, 1977, 00:00:00

                \param seconds_since_epoch seconds since clock epoch
                \param clock_type what kind of clock
            */
            TimePoint(const double seconds_since_epoch, const TimePoint::ClockType clock_type=ClockType::utc);

            /*!
                \brief Construct a new TimePoint given date, time, and optional input clock_type

                \param date date
                \param time time
            */
            TimePoint(const boost::gregorian::date& date, const boost::posix_time::time_duration& time) : ptime_{date, time} {}

            /*!
                \brief Construct a new TimePoint given date, time plus miliseconds

                \param date date
                \param time time
                \param milli milliseconds
            */
            TimePoint(const boost::gregorian::date date, const boost::posix_time::time_duration& time, const int milli) : ptime_{date, time + boost::posix_time::milliseconds{milli}} {}

            /*!
                \brief Construct a new TimePoint given a boost ptime

                \param ptime time point to represent
            */
            explicit TimePoint(const boost::posix_time::ptime& ptime) : ptime_{ptime} {}

            /*!
                \brief Construct a new TimePoint given a string

                \param datetime string representation of date (see example section for acceptable strings)
            */
            explicit TimePoint(const std::string& datetime);

            // TimePoint(const TimePoint& other){ptime_ = other.ptime_; clock_type_ = other.clock_type_;}

            /*!
                \brief Return time point as string, with the given format.

                The format is as described in the following link:

                https://www.boost.org/doc/libs/1_73_0/doc/html/date_time/date_time_io.html

                \param format format string

                \return string representing this time point
            */
            const std::string to_string(const std::string &format="%Y-%m-%dT%H:%M:%S%F") const;

            /*!
                \brief Get the native representation of this date time.

                This shouldn't really be used, as the internals might change to
                not even use a ptime if a limitation is discovered.

                \return the raw, underlying ptime used by this object
            */
            const boost::posix_time::ptime& native() const {return ptime_;} // shouldn't be relied on, but I won't wrap everything

            /*!
                \brief Output time point using the to_string() function with default parameters.

                \param out output stream to which to print
                \param t TimePoint object to print

                \return output stream passed in, for chaining
            */
            inline friend std::ostream& operator<<(std::ostream &out, const TimePoint& t){return out << t.to_string();}

            /*!
                \brief Return duration since given time.

                \param time epoch time to check

                \return duration since given epoch
            */
            // boost::posix_time::time_duration duration_since(const boost::posix_time::ptime& time) const;

            // /*! \brief Return duration since 1970-01-01 00:00:00.
            //     \return duration since 1970-01-01 00:00:00.  */
            // boost::posix_time::time_duration duration_since_1970() const {return duration_since(boost::posix_time::ptime{{1970,1,1}});}
            // /*! \brief Return seconds since 1970-01-01 00:00:00.
            //     \return seconds since 1970-01-01 00:00:00.  */
            // double seconds_since_1970() const {return (double)duration_since_1970().total_microseconds() / 1'000'000.0;}
            // /*! \brief Return duration since 1993-01-01 00:00:00.  
            //     \return duration since 1993-01-01 00:00:00.  */
            // boost::posix_time::time_duration duration_since_1993() const {return duration_since(boost::posix_time::ptime{{1993,1,1}});}
            // /*! \brief Return seconds since 1993-01-01 00:00:00.  
            //     \return seconds since 1993-01-01 00:00:00.  */
            // double seconds_since_1993() const {return (double)duration_since_1993().total_microseconds() / 1'000'000.0;}
            // /*! \brief Return duration since 1958-01-01 00:00:00.  
            //     \return duration since 1958-01-01 00:00:00.  */
            // boost::posix_time::time_duration duration_since_1958() const {return duration_since(boost::posix_time::ptime{{1958,1,1}});}
            // /*! \brief Return seconds since 1958-01-01 00:00:00.  
            //     \return seconds since 1958-01-01 00:00:00.  */
            // double seconds_since_1958() const {return (double)duration_since_1958().total_microseconds() / 1'000'000.0;}

            /*! \brief Return seconds since clock epoch adding leap seconds for TAI based clocks.
                \return seconds since epoch. 
                
                This is the epic for each clock:

                Clock | Epoch
                ----- | -----
                local | Jan 1, 1970, 00:00:00 (standard Unix epoch)
                utc   | Jan 1, 1970, 00:00:00 (standard Unix epoch)
                tai   | Jan 1, 1958, 00:00:00
                tai93 | Jan 1, 1993, 00:00:00
                gps   | Jan 6, 1980, 00:00:00
                tt    | Jan 1, 1977, 00:00:00
                
            */
            double seconds_since(const TimePoint::ClockType clock_type) const;

            /*! \brief Get UTC time (unix) in seconds since Jan 1, 1970.
                \return TAI seconds since 1970 */
            double utc() const { return seconds_since(ClockType::utc); }

            /*! \brief Get TAI time in seconds since Jan 1, 1958 including leap seconds.
                \return TAI seconds since 1958 */
            double tai() const { return seconds_since(ClockType::tai); }

            /*! \brief Get TAI93 time in seconds since Jan 1, 1993 including leap seconds.
                \return TAI seconds since 1993 */
            double tai93() const { return seconds_since(ClockType::tai93); }

            /*! \brief Get GPS time in seconds since Jan 6, 1980 including leap seconds.
                \return GPS seconds since 1980 Jan 6 */
            double gps() const { return seconds_since(ClockType::gps); }

            /*! \brief Get terrestrial time in seconds since Jan 1, 1977.
                \return Terrestrial seconds since 1977 */
            double tt() const { return seconds_since(ClockType::tt); }


            // Developer note for the following Julian functions: This is
            // re-implementing ptime's julian_date(), which returns a long,
            // losing fractional seconds.

            /*! \brief Get Julian date
                \return Julian date  */
            double julian_date() const {return (utc() / 86400.0) + 2440587.5;}
            /*! \brief Get modified Julian date
                \return modified Julian date  */
            double modified_julian_date() const {return julian_date() - 2400000.5;}
            /*! \brief Get reduced Julian date
                \return reduced Julian date  */
            double reduced_julian_date() const {return julian_date() - 2400000.0;}
            /*! \brief Get truncated Julian date
                \return truncated Julian date  */
            double truncated_julian_date() const {return std::floor(julian_date() - 2440000.5);}

            /*!
                \brief Check if two time points are equal

                \param lhs Left hand side argument
                \param rhs Right hand side argument

                \return true if the two given objects are equal
            */
            friend inline bool operator==(const TimePoint& lhs, const TimePoint& rhs){
                return lhs.ptime_ == rhs.ptime_;
            }
            /*!
                \brief Check if two time points are not equal (clock type or time differ)

                \param lhs Left hand side argument
                \param rhs Right hand side argument

                \return true if the two given objects are not equal
            */
            friend inline bool operator!=(const TimePoint& lhs, const TimePoint& rhs){ return !(lhs == rhs); }

            /*!
                \brief Return true if the left TimePoint is less than the right

                \param lhs Left hand side argument
                \param rhs Right hand side argument

                \return true if the left hand side is less
            */
            friend inline bool operator<(const TimePoint& lhs, const TimePoint& rhs){
                return lhs.ptime_ < rhs.ptime_;
            }

            /*!
                \brief Return true if the left TimePoint is greater than the right

                \param lhs Left hand side argument
                \param rhs Right hand side argument

                \return true if the left hand side is greater
            */
            friend inline bool operator> (const TimePoint& lhs, const TimePoint& rhs){ return rhs < lhs; }

            /*!
                \brief Return true if the left TimePoint is less than or equal to the right

                \param lhs Left hand side argument
                \param rhs Right hand side argument

                \return true if the left hand side is less than or equal to the right
            */
            friend inline bool operator<=(const TimePoint& lhs, const TimePoint& rhs){ return !(lhs > rhs); }

            /*!
                \brief Return true if the left TimePoint is greater than or equal to the right

                \param lhs Left hand side argument
                \param rhs Right hand side argument

                \return true if the left hand side is greater than or equal to the right
            */
            friend inline bool operator>=(const TimePoint& lhs, const TimePoint& rhs){ return !(lhs < rhs); }


            // ************************************************************************************
            // Start simple boost::posix_time wrappers

            /*! \brief Get year part of the date
                \return year */
            int16_t year() const {return static_cast<int16_t>(ptime_.date().year());}

            /*! \brief Get month part of the date (one-based, Jan = 1)
                \return month */
            uint8_t month() const {return static_cast<uint8_t>(ptime_.date().month());}

            /*! \brief Get month-based day part of the date (one-based, first of the month = 1)
                \return day */
            uint8_t day() const {return static_cast<uint8_t>(ptime_.date().day());}

            /*! \brief Get day of year part of the date (one-based, Jan 1 = 1)
                \return day of year */
            uint16_t day_of_year() const {return static_cast<uint8_t>(ptime_.date().day_of_year());}

            /*! \brief Get hour part of the time
                \return hour */
            uint8_t hour() const {return static_cast<uint8_t>(ptime_.time_of_day().hours());}

            /*! \brief Get minute part of the time
                \return minute */
            uint8_t minute() const {return static_cast<uint8_t>(ptime_.time_of_day().minutes());}

            /*! \brief Get second part of the time
                \return second */
            uint8_t second() const {return static_cast<uint8_t>(ptime_.time_of_day().seconds());}

            // End simple boost::posix_time wrappers
            // ************************************************************************************
        private:
            boost::posix_time::ptime ptime_;
    };
} // namespace focs

#endif // FOCS_TIMEPOINT

