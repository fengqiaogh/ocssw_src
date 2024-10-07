#ifndef FOCS_LEAPSECONDDATABASE
#define FOCS_LEAPSECONDDATABASE

#include <string>
#include <vector>

#include <boost/date_time/posix_time/posix_time.hpp>

namespace focs {
    /*!
        \class focs::LeapSecond
    
        \brief Class representing the instant in time a leap second was/will be added
    */
    class LeapSecond {
        public:
            /*!
                \brief Construct a with a date and the number of leap seconds in effect starting at the input date

                \param year year
                \param month month
                \param day day
                \param leap_seconds_ Leap seconds starting here
            
            */
            LeapSecond(unsigned short year, unsigned short month, unsigned short day, float leap_seconds_) : time{boost::gregorian::date{year, month, day}}, leap_seconds{leap_seconds_}{}
            boost::posix_time::ptime time; //!< Instance in time of leap second addition
            float leap_seconds; //!< New count of leap seconds
    };
    /*!
        \class focs::LeapSecondDatabase
    
        \brief Reader of leap second databases
    
        \section Example
        \snippet tests/leap_second_database.cpp focs::LeapSecond example
    */
    class LeapSecondDatabase {
        public:
            /*!
                \brief Default constructor, using default locations for database file
            
                The default locations are checked are:
                + The environment variable $LEAPSECOND_DAT
                + $OCVARROOT/viirsn/IETTime.dat
                + $OCVARROOT/modis/leapsec.dat
            */
            LeapSecondDatabase();
            /*!
                \brief Read the database at the given location
            */
            explicit LeapSecondDatabase(const std::string& file){load_file(file);}
            /*!
                \brief Obtain a singleton of a default database
            */
            inline static LeapSecondDatabase& get_default(){
                static LeapSecondDatabase _instance{};
                return _instance;
            }

            /*!
                \brief Look up how many leap seconds are in effect on a given date
            
                \param time The date/time to check
            
                \return the number of leap seconds in effect
            */
            float leap_seconds_since(const boost::posix_time::ptime& time);

            /*!
                \brief Check if a given time is a leap second
            
                \param time The date/time to check
            
                \return true if time is a leap second
            */
            bool is_leap_second(const boost::posix_time::ptime& time);
        private:
            std::string path_;
            std::vector<LeapSecond> leap_seconds_{};

            void load_file(const std::string& file);
    };
    
} // namespace focs

#endif // FOCS_LEAPSECONDDATABASE

