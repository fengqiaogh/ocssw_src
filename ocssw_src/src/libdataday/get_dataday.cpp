/*
 * Given a L2 file, determine the dataday(s) it contains
 */

#include "get_dataday.h"
#include "get_dataday.hpp"
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/optional.hpp>
#include <get_geospatial.hpp>
#include <cstdlib>
#include <iostream>
#include <netcdf>

#include "sensorDefs.h"
#include "sensorInfo.h"
#include "timeutils.h"
#include "version.h"
#define _PRINT_(...) if(ddvals::verbosity) printf(__VA_ARGS__)

namespace ddvals {
boost::optional<float> equatorialCrossingTime;
int32_t plus_day{0};
float deltaeqcross{0.0};
bool verbosity = false;
}  // namespace ddvals
using namespace std;
using namespace netCDF;
using namespace netCDF::exceptions;
namespace bg = boost::geometry;

typedef bg::model::point<double, 2, bg::cs::spherical_equatorial<bg::degree>> Point_t;
typedef bg::model::linestring<Point_t> Linestring_t;
typedef bg::model::polygon<Point_t> Polygon_t;
typedef bg::model::multi_point<Point_t> MultiPoint_t;

//void print(std::ostream& stream, const char* format) {
//    if (!ddvals::verbosity)
//        return;
//    stream << format;
//}
// h
/**
 * @brief
 * ttps://www.vertica.com/docs/11.0.x/HTML/Content/Authoring/AnalyzingData/Geospatial/Spatial_Definitions/WellknownTextWKT.htm
 * @param wkt_string
 */
Polygon_t parse_wkr_polygon(const std::string& wkt_string) {
    Polygon_t poly;
    boost::geometry::read_wkt(wkt_string, poly);
    return poly;
}

//template <typename T, typename... Targs>
//void print(std::ostream& stream, const char* format, T&& value, Targs&&... fargs) {
//    if (!ddvals::verbosity)
//        return;
//    for (; *format != '\0'; format++) {
//        if (*format == '%') {
//            stream << value;
//            print(stream, format + 1, fargs...);
//            return;
//        }
//        stream << *format;
//    }
//}
//
//void print(const char* format) {
//    if (!ddvals::verbosity)
//        return;
//    cout << format;
//}
//
//template <typename T, typename... Targs>
//void print(const char* format, T&& value, Targs&&... fargs) {
//    print(std::cout, format, value, fargs...);
//}

enum class DL { NOT_CROSSED, CROSSED, CROSSED_NORTH_POLE, CROSSED_SOUTH_POLE, CROSSED_IRREGULAR };

Polygon_t gRing_to_gPolygon(float* gRingLats, float* gRingLons, size_t length) {
    Polygon_t gPolygon;

    for (size_t i = 0; i < length; i++) {
        double lat = boost::lexical_cast<double>(gRingLats[i]);
        double lon = boost::lexical_cast<double>(gRingLons[i]);
        gPolygon.outer().push_back(Point_t{lon, lat});
    }
    gPolygon.outer().push_back(
        Point_t{boost::lexical_cast<double>(gRingLons[0]), boost::lexical_cast<double>(gRingLats[0])});

    return gPolygon;
}

void printUsage(int32_t exitStatus) {
    std::string softwareVersion;
    softwareVersion += std::to_string(VERSION_MAJOR);
    softwareVersion += ".";
    softwareVersion += std::to_string(VERSION_MINOR);
    softwareVersion += ".";
    softwareVersion += std::to_string(VERSION_PATCH);
    softwareVersion += "-";
    softwareVersion += GITSHA;

    cout << "get_dataday " << softwareVersion << endl;
    cout << "\nUsage: get_dataday ifile" << endl;
    cout << "where:\n\tifile is an L2 file generated by l2gen\n" << endl;
    cout << "Options:\n\t-h, --help: get this clever little usage statement" << endl;
    exit(exitStatus);
}

void get_datadays(time_t starttime,             /* (in)  swath start time
                                           (seconds since 1-Jan-1970 00:00:00 GMT) */
                  float equatorialCrossingTime, /* (in)  sensor's nominal local equator
                                        crossing time (expressed in hours) */
                  DL dateLineCrossed,           /* (in)  indicates whether swath crosses dateline */
                  float west,                   /* (in)  westernmost longitude of swath */
                  float east,                   /* (in)  easternmost longitude of swath */
                  int32_t* dataday0,            /* (out) dataday of swath (days since 1-Jan-1970) */
                  int32_t* dataday1             /* (out) 2nd dataday for datline-spanning swaths */
) {
    time_t referenceTime;
    int referenceDay;
    float referenceHour;

    referenceTime = (time_t)(starttime + (12 - (double)equatorialCrossingTime) * 3600);
    referenceDay = referenceTime / 86400;
    referenceHour = (referenceTime % 86400) / 3600.0;

    if (dateLineCrossed == DL::NOT_CROSSED) {
        *dataday1 = *dataday0 = referenceDay;
    } else if (referenceHour < 6) {
        *dataday0 = referenceDay - 1;
        *dataday1 = referenceDay;
    } else if (referenceHour > 18) {
        *dataday0 = referenceDay;
        *dataday1 = referenceDay + 1;
    } else if (dateLineCrossed == DL::CROSSED_NORTH_POLE || dateLineCrossed == DL::CROSSED_SOUTH_POLE) {
        *dataday1 = *dataday0 = referenceDay;
    } else {
        float westOfDateline = 180 - west; /* number of degrees west of dateline */
        float eastOfDateline = east + 180; /* number of degrees east of dateline */

        if (westOfDateline > eastOfDateline) {
            *dataday0 = referenceDay - 1;
            *dataday1 = referenceDay;
        } else {
            *dataday0 = referenceDay;
            *dataday1 = referenceDay + 1;
        }
    }
}

void getEquatorCrossingTime(int32_t sensorID, bool dayNight, time_t starttime, float* equatorialCrossingTime,
                            int32_t* plusDay) {
    if(std::abs(ddvals::deltaeqcross) > 1e-5)
        *equatorialCrossingTime = 12.0 + ddvals::deltaeqcross/ 60.0;
    else
    {
        switch (sensorID) {
            case MERIS:
            case OLCIS3A:
            case OLCIS3B:
                *equatorialCrossingTime = 10.0;
                break;
            case MODIST:
            case OCTS:
                *equatorialCrossingTime = 10.5;
                break;
            case HAWKEYE:
            case CZCS:
            case HICO:
                *equatorialCrossingTime = 12.0;
                break;
            case OCI:
            case SPEXONE:
            case HARP2:
                *equatorialCrossingTime = 13.0;
                break;
            case MODISA:
            case VIIRSN:
            case VIIRSJ1:
            case VIIRSJ2:
                *equatorialCrossingTime = 13.5;
                break;
            case SEAWIFS:
                *equatorialCrossingTime = 12.0;
                int16_t year;
                int16_t day;
                double secs;
                unix2yds(starttime, &year, &day, &secs);
                if (year > 2002) {
                    // The constant 10957 makes d the number of days since 1 Jan.
                    // 2000
                    int32_t d = starttime / 86400.0 - 10957.0;
                    /*
                     * On 10 July 2010 (d=3843) OrbView-2/SeaWiFS was nearing the
                     * end of its orbit-raising maneuvers.  Before these maneuvers
                     * the node-crossing time was drifting further into the
                     afternoon
                     * (first equation below).  After the orbit raising, the
                     node-crossing
                     * time was drifting back towards noon (second equation).

                     * Correction equations provided by Fred Patt.
                     */
                    double deg;
                    if (d < 3843) {
                        deg = 7.7517951e-10 * d * d * d - 2.1692192e-06 * d * d + 0.0070669241 * d - 4.1300585;
                    } else {
                        /*
                         * This equation may need to be replaced with a polynomial
                         * once we get more orbit data under our belts. (16-Jul-2010
                         * N.Kuring)
                         * ...note...the above comment was written before the demise
                         * of SeaWiFS on 11-Dec-2010 ...no further modifications
                         * necessary...
                         */
                        deg = -0.024285181 * d + 128.86093;
                    }

                    // The above polynomials yield degrees; convert to hours.
                    float hours = (float)(deg / 15.);

                    *equatorialCrossingTime = 12.0 + hours;
                }

                break;
            default:
                cerr << "-W- Unknown equator crossing time for sensorID=" << sensorID << endl;
                cerr << "-W- Assuming local noon crossing" << endl;
                *equatorialCrossingTime = 12.0;
                break;
        }
    }
    // Shift the equatorial crossing time by 12 hours if doing night time
    // binning This is done so that if crossing the dateline the reference hour
    // is being used to move into a dataday in the same direction in function
    // get_datadays
    if (dayNight) {
        *equatorialCrossingTime = *equatorialCrossingTime - 12;
        if (*equatorialCrossingTime < 0) {
            *equatorialCrossingTime = *equatorialCrossingTime + 24;
            *plusDay = 1;  // if we do this, we're effectively going back
            // a day, and we so we need to add a day to the
            // output of get_datadays...
        }
    }
    // setting equatorialCrossingTime
    ddvals::equatorialCrossingTime = *equatorialCrossingTime;
    ddvals::plus_day = *plusDay;
}

std::string dataday_to_isodate(int32_t dataday) {
    double datatime = dataday * 86400.0;
    std::string isoDate(unix2isodate(datatime, 'G'));
    return isoDate;
}

int32_t get_plusday()
{
    return ddvals::plus_day;
}

float get_equatorial_crossingTime() {
    if (ddvals::equatorialCrossingTime == boost::none) {
        std::cerr << "Equatorial Crossing Time is not set. Exiting ...";
        exit(EXIT_FAILURE);
    } else {
        return ddvals::equatorialCrossingTime.get();
    }
}

std::pair<int32_t,int32_t> get_datadays(const std::string &filepath)
{
        std::pair<int32_t,int32_t> res;
        try {
            NcFile nc_input(filepath, NcFile::read);
            res =  get_datadays(nc_input);
            nc_input.close();

        } catch (NcException& e) {
            e.what();
            cerr << "\nFailure opening  input file: " + filepath + "\n" << endl;
            printUsage(NC2_ERR);
        }
        return res;
}

std::pair<int32_t, int32_t> get_datadays(const NcFile& nc_input, float deltaeqcross, int night_flag) {
    // double const earth_radius = 6371.229;

    // Define the dateline, north and south poles
    Linestring_t dateline;
    bg::append(dateline, Point_t(180.0, 89.9999));
    bg::append(dateline, Point_t(180.0, -89.9999));
    Point_t northPole = {0, 90.0};
    Point_t southPole = {0, -90.0};
    std::string instrument;
    std::string platform;
    std::string day_night_flag;
    float maxEast;
    float maxWest;
    std::string sceneStartTimeISO;
    Polygon_t gPolygon;
    int32_t sensorID;
    double sceneStartTime;
    // starting reading the file
    

        Geospatialbounds geo_bounds(nc_input);
        maxEast = geo_bounds.get_max_lon();
        maxWest = geo_bounds.get_min_lon();
        sceneStartTimeISO = geo_bounds.get_time_coverage_start();
        platform = geo_bounds.get_platform();
        instrument = geo_bounds.get_instrument();
        if (platform.empty())
            std::cerr << "Warning: platform is not found" << std::endl;
        if (instrument.empty())
            std::cerr << "Warning: instrument is not found" << std::endl;
        if (sceneStartTimeISO.empty()) {
            std::cerr << "Error: start time is not found. Exiting" << std::endl;
            exit(EXIT_FAILURE);
        }
        sensorID = instrumentPlatform2SensorId(instrument.c_str(), platform.c_str());
        sceneStartTime = isodate2unix(sceneStartTimeISO.c_str());
        std::string geospatial_bounds = geo_bounds.get_bounds_wkt();
        gPolygon = parse_wkr_polygon(geospatial_bounds);
        day_night_flag = geo_bounds.get_day_night_flag();
    ddvals::deltaeqcross = deltaeqcross;
    
    /// ending reading the file
    DL dateLineCrossed = DL::NOT_CROSSED;
    std::vector<Point_t> result;
    bg::intersection(dateline, gPolygon, result);
    if (result.size()) {
        if (bg::within(northPole, gPolygon)) {
            dateLineCrossed = DL::CROSSED_NORTH_POLE;
            _PRINT_("# Scene crossed North pole\n");
        } else if (bg::within(southPole, gPolygon)) {
            _PRINT_("# Scene crossed South pole\n");
            dateLineCrossed = DL::CROSSED_SOUTH_POLE;
        } else {
            _PRINT_("# Scene crossed dateline\n");
            dateLineCrossed = DL::CROSSED;
        }
    }

    int32_t dataday0, dataday1;
    int32_t plusDay = 0;
    float equatorialCrossingTime;
    bool nightScene = night_flag;
    if (day_night_flag.c_str() == std::string("Night")) {
        nightScene = true;
    }
    getEquatorCrossingTime(sensorID, nightScene, sceneStartTime, &equatorialCrossingTime, &plusDay);

    get_datadays(sceneStartTime, equatorialCrossingTime, dateLineCrossed, maxWest, maxEast, &dataday0,
                 &dataday1);

    std::string startDataDay = dataday_to_isodate(dataday0 + plusDay);
    std::string stopDataDay = dataday_to_isodate(dataday1 + plusDay);
    // Print out day/night flag
    _PRINT_("Day_Night_Flag=%s\n", day_night_flag.c_str());
    // Print out the derived dataday(s)
    if (dataday0 == dataday1) {
        _PRINT_("DataDay0=%s\n", startDataDay.substr(0, 10).c_str());

    } else {
        _PRINT_("DataDay0=%s\nDataDay1=%s\n", startDataDay.substr(0, 10).c_str(), stopDataDay.substr(0, 10).c_str());
    }
    return {dataday0 + plusDay, dataday1 + plusDay};
}

int get_datadays(const char* path, int32_t* day0, int32_t* day1) {
    const auto data = get_datadays(path);
    *day0 = data.first;
    *day1 = data.second;
    return EXIT_SUCCESS;
}

void set_verbosity(int val) {
    if (val != 0) {
        ddvals::verbosity = true;
    }
}
