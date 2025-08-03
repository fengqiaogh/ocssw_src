#include "get_geospatial.hpp"
#include "find_variable.hpp"
#include "boost/algorithm/string.hpp"
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <algorithm>
#include <stack>
#include "scaled_nc_var.hpp"
namespace bg = boost::geometry;
typedef bg::model::point<double, 2, bg::cs::spherical_equatorial<bg::degree>> Point_t;
typedef bg::model::linestring<Point_t> Linestring_t;
typedef bg::model::polygon<Point_t> Polygon_t;
typedef bg::model::multi_point<Point_t> MultiPoint_t;

Geospatialbounds::Geospatialbounds() {
    scan_line_vars["elat"] = std::vector<float>();
    scan_line_vars["slat"] = std::vector<float>();
    scan_line_vars["clat"] = std::vector<float>();
    scan_line_vars["elon"] = std::vector<float>();
    scan_line_vars["slon"] = std::vector<float>();
}

std::string Geospatialbounds::get_day_night_flag(const netCDF::NcFile &nc_file) {
    ScaledNcVar solzen_nc;
    for (const auto &name : possible_sozlen_names) {
        solzen_nc = find_nc_variable_cpp<const netCDF::NcFile, ScaledNcVar>(name, nc_file);
        if (!solzen_nc.isNull())
            break;
    }
    if (solzen_nc.isNull()) {
        return "Undefined";
    }
    std::vector<float> solzen;
    allocate_var(solzen, solzen_nc);
    const float night_angle = 90.0f;
    std::stack<bool> state;
    for (const auto &angle : solzen) {
        if (angle > 180.0 || angle < -180 || std::isnan(angle) || !std::isfinite(angle))
            continue;
        bool flag = angle < night_angle;
        if (state.empty()) {
            state.push(flag);
        } else if (state.top() != flag) {
            return "Mixed";
        }
    }
    // release memory
    free_vector(solzen);
    if (state.empty()) {
        return "Undefined";
    } else {
        if (state.top())
            return "Day";
        else
            return "Night";
    }
}

void Geospatialbounds::get_lat_lon(const netCDF::NcFile &nc_input) {
    netCDF::NcVar latnc;
    netCDF::NcVar lonnc;
    find_lat_lon<const netCDF::NcFile,netCDF::NcVar>(nc_input, latnc, lonnc);
    allocate_var(lat, latnc);
    allocate_var(lon, lonnc);
    std::vector<netCDF::NcDim> dims = latnc.getDims();
    number_of_lines = dims.at(0).getSize();
    pixels_per_line = dims.at(1).getSize();
}

void Geospatialbounds::allocate_var(std::vector<float> &data, const netCDF::NcVar &var) {
    if (!var.isNull()) {
        std::vector<netCDF::NcDim> dims = var.getDims();
        size_t size = 1;
        for (const auto &dim : dims) {
            size *= dim.getSize();
        }
        data.resize(size);
        var.getVar(data.data());
    }
}

void Geospatialbounds::allocate_attr(std::vector<float> &data, const netCDF::NcAtt &att) {
    if (!att.isNull()) {
        size_t size = att.getAttLength();
        if (size == 0)
            return;
        if (att.getType() == NC_STRING || att.getType() == NC_CHAR)
            return;
        data.resize(size);
        att.getValues(data.data());
    }
}

Geospatialbounds::Geospatialbounds(const std::string &path)
    : Geospatialbounds(netCDF::NcFile(path, netCDF::NcFile::read)) {
}

Geospatialbounds::Geospatialbounds(const netCDF::NcFile &nc_input) : Geospatialbounds() {
    netCDF::NcGroupAtt attr = nc_input.getAtt("instrument");
    if (attr.isNull()) {
        std::cerr << "Warning: Instrument not found. Returning an empty string" << std::endl;
    } else {
        attr.getValues(instrument);
    }
    attr = nc_input.getAtt("platform");
    if (attr.isNull()) {
        std::cerr << "Warning: Platform not found. Returning an empty string" << std::endl;
    } else {
        attr.getValues(platform);
    }
    // search for gring lat/lon or polygon in global attributes

    const auto attrs = nc_input.getAtts();
    for (const auto &_attr : attrs) {
        std::string attr_name = _attr.first;
        boost::algorithm::to_lower(attr_name);
        if (attr_name == "gringpointlatitude") {
            allocate_attr(gringpointlatitude, _attr.second);
        }
        if (attr_name == "gringpointlongitude") {
            allocate_attr(gringpointlongitude, _attr.second);
        }
        if (attr_name == "geospatial_bounds") {
            _attr.second.getValues(geospatial_bounds_wkt);
        }
    }
    // search for group attributes
    if (gringpointlatitude.empty()) {
        auto gringpointlatitude_attr = find_nc_grp_attribute_cpp("gringpointlatitude", nc_input);
        allocate_attr(gringpointlatitude, gringpointlatitude_attr);
    }
    if (gringpointlongitude.empty()) {
        auto gringpointlongitude_attr = find_nc_grp_attribute_cpp("gringpointlongitude", nc_input);
        allocate_attr(gringpointlongitude, gringpointlongitude_attr);
    }
    if (geospatial_bounds_wkt.empty()) {
        auto geospatial_bounds_wkt_attr = find_nc_grp_attribute_cpp("geospatial_bounds", nc_input);
        if (!geospatial_bounds_wkt_attr.isNull())
            geospatial_bounds_wkt_attr.getValues(geospatial_bounds_wkt);
    }
    auto crs = find_nc_grp_attribute_cpp("geospatial_bounds_crs", nc_input);
    if (!crs.isNull()) {
        Polygon_t poly;

        // get rid of 0 at end of string
        if(geospatial_bounds_wkt.back() == 0)
            geospatial_bounds_wkt.pop_back();

        boost::geometry::read_wkt(geospatial_bounds_wkt, poly);
        geospatial_bounds_wkt = "POLYGON((";
        double first_ini = true;
        double x_first, y_first;
        double x, y;
        // getting the vertices back
        for (auto it = boost::begin(boost::geometry::exterior_ring(poly));
             it != boost::end(boost::geometry::exterior_ring(poly)); ++it) {
            x = bg::get<0>(*it);
            y = bg::get<1>(*it);
            if (first_ini) {
                x_first = x;
                y_first = y;
                first_ini = false;
            }
            geospatial_bounds_wkt += std::to_string(y);
            geospatial_bounds_wkt += " ";
            geospatial_bounds_wkt += std::to_string(x);
            geospatial_bounds_wkt += ", ";
            // use the coordinates...
        }
        // Polygon must be closed though it does not affect dateline crossing
        if(x != x_first || y != y_first) {
            geospatial_bounds_wkt += std::to_string(y_first);
            geospatial_bounds_wkt += " ";
            geospatial_bounds_wkt += std::to_string(x_first);
        } else {
            geospatial_bounds_wkt.pop_back();
            geospatial_bounds_wkt.pop_back();
        }
        geospatial_bounds_wkt += "))";
    }
    if (std::isnan(max_lon)) {
        auto geospatial_lon_max_attr = find_nc_grp_attribute_cpp("geospatial_lon_max", nc_input);
        if (!geospatial_lon_max_attr.isNull())
            geospatial_lon_max_attr.getValues(&max_lon);
    }
    if (std::isnan(min_lon)) {
        auto geospatial_lon_min_attr = find_nc_grp_attribute_cpp("geospatial_lon_min", nc_input);
        if (!geospatial_lon_min_attr.isNull())
            geospatial_lon_min_attr.getValues(&min_lon);
    }
    if (time_coverage_start.empty()) {
        auto attr_time_coverage = find_nc_grp_attribute_cpp("time_coverage_start", nc_input);
        if (!attr_time_coverage.isNull())
            attr_time_coverage.getValues(time_coverage_start);
    }
    if (time_coverage_end.empty()) {
        auto attr_time_coverage = find_nc_grp_attribute_cpp("time_coverage_end", nc_input);
        if (!attr_time_coverage.isNull())
            attr_time_coverage.getValues(time_coverage_end);
    }
    if (day_night_flag.empty()) {
        auto attr_day_night_flag = find_nc_grp_attribute_cpp("day_night_flag", nc_input);
        if (!attr_day_night_flag.isNull())
            attr_day_night_flag.getValues(day_night_flag);
        else
            day_night_flag = get_day_night_flag(nc_input);
    }
    for (auto &sc_ln_attr : scan_line_vars) {
        std::string var_name = sc_ln_attr.first;
        auto &var_val = sc_ln_attr.second;
        netCDF::NcVar var = find_nc_variable_cpp(var_name, nc_input);
        allocate_var(var_val, var);
    }
    for (const auto &sc_ln_attr : scan_line_vars) {
        auto &var_val = sc_ln_attr.second;
        if (var_val.empty()) {
            get_lat_lon(nc_input);
            break;
        }
        number_of_lines = var_val.size();
    }
}

void Geospatialbounds::calc_scan_line(const std::vector<float> &lat, const std::vector<float> &lon,
                                      size_t numScans, size_t numPixels,
                                      std::unordered_map<std::string, std::vector<float>> &scan_line_vars) {
    assert(lat.size() == numPixels * numScans);
    assert(lon.size() == numPixels * numScans);
    scan_line_vars.at("elat") = std::vector<float>(numScans, NAN);
    scan_line_vars.at("elon") = std::vector<float>(numScans, NAN);
    scan_line_vars.at("slon") = std::vector<float>(numScans, NAN);
    scan_line_vars.at("slat") = std::vector<float>(numScans, NAN);
    scan_line_vars.at("clat") = std::vector<float>(numScans, NAN);
    std::vector<int8_t> mask(lat.size(), 1);
    auto f_lat = [](const float &v) {
        if (v > 90.0 || v < -90 || std::isnan(v) || !std::isfinite(v))
            return 0;
        else
            return 1;
    };
    auto f_lon = [](const float &v, const int8_t &mask) {
        if (v > 180.0 || v < -180 || std::isnan(v) || !std::isfinite(v) || !mask)
            return 0;
        else
            return 1;
    };
    std::transform(lat.begin(), lat.end(), mask.begin(), f_lat);
    std::transform(lat.begin(), lat.end(), mask.begin(), mask.begin(), f_lon);
    auto get_index = [&](size_t is, size_t ip) { return numPixels * is + ip; };
    float last_lat = NAN;
    for (size_t is = 0; is < numScans; is++) {
        size_t spix = 0;
        size_t cpix = numPixels / 2;
        size_t epix = numPixels - 1;
        // find good geolocation for spix and epix
        size_t good_spix = spix;
        size_t good_epix = epix;

        for (size_t ip = spix; ip < epix; ip++) {
            size_t i = get_index(is, ip);
            if (mask.at(i)) {
                good_spix = ip;
                break;
            }
        }
        for (size_t ip = epix; ip >= spix; ip--) {
            size_t i = get_index(is, ip);
            if (mask.at(i)) {
                good_epix = ip;
                break;
            }
        }
        size_t i = get_index(is, good_epix);
        scan_line_vars.at("elat")[is] = lat.at(i);
        scan_line_vars.at("elon")[is] = lon.at(i);
        i = get_index(is, good_spix);
        scan_line_vars.at("slat")[is] = lat.at(i);
        scan_line_vars.at("slon")[is] = lon.at(i);
        if (std::isnan(last_lat))
            last_lat = (scan_line_vars.at("elat")[is] + scan_line_vars.at("slat")[is]) / 2;
        i = get_index(is, cpix);
        if (f_lat(lat.at(i))) {
            last_lat = lat.at(i);
        }
        scan_line_vars.at("clat")[is] = last_lat;
    }
};

std::pair<std::vector<float>, std::vector<float>> Geospatialbounds::calc_gring(
    const std::vector<float> &slat, const std::vector<float> &elat, const std::vector<float> &clat,
    const std::vector<float> &slon, const std::vector<float> &elon) {
    const float geobox_bounds{20.0f};
    std::vector<std::vector<float>> geobox{{}, {}, {}, {}};
    // find last valid scan
    size_t nscan = clat.size();
    // find first valid scan
    size_t first_scan = nscan;
    size_t last_scan = nscan;
    for (size_t is = 0; is < nscan; is++) {
        if (!std::isnan(elat.at(is))) {
            if (first_scan == nscan)
                first_scan = is;
            last_scan = is;
        }
    }
    if (first_scan == nscan) {
        std::cerr << "-Error- : not valid elat is supplied. Exiting";
        exit(EXIT_FAILURE);
    }
    geobox.at(0).push_back(slon.at(first_scan));
    geobox.at(1).push_back(slat.at(first_scan));
    geobox.at(2).push_back(elon.at(first_scan));
    geobox.at(3).push_back(elat.at(first_scan));
    //
    float last_lat = clat.at(first_scan);
    for (size_t is = first_scan + 1; is < last_scan; is++) {
        if (std::isnan(clat.at(is)))
            continue;
        if (std::isnan(elat.at(is)))
            continue;
        if (std::abs(last_lat - clat.at(is)) > geobox_bounds) {
            geobox.at(0).push_back(slon.at(is));
            geobox.at(1).push_back(slat.at(is));
            geobox.at(2).push_back(elon.at(is));
            geobox.at(3).push_back(elat.at(is));
            last_lat = clat.at(is);
        }
    }
    geobox.at(0).push_back(slon.at(last_scan));
    geobox.at(1).push_back(slat.at(last_scan));
    geobox.at(2).push_back(elon.at(last_scan));
    geobox.at(3).push_back(elat.at(last_scan));
    std::vector<float> lat_gring;
    std::vector<float> lon_gring;
    lat_gring.push_back(geobox.at(1).at(0));
    for (const auto &val : geobox.at(3)) {
        lat_gring.push_back(val);
    }
    while (geobox.at(1).size() > 1) {
        lat_gring.push_back(geobox.at(1).back());
        geobox.at(1).pop_back();
    }

    lon_gring.push_back(geobox.at(0).at(0));
    for (const auto &val : geobox.at(2)) {
        lon_gring.push_back(val);
    }
    while (geobox.at(0).size() > 1) {
        lon_gring.push_back(geobox.at(0).back());
        geobox.at(0).pop_back();
    }
    return {lat_gring, lon_gring};
}

std::string Geospatialbounds::calc_gring(const std::pair<std::vector<float>, std::vector<float>> &gring) {
    const std::vector<float> &lat_gring = gring.first;
    const std::vector<float> &lon_gring = gring.second;
    std::string wkt = "POLYGON((";
    size_t len = lat_gring.size();
    for (size_t i = 0; i < len; i++) {
        wkt += std::to_string(lon_gring.at(i));
        wkt += " ";
        wkt += std::to_string(lat_gring.at(i));
        wkt += ", ";
    }
    // Polygon must be closed though it does not affect dateline crossing
    wkt += std::to_string(lon_gring.at(0));
    wkt += " ";
    wkt += std::to_string(lat_gring.at(0));
    wkt += "))";
    return wkt;
};

const std::unordered_map<std::string, std::vector<float>> &Geospatialbounds::get_bounds() {
    if (scan_line_vars.at("elat").empty()) {
        calc_scan_line(lat, lon, number_of_lines, pixels_per_line, scan_line_vars);
        // release memory
        free_vector(lat);
        free_vector(lon);
    }
    return scan_line_vars;
}

std::pair<std::vector<float>, std::vector<float>> Geospatialbounds::get_gring() {
    if (gringpointlatitude.empty()) {
        scan_line_vars = get_bounds();
        const auto &slat = scan_line_vars.at("slat");
        const auto &elat = scan_line_vars.at("elat");
        const auto &clat = scan_line_vars.at("clat");
        const auto &slon = scan_line_vars.at("slon");
        const auto &elon = scan_line_vars.at("elon");
        auto res = calc_gring(slat, elat, clat, slon, elon);
        gringpointlatitude = res.first;
        gringpointlongitude = res.second;
    }
    return {gringpointlatitude, gringpointlongitude};
}

std::string Geospatialbounds::get_bounds_wkt() {
    if (geospatial_bounds_wkt.empty()) {
        geospatial_bounds_wkt = calc_gring(get_gring());
    }
    return geospatial_bounds_wkt;
}

float Geospatialbounds::get_max_lon() {
    if (std::isnan(max_lon)) {
        auto res = get_gring();
        const auto &longring = res.second;
        max_lon = *std::max_element(longring.begin(), longring.end());
        min_lon = *std::min_element(longring.begin(), longring.end());
    }
    return max_lon;
}

float Geospatialbounds::get_min_lon() {
    get_max_lon();
    return min_lon;
}

const float *Geospatialbounds::get_elat() {
    if (scan_line_vars.at("elat").empty())
        scan_line_vars = get_bounds();
    return scan_line_vars.at("elat").data();
}

const float *Geospatialbounds::get_slat() {
    if (scan_line_vars.at("slat").empty())
        scan_line_vars = get_bounds();
    return scan_line_vars.at("slat").data();
}