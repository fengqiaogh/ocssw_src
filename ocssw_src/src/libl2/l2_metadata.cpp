#include "l2_metadata.hpp"


     void MetaL2::read(netCDF::NcFile& file) {
        // Read global attributes
        auto attrs = file.getAtts();
        for (auto& [name, attr] : attrs) {
            if (attr.getType() == netCDF::NcType::nc_INT) {
                if (attr.getAttLength() == 1) {
                    int value;
                    attr.getValues(&value);
                    attributes[attr.getName()] = value;
                } else {
                    std::vector<int> value(attr.getAttLength());
                    attr.getValues(value.data());
                    attributes[attr.getName()] = value;
                }
            } else if (attr.getType() == netCDF::NcType::nc_FLOAT) {
                if (attr.getAttLength() == 1) {
                    float value;
                    attr.getValues(&value);
                    attributes[attr.getName()] = value;
                } else {
                    std::vector<float> value(attr.getAttLength());
                    attr.getValues(value.data());
                    attributes[attr.getName()] = value;
                }
            } else if (attr.getType() == netCDF::NcType::nc_STRING) {
                std::string value;
                attr.getValues(value);
                attributes[attr.getName()] = value;
            } else if (attr.getType() == netCDF::NcType::nc_DOUBLE) {
                if (attr.getAttLength() == 1) {
                    double value;
                    attr.getValues(&value);
                    attributes[attr.getName()] = value;
                } else {
                    std::vector<double> value(attr.getAttLength());
                    attr.getValues(value.data());
                    attributes[attr.getName()] = value;
                }
            } else if (attr.getType() == netCDF::NcType::nc_SHORT) {
                if (attr.getAttLength() == 1) {
                    short value;
                    attr.getValues(&value);
                    attributes[attr.getName()] = value;
                } else {
                    std::vector<short> value(attr.getAttLength());
                    attr.getValues(value.data());
                    attributes[attr.getName()] = value;
                }
            } else if (attr.getType() == netCDF::NcType::nc_INT64) {
                if (attr.getAttLength() == 1) {
                    long long value;
                    attr.getValues(&value);
                    attributes[attr.getName()] = value;
                } else {
                    std::vector<long long> value(attr.getAttLength());
                    attr.getValues(value.data());
                    attributes[attr.getName()] = value;
                }
            } else if (attr.getType() == netCDF::NcType::nc_CHAR) {
                if (attr.getAttLength() == 1) {
                    char value;
                    attr.getValues(&value);
                    attributes[attr.getName()] = value;
                } else {
                    std::string value;
                    attr.getValues(value);
                    attributes[attr.getName()] = value;
                }
            }
        }
        // read group specific attributes and variables
        netCDF::NcGroup scan_time_group = file.getGroup("scan_line_attributes");
        if (!scan_time_group.isNull()) {
            netCDF::NcVar scan_time_var = scan_time_group.getVar("scan_time");
            if (scan_time_var.isNull()) {
                scan_time_var = scan_time_group.getVar("time");
            }
            if (scan_time_var.isNull()) {
                scan_time_var = scan_time_group.getVar("scantime");
            }
            if (!scan_time_var.isNull()) {
                std::vector<double> scan_time_values(scan_time_var.getDim(0).getSize());
                scan_time_var.getVar(scan_time_values.data());
                attributes["scan_time"] = scan_time_values;
            }
            netCDF::NcVar year_var = scan_time_group.getVar("year");
            if (!year_var.isNull()) {
                std::vector<int> year_values(year_var.getDim(0).getSize());
                year_var.getVar(year_values.data());
                attributes["year"] = year_values;
            }
            netCDF::NcVar day_var = scan_time_group.getVar("day");
            if (!day_var.isNull()) {
                std::vector<int> day_values(day_var.getDim(0).getSize());
                day_var.getVar(day_values.data());
                attributes["day"] = day_values;
            }
            netCDF::NcVar msec_var = scan_time_group.getVar("msec");
            if (!msec_var.isNull()) {
                std::vector<int> msec_values(msec_var.getDim(0).getSize());
                msec_var.getVar(msec_values.data());
                attributes["msec"] = msec_values;
            }
        }
        netCDF::NcGroup processing_control = file.getGroup("processing_control");
        if (!processing_control.isNull()) {
            netCDF::NcGroupAtt sourc_nc = processing_control.getAtt("source");
            if (!sourc_nc.isNull()) {
                std::string source;
                sourc_nc.getValues(source);
                attributes["source"] = source;
            }
        }
    }

    MetaL2::MetaL2(netCDF::NcFile& file) {
        read(file);
        if (get_attribute<std::string>("time_coverage_start").has_value() &&
            get_attribute<std::string>("time_coverage_end").has_value()) {
            time_coverage_start = get_attribute<std::string>("time_coverage_start").value();
            time_coverage_end = get_attribute<std::string>("time_coverage_end").value();
            unix2yds(isodate2unix(time_coverage_start.c_str()), &syear, &sday, &smsec);
            unix2yds(isodate2unix(time_coverage_end.c_str()), &eyear, &eday, &emsec);
            smsec = std::roundl(smsec * 1000);
            emsec = std::roundl(emsec * 1000);
        }
        if (get_attribute<float>("westernmost_longitude").has_value())
            westlon = get_attribute<float>("westernmost_longitude").value();
        if (get_attribute<float>("easternmost_longitude").has_value())
            eastlon = get_attribute<float>("easternmost_longitude").value();
        if (get_attribute<float>("southernmost_longitude").has_value())
            southlat = get_attribute<float>("southernmost_latitude").value();
        if (get_attribute<float>("northernmost_longitude").has_value())
            northlat = get_attribute<float>("northernmost_latitude").value();
        if (get_attribute<std::string>("source").has_value())
            source = get_attribute<std::string>("source").value();
        if (get_attribute<std::string>("instrument").has_value())
            sensor_name = get_attribute<std::string>("instrument").value();
        if (get_attribute<std::string>("platform").has_value())
            mission = get_attribute<std::string>("platform").value();
        if (get_attribute<std::string>("flag_meanings").has_value())
            flag_names = get_attribute<std::string>("flag_meanings").value();
        if (get_attribute<std::string>("title").has_value())
            title = get_attribute<std::string>("title").value();
        if (get_attribute<std::vector<int>>("flag_masks").has_value())
            bits = get_attribute<std::vector<int>>("flag_masks").value();
        if (get_attribute<std::vector<double>>("scan_time").has_value())
            scan_time = get_attribute<std::vector<double>>("scan_time").value();
        if (get_attribute<std::vector<int>>("year").has_value())
            year = get_attribute<std::vector<int>>("year").value();
        if (get_attribute<std::vector<int>>("day").has_value())
            day = get_attribute<std::vector<int>>("day").value();
        if (get_attribute<std::vector<int>>("msec").has_value())
            msec = get_attribute<std::vector<int>>("msec").value();
    }