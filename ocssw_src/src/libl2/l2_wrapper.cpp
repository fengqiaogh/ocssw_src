#include <stdlib.h>
#include <stdio.h>
#include <netcdf.h>
#include <genutils.h>
#include "l2_wrapper.h"
#include <sensorInfo.h>
#include <timeutils.h>
#include <nc4utils.h>
#include <hdf.h>
#include "l2_reader.hpp"

namespace {
// Array of L2_Reader objects to handle multiple L2 files (up to MAXNFILES)
std::array<L2_Reader, MAXNFILES> l2_reader;

// Arrays to store L2 data for each file
std::array<std::vector<float *>, MAXNFILES> l2_data;  // Product data values

// Arrays to store temporal information
std::array<std::vector<double>, MAXNFILES> l2_scan_time;  // Scan timestamps
std::array<std::vector<int>, MAXNFILES> l2_year;          // Year values
std::array<std::vector<int>, MAXNFILES> l2_day;           // Day of year values
std::array<std::vector<int>, MAXNFILES> l2_msec;          // Milliseconds of day

// Arrays to store geolocation data
std::array<std::vector<float>, MAXNFILES> l2_latitude;   // Latitude values
std::array<std::vector<float>, MAXNFILES> l2_longitude;  // Longitude values

// Arrays to store sensor-specific information
std::array<std::vector<uint8_t>, MAXNFILES> mside;   // Mirror side
std::array<std::vector<uint8_t>, MAXNFILES> detnum;  // Detector number
std::array<std::vector<int32_t>, MAXNFILES> pixnum;  // Pixel number

// Index tracking which file is currently being processed
size_t current_file_index = 0;
}  // namespace

int32_t get_dtype(int32_t dtype, ds_format_t fileformat) {
    if (fileformat == DS_NCDF) {
        if (dtype == DFNT_INT8)
            dtype = NC_BYTE;
        else if (dtype == DFNT_UINT8)
            dtype = NC_UBYTE;
        else if (dtype == DFNT_INT16)
            dtype = NC_SHORT;
        else if (dtype == DFNT_UINT16)
            dtype = NC_USHORT;
        else if (dtype == DFNT_INT32)
            dtype = NC_INT;
        else if (dtype == DFNT_UINT32)
            dtype = NC_UINT;
        else if (dtype == DFNT_FLOAT32)
            dtype = NC_FLOAT;
        else if (dtype == DFNT_FLOAT64)
            dtype = NC_DOUBLE;
        else {
            printf("-E- %s:%d - unknown data type.\n", __FILE__, __LINE__);
            exit(EXIT_FAILURE);
        }
    }
    return dtype;
}
//----------------------------------------------------------------

int32_t openL2(const char *fname, const char *plist, l2_prod *l2_str) {
    l2_reader[current_file_index] = L2_Reader(fname);
    // put geolocation check for set Variables.
    l2_reader[current_file_index].iniGeolocation();
    std::vector<std::string> products_requested_l3b;
    std::string product_list = "ALL";
    if (plist)
        product_list = plist;
    std::string wave_list = "ALL";
    l2_reader[current_file_index].setVariables(product_list, wave_list, products_requested_l3b);
    int status = l2_reader[current_file_index].ini_l2_flags();
    l2_reader[current_file_index].set_l2_metadata();
    size_t nsamp, nrec;
    l2_reader[current_file_index].getDimensions(nrec, nsamp);
    l2_str->nrec = nrec;
    l2_str->nsamp = nsamp;
    l2_str->sday = l2_reader[current_file_index].get_l2_metadata().sday;
    l2_str->syear = l2_reader[current_file_index].get_l2_metadata().syear;
    l2_str->smsec = l2_reader[current_file_index].get_l2_metadata().smsec;
    l2_str->eday = l2_reader[current_file_index].get_l2_metadata().eday;
    l2_str->eyear = l2_reader[current_file_index].get_l2_metadata().eyear;
    l2_str->emsec = l2_reader[current_file_index].get_l2_metadata().emsec;
    l2_str->year = l2_str->syear;
    l2_str->day = l2_str->sday;
    l2_str->msec = l2_str->smsec;
    l2_year[current_file_index] = l2_reader[current_file_index].get_l2_metadata().year;
    l2_day[current_file_index] = l2_reader[current_file_index].get_l2_metadata().day;
    l2_msec[current_file_index] = l2_reader[current_file_index].get_l2_metadata().msec;
    l2_str->flagnames = nullptr;
    if (status == 0) {
        constexpr size_t flagnames_size = 1024;
        l2_str->flagnames = new char[flagnames_size];
        auto bits_flags = l2_reader[current_file_index].get_l2_meaning_bit_dict();
        size_t index = 0;
        for (const auto &[meaning, flag] : bits_flags) {
            if (index == 0) {
                strncpy(l2_str->flagnames, meaning.c_str(), flagnames_size - 1);
            } else {
                strncat(l2_str->flagnames, ("," + meaning).c_str(),
                        flagnames_size - strlen(l2_str->flagnames) - 1);
            }
            index++;
        }
    }
    l2_data[current_file_index].resize(products_requested_l3b.size());
    l2_scan_time[current_file_index] = l2_reader[current_file_index].get_l2_metadata().scan_time;
    l2_str->l2_data = new float *[products_requested_l3b.size()];
    l2_str->nprod = products_requested_l3b.size();
    for (int i = 0; i < l2_str->nprod; i++) {
        l2_str->prodname[i] = strdup(products_requested_l3b[i].c_str());
    }
    l2_str->ntilts = 1;
    l2_str->tilt_flags[0] = 0;
    l2_str->tilt_ranges[0][0] = 0;
    l2_str->tilt_ranges[1][0] = 0;
    l2_latitude[current_file_index].resize(nsamp);
    l2_longitude[current_file_index].resize(nsamp);
    strncpy(l2_str->filename, l2_reader[current_file_index].get_filename().c_str(),
            sizeof(l2_str->filename) - 1);
    l2_str->fileindex = current_file_index;

    netCDF::NcVar var = l2_reader[current_file_index].get_variable("mside");
    if (!var.isNull()) {
        mside[current_file_index].resize(nrec);
        var.getVar(mside[current_file_index].data());
        l2_str->mside = mside[current_file_index].data();
    }
    var = l2_reader[current_file_index].get_variable("pixnum");
    if (!var.isNull()) {
        pixnum[current_file_index].resize(nrec * nsamp);
        var.getVar(pixnum[current_file_index].data());
        l2_str->pixnum = pixnum[current_file_index].data();
    }
    var = l2_reader[current_file_index].get_variable("detnum");
    if (!var.isNull()) {
        detnum[current_file_index].resize(nrec);
        var.getVar(detnum[current_file_index].data());
        l2_str->detnum = detnum[current_file_index].data();
    }
    // increment current file index
    current_file_index++;
    return 0;
}

int32_t reopenL2(int32_t fileindex, l2_prod *l2_str) {
    l2_reader[fileindex].reopenL2();
    return 0;
}

int32_t readL2(l2_prod *l2_str, int32_t ifile, int32_t recnum, int32_t iprod,
               unsigned char *scan_in_rowgroup) {
    std::vector<u_char> mask;
    if (scan_in_rowgroup) {
        mask.resize(l2_str->nrec);
        std::copy(scan_in_rowgroup, scan_in_rowgroup + l2_str->nrec, mask.begin());
        l2_reader[ifile].set_cache_flag(true);
    } else {
        // don't use cache, the user is inteded to read random chuncks or one scan at a time
        l2_reader[ifile].set_cache_flag(false);
    }
    if (iprod != -1) {
        std::vector<size_t> start = {(size_t)recnum, 0};
        std::vector<size_t> count = {1, (size_t)l2_str->nsamp};
        l2_reader[ifile].readL2data(&l2_str->l2_data[iprod], start, count, iprod, mask, false);
    } else {
        l2_reader[ifile].readL2dataScan(l2_data[ifile], recnum, mask);
        for (size_t i = 0; i < l2_data[ifile].size(); i++) {
            l2_str->l2_data[i] = l2_data[ifile][i];
        }
    }
    l2_reader[ifile].readL2FlagsScan(&l2_str->l2_flags, recnum, mask);
    l2_reader[ifile].readLatitudeScan(&l2_str->latitude, recnum, mask);
    l2_reader[ifile].readLongitudeScan(&l2_str->longitude, recnum, mask);
    if (!l2_scan_time[ifile].empty()) {
        double scan_time = l2_scan_time[ifile][recnum];
        int16_t year, day;
        double sec;
        unix2yds(scan_time, &year, &day, &sec);
        l2_str->year = year;
        l2_str->day = day;
        l2_str->msec = std::roundl(sec * 1000);
    }
    if (!l2_year[ifile].empty()) {
        l2_str->year = l2_year[ifile][recnum];
    }
    if (!l2_day[ifile].empty()) {
        l2_str->day = l2_day[ifile][recnum];
    }
    if (!l2_msec[ifile].empty()) {
        l2_str->msec = l2_msec[ifile][recnum];
    }
    return 0;
}

int32_t readlonlat(l2_prod *l2_str, int32_t ifile, int32_t *start, int32_t *edges,
                   unsigned char *scan_in_rowgroup) {
    std::vector<size_t> start_nc = {(size_t)start[0], (size_t)start[1]};
    std::vector<size_t> count = {(size_t)edges[0], (size_t)edges[1]};
    l2_reader[ifile].readLatitude(l2_latitude[ifile].data(), start_nc, count);
    l2_reader[ifile].readLongitude(l2_longitude[ifile].data(), start_nc, count);
    l2_str->latitude = l2_latitude[ifile].data();
    l2_str->longitude = l2_longitude[ifile].data();
    return 0;
}

int32_t closeL2(l2_prod *l2_str, int32_t ifile) {
    l2_reader[ifile].reset_cache();
    return 0;
}

int32_t freeL2(l2_prod *l2_str) {
    if(!l2_str)
        return 1;
    delete[] l2_str->l2_data;
    delete[] l2_str->flagnames;
    l2_str->l2_data = nullptr;
    l2_str->flagnames = nullptr;
    for (int i = 0; i < l2_str->nprod; i++) {
        free(l2_str->prodname[i]);
    }
    return 0;
}

int32_t findprod(l2_prod *l2_str, char *prodname) {
    size_t index;
    int status = l2_reader[current_file_index].find_product_index(prodname, index);
    if (status != 0) {
        return -1;
    }
    return index;
}

int32_t readL2meta(meta_l2Type *meta_l2, int32_t ifile) {
    meta_l2->title = strdup(l2_reader[ifile].get_l2_metadata().title.c_str());
    meta_l2->sensor = strdup(l2_reader[ifile].get_l2_metadata().sensor_name.c_str());
    meta_l2->mission = strdup(l2_reader[ifile].get_l2_metadata().mission.c_str());
    meta_l2->sensor_name = strdup(instrumentPlatform2SensorName(meta_l2->sensor, meta_l2->mission));
    meta_l2->westlon = l2_reader[ifile].get_l2_metadata().westlon;
    meta_l2->eastlon = l2_reader[ifile].get_l2_metadata().eastlon;
    meta_l2->southlat = l2_reader[ifile].get_l2_metadata().southlat;
    meta_l2->northlat = l2_reader[ifile].get_l2_metadata().northlat;
    return 0;
}

int32_t freeL2meta(meta_l2Type *meta_l2) {
#define FREE(ptr)     \
    if ((ptr) != 0x0) \
        free(ptr);
    FREE(meta_l2->title);
    FREE(meta_l2->sensor_name);
    FREE(meta_l2->mission);
    FREE(meta_l2->sensor);
    return 0;
}

int32_t getL3units(l2_prod *l2_str, int32_t ifile, char *l3b_prodname, char *units) {
    size_t product_index;
    int status = l2_reader[ifile].find_product_index(l3b_prodname, product_index);
    if (status != 0) {
        return -1;
    }
    units = strdup(l2_reader[ifile].get_units()[product_index].c_str());
    return 0;
}
