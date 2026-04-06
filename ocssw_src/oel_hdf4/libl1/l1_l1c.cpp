/* ============================================================================ */
/* module l1_l1c.c - functions to read L1C for l2gen (l1info really)            */
/* Written By:  Don Shea                                                        */
/*                                                                              */
/* ============================================================================ */

#include "l1_l1c.h"
#include <netcdf>
#include <nc4utils.h>
#include "libnav.h"
#include <cstdio>
#include <cmath>
#include <memory>
#include <string>
#include <fstream>
#include <scaled_nc_var.hpp>
#include <regex>
#include <algorithm>

static bool fileExists(const std::string &filename) {
    std::ifstream file(filename);
    return file.good();
}

void exit_verbose(const char *msg, int line, const char *file) {
    fprintf(stderr, "-E- %s:%d: %s\n", file, line, msg);
    exit(EXIT_FAILURE);
}

class L1C_FileReader {
   protected:
    bool file_open = false;
    netCDF::NcFile file;
    ScaledNcVar latitude;
    ScaledNcVar longitude;
    ScaledNcVar sensor_azimuth_angle;
    ScaledNcVar sensor_zenith_angle;
    ScaledNcVar solar_azimuth_angle;
    ScaledNcVar solar_zenith_angle;
    ScaledNcVar sensor_view_angle;
    ScaledNcVar scan_quality_flags;
    double start_time;
    ScaledNcVar i;
    ScaledNcVar qc;
    ScaledNcVar f0;
    ScaledNcVar height;
    ScaledNcVar wavelengths;
    netCDF::NcDim first_dimension;
    netCDF::NcDim second_dimension;
    netCDF::NcDim bands_dimension;
    netCDF::NcDim view_dimension;
    size_t num_scans, num_pixels;
    size_t num_bands, num_views;
    std::vector<double> scan_time;
    std::vector<float> tilt;
    std::vector<unsigned char> quality_geolocation;
   public:
    L1C_FileReader() = default;
    explicit L1C_FileReader(const std::string &filepath) {
        char msg[1024];
        if (!fileExists(filepath)) {
            sprintf(msg, "File %s does not exist\n", filepath.c_str());
            exit_verbose(msg, __LINE__, __FILE__);
        }
        try {
            file.open(filepath, netCDF::NcFile::read);
            file_open = true;
        } catch (netCDF::exceptions::NcException &e) {
            sprintf(msg, "File %s is not a netCDF file, %s \n", filepath.c_str(), e.what());
            exit_verbose(msg, __LINE__, __FILE__);
        }
        // find latitude and longitude
        // first find geolocation_data group
        netCDF::NcGroup geolocation_data = file.getGroup("geolocation_data");
        if (geolocation_data.isNull()) {
            sprintf(msg, "File %s does not contain geolocation_data group\n", filepath.c_str());
            exit_verbose(msg, __LINE__, __FILE__);
        }
        latitude = geolocation_data.getVar("latitude");
        if (latitude.isNull()) {
            sprintf(msg, "File %s does not contain latitude variable\n", filepath.c_str());
            exit_verbose(msg, __LINE__, __FILE__);
        }
        longitude = geolocation_data.getVar("longitude");
        if (longitude.isNull()) {
            sprintf(msg, "File %s does not contain longitude variable\n", filepath.c_str());
            exit_verbose(msg, __LINE__, __FILE__);
        }
        // check their dimensions
        auto dimensions = latitude.getDims();
        if (dimensions.size() != 2) {
            sprintf(msg, "File %s latitude variable has %zu dimensions, expected 2\n", filepath.c_str(),
                    dimensions.size());
            exit_verbose(msg, __LINE__, __FILE__);
        }
        first_dimension = dimensions[0];
        second_dimension = dimensions[1];
        dimensions = longitude.getDims();
        if (dimensions.size() != 2) {
            sprintf(msg, "File %s longitude variable has %zu dimensions, expected 2\n", filepath.c_str(),
                    dimensions.size());
            exit_verbose(msg, __LINE__, __FILE__);
        }
        if (first_dimension != dimensions[0]) {
            sprintf(msg, "File %s latitude and longitude dimensions do not match\n", filepath.c_str());
            exit_verbose(msg, __LINE__, __FILE__);
        }
        if (second_dimension != dimensions[1]) {
            sprintf(msg, "File %s latitude and longitude dimensions do not match\n", filepath.c_str());
            exit_verbose(msg, __LINE__, __FILE__);
        }
        height = geolocation_data.getVar("height");
        if (height.isNull()) {
            sprintf(msg, "File %s does not contain height variable\n", filepath.c_str());
            exit_verbose(msg, __LINE__, __FILE__);
        }
        dimensions = height.getDims();
        if (dimensions.size() != 2) {
            sprintf(msg, "File %s height variable has %zu dimensions, expected 2\n", filepath.c_str(),
                    dimensions.size());
            exit_verbose(msg, __LINE__, __FILE__);
        }
        if (first_dimension != dimensions[0]) {
            sprintf(msg, "File %s height and geolocation_data dimensions do not match\n", filepath.c_str());
            exit_verbose(msg, __LINE__, __FILE__);
        }
        if (second_dimension != dimensions[1]) {
            sprintf(msg, "File %s height and geolocation_data dimensions do not match\n", filepath.c_str());
            exit_verbose(msg, __LINE__, __FILE__);
        }

        // getting sensor_azimuth_angle
        sensor_azimuth_angle = geolocation_data.getVar("sensor_azimuth_angle");
        if (sensor_azimuth_angle.isNull()) {
            sprintf(msg, "File %s does not contain sensor_azimuth_angle variable\n", filepath.c_str());
            exit_verbose(msg, __LINE__, __FILE__);
        }
        dimensions = sensor_azimuth_angle.getDims();
        if (dimensions.size() != 3) {
            sprintf(msg, "File %s sensor_azimuth_angle variable has %zu dimensions, expected 3\n",
                    filepath.c_str(), dimensions.size());
            exit_verbose(msg, __LINE__, __FILE__);
        }
        if (first_dimension != dimensions[0]) {
            sprintf(msg, "File %s sensor_azimuth_angle and geolocation_data dimensions do not match\n",
                    filepath.c_str());
            exit_verbose(msg, __LINE__, __FILE__);
        }
        if (second_dimension != dimensions[1]) {
            sprintf(msg, "File %s sensor_azimuth_angle and geolocation_data dimensions do not match\n",
                    filepath.c_str());
            exit_verbose(msg, __LINE__, __FILE__);
        }
        view_dimension = dimensions[2];
        // getting sensor_zenith_angle
        sensor_zenith_angle = geolocation_data.getVar("sensor_zenith_angle");
        if (sensor_zenith_angle.isNull()) {
            sprintf(msg, "File %s does not contain sensor_zenith_angle variable\n", filepath.c_str());
            exit_verbose(msg, __LINE__, __FILE__);
        }
        if (sensor_zenith_angle.getDims().size() != 3) {
            sprintf(msg, "File %s sensor_zenith_angle variable has %lu dimensions, expected 3\n",
                    filepath.c_str(), sensor_zenith_angle.getDims().size());
            exit_verbose(msg, __LINE__, __FILE__);
        }
        for (size_t idim = 0; idim < dimensions.size(); idim++) {
            if (dimensions[idim] != sensor_zenith_angle.getDims()[idim]) {
                sprintf(msg, "File %s sensor_zenith_angle and sensor_azimuth_angle dimensions do not match\n",
                        filepath.c_str());
                exit_verbose(msg, __LINE__, __FILE__);
            }
        }
        // getting solar_azimuth_angle
        solar_azimuth_angle = geolocation_data.getVar("solar_azimuth_angle");
        if (solar_azimuth_angle.isNull()) {
            sprintf(msg, "File %s does not contain solar_azimuth_angle variable\n", filepath.c_str());
            exit_verbose(msg, __LINE__, __FILE__);
        }
        if (solar_azimuth_angle.getDims().size() != 3) {
            sprintf(msg, "File %s solar_azimuth_angle variable has %zu dimensions, expected 3\n",
                    filepath.c_str(), solar_azimuth_angle.getDims().size());
            exit_verbose(msg, __LINE__, __FILE__);
        }
        for (size_t idim = 0; idim < dimensions.size(); idim++) {
            if (dimensions[idim] != solar_azimuth_angle.getDims()[idim]) {
                sprintf(msg, "File %s solar_azimuth_angle and sensor_azimuth_angle dimensions do not match\n",
                        filepath.c_str());
                exit_verbose(msg, __LINE__, __FILE__);
            }
        }

        // getting solar_zenith_angle
        solar_zenith_angle = geolocation_data.getVar("solar_zenith_angle");
        if (solar_zenith_angle.isNull()) {
            sprintf(msg, "File %s does not contain solar_zenith_angle variable\n", filepath.c_str());
            exit_verbose(msg, __LINE__, __FILE__);
        }
        if (solar_zenith_angle.getDims().size() != 3) {
            sprintf(msg, "File %s solar_zenith_angle variable has %zu dimensions, expected 3\n",
                    filepath.c_str(), solar_zenith_angle.getDims().size());
            exit_verbose(msg, __LINE__, __FILE__);
        }
        for (size_t idim = 0; idim < dimensions.size(); idim++) {
            if (dimensions[idim] != solar_zenith_angle.getDims()[idim]) {
                sprintf(msg, "File %s solar_zenith_angle and sensor_azimuth_angle dimensions do not match\n",
                        filepath.c_str());
            }
        }

        // find observation_data group
        netCDF::NcGroup observation_data = file.getGroup("observation_data");
        if (observation_data.isNull()) {
            sprintf(msg, "File %s does not contain observation_data group\n", filepath.c_str());
            exit_verbose(msg, __LINE__, __FILE__);
        }
        i = observation_data.getVar("i");
        dimensions = i.getDims();
        if (dimensions.size() != 4) {
            sprintf(msg, "File %s i variable has %zu dimensions, expected 4\n", filepath.c_str(),
                    dimensions.size());
            exit_verbose(msg, __LINE__, __FILE__);
        }
        if (first_dimension.getSize() != dimensions[0].getSize()) {
            sprintf(msg, "File %s i and observation_data dimensions do not match\n", filepath.c_str());
            exit_verbose(msg, __LINE__, __FILE__);
        }
        if (second_dimension.getSize() != dimensions[1].getSize()) {
            sprintf(msg, "File %s i and observation_data dimensions do not match\n", filepath.c_str());
            exit_verbose(msg, __LINE__, __FILE__);
        }
        if (view_dimension.getSize() != dimensions[2].getSize()) {
            sprintf(msg, "File %s i and observation_data dimensions do not match\n", filepath.c_str());
            exit_verbose(msg, __LINE__, __FILE__);
        }
        bands_dimension = dimensions[3];
        if (bands_dimension.getName().find("band") == std::string::npos) {
            sprintf(msg, "File %s bands_dimension is not named band\n", filepath.c_str());
            exit_verbose(msg, __LINE__, __FILE__);
        }
        // find qc
        qc = observation_data.getVar("qc");
        if (qc.isNull()) {
            sprintf(msg, "File %s does not contain qc variable\n", filepath.c_str());
            exit_verbose(msg, __LINE__, __FILE__);
        }
        dimensions = qc.getDims();
        if (dimensions.size() != 4) {
            sprintf(msg, "File %s qc variable has %zu dimensions, expected 4\n", filepath.c_str(),
                    dimensions.size());
            exit_verbose(msg, __LINE__, __FILE__);
        }
        if (first_dimension.getSize() != dimensions[0].getSize()) {
            sprintf(msg, "File %s qc and observation_data dimensions do not match\n", filepath.c_str());
            exit_verbose(msg, __LINE__, __FILE__);
        }
        if (second_dimension.getSize() != dimensions[1].getSize()) {
            sprintf(msg, "File %s qc and observation_data dimensions do not match\n", filepath.c_str());
            exit_verbose(msg, __LINE__, __FILE__);
        }
        if (view_dimension.getSize() != dimensions[2].getSize()) {
            sprintf(msg, "File %s qc and observation_data dimensions do not match\n", filepath.c_str());
            exit_verbose(msg, __LINE__, __FILE__);
        }
        if (bands_dimension.getSize() != dimensions[3].getSize()) {
            sprintf(msg, "File %s qc and observation_data dimensions do not match\n", filepath.c_str());
            exit_verbose(msg, __LINE__, __FILE__);
        }

        // find sensor_views_bands group
        netCDF::NcGroup sensor_views_bands = file.getGroup("sensor_views_bands");
        if (sensor_views_bands.isNull()) {
            sprintf(msg, "File %s does not contain sensor_views_bands group\n", filepath.c_str());
            exit_verbose(msg, __LINE__, __FILE__);
        }

        // find intensity_f0
        f0 = sensor_views_bands.getVar("intensity_f0");

        dimensions = f0.getDims();
        if (dimensions.size() != 2) {
            sprintf(msg, "File %s f0 variable has %zu dimensions, expected 2\n", filepath.c_str(),
                    dimensions.size());
            exit_verbose(msg, __LINE__, __FILE__);
        }
        if (view_dimension.getSize() != dimensions[0].getSize()) {
            sprintf(msg, "File %s f0 and sensor_views_bands dimensions do not match\n", filepath.c_str());
            exit_verbose(msg, __LINE__, __FILE__);
        }
        if (bands_dimension.getSize() != dimensions[1].getSize()) {
            sprintf(msg, "File %s f0 and sensor_views_bands dimensions do not match\n", filepath.c_str());
            exit_verbose(msg, __LINE__, __FILE__);
        }

        // find intensity_wavelength
        wavelengths = sensor_views_bands.getVar("intensity_wavelength");
        dimensions = wavelengths.getDims();
        if (dimensions.size() != 2) {
            sprintf(msg, "File %s wavelengths variable has %zu dimensions, expected 2\n", filepath.c_str(),
                    dimensions.size());
            exit_verbose(msg, __LINE__, __FILE__);
        }
        if (view_dimension.getSize() != dimensions[0].getSize()) {
            sprintf(msg, "File %s wavelengths and sensor_views_bands dimensions do not match\n",
                    filepath.c_str());
            exit_verbose(msg, __LINE__, __FILE__);
        }
        if (bands_dimension.getSize() != dimensions[1].getSize()) {
            sprintf(msg, "File %s wavelengths and sensor_views_bands dimensions do not match\n",
                    filepath.c_str());
            exit_verbose(msg, __LINE__, __FILE__);
        }
        // reading group bin_attributes
        netCDF::NcGroup bin_attributes = file.getGroup("bin_attributes");
        if (bin_attributes.isNull()) {
            sprintf(msg, "File %s does not contain bin_attributes group\n", filepath.c_str());
            exit_verbose(msg, __LINE__, __FILE__);
        }
        // get var nadir_view_time
        netCDF::NcVar nadir_view_time = bin_attributes.getVar("nadir_view_time");
        if (nadir_view_time.isNull()) {
            sprintf(msg, "File %s does not contain nadir_view_time variable\n", filepath.c_str());
            exit_verbose(msg, __LINE__, __FILE__);
        }
        dimensions = nadir_view_time.getDims();
        if (dimensions.size() != 1) {
            sprintf(msg, "File %s nadir_view_time variable has %zu dimensions, expected 1\n",
                    filepath.c_str(), dimensions.size());
            exit_verbose(msg, __LINE__, __FILE__);
        }
        if (first_dimension.getSize() != dimensions[0].getSize()) {
            sprintf(msg, "File %s nadir_view_time and bin_attributes dimensions do not match\n",
                    filepath.c_str());
            exit_verbose(msg, __LINE__, __FILE__);
        }
        // read units string attribute from nadir_view_time
        std::string units;
        netCDF::NcVarAtt units_nc = nadir_view_time.getAtt("units");
        if (units_nc.isNull()) {
            sprintf(msg, "File %s does not contain units attribute for nadir_view_time variable\n",
                    filepath.c_str());
            exit_verbose(msg, __LINE__, __FILE__);
        }
        units_nc.getValues(units);

        std::regex pattern(R"((\d{4})-(\d{2})-(\d{2}))");
        std::smatch matches;
        if (std::regex_search(units, matches, pattern)) {
            int16_t year, month, day;
            year = std::stoi(matches[1]);
            month = std::stoi(matches[2]);
            day = std::stoi(matches[3]);
            start_time = ymds2unix(year, month, day, 0.0);
        } else {
            // get time_coverage_start global attribute
            std::string time_coverage_start;
            netCDF::NcGroupAtt time_coverage_start_nc = file.getAtt("time_coverage_start");
            if (time_coverage_start_nc.isNull()) {
                sprintf(msg, "File %s does not contain time_coverage_start attribute\n", filepath.c_str());
                exit_verbose(msg, __LINE__, __FILE__);
            }
            time_coverage_start_nc.getValues(time_coverage_start);
            double starTtime = isodate2unix(time_coverage_start.c_str());
            int16_t syear, smon, sday;
            double secs;
            unix2ymds(starTtime, &syear, &smon, &sday, &secs);
            start_time = ymds2unix(syear, smon, sday, 0.0);
        }

        // OCI specific
        num_scans = first_dimension.getSize();
        num_pixels = second_dimension.getSize();
        num_bands = bands_dimension.getSize();
        num_views = view_dimension.getSize();
        scan_time.resize(num_scans);
        nadir_view_time.getVar(scan_time.data());
        // optional
        scan_quality_flags = bin_attributes.getVar("scan_quality_flags");
        if (!scan_quality_flags.isNull()) {
            if (scan_quality_flags.getDimCount() !=1) {
                sprintf(msg, "File %s scan_quality_flags variable has %d dimensions, expected 1\n",
                        filepath.c_str(), scan_quality_flags.getDimCount());
                exit_verbose(msg, __LINE__, __FILE__);
            }
            if (scan_quality_flags.getDims()[0].getSize() != num_scans) {
                sprintf(msg, "File %s scan_quality_flags and bin_attributes  dimensions do not match\n",
                        filepath.c_str());
            }
            quality_geolocation.resize(num_scans);
            scan_quality_flags.getVar(quality_geolocation.data());
        }
        sensor_view_angle = sensor_views_bands.getVar("sensor_view_angle");
        if (!sensor_view_angle.isNull()) {
            if (sensor_view_angle.getDimCount() != 1) {
                sprintf(msg, "File %s sensor_view_angle variable has %d dimensions, expected 1\n",filepath.c_str(),sensor_view_angle.getDimCount());
                exit_verbose(msg, __LINE__, __FILE__);
            }
            if (sensor_view_angle.getDims()[0].getSize() != num_scans) {
                sprintf(msg, "File %s sensor_view_angle and sensor_views_bands dimensions do not match\n",
                        filepath.c_str());
            }
            tilt.resize(num_scans);
            sensor_view_angle.getVar(tilt.data());
        }



    }
    explicit L1C_FileReader(const std::string &filepath, filehandle *file) : L1C_FileReader(filepath) {
        setl1cfile(file);
    }
    virtual void setl1cfile(filehandle *file_l2) {
        file_l2->sd_id = file.getId();
        file_l2->nbands = num_bands;
        file_l2->npix = num_pixels;
        file_l2->nscan = num_scans;
        file_l2->ndets = 1;
        file_l2->terrain_corrected = 1;  // presumed.
        strcpy(file_l2->spatialResolution, "5.2km");
    }
    virtual void closeFile() {
        if (file_open)
            file.close();
        file_open = false;
    }
    // OCI specific
    virtual void readScan(size_t scan, l1str *l1rec) {
        for (size_t ip = 0; ip < num_pixels; ip++) {
            l1rec->pixnum[ip] = ip;
        }
        if (!tilt.empty()) {
            l1rec->tilt = tilt[scan];
        }
        if (!quality_geolocation.empty()) {
            if (quality_geolocation[scan] & 1)
                for (size_t i = 0; i < num_pixels; i++)
                    l1rec->navwarn[i] = 1;
        }
        l1rec->npix = num_pixels;
        if (scan_time[scan] == BAD_FLT) {
            l1rec->scantime = BAD_FLT;
            l1rec->fsol = BAD_FLT;
        } else {
            l1rec->scantime = start_time + scan_time[scan];

            int16_t syear, sday;
            double secs;
            unix2yds(l1rec->scantime, &syear, &sday, &secs);

            int32_t yr = syear;
            int32_t dy = sday;
            int32_t msec = (int32_t)(secs * 1000.0);
            double esdist = esdist_(&yr, &dy, &msec);

            l1rec->fsol = pow(1.0 / esdist, 2);
        }
        std::vector<size_t> start = {scan, 0};
        std::vector<size_t> count = {1, num_pixels};
        try {
            latitude.getVar(start, count, l1rec->lat);
            longitude.getVar(start, count, l1rec->lon);
            height.getVar(start, count, l1rec->height);
        } catch (netCDF::exceptions::NcException &e) {
            fprintf(stderr, "Error reading latitude/longitude/height variables: %s\n", e.what());
            exit(1);
        }
        size_t valid_view = 0;

        start = {scan, 0, valid_view};
        count = {1, num_pixels, 1};
        sensor_zenith_angle.getVar(start, count, l1rec->senz);

        bool all_fill =
            std::all_of(l1rec->senz, l1rec->senz + num_pixels, [](float x) { return x == BAD_FLT; });
        if (all_fill)
            valid_view = std::min(1lu,num_views - 1);

        start = {scan, 0, valid_view, 0};
        count = {1, num_pixels, 1, num_bands};
        try {
            i.getVar(start, count, l1rec->Lt);
        } catch (netCDF::exceptions::NcException &e) {
            fprintf(stderr, "Error reading i variable: %s\n", e.what());
            exit(1);
        }

        for (size_t ip = 0; ip < num_pixels; ip++) {
            for (size_t ib = 0; ib < num_bands; ib++) {
                size_t inb = num_bands * ip + ib;
                if (l1rec->Lt[inb]!=BAD_FLT)
                    l1rec->Lt[inb] = l1rec->Lt[inb] / 10;
            }
        }

        start = {scan, 0, valid_view};
        count = {1, num_pixels, 1};
        try {
            sensor_azimuth_angle.getVar(start, count, l1rec->sena);
            solar_azimuth_angle.getVar(start, count, l1rec->sola);
            sensor_zenith_angle.getVar(start, count, l1rec->senz);
            solar_zenith_angle.getVar(start, count, l1rec->solz);
        } catch (netCDF::exceptions::NcException &e) {
            fprintf(stderr,
                    "Error reading "
                    "sensor_azimuth_angle/solar_azimuth_angle/sensor_zenith_angle/solar_zenith_angle "
                    "variables: %s\n",
                    e.what());
            exit(1);
        }
        // set navigation fail
        for (size_t ip = 0; ip < num_pixels; ip++) {
            l1rec->navfail[ip] = 0;
            if (l1rec->sena[ip] == BAD_FLT) {
                l1rec->navfail[ip] = 1;
            }
            if (l1rec->sola[ip] == BAD_FLT) {
                l1rec->navfail[ip] = 1;
            }
            if (l1rec->senz[ip] == BAD_FLT) {
                l1rec->navfail[ip] = 1;
            }
            if (l1rec->solz[ip] == BAD_FLT) {
                l1rec->navfail[ip] = 1;
            }
            if (l1rec->lat[ip] == BAD_FLT) {
                l1rec->navfail[ip] = 1;
            }
            if (l1rec->lon[ip] == BAD_FLT) {
                l1rec->navfail[ip] = 1;
            }
            if (l1rec->navfail[ip]) {
                l1rec->sena[ip] = BAD_FLT;
                l1rec->sola[ip] = BAD_FLT;
                l1rec->senz[ip] = BAD_FLT;
                l1rec->solz[ip] = BAD_FLT;
                l1rec->lat[ip] = BAD_FLT;
                l1rec->lon[ip] = BAD_FLT;
            }
        }
    }
    virtual ~L1C_FileReader() {
        if (file_open)
            file.close();
    }
};

std::unique_ptr<L1C_FileReader> l1file;

/**
 * Get
 * @param file input l1c oci file
 * @return
 */
int openl1_l1c(filehandle *file) {
    l1file = std::make_unique<L1C_FileReader>(file->name, file);
    return LIFE_IS_GOOD;
}

/**
 *
 * @param file input l1c oci file
 * @param line scan
 * @param l1rec l1 record to fill
 * @return
 */
int readl1_l1c(filehandle *file, int32_t line, l1str *l1rec) {
    l1file->readScan(line, l1rec);
    return LIFE_IS_GOOD;
}

/**
 * @param file closes l1 file handle
 * @return
 */
int closel1_l1c(filehandle *file) {
    l1file->closeFile();

    return LIFE_IS_GOOD;
}
