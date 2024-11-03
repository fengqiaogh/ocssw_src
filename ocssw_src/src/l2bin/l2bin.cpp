#include <cstdlib>
#include <cstdint>
#include <libgen.h>
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <netcdf>
#include <boost/algorithm/string.hpp>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/foreach.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <get_dataday.hpp>
#include <meta_l3b.h>
#include <seaproto.h>
#include <readL2scan.h>
#include <timeutils.h>
#include <genutils.h>
#include <setupflags.h>
#include <sensorInfo.h>
#include <nc4utils.h>
#include <ncdfbin_utils.h>
#include <chrono>
#include <hdf.h>
#include "l2bin_input.h"
#include "L3Shape.h"
#include "L3ShapeIsine.h"
#include <get_geospatial.hpp>
#include "expand3D.hpp"
#include "find_variable.hpp"

// using namespace std;
// using namespace netCDF;
// using namespace netCDF::exceptions;

namespace bg = boost::geometry;

// #define DEBUG_PRINT
// #define DEBUG_PRINT2
// extern void cdata_();
#ifdef DEBUG_PRINT
bool enableDebugPrint = true;
#endif

#define MTILT_DIMS_2 20
#define LTILT_DIMS_2 2
#define MAXALLOCPERBIN 20


static constexpr int max_l3b_products = MAXNUMBERPRODUCTS;

/* Global variables */
static instr input;
static l2_prod l2_str[MAXNFILES];
// bin file shape
static int32_t nrows = -1;
static l3::L3Shape *shape;
static int64_t *basebin;

// global variables for the group
int32_t n_allocperbin;    // number of obs to allocate to make another chunk
int16_t *allocated_space; // number of obs allocated for each bin
static int16_t *nobs;     // number of observations in each bin

static float32 **data_values;
static double **data_areas; // in fractions of a bin
static int16_t **file_index;
static uint8_t **data_quality;
static float64 **time_value;
static double time_tai93; // TAI 93 time for the current L2 line

int32_t *bin_flag;
int16_t *tilt, *qual, *nscenes, *lastfile;

int32_t n_bins_in_group;      // number of bins in group
int32_t n_rows_in_group = -1; // number of rows in group
int32_t krow;                 // first bin row number of group

// other globals
int32_t tiltstate = 0;
int32_t l3b_nprod; // number of products in bin file
static char buf[LG_ATTRSZ];

int32_t sensorID[MAXNFILES];

int16_t *numer[MAXNFILES], *denom[MAXNFILES];
int16_t qual_prod_index[MAXNFILES];
int16_t composite_prod_index[MAXNFILES];
int16_t composite_l3prod_index = -1;

float32 f32;
std::vector<int32_t> thirdDimId;
std::vector<float32> min_value, max_value;

typedef bg::model::point<double, 2, bg::cs::geographic<bg::degree>> Point_t;
typedef bg::model::polygon<Point_t> Polygon_t;
typedef bg::model::box<Point_t> Box_t;

#define VERSION "7.1.0"
#define PROGRAM "l2bin"

int32_t get_l2prod_index(const l2_prod &l2, const char *prodname)
{
    int32_t index;
    for (index = 0; index < l2.nprod; index++)
        if (strcmp(prodname, l2.prodname[index]) == 0)
            break;
    if (index == l2.nprod)
        index = -1;
    return index;
}

void addPixelToBin(int32_t ifile, int32_t ipixl, uint64_t bin, bool is_l2_flags_defined = true, double areaFrac = 1.0)
{

    int32_t ibin = bin - basebin[krow];

    /* if bin not within bin row group then continue */
    /* --------------------------------------------- */
    if (ibin < 0 ||
        ibin >= n_bins_in_group ||
        (int64_t)bin >= basebin[krow + n_rows_in_group] ||
        (int64_t)bin < basebin[krow])
        return;

        /* GOOD OBSERVATION FOUND */
        /* ---------------------- */
#ifdef MALLINFO
    if (input.dcinfo)
    {
        if ((l2_str[ifile].longitude[ipixl] <= -160) ||
            (l2_str[ifile].longitude[ipixl] >= +160))
        {
            printf("DC: %10ld %12d %8.2f %8.2f\n",
                   (long int)bin,
                   (int)dbldate,
                   l2_str[ifile].longitude[ipixl],
                   l2_str[ifile].latitude[ipixl]);
        }
    }
#endif

    /* "OR" flags in swath pixel & set tilt & increment nscenes */
    /* -------------------------------------------------------- */
    if (is_l2_flags_defined)
        bin_flag[ibin] = bin_flag[ibin] | l2_str[ifile].l2_flags[ipixl];
    tilt[ibin] = tiltstate;
    if (ifile != lastfile[ibin])
    {
        nscenes[ibin]++;
        lastfile[ibin] = ifile;
    }

    /* Allocate space for file index & bin data values */
    /* ----------------------------------------------- */
    if (file_index[ibin] == NULL)
    {
        file_index[ibin] = (int16_t *)
            calloc(n_allocperbin, sizeof(int16_t));

        data_values[ibin] = (float32 *)
            calloc(n_allocperbin * l3b_nprod, sizeof(float32));
        if (data_values[ibin] == 0x0)
        {
            perror(buf);
            printf("Allocation failed for data_values[ibin]: %d %s\n",
                   ibin, buf);
            exit(EXIT_FAILURE);
        }

        data_areas[ibin] = (double *)calloc(n_allocperbin, sizeof(double));
        if (data_areas[ibin] == 0x0)
        {
            perror(buf);
            printf("Allocation failed for data_areas[ibin]: %d %s\n",
                   ibin, buf);
            exit(EXIT_FAILURE);
        }

        data_quality[ibin] = (uint8_t *)
            calloc(n_allocperbin, sizeof(uint8_t));
        if (data_quality[ibin] == 0x0)
        {
            perror(buf);
            printf("Allocation failed for data_quality[ibin]: %d %s\n",
                   ibin, buf);
            exit(EXIT_FAILURE);
        }

        time_value[ibin] = (float64 *)
            calloc(n_allocperbin, sizeof(float64));
        if (time_value[ibin] == 0x0)
        {
            perror(buf);
            printf("Allocation failed for time_value[ibin]: %d %s\n",
                   ibin, buf);
            exit(EXIT_FAILURE);
        }

        allocated_space[ibin] = n_allocperbin;
    }

    /* Set file_index for each observation */
    /* ----------------------------------- */
    file_index[ibin][nobs[ibin]] = ifile;

    /* Get data quality */
    /* ---------------- */
    if (input.qual_prod[0] != 0)
    {
        data_quality[ibin][nobs[ibin]] = l2_str[ifile].l2_data[qual_prod_index[ifile]][ipixl];
        // a temporary fix to mask the last 4 pixels of a MODIS Terra scan
        // for SST product
        if (strncmp(input.suite, "SST", 3) == 0 && sensorID[ifile] == MODIST && ipixl > 1349)
        {
            data_quality[ibin][nobs[ibin]] = 3;
        }
    }

    /* Store time_value (TAI93) */
    /* ---------------------- */
    time_value[ibin][nobs[ibin]] = time_tai93;

    /* Get composite data */
    /* ------------------- */
    if (input.composite_prod[0] != 0)
    {
        int idx = composite_prod_index[ifile];
        f32 = l2_str[ifile].l2_data[idx][ipixl];
        if (f32 == -32767)
            return;

        if (nobs[ibin] != 0)
        {
            if (strcmp(input.composite_scheme, "max") == 0)
            {
                if (f32 < data_values[ibin][composite_l3prod_index])
                    return;
            }
            else
            {
                if (f32 > data_values[ibin][composite_l3prod_index])
                    return;
            }
            nobs[ibin] = 0;
        }
    }

    /* Get data area for pixel */
    /* ----------------------------------- */
    data_areas[ibin][nobs[ibin]] = areaFrac;

    /* Get data values for all L3 products */
    /* ----------------------------------- */
    for (int32_t l3_iprod = 0; l3_iprod < l3b_nprod; l3_iprod++)
    {

        int32_t l2_iprod = numer[ifile][l3_iprod];
        f32 = l2_str[ifile].l2_data[l2_iprod][ipixl * l2_str[ifile].thirdDim[l2_iprod] + thirdDimId[l3_iprod]]; // probably the best way is to pass wavelength

        /* Set -32767 value to "bad" quality */
        if (f32 == -32767)
            if (input.qual_prod[0] != 0)
                data_quality[ibin][nobs[ibin]] = 4;

        if (input.composite_prod[0] != 0)
        {
            data_values[ibin][l3b_nprod * nobs[ibin] + l3_iprod] = f32;
        }
        else
        {
            if (denom[ifile][l3_iprod] == -1 && f32 >= min_value[l3_iprod] && f32 <= max_value[l3_iprod])
                data_values[ibin][l3b_nprod * nobs[ibin] + l3_iprod] = f32;

            if (denom[ifile][l3_iprod] == -1 && f32 < min_value[l3_iprod])
                data_values[ibin][l3b_nprod * nobs[ibin] + l3_iprod] = min_value[l3_iprod];

            if (denom[ifile][l3_iprod] == -1 && f32 > max_value[l3_iprod])
                data_values[ibin][l3b_nprod * nobs[ibin] + l3_iprod] = max_value[l3_iprod];

            if (denom[ifile][l3_iprod] == -2 && is_l2_flags_defined)
                data_values[ibin][l3b_nprod * nobs[ibin] + l3_iprod] =
                    (l2_str[ifile].l2_flags[ipixl] >>
                     numer[ifile][l3_iprod]) &
                    1;
        }

        /* ratio product */
        /* ------------- */
        if (denom[ifile][l3_iprod] >= 0)
        {

            data_values[ibin][l3b_nprod * nobs[ibin] + l3_iprod] = f32;

            l2_iprod = denom[ifile][l3_iprod];
            f32 = l2_str[ifile].l2_data[l2_iprod][ipixl * l2_str[ifile].thirdDim[l2_iprod] + thirdDimId[l3_iprod]];

            if (f32 >= min_value[l3_iprod] && f32 <= max_value[l3_iprod])
                data_values[ibin][l3b_nprod * nobs[ibin] + l3_iprod] /= f32;
            else if (f32 < min_value[l3_iprod])
                data_values[ibin][l3b_nprod * nobs[ibin] + l3_iprod] /= min_value[l3_iprod];
            else
                data_values[ibin][l3b_nprod * nobs[ibin] + l3_iprod] /= max_value[l3_iprod];
        }

    } /* iprod loop */

    /* Increment number of observations in bin */
    /* --------------------------------------- */
    nobs[ibin]++;
    if (input.composite_prod[0] != 0)
        nobs[ibin] = 1;

    /* Reallocate if necessary */
    /* ----------------------- */
    if (nobs[ibin] == allocated_space[ibin])
    {

        file_index[ibin] = (int16_t *)
            realloc(file_index[ibin],
                    (nobs[ibin] + n_allocperbin) * sizeof(int16_t));

        data_values[ibin] = (float32 *)
            realloc(data_values[ibin],
                    (nobs[ibin] + n_allocperbin) * l3b_nprod * sizeof(float32));
        if (data_values[ibin] == 0x0)
        {
            perror(buf);
            printf("Reallocation failed for data_values[ibin]: %d %s\n",
                   ibin, buf);
            exit(EXIT_FAILURE);
        }

        data_areas[ibin] = (double *)realloc(data_areas[ibin],
                                             (nobs[ibin] + n_allocperbin) * sizeof(double));
        if (data_areas[ibin] == 0x0)
        {
            perror(buf);
            printf("Reallocation failed for data_areas[ibin]: %d %s\n",
                   ibin, buf);
            exit(EXIT_FAILURE);
        }

        data_quality[ibin] = (uint8_t *)
            realloc(data_quality[ibin],
                    (nobs[ibin] + n_allocperbin) * sizeof(uint8_t));
        if (data_quality[ibin] == 0x0)
        {
            perror(buf);
            printf("Reallocation failed for data_quality[ibin]: %d %s\n",
                   ibin, buf);
            exit(EXIT_FAILURE);
        }

        time_value[ibin] = (float64 *)
            realloc(time_value[ibin],
                    (nobs[ibin] + n_allocperbin) * sizeof(float64));
        if (time_value[ibin] == 0x0)
        {
            perror(buf);
            printf("Reallocation failed for time_value[ibin]: %d %s\n",
                   ibin, buf);
            exit(EXIT_FAILURE);
        }

        allocated_space[ibin] += n_allocperbin;
    } /* end reallocate */
}

bool binIntersectsPixel(int32_t row, int32_t col, Box_t &pixelBox, double &areaFrac)
{
    bool result = false;
    double n, s, e, w;
    shape->rowcol2bounds(row, col, n, s, e, w);

#ifdef DEBUG_PRINT2
    // debug plot the bin
    if (enableDebugPrint)
    {
        printf("lat = [%f, %f, %f, %f, %f]\n", s, n, n, s, s);
        printf("lon = [%f, %f, %f, %f, %f]\n", w, w, e, e, w);
        printf("plt.plot(lon, lat)\n\n");
    }
#endif

    Box_t box(Point_t(w, s), Point_t(e, n));
    areaFrac = 0;

    if (!bg::disjoint(box, pixelBox))
    {
        Box_t output;
        if (bg::intersection(box, pixelBox, output))
        {
            double intersectArea = bg::area(output);
            if (intersectArea > 0)
            {
                result = true;
                double binArea = (n - s) * (e - w);
                areaFrac = intersectArea / binArea;
            }
        }
    }
    return result;
}

bool binIntersectsPixel(int32_t row, int32_t col, Polygon_t &pixelPoly, double &areaFrac)
{
    bool result = false;
    double n, s, e, w;
    shape->rowcol2bounds(row, col, n, s, e, w);

#ifdef DEBUG_PRINT2
    // debug plot the bin
    if (enableDebugPrint)
    {
        printf("lat = [%f, %f, %f, %f, %f]\n", s, n, n, s, s);
        printf("lon = [%f, %f, %f, %f, %f]\n", w, w, e, e, w);
        printf("plt.plot(lon, lat)\n\n");
    }
#endif

    Box_t box(Point_t(w, s), Point_t(e, n));
    areaFrac = 0;

    if (!bg::disjoint(box, pixelPoly))
    {
        std::deque<Polygon_t> output;
        if (bg::intersection(box, pixelPoly, output))
        {
            double binArea = (n - s) * (e - w);
            BOOST_FOREACH (Polygon_t const &p, output)
            {
                double intersectArea = bg::area(p);
                if (intersectArea > 0.0)
                {
                    result = true;
                    areaFrac += intersectArea / binArea;
                }
            }
        }
    }
    return result;
}

template <class T>
bool getBinsFromRow(double lat, double lon, T &pixelPoly, std::map<uint64_t, double> &areas)
{
    int32_t row0, col0;
    int32_t col;
    bool result = false;
    double areaFrac;
    uint64_t bin;

    shape->latlon2rowcol(lat, lon, row0, col0);

    // look right
    col = col0;
    while (binIntersectsPixel(row0, col, pixelPoly, areaFrac))
    {
        result = true;
        bin = shape->rowcol2bin(row0, col);
        areas.emplace(bin, areaFrac);
        col++;
        shape->constrainRowCol(row0, col);
        if (col == col0)
            break;
    }

    // look left
    col = col0 - 1;
    while (binIntersectsPixel(row0, col, pixelPoly, areaFrac))
    {
        result = true;
        bin = shape->rowcol2bin(row0, col);
        areas.emplace(bin, areaFrac);
        col--;
        shape->constrainRowCol(row0, col);
        if (col == col0 - 1)
            break;
    }

    return result;
}

void getBins(int32_t ifile, int32_t ipixl, std::map<uint64_t, double> &areas)
{

    areas.clear();

    double lat0 = l2_str[ifile].latitude[ipixl];
    double lon0 = l2_str[ifile].longitude[ipixl];
    int32_t row0;
    double lat;

#ifdef DEBUG_PRINT
    // debug plotting
    // plot center point
    if (enableDebugPrint)
    {
        printf("lat=%f\n", lat0);
        printf("lon=%f\n", lon0);
        printf("plt.plot(lon, lat, \"ro\")\n\n");
    }
#endif

    // use the center of the bin lat for search
    row0 = shape->lat2row(lat0);
    lat0 = shape->row2lat(row0);

    double deltaLat = 180.0 / shape->getNumRows();

    if (input.area_weighting == 1)
    {
        Box_t pixelBox;

        // use the pixel box
        // lat1/lon1 contain pixel delta
        pixelBox.min_corner().set<0>(l2_str[ifile].longitude[ipixl] - l2_str[ifile].lon1[ipixl]);
        pixelBox.min_corner().set<1>(l2_str[ifile].latitude[ipixl] - l2_str[ifile].lat1[ipixl]);
        pixelBox.max_corner().set<0>(l2_str[ifile].longitude[ipixl] + l2_str[ifile].lon1[ipixl]);
        pixelBox.max_corner().set<1>(l2_str[ifile].latitude[ipixl] + l2_str[ifile].lat1[ipixl]);

#ifdef DEBUG_PRINT
        // debug
        // plot pixel box
        if (enableDebugPrint)
        {
            printf("lat = [%f, %f, %f, %f, %f]\n",
                   pixelBox.min_corner().y(), pixelBox.max_corner().y(), pixelBox.max_corner().y(),
                   pixelBox.min_corner().y(), pixelBox.min_corner().y());
            printf("lon = [%f, %f, %f, %f, %f]\n",
                   pixelBox.min_corner().x(), pixelBox.min_corner().x(), pixelBox.max_corner().x(),
                   pixelBox.max_corner().x(), pixelBox.min_corner().x());
            printf("plt.plot(lon, lat)\n\n");
        }
#endif

        // look up
        lat = lat0;
        while (getBinsFromRow(lat, lon0, pixelBox, areas))
        {
            lat += deltaLat;
            if (lat > 90.0)
                break;
        }

        // look down
        lat = lat0 - deltaLat;
        while (getBinsFromRow(lat, lon0, pixelBox, areas))
        {
            lat -= deltaLat;
            if (lat < -90.0)
                break;
        }
    }
    else if (input.area_weighting == 2)
    {
        Box_t pixelBox;

        // use the pixel bounding box
        float latMax = std::max(l2_str[ifile].lat1[ipixl], l2_str[ifile].lat1[ipixl + 1]);
        latMax = std::max(latMax, l2_str[ifile].lat2[ipixl]);
        latMax = std::max(latMax, l2_str[ifile].lat2[ipixl + 1]);

        float latMin = std::min(l2_str[ifile].lat1[ipixl], l2_str[ifile].lat1[ipixl + 1]);
        latMin = std::min(latMin, l2_str[ifile].lat2[ipixl]);
        latMin = std::min(latMin, l2_str[ifile].lat2[ipixl + 1]);

        float lonMax = std::max(l2_str[ifile].lon1[ipixl], l2_str[ifile].lon1[ipixl + 1]);
        lonMax = std::max(lonMax, l2_str[ifile].lon2[ipixl]);
        lonMax = std::max(lonMax, l2_str[ifile].lon2[ipixl + 1]);

        float lonMin = std::min(l2_str[ifile].lon1[ipixl], l2_str[ifile].lon1[ipixl + 1]);
        lonMin = std::min(lonMin, l2_str[ifile].lon2[ipixl]);
        lonMin = std::min(lonMin, l2_str[ifile].lon2[ipixl + 1]);

        if ((lonMax - lonMin) > 180)
        {
            float tmpF = lonMin;
            lonMin = lonMax;
            lonMax = tmpF + 360.0;
        }

        pixelBox.min_corner().set<0>(lonMin);
        pixelBox.min_corner().set<1>(latMin);
        pixelBox.max_corner().set<0>(lonMax);
        pixelBox.max_corner().set<1>(latMax);

#ifdef DEBUG_PRINT
        // debug
        // plot pixel box
        if (enableDebugPrint)
        {
            printf("lat = [%f, %f, %f, %f, %f]\n",
                   pixelBox.min_corner().y(), pixelBox.max_corner().y(), pixelBox.max_corner().y(),
                   pixelBox.min_corner().y(), pixelBox.min_corner().y());
            printf("lon = [%f, %f, %f, %f, %f]\n",
                   pixelBox.min_corner().x(), pixelBox.min_corner().x(), pixelBox.max_corner().x(),
                   pixelBox.max_corner().x(), pixelBox.min_corner().x());
            printf("plt.plot(lon, lat)\n\n");
        }
#endif

        // look up
        lat = lat0;
        while (getBinsFromRow(lat, lon0, pixelBox, areas))
        {
            lat += deltaLat;
            if (lat > 90.0)
                break;
        }

        // look down
        lat = lat0 - deltaLat;
        while (getBinsFromRow(lat, lon0, pixelBox, areas))
        {
            lat -= deltaLat;
            if (lat < -90.0)
                break;
        }
    }
    else
    {
        // use the exact pixel polygon
        Polygon_t pixelPoly;

#ifdef DEBUG_PRINT
        // debug
        // plot pixel polygon
        if (enableDebugPrint)
        {
            printf("lat = [%f, %f, %f, %f, %f]\n",
                   l2_str[ifile].lat1[ipixl],
                   l2_str[ifile].lat2[ipixl],
                   l2_str[ifile].lat2[ipixl + 1],
                   l2_str[ifile].lat1[ipixl + 1],
                   l2_str[ifile].lat1[ipixl]);
            printf("lon = [%f, %f, %f, %f, %f]\n",
                   l2_str[ifile].lon1[ipixl],
                   l2_str[ifile].lon2[ipixl],
                   l2_str[ifile].lon2[ipixl + 1],
                   l2_str[ifile].lon1[ipixl + 1],
                   l2_str[ifile].lon1[ipixl]);
            printf("plt.plot(lon, lat)\n\n");
        }
#endif

        bg::append(pixelPoly.outer(), Point_t(l2_str[ifile].lon1[ipixl], l2_str[ifile].lat1[ipixl]));
        bg::append(pixelPoly.outer(), Point_t(l2_str[ifile].lon2[ipixl], l2_str[ifile].lat2[ipixl]));
        bg::append(pixelPoly.outer(), Point_t(l2_str[ifile].lon2[ipixl + 1], l2_str[ifile].lat2[ipixl + 1]));
        bg::append(pixelPoly.outer(), Point_t(l2_str[ifile].lon1[ipixl + 1], l2_str[ifile].lat1[ipixl + 1]));
        bg::append(pixelPoly.outer(), Point_t(l2_str[ifile].lon1[ipixl], l2_str[ifile].lat1[ipixl]));

        // make sure the polygon is defined properly
        bg::correct(pixelPoly);

        // look up
        lat = lat0;
        while (getBinsFromRow(lat, lon0, pixelPoly, areas))
        {
            lat += deltaLat;
            if (lat > 90.0)
                break;
        }

        // look down
        lat = lat0 - deltaLat;
        while (getBinsFromRow(lat, lon0, pixelPoly, areas))
        {
            lat -= deltaLat;
            if (lat < -90.0)
                break;
        }
    }
}

/**
* @brief checks if lon should be skipped
* @param lon input longitude
* @param side side (west or east) of the dateline should be included
* @param night_flag night
* @param end_day time difference between the start scan time and the end day
* @param beg_day time difference between the start scan time and the beginning day
*/
bool skip_DL(float lon, int side, int night_flag, time_t end_day, time_t beg_day) {
    if (night_flag == 1) {
        if ((side == -1) && (beg_day == -1) && (lon < 0))
            return true;
        if ((side == +1) && (end_day == 0) && (lon > 0))
            return true;
    } else {
        if ((side == -1) && (beg_day <= 0) && (lon < 0))
            return true;
        if ((side == +1) && (end_day >= 0) && (lon > 0))
            return true;
    }
    return false;
}

int main(int argc, char **argv)
{
    int i, j, k;
    int status;
    intn ret_status = 0;

    int32_t index;
    int32_t ifile, jsrow, ipixl;
    uint64_t bin;
    int32_t ibin;
    int32_t nfiles;
    int *fileused;
    int32_t n_active_files;

    int32_t within_flag;

    int32_t n_filled_bins;
    int64_t total_filled_bins = 0;
    int32_t date;

    int16_t brk_scan[MAXNFILES];

    float32 *slat = NULL;
    float32 *elat = NULL;
    float32 dlat;
    double latbin = 0.0;
    double lonbin = 0.0;

    int32_t igroup, ngroup;

    static int32_t *bscan_row[MAXNFILES];
    static int32_t *escan_row[MAXNFILES];
    static unsigned char *scan_in_rowgroup[MAXNFILES];

    int32_t fileid_w;
    int32_t out_grp;

    int64_t *beg;
    int32_t *ext;
    int64_t *binnum_data;
    int32_t i32;
    int32_t iprod;

    time_t diffday_beg, diffday_end;
    int32_t syear, sday, eyear, eday;

    int32_t ntilts;
    int16_t tilt_flags[MTILT_DIMS_2];
    int16_t tilt_ranges[LTILT_DIMS_2][MTILT_DIMS_2];

    uint32_t flagusemask;
    uint32_t required;
    uint32_t flagcheck;

    int32_t proc_day_beg, proc_day_end;

    uint8_t qual_max_allowed;
    std::vector<uint8_t> best_qual;
    double *sum_bin;
    double *sum2_bin;
    double wgt;
    float32 northmost = -90.0, southmost = 90.0, eastmost = -180.0, westmost = 180.0;


    double time_rec = 0;
    time_t tnow;
    int32_t dataday0, dataday1, startdate, enddate;
    double  dbldate;
    struct tm *tmnow;
    static meta_l2Type meta_l2;
    static meta_l3bType meta_l3b;

    static float lonRanges[4];
    lonRanges[0] = 0; // min, max in Western hemisphere
    lonRanges[1] = -180;
    lonRanges[2] = 180; // min, max in Eastern hemisphere
    lonRanges[3] = 0;

    char units[1024];
    std::vector<std::string> l2_prodname;
    std::vector<std::string> l3_prodname;

    /* Function Prototypes */
    int64_t getbinnum(int32_t, int32_t);

    FILE *fp2 = NULL;

    binListStruct_nc *binList32nc = NULL;
    binListStruct64_nc *binList64nc = NULL;

    binIndexStruct_nc *binIndex32nc = NULL;
    binIndexStruct64_nc *binIndex64nc = NULL;

#ifdef MALLINFO
    struct mallinfo minfo;
#endif

    init_rowgroup_cache();

    /* From Fred Patt

        sum(data)    sum(data)*sqrt(n)
   s =  --------- =  -----------------  =  avg(data)*sqrt(n)
         sqrt(n)            n

    */

    setlinebuf(stdout);
    printf("%s %s (%s %s)\n", PROGRAM, VERSION, __DATE__, __TIME__);

    l2bin_input(argc, argv, &input, PROGRAM, VERSION);

    // set dictionary for file names
    std::map<std::string,std::string> output_l3_filenames;
    if(input.output_product_names[0]){
        std::vector<std::string> output_pairs_l2_l3_names;
        boost::split(output_pairs_l2_l3_names,std::string(input.output_product_names),boost::is_any_of(","));
        for(const auto & pair_l2_l3 : output_pairs_l2_l3_names){
            std::vector<std::string> l2_l3_name;
            boost::split(l2_l3_name,pair_l2_l3,boost::is_any_of(":"));
            boost::trim(l2_l3_name[0]);
            boost::trim(l2_l3_name[1]);
            output_l3_filenames[l2_l3_name[0]] = l2_l3_name[1];
        }
    }

    // check to make sure area weighting AND composite_prod are not both set
    if(input.composite_prod[0] != 0 && input.area_weighting != 0) {
        printf("-E- Arguments \"composite_prod\" and \"area_weighting\" are incompatable.  Do not set both.\n");
        exit(EXIT_FAILURE);
    }

    switch (input.area_weighting)
    {
    case 0:
        enableL2PixelArea(L2PixelOff); // no corner calc
        break;
    case 1:
        enableL2PixelArea(L2PixelDelta); // pixel deltas
        break;
    default:
        enableL2PixelArea(L2PixelCorner); // pixel corners
        break;
    }

    // Get times for each file
    nfiles = input.files.size();
    printf("%d input files\n", nfiles);

    bool flags_l2_use = set_l2_flags_use(input.flaguse);
    std::string products_requested_l3_temp,products_requested_l3;
    std::string products_requested;
    for (ifile = 0; ifile < nfiles; ifile++)
    {
        std::string product_list_temp;
        std::string product_list;
        products_requested = input.l3bprod;
        if (products_requested == "ALL")
        {
            products_requested = "";
            find_variables_geo_physical(input.files[ifile],products_requested);
        }

        // check for quality products and composite products
        if (input.composite_prod[0])
        {
            if (products_requested.find(input.composite_prod) == std::string::npos)
            {
                products_requested.append(",");
                products_requested.append(input.composite_prod);
            }
        }
        if (input.qual_prod[0])
        {
            if (products_requested.find(input.qual_prod) == std::string::npos)
            {
                products_requested.append(",");
                products_requested.append(input.qual_prod);
            }
        }

        // set the list of requested products
        std::vector<std::string> products_requested_separated;
        set_mapper(input.files[ifile], products_requested, products_requested_separated,input.output_wavelengths);

        // sanitize the input product list: remove duplicates
        std::unordered_set<std::string> products_l2_unique;
        for (size_t i = 0; i < products_requested_separated.size(); i++)
        {
            const auto &prod = products_requested_separated.at(i);
            const auto name_to_pass = prod; // return_l2_name(prod);
            if (products_l2_unique.count(name_to_pass) > 0)
                continue;
            if (!product_list_temp.empty())
                product_list_temp += "," + name_to_pass;
            else
                product_list_temp += name_to_pass; //
            if(name_to_pass!=input.qual_prod)
            {
                if (!product_list.empty())
                product_list += "," + name_to_pass;
                    else
                product_list += name_to_pass;
            }//

            products_l2_unique.insert(name_to_pass);
        }
        products_requested_l3_temp = product_list_temp;
        products_requested_l3 = product_list;
        status = openL2(input.files[ifile].c_str(), products_requested_l3_temp.c_str(), &l2_str[ifile]);
        // l2_stimes[ifile] = (time_t)yds2unix(l2_str[ifile].syear,
        //                                     l2_str[ifile].sday,
        //                                     l2_str[ifile].smsec / 1000); // 86400;
        closeL2(&l2_str[ifile], ifile);
    }

    if ((fileused = (int *)calloc(nfiles, sizeof(int))) == NULL)
    {
        printf("-E- Problem allocating memory for fileused element of metadata structure\n");
        exit(EXIT_FAILURE);
    }

    proc_day_beg = input.sday;
    proc_day_end = input.eday;
    syear = (int32_t)input.sday / 1000.;
    sday = input.sday - 1000 * syear;
    startdate = (int32_t)(yds2unix(syear, sday, 0) / 86400);
    eyear = (int32_t)input.eday / 1000.;
    eday = input.eday - 1000 * eyear;
    enddate = (int32_t)(yds2unix(eyear, eday, 0) / 86400);

    resolve2binRows(input.resolve, nrows, meta_l3b.resolution);

    if (nrows == -1)
    {
        printf("Grid resolution not defined.\n");
        exit(EXIT_FAILURE);
    }
    dlat = 180. / nrows;
    bool is64bit = (nrows > 2160 * 16);

    qual_max_allowed = input.qual_max;

    printf("Resolution: %s\n", input.resolve);
    printf("Max Qual Allowed: %d\n", input.qual_max);

    /* Setup flag mask */
    /* --------------- */
    if (flags_l2_use)
    {
        strcpy(buf, l2_str[0].flagnames);
        setupflags(buf, input.flaguse, &flagusemask, &required, &status,l2_str[0].l2_bits);
    }

    // Parse L3 Product list
    std::string delim1 = ", "; // product delimiters
    std::string l3bprod = products_requested_l3;
    boost::trim_if(l3bprod, boost::is_any_of(delim1));
    std::vector<std::string> prodparam;
    boost::algorithm::split(prodparam, l3bprod,
                            boost::is_any_of(delim1));

    // set min, max and thirdim
    set_prodname_3d_to_l2(prodparam, l2_str[0], l2_prodname, l3_prodname, thirdDimId, min_value, max_value);
    l3b_nprod = l3_prodname.size();
    if (l3b_nprod > max_l3b_products)
    {
        printf("-E- Number of output products is greater than %d\n", max_l3b_products);
        exit(EXIT_FAILURE);
    }

    /* Initialize bscan_row, escan_row, numer, denom */
    /* --------------------------------------------- */
    for (i = 0; i < MAXNFILES; i++)
    {
        bscan_row[i] = NULL;
        escan_row[i] = NULL;
        numer[i] = NULL;
        denom[i] = NULL;
        scan_in_rowgroup[i] = NULL;
    }

    /* Check each L2 file for required products */
    /* ---------------------------------------- */
    for (ifile = 0; ifile < nfiles; ifile++)
    {

        numer[ifile] = (int16_t *)calloc(l3b_nprod, sizeof(int16_t));
        denom[ifile] = (int16_t *)calloc(l3b_nprod, sizeof(int16_t));

        // Check whether L3 products exist in L2
        for (iprod = 0; iprod < l3b_nprod; iprod++)
        {

            std::vector<std::string> operands;
            boost::algorithm::split(operands, l2_prodname[iprod],
                                    boost::is_any_of("/"));

            // find numerator, or solo product
            index = get_l2prod_index(l2_str[ifile], operands[0].c_str());
            if (index < 0)
            {
                printf("L3 product: \"%s\" not found in L2 file \"%s\".\n",
                       operands[0].c_str(), l2_str[ifile].filename);
                exit(EXIT_FAILURE);
            }
            else
                numer[ifile][iprod] = index;

            // find denominator product (as needed)
            denom[ifile][iprod] = -1;
            if (operands.size() > 1)
            {
                boost::replace_all(l3_prodname[iprod], "/", "_");
                index = get_l2prod_index(l2_str[ifile], operands[1].c_str());
                if (index < 0)
                {
                    printf("L3 product: \"%s\" not found in L2 file \"%s\".\n",
                           operands[1].c_str(), l2_str[ifile].filename);
                    exit(EXIT_FAILURE);
                }
                else
                    denom[ifile][iprod] = index;
            }

        } // iprod loop

        // Check whether Quality product exists in L2
        if (input.qual_prod[0] != 0)
        {
            index = get_l2prod_index(l2_str[ifile], input.qual_prod);
            if (index < 0)
            {
                printf("Quality product: \"%s\" not found in L2 file \"%s\".\n.See line %d in %s\n",
                       input.qual_prod, l2_str[ifile].filename, __LINE__, __FILE__);
                exit(EXIT_FAILURE);
            }
            else
                qual_prod_index[ifile] = index;
        }

        // Check whether Composite product exists in L2
        if (input.composite_prod[0] != 0)
        {
            index = get_l2prod_index(l2_str[ifile], input.composite_prod);
            if (index < 0)
            {
                printf("Composite product: \"%s\" not found in L2 file \"%s\".\n",
                       input.composite_prod, l2_str[ifile].filename);
                exit(EXIT_FAILURE);
            }
            else
                composite_prod_index[ifile] = index;
        }

    } // ifile loop

    /* Check whether composite product exists in L3 product list */
    /* --------------------------------------------------------- */
    if (input.composite_prod[0] != 0)
    {
        for (iprod = 0; iprod < l3b_nprod; iprod++)
        {
            if (l3_prodname[iprod] == input.composite_prod)
            {
                composite_l3prod_index = iprod;
                break;
            }
        }
        if (composite_l3prod_index == -1)
        {
            printf("Composite product: \"%s\" not found in L3 product list.\n",
                   input.composite_prod);
            exit(EXIT_FAILURE);
        }
    }

    /* Find begin and end scan latitudes for each swath row */
    /* ---------------------------------------------------- */
    for (ifile = 0; ifile < nfiles; ifile++)
    {
        size_t nrec = l2_str[ifile].nrec;
        slat = (float32 *)calloc(l2_str[ifile].nrec, sizeof(float32));
        elat = (float32 *)calloc(l2_str[ifile].nrec, sizeof(float32));
        bscan_row[ifile] = (int32_t *)calloc(l2_str[ifile].nrec, sizeof(int32_t));
        escan_row[ifile] = (int32_t *)calloc(l2_str[ifile].nrec, sizeof(int32_t));
        std::string instrument, platform;
        try {
            netCDF::NcFile nc_input(input.files[ifile], netCDF::NcFile::read);
            Geospatialbounds geo_bounds(nc_input);
            instrument = geo_bounds.get_instrument();
            platform = geo_bounds.get_platform();
            std::copy(geo_bounds.get_slat(),geo_bounds.get_slat() + nrec,slat);
            std::copy(geo_bounds.get_elat(),geo_bounds.get_elat() + nrec,elat);
            nc_input.close();
        }
        catch (netCDF::exceptions::NcException const &e)
        {
            std::cout << e.what() << std::endl;
            exit(EXIT_FAILURE);
        }
        catch (std::exception const &e)
        {
            std::cerr << "-X- " << __FILE__ " line " << __LINE__ << ": "
                      << e.what() << std::endl;
            exit(EXIT_FAILURE);
        }

        /* Get minimum/maximu value from product.xml */
        /* ---------------------------------- */
        sensorID[ifile] = instrumentPlatform2SensorId(instrument.c_str(),
                                                      platform.c_str());
        /* Note: bscan > escan */

        for (jsrow = 0; jsrow < l2_str[ifile].nrec; jsrow++)
        {
            escan_row[ifile][jsrow] = (int32_t)((90 + elat[jsrow]) / dlat);
            bscan_row[ifile][jsrow] = (int32_t)((90 + slat[jsrow]) / dlat);

            if (escan_row[ifile][jsrow] > bscan_row[ifile][jsrow])
            {
                k = escan_row[ifile][jsrow];
                escan_row[ifile][jsrow] = bscan_row[ifile][jsrow];
                bscan_row[ifile][jsrow] = k;
            }
            escan_row[ifile][jsrow] -= 10;
            bscan_row[ifile][jsrow] += 10;
        }

        free(slat);
        free(elat);

    } /* ifile loop */

    /* Find begin & end scans for each input file */
    /* ------------------------------------------ */
    n_active_files = nfiles;
    for (ifile = 0; ifile < nfiles; ifile++)
    {

        /* Determine brk_scan value */
        /* ------------------------ */
        brk_scan[ifile] = 0;

        /* Regional Product */
        /* ---------------- */
        if (strcasecmp(input.prodtype, "regional") == 0)
        {
            // printf("%s   brk:%5d  %5d %3d %6d\n",
            //        buf, brk_scan[ifile],
            //        l2_str[ifile].nrec, l2_str[ifile].sday,
            //        l2_str[ifile].smsec / 1000);
            continue;
        }
        // get data days

        elat = (float32 *)calloc(l2_str[ifile].nrec, sizeof(float32));
        slat = (float32 *)calloc(l2_str[ifile].nrec, sizeof(float32));
        size_t nrec = l2_str[ifile].nrec;
        try {
            netCDF::NcFile nc_input(input.files[ifile], netCDF::NcFile::read);
            Geospatialbounds geo_bounds(nc_input);
            std::copy(geo_bounds.get_slat(),geo_bounds.get_slat() + nrec,slat);
            std::copy(geo_bounds.get_elat(),geo_bounds.get_elat() + nrec,elat);
            std::pair<int32_t,int32_t> dates_0_1 = get_datadays(nc_input,input.deltaeqcross,input.night);
            dataday0 = dates_0_1.first;
            dataday1 = dates_0_1.second;
            nc_input.close();
        }
        catch (netCDF::exceptions::NcException const &e)
        {
            std::cout << e.what();
            exit(EXIT_FAILURE);
        }
        catch (std::exception const &e)
        {
            std::cerr << "-X- " << __FILE__ " line " << __LINE__ << ": "
                      << e.what() << std::endl;
            exit(EXIT_FAILURE);
        }

        date = (time_t)yds2unix(l2_str[ifile].syear,
                                l2_str[ifile].sday,
                                l2_str[ifile].smsec / 1000) /
               86400;
        diffday_beg = date - startdate;
        diffday_end = date - enddate;
        if (dataday1 == dataday0)
        {
            if (dataday0 < startdate || dataday1 > enddate)
                brk_scan[ifile] = -9999;
            else
                brk_scan[ifile] = 0; // -9999;
        }
        else
        {
            if (dataday1 < startdate) // startdate is dataday conversion of input sday
                brk_scan[ifile] = -9999;
            else if (dataday0 > enddate)
            { // enddate is dataday conversion of input eday
                brk_scan[ifile] = -9999;
            }
            else
            {
                if (dataday1 > enddate)
                    brk_scan[ifile] = 1;
                else
                    brk_scan[ifile] = -1;
            }
        }
        if (brk_scan[ifile] == -9999)
            n_active_files--;
        free(elat);
        free(slat);
    } /* ifile loop */


    shape = new l3::L3ShapeIsine(nrows);

    /* Compute basebin array (Starting bin of each row [1-based]) */
    /* ---------------------------------------------------------- */
    basebin = (int64_t *)calloc(nrows + 1, sizeof(int64_t));
    for (i = 0; i <= nrows; i++)
    {
        basebin[i] = shape->getBaseBin(i);
    }
    basebin[nrows] = shape->getNumBins() + 1;

    printf("total number of bins: %ld\n", (long int)basebin[nrows] - 1);

    /*
     * Create output file
     */
    strcpy(buf, input.ofile);
    status = nc_create(buf, NC_NETCDF4, &fileid_w);
    check_err(status, __LINE__, __FILE__);

    status = nc_def_grp(fileid_w, "level-3_binned_data", &out_grp);
    check_err(status, __LINE__, __FILE__);

    idDS ds_id;
    ds_id.fid = fileid_w;
    ds_id.sid = -1;
    ds_id.fftype = DS_NCDF; // FMT_L2NCDF

    if (is64bit)
        status = defineBinList64_nc(input.deflate, out_grp);
    else
        status = defineBinList_nc(input.deflate, out_grp);
    if (status)
    {
        printf("-E- %s:%d Could not define binList variable in output file\n",
               __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }

    char **prodnames = (char **)malloc(l3_prodname.size() * sizeof(char *));
    for (size_t size = 0; size < l3_prodname.size(); size++)
    {
        // look for user requested product names
        if(output_l3_filenames.find(l3_prodname[size]) != output_l3_filenames.end())
            prodnames[size] = strdup(output_l3_filenames[l3_prodname[size]].c_str());
        else
            prodnames[size] = strdup(l3_prodname[size].c_str());
    }
    status = defineBinData_nc(input.deflate, out_grp, l3b_nprod, prodnames);
    // set  wave float attribute here
    for (size_t size = 0; size < l3_prodname.size(); size++) {
        float wave = BAD_FLT;
        wave = l3_attr(l3_prodname[size]);
        if(wave !=BAD_FLT) {
            int varid;
            nc_inq_varid(out_grp, prodnames[size] , &varid);
            nc_put_att_float(out_grp,varid,"wavelength",NC_FLOAT,1,&wave);
        }
    }

    if (status)
    {
        printf("-E- %s:%d Could not define binData variable in output file\n",
               __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    for (size_t i = 0; i < l3_prodname.size(); i++)
    {
        free(prodnames[i]);
    }
    free(prodnames);

    if (is64bit)
        status = defineBinIndex64_nc(input.deflate, out_grp);
    else
        status = defineBinIndex_nc(input.deflate, out_grp);
    if (status)
    {
        printf("-E- %s:%d Could not define binIndex variable in output file\n",
               __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }

    if (input.qual_prod[0] != 0)
    {
        status = defineQuality_nc(input.deflate, out_grp);
        if (status)
        {
            printf("-E- %s:%d Could not define quality variable in output file\n",
                   __FILE__, __LINE__);
            exit(EXIT_FAILURE);
        }
    }

    /* Allocate Arrays for Bin Index */
    /* ----------------------------- */
    beg = (int64_t *)calloc(nrows, sizeof(int64_t));
    ext = (int32_t *)calloc(nrows, sizeof(int32_t));

    if (is64bit)
        binIndex64nc = (binIndexStruct64_nc *)calloc(nrows, sizeof(binIndexStruct64_nc));
    else
        binIndex32nc = (binIndexStruct_nc *)calloc(nrows, sizeof(binIndexStruct_nc));

    /* Initialize bin_indx array */
    /* ------------------------- */
    for (i = 0; i < nrows; i++)
    {

        i32 = i;

        if (i32 < 0 || i32 >= nrows)
        {
            printf("%d %d\n", i, nrows);
            exit(-1);
        }

        if (is64bit)
        {
            binIndex64nc[i].begin = beg[i32];
            binIndex64nc[i].extent = ext[i32];
            binIndex64nc[i].max = shape->getNumCols(i32);
        }
        else
        {
            binIndex32nc[i].begin = beg[i32];
            binIndex32nc[i].extent = ext[i32];
            binIndex32nc[i].max = shape->getNumCols(i32);
        }
    }

    // Row Group
    n_rows_in_group = input.rowgroup;
    if (n_rows_in_group <= 0)
    {
        printf("row_group not defined, using 270.\n");
        n_rows_in_group = 270;
    }
    if (input.verbose)
        printf("%d %d %d\n", proc_day_beg, proc_day_end, n_rows_in_group);

    /* Find row_group that divides nrows */
    /* --------------------------------- */
    for (i = nrows; i > 0; i--)
    {
        if ((nrows % i) == 0)
        {
            if (i <= n_rows_in_group)
            {
                n_rows_in_group = i;
                break;
            }
        }
    }
    if (input.rowgroup != n_rows_in_group)
    {
        printf("Input row_group: %d   Actual row_group: %d\n",
               input.rowgroup, n_rows_in_group);
    }
    ngroup = nrows / n_rows_in_group;

    /* Process each group of bin rows (Main Loop) */
    /* ========================================== */
    int64_t total_count = 0;
    for (krow = 0, igroup = 0; igroup < ngroup; igroup++)
    {

        if (shape->row2lat(krow + n_rows_in_group) < input.latsouth)
        {
            krow += n_rows_in_group;
            continue;
        }
        if (shape->row2lat(krow) > input.latnorth)
        {
            krow += n_rows_in_group;
            continue;
        }

        /* Print info on rowgroup */
        /* ---------------------- */
        time(&tnow);
        tmnow = localtime(&tnow);
        printf("krow:%6d out of %6d  (%6.2f to %6.2f) ",
               krow, nrows,
               ((float32)(krow) / nrows) * 180 - 90,
               ((float32)(krow + n_rows_in_group) / nrows) * 180 - 90);
        printf("%s", asctime(tmnow));

        n_bins_in_group = basebin[krow + n_rows_in_group] - basebin[krow];
        within_flag = 0;

        /* Determine relevant swath rows for this bin row group for each file */
        /* ------------------------------------------------------------------ */
        for (ifile = 0; ifile < nfiles; ifile++)
        {

            /* add an extra 0 to the end of scan_in_rowgroup so the caching
             * code never reads past the end of the file */
            scan_in_rowgroup[ifile] = (unsigned char *)
                calloc(l2_str[ifile].nrec + 1, sizeof(unsigned char));

            for (jsrow = 0; jsrow < l2_str[ifile].nrec; jsrow++)
            {
                scan_in_rowgroup[ifile][jsrow] = 1;
                if (bscan_row[ifile][jsrow] < krow ||
                    escan_row[ifile][jsrow] >= (krow + n_rows_in_group - 1))
                {
                    scan_in_rowgroup[ifile][jsrow] = 255;
                }
            } /* jsrow loop */

            /* Determine if within bin row group */
            /* --------------------------------- */
            if (!within_flag)
            {
                for (jsrow = 0; jsrow < l2_str[ifile].nrec; jsrow++)
                {
                    if (scan_in_rowgroup[ifile][jsrow] == 1)
                    {
                        within_flag = 1;
                        break;
                    }
                } /* scan row loop */
            }

        } /* ifile loop */

        /* If no swath rows within group then continue to next group */
        /* --------------------------------------------------------- */
        if (within_flag == 0)
        {
            for (ifile = 0; ifile < nfiles; ifile++)
            {
                if (scan_in_rowgroup[ifile] != NULL)
                {
                    free(scan_in_rowgroup[ifile]);
                    scan_in_rowgroup[ifile] = NULL;
                }
            }
            krow += n_rows_in_group;
            continue;
        }

        /* Allocate # pixels in bin, bin_flag, tilt, qual, & nscenes arrays */
        /* ---------------------------------------------------------------- */
        n_filled_bins = 0;
        bin_flag = (int32_t *)calloc(n_bins_in_group, sizeof(int32_t));
        tilt = (int16_t *)calloc(n_bins_in_group, sizeof(int16_t));
        qual = (int16_t *)calloc(n_bins_in_group, sizeof(int16_t));
        nscenes = (int16_t *)calloc(n_bins_in_group, sizeof(int16_t));
        lastfile = (int16_t *)calloc(n_bins_in_group, sizeof(int16_t));

        for (i = 0; i < n_bins_in_group; i++)
        {
            tilt[i] = -1;
            qual[i] = 3;
            lastfile[i] = -1;
        }

        /* Allocate bin accumulator & data value arrays */
        /* -------------------------------------------- */
        nobs = (int16_t *)calloc(n_bins_in_group, sizeof(int16_t));
        allocated_space = (int16_t *)calloc(n_bins_in_group, sizeof(int16_t));
        data_values = (float32 **)calloc(n_bins_in_group, sizeof(float32 *));
        data_areas = (double **)calloc(n_bins_in_group, sizeof(double *));
        file_index = (int16_t **)calloc(n_bins_in_group, sizeof(int16_t *));
        data_quality = (uint8_t **)calloc(n_bins_in_group, sizeof(uint8_t *));
        time_value = (float64 **)calloc(n_bins_in_group, sizeof(float64 *));

        /* Initialize bin counters */
        /* ----------------------- */
        n_allocperbin =
            n_active_files * l2_str[0].nrec * l2_str[0].nsamp / 50000000;

        if (n_allocperbin < 2)
            n_allocperbin = 2;
        if (n_allocperbin > MAXALLOCPERBIN)
            n_allocperbin = MAXALLOCPERBIN;

        if (input.verbose)
            printf("%-20s:%8d\n\n", "# allocated per bin", n_allocperbin);

        for (i = 0; i < n_bins_in_group; i++)
        {
            nobs[i] = 0;
            allocated_space[i] = 0;
            lastfile[i] = -1;
        }

        /* Read L2 files and fill data_values (L3b) array */
        /* ++++++++++++++++++++++++++++++++++++++++++++++ */
        for (ifile = 0; ifile < nfiles; ifile++)
        {

            free_rowgroup_cache();

            /* if "early" or "late" input file then skip */
            /* ----------------------------------------- */
            if (brk_scan[ifile] == -9999)
                continue;

            /* if no scans in rowgroup for this file then skip */
            /* ----------------------------------------------- */
            for (jsrow = 0; jsrow < l2_str[ifile].nrec; jsrow++)
            {
                if (scan_in_rowgroup[ifile][jsrow] == 1)
                {
                    break;
                }
            }
            if (jsrow == l2_str[ifile].nrec)
            {
                continue;
            }

            reopenL2(ifile, &l2_str[ifile]);

            /* Get tilt flags & ranges */
            /* ----------------------- */
            ntilts = l2_str[ifile].ntilts;
            for (i = 0; i < ntilts; i++)
            {
                tilt_flags[i] = l2_str[ifile].tilt_flags[i];
                tilt_ranges[0][i] = l2_str[ifile].tilt_ranges[0][i];
                tilt_ranges[1][i] = l2_str[ifile].tilt_ranges[1][i];
            }

            /* Get date stuff */
            /* -------------- */
            date = (time_t)yds2unix(l2_str[ifile].syear, l2_str[ifile].sday,
                                    l2_str[ifile].smsec / 1000) /
                   86400;
            diffday_beg = date - startdate;
            diffday_end = date - enddate;

            /* Loop over swath rows */
            /* ^^^^^^^^^^^^^^^^^^^^ */
            for (jsrow = 0; jsrow < l2_str[ifile].nrec; jsrow++)
            {

                /* if swath row not within group then continue */
                /* ------------------------------------------- */
                if (scan_in_rowgroup[ifile][jsrow] != 1)
                    continue;

                /* Read swath record from L2 */
                /* ------------------------- */
                 auto start = std::chrono::high_resolution_clock::now();
                status = readL2(&l2_str[ifile], ifile, jsrow, -1,
                                scan_in_rowgroup[ifile]);
                auto stop = std::chrono::high_resolution_clock::now();
                auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
                auto count_ = duration.count();
                total_count += count_;
                if (status == 5)
                    continue;
                int scan_crosses_dateline = 0;
                float scan_lon = BAD_FLT;
                if (brk_scan[ifile] != 0) {
                    int npixls = l2_str[ifile].nsamp - 1;
                    float *ptr_lon = l2_str[ifile].longitude;
                    float slon = BAD_FLT;
                    float elon = BAD_FLT;

                    for (int i_p = 0; i_p <= npixls; i_p++) {
                        if (std::abs(slon) <= 180.0)
                            break;
                        slon = ptr_lon[i_p];
                    }
                    for (int i_p = npixls; i_p >= 0; i_p--) {
                        if (std::abs(elon) <= 180.0)
                            break;
                        elon = ptr_lon[i_p];
                    }
                    if (std::abs(elon) > 180.0 || std::abs(slon) > 180.0)
                        continue;
                    scan_lon = elon;
                    if (slon * elon < 0)
                        scan_crosses_dateline = 1;
                }
                if (scan_crosses_dateline == 0 && (strcasecmp(input.prodtype, "regional") != 0)) {
                    if (skip_DL(scan_lon, brk_scan[ifile], input.night, diffday_end, diffday_beg))
                        continue;
                }

                /* Check tilt state */
                /* ---------------- */
                for (i = 0; i < ntilts; i++)
                {
                    if ((jsrow + 1) <= tilt_ranges[1][i])
                    {
                        tiltstate = (tilt_flags[i] & 0xFF);
                        break;
                    }
                }

                if (l2_str[ifile].nsamp == 0)
                    continue;

                if ((jsrow % 100) == 0 && input.verbose)
                {
                    printf("ifile:%4d  jsrow:%6d  nsamp:%8d\n", ifile, jsrow, l2_str[ifile].nsamp);
                }

                // save the tai93 date for this line, so is can be saved in the bins
                dbldate = yds2unix(l2_str[ifile].year, l2_str[ifile].day, l2_str[ifile].msec / 1000.0);
                time_tai93 = unix_to_tai93(dbldate);

                /* ##### Loop over L2 pixels ##### */
                /* ------------------------------- */
                for (ipixl = 0; ipixl < l2_str[ifile].nsamp; ipixl++)
                {

#ifdef DEBUG_PRINT
                    if (ipixl >= 510 && ipixl < 512)
                        enableDebugPrint = true;
                    else
                        enableDebugPrint = false;
#endif

                    /* if bin flagged then continue */
                    /* ---------------------------- */
                    if (flags_l2_use)
                    {
                        flagcheck = (l2_str[ifile].l2_flags[ipixl]);
                        if ((flagcheck & flagusemask) != 0)
                            continue;
                        if ((flagcheck & required) != required)
                            continue;
                    }

                    /* Check for dateline crossing */
                    /* --------------------------- */
                    if (scan_crosses_dateline == 1) {
                        if (skip_DL(l2_str[ifile].longitude[ipixl], brk_scan[ifile], input.night, diffday_end,
                                    diffday_beg))
                            continue;
                    }
                    /* Check for bad value in any of the products */
                    /* ------------------------------------------ */
                    int bad_value = 0;
                    for (iprod = 0; iprod < l3b_nprod; iprod++) {
                        int32_t l2_iprod = numer[ifile][iprod];
                        for (int i = 0; i < l2_str[ifile].thirdDim[l2_iprod]; i++) {
                            f32 =
                                l2_str[ifile].l2_data[l2_iprod][ipixl * l2_str[ifile].thirdDim[l2_iprod] + i];
                            if (f32 == BAD_FLT) {
                                bad_value = 1;
                                break;
                            }
                        }
                        if (bad_value)
                            break;
                    }
                    if (bad_value == 1)
                        continue;

                    /* Check if within longitude boundaries */
                    /* ------------------------------------ */
                    if (input.lonwest != 0.0 || input.loneast != 0.0)
                    {
                        if (l2_str[ifile].longitude[ipixl] < input.lonwest)
                            continue;
                        if (l2_str[ifile].longitude[ipixl] > input.loneast)
                            continue;
                    }
                    /* Check if within latitude boundaries */
                    /* ------------------------------------ */
                    if (l2_str[ifile].latitude[ipixl] < input.latsouth)
                        continue;
                    if (l2_str[ifile].latitude[ipixl] > input.latnorth)
                        continue;
                    if(std::abs(l2_str[ifile].latitude[ipixl]) > 90)
                        continue;
                    if(std::abs(l2_str[ifile].longitude[ipixl]) > 180)
                        continue;
                    if (input.area_weighting) {
                        if (std::abs(l2_str[ifile].lat1[ipixl]) > 90)
                            continue;
                        if (std::abs(l2_str[ifile].lon1[ipixl]) > 180)
                            continue;
                        if (std::abs(l2_str[ifile].lat1[ipixl + 1]) > 90)
                            continue;
                        if (std::abs(l2_str[ifile].lon1[ipixl + 1]) > 180)
                            continue;
                    }
                    if (input.area_weighting >= 2) {
                        if (std::abs(l2_str[ifile].lat2[ipixl]) > 90)
                            continue;
                        if (std::abs(l2_str[ifile].lon2[ipixl]) > 180)
                            continue;
                        if (std::abs(l2_str[ifile].lat2[ipixl + 1]) > 90)
                            continue;
                        if (std::abs(l2_str[ifile].lon2[ipixl + 1]) > 180)
                            continue;
                    }
                    if (input.area_weighting)
                    {
                        /* Get map of bin numbers and intersection area that intersect pixel */
                        std::map<uint64_t, double> areas;
                        getBins(ifile, ipixl, areas);
                        for (auto const &val : areas)
                        {
                            addPixelToBin(ifile, ipixl, val.first, flags_l2_use, val.second);
                        }
                    }
                    else
                    {
                        /* Get Bin Number for Pixel */
                        /* ------------------------ */
                        bin = getbinnum(ifile, ipixl); // bin is 1-based
                        addPixelToBin(ifile, ipixl, bin, flags_l2_use);
                    }

                } /* ##### i (ipixl) loop ##### */

            } /* ^^^^^^^^^^ jsrow loop ^^^^^^^^^^ */

            closeL2(&l2_str[ifile], ifile);

#ifdef MALLINFO
            if (input.meminfo)
            {
                int32_t total_alloc_space;
                /*      malloc_stats();*/
                minfo = mallinfo();
                total_alloc_space = 0;
                for (i = 0; i < n_bins_in_group; i++)
                {
                    total_alloc_space += allocated_space[i];
                }
                printf("Used space: %10d\n", minfo.uordblks);
                printf("Allo space: %10d\n", total_alloc_space * (2 + l3b_nprod * 4));
            }
#endif

        } /* ++++++++++ ifile loop ++++++++++ */

        if (input.verbose)
        {
            time(&tnow);
            tmnow = localtime(&tnow);
            printf("krow:%5d After data_value fill: %s\n", krow, asctime(tmnow));
        }

        /* Compute Total # of filled bins */
        /* ------------------------------ */
        for (ibin = 0; ibin < n_bins_in_group; ibin++)
        {
            if (nobs[ibin] > 0 && nobs[ibin] < input.minobs)
                nobs[ibin] = 0;

            if (nobs[ibin] != 0)
                n_filled_bins++;
        } /* ibin loop */

        best_qual = std::vector<uint8_t>(n_bins_in_group, 255);
        // memset(best_qual, 255, n_bins_in_group * sizeof(uint8_t));

        /* ********** If filled bins ********** */
        /* ------------------------------------ */
        if (n_filled_bins > 0)
        {

            /* Fill "Bin List" vdata array */
            /* --------------------------- */
            if (is64bit)
                binList64nc = (binListStruct64_nc *)
                    calloc(n_filled_bins, sizeof(binListStruct64_nc));
            else
                binList32nc = (binListStruct_nc *)
                    calloc(n_filled_bins, sizeof(binListStruct_nc));

            i = 0;
            for (ibin = 0; ibin < n_bins_in_group; ibin++)
            {

                /* Adjust for bins with "bad" quality values */
                /* ----------------------------------------- */
                if (input.qual_prod[0] != 0 && nobs[ibin] > 0)
                {
                    best_qual[ibin] = 255;
                    for (j = 0; j < nobs[ibin]; j++)
                        if (data_quality[ibin][j] < best_qual[ibin])
                            best_qual[ibin] = data_quality[ibin][j];

                    k = 0;
                    for (j = 0; j < nobs[ibin]; j++)
                    {
                        if ((data_quality[ibin][j] <= best_qual[ibin]) &&
                            (data_quality[ibin][j] <= qual_max_allowed))
                        {
                            if (k < j)
                            {
                                data_areas[ibin][k] = data_areas[ibin][j];
                                for (iprod = 0; iprod < l3b_nprod; iprod++)
                                {
                                    data_values[ibin][k * l3b_nprod + iprod] =
                                        data_values[ibin][j * l3b_nprod + iprod];
                                }
                            }
                            k++;
                        }
                    }
                    nobs[ibin] = k;

                    if (nobs[ibin] == 0)
                        n_filled_bins--;
                }

                if (nobs[ibin] != 0)
                {
                    bin = (uint64_t)ibin + basebin[krow];

                    if (is64bit)
                    {
                        binList64nc[i].binnum = bin;
                        binList64nc[i].nobs = nobs[ibin];
                        binList64nc[i].nscenes = nscenes[ibin];
                    }
                    else
                    {
                        binList32nc[i].binnum = bin;
                        binList32nc[i].nobs = nobs[ibin];
                        binList32nc[i].nscenes = nscenes[ibin];
                    }

                    /* weights {=sqrt(# of L2 files in given bin)} */
                    /* ------------------------------------------- */
                    wgt = 0.0;
                    for (ifile = 0; ifile <= nfiles; ifile++)
                    {
                        double area = 0.0;
                        for (j = 0; j < nobs[ibin]; j++)
                        {
                            if (file_index[ibin][j] == ifile)
                                area += data_areas[ibin][j];
                        }
                        wgt += sqrt(area);
                    }

                    time_rec = 0.0;
                    for (j = 0; j < nobs[ibin]; j++)
                        time_rec += time_value[ibin][j];

                    if (is64bit)
                    {
                        binList64nc[i].weights = wgt;
                        binList64nc[i].time_rec = time_rec;
                    }
                    else
                    {
                        binList32nc[i].weights = wgt;
                        binList32nc[i].time_rec = time_rec;
                    }

                    i++;

                    /* Update Max/Min Lon/Lat */
                    /* ---------------------- */
                    shape->bin2latlon(bin, latbin, lonbin);

                    if (latbin > northmost)
                        northmost = latbin;
                    if (latbin < southmost)
                        southmost = latbin;

                    float minNeg = lonRanges[0]; // min, max in Western hemisphere
                    float maxNeg = lonRanges[1];
                    float minPos = lonRanges[2]; // min, max in Eastern hemisphere
                    float maxPos = lonRanges[3];
                    if (lonbin < 0)
                    { // Western hemisphere
                        minNeg = fmin(minNeg, lonbin);
                        maxNeg = fmax(maxNeg, lonbin);
                    }
                    else
                    { // Eastern hemisphere
                        minPos = fmin(minPos, lonbin);
                        maxPos = fmax(maxPos, lonbin);
                    }
                    // Adjust east and west granule bounding coordinates
                    float max_angle = 90.0;
                    int8_t lon000 = ((minPos - maxNeg) < max_angle);
                    int8_t lon180 = ((maxPos - minNeg) > (360.0 - max_angle));
                    int8_t polar = (lon000 && lon180) ||
                                   (90.0 - fmax(northmost, -1 * southmost) <= FLT_EPSILON);

                    if (minNeg >= maxNeg)
                    { // Entirely in Eastern hemisphere
                        eastmost = maxPos;
                        westmost = minPos;
                    }
                    else if (minPos >= maxPos)
                    { // Entirely in Western hemisphere
                        eastmost = maxNeg;
                        westmost = minNeg;
                    }
                    else if (polar)
                    { // Polar granule
                        eastmost = 180;
                        westmost = -180;
                        if (northmost > 0)
                            northmost = 90; // North pole
                        else
                            southmost = -90; // South pole
                    }
                    else if (lon000)
                    { // Granule crosses Longitude 0
                        eastmost = maxPos;
                        westmost = minNeg;
                    }
                    else if (lon180)
                    { // Granule crosses Longitude 180
                        eastmost = maxNeg;
                        westmost = minPos;
                    }
                    lonRanges[0] = minNeg;
                    lonRanges[1] = maxNeg;
                    lonRanges[2] = minPos;
                    lonRanges[3] = maxPos;
                } /* nobs[ibin] != 0 */
            }     /* ibin loop */

            /* if no good obs left than bail */
            /* ----------------------------- */
            if (n_filled_bins == 0)
                goto freemem;

            /* Print info on filled row group */
            /* ------------------------------ */
            printf("%-20s:%8d\n", "# bins in row group", n_bins_in_group);
            printf("%-20s:%8d\n", "# filled bins", n_filled_bins);

            /* Write "Bin List" vdata */
            /* ---------------------- */
            if (is64bit)
            {
                writeBinList_nc(out_grp, n_filled_bins, (VOIDP)binList64nc);
                free(binList64nc);
            }
            else
            {
                writeBinList_nc(out_grp, n_filled_bins, (VOIDP)binList32nc);
                free(binList32nc);
            }

            /* Allocate sum & sum-squared arrays */
            /* --------------------------------- */
            sum_bin = (double *)calloc(n_bins_in_group, sizeof(double));
            sum2_bin = (double *)calloc(n_bins_in_group, sizeof(double));

            /* Loop over all L3 products to fill sum arrays */
            /* -------------------------------------------- */
            for (iprod = 0; iprod < l3b_nprod; iprod++)
            {
                memset(sum_bin, 0, n_bins_in_group * sizeof(double));
                memset(sum2_bin, 0, n_bins_in_group * sizeof(double));

                /* Process bins */
                /* ------------ */
                for (ibin = 0; ibin < n_bins_in_group; ibin++)
                {
                    if (nobs[ibin] == 0)
                        continue;

                    /* Process data file by file */
                    /* ------------------------- */
                    int32_t npix_file; // for weighting spatial average
                    double sum_file;   // accumulators for each file
                    double sum2_file;
                    double area_file;
                    int16_t prev_file; // index of previous file considered

                    // initialize sum for this bin with first file
                    j = 0;
                    float32 pixval = data_values[ibin][j * l3b_nprod + iprod];
                    double pixarea = data_areas[ibin][j];
                    npix_file = 1;
                    area_file = pixarea;
                    sum_file = pixval * pixarea;
                    sum2_file = pixval * pixarea * pixval * pixarea;
                    prev_file = file_index[ibin][j];

                    // add weighted data for rest of observations (files)
                    for (j = 1; j < nobs[ibin]; j++)
                    {
                        pixval = data_values[ibin][j * l3b_nprod + iprod];
                        pixarea = data_areas[ibin][j];

                        if (file_index[ibin][j] == prev_file)
                        { // same file
                            npix_file += 1;
                            area_file += pixarea;
                            sum_file += pixval * pixarea;
                            sum2_file += pixval * pixarea * pixval * pixarea;
                        }
                        else
                        { // new file

                            // finalize weighted sum for previous file
                            sum_bin[ibin] += (sum_file / sqrt(area_file));
                            sum2_bin[ibin] += (sum2_file / sqrt(area_file));

                            // initialize sum for current file
                            npix_file = 1;
                            area_file = pixarea;
                            sum_file = pixval * pixarea;
                            sum2_file = pixval * pixarea * pixval * pixarea;
                            prev_file = file_index[ibin][j];
                        }
                    } /* observation loop */

                    // add data from final file
                    sum_bin[ibin] += (sum_file / sqrt(area_file));
                    sum2_bin[ibin] += (sum2_file / sqrt(area_file));

                } /* ibin loop */

                /* Write Product Vdatas */
                /* -------------------- */
                float *productData = (float *)calloc(2 * n_filled_bins, sizeof(float));

                /* Fill bin data array */
                /* ------------------- */
                i = 0;
                for (ibin = 0; ibin < n_bins_in_group; ibin++)
                {
                    if (nobs[ibin] != 0)
                    {
                        productData[i * 2] = sum_bin[ibin];
                        productData[i * 2 + 1] = sum2_bin[ibin];
                        // memcpy(&productData[i * 8], &sum_bin[ibin], 4);
                        // memcpy(&productData[i * 8 + 4], &sum2_bin[ibin], 4);
                        i++;
                    }
                }

                // char_ptr1 = strchr(prodname[iprod], '/');
                // if (char_ptr1 != NULL) *char_ptr1 = '_';
                writeBinData_nc(out_grp, i, iprod, productData);
                // if (char_ptr1 != NULL) *char_ptr1 = '/';

                free(productData);

            } /* iprod loop */

            /* Write Quality vdata */
            /* ------------------- */
            if (input.qual_prod[0] != 0)
            {
                uint8_t *qualityData = (uint8_t *)calloc(n_filled_bins, 1);

                i = 0;
                for (ibin = 0; ibin < n_bins_in_group; ibin++)
                {
                    if (nobs[ibin] != 0)
                    {
                        memcpy(&qualityData[i], &best_qual[ibin], 1);
                        i++;
                    }
                }

                writeQuality_nc(out_grp, n_filled_bins, (VOIDP)qualityData);
                free(qualityData);
            }

            /* Free dynamic memory */
            /* ------------------- */
            if (sum_bin != NULL)
                free(sum_bin);
            if (sum2_bin != NULL)
                free(sum2_bin);

            /* Compute "begin" & "extent" vdata entries */
            /* ---------------------------------------- */
            binnum_data = (int64_t *)calloc(n_filled_bins, sizeof(int64_t));

            i = 0;
            for (ibin = 0; ibin < n_bins_in_group; ibin++)
            {
                if (nobs[ibin] != 0)
                {
                    binnum_data[i] = (int64_t)ibin + basebin[krow];

                    if (i < 0 || i >= n_filled_bins)
                    {
                        printf("Error: %d %d %d %d\n",
                               i, ibin, n_filled_bins, n_bins_in_group);
                    }
                    i++;
                }
            }

            get_beg_ext(n_filled_bins, binnum_data, basebin, nrows, beg, ext);

            free(binnum_data);

            total_filled_bins += n_filled_bins;

        } /* ********** n_filled_bin > 0 ********** */

        // free(best_qual);

        if (input.verbose)
        {
            time(&tnow);
            tmnow = localtime(&tnow);
            printf("krow:%5d After bin processing:  %s", krow, asctime(tmnow));
        }

        /* Fill BinIndex Vdata */
        /* ------------------- */
        for (i = 0; i < n_rows_in_group; i++)
        {

            i32 = i + krow;
            if (i32 < 0 || i32 >= nrows)
            {
                printf("Error: %d %d\n", i, krow);
                exit(-1);
            }

            if (is64bit)
            {
                binIndex64nc[i + krow].start_num = basebin[i32];
                binIndex64nc[i + krow].begin = beg[i32];
                binIndex64nc[i + krow].extent = ext[i32];
                binIndex64nc[i + krow].max = shape->getNumCols(i32);
            }
            else
            {
                binIndex32nc[i + krow].start_num = basebin[i32];
                binIndex32nc[i + krow].begin = beg[i32];
                binIndex32nc[i + krow].extent = ext[i32];
                binIndex32nc[i + krow].max = shape->getNumCols(i32);
            }

        } /* row_group loop */

        /* End Bin Index Fill */

    freemem:
        /* Free Dynamic Memory */
        /* ------------------- */
        if (bin_flag != NULL)
            free(bin_flag);
        if (tilt != NULL)
            free(tilt);
        if (qual != NULL)
            free(qual);
        if (nscenes != NULL)
            free(nscenes);
        if (lastfile != NULL)
            free(lastfile);

        for (ifile = 0; ifile < nfiles; ifile++)
        {
            if (scan_in_rowgroup[ifile] != NULL)
                free(scan_in_rowgroup[ifile]);
        }

        if (input.verbose)
        {
            time(&tnow);
            tmnow = localtime(&tnow);
            printf("krow:%5d Befre free dynic mem:  %s", krow, asctime(tmnow));
        }

        for (i = 0; i < n_bins_in_group; i++)
        {
            if (file_index[i] != NULL)
                free(file_index[i]);
            if (data_values[i] != NULL)
                free(data_values[i]);
            if (data_areas[i] != NULL)
                free(data_areas[i]);
            if (data_quality[i] != NULL)
                free(data_quality[i]);
            if (time_value[i] != NULL)
                free(time_value[i]);
        }

        if (input.verbose)
        {
            time(&tnow);
            tmnow = localtime(&tnow);
            printf("krow:%5d After free dynic mem:  %s", krow, asctime(tmnow));
        }

        free(data_values);
        free(data_areas);
        free(data_quality);
        free(time_value);
        free(nobs);
        free(allocated_space);
        free(file_index);

        krow += n_rows_in_group;

    } /* ========== End krow (Main) loop ========== */

    if (input.verbose)
    {
        time(&tnow);
        tmnow = localtime(&tnow);
        printf("krow:%5d %s", krow, asctime(tmnow));
    }
    printf("total_filled_bins: %ld\n", (long int)total_filled_bins);

    // {
    //     std::cout << "I/O (readL2) takes  " << total_count/1000.0 << " seconds" << std::endl;
    // }

    if (total_filled_bins == 0)
    {
        strcpy(buf, "rm -f ");
        strcat(buf, input.ofile);
        strcat(buf, "*");
        printf("%s\n", buf);
        system(buf);
        ret_status = 110;
        exit(ret_status);
    }

    /* Write and Close BinIndex Vdata */
    /* ------------------------------ */
    if (is64bit)
        writeBinIndex_nc(out_grp, nrows, binIndex64nc);
    else
        writeBinIndex_nc(out_grp, nrows, binIndex32nc);

    /*
     * determine list of files actually used in the bin output
     */
    if (input.fileuse[0] != 0)
    {
        fp2 = fopen(input.fileuse, "w");
        for (ifile = 0; ifile < nfiles; ifile++)
        {
            if (brk_scan[ifile] != -9999)
            {
                fileused[ifile] = 1;
                fprintf(fp2, "%s\n", l2_str[ifile].filename);
            }
        }
        fclose(fp2);
    }
    else
    {
        for (ifile = 0; ifile < nfiles; ifile++)
        {
            if (brk_scan[ifile] != -9999)
            {
                fileused[ifile] = 1;
            }
        }
    }

    /* Read and write global attributes */
    /* -------------------------------- */
    if (input.verbose)
        printf("Writing Global Attributes\n");

    status = reopenL2(0, &l2_str[0]);
    status = readL2meta(&meta_l2, 0);

    strcpy(meta_l3b.product_name, input.ofile);
    if (meta_l2.sensor_name != NULL)
        strcpy(meta_l3b.sensor_name, meta_l2.sensor_name);
    strcpy(meta_l3b.sensor, meta_l2.sensor);
    meta_l3b.sensorID = sensorID[0];
    strcpy(meta_l3b.title, meta_l3b.sensor_name);
    strcat(meta_l3b.title, " Level-3 Binned Data");

    if (input.suite) {
        strcat(meta_l3b.title, " ");
        strncat(meta_l3b.title, input.suite, sizeof(meta_l3b.title) - strlen(meta_l3b.title) - strlen(input.suite));
    }

    if (meta_l2.mission != NULL)
        strcpy(meta_l3b.mission, meta_l2.mission);
    strcpy(meta_l3b.prod_type, input.prodtype);

    strcat(meta_l3b.pversion, input.pversion);
    strcat(meta_l3b.soft_name, "l2bin");
    strcat(meta_l3b.soft_ver, VERSION);

    /*
     * loop through the input files
     *
     * set some ridiculous start and end times to get the ball rolling...
     * if this code is still in use in the 22nd century, damn we're good!
     * ...but this line will need tweaking.
     *
     */
    meta_l3b.end_orb = 0;
    meta_l3b.start_orb = 1000000;
    double startTime = yds2unix(2100, 1, 0);
    double endTime = yds2unix(1900, 1, 0);
    double tmpTime;
    strcpy(meta_l3b.infiles, basename(l2_str[0].filename));
    for (i = 0; i < nfiles; i++)
    {
        if (fileused[i] == 1)
        {
            // update orbit numbers
            if (l2_str[i].orbit > meta_l3b.end_orb)
                meta_l3b.end_orb = l2_str[i].orbit;
            if (l2_str[i].orbit < meta_l3b.start_orb)
                meta_l3b.start_orb = l2_str[i].orbit;

            // update start/end times
            tmpTime = yds2unix(l2_str[i].syear,
                               l2_str[i].sday,
                               (float64)(l2_str[i].smsec / 1000.0));
            if (tmpTime < startTime)
                startTime = tmpTime;
            tmpTime = yds2unix(l2_str[i].eyear,
                               l2_str[i].eday,
                               (float64)(l2_str[i].emsec / 1000.0));
            if (tmpTime > endTime)
                endTime = tmpTime;
        }
        // append to the file list - all files input, not just those used in the binning.
        if (i > 0)
        {
            strcat(meta_l3b.infiles, ",");
            strcat(meta_l3b.infiles, basename(l2_str[i].filename));
        }
    }

    meta_l3b.startTime = startTime;
    meta_l3b.endTime = endTime;

    strcpy(meta_l3b.binning_scheme, "Integerized Sinusoidal Grid");

    meta_l3b.north = northmost;
    meta_l3b.south = southmost;
    meta_l3b.east = eastmost;
    meta_l3b.west = westmost;

    strcpy(buf, argv[0]);
    for (i = 1; i < argc; i++)
    {
        strcat(buf, " ");
        strcat(buf, argv[i]);
    }
    strcpy(meta_l3b.proc_con, buf);
    strcpy(meta_l3b.input_parms, input.parms);
    strcpy(meta_l3b.flag_names, input.flaguse);

    buf[0] = 0;
    for (iprod = 0; iprod < l3b_nprod; iprod++)
    {
        char *prodName = strdup(l3_prodname[iprod].c_str());
        getL3units(&l2_str[0], 0, prodName, units);
        strcat(buf, prodName);
        free(prodName);
        strcat(buf, ":");
        strcat(buf, units);
        strcat(buf, ",");
        if (strlen(buf) >= MD_ATTRSZ)
        {
            printf("-E- units metadata is too long.  Remove some products from product list.\n");
            exit(EXIT_FAILURE);
        }
    }
    buf[strlen(buf) - 1] = 0;
    strcpy(meta_l3b.units, buf);

    meta_l3b.data_bins = total_filled_bins;
    int64_t totalNumBins;
    totalNumBins = basebin[nrows] - 1;
    meta_l3b.pct_databins = (float)(total_filled_bins) / totalNumBins * 100.0;

    strcpy(meta_l3b.doi, input.doi);

    const char *keywordStr = getGCMDKeywords(input.suite);
    if (keywordStr)
        strcpy(meta_l3b.keywords, keywordStr);

    time(&tnow);
    strcpy(meta_l3b.ptime, unix2isodate(tnow, 'G'));
    strcpy(meta_l3b.mission, sensorId2PlatformName(sensorID[0]));
    write_l3b_meta_netcdf4(ds_id, &meta_l3b, is64bit);

    freeL2meta(&meta_l2);

    /* Free Dynamic Memory */
    /* ------------------- */
    if (input.verbose)
        printf("Freeing Dynamic Memory\n");
    free(fileused);

    if (is64bit)
        free(binIndex64nc);
    else
        free(binIndex32nc);

    free(ext);
    free(beg);

    free(basebin);

    for (ifile = 0; ifile < nfiles; ifile++)
    {
        if (bscan_row[ifile] != NULL)
            free(bscan_row[ifile]);
        if (escan_row[ifile] != NULL)
            free(escan_row[ifile]);

        if (numer[ifile] != NULL)
            free(numer[ifile]);
        if (denom[ifile] != NULL)
            free(denom[ifile]);
    }

    /* Close L2 input files and free allocated L2 memory */
    /* ------------------------------------------------- */
    if (input.verbose)
        printf("Freeing L2 arrays\n");
    closeL2(&l2_str[0], 0);
    for (ifile = 0; ifile < nfiles; ifile++)
    {
        freeL2(&l2_str[ifile]);
    }
    freeL2(NULL);

    /* Close L3B output file */
    /* --------------------- */
    status = nc_close(fileid_w);

    return ret_status;
}

int64_t getbinnum(int32_t ifile, int32_t ipixl)
{

    return shape->latlon2bin(l2_str[ifile].latitude[ipixl], l2_str[ifile].longitude[ipixl]);
}
