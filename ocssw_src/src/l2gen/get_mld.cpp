#include <cstdio>
#include "l12_proto.h"
#include <sys/types.h>
#include <netcdf>
#include <vector>
#include <algorithm>

// module get_mld.c - retrieve mixed layer depth from ancillary data

static float mldbad = BAD_FLT;

#define MLDXANC 360
#define MLDYANC 180

char* mldclimatology_path() {
    char *ocdataroot;
    char *filepath = (char*) malloc(sizeof(char) * FILENAME_MAX);
    if ((ocdataroot = getenv("OCDATAROOT")) == NULL) {
        printf("-E- OCDATAROOT environment variable is not defined.\n");
        exit(EXIT_FAILURE);
    }
    sprintf(filepath, "%s/common/mld_climatology_woa1994.hdf", ocdataroot);
    return filepath;
}

/* This program opens the World Ocean Atlas Climatology hdf4 file
 and extracts mixed layer depth information. */

/*  :Title = "WOA Mixed Layer Depth Monthly Climatology"
    :Description = "World Ocean Atlas 1994, Mixed Layer Depth, Monthly"
    :Reference = "http://www.nodc.noaa.gov/OC5/WOA94/mix.html"
 */

float get_mld_hdf(float lon, float lat, int day) {
    static int firstCall = 1;
    static int mldx = MLDXANC;
    static int mldy = MLDYANC;
    static float dx = 360.0 / MLDXANC;
    static float dy = 180.0 / MLDYANC;

    typedef float ref_t[MLDXANC + 2];
    static ref_t* mldref;
    static ref_t* mld_sav;

    float mld = mldbad;
    int i, j, ii;
    float xx, yy;
    float t, u;

    if (firstCall) {

        float mldtmp[MLDYANC][MLDXANC];
        char name[H4_MAX_NC_NAME];
        char sdsname[H4_MAX_NC_NAME];
        int32 sd_id;
        int32 sds_id;
        int32 rank;
        int32 status;
        int32 sds_index;
        int32 nt;
        int32 dims[H4_MAX_VAR_DIMS];
        int32 nattrs;
        int32 start[2] = {0, 0};
        int32 ix, iy, ct;
        int32 lymin, lxmin, lymax, lxmax;
        float sum;

        char mldclimatology_file[FILENAME_MAX];
        char mldmonth[7];

        firstCall = 0;

        // allocate data
        mldref = (ref_t*) allocateMemory((MLDYANC + 2) * sizeof (ref_t), "mldref");
        mld_sav = (ref_t*) allocateMemory((MLDYANC + 2) * sizeof (ref_t), "mld_sav");

        strcpy(mldclimatology_file, mldclimatology_path());
        sd_id = SDstart(mldclimatology_file, DFACC_RDONLY);
        if (sd_id == -1) {
            printf("-E- %s:  Error opening MLD climatology file %s.\n",
                    __FILE__, mldclimatology_file);
            exit(EXIT_FAILURE);
        }
        printf("Opening MLD climatology file %s\n\n", mldclimatology_file);

        /* calculate the correct date and day of the month using yd2md utility. */
        int16 mon = (int) day / 31 + 1; // month of year
        // (no need for perfection..at least according to the sea surface salinity reference algorithm)

        sprintf(mldmonth, "mld%02i", mon);

        /* Get the SDS index */
        strcpy(sdsname, mldmonth);
        sds_index = SDnametoindex(sd_id, sdsname);
        if (sds_index == -1) {
            printf("-E- %s:  Error seeking %s SDS from %s.\n", __FILE__,
                    sdsname, mldclimatology_file);
            exit(EXIT_FAILURE);
        }
        sds_id = SDselect(sd_id, sds_index);

        status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
        if (status != 0) {
            printf("-E- %s:  Error reading SDS info for %s from %s.\n",
                    __FILE__, sdsname, mldclimatology_file);
            exit(EXIT_FAILURE);
        }
        status = SDreaddata(sds_id, start, NULL, dims, (VOIDP) & mldtmp[0][0]);
        if (status != 0) {
            printf("-E- %s:  Error reading SDS %s from %s.\n", __FILE__,
                    sdsname, mldclimatology_file);
            exit(EXIT_FAILURE);
        }

        status = SDendaccess(sds_id);
        status = SDend(sd_id);

        /* rotate 180-deg and add wrapping border to simplify interpolation */
        /* new grid is -180.5,180.5 in i=0,361 and -90.5,90.5 in j=0,181    */
        for (j = 0; j < mldy; j++) {
            for (i = 0; i < mldx; i++) {
                ii = (i < mldx / 2) ? i + mldx / 2 : i - mldx / 2;
                if (mldtmp[j][i] > 0)
                    mldref[j + 1][ii + 1] = mldtmp[j][i];
                else
                    mldref[j + 1][ii + 1] = mldbad;
            }
            mldref[j + 1][0] = mldref[j + 1][mldx];
            mldref[j + 1][mldx + 1] = mldref[j + 1][1];
        }
        for (i = 0; i < mldx + 2; i++) {
            mldref[0][i] = mldref[1][i];
            mldref[mldy + 1][i] = mldref[mldy][i];
        }

        /*
         *  Extend the grid by 1 point (use only full grid boxes later)
         */
        memcpy(mld_sav, mldref,
                (MLDYANC + 2) * (MLDXANC + 2) * sizeof (float));
        for (iy = 0; iy < mldy + 2; iy++) {
            lymin = (iy == 0) ? 0 : iy - 1;
            lymax = (iy == mldy + 1) ? mldy + 1 : iy + 1;

            for (ix = 0; ix < mldx + 2; ix++) {
                if (mldref[iy][ix] < mldbad + 1) {
                    lxmin = (ix == 0) ? 0 : ix - 1;
                    lxmax = (ix == mldx + 1) ? mldx + 1 : ix + 1;

                    sum = 0.;
                    ct = 0;
                    for (j = lymin; j < lymax + 1; j++)
                        for (i = lxmin; i < lxmax + 1; i++) {
                            if (mld_sav[j][i] > mldbad + 1) {
                                sum += mld_sav[j][i];
                                ct++;
                            }
                        }
                    if (ct > 0)
                        mldref[iy][ix] = sum / ct;
                }
            }
        }

    } /* end 1-time grid set up */

    /* locate LL position within reference grid */
    i = MAX(MIN((int) ((lon + 180.0 + dx / 2) / dx), MLDXANC + 1), 0);
    j = MAX(MIN((int) ((lat + 90.0 + dy / 2) / dy), MLDYANC + 1), 0);

    /* compute longitude and latitude of that grid element */
    xx = i * dx - 180.0 - dx / 2;
    yy = j * dy - 90.0 - dy / 2;

    /* Bilinearly interpolate, only using full grid boxes */
    t = (lon - xx) / dx;
    u = (lat - yy) / dy;

    if ((mldref[j][i] > mldbad + 1) && (mldref[j][i + 1] > mldbad + 1)
            && (mldref[j + 1][i] > mldbad + 1)
            && (mldref[j + 1][i + 1] > mldbad + 1)) {

        mld = (1 - t) * (1 - u) * mldref[j][i] + t * (1 - u) * mldref[j][i + 1]
                + t * u * mldref[j + 1][i + 1] + (1 - t) * u * mldref[j + 1][i];

    } else
        mld = mldbad;

    return (mld);
}

/*******************************************************************************/

using namespace std;
using namespace netCDF;
using namespace netCDF::exceptions;

int closest_index(vector<float> const& refvals, const float value) {
    // assumes refvals is sorted

    if (value < refvals.front() ||
        value > refvals.back())
        return -1;
    auto head = refvals.begin();
    auto tail = refvals.end();

    // iterator to first element with value not less than input (<=)
    auto next = lower_bound(head, tail, value);
    if (next == head) return next - head;

    // iterator to last element with value less than input (<)
    auto prev = next;
    --prev;
    if (prev == tail) return prev - head;

    // return index of element closer to value
    return (*next - value) < (value - *prev) ? next - head : prev - head;
}

float wrap(const float val, const float lo, const float hi) {
    const float width = hi - lo;
    const float offset = val - lo;
    return (offset - (floor(offset / width) * width)) + lo;
}

typedef struct mld_info_struct {
    NcVar mldvar;
    vector<float> lon;
    vector<float> lat;
} mld_info;

vector<float> get_coords(NcFile* ncid, vector<string> varnames) {
    NcDim Dim;
    NcVar Var;
    for (auto name = varnames.begin(); name != varnames.end(); ++name) {
        ncid->getCoordVar(*name, Dim, Var);
        if (!Dim.isNull()) {
            uint nvar = Dim.getSize();
            vector<float> var(nvar);
            Var.getVar(&var[0]);
            return var;
        }
    }
    // if we got this far, the coordinate variable wasn't found.
    return {};
}

int init_mld_nc(char* mldfile, mld_info* mldinfo) {
    try {
        NcFile* ncid = new NcFile(mldfile, NcFile::read);
        string varname;

        // get MLD info
        varname = "mixed_layer_thickness";
        mldinfo->mldvar = ncid->getVar(varname);

        // read longitude and latitude
        vector<string> lonnames { "lon", "Longitude", "longitude" };
        mldinfo->lon = get_coords(ncid, lonnames);
        if (mldinfo->lon.size() == 0) {
            cerr << "-X- " << __FILE__ " line " << __LINE__ << ": "
                 << "No Longitude coordinate variable in " <<  mldfile << endl;
            exit(EXIT_FAILURE);
        }

        vector<string> latnames { "lat", "Latitude", "latitude" };
        mldinfo->lat = get_coords(ncid, latnames);
        if (mldinfo->lat.size() == 0) {
            cerr << "-X- " << __FILE__ " line " << __LINE__ << ": "
                 << "No Latitude coordinate variable in " <<  mldfile << endl;
            exit(EXIT_FAILURE);
        }

        // check that lat, lon arrays are sorted
        vector<float> sorted;
        sorted = mldinfo->lon;
        sort(sorted.begin(), sorted.end());
        if (mldinfo->lon != sorted) {
            cerr << "-X- " << __FILE__ " line " << __LINE__ << ": "
                 << "Longitude array is not monotonically increasing\n";
            exit(EXIT_FAILURE);
        }
        sorted = mldinfo->lat;
        sort(sorted.begin(), sorted.end());
        if (mldinfo->lat != sorted) {
            cerr << "-X- " << __FILE__ " line " << __LINE__ << ": "
                 << "Latitude array is not monotonically increasing\n";
            exit(EXIT_FAILURE);
        }

    } catch (NcException const & e) {
        return 1;
    } catch (exception const & e) {
        e.what();
        cerr << "-X- " << __FILE__ " line " << __LINE__ << ": "
             << e.what() << endl;
        exit(EXIT_FAILURE);
    }

    return 0;
}

float get_mld_nc(char* mldfile, float lon, float lat) {
    static int firstCall = 1;
    static mld_info mldinfo;
    float mld;

    // initial file setups
    if (firstCall) {
        firstCall = 0;
        if (init_mld_nc(mldfile, &mldinfo)) {
            return (BAD_FLT);
        }
    }

    // normalize input coords and find closest index
    float keeplon = lon;
    lon = wrap(lon,  // put into same range as input lon
               mldinfo.lon.front(),
               mldinfo.lon.front() + 360);
    int ilon = closest_index(mldinfo.lon, lon);
    int ilat = closest_index(mldinfo.lat, lat);
    lon = keeplon;

    // look up mld value
    if ((ilon < 0) || (ilat < 0)) {
        mld = BAD_FLT;
    } else {
        try {
            vector<size_t> index = { 0, (size_t) ilat, (size_t) ilon };
            mldinfo.mldvar.getVar(index, &mld);

        } catch (exception const & e) {
            e.what();
            cerr << "-X- " << __FILE__ " line " << __LINE__ << ": "
                 << e.what()
                 << endl;
            exit(EXIT_FAILURE);
        }
    }

    return (mld);
}

/*******************************************************************************/

float get_mld(char* mldfile, float lon, float lat, int day) {
    static int netcdf_mld;
    static int firstCall = 1;
    int badval = 0;
    float mld;

    // initial file tests
    if (firstCall) {
        firstCall = 0;

        if (mldfile != NULL && mldfile[0] != 0) {
            // if the file is an HDF4 file, we will assume for now it is the climatology
            if (!Hishdf(mldfile)) {
                mld = get_mld_nc(mldfile, lon, lat);
                netcdf_mld = (mld != BAD_FLT); // is it valid format?
                if (netcdf_mld) {
                    printf("Opening mixed layer depth file: %s\n", mldfile);
                } else {
                    printf("Failed to read the specified mixed layer depth file: %s\n", mldfile);
                    exit(EXIT_FAILURE);
                }
            }
        } else {
            printf("Reading mixed layer depth info from default climatology\n");
            strcpy(mldfile, mldclimatology_path());
            netcdf_mld = 0;
        }
    } // firstCall done

    // try reading from netcdf file; check result
    if (netcdf_mld) {
        mld = get_mld_nc(mldfile, lon, lat);
        badval = (mld == BAD_FLT);
    }

    // read from climatology file as needed
    if (!netcdf_mld || badval) {
        mld = get_mld_hdf(lon, lat, day);
    }

    return (mld);
}
