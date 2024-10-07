//-*- Mode: C++; -*-

#ifndef HDF5UTILS_H
#define HDF5UTILS_H

#include <sstream>

using namespace std;

#include "hdf5.h"
//#define VOIDP (void*)

#define DIFF1980_1970 315532800
#define DIFFJAN0680_1970 (315532800 + 432000)
#define DIFF2000_1980 631152000
#define DIFF2000_JAN0680 (631152000 - 432000)

int32_t get_millisec(string *ydhmsf_str);
double get_tai(int32_t year, int32_t doy, double millisec);
double get_tai(char *orbString);
double gpstai2utc2000(double gpstai);
double gpstai2unix(double gpstai);
double unix2gpstai(double unixtime);


namespace Hdf {

int CreateH5D(
        hid_t grp, /* group id */
        const char *sname, /* short name */
        const char *lname, /* long name */
        const char *units, /* units (not set if passed NULL or "") */
        double low, /* low end of valid range */
        double high, /* high end of range (no range set if low >= high) */
        float slope, /* scale factor (not set if 0)  */
        float offset, /* scaling offset (not set if 0)  */
        hid_t nt, /* HDF number type */
        int rank, /* number of dimensions (must be <= 3) */
        int32_t d0, /* size of 1st dimension */
        int32_t d1, /* size of 2nd dimension */
        int32_t d2, /* size of 3rd dimension */
        int32_t d3, /* size of 4th dimension */
        int32_t d4, /* size of 5th dimension */
        int32_t d5, /* size of 6t dimension  */
        const char *dn0, /* name of 1st dimension */
        const char *dn1, /* name of 2nd dimension */
        const char *dn2, /* name of 3rd dimension */
        const char *dn3, /* name of 4th dimension */
        const char *dn4, /* name of 5th dimension */
        const char *dn5, /* name of 6th dimension */
        hid_t plist /* Dataset property list */
        );

int CreateH5D(
        hid_t grp, /* group id */
        const char *sname, /* short name */
        const char *lname, /* long name */
        const char *stdname, /* standard name */
        const char *units, /* units (not set if passed NULL or "") */
        double low, /* low end of valid range */
        double high, /* high end of range (no range set if low >= high) */
        float slope, /* scale factor (not set if 0)  */
        float offset, /* scaling offset (not set if 0)  */
        float fillvalue, /* fill value */
        hid_t nt, /* HDF number type */
        int rank, /* number of dimensions (must be <= 3) */
        int32_t d0, /* size of 1st dimension */
        int32_t d1, /* size of 2nd dimension */
        int32_t d2, /* size of 3rd dimension */
        int32_t d3, /* size of 4th dimension */
        int32_t d4, /* size of 5th dimension */
        int32_t d5, /* size of 6th dimension */
        const char *dn0, /* name of 1st dimension */
        const char *dn1, /* name of 2nd dimension */
        const char *dn2, /* name of 3rd dimension */
        const char *dn3, /* name of 4th dimension */
        const char *dn4, /* name of 5th dimension */
        const char *dn5, /* name of 6th dimension */
        hid_t plist /* Dataset property list */
        );


int h5a_set(hid_t dataset, const char *nam, hid_t typ,
        hid_t cnt, void* data);

hid_t h5d_create(hid_t grp, const char *nam, hid_t typ,
        int rank,
        hsize_t d0, hsize_t d1, hsize_t d2,
        hsize_t d3, hsize_t d4, hsize_t d5,
        hid_t *dataset,
        hid_t *dataspace,
        hid_t plist);

herr_t h5d_read(hid_t id, const char *name, void* data, hsize_t rank,
        hsize_t s[6], hsize_t e[6]);

herr_t h5d_write(hid_t id, const char *name, void* data, hsize_t rank,
        hsize_t s[6], hsize_t e[6]);

int SetScalarH5A(hid_t id, const char *name, hid_t type, const void *value);

}
#endif

