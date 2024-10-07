#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>

#include "hdf5utils.h"
#include "passthebuck.h"

#define FAIL -1

#define DATACENTER "NASA/GSFC Aquarius Data Processing Center"
#define MISSIONCHARACTERISTICS "Nominal orbit: inclination=98.0 (Sun-synchronous); node=6PM (ascending); eccentricity=<0.002; altitude=657 km; ground speed=6.825 km/sec"

#define SENSORCHARACTERISTICS "Number of beams=3; channels per receiver=4; frequency 1.413 GHz; bits per sample=16; instatntaneous field of view=6.5 degrees; science data block period=1.44 sec."

// Millisecs since start of day

int32_t get_millisec(string *ydhmsf_str) {
    istringstream istr;
    int32_t itemp;

    istr.clear();
    istr.str(ydhmsf_str->substr(7, 2));
    istr >> itemp;
    int32_t millisec = itemp * 3600000;
    istr.clear();
    istr.str(ydhmsf_str->substr(9, 2));
    istr >> itemp;
    millisec += itemp * 60000;
    istr.clear();
    istr.str(ydhmsf_str->substr(11, 5));
    istr >> itemp;
    millisec += itemp;

    return millisec;
}

namespace Hdf {

int h5a_set(hid_t dataset, const char *nam, hid_t typ,
        hid_t cnt, void* data) {

    hid_t attr, atype, wr_typ, aid;
    hsize_t fdim[1];

    /* Save old error handler */
    H5E_auto_t old_func;
    void *old_client_data;
    H5Eget_auto(H5E_DEFAULT, &old_func, &old_client_data);

    /* Turn off error handling */
    H5Eset_auto(H5E_DEFAULT, NULL, NULL);

    /* Probe. Likely to fail, but that's okay */
    attr = H5Aopen_name(dataset, nam);

    /* Restore previous error handler */
    H5Eset_auto(H5E_DEFAULT, old_func, old_client_data);

    // Write attribute if it exists and return
    if (attr > 0) {
        if (H5Awrite(attr, typ, data)) {
            fprintf(stderr, "-E- %s line %d: ", __FILE__, __LINE__);
            fprintf(stderr, "H5Awrite(%ld,\"%s\",%ld,%ld,data) failed.  ",
                    (long) dataset, nam, (long) typ, (long) cnt);
            return (HDF_FUNCTION_ERROR);
        }
        H5Aclose(attr);
        return (LIFE_IS_GOOD);
    }


    if (typ == H5T_STRING) {
        /*
         * Create string attribute.
         */
        aid = H5Screate(H5S_SCALAR);

        atype = H5Tcopy(H5T_C_S1);
        H5Tset_size(atype, cnt);

        attr = H5Acreate1(dataset, nam, atype, aid, H5P_DEFAULT);
        wr_typ = atype;

    } else if (cnt == 1) {
        /*
         * Create scalar attribute.
         */
        aid = H5Screate(H5S_SCALAR);
        attr = H5Acreate1(dataset, nam, typ, aid, H5P_DEFAULT);
        wr_typ = typ;

    } else {
        aid = H5Screate(H5S_SIMPLE);
        fdim[0] = cnt;
        H5Sset_extent_simple(aid, 1, fdim, NULL);

        /*
         * Create array attribute.
         */
        attr = H5Acreate1(dataset, nam, typ, aid, H5P_DEFAULT);
        wr_typ = typ;
    }

    /*
     * Write attribute.
     */
    if (H5Awrite(attr, wr_typ, data)) {
        fprintf(stderr, "-E- %s line %d: ", __FILE__, __LINE__);
        fprintf(stderr, "H5Awrite(%ld,\"%s\",%ld,%ld,data) failed.  ",
                (long) dataset, nam, (long) typ, (long) cnt);
        return (HDF_FUNCTION_ERROR);
    }

    /*
     * Close attribute and file dataspaces.
     */
    H5Sclose(aid);

    /*
     * Close the attributes.
     */
    H5Aclose(attr);

    return (LIFE_IS_GOOD);
}

hid_t h5d_create(hid_t grp, const char *nam, hid_t typ,
        int rank,
        hsize_t d0, hsize_t d1, hsize_t d2,
        hsize_t d3, hsize_t d4, hsize_t d5,
        hid_t *dataset,
        hid_t *dataspace,
        hid_t plist) {

    hsize_t dimsizes[6];

    if (rank > 6) {
        fprintf(stderr, "-E- %s line %d: ", __FILE__, __LINE__);
        fprintf(stderr, "sd_create() expects to be passed a rank <= 6.  ");
        return (PROGRAMMER_BOOBOO);
    }
    dimsizes[0] = d0;
    dimsizes[1] = d1;
    dimsizes[2] = d2;
    dimsizes[3] = d3;
    dimsizes[4] = d4;
    dimsizes[5] = d5;

    *dataspace = H5Screate(H5S_SIMPLE);
    H5Sset_extent_simple(*dataspace, rank, dimsizes, NULL);

    if ((*dataset = H5Dcreate1(grp, nam, typ,
            *dataspace, plist)) == FAIL) {
        fprintf(stderr, "-E- %s line %d: ", __FILE__, __LINE__);
        fprintf(stderr, "H5Dcreate(%ld,\"%s\",%ld,%ld,[%ld,%ld,%ld]) failed.  ",
                (long) grp, nam, (long) typ, (long) rank,
                (long) d0, (long) d1, (long) d2);
        return (HDF_FUNCTION_ERROR);
    }

    return (LIFE_IS_GOOD);
}

int SetScalarH5A(hid_t id, const char *name, hid_t type, const void *value) {
    if (type == H5T_STRING)
        h5a_set(id, name, type, strlen((char *) value) + 1, (void*) value);
    else
        h5a_set(id, name, type, 1, (void*) value);

    return (LIFE_IS_GOOD);
}



/* ----------------------------------------------------- */
/* Create an H5D using the wrappers for the HDF routines */

/* ----------------------------------------------------- */
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
        int32_t d5, /* size of 6th dimension */
        const char *dn0, /* name of 1st dimension */
        const char *dn1, /* name of 2nd dimension */
        const char *dn2, /* name of 3rd dimension */
        const char *dn3, /* name of 4th dimension */
        const char *dn4, /* name of 5th dimension */
        const char *dn5, /* name of 6th dimension */
        hid_t plist /* Dataset property list */
        ) {

    hid_t dataset, dataspace;

    /* Create the SDS */
    PTB(h5d_create(grp, sname, nt, rank, d0, d1, d2, d3, d4, d5,
            &dataset, &dataspace, plist));

    /* Add a "long_name" attribute */
    PTB(SetScalarH5A(dataset, "long_name", (hid_t) H5T_STRING, lname));
    /* Add a "valid_range" attribute if one is specified */
    if (nt == H5T_NATIVE_UCHAR) {
        unsigned char vr[2];
        vr[0] = (unsigned char) low;
        vr[1] = (unsigned char) high;
        PTB(h5a_set(dataset, "valid_min", H5T_NATIVE_UCHAR, 1, (void*) & vr[0]));
        PTB(h5a_set(dataset, "valid_max", H5T_NATIVE_UCHAR, 1, (void*) & vr[1]));
    } else if (nt == H5T_NATIVE_USHORT) {
        short vr[2];
        vr[0] = (short) low;
        vr[1] = (short) high;
        PTB(h5a_set(dataset, "valid_min", H5T_NATIVE_USHORT, 1, (void*) & vr[0]));
        PTB(h5a_set(dataset, "valid_max", H5T_NATIVE_USHORT, 1, (void*) & vr[1]));
    } else if (nt == H5T_STD_I32LE) {
        int32_t vr[2];
        vr[0] = (int32_t) low;
        vr[1] = (int32_t) high;
        PTB(h5a_set(dataset, "valid_min", H5T_STD_I32LE, 1, (void*) & vr[0]));
        PTB(h5a_set(dataset, "valid_max", H5T_STD_I32LE, 1, (void*) & vr[1]));
    } else if (nt == H5T_STD_U32LE) {
        int32_t vr[2];
        vr[0] = (int32_t) low;
        vr[1] = (int32_t) high;
        PTB(h5a_set(dataset, "valid_min", H5T_STD_U32LE, 1, (void*) & vr[0]));
        PTB(h5a_set(dataset, "valid_max", H5T_STD_U32LE, 1, (void*) & vr[1]));
    } else if (nt == H5T_NATIVE_FLOAT) {
        float vr[2];
        vr[0] = (float) low;
        vr[1] = (float) high;
        PTB(h5a_set(dataset, "valid_min", H5T_NATIVE_FLOAT, 1, (void*) & vr[0]));
        PTB(h5a_set(dataset, "valid_max", H5T_NATIVE_FLOAT, 1, (void*) & vr[1]));
    } else if (nt == H5T_NATIVE_DOUBLE) {
        double vr[2];
        vr[0] = (float) low;
        vr[1] = (float) high;
        PTB(h5a_set(dataset, "valid_min", H5T_NATIVE_DOUBLE, 1, (void*) & vr[0]));
        PTB(h5a_set(dataset, "valid_max", H5T_NATIVE_DOUBLE, 1, (void*) & vr[1]));
    } else {
        fprintf(stderr, "-E- %s line %d: ", __FILE__, __LINE__);
        fprintf(stderr, "Got unsupported number type (%ld) ", (long) nt);
        fprintf(stderr, "while trying to create SDS, \"%s\", ", sname);
        exit(1);
    }


    /* Add a "slope" attribute if one is specified, and also
       an intercept attribute */
    if (slope != 0) {
        PTB(h5a_set(dataset, "slope", H5T_NATIVE_FLOAT, 1, (void*) & slope));
    }
    if (offset != 0) {
        PTB(h5a_set(dataset, "intercept", H5T_NATIVE_FLOAT, 1, (void*) & offset));
    }

    /* Add a "units" attribute if one is specified */
    if (units != NULL && *units != 0) {
        PTB(SetScalarH5A(dataset, "units", (hid_t) H5T_STRING, units));
    }


    /* Release this SDS */
    PTB(H5Sclose(dataspace));
    PTB(H5Dclose(dataset));

    return (LIFE_IS_GOOD);
}



/* ----------------------------------------------------- */
/* Create an H5D using the wrappers for the HDF routines */

/* ----------------------------------------------------- */
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
        ) {

    hid_t dataset;

    CreateH5D(grp, sname, lname, units, low, high, slope, offset,
            nt, rank, d0, d1, d2, d3, d4, d5,
            dn0, dn1, dn2, dn3, dn4, dn5, plist);

    dataset = H5Dopen1(grp, sname);

    if (nt == H5T_NATIVE_SHORT) {
        short _FillValue;
        _FillValue = (short) fillvalue;
        PTB(h5a_set(dataset, "_FillValue", H5T_NATIVE_SHORT, 1, (void*) & _FillValue));
    } else if (nt == H5T_STD_I32LE) {
        int32_t _FillValue;
        _FillValue = (int32_t) fillvalue;
        PTB(h5a_set(dataset, "_FillValue", H5T_STD_I32LE, 1, (void*) & _FillValue));
    } else if (nt == H5T_NATIVE_FLOAT) {
        float _FillValue;
        _FillValue = (float) fillvalue;
        PTB(h5a_set(dataset, "_FillValue", H5T_NATIVE_FLOAT, 1, (void*) & _FillValue));
    } else if (nt == H5T_NATIVE_DOUBLE) {
        double _FillValue;
        _FillValue = (float) fillvalue;
        PTB(h5a_set(dataset, "_FillValue", H5T_NATIVE_DOUBLE, 1, (void*) & _FillValue));
    } else {
        fprintf(stderr, "-E- %s line %d: ", __FILE__, __LINE__);
        fprintf(stderr, "Got unsupported number type (%ld) ", (long) nt);
        fprintf(stderr, "while trying to create SDS, \"%s\", ", sname);
        return (PROGRAMMER_BOOBOO);
    }


    // Standard Name
    if (strlen(stdname) != 0) {
        PTB(SetScalarH5A(dataset, "standard_name", (hid_t) H5T_STRING, stdname));
    }

    /* Release this SDS */
    PTB(H5Dclose(dataset));

    return (LIFE_IS_GOOD);
}



/*---------------------- */
/* Write to HDF5 dataset */

/* --------------------- */
herr_t h5d_write(hid_t id, const char *name, void* data, hsize_t rank,
        hsize_t s[6], hsize_t e[6]) {

    hid_t dataset, dataspace, filespace, datatype;
    hsize_t start[6], edge[6], dims[6];

    for (size_t i = 0; i < 6; i++) {
        start[i] = s[i];
        edge[i] = e[i];
    }


    dataset = H5Dopen1(id, name);
    if (dataset == -1) {
        printf("Dataset: %s not found (%d)\n", name, (int)id);
        exit(1);
    }

    /* Get datatype */
    /* ------------ */
    datatype = H5Dget_type(dataset);

    /*
     * Select a hyperslab.
     */
    filespace = H5Dget_space(dataset);


    // Check for out-of-bounds write
    H5Sget_simple_extent_dims(filespace, dims, NULL);
    for (size_t i = 0; i < rank; i++) {
        if (start[i] + edge[i] > dims[i]) {
            printf("Write to dim (0-based): %d of \"%s\" out of bounds (%d %d)\n",
                    (int) i, name, (int) (start[i] + edge[i]), (int) dims[i]);
            exit(110);
        }
    }

    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start, NULL,
            edge, NULL);

    dataspace = H5Screate_simple(rank, edge, NULL);

    /*
     * Write the data to the hyperslab.
     */
    H5Dwrite(dataset, datatype, dataspace, filespace,
            H5P_DEFAULT, data);

    H5Dclose(dataset);
    H5Sclose(dataspace);
    H5Sclose(filespace);

    return (LIFE_IS_GOOD);
}



/*----------------------- */
/* Read from HDF5 dataset */

/* ---------------------- */
herr_t h5d_read(hid_t id, const char *name, void* data, hsize_t rank,
        hsize_t s[6], hsize_t e[6]) {

    hid_t dataset, dataspace, filespace, datatype;
    hsize_t start[6], edge[6], dims[6];

    dataset = H5Dopen1(id, name);
    if (dataset == -1) {
        printf("Dataset: %s not found (%d)\n", name, (int)id);
        exit(1);
    }

    // Get datatype
    datatype = H5Dget_type(dataset);

    // Select hyperslab
    filespace = H5Dget_space(dataset);

    // Return dimensions if data = NULL
    if (data == NULL) {
        for (size_t i = 0; i < 6; i++) e[i] = 0;
        H5Sget_simple_extent_dims(filespace, e, NULL);
        H5Dclose(dataset);
        H5Sclose(filespace);
        return (LIFE_IS_GOOD);
    }

    // Copy start and edge
    for (size_t i = 0; i < 6; i++) {
        start[i] = s[i];
        edge[i] = e[i];
    }

    // Check for out-of-bounds read
    H5Sget_simple_extent_dims(filespace, dims, NULL);
    for (size_t i = 0; i < rank; i++) {
        if (start[i] + edge[i] > dims[i]) {
            printf("Read to dim (0-based): %d of \"%s\" out of bounds (%d %d)\n",
                    (int) i, name, (int) (start[i] + edge[i]), (int) dims[i]);
            exit(110);
        }
    }

    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start, NULL,
            edge, NULL);

    dataspace = H5Screate_simple(rank, edge, NULL);

    /*
     * Read the data from the hyperslab.
     */
    H5Dread(dataset, datatype, dataspace, filespace,
            H5P_DEFAULT, data);

    H5Dclose(dataset);
    H5Sclose(dataspace);
    H5Sclose(filespace);
    H5Tclose(datatype);

    return (LIFE_IS_GOOD);
}

} // namespace Hdf
