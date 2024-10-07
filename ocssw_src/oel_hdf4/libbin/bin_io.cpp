#include "netcdf.h"
#include "hdf_bin.h"
#include "hdf5utils.h"
#include "passthebuck.h"
#include <seaproto.h>
#include <sensorInfo.h>
#include <timeutils.h>
#include <genutils.h>
#include <nc4utils.h>
#include <ncdfbin_utils.h>
#include "L3ShapeIsine.h"

#include <hdf.h>
#include <mfhdf.h>

// max # rows where the ISIN grid bin num fits into a 32bit int
#define MAX_ROWS_32BIT  41068

using namespace std;

//bool hdf5Type;

namespace Hdf {

/**
 * create an new bin file object and open it.
 * Will return an hdf4, hdf5 or netCDF bin file object.
 *
 * @param binFileName name of the bin file to open
 * @return new bin object that needs to be deleted.  NULL on failure.
 */
hdf_bin* openBinObject(const char* binFileName) {
    hdf_bin* binFile = NULL;
    int ncid;
    int status;

    if (Hishdf(binFileName) == TRUE) {
        binFile = new hdf4_bin;
    } else if (H5Fis_hdf5(binFileName) == TRUE) {
        char nam_buf[80];

        /* Turn off HDF5 warning/error handling - it's ugly stuff anyway */
        H5Eset_auto(H5E_DEFAULT, NULL, NULL);

        if (nc_open(binFileName, NC_NOWRITE, &ncid) == NC_NOERR) {
            status = nc_get_att(ncid, NC_GLOBAL, "Mission", nam_buf);
            if (status == NC_NOERR) {
                nc_close(ncid);
                if (strcmp(nam_buf, "SAC-D Aquarius") == 0) {
                    binFile = new hdf5_bin;
                }
            } else {
                binFile = new cdf4_bin;
            }
        }
    }

    if (binFile) {
        // open the bin file
        status = binFile->open(binFileName);
        if (status) {
            delete binFile;
            return NULL;
        }

        int numProducts = binFile->nprod();
        if (want_verbose) {
            printf("numProducts = %d\n", numProducts);
            char prodName[1024];
            for (int i = 0; i < numProducts; i++) {
                binFile->get_prodname(i, prodName);
                printf("  %2d - %s\n", i, prodName);
            }
        }
    }

    return binFile;
}

hdf_bin::hdf_bin() {
    binShape = NULL;
    totbins = 0;
    nrows = 0;

    binListPtr = 0;
    lastBinListPtr = -1;
    lastNumBins = 0;

    product_array[0] = NULL;
    n_data_prod = 0;
    n_data_records = 0;
    n_active_prod = 0;
    isHDF5 = false;
    isCDF4 = false;
    hasQual = false;
    hasNoext = true;
    deflate = 0;
}

hdf_bin::~hdf_bin() {
    if(binShape)
        delete binShape;
}

hdf4_bin::hdf4_bin() {
    file_id = -1;
    bin_ptr = 0;
    noext = 0;
    totbins = 0;
    binqual_idx = -1;

    //    cout << "in hdf4_bin constructor" << endl;
    //    for (size_t i=0; i<MAXNPROD; i++) product_array[i] = NULL;
    for (size_t i = 0; i < MAXNPROD; i++) aid[i] = -1;
    for (size_t i = 0; i < MAXNVDATA; i++) vdata_id[i] = -1;
}

hdf4_bin::~hdf4_bin() {
}

hdf5_bin::hdf5_bin() {
    isHDF5 = true;
    h5fid = -1;
    n_datasets = 0;
    binDataPtr = 0;

    for (size_t i = 0; i < MAXNPROD; i++) h5table_id[0][i] = -1;
    for (size_t i = 0; i < MAXNPROD; i++) h5table_id[1][i] = -1;
    for (size_t i = 0; i < MAXNPROD; i++) h5table_id[2][i] = -1;
}

hdf5_bin::~hdf5_bin() {
}

cdf4_bin::cdf4_bin() {
    isCDF4 = true;
    is64bit = false;
    ncid = -1;
    grp0 = -1;
    grp1 = -1;
    access_mode = NC_NOWRITE;
    n_datasets = 0;

    bindata_idx = -1;
    binqual_idx = -1;
    binlist_idx = -1;
    binindex_idx = -1;

    binDataPtr = 0;
    binQualityPtr = 0;
    binList = NULL;
    binList64 = NULL;
}

cdf4_bin::~cdf4_bin() {
}

// Create Bin File

int hdf4_bin::create(const char* l3b_filename, int32_t nrows) {
    uint8_t *a;
    int32_t i32;
    int32_t type[16];

    const char* fldname1[] = {"registration", "straddle", "bins", "radius",
        "max_north", "max_south", "seam_lon"};

    file_id = Hopen(l3b_filename, DFACC_CREATE, 0);
    sd_id = SDstart(l3b_filename, DFACC_RDWR);
    access_mode = DFACC_CREATE;
    this->nrows = nrows;
    if (this->hasNoext == true)
        noext = 1;

    Vstart(file_id);

    vg_id = Vattach(file_id, -1, "w");

    Vsetname(vg_id, "Level-3 Binned Data");
    Vsetclass(vg_id, "PlanetaryGrid");

    // Write "SEAGrid"
    // ---------------
    a = (uint8_t *) malloc(44);

    int32_t zero = 0;
    int32_t five = 5;

    double radius = 6378.137000;
    double north = 90.0;
    double south = -90.0;
    double seam_lon = -180.0;

    memcpy(&a[0], &five, 4);
    memcpy(&a[4], &zero, 4);
    i32 = 2 * nrows;
    memcpy(&a[8], &i32, 4);
    memcpy(&a[12], &radius, 8);
    memcpy(&a[20], &north, 8);
    memcpy(&a[28], &south, 8);
    memcpy(&a[36], &seam_lon, 8);

    type[0] = DFNT_INT32;
    type[1] = DFNT_INT32;
    type[2] = DFNT_INT32;
    type[3] = DFNT_FLOAT64;
    type[4] = DFNT_FLOAT64;
    type[5] = DFNT_FLOAT64;
    type[6] = DFNT_FLOAT64;

    seagrid_idx = 0;
    create_vdata(file_id, vg_id, &vdata_id[seagrid_idx], "SEAGrid",
            "Geometry", 7, fldname1, type, noext, NULL);

    write_vdata(vdata_id[seagrid_idx], 1, a);

    free(a);

    n_data_records = 0;

    binShape = new l3::L3ShapeIsine(nrows);
    totbins = binShape->getNumBins();
    
    binList = (binListStruct *) calloc(2 * nrows, sizeof (binListStruct));

    binqual_idx = -1;

    return 0;
}

// Open Bin File

int hdf4_bin::open(const char* l3b_filename) {
    int32_t vg_ref;
    int32_t tag;
    int32_t ref;
    char nam_buf[80];
    char cls_buf[80];

    hasQual = this->has_qual();
    n_data_prod = this->nprod();
    for (int32_t i = 0; i < n_data_prod; i++)
        this->get_prodname(i, proddata_name[i]);

    binqual_idx = -1;

    if ((file_id = Hopen(l3b_filename, DFACC_RDONLY, 0)) < 0) {
        fprintf(stderr, "Error: Cannot open input HDF file on Hopen - %s\n",
                l3b_filename);
        return FAIL;
    }

    Vstart(file_id);

    if ((sd_id = SDstart(l3b_filename, DFACC_RDONLY)) < 0) {
        fprintf(stderr, "Error: Cannot open input HDF file on SDstart- %s\n",
                l3b_filename);
        return FAIL;
    }

    access_mode = DFACC_RDONLY;

    // Read metadata
    if (read_l3b_meta_hdf4(sd_id, &meta_l3b)) {
        fprintf(stderr, "Error: Cannot read metadata from file - %s\n",
                l3b_filename);
        return FAIL;
    }

    vg_id = -1;
    for (size_t i = 0; i < MAXNVDATA; i++) vdata_id[i] = -1;

    vg_ref = Vfind(file_id, "Level-3 Binned Data");
    if (vg_ref > 0) {
        vg_id = Vattach(file_id, vg_ref, "r");

        //      status = HXsetdir((const char*) dirname(l3b_filename));

        // Loop through Vdatas
        bool first = true;
        for (int32_t i = 0; i < Vntagrefs(vg_id); i++) {
            Vgettagref(vg_id, i, &tag, &ref);
            vdata_id[i] = VSattach(file_id, ref, "r");

            if (vdata_id[i] == -1) {
                printf("Error opening Vdata (reference # %d)\n", ref);
                exit(-1);
            }

            VSgetname(vdata_id[i], nam_buf);
            VSgetclass(vdata_id[i], cls_buf);

            if (strcmp(cls_buf, "Geometry") == 0) {
                seagrid_idx = i;
            }

            if (strcmp(cls_buf, "Index") == 0) {
                nrows = VSelts(vdata_id[i]);
                binindex_idx = i;
            }

            if (strcmp(cls_buf, "DataMain") == 0) {
                n_data_records = VSelts(vdata_id[i]);
                binlist_idx = i;
            }

            if (strcmp(cls_buf, "DataSubordinate") == 0) {
                strcpy(proddata_name[n_data_prod++], nam_buf);
                if (first) {
                    bindata_idx = i;
                    first = false;
                }
                // Don't include quality vdata as product vdata
                //      ???  if (strcmp("qual_l3", nam_buf) == 0) n_data_prod--;
            }

            if (strcmp(cls_buf, "DataQuality") == 0) {
                binqual_idx = i;
            }
        }
    }

    binShape = new l3::L3ShapeIsine(nrows);
    totbins = binShape->getNumBins();

    binList = (binListStruct *) calloc(2 * nrows, sizeof (binListStruct));

    return SUCCEED;
}

// Read BinIndex

int hdf4_bin::readBinIndex(int row_num_to_read) {
    const char *binIndexFields =
            "row_num,vsize,hsize,start_num,begin,extent,max";

    // row_num_to_read is 0-based
    VSseek(vdata_id[binindex_idx], row_num_to_read);

    VSsetfields(vdata_id[binindex_idx], binIndexFields);

    int records_read = VSread(vdata_id[binindex_idx], (uint8_t *) & binIndex,
            1, NO_INTERLACE);

    // Correction for 64 bit
    if (sizeof (binIndex) == 40) {
        memmove(((char *) &binIndex) + 8,
                ((char *) &binIndex) + 4,
                32);
        memset(((char *) &binIndex) + 4, 0, 4);
    }

    return records_read;
}

// Read Bin File Product Data Size

int hdf_bin::read(char* product_list) {
    int n = strlen(product_list);
    char *buf = (char *) malloc(n + 1);
    strcpy(buf, product_list);

    for (int32_t i = 0; i < n_data_prod; i++)
        active_data_prod[i] = false;

    char *pch = strtok(buf, ",");
    while (pch != NULL) {
        if (strcmp(pch, "qual_l3") == 0) {
            pch = strtok(NULL, ",");
            continue;
        }

        bool found = false;
        for (int32_t i = 0; i < n_data_prod; i++) {
            if (strcmp(pch, proddata_name[i]) == 0) {
                found = true;
                active_data_prod[i] = true;
                break;
            }
        }
        if (!found) {
            cout << "Product: " << pch << " not found.";
            free(buf);
            return FAIL;
        }
        pch = strtok(NULL, ",");
    }

    n_active_prod = 0;
    for (int32_t i = 0; i < n_data_prod; i++)
        if (active_data_prod[i])
            n_active_prod++;

    free(buf);
    return (n_data_records * n_active_prod);
}

// Read Bin File Product Data

int hdf4_bin::read(float* data, binListStruct* binList, int nbins_to_read) {
    char buf2[512];
    int32_t records_read;

    int32_t records_remaining = VSelts(vdata_id[binlist_idx]) - bin_ptr;
    if (records_remaining < nbins_to_read) {
        cout << nbins_to_read << " requested, ";
        cout << records_remaining << " to be read at bin_ptr: ";
        cout << bin_ptr << endl;
        nbins_to_read = records_remaining;
    }

    VSsetfields(vdata_id[binlist_idx],
            "bin_num,nobs,nscenes,time_rec,weights,sel_cat,flags_set");

    read_binList(nbins_to_read, vdata_id[binlist_idx], binList, binShape);

    binListStruct binListRec;
    int k = 0;
    bool first = true;
    int bin_ptr0 = bin_ptr;
    for (int32_t i = 0; i < n_data_prod; i++) {
        if (active_data_prod[i] == true) {
            strcpy(buf2, proddata_name[i]);
            strcat(buf2, "_sum");
            VSsetfields(vdata_id[bindata_idx + i], buf2);

            VSseek(vdata_id[bindata_idx + i], bin_ptr0);

            records_read = VSread(vdata_id[bindata_idx + i],
                    (uint8_t *) & data[k * nbins_to_read],
                    nbins_to_read, FULL_INTERLACE);

            if (records_read == -1) {
                cout << "Read failure at bin_ptr: " << bin_ptr << endl;
                exit(4);
            }

            if (first) {
                bin_ptr += records_read;
                first = false;
            }

            for (int32_t j = 0; j < nbins_to_read; j++) {
                memcpy((void *) &binListRec, &binList[j],
                        sizeof (binListStruct));
                //        cout << binListRec.weights << endl;
                data[k * nbins_to_read + j] /= binListRec.weights;
            }
            k++;
        } // if active vdata
    } // vdata loop

    return records_read;
}

// Read Bin File Product Data

int hdf4_bin::read(float* data, float* var,
        binListStruct* binList, int nbins_to_read) {
    char buf2[512 * 2];
    int32_t records_read;

    int32_t records_remaining = VSelts(vdata_id[binlist_idx]) - bin_ptr;
    if (records_remaining < nbins_to_read) {
        cout << nbins_to_read << " requested, ";
        cout << records_remaining << " to be read at bin_ptr: ";
        cout << bin_ptr << endl;
        nbins_to_read = records_remaining;
    }

    VSsetfields(vdata_id[binlist_idx],
            "bin_num,nobs,nscenes,time_rec,weights,sel_cat,flags_set");

    read_binList(nbins_to_read, vdata_id[binlist_idx], binList, binShape);

    binListStruct binListRec;
    int k = 0;
    bool first = true;
    int bin_ptr0 = bin_ptr;
    for (int32_t i = 0; i < n_data_prod; i++) {
        if (active_data_prod[i] == true) {

            // Read sum
            strcpy(buf2, proddata_name[i]);
            strcat(buf2, "_sum");
            VSsetfields(vdata_id[bindata_idx + i], buf2);

            VSseek(vdata_id[bindata_idx + i], bin_ptr0);
            records_read = VSread(vdata_id[bindata_idx + i],
                    (uint8_t *) & data[k * nbins_to_read],
                    nbins_to_read, FULL_INTERLACE);

            if (records_read == -1) {
                cout << "Read failure at bin_ptr: " << bin_ptr << endl;
                exit(4);
            }

            // Read sum_sq
            strcpy(buf2, proddata_name[i]);
            strcat(buf2, "_sum_sq");
            VSsetfields(vdata_id[bindata_idx + i], buf2);

            VSseek(vdata_id[bindata_idx + i], bin_ptr0);
            records_read = VSread(vdata_id[bindata_idx + i],
                    (uint8_t *) & var[k * nbins_to_read],
                    nbins_to_read, FULL_INTERLACE);

            if (records_read == -1) {
                cout << "Read failure at bin_ptr: " << bin_ptr << endl;
                exit(4);
            }

            if (first) {
                bin_ptr += records_read;
                first = false;
            }

            for (int32_t j = 0; j < nbins_to_read; j++) {
                memcpy((void *) &binListRec, &binList[j],
                        sizeof (binListStruct));
                float sum = data[k * nbins_to_read + j];
                float wgt = binListRec.weights;
                data[k * nbins_to_read + j] /= wgt;

                if (binListRec.nscenes > 1) {
                    float tmp = var[k * nbins_to_read + j] * wgt;
                    tmp -= sum * sum;
                    tmp /= (wgt * wgt - binListRec.nobs);
                    var[k * nbins_to_read + j] = tmp;
                } else
                    var[k * nbins_to_read + j] = 1.0;
            }
            k++;
        } // if active vdata
    } // vdata loop

    return records_read;
}

// Read BinList

int hdf4_bin::readBinList(int nbins_to_read) {
    int32_t records_read;

    // return if we are reading the exact same data we already have
    if (lastBinListPtr == binListPtr && lastNumBins == nbins_to_read)
        return nbins_to_read;
    lastBinListPtr = binListPtr;
    lastNumBins = nbins_to_read;

    if (nbins_to_read == 0)
        return 0;

    const char *binListFields =
            "bin_num,nobs,nscenes,time_rec,weights,sel_cat,flags_set";

    VSsetfields(vdata_id[binlist_idx], binListFields);

    records_read = read_binList(nbins_to_read, vdata_id[binlist_idx], binList, binShape);

    binListPtr += records_read;

    return records_read;
}

// Read BinList

int hdf4_bin::readBinList() {
    int32_t records_read;

    records_read = this->readBinList(n_data_records);

    return records_read;
}

int hdf4_bin::readBinList(int nbins_to_read, int32_t list_reset_ptr) {
    int32_t records_read;

    VSseek(vdata_id[binlist_idx], list_reset_ptr);

    records_read = this->readBinList(nbins_to_read);

    binListPtr = list_reset_ptr + records_read;

    return records_read;
}

// Read Bin File Product Data

int hdf4_bin::read(float* data, binListStruct* binList) {
    int32_t records_read;

    records_read = this->read(data, binList, n_data_records);

    return records_read;
}

// Read Bin File Product Data/Var

int hdf4_bin::read(float* data, float* var, binListStruct* binList) {
    int32_t records_read;

    records_read = this->read(data, var, binList, n_data_records);

    return records_read;
}

int hdf4_bin::write(char *product_list, int32_t nwrite, float *outData,
        binListStruct* binList) {
    // Write BinList Vdata
    const char* fldname2[] = {"bin_num", "nobs", "nscenes", "time_rec",
        "weights", "sel_cat", "flags_set"};

    int32_t type[16];

    binlist_idx = 1;
    if (vdata_id[binlist_idx] == -1) {
        type[0] = DFNT_INT32;
        type[1] = DFNT_INT16;
        type[2] = DFNT_INT16;
        type[3] = DFNT_INT16;
        type[4] = DFNT_FLOAT32;
        type[5] = DFNT_UCHAR8;
        type[6] = DFNT_INT32;
        create_vdata(file_id, vg_id, &vdata_id[binlist_idx], "BinList",
                "DataMain", 7, fldname2, type, noext, NULL);
    }
    write_binList(nwrite, vdata_id[binlist_idx], binList);

    // Write Product Data Vdatas
    bindata_idx = 2; // Previous vdatas: SEAGrid, BinList
    char *fldname3[2];
    int n = strlen(product_list);
    char *buf = (char *) malloc(n + 1);
    strcpy(buf, product_list);
    char *pch = strtok(buf, ",");
    type[0] = DFNT_FLOAT32;
    type[1] = DFNT_FLOAT32;

    int i = 0;

    while (pch != NULL) {
        //      printf ("%s\n",pch);

        if (vdata_id[bindata_idx + i] == -1) {
            strcpy(proddata_name[n_data_prod++], pch);

            fldname3[0] = (char *) calloc(strlen(pch) + 5, sizeof (char));
            fldname3[1] = (char *) calloc(strlen(pch) + 8, sizeof (char));
            strcpy(fldname3[0], pch);
            strcat(fldname3[0], "_sum");
            strcpy(fldname3[1], pch);
            strcat(fldname3[1], "_sum_sq");

            create_vdata(file_id, vg_id, &vdata_id[bindata_idx + i], pch,
                    "DataSubordinate", 2, fldname3, type, noext, &aid[i]);

            free(fldname3[0]);
            free(fldname3[1]);
        }
        write_prodData(nwrite, vdata_id[bindata_idx + i],
                &outData[(bindata_idx + i - 2) * nwrite], binList);
        i++;
        pch = strtok(NULL, ",");
    }
    free(buf);

    n_data_records += nwrite;

    return 0;
}

int hdf4_bin::copy(char *product_list, int32_t nwrite, int32_t *binsToCopy,
        Hdf::binListStruct *inBinList, Hdf::hdf4_bin *input_binfile) {
    // Copy BinList Vdata
    const char* fldname2[] = {"bin_num", "nobs", "nscenes", "time_rec",
        "weights", "sel_cat", "flags_set"};

    int32_t type[16];

    binlist_idx = 1;
    if (vdata_id[binlist_idx] == -1) {
        type[0] = DFNT_INT32;
        type[1] = DFNT_INT16;
        type[2] = DFNT_INT16;
        type[3] = DFNT_INT16;
        type[4] = DFNT_FLOAT32;
        type[5] = DFNT_UCHAR8;
        type[6] = DFNT_INT32;
        create_vdata(file_id, vg_id, &vdata_id[binlist_idx], "BinList",
                "DataMain", 7, fldname2, type, noext, NULL);
    }
    for (int32_t i = 0; i < nwrite; i++) {
        write_binList(1, vdata_id[binlist_idx], &inBinList[binsToCopy[i]]);
    }

    // Copy Product Data Vdatas
    char *fldname3[2];
    int n = strlen(product_list);
    char *buf = (char *) malloc(n + 1);
    strcpy(buf, product_list);
    char *pch = strtok(buf, ",");
    type[0] = DFNT_FLOAT32;
    type[1] = DFNT_FLOAT32;

    int i = 0;

    while (pch != NULL) {
        strcpy(proddata_name[n_data_prod++], pch);

        fldname3[0] = (char *) calloc(strlen(pch) + 5, sizeof (char));
        fldname3[1] = (char *) calloc(strlen(pch) + 8, sizeof (char));
        strcpy(fldname3[0], pch);
        strcat(fldname3[0], "_sum");
        strcpy(fldname3[1], pch);
        strcat(fldname3[1], "_sum_sq");

        create_vdata(file_id, vg_id, &vdata_id[bindata_idx + i], pch,
                "DataSubordinate", 2, fldname3, type, noext, &aid[i]);

        int32_t ref = VSfind(input_binfile->file_id, pch);
        int32_t in_vdata_id = VSattach(input_binfile->file_id, ref, "r");
        if (in_vdata_id == -1) {
            printf("Error opening Vdata (reference # %d)\n", ref);
            exit(-1);
        }
        copy_prodData(nwrite, binsToCopy, fldname3,
                in_vdata_id, vdata_id[bindata_idx + i]);
        VSdetach(in_vdata_id);

        free(fldname3[0]);
        free(fldname3[1]);

        i++;
        pch = strtok(NULL, ",");
    }
    free(buf);

    n_data_records += nwrite;

    return 0;
}

// Read Quality

int hdf4_bin::readQual(uint8_t* qual, int32_t nbins_to_read) {
    int records_read = -1;
    if (binqual_idx != -1) {
        if (nbins_to_read > 0) {
            records_read = VSread(vdata_id[binqual_idx], (uint8_t *) qual,
                    nbins_to_read, FULL_INTERLACE);
        } else {
            return 0;
        }
    }
    return records_read;
}

int hdf4_bin::readQual(uint8_t* qual, int32_t nbins_to_read,
        int32_t row_num_to_read) {
    int records_read = -1;
    if (binqual_idx != -1) {
        if (nbins_to_read > 0) {
            // row_num_to_read is 0-based
            VSseek(vdata_id[binqual_idx], row_num_to_read);

            records_read = VSread(vdata_id[binqual_idx], (uint8_t *) qual,
                    nbins_to_read, FULL_INTERLACE);
        } else {
            return 0;
        }
    }
    return records_read;
}

// Read Bin File Product Data

int hdf4_bin::readSums(float* sums, int32_t nbins_to_read, int iprod) {
    char buf2[512 * 2];
    int32_t records_read;

    int k = 0;
    for (int32_t i = 0; i < n_data_prod; i++) {
        if (active_data_prod[i] == true && (i == iprod || iprod == -1)) {

            // Read sums
            strcpy(buf2, proddata_name[i]);
            strcat(buf2, "_sum,");
            strcat(buf2, proddata_name[i]);
            strcat(buf2, "_sum_sq");
            VSsetfields(vdata_id[bindata_idx + i], buf2);

            records_read = VSread(vdata_id[bindata_idx + i],
                    (uint8_t *) & sums[2 * k * nbins_to_read],
                    nbins_to_read, FULL_INTERLACE);

            if (records_read == -1) {
                //cout << "Read failure at bin_ptr: " << bin_ptr << endl;
                exit(4);
            }

            k++;
        } // if active vdata
    } // vdata loop

    return records_read;
}

// Read Bin File Product Data for a list of bin index numbers

int hdf4_bin::readSums(float* sums, int32_t* listOfBins,
        int32_t nbins_to_read, int iprod) {
    char buf2[512 * 2];
    int32_t records_read;

    int k = 0;
    for (int32_t i = 0; i < n_data_prod; i++) {
        if (active_data_prod[i] == true && (i == iprod || iprod == -1)) {

            // Read sums
            strcpy(buf2, proddata_name[i]);
            strcat(buf2, "_sum,");
            strcat(buf2, proddata_name[i]);
            strcat(buf2, "_sum_sq");
            VSsetfields(vdata_id[bindata_idx + i], buf2);

            for (int32_t j = 0; j < nbins_to_read; j++) {

                VSseek(vdata_id[bindata_idx + i], listOfBins[j]);

                records_read = VSread(vdata_id[bindata_idx + i],
                        (uint8_t *) & sums[2 * k], 1, FULL_INTERLACE);

                if (records_read == -1) {
                    //cout << "Read failure at bin_ptr: " << bin_ptr << endl;
                    exit(4);
                }

                k++;
            } // j-loop
        } // if active vdata
    } // vdata loop

    return records_read;
}

int hdf4_bin::setDataPtrAbsolute(int32_t recordNum) {
    for (int32_t i = 0; i < n_data_prod; i++) {
        if (VSseek(vdata_id[bindata_idx + i], recordNum) < 0)
            return -1;
    }
    return 0;
}

int hdf4_bin::writeBinList(int32_t nbins_to_write) {
    // Write BinList Vdata
    const char* fldname2[] = {"bin_num", "nobs", "nscenes", "time_rec",
        "weights", "sel_cat", "flags_set"};

    int32_t type[16];

    binlist_idx = 1;
    if (vdata_id[binlist_idx] == -1) {
        type[0] = DFNT_INT32;
        type[1] = DFNT_INT16;
        type[2] = DFNT_INT16;
        type[3] = DFNT_INT16;
        type[4] = DFNT_FLOAT32;
        type[5] = DFNT_UCHAR8;
        type[6] = DFNT_INT32;
        create_vdata(file_id, vg_id, &vdata_id[binlist_idx], "BinList",
                "DataMain", 7, fldname2, type, noext, NULL);
    }
    write_binList(nbins_to_write, vdata_id[binlist_idx], &binList[0]);

    n_data_records += nbins_to_write;

    return 0;
}

// Write Quality

int hdf4_bin::writeQual(uint8_t* qual, int32_t nbins_to_write) {
    const char* fldname5[] = {"qual_l3"};
    int32_t type[1];
    type[0] = DFNT_UCHAR8;

    if (binqual_idx == -1) {
        binqual_idx = bindata_idx + n_data_prod + 1;
        create_vdata(file_id, vg_id, &vdata_id[binqual_idx], "qual_l3",
                "DataQuality", 1, fldname5, type, noext, NULL);
    }

    int records_written = -1;
    if (binqual_idx != -1) {

        if (nbins_to_write > 0) {
            records_written = VSwrite(vdata_id[binqual_idx], (uint8_t *) qual,
                    nbins_to_write, FULL_INTERLACE);
        }
    }
    return records_written;
}

int hdf4_bin::writeSums(float* sums, int32_t nbins_to_write,
        const char *prodname) {
    bindata_idx = 2; // Previous vdatas: SEAGrid, BinList

    int iprod;
    bool found = false;
    for (iprod = 0; iprod < n_data_prod; iprod++) {
        if (strcmp(proddata_name[iprod], prodname) == 0) {
            found = true;
            break;
        }
    }

    char *fldname3[2];
    int32_t type[16];
    type[0] = DFNT_FLOAT32;
    type[1] = DFNT_FLOAT32;
    if (found == false) {
        strcpy(proddata_name[n_data_prod], prodname);
        iprod = n_data_prod;
        n_data_prod++;
    }

    if (vdata_id[bindata_idx + iprod] == -1) {
        fldname3[0] = (char *) calloc(strlen(prodname) + 5, sizeof (char));
        fldname3[1] = (char *) calloc(strlen(prodname) + 8, sizeof (char));

        strcpy(fldname3[0], prodname);
        strcat(fldname3[0], "_sum");
        strcpy(fldname3[1], prodname);
        strcat(fldname3[1], "_sum_sq");

        create_vdata(file_id, vg_id, &vdata_id[bindata_idx + iprod],
                prodname, "DataSubordinate", 2, fldname3, type,
                noext, &aid[iprod]);
        free(fldname3[0]);
        free(fldname3[1]);

    }
    int records_written = VSwrite(vdata_id[bindata_idx + iprod],
            (const uint8_t*) sums, nbins_to_write,
            FULL_INTERLACE);

    return records_written;
}

int hdf4_bin::close() {
    intn status;
    int32_t n_total_bins = 0;
    float slat = +90, nlat = -90, elon = -180, wlon = +180;

    if (access_mode == DFACC_CREATE) {
        int32_t *binnum_data = (int32_t *) calloc(n_data_records, sizeof (int32_t));

        status = VSsetfields(vdata_id[binlist_idx], "bin_num");
        if (status == -1) {
            printf("Unable to set \"bin_num\" in VSsetfields for id: %d.\n",
                    vdata_id[binlist_idx]);
            exit(1);
        }
        VSseek(vdata_id[binlist_idx], 0);
        int32_t nread = VSread(vdata_id[binlist_idx], (uint8_t *) binnum_data,
                n_data_records, FULL_INTERLACE);
        if (nread == -1) {
            printf("Unable to read bin number input in hdf4_bin::close()\n");
            exit(1);
        }

        // for (int32_t i=0; i<n_data_records-1; i++) {
        //     if (binnum_data[i] >= binnum_data[i+1]) {
        //         printf("Improper bin number: %d at element: %d\n",
        //                binnum_data[i], (int)i);
        //     }
        // }

        int32_t *beg = (int32_t *) calloc(nrows, sizeof (int32_t));
        int32_t *ext = (int32_t *) calloc(nrows, sizeof (int32_t));

        int32_t *numbin_out = (int32_t *) calloc(nrows, sizeof (int32_t));
        int32_t *basebin_out = (int32_t *) calloc(nrows + 1, sizeof (int32_t));

        for (int32_t i = 0; i < nrows; i++) {
            float latbin = (i + 0.5) * (180.0 / (nrows)) - 90.0;
            numbin_out[i] = (int32_t) (cos(latbin * PI / 180.0) * (2.0 * nrows) + 0.5);
            n_total_bins += numbin_out[i];
        }

        basebin_out[0] = 1;
        for (int32_t i = 1; i <= nrows; i++) {
            basebin_out[i] = basebin_out[i - 1] + numbin_out[i - 1];
        }

        /* Determine beg and extent values for each row */
        /* -------------------------------------------- */
        get_beg_ext32(n_data_records, binnum_data, basebin_out, nrows, beg, ext);

        /* Write BinIndex Vdata */
        /* -------------------- */
        uint8_t *a = (uint8_t *) calloc(36 * nrows, 1);
        int32_t i32 = 0;
        for (int32_t i = 0; i < nrows; i++) {
            double vsize = 180.0 / nrows;
            double hsize = 360.0 / numbin_out[i];
            memcpy(&a[i * 36], &i, 4);
            memcpy(&a[i * 36 + 4], &vsize, 8);
            memcpy(&a[i * 36 + 12], &hsize, 8);
            memcpy(&a[i * 36 + 20], &basebin_out[i], 4);
            memcpy(&a[i * 36 + 24], &beg[i], 4);
            memcpy(&a[i * 36 + 28], &ext[i], 4);
            memcpy(&a[i * 36 + 32], &numbin_out[i], 4);

            /* Determine NSEW boundries */
            /* ------------------------ */
            float lon, lat;
            if (beg[i] > 0) {
                lat = ((i + 0.5) / nrows) * 180.0 - 90.0;
                if (lat > nlat)
                    nlat = lat;
                if (lat < slat)
                    slat = lat;

                lon = 360.0 * (beg[i] - basebin_out[i] + 0.5)
                        / numbin_out[i] - 180.0;
                if (lon < wlon)
                    wlon = lon;

                lon = 360.0 * (binnum_data[i32 + ext[i] - 1] - basebin_out[i] + 0.5)
                        / numbin_out[i] - 180.0;
                if (lon > elon)
                    elon = lon;
            }
            i32 += ext[i];
        }

        const char* fldname4[] = {"row_num", "vsize", "hsize", "start_num",
            "begin", "extent", "max"};

        int32_t type[16];
        type[0] = DFNT_INT32;
        type[1] = DFNT_FLOAT64;
        type[2] = DFNT_FLOAT64;
        type[3] = DFNT_INT32;
        type[4] = DFNT_INT32;
        type[5] = DFNT_INT32;
        type[6] = DFNT_INT32;

        binindex_idx = 2 + n_data_prod;
        create_vdata(file_id, vg_id, &vdata_id[binindex_idx], "BinIndex",
                "Index", 7, fldname4, type, noext, NULL);

        write_vdata(vdata_id[binindex_idx], nrows, a);

        free(a);

        free(binnum_data);
        free(beg);
        free(ext);
        free(numbin_out);
        free(basebin_out);

        /* Write Global Attributes */
        /* ----------------------- */
        meta_l3b.north = nlat;
        meta_l3b.south = slat;
        meta_l3b.east = elon;
        meta_l3b.west = wlon;
        meta_l3b.data_bins = n_data_records;
        meta_l3b.pct_databins = 100 * ((float) n_data_records) / n_total_bins;

        if (nrows == 180) {
            meta_l3b.resolution = str2resolution("1deg");
        } else if (nrows == 360) {
            meta_l3b.resolution = str2resolution("0.5deg");
        } else if (nrows == 540) {
            meta_l3b.resolution = str2resolution("36km");
        } else if (nrows == 720) {
            meta_l3b.resolution = str2resolution("0.25deg");
        } else if (nrows == 1080) {
            meta_l3b.resolution = str2resolution("18km");
        } else if (nrows == 2160) {
            meta_l3b.resolution = str2resolution("9km");
        } else if (nrows == 4320) {
            meta_l3b.resolution = str2resolution("4km");
        } else if (nrows == 8640) {
            meta_l3b.resolution = str2resolution("2km");
        } else if (nrows == 17280) {
            meta_l3b.resolution = str2resolution("1km");
        } else if (nrows == 34560) {
            meta_l3b.resolution = str2resolution("hkm");
        } else if (nrows == 69120) {
            meta_l3b.resolution = str2resolution("qkm");
        } else if (nrows == 172800) {
            meta_l3b.resolution = str2resolution("hqkm");
        } else if (nrows == 345600) {
            meta_l3b.resolution = str2resolution("hhkm");
        } else {
            meta_l3b.resolution = BAD_FLT;
        }

        write_l3b_meta_hdf4(sd_id, &meta_l3b);

    } // CREATE

    free(binList);

    // Detach Vdatas
    status = VSdetach(vdata_id[seagrid_idx]);
    status = VSdetach(vdata_id[binlist_idx]);
    status = VSdetach(vdata_id[binindex_idx]);

    if (binqual_idx != -1)
        VSdetach(vdata_id[binqual_idx]);

    for (int32_t i = 0; i < n_data_prod; i++) {
        status = VSdetach(vdata_id[bindata_idx + i]);
        if (aid[i] != -1)
            Hendaccess(aid[i]);
    }

    // Detach HDF file structures
    status = Vdetach(vg_id);
    status = Vend(file_id);
    status = SDend(sd_id);
    status = Hclose(file_id);

    return 0;
}

// Check whether bin file has qual field

bool hdf4_bin::has_qual() {

    if (binqual_idx != -1) {
        hasQual = true;
    } else {
        hasQual = false;
    }
    return hasQual;
}

// Create Bin File

int hdf5_bin::create(const char *l3b_filename, int32_t nrows) {
    h5fid = H5Fcreate(l3b_filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (h5fid < 0) {
        cout << "Cannot create: " << l3b_filename << endl;
        exit(1);
    }

    this->nrows = nrows;

    grp0 = H5Gopen1(h5fid, "/");
    grp1 = H5Gcreate1(h5fid, "Level-3 Binned Data", 0);

    n_data_records = 0;

    binShape = new l3::L3ShapeIsine(nrows);
    totbins = binShape->getNumBins();

    binList = (binListStruct_hdf5 *) calloc(2 * nrows,
            sizeof (binListStruct_hdf5));

    return 0;
}

// Open Bin File

int hdf5_bin::open(const char* l3b_filename) {
    access_mode = H5F_ACC_RDONLY;
    h5fid = H5Fopen(l3b_filename, access_mode, H5P_DEFAULT);
    if (h5fid < 0) {
        cout << l3b_filename << " cannot be opened." << endl;
        exit(1);
    }

    grp0 = H5Gopen1(h5fid, "/");
    grp1 = H5Gopen1(h5fid, "Level-3 Binned Data");
    if (grp1 < 0) {
        cout << "Group \"Level-3 Binned Data\" cannot be opened in \""
                << l3b_filename << "\"" << endl;
        exit(1);
    }

    // Read metadata
    read_l3b_meta_hdf5(grp0, &meta_l3b);

    // Product Type Kludge
    strcpy(meta_l3b.prod_type, "");
    char *fstunderscore = (char *) strchr(l3b_filename, '_');
    if (fstunderscore != NULL) {
        char *sndunderscore = strchr(fstunderscore + 1, '_');
        if (sndunderscore != NULL) {
            char prodtypebuf[8];
            memset(prodtypebuf, 0, 8);
            memcpy(prodtypebuf, fstunderscore + 1,
                    sndunderscore - fstunderscore - 1);
            strcpy(meta_l3b.prod_type, prodtypebuf);
        }
    }

    hsize_t numTables;
    H5Gget_num_objs(grp1, &numTables);
    char nam_buf[80];

    bindata_idx = -1;

    bool first = true;
    for (size_t i = 0; i < numTables; i++) {
        H5Gget_objname_by_idx(grp1, i, nam_buf, 80);
        H5G_obj_t objType = H5Gget_objtype_by_idx(grp1, i);
        if (objType == H5G_DATASET) {

            // BinIndex
            if (strcmp(nam_buf, "BinIndex") == 0) {
                h5table_id[0][i] = H5Dopen1(grp1, "BinIndex");
                if (h5table_id[0][i] < 0) {
                    cout << "BinIndex cannot be opened" << endl;
                    exit(1);
                }

                binindex_idx = i;
                h5table_id[1][i] = H5Dget_space(h5table_id[0][i]);
                if (h5table_id[1][i] < 0) {
                    cout << "Filespace error for BinIndex" << endl;
                    exit(1);
                }

                h5table_id[2][i] = H5Dget_type(h5table_id[0][i]);
                if (h5table_id[2][i] < 0) {
                    cout << "Datatype error for BinIndex" << endl;
                    exit(1);
                }

                H5Sget_simple_extent_dims(h5table_id[1][i],
                        (hsize_t*) & nrows, NULL);

                n_datasets++;
                continue;
            } // BinIndex

            // BinList
            if (strcmp(nam_buf, "BinList") == 0) {
                h5table_id[0][i] = H5Dopen1(grp1, "BinList");
                if (h5table_id[0][i] < 0) {
                    cout << "BinList cannot be opened" << endl;
                    exit(1);
                }

                binlist_idx = i;
                h5table_id[1][i] = H5Dget_space(h5table_id[0][i]);
                if (h5table_id[1][i] < 0) {
                    cout << "Filespace error for BinList" << endl;
                    exit(1);
                }

                h5table_id[2][i] = H5Dget_type(h5table_id[0][i]);
                if (h5table_id[2][i] < 0) {
                    cout << "Datatype error for BinList" << endl;
                    exit(1);
                }

                H5Sget_simple_extent_dims(h5table_id[1][i],
                        (hsize_t*) & n_data_records, NULL);
                n_datasets++;
                continue;
            } // BinList

            // BinData Fields
            h5table_id[0][i] = H5Dopen1(grp1, nam_buf);
            if (h5table_id[0][i] < 0) {
                cout << nam_buf << " cannot be opened" << endl;
                exit(1);
            }

            h5table_id[1][i] = H5Dget_space(h5table_id[0][i]);
            if (h5table_id[1][i] < 0) {
                cout << "Filespace error for " << nam_buf << endl;
                exit(1);
            }

            h5table_id[2][i] = H5Dget_type(h5table_id[0][i]);
            if (h5table_id[2][i] < 0) {
                cout << "Datatype error for " << nam_buf << endl;
                exit(1);
            }

            strcpy(proddata_name[n_data_prod++], nam_buf);
            n_datasets++;
            if (first) {
                bindata_idx = i;
                first = false;
            }

        } // objType == H5G_DATASET
    }

    binShape = new l3::L3ShapeIsine(nrows);
    totbins = binShape->getNumBins();

    binList = (binListStruct_hdf5 *) calloc(2 * nrows,
            sizeof (binListStruct_hdf5));

    return SUCCEED;
}

// Read BinIndex

int hdf5_bin::readBinIndex(int row_num_to_read) {
    herr_t status;
    hsize_t one = 1;

    // row_num_to_read is 0-based
    hid_t dataset = h5table_id[0][binindex_idx];
    if (dataset < 0) {
        cout << "Improper BinIndex id" << endl;
        exit(1);
    }

    hid_t filespace = h5table_id[1][binindex_idx];
    if (filespace < 0) {
        cout << "Improper Filespace id for BinIndex" << endl;
        exit(1);
    }

    hid_t datatype = h5table_id[2][binindex_idx];
    if (datatype < 0) {
        cout << "Improper Datatype id for BinIndex" << endl;
        exit(1);
    }

    hsize_t strt = row_num_to_read;
    status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET,
            &strt, NULL, &one, NULL);
    if (status < 0) {
        cout << "H5Sselect_hyperslab error (readBinIndex)" << endl;
        exit(1);
    }

    hid_t dataspace = H5Screate_simple(1, &one, NULL);
    if (dataspace < 0) {
        cout << "Dataspace cannot be created for BinIndex." << endl;
        exit(1);
    }

    status = H5Dread(dataset, datatype, dataspace, filespace,
            H5P_DEFAULT, &binIndex);
    if (status < 0) {
        cout << "H5Dread error (readBinIndex)" << endl;
        exit(1);
    }

    H5Sclose(dataspace);

    // Correction for 64 bit
    if (sizeof (binIndex) == 40) {
        memmove(((char *) &binIndex) + 8,
                ((char *) &binIndex) + 4,
                32);
        memset(((char *) &binIndex) + 4, 0, 4);
    }

    return status;
}

// Read BinIndex

int hdf5_bin::readBinList(int nbins_to_read) {
    herr_t status;
    hsize_t n2read = nbins_to_read;

    // return good status if already read the requested data
    if (lastBinListPtr == binListPtr && lastNumBins == nbins_to_read)
        return 0;
    lastBinListPtr = binListPtr;
    lastNumBins = nbins_to_read;

    hid_t dataset = h5table_id[0][binlist_idx];
    if (dataset < 0) {
        cout << "Improper BinList id" << endl;
        exit(1);
    }

    hid_t filespace = h5table_id[1][binlist_idx];
    if (filespace < 0) {
        cout << "Improper Filespace id for BinList" << endl;
        exit(1);
    }

    hid_t datatype = h5table_id[2][binlist_idx];
    if (datatype < 0) {
        cout << "Improper Datatype id for BinList" << endl;
        exit(1);
    }

    hsize_t binListPtr_hsize;
    binListPtr_hsize = binListPtr;
    status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET,
            &binListPtr_hsize, NULL,
            &n2read, NULL);
    if (status < 0) {
        cout << "H5Sselect_hyperslab error (readBinList)" << endl;
        exit(1);
    }

    hid_t dataspace = H5Screate_simple(1, &n2read, NULL);
    if (dataspace < 0) {
        cout << "Dataspace cannot be created for BinList." << endl;
        exit(1);
    }

    status = H5Dread(dataset, datatype, dataspace, filespace,
            H5P_DEFAULT, binList);
    if (status < 0) {
        cout << "H5Dread error (readBinList)" << endl;
        exit(1);
    }

    binListPtr += nbins_to_read;

    H5Sclose(dataspace);

    return status;
}

int hdf5_bin::readBinList(int nbins_to_read, int32_t list_reset_ptr) {
    int32_t records_read;

    binListPtr = list_reset_ptr;

    records_read = this->readBinList(nbins_to_read);

    binListPtr = list_reset_ptr + records_read;

    return records_read;
}

// Read BinList

int hdf5_bin::readBinList() {
    // stub
    return 0;
}

// Read Bin File Quality Data

int hdf5_bin::readQual(unsigned char* qual, int32_t nbins_to_read) {
    // Currently no HDF5 binfiles have quality
    return -1;
}

int hdf5_bin::readQual(unsigned char* qual, int32_t nbins_to_read,
        int32_t row_num_to_read) {
    // Currently no HDF5 binfiles have quality
    return -1;
}

// Read Bin File Product Data

int hdf5_bin::readSums(float* sums, int32_t nbins_to_read, int iprod) {

    herr_t status;
    hsize_t n2read = nbins_to_read;

    int k = 0;
    for (int32_t i = 0; i < n_data_prod; i++) {
        if (active_data_prod[i] == true && (i == iprod || iprod == -1)) {

            hid_t dataset = h5table_id[0][bindata_idx + i];
            if (dataset < 0) {
                cout << "Improper BinData id for product: "
                        << iprod << endl;
                exit(1);
            }

            hid_t filespace = h5table_id[1][bindata_idx + i];
            if (filespace < 0) {
                cout << "Improper Filespace id for product: "
                        << iprod << endl;
                exit(1);
            }

            hid_t datatype = h5table_id[2][bindata_idx + i];
            if (datatype < 0) {
                cout << "Improper Datatype id for product: "
                        << iprod << endl;
                exit(1);
            }

            status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET,
                    &binDataPtr, NULL, &n2read, NULL);
            if (status < 0) {
                cout << "H5Sselect_hyperslab error for product: "
                        << iprod << endl;
                exit(1);
            }

            hid_t dataspace = H5Screate_simple(1, &n2read, NULL);
            if (dataspace < 0) {
                cout << "Dataspace cannot be created for product: "
                        << iprod << endl;
                exit(1);
            }

            status = H5Dread(dataset, datatype, dataspace, filespace,
                    H5P_DEFAULT, &sums[2 * k * nbins_to_read]);
            if (status < 0) {
                cout << "H5Dread error for product: "
                        << iprod << endl;
                exit(1);
            }

            H5Sclose(dataspace);

            k++;
        } // if active vdata
    } // vdata loop

    return status;
}

// Read Bin File Product Data for list of bins

int hdf5_bin::readSums(float* sums, int32_t* listOfBins,
        int32_t nbins_to_read, int iprod) {

    herr_t status;
    hsize_t one = 1;
    hsize_t binPtr;

    int k = 0;
    for (int32_t i = 0; i < n_data_prod; i++) {
        if (active_data_prod[i] == true && (i == iprod || iprod == -1)) {

            hid_t dataset = h5table_id[0][bindata_idx + i];
            if (dataset < 0) {
                cout << "Improper BinData id for product: "
                        << iprod << endl;
                exit(1);
            }

            hid_t filespace = h5table_id[1][bindata_idx + i];
            if (filespace < 0) {
                cout << "Improper Filespace id for product: "
                        << iprod << endl;
                exit(1);
            }

            hid_t datatype = h5table_id[2][bindata_idx + i];
            if (datatype < 0) {
                cout << "Improper Datatype id for product: "
                        << iprod << endl;
                exit(1);
            }

            hid_t dataspace = H5Screate_simple(1, &one, NULL);
            if (dataspace < 0) {
                cout << "Dataspace cannot be created for product: "
                        << iprod << endl;
                exit(1);
            }

            for (int32_t j = 0; j < nbins_to_read; j++) {
                binPtr = listOfBins[j];
                status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET,
                        &binPtr, NULL, &one, NULL);
                if (status < 0) {
                    cout << "H5Sselect_hyperslab error for product: "
                            << iprod << endl;
                    exit(1);
                }

                status = H5Dread(dataset, datatype, dataspace, filespace,
                        H5P_DEFAULT, &sums[2 * k]);
                if (status < 0) {
                    cout << "H5Dread error for product: "
                            << iprod << endl;
                    exit(1);
                }

                k++;
            }

            H5Sclose(dataspace);

        } // if active vdata
    } // vdata loop

    return status;
}

int hdf5_bin::write(const char *product_list, hsize_t nwrite, float *outData,
        binListStruct_hdf5* binList) {
    // Write BinList DataSet
    herr_t status;

    const char* fldname2[] = {"bin_num", "nobs", "nscenes",
        "weights", "flags_set"};
    binlist_idx = 0;
    hsize_t dim = nwrite;
    hsize_t unlimited = H5S_UNLIMITED;

    hid_t dataspace = H5Screate_simple(1, &dim, &unlimited);
    if (dataspace < 0) {
        cout << "Dataspace cannot be created (write)." << endl;
        exit(1);
    }

    if (binList != NULL) {
        if (h5table_id[0][binlist_idx] == -1) {
            hid_t type[5] = {
                H5T_STD_U32LE,
                H5T_NATIVE_SHORT,
                H5T_NATIVE_SHORT,
                H5T_NATIVE_FLOAT,
                H5T_STD_U64LE
            };
            size_t offset[5] = {
                HOFFSET(binListStruct_hdf5, bin_num),
                HOFFSET(binListStruct_hdf5, nobs),
                HOFFSET(binListStruct_hdf5, nscenes),
                HOFFSET(binListStruct_hdf5, weights),
                HOFFSET(binListStruct_hdf5, flags_set)
            };

            create_compound(grp1, "BinList",
                    &h5table_id[0][binlist_idx],
                    &h5table_id[2][binlist_idx],
                    sizeof (binListStruct_hdf5), 5, fldname2, offset, type,
                    &h5table_id[1][binlist_idx], dataspace);
        } else {
            hsize_t newsize = n_data_records + nwrite;
            status = H5Dextend(h5table_id[0][binlist_idx], &newsize);
            H5Sset_extent_simple(h5table_id[1][binlist_idx], 1, &newsize, NULL);
        }
        hid_t dataset = h5table_id[0][binlist_idx];
        if (dataset < 0) {
            cout << "Improper Dataset id for BinList (write)" << endl;
            exit(1);
        }

        hid_t filespace = h5table_id[1][binlist_idx];
        if (filespace < 0) {
            cout << "Improper Filespace id for BinList (write)" << endl;
            exit(1);
        }

        hid_t binlist_tid = h5table_id[2][binlist_idx];
        if (binlist_tid < 0) {
            cout << "Improper Datatype id for BinList (write)" << endl;
            exit(1);
        }

        const hsize_t start = n_data_records;
        const hsize_t count = nwrite;
        status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET,
                &start, NULL, &count, NULL);
        if (status < 0) {
            cout << "H5Sselect_hyperslab error for BinList (write)" << endl;
            exit(1);
        }

        status = H5Dwrite(dataset, binlist_tid, dataspace,
                filespace, H5P_DEFAULT, binList);
        if (status < 0) {
            cout << "Error writing BinList (write)." << endl;
            cout << "Dataset id:   " << dataset << endl;
            cout << "Datatype id:  " << binlist_tid << endl;
            cout << "Filespace id: " << filespace << endl;
            cout << "Start: " << n_data_records << endl;
            cout << "Count: " << nwrite << endl;
            exit(1);
        }
    }

    // if no products to write then return
    if (product_list == NULL) {
        H5Sclose(dataspace);
        return 0;
    }

    // Write Product Dataset
    bindata_idx = 1; // Previous datasets: BinList
    char *fldname3[2];
    int n = strlen(product_list);
    char *buf = (char *) malloc(n + 1);
    strcpy(buf, product_list);
    char *pch = strtok(buf, ",");

    while (pch != NULL) {

        char ds_name[200];
        int i = 0;
        while (h5table_id[0][bindata_idx + i] != -1) {
            H5Iget_name(h5table_id[0][bindata_idx + i], ds_name, 200);
            if (strcmp(pch, &ds_name[21]) == 0)
                break;
            i++;
        }

        if (h5table_id[0][bindata_idx + i] == -1) {

            strcpy(proddata_name[n_datasets++], pch);

            fldname3[0] = (char *) calloc(strlen(pch) + 5, sizeof (char));
            fldname3[1] = (char *) calloc(strlen(pch) + 8, sizeof (char));
            strcpy(fldname3[0], pch);
            strcat(fldname3[0], "_sum");
            strcpy(fldname3[1], pch);
            strcat(fldname3[1], "_sum_sq");

            hid_t type[2] = {H5T_NATIVE_FLOAT, H5T_NATIVE_FLOAT};
            size_t offset[2] = {0, 4};
            create_compound(grp1, pch,
                    &h5table_id[0][bindata_idx + i],
                    &h5table_id[2][bindata_idx + i],
                    8, 2, fldname3, offset, type,
                    &h5table_id[1][bindata_idx + i], dataspace);

            free(fldname3[0]);
            free(fldname3[1]);
        } else {
            hsize_t newsize = n_data_records + nwrite;
            status = H5Dextend(h5table_id[0][bindata_idx + i], &newsize);
            H5Sset_extent_simple(h5table_id[1][bindata_idx + i], 1,
                    &newsize, NULL);
        }

        hid_t dataset = h5table_id[0][bindata_idx + i];
        if (dataset < 0) {
            cout << "Improper Dataset id for: " << pch << endl;
            exit(1);
        }

        hid_t filespace = h5table_id[1][bindata_idx + i];
        if (filespace < 0) {
            cout << "Improper Filespace id for: " << pch << endl;
            exit(1);
        }

        hid_t datatype = h5table_id[2][bindata_idx + i];
        if (datatype < 0) {
            cout << "Improper Datatype id for: " << pch << endl;
            exit(1);
        }

        const hsize_t start = n_data_records;
        const hsize_t count = nwrite;
        status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET,
                &start, NULL, &count, NULL);
        if (status < 0) {
            cout << "H5Sselect_hyperslab error for: " << pch << endl;
            exit(1);
        }

        status = H5Dwrite(dataset, datatype,
                dataspace, filespace, H5P_DEFAULT, outData);
        if (status < 0) {
            cout << "Error writing: " << pch << endl;
            cout << "Dataset id:   " << dataset << endl;
            cout << "Datatype id:  " << datatype << endl;
            cout << "Filespace id: " << filespace << endl;
            cout << "Start: " << n_data_records << endl;
            cout << "Count: " << nwrite << endl;
            exit(1);
        }

        pch = strtok(NULL, ",");
    }
    free(buf);

    H5Sclose(dataspace);

    return 0;
}

int hdf5_bin::writeBinList(int32_t nbins_to_write) {
    // Write BinList DataSet
    herr_t status;

    const char* fldname2[] = {"bin_num", "nobs", "nscenes",
        "weights", "flags_set"};
    binlist_idx = 0;
    hsize_t dim = nbins_to_write;

    hsize_t unlimited = H5S_UNLIMITED;

    hid_t dataspace = H5Screate_simple(1, &dim, &unlimited);

    if (h5table_id[0][binlist_idx] == -1) {
        hid_t type[5] = {
            H5T_STD_U32LE,
            H5T_NATIVE_SHORT,
            H5T_NATIVE_SHORT,
            H5T_NATIVE_FLOAT,
            H5T_STD_U64LE
        };
        size_t offset[5] = {
            HOFFSET(binListStruct_hdf5, bin_num),
            HOFFSET(binListStruct_hdf5, nobs),
            HOFFSET(binListStruct_hdf5, nscenes),
            HOFFSET(binListStruct_hdf5, weights),
            HOFFSET(binListStruct_hdf5, flags_set)
        };

        create_compound(grp1, "BinList",
                &h5table_id[0][binlist_idx],
                &h5table_id[2][binlist_idx],
                sizeof (binListStruct_hdf5), 5, fldname2, offset, type,
                &h5table_id[1][binlist_idx], dataspace);

        n_datasets++;
    } else {
        hsize_t newsize = n_data_records + nbins_to_write;
        status = H5Dextend(h5table_id[0][binlist_idx], &newsize);
        H5Sset_extent_simple(h5table_id[1][binlist_idx], 1, &newsize, NULL);
    }
    hid_t dataset = h5table_id[0][binlist_idx];
    if (dataset < 0) {
        cout << "Improper Dataset id for BinList" << endl;
        exit(1);
    }

    hid_t filespace = h5table_id[1][binlist_idx];
    if (filespace < 0) {
        cout << "Improper Filespace id for BinList (writeBinList)." << endl;
        exit(1);
    }

    //    cout << H5Sget_simple_extent_npoints(filespace) << endl;
    //hsize_t d;
    //H5Sget_simple_extent_dims(filespace, &d, NULL);
    //cout << "d: " << d << endl;

    hid_t binlist_tid = h5table_id[2][binlist_idx];
    if (binlist_tid < 0) {
        cout << "Improper Datatype id for BinList (writeBinList)." << endl;
        exit(1);
    }

    const hsize_t start = n_data_records;
    const hsize_t count = nbins_to_write;
    status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, &start,
            NULL, &count, NULL);

    status = H5Dwrite(dataset, binlist_tid, dataspace,
            filespace, H5P_DEFAULT, binList);
    if (status < 0) {
        cout << "Error writing BinList (writeBinList)." << endl;
        cout << "Dataset id:   " << dataset << endl;
        cout << "Datatype id:  " << binlist_tid << endl;
        cout << "Filespace id: " << filespace << endl;
        cout << "Start: " << n_data_records << endl;
        cout << "Count: " << nbins_to_write << endl;
        exit(1);
    }

    H5Sclose(dataspace);

    return 0;
}

int hdf5_bin::writeSums(float* sums, int32_t nbins_to_write,
        const char *prodname) {
    int iprod;
    bool found = false;
    for (iprod = 0; iprod < n_data_prod; iprod++) {
        if (strcmp(proddata_name[iprod], prodname) == 0) {
            found = true;
            break;
        }
    }
    bindata_idx = 1; // Previous Datasets: BinList
    hsize_t dim = nbins_to_write;

    hsize_t unlimited = H5S_UNLIMITED;
    hid_t dataspace = H5Screate_simple(1, &dim, &unlimited);

    char *fldname3[2];
    herr_t status;

    //    int i = 0;

    if (found == false) {
        strcpy(proddata_name[n_data_prod], prodname);

        fldname3[0] = (char *) calloc(strlen(prodname) + 5, sizeof (char));
        fldname3[1] = (char *) calloc(strlen(prodname) + 8, sizeof (char));

        strcpy(fldname3[0], prodname);
        strcat(fldname3[0], "_sum");
        strcpy(fldname3[1], prodname);
        strcat(fldname3[1], "_sum_sq");
        hid_t type[2] = {H5T_NATIVE_FLOAT, H5T_NATIVE_FLOAT};
        size_t offset[2] = {0, 4};
        create_compound(grp1, prodname,
                &h5table_id[0][bindata_idx + n_data_prod],
                &h5table_id[2][bindata_idx + n_data_prod],
                8, 2, fldname3, offset, type,
                &h5table_id[1][bindata_idx + n_data_prod], dataspace);

        free(fldname3[0]);
        free(fldname3[1]);

        n_data_prod++;
        n_datasets++;
    } else {
        hsize_t newsize = n_data_records + nbins_to_write;
        status = H5Dextend(h5table_id[0][bindata_idx + iprod], &newsize);
        H5Sset_extent_simple(h5table_id[1][bindata_idx + iprod], 1,
                &newsize, NULL);
    }
    hid_t dataset = h5table_id[0][bindata_idx + iprod];
    hid_t filespace = h5table_id[1][bindata_idx + iprod];
    hid_t bindata_tid = h5table_id[2][bindata_idx + iprod];
    //    hsize_t nn1=H5Sget_simple_extent_npoints(filespace);

    const hsize_t start = n_data_records;
    const hsize_t count = nbins_to_write;
    status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, &start,
            NULL, &count, NULL);

    status = H5Dwrite(h5table_id[0][bindata_idx + iprod], bindata_tid,
            dataspace, filespace, H5P_DEFAULT, sums);
    if (status < 0) {
        cout << "Error writing product: " << iprod << endl;
        cout << "Dataset id:   " << dataset << endl;
        cout << "Datatype id:  " << bindata_tid << endl;
        cout << "Filespace id: " << filespace << endl;
        cout << "Start: " << n_data_records << endl;
        cout << "Count: " << nbins_to_write << endl;
        exit(1);
    }

    H5Sclose(dataspace);

    return 0;
}

int hdf5_bin::close() {
    herr_t status;
    int32_t n_total_bins = 0;
    float slat = +90, nlat = -90, elon = -180, wlon = +180;

    // Save old error handler
    H5E_auto_t old_func;
    void *old_client_data;
    H5Eget_auto(H5E_DEFAULT, &old_func, &old_client_data);

    H5Eset_auto(H5E_DEFAULT, NULL, NULL);
    hid_t dataset = H5Dopen1(grp1, "BinIndex");

    // Restore previous error handler
    H5Eset_auto(H5E_DEFAULT, old_func, old_client_data);

    if (dataset == -1) {
        int32_t *binnum_data = (int32_t *) calloc(n_data_records, sizeof (int32_t));

        hid_t tid = H5Tcreate(H5T_COMPOUND, sizeof (int32_t));
        // Note: Name of field must match original name
        status = H5Tinsert(tid, "bin_num", 0, H5T_STD_U32LE);

        status = H5Dread(h5table_id[0][binlist_idx], tid,
                H5S_ALL, H5S_ALL, H5P_DEFAULT, binnum_data);

        H5Tclose(tid);
        //      goto skip;

        for (int32_t i = 0; i < n_data_records - 1; i++) {
            if (binnum_data[i] >= binnum_data[i + 1]) {
                printf("Improper bin number: %ld at element: %d\n",
                        (long int) binnum_data[i], (int) i);
            }
        }

        int32_t *beg = (int32_t *) calloc(nrows, sizeof (int32_t));
        int32_t *ext = (int32_t *) calloc(nrows, sizeof (int32_t));

        int32_t *numbin_out = (int32_t *) calloc(nrows, sizeof (int32_t));
        int32_t *basebin_out = (int32_t *) calloc(nrows + 1, sizeof (int32_t));

        for (int32_t i = 0; i < nrows; i++) {
            float latbin = (i + 0.5) * (180.0 / (nrows)) - 90.0;
            numbin_out[i] = (int32_t) (cos(latbin * PI / 180.0) * (2.0 * nrows) + 0.5);
            n_total_bins += numbin_out[i];
        }

        basebin_out[0] = 1;
        for (int32_t i = 1; i <= nrows; i++) {
            basebin_out[i] = basebin_out[i - 1] + numbin_out[i - 1];
        }

        /* Determine beg and extent values for each row */
        /* -------------------------------------------- */
        get_beg_ext32(n_data_records, binnum_data, basebin_out, nrows, beg, ext);

        /* Write BinIndex Vdata */
        /* -------------------- */
        uint8_t *a = (uint8_t *) calloc(36 * nrows, 1);
        int32_t i32 = 0;
        for (int32_t i = 0; i < nrows; i++) {
            double vsize = 180.0 / nrows;
            double hsize = 360.0 / numbin_out[i];
            memcpy(&a[i * 36], &i, 4);
            memcpy(&a[i * 36 + 4], &vsize, 8);
            memcpy(&a[i * 36 + 12], &hsize, 8);
            memcpy(&a[i * 36 + 20], &basebin_out[i], 4);
            memcpy(&a[i * 36 + 24], &beg[i], 4);
            memcpy(&a[i * 36 + 28], &ext[i], 4);
            memcpy(&a[i * 36 + 32], &numbin_out[i], 4);

            /* Determine NSEW boundries */
            /* ------------------------ */
            float lon, lat;
            if (beg[i] > 0) {
                lat = ((i + 0.5) / nrows) * 180.0 - 90.0;
                if (lat > nlat)
                    nlat = lat;
                if (lat < slat)
                    slat = lat;

                lon = 360.0 * (beg[i] - basebin_out[i] + 0.5)
                        / numbin_out[i] - 180.0;
                if (lon < wlon)
                    wlon = lon;

                lon = 360.0 * (binnum_data[i32 + ext[i] - 1] - basebin_out[i] + 0.5)
                        / numbin_out[i] - 180.0;
                if (lon > elon)
                    elon = lon;
            }
            i32 += ext[i];
        }

        const char* fldname4[] = {"row_num", "vsize", "hsize", "start_num",
            "begin", "extent", "max"};

        hid_t type[7];
        type[0] = H5T_STD_I32LE;
        type[1] = H5T_NATIVE_DOUBLE;
        type[2] = H5T_NATIVE_DOUBLE;
        type[3] = H5T_STD_I32LE;
        type[4] = H5T_STD_I32LE;
        type[5] = H5T_STD_I32LE;
        type[6] = H5T_STD_I32LE;

        binindex_idx = n_datasets;
        size_t offset[7] = {0, 4, 12, 20, 24, 28, 32};

        hsize_t start[1] = {0};
        hsize_t dim = nrows;
        hid_t dataspace = H5Screate_simple(1, &dim, NULL);
        create_compound(grp1, "BinIndex",
                &h5table_id[0][binindex_idx],
                &h5table_id[2][binindex_idx],
                36, 7, fldname4, offset, type,
                &h5table_id[1][binindex_idx], dataspace);

        status = H5Sselect_hyperslab(h5table_id[1][binindex_idx],
                H5S_SELECT_SET, start, NULL, &dim, NULL);

        status = H5Dwrite(h5table_id[0][binindex_idx],
                h5table_id[2][binindex_idx], dataspace,
                h5table_id[1][binindex_idx], H5P_DEFAULT, a);
        if (status < 0) {
            cout << "Error writing BinIndex." << endl;
            exit(1);
        }

        H5Sclose(dataspace);

        n_datasets++;

        free(a);

        free(binnum_data);
        free(beg);
        free(ext);
        free(numbin_out);
        free(basebin_out);

        // Write Global Attributes
        // -----------------------
        char buf[50000];

        // Check if metadata already written (l2bin_aquarius)
        strcpy(buf, this->meta_l3b.product_name);
        if (strcmp(buf, "___") == 0)
            goto skip;

        PTB(SetScalarH5A(grp0, "Product Name", H5T_STRING, buf));

        strcpy(buf, this->meta_l3b.sensor_name);
        strcat(buf, " Level-3 Binned Data");
        PTB(SetScalarH5A(grp0, "Title", H5T_STRING, buf));

        strcpy(buf, this->meta_l3b.sensor);
        PTB(SetScalarH5A(grp0, "Sensor", H5T_STRING, buf));

        strcpy(buf, this->meta_l3b.sensor_char);
        PTB(SetScalarH5A(grp0, "Sensor Characteristics", H5T_STRING, buf));

        strcpy(buf, this->meta_l3b.mission);
        PTB(SetScalarH5A(grp0, "Mission", H5T_STRING, buf));

        strcpy(buf, this->meta_l3b.mission_char);
        PTB(SetScalarH5A(grp0, "Mission Characteristics", H5T_STRING, buf));

        strcpy(buf, this->meta_l3b.soft_name);
        PTB(SetScalarH5A(grp0, "Software Name", H5T_STRING, buf));

        strcpy(buf, this->meta_l3b.soft_ver);
        PTB(SetScalarH5A(grp0, "Software ID", H5T_STRING, buf));

        strcpy(buf, this->meta_l3b.ptime);
        PTB(SetScalarH5A(grp0, "Processing Time", H5T_STRING, buf));

        strcpy(buf, this->meta_l3b.proc_con);
        PTB(SetScalarH5A(grp0, "Processing Control", H5T_STRING, buf));

        strcpy(buf, this->meta_l3b.flag_names);
        PTB(SetScalarH5A(grp0, "L2 Flag Names", H5T_STRING, buf));

        strcpy(buf, this->meta_l3b.infiles);
        PTB(SetScalarH5A(grp0, "Input Files", H5T_STRING, buf));

        PTB(SetScalarH5A(grp0, "Start Orbit", H5T_STD_U32LE,
                (VOIDP) & this->meta_l3b.start_orb));
        PTB(SetScalarH5A(grp0, "End Orbit", H5T_STD_U32LE,
                (VOIDP) & this->meta_l3b.end_orb));

        strcpy(buf, ydhmsf(this->meta_l3b.startTime, 'G'));
        PTB(SetScalarH5A(grp0, "Start Time", H5T_STRING, buf));

        int16_t syear, sday, eyear, eday;
        double msec;
        int32_t smsec, emsec;
        unix2yds(meta_l3b.startTime, &syear, &sday, &msec);
        smsec = (int32_t) msec;
        unix2yds(meta_l3b.endTime, &eyear, &eday, &msec);
        emsec = (int32_t) msec;

        PTB(SetScalarH5A(grp0, "Start Year", H5T_NATIVE_SHORT, (VOIDP) & syear));
        PTB(SetScalarH5A(grp0, "Start Day", H5T_NATIVE_SHORT, (VOIDP) & sday));
        PTB(SetScalarH5A(grp0, "Start Millisec", H5T_STD_U32LE, (VOIDP) & smsec));

        strcpy(buf, ydhmsf(this->meta_l3b.endTime, 'G'));
        PTB(SetScalarH5A(grp0, "End Time", H5T_STRING, buf));

        PTB(SetScalarH5A(grp0, "End Year", H5T_NATIVE_SHORT, (VOIDP) & eyear));
        PTB(SetScalarH5A(grp0, "End Day", H5T_NATIVE_SHORT, (VOIDP) & eday));
        PTB(SetScalarH5A(grp0, "End Millisec", H5T_STD_U32LE, (VOIDP) & emsec));

        strcpy(buf, "degrees North");
        PTB(SetScalarH5A(grp0, "Latitude Units", H5T_STRING, buf));
        strcpy(buf, "degrees East");
        PTB(SetScalarH5A(grp0, "Longitude Units", H5T_STRING, buf));

        PTB(SetScalarH5A(grp0, "Northernmost Latitude", H5T_NATIVE_FLOAT, (VOIDP) & nlat));
        PTB(SetScalarH5A(grp0, "Southernmost Latitude", H5T_NATIVE_FLOAT, (VOIDP) & slat));
        PTB(SetScalarH5A(grp0, "Easternmost Longitude", H5T_NATIVE_FLOAT, (VOIDP) & elon));
        PTB(SetScalarH5A(grp0, "Westernmost Longitude", H5T_NATIVE_FLOAT, (VOIDP) & wlon));

        PTB(SetScalarH5A(grp0, "Data Bins", H5T_STD_U32LE, (VOIDP) & n_data_records));
        float f32 = 100 * ((float) n_data_records) / n_total_bins;
        PTB(SetScalarH5A(grp0, "Percent Data Bins", H5T_NATIVE_FLOAT, (VOIDP) & f32));
    } // CREATE

skip:
    for (size_t i = 0; i < MAXNPROD; i++) {
        if (h5table_id[0][i] != -1)
            H5Dclose(h5table_id[0][i]);
        if (h5table_id[1][i] != -1)
            H5Sclose(h5table_id[1][i]);
        if (h5table_id[2][i] != -1)
            H5Tclose(h5table_id[2][i]);
    }

    H5Gclose(grp0);
    H5Gclose(grp1);

    H5Fclose(h5fid);

    return 0;
}

// Create Bin File

int cdf4_bin::create(const char *l3b_filename, int32_t nrows) {
    int status;

    if(nrows > MAX_ROWS_32BIT)
        is64bit = true;
    else
        is64bit = false;

    if (n_data_prod < 1) {
        cout << "Error - " << __FILE__ << ":" << __LINE__ <<
                " - setProductList needs to be set before calling create." << endl;
        exit(1);
    }

    status = nc_create(l3b_filename, NC_NETCDF4, &ncid);
    check_err(status, __LINE__, __FILE__);
    status = nc_def_grp(ncid, "level-3_binned_data", &grp1);
    check_err(status, __LINE__, __FILE__);

    this->nrows = nrows;
    n_data_records = 0;

    
    binShape = new l3::L3ShapeIsine(nrows);
    totbins = binShape->getNumBins();


    if(is64bit) {
        binList64 = new binListStruct64_cdf4[2 * nrows];
    } else {
        binList = new binListStruct_cdf4[2 * nrows];
    }
        
    binindex_idx = -1;

    if(is64bit)
        defineBinList64_nc(deflate, grp1);
    else
        defineBinList_nc(deflate, grp1);

    defineBinData_nc(deflate, grp1, n_data_prod, product_array);

    if(is64bit)
        defineBinIndex64_nc(deflate, grp1);
    else
        defineBinIndex_nc(deflate, grp1);

    if (hasQual)
        defineQuality_nc(deflate, grp1);

    return 0;
}

// Open Bin File

int cdf4_bin::open(const char* l3b_filename) {
    access_mode = NC_NOWRITE;

    int status;
    int dimid;
    status = nc_open(l3b_filename, access_mode, &ncid);
    if (status != NC_NOERR) {
        cout << l3b_filename << " cannot be opened." << endl;
        exit(1);
    }
    status = nc_inq_ncid(ncid, "level-3_binned_data", &grp1);
    status = nc_inq_dimid(grp1, "binIndexDim", &dimid);
    status = nc_inq_dimlen(grp1, dimid, (size_t *) & nrows);

    status = nc_inq_dimid(grp1, "binListDim", &dimid);
    status = nc_inq_dimlen(grp1, dimid, (size_t *) & n_data_records);

    status = nc_inq_varid(grp1, "BinIndex", &binindex_idx);
    status = nc_inq_varid(grp1, "BinList", &binlist_idx);
    //status = nc_inq_varid(grp1, pname, bindata_id);

    status = nc_inq_varid(grp1, "qual_l3", &binqual_idx);
    if (status == NC_NOERR) {
        hasQual = true;
    } else {
        binqual_idx = -1;
        hasQual = false;
    }

    // Read metadata
    read_l3b_meta_netcdf4(ncid, &meta_l3b);

    int numTables;
    char nam_buf[80];

    // determine if the file is 32 or 64 bit bin numbers
    is64bit = false;
    nc_type xtype;
    status = nc_inq_vartype(grp1, binlist_idx, &xtype);
    if(status == NC_NOERR) {
        int fieldId;
        status = nc_inq_compound_fieldindex(grp1, xtype, "bin_num", &fieldId);
        if(status == NC_NOERR) {
            nc_type fieldType;
            status = nc_inq_compound_fieldtype(grp1, xtype, fieldId, &fieldType);
            if(status == NC_NOERR) {
                if(fieldType == NC_UINT64)
                    is64bit = true;
            }
        }
    }
    
    status = nc_inq_varids(grp1, &numTables, NULL);
    int *varids = (int *) calloc(numTables, sizeof (int));
    status = nc_inq_varids(grp1, &numTables, varids);

    bindata_idx = -1;

    bool first = true;
    for (int i = 0; i < numTables; i++) {
        status = nc_inq_varname(grp1, varids[i], nam_buf);

        if (strcmp(nam_buf, "BinIndex") == 0)
            continue;
        if (strcmp(nam_buf, "BinList") == 0)
            continue;

        strcpy(proddata_name[n_data_prod++], nam_buf);
        //      n_datasets++;
        if (first) {
            bindata_idx = i;
            first = false;
        }
    }
    free(varids);
    
    status = nc_get_att(ncid, NC_GLOBAL, "binning_scheme", nam_buf);
    if (status == NC_NOERR) {
        strcpy(this->meta_l3b.binning_scheme, nam_buf);
    } else {
        cout << "Binning Scheme metadata not found" << endl;
        exit(1);
    }

    binShape = new l3::L3ShapeIsine(nrows);
    totbins = binShape->getNumBins();

    if(is64bit)
        binList64 = new binListStruct64_cdf4[2 * nrows];
    else
        binList = new binListStruct_cdf4[2 * nrows];
        
    return SUCCEED;
}

// Read BinIndex

int cdf4_bin::readBinIndex(int row_num_to_read) {
    int status;
    size_t start = row_num_to_read;

    if(is64bit)
        status = nc_get_var1(grp1, binindex_idx, &start, (void *) &binIndex64);
    else
        status = nc_get_var1(grp1, binindex_idx, &start, (void *) &binIndex);
    if (status != NC_NOERR) {
        printf("-E- %s %d: %s for %s\n",
                __FILE__, __LINE__, nc_strerror(status), "binIndex");
        exit(1);
    }

    return status;
}

// Read BinList

int cdf4_bin::readBinList() {
    // stub
    return 0;
}

// Read BinList

int cdf4_bin::readBinList(int nbins_to_read) {
    int status;
    size_t n2read = nbins_to_read;

    if (lastBinListPtr == binListPtr && lastNumBins == nbins_to_read) {
        binListPtr += nbins_to_read;
        return nbins_to_read;
    }
    lastBinListPtr = binListPtr;
    lastNumBins = nbins_to_read;

    if(nbins_to_read > 0) {
        if(is64bit)
            status = nc_get_vara(grp1, binlist_idx, &binListPtr, &n2read, (void *) binList64);
        else
            status = nc_get_vara(grp1, binlist_idx, &binListPtr, &n2read, (void *) binList);

        if (status != NC_NOERR) {
            printf("-E- %s %d: %s for %s\n",
                    __FILE__, __LINE__, nc_strerror(status), "binList");
            exit(1);
        }

        binListPtr += nbins_to_read;
    }
        
    return nbins_to_read;
}

int cdf4_bin::readBinList(int nbins_to_read, int32_t list_reset_ptr) {
    binListPtr = list_reset_ptr;
    this->readBinList(nbins_to_read);

    return nbins_to_read;
}

// Read Bin File Quality Data

int cdf4_bin::readQual(unsigned char* qual, int32_t nbins_to_read) {
    int status = NC_NOERR;
    size_t n2read = nbins_to_read;

    if(nbins_to_read > 0) {
        status = nc_get_vara(grp1, binqual_idx, &binQualityPtr, &n2read, (void *) qual);
        if (status != NC_NOERR) {
            printf("-E- %s %d: %s in readQual\n",
                __FILE__, __LINE__, nc_strerror(status));
            exit(1);
        }
        binQualityPtr += nbins_to_read;
    }

    return status;
}

int cdf4_bin::readQual(unsigned char* qual, int32_t nbins_to_read,
        int32_t row_num_to_read) {

    int status = NC_NOERR;
    size_t n2read = nbins_to_read;
    size_t start = row_num_to_read;

    if(nbins_to_read > 0) {
        status = nc_get_vara(grp1, binqual_idx, &start, &n2read, (void *) qual);
        if (status != NC_NOERR) {
            printf("-E- %s %d: %s in readQual\n",
                    __FILE__, __LINE__, nc_strerror(status));
            exit(1);
        }
    }

    binQualityPtr = row_num_to_read + nbins_to_read;
    
    return status;
}

// Read Bin File Product Data

int cdf4_bin::readSums(float* sums, int32_t nbins_to_read, int iprod) {
    int status;
    size_t n2read = nbins_to_read;

    int k = 0;
    for (int32_t i = 0; i < n_data_prod; i++) {
        if (active_data_prod[i] == true && (i == iprod || iprod == -1)) {

            status = nc_get_vara(grp1, bindata_idx + i,
                    (size_t *) & binDataPtr, &n2read,
                    (void *) &sums[2 * k * nbins_to_read]);
            if (status != NC_NOERR) {
                printf("-E- %s %d: %s for prod# %d\n",
                        __FILE__, __LINE__, nc_strerror(status), iprod);
                exit(1);
            }

            k++;
        } // if active vdata
    } // vdata loop

    return status;
}

// Read Bin File Product Data for specified bins

int cdf4_bin::readSums(float* sums, int32_t* listOfBins,
        int32_t nbins_to_read, int iprod) {
    int status;
    size_t one = 1;

    int k = 0;
    size_t binPtr;
    for (int32_t i = 0; i < n_data_prod; i++) {
        if (active_data_prod[i] == true && (i == iprod || iprod == -1)) {

            for (int32_t j = 0; j < nbins_to_read; j++) {
                binPtr = listOfBins[j];
                status = nc_get_vara(grp1, bindata_idx + i,
                        (size_t *) & binPtr, &one, (void *) &sums[2 * k]);
                if (status != NC_NOERR) {
                    printf("-E- %s %d: %s for prod# %d\n",
                            __FILE__, __LINE__, nc_strerror(status), iprod);
                    exit(1);
                }
                k++;
            } // j-loop
        } // if active vdata
    } // vdata loop

    return status;
}

int cdf4_bin::write(const char *product_list, hsize_t nwrite, float *outData,
        binListStruct_cdf4* binList) {
    return 0;
}

int cdf4_bin::writeBinList(int32_t nbins_to_write) {
    // Write BinList DataSet
    if(is64bit)
        writeBinList_nc(grp1, nbins_to_write, (VOIDP) binList64);
    else
        writeBinList_nc(grp1, nbins_to_write, (VOIDP) binList);
    
    return 0;
}

int cdf4_bin::writeQual(uint8_t* qual, int32_t nbins_to_write) {
    // Write Quality
    writeQuality_nc(grp1, nbins_to_write, qual);

    return 0;
}

int cdf4_bin::writeSums(float* sums, int32_t nbins_to_write,
        const char *prodname) {
    int iprod;
    bool found = false;

    strcpy(proddata_name[n_data_prod], prodname);
    for (iprod = 0; iprod < n_data_prod; iprod++) {
        if (strcmp(proddata_name[iprod], prodname) == 0) {
            found = true;
            break;
        }
    }
    if (found == false)
        n_data_prod++;

    writeBinData_nc(grp1, nbins_to_write, iprod, sums);

    return 0;
}

int cdf4_bin::close() {
    int retval = 0;
    int status;
    int64_t n_total_bins = 0;
    float slat = +90, nlat = -90, elon = -180, wlon = +180;
    size_t start;
    size_t count;
    
    if (binindex_idx == -1) {
        if(n_data_records > 0) {
            int64_t *binnum_data = (int64_t *) calloc(n_data_records, sizeof (int64_t));
            status = nc_inq_varid(grp1, "BinList", &binlist_idx);
            if (status != NC_NOERR) {
                cout << "Error looking for BinList in output file" << endl;
                exit(1);
            }

            if(is64bit) {
                binListStruct64_cdf4 *binList0;
                binList0 = new binListStruct64_cdf4[n_data_records];
                start = 0;
                count = n_data_records;
                status = nc_get_vara(grp1, binlist_idx, &start, &count, binList0);
                if (status != NC_NOERR) {
                    cout << "Error reading output file binlist records" << endl;
                    exit(1);
                }
                for (size_t i = 0; (int) i < n_data_records; i++) {
                    binnum_data[i] = binList0[i].bin_num;
                }
                delete[] binList0;
            } else {
                binListStruct_cdf4 *binList0;
                binList0 = new binListStruct_cdf4[n_data_records];
                start = 0;
                count = n_data_records;
                status = nc_get_vara(grp1, binlist_idx, &start, &count, binList0);
                if (status != NC_NOERR) {
                    cout << "Error reading output file binlist records" << endl;
                    exit(1);
                }
                for (size_t i = 0; (int) i < n_data_records; i++) {
                    binnum_data[i] = binList0[i].bin_num;
                }
                delete[] binList0;
            }

            for (int32_t i = 0; i < n_data_records - 1; i++) {
                if (binnum_data[i] >= binnum_data[i + 1]) {
                    printf("Improper bin number: %ld at element: %d\n",
                            (long int) binnum_data[i], (int) i);
                }
            }

            int64_t *beg = (int64_t *) calloc(nrows, sizeof (int64_t));
            int32_t *ext = (int32_t *) calloc(nrows, sizeof (int32_t));

            int64_t *numbin_out = (int64_t *) calloc(nrows, sizeof (int64_t));
            int64_t *basebin_out = (int64_t *) calloc(nrows + 1, sizeof (int64_t));

            if ((nrows % 2) == 0) {
                for (int32_t i = 0; i < nrows; i++) {
                    double latbin = (i + 0.5) * (180.0 / (nrows)) - 90.0;
                    numbin_out[i] = (int64_t) (cos(latbin * PI / 180.0)
                            * (2.0 * nrows) + 0.5);
                    n_total_bins += numbin_out[i];
                }
                basebin_out[0] = 1;
            } else {
                int nside = (nrows + 1) / 4;
                totbins = 12 * nside * nside;
                for (int32_t i = 0; i <= nrows / 2; i++) {
                    numbin_out[i] = 4 * (i + 1);
                    if (numbin_out[i] > 4 * nside)
                        numbin_out[i] = 4 * nside;
                    numbin_out[nrows - 1 - i] = numbin_out[i];
                }
                basebin_out[0] = 0;
            }

            for (int32_t i = 1; i <= nrows; i++) {
                basebin_out[i] = basebin_out[i - 1] + numbin_out[i - 1];
            }

            /* Determine beg and extent values for each row */
            /* -------------------------------------------- */
            get_beg_ext(n_data_records, binnum_data, basebin_out, nrows, beg, ext);

            /* Write BinIndex Vdata */
            /* -------------------- */
            if(is64bit) {
                binIndexStruct64_nc* binIndex0 = new binIndexStruct64_nc[nrows];
                for (int32_t i = 0; i < nrows; i++) {
                    binIndex0[i].start_num = basebin_out[i];
                    binIndex0[i].begin = beg[i];
                    binIndex0[i].extent = ext[i];
                    binIndex0[i].max = numbin_out[i];
                }
                writeBinIndex_nc(grp1, nrows, (VOIDP) binIndex0);
                delete [] binIndex0;
            } else {
                binIndexStruct_nc* binIndex0 = new binIndexStruct_nc[nrows];
                for (int32_t i = 0; i < nrows; i++) {
                    binIndex0[i].start_num = basebin_out[i];
                    binIndex0[i].begin = beg[i];
                    binIndex0[i].extent = ext[i];
                    binIndex0[i].max = numbin_out[i];
                }
                writeBinIndex_nc(grp1, nrows, (VOIDP) binIndex0);
                delete [] binIndex0;
            }

            /* Determine NSEW boundries */
            /* ------------------------ */
            int32_t i32 = 0;
            for (int32_t i = 0; i < nrows; i++) {
                float lon, lat;
                if (beg[i] > 0) {
                    lat = ((i + 0.5) / nrows) * 180.0 - 90.0;
                    if (lat > nlat)
                        nlat = lat;
                    if (lat < slat)
                        slat = lat;

                    lon = 360.0 * (beg[i] - basebin_out[i] + 0.5)
                            / numbin_out[i] - 180.0;
                    if (lon < wlon)
                        wlon = lon;

                    lon = 360.0 * (binnum_data[i32 + ext[i] - 1] - basebin_out[i] + 0.5)
                            / numbin_out[i] - 180.0;
                    if (lon > elon)
                        elon = lon;
                }
                i32 += ext[i];
            }

            free(binnum_data);
            free(beg);
            free(ext);
            free(numbin_out);
            free(basebin_out);

            // Write Global Attributes
            // -----------------------
            idDS ds_id;
            ds_id.fid = ncid;
            ds_id.sid = -1;
            ds_id.fftype = DS_NCDF; // FMT_L2NCDF

            meta_l3b.north = nlat;
            meta_l3b.south = slat;
            meta_l3b.east = elon;
            meta_l3b.west = wlon;
            meta_l3b.data_bins = n_data_records;
            meta_l3b.pct_databins = 100 * ((float) n_data_records) / n_total_bins;

            if (nrows == 180) {
                meta_l3b.resolution = str2resolution("1deg");
            } else if (nrows == 360) {
                meta_l3b.resolution = str2resolution("0.5deg");
            } else if (nrows == 540) {
                meta_l3b.resolution = str2resolution("36km");
            } else if (nrows == 720) {
                meta_l3b.resolution = str2resolution("0.25deg");
            } else if (nrows == 1080) {
                meta_l3b.resolution = str2resolution("18km");
            } else if (nrows == 2160) {
                meta_l3b.resolution = str2resolution("9km");
            } else if (nrows == 4320) {
                meta_l3b.resolution = str2resolution("4km");
            } else if (nrows == 8640) {
                meta_l3b.resolution = str2resolution("2km");
            } else if (nrows == 17280) {
                meta_l3b.resolution = str2resolution("1km");
            } else if (nrows == 34560) {
                meta_l3b.resolution = str2resolution("hkm");
            } else if (nrows == 69120) {
                meta_l3b.resolution = str2resolution("qkm");
            } else if (nrows == 172800) {
                meta_l3b.resolution = str2resolution("hqkm");
            } else if (nrows == 345600) {
                meta_l3b.resolution = str2resolution("hhkm");
            } else {
                meta_l3b.resolution = BAD_FLT;
            }

            write_l3b_meta_netcdf4(ds_id, &meta_l3b, is64bit);

        } else { // any data?
            retval = 110; // no bins were filled
        }
    } // CREATE
        
    nc_close(ncid);

    if(is64bit)
        delete[] binList64;
    else
        delete[] binList;

    return retval;
}

// Check whether bin file has qual field

bool cdf4_bin::has_qual() {

    if (binqual_idx != -1) {
        hasQual = true;
    } else {
        hasQual = false;
    }
    return hasQual;
}

/**
 * Get the size of the string needed to hold the names of all the products
 * separated by commas.
 * @return length of product string
 */
int hdf_bin::query() {
    int nchar = 0;
    for (size_t i = 0; (int) i < n_data_prod; i++) {
        nchar += strlen(proddata_name[i]) + 1;
    }
    return nchar;
}

/**
 * Make a string of all the product names separated by commas.
 * @param product_list pointer where string is copied to.
 * @return 0 if all good
 */
int hdf_bin::query(char* product_list) {
    product_list[0] = 0;
    for (int32_t i = 0; i < n_data_prod; i++) {
        strcat(product_list, proddata_name[i]);
        if (i != n_data_prod - 1)
            strcat(product_list, ",");
    }
    return 0;
}

// Query Bin File Product List

int hdf_bin::query(char ***prod_array) {
    // only copy pointers once
    if (product_array[0] == NULL) {
        for (int32_t i = 0; i < n_data_prod; i++) {
            product_array[i] = proddata_name[i];
        }
    }
    if (prod_array != NULL)
        *prod_array = &product_array[0];

    return n_data_prod;
}

int hdf_bin::get_prodname(int iprod, char *prodname) {
    strcpy(prodname, proddata_name[iprod]);
    return 0;
}

void hdf_bin::setProductList(int numProducts, char* prodNames[]) {
    if (numProducts > MAXNPROD) {
        printf("Error - %s:%d - setProductList - Maximum number of products exceeded",
                __FILE__, __LINE__);
        exit(1);
    }

    n_data_prod = numProducts;
    for (int i = 0; i < numProducts; i++) {
        if (strlen(prodNames[i]) >= 80) {
            printf("Error - %s:%d - setProductList - Product name too long \"%s\"",
                    __FILE__, __LINE__, prodNames[i]);
            exit(1);
        }
        strcpy(proddata_name[i], prodNames[i]);
        product_array[i] = proddata_name[i];
    }
}

/**
 * Get the name of a product given the index.
 * @param prodIndex index of the product
 * @return internal pointer to the product name
 */
const char* hdf_bin::getProdName(int prodIndex) const {
    if (prodIndex < 0 || prodIndex >= MAXNPROD)
        return NULL;
    return proddata_name[prodIndex];
}

/**
 * Get the index for the product name
 * @param prodname name of the product desired
 * @return product index, or -1 if not found
 */
int hdf_bin::getProdIndex(const char *prodname) const {
    int iprod = -1;
    for (int i = 0; i < n_data_prod; i++) {
        if (strcmp(proddata_name[i], prodname) == 0) {
            iprod = i;
            break;
        }
    }
    return iprod;
}

/**
 * Get the name of a active product given the index.
 * @param prodIndex index of the product
 * @return internal pointer to the product name
 */
const char* hdf_bin::getActiveProdName(int prodIndex) const {
    if (prodIndex < 0 || prodIndex >= n_active_prod)
        return NULL;
    int activeProd = 0;
    int i;
    for (i = 0; i < n_data_prod; i++) {
        if (active_data_prod[i]) {
            if (activeProd == prodIndex) {
                return proddata_name[i];
            }
            activeProd++;
        }
    }
    return NULL;
}

int hdf_bin::copymeta(int32_t nfiles, Hdf::hdf_bin *input_binfile[]) {
    strcpy(this->meta_l3b.sensor, input_binfile[0]->meta_l3b.sensor);
    strcpy(this->meta_l3b.sensor_name, input_binfile[0]->meta_l3b.sensor_name);
    strcpy(this->meta_l3b.data_center, input_binfile[0]->meta_l3b.data_center);
    strcpy(this->meta_l3b.mission, input_binfile[0]->meta_l3b.mission);
    strcpy(this->meta_l3b.mission_char,
            input_binfile[0]->meta_l3b.mission_char);
    strcpy(this->meta_l3b.sensor_char, input_binfile[0]->meta_l3b.sensor_char);
    strcpy(this->meta_l3b.pversion, input_binfile[0]->meta_l3b.pversion);
    strcpy(this->meta_l3b.prod_type, input_binfile[0]->meta_l3b.prod_type);

    int32_t sorbit = input_binfile[0]->meta_l3b.start_orb;
    for (int32_t ifile = 1; ifile < nfiles; ifile++) {
        int32_t i32 = input_binfile[ifile]->meta_l3b.start_orb;
        if (i32 < sorbit)
            sorbit = i32;
    }
    this->meta_l3b.start_orb = sorbit;

    int32_t eorbit = input_binfile[0]->meta_l3b.end_orb;
    for (int32_t ifile = 1; ifile < nfiles; ifile++) {
        int32_t i32 = input_binfile[ifile]->meta_l3b.end_orb;
        if (i32 > eorbit)
            eorbit = i32;
    }
    this->meta_l3b.end_orb = eorbit;

    strcpy(this->meta_l3b.flag_names, input_binfile[0]->meta_l3b.flag_names);

    double stime;
    double etime;

    int32_t sfile = 0;
    stime = input_binfile[0]->meta_l3b.startTime;
    for (int32_t ifile = 1; ifile < nfiles; ifile++) {
        if (input_binfile[ifile]->meta_l3b.startTime < stime) {
            stime = input_binfile[ifile]->meta_l3b.startTime;
            sfile = ifile;
        }
    }
    this->meta_l3b.startTime = input_binfile[sfile]->meta_l3b.startTime;

    int32_t efile = 0;
    etime = input_binfile[0]->meta_l3b.endTime;
    for (int32_t ifile = 1; ifile < nfiles; ifile++) {
        if (input_binfile[ifile]->meta_l3b.endTime > etime) {
            etime = input_binfile[ifile]->meta_l3b.endTime;
            efile = ifile;
        }
    }
    this->meta_l3b.endTime = input_binfile[efile]->meta_l3b.endTime;
#if(0)
    this->meta_l3b.bin_syear = input_binfile[sfile]->meta_l3b.bin_syear;
    this->meta_l3b.bin_sday = input_binfile[sfile]->meta_l3b.bin_sday;
    this->meta_l3b.bin_eyear = input_binfile[efile]->meta_l3b.bin_eyear;
    this->meta_l3b.bin_eday = input_binfile[efile]->meta_l3b.bin_eday;

    this->meta_l3b.syear = input_binfile[sfile]->meta_l3b.syear;
    this->meta_l3b.sday = input_binfile[sfile]->meta_l3b.sday;
    this->meta_l3b.smsec = input_binfile[sfile]->meta_l3b.smsec;

    this->meta_l3b.eyear = input_binfile[efile]->meta_l3b.eyear;
    this->meta_l3b.eday = input_binfile[efile]->meta_l3b.eday;
    this->meta_l3b.emsec = input_binfile[efile]->meta_l3b.emsec;
#endif
    strcpy(this->meta_l3b.station, input_binfile[0]->meta_l3b.station);
    strcpy(this->meta_l3b.units, input_binfile[0]->meta_l3b.units);

    this->meta_l3b.station_lat = input_binfile[0]->meta_l3b.station_lat;
    this->meta_l3b.station_lon = input_binfile[0]->meta_l3b.station_lon;

    strcpy(this->meta_l3b.binning_scheme,
            input_binfile[sfile]->meta_l3b.binning_scheme);

    strcpy(this->meta_l3b.doi, input_binfile[0]->meta_l3b.doi);
    strcpy(this->meta_l3b.keywords, input_binfile[0]->meta_l3b.keywords);
    
    return 0;
}

}
