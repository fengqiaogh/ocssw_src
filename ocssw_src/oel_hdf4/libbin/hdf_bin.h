#ifndef hdf_bin_h
#define hdf_bin_h

#pragma GCC diagnostic ignored "-Wpadded"
#include "bin_util.h"
#include "meta_l3b.h"
#include "L3Shape.h"

namespace Hdf {
class hdf_bin;

hdf_bin* openBinObject(const char* binFileName);

class hdf_bin {
protected:
    int32_t n_data_prod;
    char proddata_name[MAXNPROD][80];
    char *product_array[MAXNPROD];

    l3::L3Shape* binShape;
    
    size_t binListPtr;
    size_t lastBinListPtr;
    int lastNumBins;

public:
    hdf_bin();
    virtual ~hdf_bin();

    virtual int query();
    virtual int query(char* product_list);
    virtual int query(char ***prod_array);
    virtual int get_prodname(int iprod, char *prodname);
    virtual void setProductList(int numProducts, char* prodNames[]);

    virtual const char* getProdName(int prodNum) const;
    virtual int getProdIndex(const char *prodname) const;
    virtual const char* getActiveProdName(int prodNum) const;

    virtual int read(char* product_list);

    virtual int64_t get_beg() = 0;
    virtual int get_ext() = 0;

    virtual int open(const char* l3b_filename) = 0;
    virtual int create(const char* l3b_filename, int32_t nrows) = 0;

    virtual int readBinIndex(int row_num_to_read) = 0;
    virtual int readBinList(int nbins_to_read) = 0;
    virtual int readBinList(int nbins_to_read, int32_t list_reset_ptr) = 0;
    virtual int readBinList() = 0;

    virtual int writeBinList(int32_t nbins_to_write) = 0;
    virtual int readQual(uint8_t* qual, int32_t nbins_to_read) = 0;
    virtual int readQual(uint8_t* qual, int32_t nbins_to_read,
            int32_t row_num_to_read) = 0;

    /**
     * Read Bin File Product Data
     * @param sums array to place sum and sumSquares for each product
     * @param nbins_to_read number of consecutive bins to read
     * @param iprod product index or -1 to read all active products set by read()
     * @return 0 if OK
     */
    virtual int readSums(float* sums, int32_t nbins_to_read, int iprod) = 0;
    virtual int readSums(float* sums, int32_t* listOfBins,
            int32_t nbins_to_read, int iprod) = 0;
    virtual int writeQual(uint8_t* qual, int32_t nbins_to_write) = 0;
    virtual int writeSums(float* sums, int32_t nbins_to_write,
            const char *prodname) = 0;

    virtual int64_t get_numbin(int irow) {
        return binShape->getNumCols(irow);
    };
    
    virtual int64_t get_basebin(int irow) {
        return binShape->getBaseBin(irow);
    };

    virtual void bin2latlon(int64_t bin_num, float &lat, float &lon) {
        binShape->bin2latlon(bin_num, lat, lon);
    }
    
    virtual int64_t get_bin_num(int kbin) = 0;
    virtual int get_nobs(int kbin) = 0;
    virtual int get_nscenes(int kbin) = 0;
    virtual float get_weights(int kbin) = 0;
    virtual float get_time_rec(int kbin) = 0;

    virtual hid_t get_index_table() = 0;
    virtual hid_t get_list_table() = 0;
    virtual hid_t get_data_table(int i) = 0;

    virtual int clear_binlist() = 0;
    virtual int copy_binlist(int src, int dest) = 0;
    virtual int set_bin_num(int offset, int64_t bin_num) = 0;
    virtual int inc_nobs(int offset, int nobs) = 0;
    virtual int set_nobs(int offset, int nobs) = 0;
    virtual int inc_nscenes(int offset, int nscenes) = 0;
    virtual int set_nscenes(int offset, int nscenes) = 0;
    virtual int inc_weights(int offset, float weights) = 0;
    virtual int set_weights(int offset, float weights) = 0;
    virtual int inc_time_rec(int offset, float time_rec) = 0;
    virtual bool has_qual() = 0;
    virtual int setDataPtr(int nbins_to_read) = 0;
    virtual int setDataPtrAbsolute(int32_t recordNum) = 0;
    virtual int incNumRec(int n_write) = 0;
    virtual int close() = 0;

    virtual int32_t nprod() {
        return n_data_prod;
    }

    virtual int32_t get_list_ptr() {
        return binListPtr;
    }

    virtual int copymeta(int32_t nfiles, Hdf::hdf_bin *input_binfile[]);

    int64_t totbins;
    int32_t nrows;

    bool active_data_prod[MAXNVDATA];
    int32_t n_data_records;
    int32_t n_active_prod;
    bool isHDF5;
    bool isCDF4;
    bool hasQual;
    bool hasNoext;

    uint32_t deflate;

    meta_l3bType meta_l3b;
};

class hdf4_bin : public hdf_bin {
    int32_t file_id;
    int32_t sd_id;
    int32_t vg_id;
    int32_t access_mode;

    int32_t vdata_id[MAXNVDATA];
    int32_t seagrid_idx;
    int32_t binindex_idx;
    int32_t binlist_idx;
    int32_t bin_ptr;
    //int32_t list_reset_ptr;
    int32_t bindata_idx;
    int32_t binqual_idx;

    int32_t aid[MAXNPROD];

    binIndexStruct binIndex;
    binListStruct *binList;

public:
    hdf4_bin();
    virtual ~hdf4_bin();

    int64_t get_beg() {
        return binIndex.beg;
    }

    int get_ext() {
        return binIndex.ext;
    }

    int create(const char* l3b_filename, int32_t nrows);
    int open(const char* l3b_filename);

    int readBinIndex(int row_num_to_read);
    
    using hdf_bin::read;
    int read(float* data, binListStruct* binList);
    int read(float* data, float* var, binListStruct* binList);
    int read(float* data, binListStruct* binList, int nbins_to_read);
    int read(float* data, float* var, binListStruct* binList,
            int nbins_to_read);

    int readBinList(int nbins_to_read);
    int readBinList();
    int readBinList(int nbins_to_read, int32_t list_reset_ptr);
    int readQual(uint8_t* qual, int32_t nbins_to_read);
    int readQual(uint8_t* qual, int32_t nbins_to_read, int32_t row_num_to_read);

    int readSums(float* sums, int32_t nbins_to_read, int iprod);
    int readSums(float* sums, int32_t* listOfBins, int32_t nbins_to_read,
            int iprod);

    int64_t get_bin_num(int kbin) {
        return binList[kbin].bin_num;
    }

    int get_nobs(int kbin) {
        return binList[kbin].nobs;
    }

    int get_nscenes(int kbin) {
        return binList[kbin].nscenes;
    }

    float get_weights(int kbin) {
        return binList[kbin].weights;
    }

    float get_time_rec(int kbin) {
        return binList[kbin].time_rec;
    }

    int set_bin_num(int offset, int64_t bin_num) {
        binList[offset].bin_num = bin_num;
        return 0;
    }

    int inc_nobs(int offset, int nobs) {
        binList[offset].nobs += nobs;
        return 0;
    }

    int set_nobs(int offset, int nobs) {
        binList[offset].nobs = nobs;
        return 0;
    }

    int inc_nscenes(int offset, int nscenes) {
        binList[offset].nscenes += nscenes;
        return 0;
    }

    int set_nscenes(int offset, int nscenes) {
        binList[offset].nscenes = nscenes;
        return 0;
    }

    int inc_weights(int offset, float weights) {
        binList[offset].weights += weights;
        return 0;
    }

    int set_weights(int offset, float weights) {
        binList[offset].weights = weights;
        return 0;
    }

    int inc_time_rec(int offset, float time_rec) {
        binList[offset].time_rec += time_rec;
        return 0;
    }

    int clear_binlist() {
        memset(binList, 0, 2 * nrows * sizeof (binListStruct));
        return 0;
    }

    int copy_binlist(int src, int dest) {
        memcpy(&binList[dest], &binList[src], sizeof (binListStruct));
        return 0;
    }

    int write(char *product_list, int32_t nwrite, float *data,
            binListStruct* binList);
    int writeBinList(int32_t nbins_to_write);
    int writeQual(uint8_t* qual, int32_t nbins_to_write);
    int writeSums(float* sums, int32_t nbins_to_write, const char *prodname);

    int copy(char *product_list, int32_t nwrite, int32_t*binsToCopy,
            Hdf::binListStruct *inBinList,
            Hdf::hdf4_bin *input_binfile);

    int close();

    bool has_qual();

    int setDataPtr(int nbins_to_read) {
        return 0;
    }
    int setDataPtrAbsolute(int32_t recordNum);

    int incNumRec(int n_write) {
        return 0;
    }

    hid_t get_index_table() {
        return 0;
    }

    hid_t get_list_table() {
        return 0;
    }

    hid_t get_data_table(int i) {
        return 0;
    }

    int32_t noext;

};

class hdf5_bin : public hdf_bin {
    hid_t h5fid;
    hid_t grp0, grp1;
    hid_t access_mode;

    int32_t n_datasets;
    hid_t h5table_id[3][MAXNPROD];

    int32_t bindata_idx;
    int32_t binlist_idx;
    int32_t binindex_idx;

    binIndexStruct binIndex;
    binListStruct_hdf5 *binList;

public:
    hdf5_bin();
    ~hdf5_bin();

    int64_t get_beg() {
        return binIndex.beg;
    }

    int get_ext() {
        return binIndex.ext;
    }

    int open(const char* l3b_filename);
    int create(const char* l3b_filename, int32_t nrows);

    int readBinIndex(int row_num_to_read);
    int readBinList(int nbins_to_read);
    int readBinList(int nbins_to_read, int32_t list_reset_ptr);
    int readBinList();

    int readQual(unsigned char* qual, int32_t nbins_to_read);
    int readQual(unsigned char* qual, int32_t nbins_to_read,
            int32_t row_num_to_read);
    int readSums(float* sums, int32_t nbins_to_read, int iprod);
    int readSums(float* sums, int32_t* listOfBins, int32_t nbins_to_read,
            int iprod);

    int writeBinList(int32_t nbins_to_write);

    int writeSums(float* sums, int32_t nbins_to_write, const char *prodname);

    int writeQual(uint8_t* qual, int32_t nbins_to_write) {
        return 0;
    }

    int write(const char *product_list, hsize_t nwrite, float *data,
            binListStruct_hdf5* binList);
    int close();

    bool has_qual() {
        return false;
    }

    hid_t get_index_table() {
        return h5table_id[0][binindex_idx];
    }

    hid_t get_list_table() {
        return h5table_id[0][binlist_idx];
    }

    hid_t get_data_table(int i) {
        return h5table_id[0][bindata_idx + i];
    }

    hid_t get_grp0() {
        return grp0;
    }

    int64_t get_bin_num(int kbin) {
        return binList[kbin].bin_num;
    }

    int get_nobs(int kbin) {
        return binList[kbin].nobs;
    }

    int get_nscenes(int kbin) {
        return binList[kbin].nscenes;
    }

    float get_weights(int kbin) {
        return binList[kbin].weights;
    }

    float get_time_rec(int kbin) {
        return 0.0;
    }

    int set_bin_num(int offset, int64_t bin_num) {
        binList[offset].bin_num = bin_num;
        return 0;
    }

    int inc_nobs(int offset, int nobs) {
        binList[offset].nobs += nobs;
        return 0;
    }

    int set_nobs(int offset, int nobs) {
        binList[offset].nobs = nobs;
        return 0;
    }

    int inc_nscenes(int offset, int nscenes) {
        binList[offset].nscenes += nscenes;
        return 0;
    }

    int set_nscenes(int offset, int nscenes) {
        binList[offset].nscenes = nscenes;
        return 0;
    }

    int inc_weights(int offset, float weights) {
        binList[offset].weights += weights;
        return 0;
    }

    int set_weights(int offset, float weights) {
        binList[offset].weights = weights;
        return 0;
    }

    int inc_time_rec(int offset, float time_rec) {
        return 0;
    }

    int clear_binlist() {
        memset(binList, 0, 2 * nrows * sizeof (binListStruct_hdf5));
        return 0;
    }

    int copy_binlist(int src, int dest) {
        memcpy(&binList[dest], &binList[src], sizeof (binListStruct_hdf5));
        return 0;
    }

    int setDataPtr(int nbins_to_read) {
        binDataPtr += nbins_to_read;
        return 0;
    }

    int setDataPtrAbsolute(int32_t recordNum) {
        binDataPtr = recordNum;
        return 0;
    }

    int incNumRec(int n_write) {
        n_data_records += n_write;
        return 0;
    }

    hsize_t binDataPtr;
};

class cdf4_bin : public hdf_bin {
    bool is64bit;
    int ncid;
    int grp0, grp1;
    int access_mode;

    int n_datasets;

    int bindata_idx;
    int binqual_idx;
    int binlist_idx;
    int binindex_idx;

    size_t binDataPtr;
    size_t binQualityPtr;

    binIndexStruct_cdf4 binIndex;
    binListStruct_cdf4 *binList;

    binIndexStruct64_cdf4 binIndex64;
    binListStruct64_cdf4 *binList64;

public:
    cdf4_bin();
    ~cdf4_bin();

    int64_t get_beg() {
        if(is64bit)
            return binIndex64.beg;
        else
            return binIndex.beg;
    }

    int get_ext() {
        if(is64bit)
            return binIndex64.ext;
        else
            return binIndex.ext;
    }

    int open(const char* l3b_filename);
    int create(const char* l3b_filename, int32_t nrows);

    int readBinIndex(int row_num_to_read);
    int readBinList(int nbins_to_read);
    int readBinList(int nbins_to_read, int32_t list_reset_ptr);
    int readBinList();

    int readQual(unsigned char* qual, int32_t nbins_to_read);
    int readQual(unsigned char* qual, int32_t nbins_to_read,
            int32_t row_num_to_read);
    int readSums(float* sums, int32_t nbins_to_read, int iprod);
    int readSums(float* sums, int32_t* listOfBins, int32_t nbins_to_read,
            int iprod);

    int writeBinList(int32_t nbins_to_write);

    int writeSums(float* sums, int32_t nbins_to_write, const char *prodname);
    int writeQual(uint8_t* qual, int32_t nbins_to_write);

    int write(const char *product_list, hsize_t nwrite, float *data,
            binListStruct_cdf4* binList);
    int close();

    bool has_qual();

    hid_t get_grp0() {
        return grp0;
    }

    hid_t get_index_table() {
        return 0;
    }

    hid_t get_list_table() {
        return 0;
    }

    hid_t get_data_table(int i) {
        return 0;
    }

    int64_t get_bin_num(int kbin) {
        if(is64bit)
            return binList64[kbin].bin_num;
        else
            return binList[kbin].bin_num;
    }

    int get_nobs(int kbin) {
        if(is64bit)
            return binList64[kbin].nobs;
        else
            return binList[kbin].nobs;
    }

    int get_nscenes(int kbin) {
        if(is64bit)
            return binList64[kbin].nscenes;
        else
            return binList[kbin].nscenes;
    }

    float get_weights(int kbin) {
        if(is64bit)
            return binList64[kbin].weights;
        else
            return binList[kbin].weights;
    }

    float get_time_rec(int kbin) {
        if(is64bit)
            return binList64[kbin].time_rec;
        else
            return binList[kbin].time_rec;
    }

    int set_bin_num(int offset, int64_t bin_num) {
        if(is64bit)
            binList64[offset].bin_num = bin_num;
        else
            binList[offset].bin_num = bin_num;
        return 0;
    }

    int inc_nobs(int offset, int nobs) {
        if(is64bit)
            binList64[offset].nobs += nobs;
        else
            binList[offset].nobs += nobs;
        return 0;
    }

    int set_nobs(int offset, int nobs) {
        if(is64bit)
            binList64[offset].nobs = nobs;
        else
            binList[offset].nobs = nobs;
        return 0;
    }

    int inc_nscenes(int offset, int nscenes) {
        if(is64bit)
            binList64[offset].nscenes += nscenes;
        else
            binList[offset].nscenes += nscenes;
        return 0;
    }

    int set_nscenes(int offset, int nscenes) {
        if(is64bit)
            binList64[offset].nscenes = nscenes;
        else
            binList[offset].nscenes = nscenes;
        return 0;
    }

    int inc_weights(int offset, float weights) {
        if(is64bit)
            binList64[offset].weights += weights;
        else
            binList[offset].weights += weights;
        return 0;
    }

    int set_weights(int offset, float weights) {
        if(is64bit)
            binList64[offset].weights = weights;
        else
            binList[offset].weights = weights;
        return 0;
    }

    int inc_time_rec(int offset, float time_rec) {
        if(is64bit)
            binList64[offset].time_rec += time_rec;
        else
            binList[offset].time_rec += time_rec;
        return 0;
    }

    int clear_binlist() {
        if(is64bit)
            memset(binList64, 0, 2 * nrows * sizeof (binListStruct64_cdf4));
        else
            memset(binList, 0, 2 * nrows * sizeof (binListStruct_cdf4));
        return 0;
    }

    int copy_binlist(int src, int dest) {
        if(is64bit)
            memcpy(&binList64[dest], &binList64[src], sizeof (binListStruct64_cdf4));
        else
            memcpy(&binList[dest], &binList[src], sizeof (binListStruct_cdf4));
        return 0;
    }

    int setDataPtr(int nbins_to_read) {
        binDataPtr += nbins_to_read;
        return 0;
    }

    int setDataPtrAbsolute(int32_t recordNum) {
        binDataPtr = recordNum;
        return 0;
    }

    int incNumRec(int n_write) {
        n_data_records += n_write;
        return 0;
    }

};

}

#endif
