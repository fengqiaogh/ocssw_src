#include <L3File.h>
#include <genutils.h>
#include <boost/algorithm/string.hpp>

#include "L3ShapeIsine.h"

using namespace std;
//using namespace boost::geometry;

namespace l3 {

//----------------------------------------------------------------------------
// L3Bin
//----------------------------------------------------------------------------

L3Bin::L3Bin(int32_t numProducts) {
    internalInit(numProducts);
}

L3Bin::L3Bin(const L3Bin &bin) {
    internalCopy(bin);
}

L3Bin::~L3Bin() {
    free(sums);
    free(sumSquares);
}

void L3Bin::clear() {
    binNum = -1;
    recordNum = -1;
    nobs = 0;
    nscenes = 0;
    timeRec = 0.0;
    weights = 0.0;
    quality = qualityUnused;
    for (int i = 0; i < numProducts; i++) {
        sums[i] = 0.0;
        sumSquares[i] = 0.0;
    }
}

/**
 * allocate memory for for the given number of products.  Frees the old
 * memory if needed.
 * @param numProducts
 */
void L3Bin::internalAllocate(int32_t numProducts) {
    if (numProducts != this->numProducts) {
        if (sums)
            free(sums);
        if (sumSquares)
            free(sumSquares);
        sums = (float*) allocateMemory(numProducts * sizeof (float),
                "sums for L3Bin");
        sumSquares = (float*) allocateMemory(numProducts * sizeof (float),
                "sumSquares for L3Bin");
        this->numProducts = numProducts;
    }
}

/**
 * init all the values in this instance.
 * @param numProducts number of products to allocate space for
 */
void L3Bin::internalInit(int32_t numProducts) {
    this->numProducts = 0;
    sums = sumSquares = NULL;
    internalAllocate(numProducts);
    clear();
}

/**
 * Copy everything form another instance assuming this instance has not allocated
 * any memory yet.
 * @param bin instance to copy from
 */
void L3Bin::internalCopy(const L3Bin &bin) {
    binNum = bin.binNum;
    recordNum = bin.recordNum;
    nobs = bin.nobs;
    nscenes = bin.nscenes;
    quality = bin.quality;
    timeRec = bin.timeRec;
    weights = bin.weights;

    internalAllocate(bin.numProducts);
    memcpy(sums, bin.sums, numProducts * sizeof (float));
    memcpy(sumSquares, bin.sumSquares, numProducts * sizeof (float));
}

/**
 * Add everything form another instance into this instance
 * @param bin instance to add
 */
void L3Bin::internalAdd(const L3Bin &bin) {
    if (numProducts != bin.numProducts) {
        setNumProducts(bin.numProducts);
    }
    if (bin.quality != qualityUnused) {
        if (bin.quality < quality) {
            clear();
            quality = bin.quality;
        }
    }
    binNum = bin.binNum; // just copy the binNum
    nobs += bin.nobs;
    nscenes += bin.nscenes;
    timeRec += bin.timeRec;
    weights += bin.weights;
    for (int i = 0; i < numProducts; i++) {
        sums[i] += bin.sums[i];
        sumSquares[i] += bin.sumSquares[i];
    }
}

void L3Bin::checkProdId(int32_t prodId) const {
    if (prodId < 0) {
        fprintf(stderr,
                "-E- L3Bin::checkProdId - prodId (%d) must be positive\n",
                prodId);
        exit(EXIT_FAILURE);
    }
    if (prodId >= numProducts) {
        fprintf(stderr,
                "-E- L3Bin::checkProdId - prodId (%d) is greater than the numProducts (%d)\n",
                prodId, numProducts);
        exit(EXIT_FAILURE);
    }
}

L3Bin& L3Bin::operator=(const L3Bin &bin) {
    internalCopy(bin);
    return *this;
}

L3Bin& L3Bin::operator+=(const L3Bin &bin) {
    internalAdd(bin);
    return *this;
}

void L3Bin::addWeighted(const L3Bin &bin, float weighting) {
    if (bin.quality != qualityUnused) {
        if (bin.quality < quality) {
            clear();
            quality = bin.quality;
        }
    }
    nobs += bin.nobs;
    nscenes += bin.nscenes;
    timeRec += bin.timeRec;
    weights += bin.weights * weighting;
    for (int i = 0; i < numProducts; i++) {
        sums[i] += bin.sums[i] * weighting;
        sumSquares[i] += bin.sumSquares[i] * weighting;
    }
}

/**
 * set the number of products this bin can hold
 * @param numProducts
 */
void L3Bin::setNumProducts(int32_t numProducts) {
    if (numProducts < 1) {
        fprintf(stderr,
                "-E- L3Bin::setNumProducts - numProducts (%d) must be at least 1\n",
                numProducts);
        exit(EXIT_FAILURE);
    }
    internalAllocate(numProducts);
    clear();
}

/**
 * get the number of products in this bin
 * @return number of products
 */
int32_t L3Bin::getNumProducts() const {
    return numProducts;
}

/**
 * get the bin number for this bin
 * @return bin number
 */
int64_t L3Bin::getBinNum() const {
    return binNum;
}

/**
 * get the number of observations for this bin
 * @return number observations
 */
int32_t L3Bin::getNobs() const {
    return nobs;
}

/**
 * get the numbers of scenes contributing to this bin
 * @return number of scenes
 */
int32_t L3Bin::getNscenes() const {
    return nscenes;
}

/**
 * get the record number in the file for this bin
 * @return record number
 */
int64_t L3Bin::getRecordNum() const {
    return recordNum;
}

/**
 * get the observation time.  sec since 1993 (TAI93)
 * @return TAI93 time
 */
float L3Bin::getObsTime() const {
    return timeRec / nobs;
}

/**
 * get the weights for this bin.  Used to calculate the mean.
 * @return weights
 */
float L3Bin::getWeights() const {
    return weights;
}

/**
 * get the SST quality factor.  255 is unused, 0 is best, 2 is poor
 * @return SST quality
 */
uint8_t L3Bin::getQuality() const {
    return quality;
}

/**
 * get this bin's sum
 * @param prodId index of the product you want
 * @return the sum
 */
float L3Bin::getSum(int32_t prodId) const {
    checkProdId(prodId);
    return sums[prodId];
}

/**
 * get the sum of the squares for this product
 * @param prodId index of the product you want
 * @return sum of squares
 */
float L3Bin::getSumSquares(int32_t prodId) const {
    checkProdId(prodId);
    return sumSquares[prodId];
}

/**
 * calculate the mean for this product
 * @param prodId index of the product you want
 * @return the average value
 */
float L3Bin::getMean(int32_t prodId) const {
    checkProdId(prodId);
    return sums[prodId] / weights;
}

/**
 * calculate the variance for this product
 * @param prodId index of the product you want
 * @return the variance value
 */
float L3Bin::getVariance(int32_t prodId) const {
    checkProdId(prodId);
    if (weights * weights > nscenes) {
        double val = sums[prodId] / weights;
        val = (sumSquares[prodId] / weights) - (val * val);
        return val * weights * weights / (weights * weights - nscenes);
    }
    return 0.0;
}

/**
 * calculate the standard deviation for this product
 * @param prodId index of the product you want
 * @return the standard deviation
 */
float L3Bin::getStdev(int32_t prodId) const {
    return sqrt(getVariance(prodId));
}


//----------------------------------------------------------------------------
// L3Row
//----------------------------------------------------------------------------

L3Row::L3Row(int32_t row, int64_t numBins, int32_t numProducts) {
    this->row = row;
    this->numBins = numBins;
    this->numProducts = numProducts;
    binArray.resize(numBins);
    for (int64_t i = 0; i < numBins; i++) {
        binArray[i] = new L3Bin(numProducts);
    }
    lastBin = numBins / 2;
}

L3Row::~L3Row() {
    setNumBins(0);
}

int32_t L3Row::getRow() const {
    return row;
}

int32_t L3Row::getNumProducts() const {
    return numProducts;
}

int64_t L3Row::getNumBins() const {
    return numBins;
}

void L3Row::setNumBins(int64_t numBins) {
    this->numBins = numBins;
    while (binArray.size() > (uint64_t) numBins) {
        delete binArray.back();
        binArray.pop_back();
    }
    while (binArray.size() < (uint64_t) numBins) {
        binArray.push_back(new L3Bin(numProducts));
    }
    lastBin = numBins / 2;
}

L3Bin* L3Row::getBin(int64_t binNum) {
    L3Bin* l3Bin;
    int64_t hi, mid, lo;
    int64_t tmpBinNum;
    int64_t tmp;

    if (numBins == 0)
        return NULL;

    // start at the bin we found last time
    lo = 0;
    hi = numBins - 1;
    mid = lastBin;

    while (lo <= hi) {
        l3Bin = binArray[mid];
        tmpBinNum = l3Bin->getBinNum();
        if (tmpBinNum == binNum) {
            lastBin = mid + 1;
            if (lastBin >= numBins)
                lastBin--;
            if (l3Bin->nscenes == 0)
                return NULL;
            else
                return l3Bin;
        }
        tmp = mid + (binNum - tmpBinNum);
        if (binNum > tmpBinNum) {
            lo = mid + 1;
            if (hi > tmp)
                hi = tmp;
        } else {
            hi = mid - 1;
            if (lo < tmp)
                lo = tmp;
        }
        mid = (lo + hi) / 2;
    }
    return NULL;
}

/**
 * Get the bin in this row by index.
 * @param index to bin in this row that was actually in the file.
 * @return pointer to bin object at given index or NULL if index out of range
 */
L3Bin* L3Row::getBinByIndex(int32_t index) {
    if (index < 0 || index >= numBins)
        return NULL;
    return binArray[index];
}
/**
 * @brief Return the bin number for the first valid bin in the row
 *
 * @return int32_t bin_number
 */
int64_t L3Row::getFirstBinNum() {
    if (binArray.empty()){
        return -1;
    }
    return binArray[0]->getBinNum();
}

/**
 * @brief Return the bin number for the last valid bin in the row
 *
 * @return int32_t bin_number
 */
int64_t L3Row::getLastBinNum() {
    if (binArray.empty()){
        return -1;
    }
    return binArray.back()->getBinNum();
}

//----------------------------------------------------------------------------
// L3File
//----------------------------------------------------------------------------

L3File::L3File() {
    shape = NULL;
    binObj = NULL;
    numCacheRows = 500;
    baseRecord = NULL;
    extentbin = NULL;
    sumBuffer = NULL;
    qualityBuffer = NULL;
    prodMap = NULL;
    numProds = 0;
}

L3File::~L3File() {
    close();
    if (shape)
        delete shape;
    if (binObj)
        delete binObj;
}

void L3File::clearCache() {
    L3Row* l3Row;
    while (!rowList.empty()) {
        l3Row = rowList.back();
        delete l3Row;
        rowList.pop_back();
    }
}

void L3File::setNumCacheRows(int32_t numRows) {
    numCacheRows = numRows;
    L3Row* l3Row;
    while (rowList.size() > (uint32_t) numRows) {
        l3Row = rowList.back();
        delete l3Row;
        rowList.pop_back();
    }
}

bool L3File::open(const char* fileName) {
    close();
    int tmpWantVerbose = want_verbose;
    want_verbose = 0;
    binObj = Hdf::openBinObject(fileName);
    want_verbose = tmpWantVerbose;
    if (binObj == NULL)
        return false;
    shape = new L3ShapeIsine(binObj->nrows);
    initRecordLookup();
    return true;
}

void L3File::close() {
    if (shape) {
        delete shape;
        shape = NULL;
    }
    if (binObj) {
        delete binObj;
        binObj = NULL;
    }
    if (sumBuffer) {
        free(sumBuffer);
        sumBuffer = NULL;
    }
    if (qualityBuffer) {
        free(qualityBuffer);
        qualityBuffer = NULL;
    }
    if (baseRecord) {
        free(baseRecord);
        baseRecord = NULL;
    }
    if (extentbin) {
        free(extentbin);
        extentbin = NULL;
    }
    clearCache();
    if (prodMap) {
        free(prodMap);
        prodMap = NULL;
    }
}

meta_l3bType* L3File::getMetaData() {
    if (binObj) {
        return &binObj->meta_l3b;
    }
    return NULL;
}

/**
 * init the internal data structures to lookup the record index
 * @return status 0 if good, -1 if bad
 */
int L3File::initRecordLookup() {

    if (binObj->nrows == 0) {
        fprintf(stderr,
                "-E- %s %d: initRecordLookup: Number of rows needs to be set and bin grid must be initialized first.\n",
                __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }

    // remove memory if it exists
    if (baseRecord) {
        free(baseRecord);
    }
    baseRecord = (int64_t*) allocateMemory(binObj->nrows * sizeof (int64_t),
            "baseRecord");
    if (extentbin) {
        free(extentbin);
    }
    extentbin = (int32_t*) allocateMemory(binObj->nrows * sizeof (int32_t),
            "extbin");

    int64_t index = 0;
    for (int row = 0; row < binObj->nrows; row++) {
        binObj->readBinIndex(row);
        int32_t extent = binObj->get_ext();
        extentbin[row] = extent;
        if (extent > 0) {
            baseRecord[row] = index;
            index += extent;
        } else {
            baseRecord[row] = -1;
        }
    }
    return 0;
}

/**
 * Given row, col and bin lookup the record number.  Note that row, col and bin
 * are NOT checked for out of bounds.  Mostly an internal function.
 * @param row row number of the bin (0-relative)
 * @param bin bin number (1-relative)
 * @return record number in the file, or -1 if bin is not in the file.
 */
int64_t L3File::rowbin2record(int32_t row, int64_t bin) {
    int64_t recordIndex;

    recordIndex = baseRecord[row];
    if (recordIndex == -1) // there are no records in this row
        return -1;
    L3Row* l3Row = getRow(row);
    if (l3Row == NULL || l3Row->getNumBins() == 0)
        return -1;
    L3Bin* l3Bin = l3Row->getBin(bin);
    if (l3Bin == NULL)
        return -1;
    return l3Bin->getRecordNum();
}

/**
 * Get the record number in the file for the given row and col
 * @param row row of the bin (0-relative)
 * @param col col of the bin (0-relative)
 * @return record number in the file or -1 if not found
 */
int64_t L3File::rowcol2record(int32_t row, int32_t col) {
    int64_t bin;
    bin = shape->rowcol2bin(row, col);
    return rowbin2record(row, bin);
}

/**
 * Get the record number in the file for the closest lat, lon
 * @param lat latitude of the bin
 * @param lon longitude of the bin
 * @return record number in the file or -1 if not found
 */
int64_t L3File::latlon2record(float lat, float lon) {
    int32_t row, col;
    int64_t bin;
    shape->latlon2rowcol(lat, lon, row, col);
    bin = shape->rowcol2bin(row, col);
    return rowbin2record(row, bin);
}

/**
 * Get the record number in the file for the given bin
 * @param bin bin number (1-relative)
 * @return record number in the file or -1 if not found
 */
int64_t L3File::bin2record(int64_t bin) {
    int32_t row;
    row = shape->bin2row(bin);
    return rowbin2record(row, bin);
}

int32_t L3File::getNumProducts() {
    if (binObj == NULL)
        return 0;
    return binObj->nprod();
}

string L3File::getProductName(size_t index) {
    return binObj->getProdName(index);
}

bool L3File::setActiveProductList(const char* prodStr) {
    //free the sum buffer since it depends on the number of products
    if (sumBuffer) {
        free(sumBuffer);
        sumBuffer = NULL;
    }
    clearCache();
    int status = binObj->read((char*) prodStr);
    if (status == -1)
        return false;

    // find number of products
    string tmpStr(prodStr);
    activeProdNameList.clear();
    boost::split(activeProdNameList, tmpStr, boost::is_any_of(","));

    // populate the product map
    prodMap = (size_t*) allocateMemory(
            activeProdNameList.size() * sizeof (size_t), "prodMap");
    numProds = activeProdNameList.size();
    for (size_t myIndex = 0; myIndex < numProds; myIndex++) {
        boost::trim(activeProdNameList[myIndex]);
        for (int32_t objIndex = 0; objIndex < binObj->n_active_prod;
                objIndex++) {
            if (activeProdNameList[myIndex].compare(
                    binObj->getActiveProdName(objIndex)) == 0) {
                prodMap[myIndex] = objIndex;
                break;
            }
        }
    }

    outBin.setNumProducts(activeProdNameList.size());

    return true;
}

int32_t L3File::getNumActiveProducts() {
    return activeProdNameList.size();
}

string L3File::getActiveProductName(size_t index) {
    return activeProdNameList[index];
}

L3Row* L3File::readRow(int32_t row) {
    L3Row* l3Row;

    if (row < 0 || row >= binObj->nrows)
        return NULL;

    int32_t extent = extentbin[row];
    int64_t base = baseRecord[row];
    if (rowList.size() < (uint32_t) numCacheRows) {
        l3Row = new L3Row(row, extent, numProds);
    } else {
        l3Row = rowList.back();
        rowList.pop_back();
        l3Row->setNumBins(extent);
        l3Row->row = row;
    }
    rowList.push_front(l3Row);

    if (extent == 0)
        return l3Row;

    // read the bin list
    binObj->readBinList(extent, base);

    // allocate the sum buffer if necessary
    if (sumBuffer == NULL)
        sumBuffer = (float*) allocateMemory(
            4 * binObj->nrows * binObj->n_active_prod * sizeof (float),
            "L3File::sumBuffer");

    binObj->setDataPtrAbsolute(base);
    binObj->readSums(sumBuffer, extent, -1);

    if (qualityBuffer)
        binObj->readQual((uint8_t*) qualityBuffer, extent, base);

    L3Bin* l3Bin;
    size_t sumOffset;
    for (int i = 0; i < extent; i++) {
        l3Bin = l3Row->binArray[i];
        l3Bin->binNum = binObj->get_bin_num(i);
        l3Bin->recordNum = base + i;
        l3Bin->nobs = binObj->get_nobs(i);
        l3Bin->nscenes = binObj->get_nscenes(i);
        l3Bin->timeRec = binObj->get_time_rec(i);
        l3Bin->weights = binObj->get_weights(i);
        for (size_t prod = 0; prod < numProds; prod++) {
            sumOffset = prodMap[prod] * extent * 2 + i * 2;
            l3Bin->sums[prod] = sumBuffer[sumOffset];
            l3Bin->sumSquares[prod] = sumBuffer[sumOffset + 1];
        }
        if (qualityBuffer)
            l3Bin->quality = qualityBuffer[i];
        else
            l3Bin->quality = L3Bin::qualityUnused;
    }
    return l3Row;
}

L3Row* L3File::getRow(int32_t row) {
    L3Row* l3Row;

    // check first cache position
    list<L3Row*>::iterator it = rowList.begin();
    if (it != rowList.end()) {
        l3Row = *it;
        if (l3Row->getRow() == row)
            return l3Row;
    }

    // check the rest of the cache
    for (it++; it != rowList.end(); it++) {
        l3Row = *it;
        if (l3Row->getRow() == row) {
            if (it != rowList.begin()) {
                rowList.erase(it);
                rowList.push_front(l3Row);
            }
            return l3Row;
        }
    }

    // not found in cache, so read the row
    return readRow(row);
}

/**
 * get the number of rows in the bin file
 * @return number of rows or -1 if file not opened yet
 */
int32_t L3File::getNumRows() {
    if (shape) {
        return shape->getNumRows();
    }
    return -1;
}

/**
 * return the bin closest to lat, lon
 * @param lat latitude of desired bin
 * @param lon longitude of desired bin
 * @return pointer to the bin object or NULL if bin not in file
 */
L3Bin* L3File::getClosestBin(float lat, float lon) {
    if (shape == NULL) {
        fprintf(stderr,
                "-E- %s %d: need to open L3File before looking for a bin.\n",
                __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    int32_t row, col;
    shape->latlon2rowcol(lat, lon, row, col);
    L3Row* l3Row = getRow(row);
    if (l3Row == NULL)
        return NULL;
    int64_t binNum = shape->rowcol2bin(row, col);
    return l3Row->getBin(binNum);
}

/**
 * return the bin at row/col
 * @param row of desired bin
 * @param col of desired bin
 * @return pointer to the bin object or NULL if bin not in file
 */
L3Bin* L3File::getBin(int32_t row, int32_t col) {
    if (shape == NULL) {
        fprintf(stderr,
                "-E- %s %d: need to open L3File before looking for a bin.\n",
                __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    if (row < 0 || row >= shape->getNumRows())
        return NULL;
    if (col < 0 || col >= shape->getNumCols(row))
        return NULL;

    L3Row* l3Row = getRow(row);
    if (l3Row == NULL)
        return NULL;
    int64_t binNum = shape->rowcol2bin(row, col);
    return l3Row->getBin(binNum);
}

/**
 * get the underlying hdf_bin object
 * @return the hdf_bin obj or NULL if file is not opened yet.
 */
Hdf::hdf_bin* L3File::getHdfBinObject() const {
    return binObj;
}

bool L3File::hasQuality() {
    if (binObj) {
        return binObj->has_qual();
    }
    return false;
}

void L3File::setQualityProcessing(bool val) {
    if (qualityBuffer)
        free(qualityBuffer);
    qualityBuffer = NULL;
    if (val) {
        if (shape == NULL) {
            fprintf(stderr,
                    "-E- %s %d: need to open L3File before setting qualityProcessing.\n",
                    __FILE__, __LINE__);
            exit(EXIT_FAILURE);
        }
        qualityBuffer = (uint8_t*) allocateMemory(2 * shape->getNumRows(),
                "L3File::qualityBuffer");
    }
}

bool L3File::getQualityProcessing() const {
    if (qualityBuffer)
        return true;
    else
        return false;
}
/**
 * @brief return the base bin number for a row
 *
 * @param row
 * @return int64_t bin_number
 */
int64_t L3File::getBaseBin(int32_t row) const {
    return baseRecord[row];
}

/**
 * @brief Get the extent of bins in a row
 *
 * @param row
 * @return int32_t extent
 */
int32_t L3File::getRowExtent(int32_t row) const {
    return extentbin[row];
}
} // namespace l3
