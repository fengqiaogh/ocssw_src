#include <L3FileSMI.h>
#include <L3ShapeSMI.h>
#include <genutils.h>
#include <boost/algorithm/string.hpp>

#include <hdf.h>
#include <mfhdf.h>

using namespace std;
using namespace netCDF;
using namespace netCDF::exceptions;

namespace l3 {

//----------------------------------------------------------------------------
// L3FileSMI
//----------------------------------------------------------------------------

L3FileSMI::L3FileSMI() {
    ncFile = NULL;
    bzero(&metaData, sizeof (meta_l3bType));
    metaData.north = metaData.south = metaData.east = metaData.west = -999.0;
}

L3FileSMI::~L3FileSMI() {
    close();
}

bool L3FileSMI::readVarCF(NcVar &var, const std::vector<size_t> &start,
        const std::vector<size_t> &count, float* data) {
    if (var.isNull()) {
        return false;
    }
    size_t num = 1;
    for (size_t i = 0; i < count.size(); i++)
        num *= count[i];

    // get scale and offset if not float
    float scale = 1;
    float offset = 0;
    bool needScale = false;
    if ((var.getType().getId() != NC_FLOAT)
            && (var.getType().getId() != NC_DOUBLE)) {
        if (readAttribute(var, "scale_factor", &scale))
            needScale = true;
        if (readAttribute(var, "add_offset", &offset))
            needScale = true;
    }

    switch (var.getType().getId()) {
    case NC_FLOAT:
    {
        float badVal = missingPixelValue;
        readAttribute(var, "_FillValue", &badVal);
        var.getVar(start, count, data);
        for (size_t i = 0; i < num; i++)
            if (data[i] == badVal)
                data[i] = missingPixelValue;
        break;
    }
    case NC_DOUBLE:
    {
        double badVal = missingPixelValue;
        readAttribute(var, "_FillValue", &badVal);
        double* tmp = new double[num];
        var.getVar(start, count, tmp);
        for (size_t i = 0; i < num; i++)
            if (tmp[i] == badVal)
                data[i] = missingPixelValue;
            else
                data[i] = tmp[i];
        delete[] tmp;
        break;
    }
    case NC_BYTE:
    {
        int8_t badVal = -128;
        readAttribute(var, "_FillValue", &badVal);
        int8_t* tmp = new int8_t[num];
        var.getVar(start, count, tmp);
        if (needScale) {
            for (size_t i = 0; i < num; i++)
                if (tmp[i] == badVal)
                    data[i] = missingPixelValue;
                else
                    data[i] = tmp[i] * scale + offset;
        } else {
            for (size_t i = 0; i < num; i++)
                if (tmp[i] == badVal)
                    data[i] = missingPixelValue;
                else
                    data[i] = tmp[i];
        }
        delete[] tmp;
        break;
    }
    case NC_UBYTE:
    {
        uint8_t badVal = 255;
        readAttribute(var, "_FillValue", &badVal);
        uint8_t* tmp = new uint8_t[num];
        var.getVar(start, count, tmp);
        if (needScale) {
            for (size_t i = 0; i < num; i++)
                if (tmp[i] == badVal)
                    data[i] = missingPixelValue;
                else
                    data[i] = tmp[i] * scale + offset;
        } else {
            for (size_t i = 0; i < num; i++)
                if (tmp[i] == badVal)
                    data[i] = missingPixelValue;
                else
                    data[i] = tmp[i];
        }
        delete[] tmp;
        break;
    }
    case NC_SHORT:
    {
        short badVal = -32767;
        readAttribute(var, "_FillValue", &badVal);
        short* tmp = new short[num];
        var.getVar(start, count, tmp);
        if (needScale) {
            for (size_t i = 0; i < num; i++)
                if (tmp[i] == badVal)
                    data[i] = missingPixelValue;
                else
                    data[i] = tmp[i] * scale + offset;
        } else {
            for (size_t i = 0; i < num; i++)
                if (tmp[i] == badVal)
                    data[i] = missingPixelValue;
                else
                    data[i] = tmp[i];
        }
        delete[] tmp;
        break;
    }
    case NC_USHORT:
    {
        uint16_t badVal = 65535;
        readAttribute(var, "_FillValue", &badVal);
        uint16_t* tmp = new uint16_t[num];
        var.getVar(start, count, tmp);
        if (needScale) {
            for (size_t i = 0; i < num; i++)
                if (tmp[i] == badVal)
                    data[i] = missingPixelValue;
                else
                    data[i] = tmp[i] * scale + offset;
        } else {
            for (size_t i = 0; i < num; i++)
                if (tmp[i] == badVal)
                    data[i] = missingPixelValue;
                else
                    data[i] = tmp[i];
        }
        delete[] tmp;
        break;
    }
    case NC_INT:
    {
        int badVal = -32767;
        readAttribute(var, "_FillValue", &badVal);
        int* tmp = new int[num];
        var.getVar(start, count, tmp);
        if (needScale) {
            for (size_t i = 0; i < num; i++)
                if (tmp[i] == badVal)
                    data[i] = missingPixelValue;
                else
                    data[i] = tmp[i] * scale + offset;
        } else {
            for (size_t i = 0; i < num; i++)
                if (tmp[i] == badVal)
                    data[i] = missingPixelValue;
                else
                    data[i] = tmp[i];
        }
        delete[] tmp;
        break;
    }
    case NC_UINT:
    {
        uint32_t badVal = 4294967295;
        readAttribute(var, "_FillValue", &badVal);
        uint32_t* tmp = new uint32_t[num];
        var.getVar(start, count, tmp);
        if (needScale) {
            for (size_t i = 0; i < num; i++)
                if (tmp[i] == badVal)
                    data[i] = missingPixelValue;
                else
                    data[i] = tmp[i] * scale + offset;
        } else {
            for (size_t i = 0; i < num; i++)
                if (tmp[i] == badVal)
                    data[i] = missingPixelValue;
                else
                    data[i] = tmp[i];
        }
        delete[] tmp;
        break;
    }
    case NC_UINT64:
    {
        uint64_t badVal = 1.844674407 * pow(10, 19);
        readAttribute(var, "_FillValue", &badVal);
        uint64_t* tmp = new uint64_t[num];
        var.getVar(start, count, tmp);
        if (needScale) {
            for (size_t i = 0; i < num; i++)
                if (tmp[i] == badVal)
                    data[i] = missingPixelValue;
                else
                    data[i] = tmp[i] * scale + offset;
        } else {
            for (size_t i = 0; i < num; i++)
                if (tmp[i] == badVal)
                    data[i] = missingPixelValue;
                else
                    data[i] = tmp[i];
        }
        delete[] tmp;
        break;
    }
    default:
        if (want_verbose)
            printf("Error - L3FileSMI::readVarCF - data type %s not supported",
                var.getType().getName().c_str());
        return false;
    }
    return true;
}

bool L3FileSMI::open(const char* fileName) {
    close();

    // bail if this is an hdf file
    if (Hishdf(fileName))
        return false;

    try {
        ncFile = new NcFile(fileName, NcFile::read);
    } catch (NcException& e) {
        e.what();
        if (want_verbose)
            printf("Error - %s is not a valid netCDF file\n", fileName);
        ncFile = NULL;
        return false;
    }

    if (ncFile->isNull()) {
        if (want_verbose)
            printf("Error - %s is not a valid netCDF file\n", fileName);
        close();
        return false;
    }

    string tmpStr = "";
    readAttribute("title", tmpStr);
    if (tmpStr.find("Level-3 Standard Mapped Image") == string::npos) {
        if (want_verbose)
            printf("Error - %s is not a Level-3 SMI file\n", fileName);
        close();
        return false;
    }

    NcDim latDim = ncFile->getDim("lat");
    if (latDim.isNull()) {
        if (want_verbose)
            printf("Error - Could not find lat dimension\n");
        close();
        return false;
    }
    size_t numRows = latDim.getSize();

    NcDim lonDim = ncFile->getDim("lon");
    if (latDim.isNull()) {
        if (want_verbose)
            printf("Error - Could not find lon dimension\n");
        close();
        return false;
    }
    size_t numCols = lonDim.getSize();

    read_l3b_meta_netcdf4(ncFile->getId(), &metaData);

    shape = new L3ShapeSMI(numRows, numCols, metaData.north, metaData.south,
            metaData.east, metaData.west);

    initRecordLookup();

    // fill up the name list
    prodNameList.clear();
    multimap<string, NcVar> varMap = ncFile->getVars();
    multimap<string, NcVar>::iterator iter = varMap.begin();
    while (iter != varMap.end()) {
        string varName = iter->first;
        if (iter->second.getDimCount() == 2) {
            size_t varRows = iter->second.getDim(0).getSize();
            size_t varCols = iter->second.getDim(1).getSize();
            size_t fileRows = shape->getNumRows();
            size_t fileCols = shape->getNumCols(0);
            if (varRows == fileRows && varCols == fileCols)
                prodNameList.push_back(varName);
        }
        iter++;
    }

    return true;
}

void L3FileSMI::close() {
    if (ncFile) {
        delete ncFile;
        ncFile = NULL;
    }
    metaData.north = metaData.south = metaData.east = metaData.west = -999.0;
    prodNameList.clear();
    activeProdNameList.clear();
    prodVarList.clear();
    L3File::close();
}

meta_l3bType* L3FileSMI::getMetaData() {
    return &metaData;
}

/**
 * init the internal data structures to lookup the record index
 * @return status 0 if good, -1 if bad
 */
int L3FileSMI::initRecordLookup() {

    if (shape->getNumRows() == 0) {
        fprintf(stderr,
                "-E- %s %d: initRecordLookup: Number of rows needs to be set and bin grid must be initialized first.\n",
                __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }

    // remove memory if it exists
    if (baseRecord) {
        free(baseRecord);
    }
    baseRecord = (int64_t*) allocateMemory(
            shape->getNumRows() * sizeof (int64_t), "baseRecord");
    if (extentbin) {
        free(extentbin);
    }
    extentbin = (int32_t*) allocateMemory(
            shape->getNumRows() * sizeof (int32_t), "extbin");

    int32_t index = 0;
    for (int row = 0; row < shape->getNumRows(); row++) {
        extentbin[row] = shape->getNumCols(row);
        baseRecord[row] = index;
        index += shape->getNumCols(row);
    }
    return 0;
}

int32_t L3FileSMI::getNumProducts() {
    return prodNameList.size();
}

string L3FileSMI::getProductName(size_t index) {
    return prodNameList[index];
}

bool L3FileSMI::setActiveProductList(const char* prodStr) {
    //free the sum buffer since it depends on the number of products
    if (sumBuffer) {
        free(sumBuffer);
        sumBuffer = NULL;
    }
    clearCache();
    outBin.clear();

    activeProdNameList.clear();
    prodVarList.clear();

    if (ncFile == NULL || ncFile->isNull()) {
        return false;
    }

    boost::split(activeProdNameList, prodStr, boost::is_any_of(","));
    for (size_t i = 0; i < activeProdNameList.size(); i++) {
        boost::trim(activeProdNameList[i]);
        NcVar tmpVar = ncFile->getVar(activeProdNameList[i]);
        if (tmpVar.isNull()) {
            activeProdNameList.clear();
            prodVarList.clear();
            if (want_verbose)
                printf("Error - product %s not found in file\n",
                    activeProdNameList[i].c_str());
            return false;
        }
        prodVarList.push_back(tmpVar);
    }
    return true;
}

L3Row* L3FileSMI::readRow(int32_t row) {
    L3Row* l3Row;

    if (row < 0 || row >= shape->getNumRows())
        return NULL;

    if (rowList.size() < (uint32_t) numCacheRows) {
        l3Row = new L3Row(row, shape->getNumCols(row),
                activeProdNameList.size());
    } else {
        l3Row = rowList.back();
        rowList.pop_back();
        l3Row->row = row;
    }
    rowList.push_front(l3Row);

    // allocate the sum buffer if necessary
    if (sumBuffer == NULL)
        sumBuffer = (float*) allocateMemory(
            shape->getNumCols(row) * sizeof (float), "L3FileSMI::sumBuffer");

    vector<size_t> start, count;
    start.push_back(shape->getNumRows() - row - 1);
    start.push_back(0);
    count.push_back(1);
    count.push_back(shape->getNumCols(row));

    for (size_t prodIndex = 0; prodIndex < prodVarList.size(); prodIndex++) {
        readVarCF(prodVarList[prodIndex], start, count, sumBuffer);

        L3Bin* l3Bin;
        for (int32_t i = 0; i < shape->getNumCols(row); i++) {
            l3Bin = l3Row->binArray[i];
            l3Bin->binNum = shape->rowcol2bin(row, i);
            l3Bin->recordNum = l3Bin->binNum;
            l3Bin->quality = L3Bin::qualityUnused;

            if (sumBuffer[i] == missingPixelValue) {
                l3Bin->nobs = 0;
                l3Bin->nscenes = 0;
                l3Bin->weights = 0;
                l3Bin->sums[prodIndex] = missingPixelValue;
                l3Bin->sumSquares[prodIndex] = missingPixelValue;
            } else {
                l3Bin->nobs = 1;
                l3Bin->nscenes = 1;
                l3Bin->weights = 1;
                l3Bin->sums[prodIndex] = sumBuffer[i];
                l3Bin->sumSquares[prodIndex] = sumBuffer[i] * sumBuffer[i];
            }
        }
    }
    return l3Row;
}

bool L3FileSMI::hasQuality() {
    return false;
}

} // namespace l3
