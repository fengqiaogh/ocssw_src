#ifndef L3FileSMI_h
#define L3FileSMI_h

#include <netcdf>
#include <L3File.h>
#include <vector>
#include <deque>

namespace l3 {

/**
 * class to read the information out of a L3 SMI file
 */
class L3FileSMI : public L3File {
protected:
    std::vector<std::string> prodNameList; ///< array of all product names in file
    netCDF::NcFile *ncFile;
    meta_l3bType metaData;
    std::vector<netCDF::NcVar> prodVarList;
    size_t wvIdx = 0;
    std::deque<size_t> wvIdxList;

    virtual int initRecordLookup();

    virtual L3Row* readRow(int32_t row);

public:
    L3FileSMI();
    virtual ~L3FileSMI();

    template <typename T>
    bool getProductAttribute(const std::string& prodName, const std::string& attributeName, T* attrValue) {
        return false;
    }

    template<typename TheType>
    bool readAttribute(std::string name, TheType &val) {
        if (ncFile == NULL || ncFile->isNull()) {
            return false;
        }
        std::multimap<std::string, netCDF::NcGroupAtt> attributeList = ncFile->getAtts();
        std::multimap<std::string, netCDF::NcGroupAtt>::iterator myIter;
        myIter = attributeList.find(name);
        if (myIter == attributeList.end())
            return false;
        myIter->second.getValues(val);
        return true;
    }

    template<typename TheType>
    bool readAttribute(const netCDF::NcGroup &group,
            std::string name, TheType &val) {
        if (group.isNull()) {
            return false;
        }
        std::multimap<std::string, netCDF::NcGroupAtt> attributeList = ncFile->getAtts();
        std::multimap<std::string, netCDF::NcGroupAtt>::iterator myIter;
        myIter = attributeList.find(name);
        if (myIter == attributeList.end())
            return false;
        myIter->second.getValues(val);
        return true;
    }

    template<typename TheType>
    bool readAttribute(const netCDF::NcVar &var,
            std::string name, TheType &val) {
        if (var.isNull()) {
            return false;
        }
        std::map<std::string, netCDF::NcVarAtt> attributeList = var.getAtts();
        std::map<std::string, netCDF::NcVarAtt>::iterator myIter;
        myIter = attributeList.find(name);
        if (myIter == attributeList.end())
            return false;
        myIter->second.getValues(val);
        return true;
    }

    template<typename TheType>
    bool readAttribute(const netCDF::NcVar &var,
            std::string name, TheType *val) {
        if (var.isNull()) {
            return false;
        }
        std::map<std::string, netCDF::NcVarAtt> attributeList = var.getAtts();
        std::map<std::string, netCDF::NcVarAtt>::iterator myIter;
        myIter = attributeList.find(name);
        if (myIter == attributeList.end())
            return false;
        myIter->second.getValues(val);
        return true;
    }

    virtual bool readVarCF(netCDF::NcVar &var, const std::vector<size_t> &start,
            const std::vector<size_t> &count, float* data);

    virtual bool open(const char* fileName);
    virtual void close();
    virtual meta_l3bType* getMetaData();
    virtual int32_t getNumProducts();
    virtual std::string getProductName(size_t index = 0);
    virtual bool setActiveProductList(const char* prodStr);
    virtual bool hasQuality();

    virtual bool getWavelength(std::string name, float& wavelength);
    virtual bool checkWavelength(const char* wvl);
    virtual bool getWavelengthList(std::vector<std::string>& wvlist);

};

} // namespace l3

#endif
