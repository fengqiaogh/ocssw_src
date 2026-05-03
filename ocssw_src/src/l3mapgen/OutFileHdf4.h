#ifndef OUTFILEHDF4_H
#define OUTFILEHDF4_H

#include "OutFile.h"

class OutFileHdf4 : public OutFile {
    void* fileData;
    int32_t sdfid;          ///< HDF SD file ID
    int32_t sdsid;          ///< HDF SDS dataset ID
    int32_t quality_sdsid;  ///< HDF SDS dataset ID for SST quality factor
    int32_t hdfDataType;

   public:
    /**
     * @brief Constructor.
     */
    OutFileHdf4();

    /**
     * @brief Virtual destructor.
     */
    virtual ~OutFileHdf4();

    /**
     * @brief Set the size of the HDF4 output file.
     * @param width  Width of the output file.
     * @param height Height of the output file.
     */
    virtual void setSize(int32_t width, int32_t height);

    /**
     * @brief Open the HDF4 output file and set its metadata attributes.
     * 
     * This function initializes and opens an HDF4 output file, setting various 
     * metadata attributes (eg: product name, sensor information, mission details, 
     * and spatial/temporal information). It also creates the necessary datasets for 
     * storing the output data and any associated quality data.
     * 
     * @return true if the file is successfully opened, exit the program if an errors are encountered. 
     */
    virtual bool open();

    /**
     * @brief Closes the HDF4 file and writes remaining metadata. 
     * @return true if the file is successfully closed, else false.
     */
    virtual bool close();

    /**
     * @brief Writes a line of data to the HDF4 file.
     * 
     * This function calculates the output line values for the current line and then
     * writes the data to the HDF4 file. If quality data is available, it writes that too.
     */
    virtual void writeLine();

    /**
     * @brief Adds a product to the HDF4 file without displaying it.
     * @param productInfo Pointer to the product information.
     * @param productAttr The ProductAttribute object containing the mean method information.
     * @return int32_t The index of the added product. // TODO: Confirm if this is the index
     */
    virtual int addProduct(productInfo_t* productInfo, bool applyMask, const ProductL3Attributes & productAttr);
};

#endif // OUTFILEHDF4_H