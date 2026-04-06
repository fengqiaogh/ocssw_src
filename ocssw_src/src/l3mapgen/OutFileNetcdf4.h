#ifndef OUTFILENETCDF4_H
#define OUTFILENETCDF4_H

#include "OutFile.h"

class OutFileNetcdf4 : public OutFile {
    void* fileData;
    netCDF::NcFile* ncFile;               ///< NetCDF file handle
    std::vector<netCDF::NcVar> prodVars;  ///< NetCDF variables
    netCDF::NcVar qualVar;                ///< NetCDF variable for quality factor
    netCDF::NcVar latVar, lonVar;         ///< NetCDF variable for latitude/longitude

    /**
     * @brief Removes the wavelength information from a long name string.
     * @param longNameWithWv The input string containing wavelength information.
     * @return A string with the wavelength removed.
     */
    std::string removeWvFromLongName(const std::string& longNameWithWv);

   public:
    /**
     * @brief Constructor for the OutFileNetcdf4 class.
     * Initializes an instance of the OutFileNetcdf4 class, setting
     * the ncFile and fileData pointers to NULL.
     */
    OutFileNetcdf4();

    /**
     * @brief Virtual destructor.
     */
    virtual ~OutFileNetcdf4();

    /**
     * @brief Sets the size of the NetCDF file.
     * @param width width of the file.
     * @param height height of the file.
     */
    virtual void setSize(int32_t width, int32_t height);

    // modifications must be done here (TODO: What did the author of this code mean by this?)
    /**
     * @brief Opens the NetCDF file and sets its metadata attributes.
     * @return true if the file is successfully opened and attributes are set, false otherwise.
     */
    virtual bool open();

    /**
     * @brief Closes the NetCDF file and writes remaining metadata.
     * @return true if the file is successfully closed.
     */
    virtual bool close();

    //  modifications must be done here (TODO: Again, what modifications?)
    /**
     * @brief Writes a line of data to the NetCDF file.
     * 
     * Calculates the output line values for the current line, then writes the data to
     * the NetCDF file. If quality data is available, it writes that too, as well as
     * latitude and longitude data if required.
     */
    virtual void writeLine();

    /**
     * @brief Adds a product to the NetCDF file.
     * @param productInfo Pointer to the product information structure.
     * @return int32_t The index of the added product.
     */
    virtual int addProduct(productInfo_t* productInfo, bool applyMask);

   private:

    /**
     * @brief Initializes compression for a NetCDF variable.
     * 
     * Sets the chunking and compression options for a given NetCDF variable
     * based on the defined deflate level.
     * 
     * @param var The NetCDF variable to be compressed.
     */
    virtual void initCompression(netCDF::NcVar var);

    /**
     * @brief Creates a product variable in the NetCDF file.
     * 
     * Creates a variable for the given product information in the NetCDF file
     * and sets the standard metadata attributes.
     * 
     * @param pInfo Pointer to the product information structure.
     * @param ncType The NetCDF data type of the variable.
     * @param ncDim The dimensions of the variable.
     * @return NcVar The created NetCDF variable.
     */
    virtual netCDF::NcVar createProduct(productInfo_t* productInfo, const netCDF::NcType& ncType,
                                        const std::vector<netCDF::NcDim> ncDim);
                                    
    /**
     * @brief Gets the NetCDF data type for the given data storage type.
     * @param dataStorage The data storage type.
     * @return NcType The corresponding NetCDF data type.
     */
    virtual netCDF::NcType getDataType(DataStorage dataStorage);
};

#endif // OUTFILENETCDF4_H
