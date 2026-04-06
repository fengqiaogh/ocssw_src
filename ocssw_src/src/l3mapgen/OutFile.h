#ifndef OutFile_h
#define OutFile_h

#include <stdio.h>
#include <stdint.h>

#include <png.h>
#include <xtiffio.h>
#include <geotiffio.h>
#include <netcdf>
#include <unordered_map>
#include <meta_l3b.h>
#include <productInfo.h>
#include <vector>
#include <string>
#include <clo.h>
#include <genutils.h>

// earth radius in meters from WGS84 equatorial radius
#define EARTH_RADIUS 6378137.0
// #define EARTH_CIRCUMFERENCE 40075016.6856

// Metadata standard strings
#define INSTITUTION \
    "NASA Goddard Space Flight Center, Ocean Ecology Laboratory, Ocean Biology Processing Group"
#define LICENSE "https://science.nasa.gov/earth-science/earth-science-data/data-information-policy/"
#define NAMING_AUTHORITY "gov.nasa.gsfc.sci.oceandata"
#define KEYWORDS_VOCABULARY "NASA Global Change Master Directory (GCMD) Science Keywords"
#define KEYWORDS_OC "Oceans > Ocean Chemistry > Chlorophyll; Oceans > Ocean Optics > Ocean Color"
#define KEYWORDS_IOP "Oceans > Ocean Optics > Ocean Color"
#define KEYWORDS_SST "Oceans > Ocean Temperature > Sea Surface Temperature"
#define STDNAME_VOCABULARY "CF Standard Name Table v36"
#define CREATOR_NAME "NASA/GSFC/OBPG"
#define CREATOR_EMAIL "data@oceancolor.gsfc.nasa.gov"
#define CREATOR_URL "https://oceandata.sci.gsfc.nasa.gov"
#define PROJECT "Ocean Biology Processing Group (NASA/GSFC/OBPG)"
#define PUBLISHER_NAME "NASA/GSFC/OBPG"
#define PUBLISHER_EMAIL "data@oceancolor.gsfc.nasa.gov"
#define PUBLISHER_URL "https://oceandata.sci.gsfc.nasa.gov"
#define DOIAUTHORITY "https://dx.doi.org"

extern clo_optionList_t* optionList;

//---------------------------------------------------

class OutFile {
   public:
    enum ColorType { GRAYSCALE, COLOR_INDEX, RGB };

    enum ScaleType { LINEAR, LOG, ARCTAN };

    enum DataStorage { BYTE_DS, UBYTE_DS, SHORT_DS, USHORT_DS, INT_DS, UINT_DS, FLOAT_DS, DOUBLE_DS };

    enum PixValues { LAND_PIX = 254, FILL_PIX = 255 };  //< special output pixel values

    enum LatLonType { LAT_LON_OFF, LAT_LON_1D, LAT_LON_2D };  //< special output pixel values

    class ProductStuff {
       public:
        int32_t width;
        productInfo_t* productInfo;
        DataStorage dataStorage;  //< data type for the file storage
        ScaleType scaleType;      //< display LINEAR, LOG, ATan scaling
        double scale;             //< display slope for scaling.  If log it is actually log10(slope)
        double offset;            //< display offset for scaling.
        double minOutputVal;      //< display min output value
        double maxOutputVal;      //< display max output value
        double minVal;            //< display min physical value
        double maxVal;            //< display max physical value
        double missingValue;      //< missing value from product XML (val stored in file)
        double* lineData;
        double landPixelValue;

        /**
         * @brief Constructs a ProductStuff object with the specified width and product information; derived from OutFile class.
         * @param width The width of the product.
         * @param productInfo Pointer to the product information.
         */
        ProductStuff(int32_t width, const productInfo_t* productInfo, double farts);
    
        /**
         * @brief Copy constructor for ProductStuff.
         * @param productStuff The ProductStuff object to copy.
         */
        ProductStuff(const OutFile::ProductStuff& productStuff);
    
        /**
         * @brief Destructor for ProductStuff.
         */
        ~ProductStuff();

        int32_t getWidth() {
            return width;
        }

        /**
         * @brief Sets the scale factors (note that default minOutputVal=0, maxOutputVal=255)
         * @param min min geophysical value
         * @param max max geophysical value
         * @param log do you want log10 scaling
         */
        void setScale(double min, double max, ScaleType scaleType);

        /**
         * @brief Sets the scale factors with custom output values.
         * @param min Min geophysical value.
         * @param max Max geophysical value.
         * @param scaleType Type of scaling to use.
         * @param minOutput Minimum output value.
         * @param maxOutput Maximum output value.
         */
        void setScale(double min, double max, ScaleType scaleType, double minOutput, double maxOutput);

        /**
         * @brief Sets the scale factors (note that default minOutputVal=0, maxOutputVal=255).
         * @param scale Slope.
         * @param offset Intercept.
         * @param scaleType Type of scaling to calculate.
         */
        void setScaleOffset(double scale, double offset, ScaleType scaleType);

        /**
         * @brief Sets the scale offset values with custom output values.
         * @param scale Slope.
         * @param offset Intercept.
         * @param scaleType Type of scaling to use.
         * @param minOutput Minimum output value.
         * @param maxOutput Maximum output value.
         */
        void setScaleOffset(double scale, double offset, ScaleType scaleType, double minOutput, double maxOutput);

        /**
         * @brief Calculates the output value from a given input value. // TODO: Clarify what kind of output value
         * @param val The input value.
         * @return The calculated output value.
         */
        double calcOutputVal(double val) const;

        /**
         * @brief Calculates the physical value from a given output value.
         * @param val The output value.
         * @return The calculated physical value.
         */
        double calcPhysicalVal(double val) const;

        /**
         * @brief Calculates the output line values and stores them in the provided buffer.
         * @param lineBuffer The buffer to store the calculated output line values.
         */
        void calcOutputLineVals(void* lineBuffer) const;
    };

   protected:
    static constexpr uint8_t qualityUnused = 255;
    static constexpr double badPixelValue = -32767.0;
    static constexpr double landPixelValue = -32766.0;
    std::unordered_map<size_t, size_t> index2d3d;
    std::unordered_map<size_t, size_t> slice2dInWv3d;
    std::unordered_map<std::string, size_t> product3dAlreadySet, product2dIndexesLastIndex;
    std::string fileName;
    int32_t width;
    int32_t height;
    int32_t samplesperpixel;
    std::string qualityName;
    uint8_t* qualityData;
    uint32_t currentLine;  //< current line number (0 based)

    ColorType colorType;

    double fileMinVal;  //< current min value of data written to file
    double fileMaxVal;  //< current max value of data written to file
    double resolution;  //< geospatial resolution of a pixel in the center of the scene (meters)
    int deflate;        //< compression setting for netCDF files

    LatLonType fullLatLon;  //< should we write full-resolution lat/lon arrays
    double* latData;        //< pointer to row of latitude data
    double* lonData;        //< pointer to row of longitude data

    uint8_t* red;
    uint8_t* green;
    uint8_t* blue;
    uint8_t* rgbLand;
    bool transparent;

    meta_l3bType* metaData;
    std::string mapProjection;

    std::string proj4String;
    double tiepoints[6];
    double pixscale[3];

    std::vector<ProductStuff*> productStuff;

    /**
     * @brief Default constructor for OutFile.
     */
    OutFile();

    /**
     * @brief Adds a product for non-display type output files.
     * @param productInfo Pointer to the product information.
     * @return The index for the new product.
     */
    virtual int addProductNonDisplay(productInfo_t* productInfo);

   public:
    /**
     * @brief Destructor for OutFile.
     */
    virtual ~OutFile();

    /**
     * @brief Sets the size of the output file.
     * @param width The width of the file.
     * @param height The height of the file.
     */
    virtual void setSize(int32_t width, int32_t height);

    /**
     * @brief Gets the width of the output file.
     * @return The width of the output file.
     */
    virtual int32_t getWidth() const;

    /**
     * @brief Get the height of the output file.
     * @return Height of the output file.
     */
    virtual int32_t getHeight() const;

    /**
     * @brief Sets the file name of the output file.
     * @param fileName The file name to set.
     */
    virtual void setFileName(std::string fileName);


    virtual std::string getFileName() {
        return fileName;
    }

    /**
     * Open the output file. Make sure you have set these functions first:
     *      setSize, setScale, setPalette, setMetaData, setProductInfo, setQualityProcessing.
     * @param fileName
     * @return true if successful.
     */
    virtual bool open() = 0;
    virtual bool close() = 0;

    virtual double getMinValue(int32_t prod = 0) {
        return productStuff[prod]->minVal;
    }
    virtual double getMaxValue(int32_t prod = 0) {
        return productStuff[prod]->maxVal;
    }

    /**
     * @brief Retrieves the scale type as a string for the specified product.
     * @param prod The product index.
     * @return The scale type as a string.
     */
    virtual std::string getScaleTypeString(int32_t prod = 0);

    /**
     * @brief Sets the pixel value for the specified product.
     * @param x The x-coordinate of the pixel.
     * @param val The value to set.
     * @param prod The product index.
     */
    virtual void setPixel(int32_t x, double val, int32_t prod = 0);

    /**
     * @brief Sets the RGB pixel value. 
     * @param x x-coordinate of the pixel.
     * @param red red component value.
     * @param green green component value.
     * @param blue blue component value.
     */
    virtual void setPixelRGB(int32_t x, float red, float green, float blue);

    /**
     * @brief Enables transparency for the output file.
     */
    virtual void setTransparency();

    /**
     * @brief Sets the quality value for the specified pixel.
     * @param x X-coordinate of the pixel.
     * @param val Quality value to set.
     */
    virtual void setQuality(int32_t x, uint8_t val);

    /**
     * @brief Marks the specified pixel as a land pixel.
     * @param x x-coordinate of the pixel.
     */
    virtual void landPixel(int32_t x);

    /**
     * @brief Marks the specified pixel as a fill pixel.
     * @param x X-coordinate of the pixel.
     */
    virtual void fillPixel(int32_t x);

    /**
     * @brief Marks the specified pixel as a missing pixel.
     * @param x X-coordinate of the pixel.
     */
    virtual void missingPixel(int32_t x);

    /**
     * @brief Sets the latitude and longitude data.
     * @param lat Pointer to the latitude data.
     * @param lon Pointer to the longitude data.
     */
    virtual void setLatLon(double* lat, double* lon);


    virtual void writeLine() = 0;

    /**
     * @brief Sets the color palette for the output file.
     * @param paletteName The name of the palette to set.
     * @param applyMask Whether to apply a mask to the palette.
     * @return True if successful, else false.
     */
    virtual bool setPalette(const char* paletteName, bool applyMask);

    /**
     * @brief Sets the RGB values for land pixels.
     * @param rgb_land_string The RGB values for land pixels as a comma-separated string.
     */
    virtual void setLandRGB(const char* rgb_land_string);

    /**
     * @brief Sets the metadata for the output file.
     * @param metaData Pointer to the metadata to set.
     */
    virtual void setMetaData(meta_l3bType* metaData);

    virtual meta_l3bType* getMetadata() {
        return metaData;
    }

    /**
     * Add a product for display type output files
     * @param productInfo info structure to copy
     * @return the index for the new product
     */
    virtual int32_t addProduct(productInfo_t* productInfo, bool applyMask);


    virtual int32_t getNumProducts() {
        return productStuff.size();
    }

    /**
     * @brief Sets the map projection for the output file.
     * @param projection The map projection to set (as a string).
     */
    virtual void setMapProjection(std::string projection);

    /**
     * @brief Sets the number of filled pixels in the output file.
     * @param num The number of filled pixels.
     */
    virtual void setNumFilledPixels(int32_t num);

    /**
     * @brief Gets the number of filled pixels in the output file.
     * @return The number of filled pixels.
     */
    virtual int32_t getNumFilledPixels();

    /**
     * @brief Gets the percentage of filled pixels in the output file.
     * @return The percentage of filled pixels.
     */
    virtual float getPercentFilledPixels();

    /**
     * @brief Resets the minimum and maximum values of the output file.
     */
    virtual void resetFileMinMax();

    virtual double getFileMinVal() {
        return fileMinVal;
    }

    virtual double getFileMaxVal() {
        return fileMaxVal;
    }

    /**
     * @brief Sets the resolution for the output file.
     * @param resolutionStr The resolution as a string.
     */
    virtual void setResolution(std::string resolutionStr);

    virtual void setResolution(double resolution) {
        this->resolution = resolution;
    }

    virtual double getResolution() {
        return resolution;
    }

    virtual void setQualityName(std::string qualName) {
        qualityName = qualName;
    }

    std::string getQualityName() {
        return qualityName;
    }

    /**
     * @brief Enables or disables quality processing.
     * @param val True to enable quality processing, false to disable.
     */
    virtual void setQualityProcessing(bool val);

    /**
     * @brief Get the current status of quality processing.
     * @return True if quality processing is enabled, else false.
     */
    virtual bool getQualityProcessing();

    virtual void setDeflate(int val) {
        deflate = val;
    }
    virtual int getDeflate() {
        return deflate;
    }
    void setFullLatLon(bool val) {
        if (val)
            fullLatLon = LAT_LON_2D;
        else
            fullLatLon = LAT_LON_OFF;
    }

    /**
     * @brief Sets the PROJ.4 information for the output file.
     * @param projStr The PROJ.4 string.
     * @param minX The minimum X coordinate.
     * @param maxY The maximum Y coordinate.
     */
    void setProj4Info(std::string projStr, double minX, double maxY);
};

#endif // OUTFILE_H