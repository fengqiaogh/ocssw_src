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

//earth radius in meters from WGS84 equatorial radius
#define EARTH_RADIUS 6378137.0
// #define EARTH_CIRCUMFERENCE 40075016.6856

// Metadata standard strings
#define INSTITUTION "NASA Goddard Space Flight Center, Ocean Ecology Laboratory, Ocean Biology Processing Group"
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

    enum ColorType {
        Grayscale, ColorIndex, RGB
    };

    enum ScaleType {
        Linear, Log, ArcTan
    };

    enum DataStorage {
        ByteDS, UByteDS, ShortDS, UShortDS, IntDS, UIntDS, FloatDS, DoubleDS
    };

    enum PixValues {
        landPix = 254, fillPix = 255
    }; //< special output pixel values

    enum LatLonType {
        LatLonOff, LatLon1D, LatLon2D
    }; //< special output pixel values

    class ProductStuff {
    public:

        int32_t width;
        productInfo_t* productInfo;
        DataStorage dataStorage; //< data type for the file storage
        ScaleType scaleType; //< display Linear, Log, ATan scaling
        double scale; //< display slope for scaling.  If log it is actually log10(slope)
        double offset; //< display offset for scaling.
        double minOutputVal; //< display min output value
        double maxOutputVal; //< display max output value
        double minVal; //< display min physical value
        double maxVal; //< display max physical value
        double missingValue; //< missing value from product XML (val stored in file)
        double* lineData;
        double landPixelValue;

        ProductStuff(int32_t width, const productInfo_t* productInfo, double landPixelValue);
        ProductStuff(const OutFile::ProductStuff& pStuff);
        ~ProductStuff();

        int32_t getWidth() {
            return width;
        }
        void setScale(double min, double max, ScaleType scaleType);
        void setScale(double min, double max, ScaleType scaleType,
                double minOutput, double maxOutput);
        void setScaleOffset(double scale, double offset, ScaleType scaleType);
        void setScaleOffset(double scale, double offset, ScaleType scaleType,
                double minOutput, double maxOutput);
        double calcOutputVal(double val) const;
        double calcPhysicalVal(double val) const;
        void calcOutputLineVals(void* lineBuffer) const;
    };

protected:

    static constexpr uint8_t qualityUnused = 255;
    static constexpr double badPixelValue = BAD_FLT;
    double landPixelValue;
    std::unordered_map<size_t, size_t> index_2d_3d;
    std::unordered_map<size_t, size_t> slice_2d_in_wv3d;
    std::unordered_map<std::string, size_t> product_3d_already_set, prod2d_indexes_last_index;
    std::string fileName;
    int32_t width;
    int32_t height;
    int32_t samplesperpixel;
    std::string qualityName;
    uint8_t* qualityData;
    uint32_t currentLine; //< current line number (0 based)

    ColorType colorType;

    double fileMinVal; //< current min value of data written to file
    double fileMaxVal; //< current max value of data written to file
    double resolution; //< geospatial resolution of a pixel in the center of the scene (meters)
    int deflate; //< compression setting for netCDF files

    LatLonType fullLatLon; //< should we write full-resolution lat/lon arrays
    double* latData; //< pointer to row of latitude data
    double* lonData; //< pointer to row of longitude data

    uint8_t* red;
    uint8_t* green;
    uint8_t* blue;
    uint8_t* rgb_land;
    bool transparent;

    meta_l3bType* metaData;
    std::string mapProjection;

    std::string proj4String;
    double tiepoints[6];
    double pixscale[3];

    std::vector<ProductStuff*> productStuff;

    OutFile();
    virtual int addProductNonDisplay(productInfo_t* productInfo);

public:
    virtual ~OutFile();

    virtual void setSize(int32_t width, int32_t height);
    virtual int32_t getWidth() const;
    virtual int32_t getHeight() const;
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
    virtual std::string getScaleTypeString(int32_t prod = 0);

    virtual void setPixel(int32_t x, double val, int32_t prod = 0);
    virtual void setPixelRGB(int32_t x, float red, float green, float blue);
    virtual void setTransparency();
    virtual void setQuality(int32_t x, uint8_t val);
    virtual void landPixel(int32_t x);
    virtual void fillPixel(int32_t x);
    virtual void missingPixel(int32_t x);
    virtual void setLatLon(double* lat, double* lon);
    virtual void writeLine() = 0;
    virtual bool setPalette(const char* paletteName, bool applyMask);
    virtual void setLandRGB(const char* rgb_land_string);
    virtual void setMetaData(meta_l3bType* metaData);

    virtual meta_l3bType* getMetadata() {
        return metaData;
    }
    virtual int32_t addProduct(productInfo_t* productInfo);

    virtual int32_t getNumProducts() {
        return productStuff.size();
    }
    virtual void setMapProjection(std::string projection);
    virtual void setNumFilledPixels(int32_t num);
    virtual int32_t getNumFilledPixels();
    virtual float getPercentFilledPixels();
    virtual void resetFileMinMax();

    virtual double getFileMinVal() {
        return fileMinVal;
    }

    virtual double getFileMaxVal() {
        return fileMaxVal;
    }
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
    virtual void setQualityProcessing(bool val);
    virtual bool getQualityProcessing();

    virtual void setDeflate(int val) {
        deflate = val;
    }
    virtual int getDeflate() {
        return deflate;
    }
    void setFullLatLon(bool val) {
        if(val)
            fullLatLon = LatLon2D;
        else
            fullLatLon = LatLonOff;
    }

    void setProj4Info(std::string projStr, double minX, double maxY);

};

//---------------------------------------------------

class OutFile_pgm : public OutFile {
protected:
    FILE *outfp;
    uint8_t* fileData;

public:
    OutFile_pgm();
    virtual ~OutFile_pgm();
    virtual void setSize(int32_t width, int32_t height);
    virtual bool open();
    virtual bool close();
    virtual void writeLine();
};

//---------------------------------------------------

class OutFile_ppm : public OutFile_pgm {
public:
    virtual void setSize(int32_t width, int32_t height);
    virtual bool open();
    virtual void writeLine();
};

//---------------------------------------------------

class OutFile_ppm_rgb : public OutFile {
    FILE *outfp;
    uint8_t* fileData;

public:
    OutFile_ppm_rgb();
    virtual ~OutFile_ppm_rgb();
    virtual void setSize(int32_t width, int32_t height);
    virtual bool open();
    virtual bool close();
    virtual void setPixel(int32_t x, double val, int32_t prod = 0);
    virtual void setPixelRGB(int32_t x, float red, float green, float blue);
    virtual void landPixel(int32_t x);
    virtual void fillPixel(int32_t x);
    virtual void missingPixel(int32_t x);
    virtual void writeLine();
};

//---------------------------------------------------

class OutFile_png : public OutFile {
    FILE *outfp;
    uint8_t* fileData;
    bool isColor;
    png_structp png_ptr;
    png_infop info_ptr;
    uint num_text;

public:
    OutFile_png(bool color);
    virtual ~OutFile_png();
    virtual void setSize(int32_t width, int32_t height);
    virtual bool open();
    virtual bool close();
    virtual void writeLine();
};

//---------------------------------------------------

class OutFile_png_rgb : public OutFile {
    FILE *outfp;
    uint8_t* fileData;
    png_structp png_ptr;
    png_infop info_ptr;
    uint num_text;

public:
    OutFile_png_rgb();
    virtual ~OutFile_png_rgb();
    virtual void setSize(int32_t width, int32_t height);
    virtual bool open();
    virtual bool close();
    virtual void setPixel(int32_t x, double val, int32_t prod = 0);
    virtual void setPixelRGB(int32_t x, float red, float green, float blue);
    virtual void landPixel(int32_t x);
    virtual void fillPixel(int32_t x);
    virtual void missingPixel(int32_t x);
    virtual void writeLine();
};

//---------------------------------------------------

class OutFile_tiff: public OutFile {
protected:
    TIFF *tiff = NULL;
    GTIF *gtif = NULL;
public:
    virtual ~OutFile_tiff();
    virtual bool open();
    virtual bool close();
    virtual void setTiffColor() = 0;
};

class OutFile_tiff_color: public OutFile_tiff {
    uint8_t* fileData = NULL;
public:
    virtual ~OutFile_tiff_color();
    virtual void setSize(int32_t width, int32_t height);
    virtual void writeLine();
    virtual void setTiffColor();
};

class OutFile_tiff_gray: public OutFile_tiff {
    float* fileData = NULL;
public:
    OutFile_tiff_gray();
    virtual ~OutFile_tiff_gray();
    virtual void setSize(int32_t width, int32_t height);
    virtual void writeLine();
    virtual void setTiffColor();
};

class OutFile_tiff_rgb: public OutFile_tiff {
    uint8_t* fileData = NULL;
public:
    virtual ~OutFile_tiff_rgb();
    virtual void setSize(int32_t width, int32_t height);
    virtual void writeLine();
    virtual void setTiffColor();
    virtual void setPixel(int32_t x, double val, int32_t prod = 0);
    virtual void setPixelRGB(int32_t x, float red, float green, float blue);
    virtual void landPixel(int32_t x);
    virtual void fillPixel(int32_t x);
    virtual void missingPixel(int32_t x);
};

//---------------------------------------------------

class OutFile_hdf4 : public OutFile {
    void* fileData;
    int32_t sdfid; ///< HDF SD file ID
    int32_t sdsid; ///< HDF SDS dataset ID
    int32_t quality_sdsid; ///< HDF SDS dataset ID for SST quality factor
    int32_t hdfDataType;

public:
    OutFile_hdf4();
    virtual ~OutFile_hdf4();
    virtual void setSize(int32_t width, int32_t height);
    virtual bool open();
    virtual bool close();
    virtual void writeLine();
    virtual int addProduct(productInfo_t* productInfo);
};

//---------------------------------------------------

class OutFile_netcdf4 : public OutFile {
    void* fileData;
    netCDF::NcFile* ncFile;  ///< NetCDF file handle
    std::vector<netCDF::NcVar> prodVars; ///< NetCDF variables
    netCDF::NcVar qualVar; ///< NetCDF variable for quality factor
    netCDF::NcVar latVar, lonVar; ///< NetCDF variable for latitude/longitude
   

public:
    OutFile_netcdf4();
    virtual ~OutFile_netcdf4();
    virtual void setSize(int32_t width, int32_t height);
    virtual bool open();
    virtual bool close();
    virtual void writeLine();
    virtual int addProduct(productInfo_t* productInfo);

private:
    virtual void initCompression(netCDF::NcVar var);
    virtual netCDF::NcVar createProduct(productInfo_t* productInfo,
                                        const netCDF::NcType& ncType,
                                        const std::vector<netCDF::NcDim> ncDim);
    virtual netCDF::NcType getDataType(DataStorage dataStorage);
};

#endif
