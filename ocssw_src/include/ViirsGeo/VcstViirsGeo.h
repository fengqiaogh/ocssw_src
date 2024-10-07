/**************************************************************************
 *
 * NAME: VcstViirsGeo.h
 *
 * DESCRIPTION: VcstViirsGeo encapsulates the geolocation algorithm. It
 * contains function declarations for functions used in the geolocation
 * routines.
 *
 *
 **************************************************************************/

#ifndef VcstViirsGeo_h
#define VcstViirsGeo_h

#include <string>
#include <VcstGeoDataStructs.h>
#include <VcstGeoParameters.h>
#include <VcstCmnGeo.h>
#include <VcstParamsReader.h>

#include <VcstGeoRctnglStruct.h>

const unsigned char GEO_SCAN_QUALITY_HAMIMP_GOODDATA = 0x00; //00000000
const unsigned char GEO_SCAN_QUALITY_HAMIMP_BADDATA = 0x04; //00000x00
const unsigned char GEO_SCAN_QUALITY_HAMIMP_DEGRADDATA = 0x08; //0000x000
const unsigned char GEO_SCAN_QUALITY_HAMIMP_MISSINGDATA = 0x0c; //0000xx00
const unsigned char GEO_SCAN_QUALITY_NOT_ABOVE_SAA = 0x00; //00000000
const unsigned char GEO_SCAN_QUALITY_ABOVE_SAA = 0x10; //000x0000
const unsigned char GEO_SCAN_QUALITY_NO_SOLAR_ECLIPSE = 0x00; //00000000
const unsigned char GEO_SCAN_QUALITY_SOLAR_ECLIPSE = 0x20; //00x00000
//for DNB only
const unsigned char GEO_SCAN_QUALITY_NO_LUNAR_ECLIPSE = 0x00; //00000000
const unsigned char GEO_SCAN_QUALITY_LUNAR_ECLIPSE = 0x40; //0x000000
const unsigned char GEO_SCAN_QUALITY_HA_MIRROR_SIDE_SHIFT = 7; //x0000000
const unsigned char GEO_SCAN_QUALITY_HA_MIRROR_SIDE_MASK = 0x80; //x0000000

//DR4767, DR4759 new flags for SCE side, non nominal ham start
const unsigned char GEO_SCAN_QUALITY_SCESIDE_A_ON = 0x00; //0000000x
const unsigned char GEO_SCAN_QUALITY_SCESIDE_B_ON = 0x01; //000000x0
const unsigned char GEO_SCAN_QUALITY_SCESIDE_INVALID = 0x02; //000000xx
const unsigned char GEO_SCAN_QUALITY_HAM_START_IS_NOMINAL = 0x00; //00000000
const unsigned char GEO_SCAN_QUALITY_HAM_START_NOT_NOMINAL = 0x04; //00000x00

const unsigned char GEO_PIXEL_QUALITY_INPUT_VALID = 0x00; //00000000
const unsigned char GEO_PIXEL_QUALITY_INPUT_INVALID = 0x01; //0000000x
const unsigned char GEO_PIXEL_QUALITY_POINTING_GOOD = 0x00; //00000000
const unsigned char GEO_PIXEL_QUALITY_POINTING_BAD = 0x02; //000000x0
const unsigned char GEO_PIXEL_QUALITY_TERRAIN_GOOD = 0x00; //00000000
const unsigned char GEO_PIXEL_QUALITY_TERRAIN_BAD = 0x04; //00000x00
const unsigned char GEO_PIXEL_QUALITY_SOLARANGLE_VALID = 0x00; //00000000
const unsigned char GEO_PIXEL_QUALITY_SOLARANGLE_INVALID = 0x08; //0000x000

typedef enum {
    OFF, ON
} switch_t;

// macros needed to pass heap arrays to functions requiring double pointers
// TBD: combine into one overloaded function
#define ptrCastFloat(array,ny,ptrarr) {  \
        ptrarr = new float *[ny];  \
        for (size_t y = 0; y < ny; y++)  \
            ptrarr[y] = array[y];  \
    }
#define ptrCastUchar(array,ny,ptrarr) {  \
        ptrarr = new unsigned char *[ny]; \
        for (size_t y = 0; y < ny; y++)  \
            ptrarr[y] = array[y];  \
    }

/***************************************************************************
NOTE: geolocation data for every pixel.

NOTE: Moon phase is in degrees and done according to standard Astronomy
usage.  That is, Moon Phase ranges from 0 to +180. Phase=0.0 means that the
center of the Sun, the center of the Earth, and the center of the Moon
are all in a perfect straight line, with the Moon and the Sun on
opposite sides of the Earth (a lunar eclipse).  A "Full Moon" is the
point in each Moon cycle where the phase reaches minimum.  Phase=+180
would mean that the center of the Moon and the center of the Sun are in
the same apparent place as viewed from the center of the Earth (a Solar
eclipse on the Earth's Equator).  A "New Moon" is the point in each
Moon cycle where the phase reaches a maximum.  We also record the
fraction of the Moon illuminated (mifrac).  mifrac=0.0 is a dark Moon.
mifrac=1.00 means that all of the visible Moon disk is illuminated.
Since phase and mifrac have very small variation during a granule, we
provide only one value of phase and mifrac for the granule.

These structures are the baseline for input to EDRs and data delivery.
These will change as RDR to SDR algorithms are delivered from the
sensor vendors.
 **********************/

/* Quality Flags */
const unsigned char VIIRS_GEO_INVALID_SCAN_START_TIME = 0x01; /* 00000001 < lsb */
const unsigned char VIIRS_GEO_CMN_GEO_FAIL = 0x02; /* 00000010 < lsb */
const unsigned char VIIRS_GEO_INVALID_MIRROR_SIDE = 0x04; /* 00000100 < lsb */
const unsigned char VIIRS_GEO_TERRAIN_CORR_FAIL = 0x08; /* 00001000 < lsb */
const unsigned char VIIRS_GEO_INVALID_TEL_HAM_ENC_DATA = 0x10; /* 00010000 < lsb */

typedef enum {
    MISS, OUTOFBND, AUTO, MAX_NUM_GEO_QUAL
} GeoQual;

#define WriteNavData -1
#define WriteGeoAttributes -2

typedef enum {
    IMG_RGEO,
    MOD_RGEO,
    IMG_RGEO_TC,
    MOD_RGEO_TC,
    DNB_GEO,
    IMG_GEO,
    MOD_GEO,
    MOD_UNAGG_GEO,
    IMG_GEO_TC,
    MOD_GEO_TC,
    MAX_NUM_GEO_OUTPUTS
} GeoOutputs;

typedef struct {
    float imgLat[4];
    float imgLon[4];
    float imgLonRanges[4]; // minNeg, maxNeg, minPos, maxPos
    float imgNSEW[4];
    float modLat[4];
    float modLon[4];
    float modLonRanges[4];
    float modNSEW[4];
    float dnbLat[4];
    float dnbLon[4];
    float dnbLonRanges[4];
    float dnbNSEW[4];
} gRingArrays;

// Structure of pointers used by geolocation

typedef struct {
    // Pointer to VCST Navigation data
    VcstCmnGeo* pCmnGeo;

    // moonArray is generated for the calibration routine, and is local
    // to SDR processing
    ViirsMoonArrays moonArray;

    GEO_param_struct* geoParams;

    // global data used in processing of geolocation
    ViirsGeoProcType* procStruct;

    ViirsGeoOutputType* sdrGeo;
    //  ViirsGeoOutputType* sdrDGeo; // Radian form products pointer(s).

    //Array will keep track of the quality metadata values needed for
    //Geolocation products
    int geoMetadataInfo[MAX_NUM_GEO_OUTPUTS][MAX_NUM_GEO_QUAL];

    int autoQual;

    // QA LUT contains the SAA threshold
    proSdrViirsCalQALUT* QALUT;

    // Servo Control Checks
    set<string> geoChecks;

} viirsSdrGeoPtrs;


class VcstCmnGeo;

class VcstViirsGeo {
public:
    // Granule data
    string platform_;
    string startDirection_;
    string endDirection_;
    string day_night_flag_;
    string time_coverage_start_;
    string time_coverage_end_;
    string pge_start_time_;
    string pge_end_time_;
    string versionid_;
    string history_;
    string source_files_;
    int filled_scans_;
    int orbit_number_;
    int format_version_;
    int instrument_number_;
    int extract_pixel_start_;
    int extract_pixel_stop_;
    int leapseconds93;

    // scan line attributes:
    double scan_time_fill_;
    double scan_time_valid_min_;
    double scan_time_valid_max_;

    static string DISABLE_SERVO_CONTROL;

public:
    /**
     * Constructor
     */
    VcstViirsGeo();

    /**
     * Destructor
     */
    virtual ~VcstViirsGeo();

    /**
     * Delete all allocated memory
     */
    void cleanupDataItems();

    /**
     * Initialize geolocation object class
     */

    int initialize();


    /**
     * Set the history class member
     */
    //    void setHistory(int argc, char* argv[]);

    /**
     * geolocate - executes geolocation algorithm
     */

    int geolocate();

    /**
     * findNSEW - find the bounding box for input coords
     */
    int findNSEW(float** lat, float** lon,
            size_t ny, size_t x0, size_t nx, float NSEW[4], float lonRanges[4]);

    /**
     * findGring - find valid corners of input coord arrays
     */
    int findGring(float** lat, float** lon,
            size_t ny, size_t x0, size_t nx,
            float gringLat[4], float gringLon[4]);

    /**
     * updateBbox - update the bounding box and gring with each scan
     */
    int updateBbox(GeoOutputs geotype, size_t x0, size_t nx);

    /**
     * write_img - write img geolocation file
     */

    int write_img(int iscan, viirsSdrGeoPtrs* geoPtrs_);

    /**
     * write_mod - write mod geolocation file
     */

    int write_mod(int iscan, viirsSdrGeoPtrs* geoPtrs_, bool modFromImg);

    /**
     * write_dnb - write dnb geolocation file
     */

    int write_dnb(int iscan, viirsSdrGeoPtrs* geoPtrs_);

    void setHistory(string history) {
        history_ = history;
    }

    string getHistory() {
        return history_;
    }

    void setSource(string source) {
        source_files_ = source;
    }

    string getSource() {
        return source_files_;
    }

protected:

    /**
     * doProcessing - IDPS/ADL doProcessing
     */

    int doProcessing();


    /**
     * Initialize L1A data
     */

    int initialize_L1A_data();


    /**
     * Initialize L1B data
     */

    int initialize_L1B_data();

    /**
     * getLandWaterMask()
     * read LWM file and populate the landwater array
     */
    int getLandWaterMask(GeoOutputs geotype, size_t x0, size_t nx);

    /**
     * write_global_attributes()
     */

    int write_global_attributes(NcFile* nc_output, string type);


    /**
     * write_scan_data()
     */

    int write_scan_data(NcFile* nc_output);


    /**
     * write_nav_data()
     */

    int write_nav_data(NcFile* nc_output);


private:

    /**
     * Removes allocated memory for temporary data
     */
    void cleanup();

    /**
     * Geolocates one granule.  It creates interpolation rectangles and
     * calculates the geolocation for points in each decimated interpolation
     * rectangle and then uses quadratic interpolation to provide geolocation
     * for all points in an interpolation rectangle. The geolocation is then
     * stored to the local 'fullGeo' buffer and later copied into the actual DMS
     * buffers.  However, note that height is put directly into the DMS buffer.
     */
    int geolocateGranule(int iscan,
            viirsSdrGeoPtrs* ptrs,
            VcstCmnGeo* geoPtr,
            ViirsGeoRctnglType* geoInterpRctngl,
            ViirsGeoDecimType* vInt,
            bool modFromImg);

    /**
     * This function stores the geolocation data to the DMS buffers.
     */
    int storeGranule(int scan,
            const ViirsGeoFullType& inFull,
            const ViirsGeoRctnglType& inRect,
            viirsSdrGeoPtrs* ptrs);

    /**
     * This function will calculate the Moderate geolocation data from the
     * Imagery geolocation data.  This method is also known as "nesting".
     *
     * @param ptrs    Geo struct data to use
     * @param geoPtr  Common GEO class to use
     */
    void calcModFromImg(int scan,
            viirsSdrGeoPtrs* ptrs,
            VcstCmnGeo *geoPtr,
            int extractPixelLimits[2]);

    /**
     * convertToDegrees() - convert arrays of radian angles to degrees.
     */
    int convertToDegrees(int iscan, char type);

    /**
     * geolocation ptrs
     */
    viirsSdrGeoPtrs* geoPtrs_;

    gRingArrays* gring;
};


/**
 * geolocatePixel
 * 
 * @param *ptrs  Structure of pointers to VIIRS SDR geolocation structs
 * @param inScan Scan number of pixel to be geolocated
 * @param inDet Detector number of pixel to be geolocated
 * @param inCol Column number of pixel to be geolocated
 * @param geoPtr visibility to a VcstCmnGeo instance
 * @param &outPix Structure containg geolocation data for a pixel
 *
 * @return PRO_FAIL or PRO_SUCCESS
 */
int geolocatePixel(viirsSdrGeoPtrs* ptrs,
        const int inScan,
        const int inDet,
        const int inCol,
        VcstCmnGeo *geoPtr,
        ViirsGeoPixelType& outPix);

/**
 * setGeoPixelQuality
 * 
 * @param *ptrs Structure of pointers to VIIRS SDR geolocation structs
 * @param inRow Row number of pixel quality to be set
 * @param inCol Column number of  pixel quality to be set
 * @param pixelQuality quality of the pixel
 *
 */
void setGeoPixelQuality(viirsSdrGeoPtrs* ptrs,
        const int inRow,
        const int inCol,
        const unsigned char pixeQuality);

/**
 * initGeoDataStructs
 * 
 * @param *ptrs  Structure of pointers to VIIRS SDR geolocation structs
 *
 * @return pro_fail or pro_success
 */
int initGeoDataStructs(int iscan, viirsSdrGeoPtrs* ptrs);

/**
 * storeGranule
 *
 * @param inFull Structure containing full interpoloation rectangle data
 * @param *ptrs  Structure of pointers to VIIRS SDR geolocation structs
 * @param env    data environment, for tasking ProCmnDataItem world
 *
 * @return void
 */
int storeGranule(const ViirsGeoFullType& inFull,
        const ViirsGeoRctnglType& inRect,
        viirsSdrGeoPtrs* ptrs);

/**
 * geolocateAllRecPix
 *
 * @param *ptrs  Structure of pointers to VIIRS SDR geolocation structs
 * @param inRec  Structure containing decimated interpolation rectangle data
 * @param geoPtr vis to a VcstCmnGeo object
 * @param ioFull Structure containing full interpoloation rectangle data
 *
 * @return void
 */
void geolocateAllRecPix(int scan,
        viirsSdrGeoPtrs* ptrs,
        ViirsGeoRctnglType* inRec,
        VcstCmnGeo* geoPtr,
        ViirsGeoFullType* ioFull,
        int extractPixelLimits[2]);

/**
 * correctFromPassThroughPI
 *
 * @param ioValue  Value to check for passing through PI and update it
 *
 * @return void
 */
void correctFromPassThroughPI(float* ioValue);

void correctFromPassThroughPI(double* ioValue);

/**
 * adjustForPassThroughPI
 *
 * @param io0 First interpolation point
 * @param io1 Second interpolation point
 * @param io2 Third interpolation point
 *
 * @return void
 */
void adjustForPassThroughPI(float* io0,
        float* io1,
        float* io2);

void adjustForPassThroughPI(double* io0,
        double* io1,
        double* io2);


/**
 * calculateNadirData
 *
 * @param inScan          Current scan number
 * @param ptrs            Pointer to viirsSdrGeoPtrs
 * @param geoPtr          vis to a VcstCmnGeo instance
 * @param outFull         Structure containing geolocation for all pixels
 *
 * @return PRO_SUCCESS, PRO_FAIL
 */
int calculateNadirData(const int inScan,
        viirsSdrGeoPtrs* ptrs,
        VcstCmnGeo* geoPtr,
        ViirsGeoFullType* outFull);


/**
 * GEO_validate_scan_encoder_data
 * 
 * @param number_of_scans  number of scans in the granule
 * @param *ptrs  Structure of pointers to VIIRS SDR geolocation structs
 *
 * @return PRO_FAIL or PRO_SUCCESS
 */
int GEO_validate_scan_encoder_data(viirsSdrGeoPtrs* ptrs);

/**
 * GEO_process_parameters
 * 
 * @param *ptrs  Structure of pointers to VIIRS SDR geolocation structs
 *
 * @return PRO_SUCCESS
 */
int GEO_process_parameters(GEO_param_struct* geoParams,
        ViirsGeoProcType* procStruct);


/**
 * GEO_absolute_limit_check
 * 
 * @param data_samples[]     array of input data samples
 * @param number_of_samples  number of samples to check
 * @param data_limits[2]     lower and upper limits
 * @param sample_flags[]     array of validation flags
 *
 * @return PRO_FAIL or PRO_SUCCESS
 */
int GEO_absolute_limit_check(double data_samples[],
        int const number_of_samples,
        double data_limits[2],
        int sample_flags[]);

/**
 * GEO_relative_limit_check
 * 
 * @param data_samples[]     array of input data samples
 * @param number_of_samples  number of samples to check
 * @param delta_limit        difference limit
 * @param sample_flags[]     validation flags indicating flagged samples
 *
 * @return PRO_FAIL or PRO_SUCCESS
 */
int GEO_relative_limit_check(double data_samples[],
        int const number_of_samples,
        double const delta_limit,
        int sample_flags[]);

/**
 * GEO_find_next_flag
 *  
 * @param sample_flags[]     array of flags for data samples
 * @param number_of_samples  array size
 * @param start_sample       start index for search
 *
 * @return index of next unflagged value or GEO_FAIL
 */
int GEO_find_next_flag(int sample_flags[],
        int const number_of_samples,
        int const start_sample);

/**
 * GEO_determine_DNB_sample_time_offsets
 * 
 * @param *ptrs  Structure of pointers to VIIRS SDR geolocation structs
 *
 * @return void
 */
void GEO_determine_DNB_sample_time_offsets(
        const focal_plane_geometry_struct& geometry_params,
        const int actScans,
        const double t_sync[MAX_SCAN_NUMBER],
        double sample_time[MAX_SCAN_NUMBER][MAX_SCAN_SAMPLE]);

/**
 * GEO_determine_sample_time_offsets 
 * 
 * This method is only used for PPC
 * The PPC algorithm uses the sample times to recalculate the view vector
 *
 * @param *ptrs  Structure of pointers to VIIRS SDR geolocation structs
 *
 * @return void
 */
void GEO_determine_sample_time_offsets(
        const focal_plane_geometry_struct& geometry_params,
        double sample_time[MAX_SCAN_SAMPLE]);

/******************************************************************
 * GEO_Upg_determine_sample_time_offsets
 *
 * @param *ptrs  Structure of pointers to VIIRS SDR geolocation structs
 *
 * @return void
 * *****************************************************************/
void GEO_Upg_determine_sample_time_offsets(
        const focal_plane_geometry_struct& geometry_params,
        const int actScans,
        const double t_sync[MAX_SCAN_NUMBER],
        double sample_time[MAX_SCAN_NUMBER][MAX_SCAN_SAMPLE]);

/**
 * GEO_determine_thermal_corrections
 * 
 * @param *ptrs  Structure of pointers to VIIRS SDR geolocation structs
 *
 * @return PRO_FAIL or PRO_SUCCESS
 */
int GEO_determine_thermal_corrections(viirsSdrGeoPtrs* ptrs);


/**
 * GEO_determine_view_vectors
 * 
 * @param scan_number              scan number
 * @param sample_number            sample number in the scan
 * @param *ptrs   Structure of pointers to VIIRS SDR geolocation structs
 * @param *flag   quality flag
 * @param u_inst[MAX_DETECTORS][3] view vectors in instrument coordinates
 *
 * @return PRO_FAIL or PRO_SUCCESS
 */
int GEO_determine_view_vectors(const int scan_number,
        const int sample_number,
        const int det,
        viirsSdrGeoPtrs* ptrs,
        unsigned char* flag,
        double u_inst[VEC_SIZE]);

/**
 * GEO_interpolate_telescope_encoder
 * 
 * @param scan_number    scan number
 * @param sample_number  sample number
 * @param *sample_enc    telescope encoder value for sample
 * @param *ptrs  Structure of pointers to VIIRS SDR geolocation structs
 *
 * @return PRO_FAIL or PRO_SUCCESS
 */
int GEO_interpolate_telescope_encoder(int const scan_number,
        int const sample_number,
        double* const sample_enc,
        viirsSdrGeoPtrs* ptrs);

/**
 * GEO_interpolate_mirror_encoder
 * 
 * @param scan_number    scan number
 * @param sample_number  sample number
 * @param *sample_enc    mirror encoder value for sample
 * @param *ptrs  Structure of pointers to VIIRS SDR geolocation structs
 *
 * @return PRO_FAIL or PRO_SUCCESS
 */
int GEO_interpolate_mirror_encoder(int const scan_number,
        int const sample_number,
        double* const sample_enc,
        viirsSdrGeoPtrs* ptrs);

/**
 * GEO_evaluate_polynomial
 * 
 * @param coef[MAX_POLY_DEGREE+1]  array of coefficients for the polynomial
 * @param in_x                     the input independent variable x
 * @param degree                   the degree of the polynomial
 * @param *out_fit                 the output value
 *
 * @return PRO_FAIL or PRO_SUCCESS
 */
int GEO_evaluate_polynomial(double coef[MAX_POLY_DEGREE + 1],
        double const in_x,
        int const degree,
        double* const out_fit);


/**
 * GEO_interp_mod_unagg
 *
 * @param scan   scan number
 * @param *ptrs  Structure of pointers to VIIRS SDR geolocation structs
 *
 * @return PRO_FAIL or PRO_SUCCESS
 */
int GEO_interp_mod_unagg(viirsSdrGeoPtrs* ptrs);


//-----------------------------------------------------------------------

inline void correctFromPassThroughPI(float* ioValue) {
    if (*ioValue > FLOAT32_FILL_TEST) {
        if (*ioValue > TWOPI) {
            *ioValue -= TWOPI;
        }

        if (*ioValue > PI) {
            *ioValue -= TWOPI;
        } else if (*ioValue < -PI) {
            *ioValue += TWOPI;
        }
    }
}

inline void correctFromPassThroughPI(double* ioValue) {
    if (*ioValue > FLOAT64_FILL_TEST) {
        if (*ioValue > TWOPI) {
            *ioValue -= TWOPI;
        }

        if (*ioValue > PI) {
            *ioValue -= TWOPI;
        } else if (*ioValue < -PI) {
            *ioValue += TWOPI;
        }
    }
}

//-----------------------------------------------------------------------

inline void adjustForPassThroughPI(float* io0,
        float* io1,
        float* io2) {
    float diff;

    if ((*io0 > FLOAT32_FILL_TEST) && (*io1 > FLOAT32_FILL_TEST) &&
            (*io2 > FLOAT32_FILL_TEST)) {
        diff = fabs(*io2 - *io0);
        if (diff > PI) {
            if (*io0 < 0.0) {
                *io0 += TWOPI;
            }

            if (*io1 < 0.0) {
                *io1 += TWOPI;
            }

            if (*io2 < 0.0) {
                *io2 += TWOPI;
            }
        }
    }
}

inline void adjustForPassThroughPI(double* io0,
        double* io1,
        double* io2) {
    double diff;

    if ((*io0 > FLOAT64_FILL_TEST) && (*io1 > FLOAT64_FILL_TEST) &&
            (*io2 > FLOAT64_FILL_TEST)) {
        diff = fabs(*io2 - *io0);
        if (diff > PI) {
            if (*io0 < 0.0) {
                *io0 += TWOPI;
            }

            if (*io1 < 0.0) {
                *io1 += TWOPI;
            }

            if (*io2 < 0.0) {
                *io2 += TWOPI;
            }
        }
    }
}

/**
 * eclipseAndSaaFlags
 *
 * @param geoPtrs
 * @param cmnGeoPtr
 *
 * @return int PRO_SUCCESS or PRO_FAIL
 */
int eclipseAndSaaFlags(viirsSdrGeoPtrs *geoPtrs,
        VcstCmnGeo *cmnGeoPtr);

/**
 *getSaa
 *
 *@param sdrPtrs
 *@param dLat 
 *@param dLon 
 *@param cmnGeoPtr
 *@param inSAA
 *
 *@return int PR0_SUCCESS or PRO_FAIL
 */
int getSaa(viirsSdrGeoPtrs* geoPtrs,
        float dLat,
        float dLon,
        VcstCmnGeo *cmnGeoPtr,
        unsigned char &inSAA);

/**
 *getSolarEclipse
 *
 *@param lat array of latitude
 *@param lon array of longitude
 *@param startTime scan start time
 *@param endTime scan end time
 *@param startRow index for the first row
 *@param endRow index for the last row
 *@param cmnGeoPtr
 *@param solEclipse out param
 *
 *@return int PRO_SUCCESS or PRO_FAIL
 */
template <int T_NUM_ROWS, int T_NUM_COLS>
int getSolarEclipse(float lat[T_NUM_ROWS][T_NUM_COLS],
        float lon[T_NUM_ROWS][T_NUM_COLS],
        long long startTime,
        long long endTime,
        int startRow,
        int endRow,
        VcstCmnGeo *cmnGeoPtr,
        unsigned char &solEclipse);


/**
 *getDnbLunarEclipse
 *
 *@param geoPtrs
 *@param startTime
 *@param endTime
 *@param startRow
 *@param endRow
 *@param cmnGeoPtr
 *@param lunEclipse
 *
 *@return int PRO_SUCCESS or PRO_FAIL
 */
int getDnbLunarEclipse(viirsSdrGeoPtrs *geoPtrs,
        double startTime,
        double endTime,
        int endRow,
        int startRow,
        VcstCmnGeo *cmnGeoPtr,
        unsigned char &lunEclipse);

/**
 * Check_Tel_Start_Not_Nominal
 *
 * @param *ptrs       structure of pointers to input and output items
 * 
 * @return PRO_FAIL or PRO_SUCCESS
 */
int Check_Tel_Start_Not_Nominal(viirsSdrGeoPtrs* ptrs);


/**
 * Check_Ham_Start_Not_Nominal
 *
 * @param *ptrs       structure of pointers to input and output items
 * 
 * @return PRO_FAIL or PRO_SUCCESS
 */
int Check_Ham_Start_Not_Nominal(viirsSdrGeoPtrs* ptrs);

/**
 * Check_Scan_Sync_Failure
 *
 * @param *ptrs       structure of pointers to input and output items
 * 
 * @return PRO_FAIL or PRO_SUCCESS
 */
int Check_Scan_Sync_Failure(viirsSdrGeoPtrs* ptrs);

/**
 * Determine the electronic side of VIIRS from OBC-IP
 *
 * @param *ptrs       structure of pointers to input and output items
 * 
 * @return PRO_FAIL or PRO_SUCCESS
 */
int calculateSCESide(viirsSdrGeoPtrs* ptrs);


/**
 * Set flag bit in Quality Flag 2 in geo struct.
 *
 * @param *ptrs    structure of pointers to input and output items
 *
 * @param flag     unsigned char flag defined in viirs_QualityFlags.h PRO_FAIL or PRO_SUCCESS
 */
int appendQualityFlag2(viirsSdrGeoPtrs* ptrs, int scan, unsigned char flag);

/**
 * Builds the data describing the interpolation rectangles
 *
 */
int VcstCreateInterpRectangles(const int bandType,
        int numZones,              //# of agg zones
        int *numPixPerZone,        //# of pixels in each agg zone
        ViirsGeoRctnglType* outRec);


#endif

