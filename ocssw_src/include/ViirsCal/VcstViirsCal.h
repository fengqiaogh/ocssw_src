/*******************************************************************************
 *
 * NAME: VcstViirsCal.h
 *
 * DESCRIPTION: Object class that provides application programming interfaces
 * for initialization, calibration, and L1B creation.
 *
 * Created on: February 25, 2015
 *     Author: Sam Anderson, VCST
 *
 *******************************************************************************/

#ifndef VcstViirsCal_h
#define VcstViirsCal_h

#include <VcstObc.h>
#include <VcstParamsReader.h>

class VcstViirsCal {
public:

    bool bRad_;
    bool bCdg_;
    bool processModByScan_;
    string history_;
    string source_files_;

    static char const* DISABLE_SERVO_CONTROL;

    static const char VIIRS_NIGHT = 0;
    static const char VIIRS_DAY = 1;
    static const char VIIRS_MIXED = 2;

    /**
     * Constructor
     */

    VcstViirsCal();

    /**
     * Destructor
     */

    virtual ~VcstViirsCal();

    /**
     * initialize()
     *
     * Create and initialize instance of granule object class
     */

    int initialize(bool bRad, bool bFilt, bool bCDG, bool bLunar);

    /**
     * calibrate()
     *
     * Sets up and executes the calibration routine by band and format
     */

    int calibrate(VIIRS_CATEGORY_ENUM bands);

    /**
     * Turn CDG generation on/off
     */

    inline const bool doCdg() {
        return bCdg_;
    };

    inline const void setCdg(const bool bCdg) {
        bCdg_ = bCdg;
    };

    /**
     * set/get history and source()
     *
     * Set source and history
     */

    inline void setHistory(std::string history) {
        history_ = history;
    }

    inline std::string getHistory() {
        return history_;
    }

    inline void setSource(std::string source) {
        source_files_ = source;
    }

    inline std::string getSource() {
        return source_files_;
    }

    int calibrateMOD(int iScan, float **l1bptrs);

private:

    /**
     * Granule object pointer.
     */

    VcstObc* pObc_;

    /**
     * Lists to store identities of bands to be calibrated
     */

    std::vector<VIIRS_BAND_ENUM> imgBands_;
    std::vector<VIIRS_BAND_ENUM> modBands_;
    std::vector<VIIRS_BAND_ENUM> dnbBands_;

    /**
     * create_band_lists()
     *
     * Create list of bands to be calibrated and reported
     */

    int create_band_lists(VIIRS_CATEGORY_ENUM bands);

    /**
     * calibrate_img_bands()
     *
     * Calibrates all imagery bands in list, and writes data out to L1B
     * files in formats specified by the enum argument
     */

    int calibrate_img_bands();

    /**
     * calibrate_mod_bands()
     *
     * Calibrates all moderate bands in list, and writes data out to L1B
     * files in formats specified by the enum argument
     */

    int calibrate_mod_bands();

    /**
     * calibrate_dnb_bands()
     *
     * Calibrates dnb if it is in list, and writes data out to L1B
     * files in format specified by the enum argument
     */

    int calibrate_dnb_bands();

    /**
     *  Write global attributes to L1B file.
     */

    int write_global_attributes(NcFile* nc_output, string type);

    /**
     *  Create imagery resolution band L1B output file and write global data.
     *  Return the output file's netCDF4 ID.
     */

    const NcFile* write_img();

    /**
     *  Create moderate resolution band L1B output file and write global data.
     *  Return the output file's netCDF4 ID.
     */

    const NcFile* write_mod();

    /**
     *  Create unaggregated dual-gain output file and write global data.
     *  Return the output file's netCDF4 ID.
     */

    const NcFile* write_cdg();

    /**
     *  Create day/night band L1B output file and write global data.
     *  Return the output file's netCDF4 ID.
     */

    const NcFile* write_dnb();

    // scan line attributes
    double scan_time_fill_value_;
    double scan_time_valid_min_;
    double scan_time_valid_max_;

    // scan state flag attributes
    static const short NUM_SCAN_STATE_FLAGS = 3;
    unsigned char scan_state_fill_value_;
    vector<unsigned char> scan_state_flag_masks_;

    // scan quality flag attributes
    static const short NUM_SCAN_QUALITY_FLAGS = 7;
    unsigned char scan_quality_fill_value_;
    vector<unsigned char> scan_quality_flag_masks_;

};

#endif

