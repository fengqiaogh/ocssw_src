constexpr double SCAN_TIME_INTERVAL = 1.4771;
constexpr double MAX_ATTERR = 0.017453293;  // radians, ~= 1 degree

// MODIS GEO variables
enum GeoVars {
    // "Geolocation Fields" stored as float
    LONGITUDE,
    LATITUDE,
    // "Data Fields" stored as scaled short
    HEIGHT,
    SOLAR_ZENITH,
    SOLAR_AZIMUTH,
    SENSOR_ZENITH,
    SENSOR_AZIMUTH,
    NUM_GEO_VARS,
};

// MODIS L1B variables
enum L1bVars {
    RSB_250,  // bands 1 & 2
    RSB_500,  // bands 3 - 7
    RSB_1KM,  // all other reflective bands
    CIR_1KM,  // band 26
    TEB_1KM,  // all thermal bands
    NUM_L1B_VARS,
};

/*
  Notes:
  All  GEO arrays are dimensioned [nscans*ndets, nframes]
  Most L1B arrays are dimensioned [nbands, nscans*ndets, nframes*nsf]
  where:
  nscans and nframes are at 1KM resolution
  nsf = #subframes = 1000/resolution
  ndets = 10*nsubframes

  Band 26 drops the 1st dimension.
  Geolocation is inherently 1KM resolution: nbands=1, ndets=10, nsf=1.

  ndets and scandim are needed for reading data one scan at a time.

  There is some redundancy between resolution, ndets, and nsf,
  but it's helpful to have them all defined.

*/

struct BandInfo {
    short varIndex;        // which variable contains the band
    short bandIndex;       // index of band within variable
    short flagIndex;       // starting index within "all-detector" arrays
    short wavelength;      // (nm)
    std::string bandName;  // for clarity
};

// Bands used by l2gen, in expected order
std::vector<BandInfo> bandList = {
    // Reflective Solar Bands
    {RSB_1KM, 0, 180, 412, "8"},     //  0
    {RSB_1KM, 1, 190, 443, "9"},     //  1
    {RSB_500, 0, 80, 469, "3"},      //  2
    {RSB_1KM, 2, 200, 488, "10"},    //  3
    {RSB_1KM, 3, 210, 531, "11"},    //  4
    {RSB_1KM, 4, 220, 547, "12"},    //  5
    {RSB_500, 1, 100, 555, "4"},     //  6
    {RSB_250, 0, 0, 645, "1"},       //  7
    {RSB_1KM, 5, 230, 667, "13lo"},  //  8
    {RSB_1KM, 7, 250, 678, "14lo"},  //  9
    {RSB_1KM, 9, 270, 748, "15"},    // 10
    {RSB_250, 1, 40, 859, "2"},      // 11
    {RSB_1KM, 10, 280, 869, "16"},   // 12
    {RSB_500, 2, 120, 1240, "5"},    // 13
    {RSB_500, 3, 140, 1640, "6"},    // 14
    {RSB_500, 4, 160, 2130, "7"},    // 15
    // Thermal Emissive Bands
    {TEB_1KM, 0, 330, 3750, "20"},    // 16
    {TEB_1KM, 2, 350, 3959, "22"},    // 17
    {TEB_1KM, 3, 360, 4050, "23"},    // 18
    {TEB_1KM, 6, 390, 6715, "27"},    // 19
    {TEB_1KM, 7, 400, 7325, "28"},    // 20
    {TEB_1KM, 8, 410, 8550, "29"},    // 21
    {TEB_1KM, 10, 430, 11000, "31"},  // 22
    {TEB_1KM, 11, 440, 12000, "32"},  // 23
    // Cirrus Band
    {CIR_1KM, 0, 320, 1375, "26"}  // 24
};

// Sentinel values for bad EV data
enum BadVals {
    L1A_SCAN_DATA_MISSING_SI = 65535,
    L1A_DN_MISSING_SI = 65534,
    SATURATED_DETECTOR_SI = 65533,
    ZERO_POINT_DN_SI = 65532,
    DEAD_DETECTOR_SI = 65531,
    RSB_DN_STAR_BELOW_MIN_SI = 65530,
    TEB_OR_RSB_GT_MAX_SI = 65529,
    AGGREGATION_FAIL_SI = 65528,
    SECTOR_ROTATION_SI = 65527,
    TEB_B1_NOT_CALCULATED = 65526,
    DEAD_SUBFRAME_SI = 65525,
    UNABLE_CALIBRATE_SI = 65524,
    UNRESCALED_HIGH_SI = 65521,
    RESCALED_L1B_SI = 65520,
    NAD_CLOSED_UPPER_SI = 65500,
};
constexpr int MIN_BAD_SI = NAD_CLOSED_UPPER_SI;
constexpr int MAX_GOOD_SI = 32767;
