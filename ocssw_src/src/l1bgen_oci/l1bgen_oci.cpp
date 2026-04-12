/**
 * @file l1bgen_oci.cpp
 * @brief OCI Level 1B data generation
 *
 * This file contains the main program for generating OCI Level 1B data from Level 1A data.
 * It performs calibration, geolocation, and other processing steps to convert raw instrument
 * data into calibrated radiances or reflectances.
 *
 * @authors Joel Gales (SAIC), Jakob Lindo (SSAI)
 * @date 02/09/2024
 *
 * @version 0.870
 *
 * @details
 * This program reads OCI Level 1A data, applies various calibration and correction factors,
 * performs geolocation, and generates Level 1B data products. It handles blue, red, and SWIR
 * bands separately, applying appropriate calibration coefficients and corrections for each.
 * The program also handles metadata, applies quality flags, and writes the processed data
 * to a new netCDF file in the Level 1B format.
 *
 * Key processing steps include:
 * - Reading calibration lookup tables
 * - Applying dark offset corrections
 * - Applying temperature and gain corrections
 * - Handling non-linearity and RVS (Response Versus Scan) corrections
 * - Aggregating hyperspectral bands
 * - Calculating solar irradiance corrections
 * - Handling saturation and quality flags
 * - Writing processed data and metadata to output file
 *
 */

#include <netcdf>
#include <math.h>
#include <string.h>
#include <algorithm>

#include <global_attrs.h>
#include <allocate2d.h>
#include <allocate3d.h>

#include "corrections.hpp"
#include "calibrations.hpp"
#include "aggregate.hpp"
#include "geolocate_oci.h"
#include "corrections.hpp"
#include "calibrations.hpp"
#include "l1b_file.hpp"
#include "processing_tracker.hpp"
#include "l1b_options.hpp"
#include "types.hpp"

using namespace std;
using namespace netCDF;
using namespace netCDF::exceptions;

template <typename T>
using matrix2D = boost::multi_array<T, 2>;

#define VERSION "2.8"
#define CHUNK_CACHE_SIZE 256 * 1024 * 1024  // 256MiB of cache memory.
#define CHUNK_CACHE_NELEMS 1033
#define CHUNK_CACHE_PREEMPTION .75
#define VARCHUNK_CACHE_SIZE 4 * 1024 * 1024  // 4Mib of cache memory.
#define VARCHUNK_CACHE_NELEMS 1033
#define VARCHUNK_CACHE_PREEMPTION .75
#define CHUNKBANDS 40
#define CHUNKPIXELS 2000
#define CHUNKLINES 16
#define EXP_DECAY_CONSTS 4

constexpr size_t PIXELS_NOMINAL = 1329;     // L1A pixels differ from L1B, b/c l1bgen_oci aggregates
constexpr size_t BLUE_BANDS_NOMINAL = 120;  // Pre-agg
constexpr size_t RED_BANDS_NOMINAL = 168;   // Pre-agg
constexpr size_t SWIR_BANDS_NOMINAL = 9;
constexpr size_t PIXELS_NBSB = 2884;
constexpr size_t BLUE_BANDS_NBSB = 60;
constexpr size_t RED_BANDS_NBSB = 60;
constexpr size_t SWIR_BANDS_NBSB = 9;

// Macros defining the views OCI has. The metadata produced from l1agen_oci differs from the OAD for
// l1bgen_oci
#define MAIN_VIEW 1  // The view from which usable science data is derived
#define DARK_VIEW 2  // Where dark data comes from

//    Modification history:
//  Programmer     Organization   Date     Ver   Description of change
//  ----------     ------------   ----     ---   ---------------------
//  Joel Gales     SAIC           04/29/20 0.01  Original development
//  Joel Gales     SAIC           01/13/21 0.02  Add support for SWIR
//  Joel Gales     SAIC           08/11/21 0.03  Add support for handling
//                                               fill values in science data
//  Joel Gales     SAIC           08/12/21 0.04  Add support for hysteresis
//                                               correction
//  Joel Gales     SAIC           09/23/21 0.045 Initialize uninitialized
//                                               variables
//  Joel Gales     SAIC           05/29/23 0.060 Implemented F.Patt code changes
//                                               (05/08/23) #272
//  Joel Gales     SAIC           06/21/23 0.061 Fix hysteresis bugs
//  Joel Gales     SAIC           06/23/23 0.062 Move get_agg_mat to common.cpp
//                                               Split off get_nib_nbb
//  Gwyn Fireman   SAIC           07/10/23 0.063  Read global metadata from json file
//  Joel Gales     SAIC           07/24/23 0.063 Add support for saturation
//  Joel Gales     SAIC           08/09/23 0.065 Add check_scan_times
//  Joel Gales     SAIC           08/25/23 0.700 Add support for reflectance
//                                               output
//  Joel Gales     SAIC           08/29/23 0.710 Convert thetap,thetas todeg
//  Joel Gales     SAIC           09/29/23 0.800 Call geolocation as function
//                                               rather than run beforehand
//  Joel Gales     SAIC           10/16/23 0.810 Add planarity correction
//  Joel Gales     SAIC           10/19/23 0.820 Fix metadata issues
//  Joel Gales     SAIC           12/06/23 0.823 Add rhot description attribute
//  Joel Gales     SAIC           12/10/23 0.840 Add tilt_home
//  Joel Gales     SAIC           01/02/24 0.860 Fix encoder interpolation bug
//  Joel Gales     SAIC           02/09/24 0.870 Add scan angle & polarization
//                                               coefficients

bool nonBaselineBandsAndPixels(NcFile *l1aFile) {
    const size_t filePixels = l1aFile->getDim("ccd_pixels").getSize();
    const size_t fileBlueBands = l1aFile->getDim("blue_bands").getSize();
    const size_t fileRedBands = l1aFile->getDim("red_bands").getSize();
    const size_t fileSwirBands = l1aFile->getDim("SWIR_bands").getSize();

    if (filePixels == PIXELS_NBSB && fileBlueBands == BLUE_BANDS_NBSB && fileRedBands == RED_BANDS_NBSB &&
        fileSwirBands == SWIR_BANDS_NBSB) {
        return true;
    } else if (filePixels == PIXELS_NOMINAL && fileBlueBands == BLUE_BANDS_NOMINAL &&
               fileRedBands == RED_BANDS_NOMINAL && fileSwirBands == SWIR_BANDS_NOMINAL) {
        return false;
    } else {
        std::cerr << "-E- Malformed pixels/bands dimension(s) in input L1A\n   "
                  << "Pixels: " << filePixels << " | Blue bands " << fileBlueBands << " | Red bands "
                  << fileRedBands << " | SWIR bands " << fileSwirBands << endl;
        exit(EXIT_FAILURE);
    }
}

int main(int argc, char *argv[]) {
    oel::L1bOptions options(argc, argv, VERSION);

    cout << clo_getVersion() << endl;

    if (argc == 1) {
        clo_printUsage(options.optionList);
        exit(EXIT_FAILURE);
    }

    // ********************************* //
    // *** Read calibration LUT file *** //
    // ********************************* //
    nc_set_chunk_cache(CHUNK_CACHE_SIZE, CHUNK_CACHE_NELEMS, CHUNK_CACHE_PREEMPTION);
    NcFile *calLutFile = new NcFile(options.calibrationLutFilename, NcFile::read);

    NcGroup calLutCommon, calLutBlue, calLutRed, calLutSwir;
    calLutCommon = calLutFile->getGroup("common");
    calLutBlue = calLutFile->getGroup("blue");
    calLutRed = calLutFile->getGroup("red");
    calLutSwir = calLutFile->getGroup("SWIR");

    vector<float> blueWavelengths(NUM_BLUE_WAVELENGTHS);       // Band centers from blue CCD
    vector<float> redWavelengths(NUM_RED_WAVELENGTHS);         // Band centers from red CCD
    vector<float> swirWavelengths(NUM_SWIR_WAVELENGTHS);       // Band centers from SWIR
    vector<double> blueSolarIrradiance(NUM_BLUE_WAVELENGTHS);  // Light in these wavelengths from the Sun
    vector<double> redSolarIrradiance(NUM_RED_WAVELENGTHS);    // Light in these wavelengths from the Sun
    vector<double> swirSolarIrradiance(NUM_SWIR_WAVELENGTHS);  // Light in these wavelengths from the Sun

    NcDim k2TimeDim = calLutFile->getDim("number_of_times");
    if (k2TimeDim.isNull()) {
        cout << "Error: could not read number_of_times from " << options.calibrationLutFilename << "\n";
        exit(EXIT_FAILURE);
    }
    size_t numK2Times = k2TimeDim.getSize();

    double *K2t = (double *)malloc(numK2Times * sizeof(double));

    calLutCommon.getVar("blue_wavelength").getVar(blueWavelengths.data());
    calLutCommon.getVar("blue_F0").getVar(blueSolarIrradiance.data());
    calLutCommon.getVar("red_wavelength").getVar(redWavelengths.data());
    calLutCommon.getVar("red_F0").getVar(redSolarIrradiance.data());
    calLutCommon.getVar("SWIR_wavelength").getVar(swirWavelengths.data());
    calLutCommon.getVar("SWIR_F0").getVar(swirSolarIrradiance.data());
    calLutCommon.getVar("K2t").getVar(K2t);

    uint32_t numBlueGainBands = 0, numRedGainBands = 0, numSwirGainBands = 0;
    NcDim K3Tdim = calLutFile->getDim("number_of_temperatures");
    float *tempCoefs = new float[K3Tdim.getSize()];  // K3T
    calLutCommon.getVar("K3T").getVar(tempCoefs);

    CalibrationLut blueCalLut = readOciCalLut(calLutFile, BLUE, calLutBlue, numBlueGainBands, 1);
    CalibrationLut redCalLut = readOciCalLut(calLutFile, RED, calLutRed, numRedGainBands, 1);
    CalibrationLut swirCalLut = readOciCalLut(calLutFile, SWIR, calLutSwir, numSwirGainBands, 2);

    // Read hysterisis parameters
    float hysterisisTimes[9][4];      // number of SWIR bands by number exponential functions
    float hysterisisAmplitude[9][4];  // number of SWIR bands by number exponential functions
    calLutSwir.getVar("hyst_time_const").getVar(&hysterisisTimes[0][0]);
    calLutSwir.getVar("hyst_amplitude").getVar(&hysterisisAmplitude[0][0]);

    calLutFile->close();
    delete (calLutFile);

    NcFile *xtkLutFile;
    size_t ncpix;  // number of influence pixels
    size_t nbands;
    size_t nxbands;
    float ***cmat_in;  // input crosstalk influence coeff matrix

    if (options.enableCrosstalk) {
        // Open and read data from Xtalk LUT file
        xtkLutFile = new NcFile(options.xtalkLutFilename, NcFile::read);

        ncpix = xtkLutFile->getDim("number_of_pixels").getSize();  // number of influence pixels
        nbands = xtkLutFile->getDim("number_of_bands").getSize();
        nxbands = xtkLutFile->getDim("number_of_xc_bands").getSize();

        cmat_in =
            allocate3d_float(ncpix, nxbands, nbands);  // blue band input crosstalk influence coeff matrix
        xtkLutFile->getVar("blue_cross_coef").getVar(&cmat_in[0][0][0]);

        xtkLutFile->close();
        delete (xtkLutFile);
    }

    static Level1bFile outfile(options.l1bFilename);
    // Append call sequence to existing history
    string history = call_sequence(argc, argv);

    // Open and read data from L1Afile
    NcFile *l1aFile = new NcFile(options.l1aFilename, NcFile::read);
    NcGroup l1aSpatSpecModes = l1aFile->getGroup("spatial_spectral_modes");
    NcGroup l1aSciData = l1aFile->getGroup("science_data");
    size_t numSpatialZones = l1aFile->getDim("spatial_zones").getSize();
    vector<short> dataTypes(numSpatialZones);
    l1aSpatSpecModes.getVar("spatial_zone_data_type").getVar(dataTypes.data());

    LocatingContext locatingContext;
    switch (dataTypes.at(MAIN_VIEW)) {
        case EARTH:
            if (nonBaselineBandsAndPixels(l1aFile)) {
                cout << "Input is high spatial resolution, turning off aggregation" << endl;
                options.aggregationOff = true;
            } else {
                cout << "Input was a normal science file" << endl;
            }
            locatingContext = GEO;
            break;

        case SOLAR_DAILY_CALIBRATION:
        case SOLAR_MONTHLY_CALIBRATION:
            cout << "Input was a solar calibration granule" << endl;
            locatingContext = HELIO;
            options.disableGeolocation = true;
            options.radianceGenerationEnabled = true;
            break;

        case LUNAR:
            cout << "Input was a lunar calibration granule" << endl;
            locatingContext = SELENO;
            break;

        default:
            cout << "-W- Unknown locating context, assuming GEO" << endl;
            locatingContext = GEO;
            break;
    }

    GeoData geoData;
    oel::ProcessingTracker processingTracker;

    try {
        geoData = locateOci(l1aFile, outfile, options, locatingContext, processingTracker);
    } catch (const exception& e) {
        cerr << "-E- Couldn't geolocate file:\n    " << e.what() << endl;
        exit(100);
    }

    cout << "Preparing to calibrate file" << endl;

    vector<int16_t> numLinesSpatZones(numSpatialZones);  // Number of lines per spatial zone
    vector<int16_t> spatialAgg(numSpatialZones);         // Spatial aggregation factors per zone
    l1aSpatSpecModes.getVar("spatial_zone_data_type").getVar(dataTypes.data());
    l1aSpatSpecModes.getVar("spatial_zone_lines").getVar(numLinesSpatZones.data());
    l1aSpatSpecModes.getVar("spatial_aggregation").getVar(spatialAgg.data());

    bool foundSpatialAggregation = false;
    size_t spatialAggregationIndex;
    // spatialAggregationIndex is looking for the index that has EARTH datatype, and not DIAGNOSTIC datatype
    // the latter always appears at a higher index value, so the loop looks for the lowest index that is
    // has datatype not equal to NO_DATA, DARK_CALIBRATIOM or EXTERNAL_SNAPSHOP_TRIGGER.
    // Presumably could simply look for index with datatype equal to EARTH instead.
    for (spatialAggregationIndex = 0; spatialAggregationIndex < dataTypes.size(); spatialAggregationIndex++) {
        int16_t dataType = dataTypes[spatialAggregationIndex];
        if (dataType != NO_DATA && dataType != DARK_CALIBRATION && dataType != EXTERNAL_SNAPSHOP_TRIGGER) {
            foundSpatialAggregation = true;
            break;
        }
    }

    if (!foundSpatialAggregation) {
        cout << "-E- " << __FILE__ << ":" << __LINE__
             << " - Could not find valid spatial aggregation type.\n";
        exit(EXIT_FAILURE);
    }

    CalibrationData blueCalData(l1aFile, geoData, BLUE, numBlueGainBands, numK2Times, K2t, spatialAgg,
                                spatialAggregationIndex, blueCalLut, blueSolarIrradiance,
                                options.aggregationOff);
    CalibrationData redCalData(l1aFile, geoData, RED, numRedGainBands, numK2Times, K2t, spatialAgg,
                               spatialAggregationIndex, redCalLut, redSolarIrradiance,
                               options.aggregationOff);
    CalibrationData swirCalData(l1aFile, geoData, SWIR, numSwirGainBands, numK2Times, K2t, spatialAgg,
                                spatialAggregationIndex, swirCalLut, swirSolarIrradiance,
                                options.aggregationOff);

    DarkData darkData;

    vector<double> cosineSolarZeniths(geoData.solarZeniths.size());

    if (!options.radianceGenerationEnabled) {
        transform(geoData.solarZeniths.begin(), geoData.solarZeniths.end(), cosineSolarZeniths.begin(),
                  [](double zenith) { return cos(zenith * OEL_DEGRAD); });
    }

    // Get date
    cout << "time_coverage_start: " << unix2isodate(geoData.unixTimeStart, 'G') << endl;
    cout << "time_coverage_end:   " << unix2isodate(geoData.unixTimeEnd, 'G') << endl;

    // ******************************************** //
    // *** Get spatial and spectral aggregation *** //
    // ******************************************** //
    uint32_t numTaps = l1aFile->getDim("number_of_taps").getSize();
    // Number of spatial zones was obtained above for locating context, as were the data types

    // numBands is assigned after aggregation
    blueCalData.saturatedPixelCounts.resize(blueCalData.numBands);
    redCalData.saturatedPixelCounts.resize(redCalData.numBands);
    swirCalData.saturatedPixelCounts.resize(swirCalData.numBands);

    if (options.enableCrosstalk) { // Crosstalk correction is only needed for blue bands right now
        blueCalData.enableCrosstalk = options.enableCrosstalk;
        makeXtalkMat(spatialAgg[spatialAggregationIndex], blueCalData, cmat_in, ncpix, nxbands, nbands);
    }

    // ********************************* //
    // *** Get dark collect location *** //
    // ********************************* //
    for (size_t i = 0; i < numSpatialZones; i++) {
        if (dataTypes[i] == DARK_CALIBRATION) {
            darkData.darkZone = (int16_t)i;
        }
    }
    if (darkData.darkZone == -1) {
        cout << "No dark collect in file: " << options.l1aFilename.c_str() << endl;
        l1aFile->close();
        return 1;
    }

    int16_t sciDarkCount = 0;   // Number of dark science lines after aggregation
    int16_t swirDarkCount = 0;  // Number of dark SWIR lines after aggregation
    for (size_t i = 0; i < (size_t)darkData.darkZone; i++) {
        if (dataTypes[i] != 0 && dataTypes[1] != 10) {
            sciDarkCount += numLinesSpatZones[i] / spatialAgg[i];
            swirDarkCount += numLinesSpatZones[i] / 8;
        }
    }
    darkData.numPix = numLinesSpatZones[darkData.darkZone] / spatialAgg[darkData.darkZone];
    int16_t numDarkSwirPix = numLinesSpatZones[darkData.darkZone] / 8;

    // *********************************************************************** //
    // *** Generate band gain structs from LUTs, date/time & gain matrices *** //
    // *********************************************************************** //

    free(K2t);

    // Read selected temperature fields and interpolate to scan times
    uint16_t numTemps = blueCalData.gains->dimensions[TEMP] + swirCalData.gains->dimensions[TEMP];

    vector<vector<double>> calibrationTemps =
        interpTemps(l1aFile, numTemps, geoData.numGoodScans, geoData.earthViewTimes.data());
    uint16_t nBCalTemps = blueCalData.gains->dimensions[TEMP];  // Number of blue calibration temps

    // Read dark collects from science data arrays
    vector<size_t> start(3), count(3);  // for use in vectorized netcdf calls
    start[0] = 0;
    start[1] = 0;
    start[2] = sciDarkCount;

    uint32_t ***blueSciPixDark =  // Blue dark sicence pixels
        (uint32_t ***)allocate3d_int(geoData.numGoodScans, blueCalData.numInsBands, darkData.numPix);
    uint32_t ***redSciPixDark =  // Red dark science pixels
        (uint32_t ***)allocate3d_int(geoData.numGoodScans, redCalData.numInsBands, darkData.numPix);
    uint32_t ***swirPixDark =
        (uint32_t ***)allocate3d_int(geoData.numGoodScans, swirCalData.numBands, numDarkSwirPix);

    if (blueCalData.numInsBands > 4) {
        count[0] = geoData.numGoodScans;
        count[1] = blueCalData.numInsBands;
        count[2] = darkData.numPix;

        l1aSciData.getVar("sci_blue").getVar(start, count, &blueSciPixDark[0][0][0]);

        filterDarkNoise(geoData.numGoodScans, blueCalData.numBands, darkData.numPix, blueCalData.fillValue,
                        blueSciPixDark);
    }
    if (redCalData.numInsBands > 4) {
        count[0] = geoData.numGoodScans;
        count[1] = redCalData.numInsBands;
        count[2] = darkData.numPix;

        l1aSciData.getVar("sci_red").getVar(start, count, &redSciPixDark[0][0][0]);

        filterDarkNoise(geoData.numGoodScans, redCalData.numBands, darkData.numPix, redCalData.fillValue,
                        redSciPixDark);
    }

    start[0] = 0;
    start[1] = 0;
    start[2] = swirDarkCount;

    count[0] = geoData.numGoodScans;
    count[1] = swirCalData.numBands;
    count[2] = numDarkSwirPix;

    l1aSciData.getVar("sci_SWIR").getVar(start, count, &swirPixDark[0][0][0]);
    filterDarkNoise(geoData.numGoodScans, swirCalData.numBands, numDarkSwirPix, swirCalData.fillValue,
                    swirPixDark);

    // Calibrated data variables
    vec2D<double> swirDigiNum(swirCalData.numBands, vector<double>(geoData.numSwirPix));
    double **swirCal = allocate2d_double(swirCalData.numBands, geoData.numSwirPix);
    uint8_t **swirQualityFlags = allocate2d_uchar(swirCalData.numBands, geoData.numSwirPix);
    vector<double> sciScanAngles(geoData.numCcdPix);
    vector<double> swirScanAngles(geoData.numSwirPix);
    vector<double> swirDarkCorrections(swirCalData.numBands);
    uint32_t **swirSciData = (uint32_t **)allocate2d_int(swirCalData.numBands, geoData.numSwirPix);
    int16_t swirAggFactor = 1;
    vector<double> swirTempCorrections(swirCalData.numBands);
    vector<double> hysteresisCorrection(geoData.numSwirPix);
    int32_t *indx = new int32_t[geoData.numSwirPix];

    double pprAngle =
        (2 * OEL_PI) *
        (geoData.mceTelem.pprOffset - geoData.geoLut.rtaNadir[geoData.mceTelem.mceBoardId % 2]) /
        MAX_ENC_COUNT;
    if (pprAngle > OEL_PI)
        pprAngle -= 2 * OEL_PI;

    ///////////////
    // Main loop //
    ///////////////

    cout << "Calibrating" << endl;

    // Read, calibrate and write science data
    for (size_t currScan = processingTracker.getStartingLine(); currScan < processingTracker.getEndingLine();
         currScan++) {
        processingTracker.update();

        // Check for valid mirror side
        if (!(geoData.hamSides[currScan] == 0 || geoData.hamSides[currScan] == 1)) {
            cout << "No mirror side index for scan: " << currScan << endl;
            continue;
        }

        // initalize quality flags
        bzero(blueCalData.qualityFlags, blueCalData.numBands * geoData.numCcdPix);
        bzero(redCalData.qualityFlags, redCalData.numBands * geoData.numCcdPix);

        // set 10 deg pixels to HAM-B_striping, if normal science data
        if (dataTypes.at(MAIN_VIEW) == EARTH) {
            if (geoData.numCcdPix > 800) {
                if (geoData.hamSides[currScan]) {
                    for (size_t band = 0; band < redCalData.numBands; band++) {
                        uint8_t *ptr = redCalData.qualityFlags + band * geoData.numCcdPix + 680;
                        for (size_t pix = 680; pix < 800; pix++) {
                            *ptr = HAM_B_STRIPING;
                            ptr++;
                        }
                    }
                }
            }
        }

        //  Get scan angle
        size_t pixelOffset = currScan * geoData.numCcdPix;
        for (size_t pix = 0; pix < geoData.numCcdPix; pix++) { 
            double angleRadians = geoData.ccdScanAngles[pixelOffset + pix];

            if (angleRadians == BAD_FLT) {
                continue;
            } else {
                sciScanAngles[pix] = angleRadians * OEL_RADEG;
            }
        }

        pixelOffset = currScan * geoData.numSwirPix;
        for (size_t pix = 0; pix < geoData.numSwirPix; pix++) {
            double angleRadians = geoData.swirScanAngles[pixelOffset + pix];

            if (angleRadians == BAD_FLT) {
                continue;
            } else {
                swirScanAngles[pix] = geoData.swirScanAngles[pixelOffset + pix] * OEL_RADEG;
            }
        }

        start[0] = currScan;
        start[1] = 0;

        count[0] = 1;
        count[1] = geoData.numCcdPix;

        outfile.navigationData.getVar("CCD_scan_angles").putVar(start, count, sciScanAngles.data());

        count[1] = geoData.numSwirPix;
        outfile.navigationData.getVar("SWIR_scan_angles").putVar(start, count, swirScanAngles.data());

        //  Blue bands
        if (blueCalData.numInsBands >= 4) {
            darkData.data = blueSciPixDark;

            calibrate(blueCalData, geoData, darkData, l1aSciData, outfile, currScan, spatialAgg,
                      spatialAggregationIndex, numTaps, blueCalLut.dimensions[4], sciScanAngles, tempCoefs,
                      calibrationTemps, cosineSolarZeniths, options.radianceGenerationEnabled);
        }

        //  Red bands
        if (redCalData.numInsBands >= 4) {
            darkData.data = redSciPixDark;

            calibrate(redCalData, geoData, darkData, l1aSciData, outfile, currScan, spatialAgg,
                      spatialAggregationIndex, numTaps, redCalLut.dimensions[4], sciScanAngles, tempCoefs,
                      calibrationTemps, cosineSolarZeniths, options.radianceGenerationEnabled);
        }

        //  SWIR bands
        start[0] = currScan;
        start[1] = 0;
        start[2] = 0;

        count[0] = 1;
        count[1] = swirCalData.numBands;
        count[2] = geoData.numSwirPix;

        l1aSciData.getVar("sci_SWIR").getVar(start, count, swirSciData[0]);

        // Compute dark offset, correct data, and apply absolute and
        // temporal gain and temperature correction

        int darkCorrRetStat = getDarkCorrection(
            currScan, geoData.numGoodScans, geoData.hamSides, darkData.numScansAvg, darkData.numPixSkip, 1, 1,
            1, {swirAggFactor}, swirCalData.fillValue, numDarkSwirPix, swirPixDark, swirDarkCorrections);

        if (darkCorrRetStat != -1) {
            // calc swirDigiNum
            for (size_t j = 0; j < swirCalData.numBands; j++) {
                for (size_t pix = 0; pix < geoData.numSwirPix; pix++) {
                    if (swirSciData[j][pix] == swirCalData.fillValue) {
                        swirDigiNum[j][pix] = BAD_FLT;
                    } else {
                        swirDigiNum[j][pix] = swirSciData[j][pix] - swirDarkCorrections[j];
                    }
                }
            }

            // When there are more CCD pixels than SWIR, there are inconsistencies in geolocation.
            // This is fixed by resampling solar zeniths to the SWIR pixels.
            vector<double> cosSolZensSwir(geoData.numSwirPix);
            if (geoData.numCcdPix > geoData.numSwirPix) {
                for (size_t i = 0; i < geoData.numSwirPix; i++) {
                    size_t j = lower_bound(sciScanAngles.begin(), sciScanAngles.end(), swirScanAngles[i]) -
                               sciScanAngles.begin();
                    if (j > 0 && j < geoData.numCcdPix) {
                        double w = (swirScanAngles[i] - sciScanAngles[j - 1]) -
                                   (sciScanAngles[j] - sciScanAngles[j - 1]);
                        cosSolZensSwir[i] =
                            cosineSolarZeniths[currScan * geoData.numCcdPix + j - 1] * (1 - w) +
                            cosineSolarZeniths[currScan * geoData.numCcdPix + j] * w;
                    } else {
                        cosSolZensSwir[i] = cosineSolarZeniths[currScan * geoData.numCcdPix + j];
                    }
                }

            } else {
                copy(cosineSolarZeniths.begin() + currScan * geoData.numSwirPix,
                     cosineSolarZeniths.begin() + (currScan + 1) * geoData.numSwirPix,
                     cosSolZensSwir.begin());
            }

            // SWIR band loop
            for (size_t band = 0; band < swirCalData.numBands; band++) {
                uint32_t numGoodPixels = 0;
                for (size_t pix = 0; pix < geoData.numSwirPix; pix++) {
                    bool notFill = swirSciData[band][pix] != swirCalData.fillValue;
                    // radiance == true, do not check geoData.qualFlag bc it does not retrieve it
                    if (options.radianceGenerationEnabled && notFill) {
                        indx[numGoodPixels++] = pix;
                    }
                    // radiance == false, check quality flag for only "Off_Earth, Solar_eclipse, Terrain_bad"
                    // (0x7)
                    else if (!options.radianceGenerationEnabled && notFill &&
                             (!(geoData.qualityFlag[currScan * geoData.numCcdPix + pix] & 0x7))) {
                        indx[numGoodPixels++] = pix;
                    } else {
                        // Handle fill value
                        swirDigiNum[band][pix] = BAD_FLT;
                        swirCal[band][pix] = BAD_FLT;
                    }
                }

                // Hysteresis correction
                double prior[EXP_DECAY_CONSTS] = {0, 0, 0, 0};
                double observation[EXP_DECAY_CONSTS] = {0, 0, 0, 0};
                hysteresisCorrection[0] = 0.0;
                array<double, EXP_DECAY_CONSTS> exponentialDecayConsts;
                for (size_t k = 0; k < EXP_DECAY_CONSTS; k++) {
                    exponentialDecayConsts[k] = exp(-1.0 / hysterisisTimes[band][k]);
                }
                vector<double> deltaScanAngles(swirScanAngles.size());  // in IDL, dtheta
                for (size_t i = 0; i < swirScanAngles.size() - 1; i++) {
                    deltaScanAngles[i] = swirScanAngles[i + 1] - swirScanAngles[i];
                }

                for (size_t pix = 1; pix < numGoodPixels; pix++) {
                    hysteresisCorrection[pix] = 0.0;

                    if (deltaScanAngles[pix - 1] < 0.1) {
                        for (size_t decayConst = 0; decayConst < EXP_DECAY_CONSTS; decayConst++) {
                            observation[decayConst] =
                                prior[decayConst] * exponentialDecayConsts[decayConst] +
                                swirDigiNum[band][indx[pix - 1]] * hysterisisAmplitude[band][decayConst];
                            hysteresisCorrection[pix] += observation[decayConst];
                            prior[decayConst] = observation[decayConst];
                        }
                    } else {
                        for (size_t decayConst = 0; decayConst < EXP_DECAY_CONSTS; decayConst++) {
                            prior[decayConst] = 0.0;
                            observation[decayConst] = 0.0;
                        }
                        hysteresisCorrection[pix] = 0.0;
                    }
                }

                double tempCorrection = getTempCorrection(
                    band, &tempCoefs[nBCalTemps], calibrationTemps[currScan], *swirCalData.gains, SWIR);
                for (size_t k = 0; k < numGoodPixels; k++) {
                    double rvsCorrection = getRvsCorrection(band, indx[k], geoData.hamSides[currScan],
                                                            *swirCalData.gains, swirScanAngles[k]);
                    double nonLinearityCorrection = getNonlinearityCorrection(
                        band, indx[k], swirCalLut.dimensions[NL], swirCalData.gains->k5Coefs, swirDigiNum);
                    swirCal[band][indx[k]] = tempCorrection *
                                          swirCalData.gains->k1K2[band][geoData.hamSides[currScan]] *
                                          (swirDigiNum[band][indx[k]] - hysteresisCorrection[k]);

                    swirCal[band][indx[k]] *= nonLinearityCorrection * rvsCorrection;

                    if (!options.radianceGenerationEnabled) {
                        // NOTE: solz is short with scale of 0.01
                        if (geoData.solarZeniths[currScan * geoData.numCcdPix + indx[k]] < MAX_SOLZ)
                            swirCal[band][indx[k]] *= OEL_PI * geoData.auCorrection /
                                                   (swirSolarIrradiance[band] * cosSolZensSwir[indx[k]]);
                        else
                            swirCal[band][indx[k]] = BAD_FLT;
                    }
                }  // (reduced) k-loop

            }  // band loop

        }  // iret != -1

        // Check for saturation
        for (size_t k = 0; k < geoData.numSwirPix; k++) {
            for (size_t band = 0; band < swirCalData.numBands; band++) {
                if (swirSciData[band][k] >= swirCalData.gains->saturationThresholds[band]){
                    swirQualityFlags[band][k] = 1;
                    swirCalData.saturatedPixelCounts[band] += 1;
                }
                else
                    swirQualityFlags[band][k] = 0;
            }
        }

        start[0] = 0;
        start[1] = currScan;
        start[2] = 0;

        count[0] = swirCalData.numBands;
        count[1] = 1;
        count[2] = geoData.numSwirPix;

        // Output to L1B file
        if (!options.radianceGenerationEnabled) {
            try {
                outfile.observationData.getVar("rhot_SWIR").putVar(start, count, &swirCal[0][0]);
            } catch (const exception &e) {
                cout << "-E- " << __FILE__ << ":" << __LINE__
                     << " - Couldn't put SWIR rhot values into L1B file" << endl;
                cout << e.what() << endl;
                exit(EXIT_FAILURE);
            }
        } else {
            try {
                outfile.observationData.getVar("Lt_SWIR").putVar(start, count, &swirCal[0][0]);
            } catch (const exception &e) {
                cout << "-E- " << __FILE__ << ":" << __LINE__
                     << " - Couldn't put SWIR Lt values into L1B file" << endl;
                cout << e.what() << endl;
                exit(EXIT_FAILURE);
            }
        }
        try {
            outfile.observationData.getVar("qual_SWIR").putVar(start, count, &swirQualityFlags[0][0]);
        } catch (const exception &e) {
            cout << "-E- " << __FILE__ << ":" << __LINE__
                 << " - Couldn't put SWIR quality flags into L1B file" << endl;
            cout << e.what() << endl;
            exit(EXIT_FAILURE);
        }

    }  // Scan loop

    // End Main loop

    cout << "Calibration finished, writing band info and granule metadata" << endl;

    // bgmat#bwave
    if (blueCalData.numInsBands >= 4) {
        blueCalData.writeBandInfo(outfile, blueCalLut, blueWavelengths);
        blueCalData.writeSaturationPercent(outfile, geoData);
    }

    if (redCalData.numInsBands >= 4) {
        redCalData.writeBandInfo(outfile, redCalLut, redWavelengths);
        redCalData.writeSaturationPercent(outfile, geoData);
    }

    swirCalData.writeBandInfo(outfile, swirCalLut, swirWavelengths);
    swirCalData.writeSaturationPercent(outfile, geoData);

    outfile.writeGranuleMetadata(unix2isodate(geoData.unixTimeStart, 'G'),
                                 unix2isodate(geoData.unixTimeEnd, 'G'), options.l1bFilename);

    set_global_attrs(outfile.l1bFile, history, options.digitalObjectId, options.processingVersion);

    outfile.writeProcessingControl(options);

    l1aFile->close();
    delete (l1aFile);

    cout << "Done" << endl;

    return 0;
}
