// Polarization Correction for OCI
// File is the same as polcor_hawkeye.cpp, but with some modifications

#include <vector>
#include <map>
#include <netcdf>
#include <math.h>
#include <allocate2d.h>
#include "l12_proto.h"
#include "polcor_oci.h"

using namespace std;
using namespace netCDF;
using namespace netCDF::exceptions;

// global variables that will be used repeately per scan
static bool firstRun = true;            // first run to grab m12, m13 coef and ham sides
static bool skipPolarization = false;   // when ccd pix != swir pix, dont do polarization

static int currScanLine = -1;
static int numPolCoef = -1;

static int numScans = -1;
static int numCcdPixels = -1;
static int numSwirPixels = -1;
static int numBands = -1; 
static int numBlueBands = -1;
static int numRedBands = -1;
static int numSwirBands = -1;

// grab m12 and m13 values on first run. The key will be the ham side.
static map<int, vector<vector<float>>> m12;
static map<int, vector<vector<float>>> m13;

// grab the SWIR quality flag for the current line. The keys
// will be the band number and gets upated every new line
static map<int, vector<unsigned char>> qualSwir;

// updates every new scan line 
static vector<unsigned char> hamSide;
static vector<float> ccdScanAngle;

/**
 *  Given an NcVar reference to a blue or red M12 or M13 Coef, 
 *  get all the coef. for the bands specified and append them 
 *  to an array that is aggrigating the blue, red and SWIR M12
 *  and M13 coefs. 
 *  @param ncVar - Ref to the band variable from NetCDF
 *  @param aggVec - vector to put the coefs
 *  @param start - start range for bands to be read
 *  @param stop - stop range for band to be read
 *  @param numPolCoef - number of polorization coef that is being read
 *  @param currHamSide - ham side that the coef is on
 */
void getPolarizationCoefs(NcVar &ncVar, vector<vector<float>> &aggVec, 
                    int start, int stop, int numPolCoef, int currHamSide) {
    
    for (int band = start; band < stop; band++) {
                
        vector<float> coefs(numPolCoef, 0.0);

        // start at the current band, current ham side for this line and first element of coefs
        vector<size_t> start = {
            (size_t)band,               // start index of band
            (size_t)currHamSide,        // start index of ham side
            0                           // start index of pol coef.
        }; 

        // want 1 line for current band, 1 line of ham side and all of the coefs (max 3)
        vector<size_t> count = {
            (size_t) 1,                          // read 1 line of band
            (size_t) 1,                          // read 1 line of ham side
            (size_t)numPolCoef          // read 3 items of coefs.
        };

        // read it into coefs
        ncVar.getVar(start, count, coefs.data());

        // add band to m12 and m13
        aggVec.push_back(coefs);

    }
}


/**
 * Function that grabs the m12 and m13 variables and set them 
 * in the static variables above for all scans.
 * It also grabs static dim values that will be used for the duration of the program.
 * @param ncFile - reference the opened NetCDF file 
 */
void initializePolCorOci(NcFile &ncFile) {

    try {
        // read in the dimensions used for m12 and m13. Also read in dims
        // that will be used to fetch other variables later on
        numBlueBands = ncFile.getDim("blue_bands").getSize();
        numRedBands = ncFile.getDim("red_bands").getSize();
        numSwirBands = ncFile.getDim("SWIR_bands").getSize();
        numPolCoef = ncFile.getDim("polarization_coefficients").getSize();
        numCcdPixels = ncFile.getDim("ccd_pixels").getSize();
        numScans = ncFile.getDim("number_of_scans").getSize();
        numSwirPixels = ncFile.getDim("SWIR_pixels").getSize();

        // OCI has 3 coefficients. Bad file if it doesn't. 
        if (numPolCoef != 3 || numCcdPixels != numSwirPixels) {
            printf("-E- %s:%d - Problem with the number of polarization coefficients. Should be 3, but got %d\n", __FILE__, __LINE__, numPolCoef);
            exit(EXIT_FAILURE);
        }
        // If ccd pixels != swir pixels, polarization will not happen
        // -- Will fix and add if we need it to work when ccd pix != swir pix ``
        if (numCcdPixels != numSwirPixels) {
            printf("-WARNING- %s:%d - Number of CCD pixels and SWIR pixels don't match. CCD: %d vs. SWIR: %d. Polarization will be 1 for all pixels.\n", __FILE__, __LINE__, numCcdPixels, numSwirPixels);
            skipPolarization = true;
            return;
        }

        /*----READING HAM SIDES----*/ 

        numScans = ncFile.getDim("number_of_scans").getSize();
        hamSide = vector<unsigned char>(numScans, 255);

        // read in HAM side
        NcVar hamVar = ncFile.getGroup("scan_line_attributes").getVar("HAM_side");
        hamVar.getVar(hamSide.data());


        /*----READING BLUE, RED, SWIR M12 AND M13----*/
        NcGroup sensorBandGroup = ncFile.getGroup("sensor_band_parameters");

        NcVar blueM12Var = sensorBandGroup.getVar("blue_m12_coef");
        NcVar blueM13Var = sensorBandGroup.getVar("blue_m13_coef"); 

        NcVar redM12Var = sensorBandGroup.getVar("red_m12_coef");
        NcVar redM13Var = sensorBandGroup.getVar("red_m13_coef");

        NcVar swirM12Var = sensorBandGroup.getVar("SWIR_m12_coef");
        NcVar swirM13Var = sensorBandGroup.getVar("SWIR_m13_coef");

        // grab m12 and m13 values for both ham sides 0 and 1 
        for (int ham = 0; ham < 2; ham++) {

            // -- BLUE -- skip last 3 because it overlaps with red
            getPolarizationCoefs(blueM12Var, m12[ham], 0, numBlueBands-3, numPolCoef, ham);
            getPolarizationCoefs(blueM13Var, m13[ham], 0, numBlueBands-3, numPolCoef, ham);

            // -- RED -- reads in all
            getPolarizationCoefs(redM12Var, m12[ham], 0, numRedBands, numPolCoef, ham);
            getPolarizationCoefs(redM13Var, m13[ham], 0, numRedBands, numPolCoef, ham);

            // -- SWIR -- read all, but 2 will be dropped during calculation
            getPolarizationCoefs(swirM12Var, m12[ham], 0, numSwirBands, numPolCoef, ham); 
            getPolarizationCoefs(swirM12Var, m13[ham], 0, numSwirBands, numPolCoef, ham);   
            
        }
        
    } catch (NcException& e) {
        printf("-E- %s:%d - Problem reading polfile %s\n", __FILE__, __LINE__, input->polfile);
        printf("            %s\n", e.what());
        exit(EXIT_FAILURE);
    }  
}

/**
 * For each new line, get the CCD scan angle for all pixels
 */
void getScanAngleForScan(NcFile &ncFile) {

    // clear the previous scan line's CCD scan angle for each pixel
    ccdScanAngle.clear(); 

    // initalize it to be the size of the number of pixels
    ccdScanAngle = vector<float>(numCcdPixels, -32767.0);

    // read in the angles
    NcVar ccdScanAngleVar = ncFile.getGroup("navigation_data").getVar("CCD_scan_angles");

    // start at scan line and the first pixel
    vector<size_t> start = {(size_t)currScanLine, 0};
    // read 1 line and all CCD pixels
    vector<size_t>count = {1, (size_t)numCcdPixels};
    // save it to vector
    ccdScanAngleVar.getVar(start, count, ccdScanAngle.data());


}

/**
 * For each new line, get the line's SWIR quality flag for low and high gain bands
 * for all pixels. The quality flag describes the bands saturation (T/F) for pix.
 * low gain bands = 2, 5
 * high gain bands = 3, 6
 * @param ncFile - ref to the current opened NC file.
 */

void getQualSwirForScan(NcFile &ncFile) {
    // clear previous scans data
    qualSwir.clear();

    NcVar qualSwirVar = ncFile.getGroup("observation_data").getVar("qual_SWIR");

    for (int band = 2; band < 7; band++) {

        // not the band we want, skip
        if (band != 2 && band != 3 && band != 5 && band != 6) {
            continue;
        }

        // everything processed here should be bands: 2, 3, 5 and 6
        // save the band's qual flag here
        vector<unsigned char> qualSwirForBand(numSwirPixels, 255);

        vector<size_t> start = {
            (size_t) band,              // start at curr band
            (size_t) currScanLine,      // start at curr scan line
            0                           // start from beginning of pixels
        };

        // how many from the start to read in
        vector<size_t> count = {
            (size_t) 1,                          // entire band
            (size_t) 1,                          // entire scan line
            (size_t) numSwirPixels      // all pixels 
        };

        // read it into the vector
        qualSwirVar.getVar(start, count, qualSwirForBand.data());

        // save it to the map
        qualSwir[band] = qualSwirForBand;

    }

}

/**
 * Grab the current lines quality flag for SWIR bands.
 * Check the high gain bands (3 and 6) for saturation.
 * @param currPix - current pixel
 * @return true if at least 1 is saturated, false if not.
 */
bool areHighGainBandsSaturated(int currPix) {

    // cast to int because reading in from NetCDF file is unsigned byte
    int band3Sat = (int) qualSwir[3][currPix];
    int band6Sat = (int) qualSwir[6][currPix];

    // true if either are saturated.
    return band3Sat || band6Sat;
}

/**
 * Calculates the polarization correction for the current pixel and sets it to index.
 * 
 * Index is the relative position of the current pixel to the polcor array, which 
 * has size numBands * numPixels for the current scan.
 * 
 * 
 * Calculate the value of m12 and m13 using the Second Degree Polynomial Function:
 *  Ax^2 + Bx + C, where:
 *      A = coef[2], 
 *      B = coef[1],
 *      C = coef[0],
 *      x = ccd scan angle
 * 
 */
void calculateAndSetPolCor(l1str* &l1rec, int currPix, int index, vector<float> m12Coefs, vector<float> m13Coefs) {

    double alpha = l1rec->alpha[currPix] / RADEG;
    double L_x = l1rec->Lt[index] / l1rec->tg_sol[index] / l1rec->tg_sen[index];
    double L_qp = l1rec->L_q[index] * cos(2 * alpha) + l1rec->L_u[index] * sin(2 * alpha);
    double L_up = l1rec->L_u[index] * cos(2 * alpha) - l1rec->L_q[index] * sin(2 * alpha);

    /*
        If L_qp and L_up are NaN, skip. polcor == 1.0.
        --
        l1rec->L_u and L_q are NaN when airmass is negative because log(negative) is undefined.
        l1rec->L_u and L_q are calculated in rayleigh.c and it uses the function ray_press_wang to
        calculate a fraction. If airmass is negative, that fraction is NaN since it uses log of airmass. 
        
        This instance happens when solar_zenith (solz) is >= 90, resulting in a negative cos after
        it is converted into radian. This value carries over when calculating airmass.
    */
    if (isnan(L_qp) || isnan(L_up)) return;

    // get the scan angle for current pixel and the coef 
    float scanAngle = ccdScanAngle[currPix];
    
    // calculate m12
    float A = m12Coefs[2];
    float B = m12Coefs[1];
    float C = m12Coefs[0];
    double m12 = (A * pow(scanAngle, 2.0)) + (B * scanAngle) + C;

    // calculate m13
    A = m13Coefs[2];
    B = m13Coefs[1];
    C = m13Coefs[0];
    double m13 = (A * pow(scanAngle, 2.0)) + (B * scanAngle) + C;

    // calculate polcor 
    l1rec->polcor[index] = 1.0 / (1.0 - m12 * L_qp / L_x - m13 * L_up / L_x);
    l1rec->dpol[index] = sqrt(pow(l1rec->L_q[index], 2.0) + pow(l1rec->L_u[index], 2.0)) / L_x;
}

/**
 * For the current scan, grab and merge the blue and red m12 and m13 bands.
 * Also grab the swir bands based on saturation (max 2).
 * Then, for every band, apply the polarization correction for currPix
 * @param l1rec - L1B file -- contains the current scan line this pixel belongs to
 * @param currPix - current pixel to apply polcor
 * 
 */
extern "C" void polcor_oci(l1str *l1rec, int32_t currPix) {

    // happens when ccd pixel != swir pixels. leave polcof as 1.000
    if (skipPolarization) {
        return;
    }

    // open NC file
    NcFile ncFile(input->ifile[0], NcFile::read);

    // num bands from the l1file is 286 after dropping 3 bands from blue
    // and 2 from SWIR 
    numBands = l1rec->l1file->nbands;

    // grab all the m12 and m13 values for blue, red and SWIR
    // also get ham side so we know which m12 or m13 side to use
    if (firstRun) {
        firstRun = false;
        initializePolCorOci(ncFile);
    }

    // new scan line, need to grab new CCD_Scan_Angles
    if (currScanLine < l1rec->iscan) {
        currScanLine = l1rec->iscan;
        getScanAngleForScan(ncFile);
        getQualSwirForScan(ncFile);
    }

    
    // ---- APPLYING POLARIZATION CORRECTION ----


    // track the number of bands that were used in calculating polcor because
    // m12 and m13 are merged
    int bandCount = 0; 
    int currHamSide = hamSide[currScanLine];


    // Process Blue Bands first, skipping last 3 bands
    // NOTE: ORDER MATTERS -- when merging the band coefs, the order is blue, red, swir
    // Skip last 3 bands because it overlaps with red. Use red over blue. m12 and m13 
    // already skips the blue bands, but not the SWIR ones.
    for (int currBand = 0; currBand < numBlueBands-3; currBand++) {

        // relative position for this pixel in the polcor array.
        // polcor = numBands * numpixels for the scan line.
        // to get the index of currPix relative to it, it will be 
        // currpix * numBands + curr band position.
        // not using currBand but bandCount bc some bands are skipped and because
        // this calculation is done for each band color
        int relativeIndex = currPix * numBands + bandCount;
        vector<float> m12Coef = m12[currHamSide][bandCount];
        vector<float> m13Coef = m13[currHamSide][bandCount];
        calculateAndSetPolCor(l1rec, currPix, relativeIndex, m12Coef, m13Coef);
        bandCount++;
    }


    // Process Red Bands
    for (int currBand = 0; currBand < numRedBands; currBand++) {

        int relativeIndex = currPix * numBands + bandCount;
        // m12 and m13 are merged, so use the current band count to get the red band
        vector<float> m12Coef = m12[currHamSide][bandCount];
        vector<float> m13Coef = m13[currHamSide][bandCount];
        calculateAndSetPolCor(l1rec, currPix, relativeIndex, m12Coef, m13Coef);
        bandCount++;
    }


    // check if the current pixels SWIR high gain bands are satured
    // if they are, we skip bands 3 and 6
    bool skipHighGain = areHighGainBandsSaturated(currPix);

    //this var will track the band for polcor, which will NOT increment when dropping
    // the high or low gain bands. bandCount will increment because it accounts for
    // all SWIR bands and it needs to update to grab the correct m12 and m13 coefs 
    int polcorBandCount = bandCount;
    
    // Process SWIR Bands
    for (int currBand = 0; currBand < numSwirBands; currBand++) {

        // if skipping high gain, skip bands 3 and 6
        if (skipHighGain && (currBand == 3 || currBand == 6)) {
            // no need to set polcor[pix] because polcor takes into account that
            // 2 of the SWIRs will be dropped for the current file
            bandCount++; // need to update to grab correct m12 and m13 coef 
            continue;
        }

        // if not skipping high gain, throw out the low gain ones
        if (!skipHighGain && (currBand == 2 || currBand == 4)) {
            bandCount++; // need to update to grab correct m12 and m13 coef
            continue;
        }

        // band count because it is used to save into polcor arr
        int relativeIndex = currPix * numBands + polcorBandCount;
        vector<float> m12Coef = m12[currHamSide][bandCount];
        vector<float> m13Coef = m13[currHamSide][bandCount];
        calculateAndSetPolCor(l1rec, currPix, relativeIndex, m12Coef, m13Coef);

        // increment bands
        polcorBandCount++;
        bandCount++;
        
    }
   
}
