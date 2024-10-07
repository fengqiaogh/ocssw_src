/*****************************************************************************
 *
 *  NAME:  VcstTerrainCorrect
 *
 *  DESCRIPTION: This is the header file for the VcstTerrainCorrect class for
 *               Terrain Correction.
 *
 *****************************************************************************/

#ifndef VCSTTERRAINCORRECT_h
#define VCSTTERRAINCORRECT_h

#include <VcstPolarTileDatabase.h>
#include <VcstPolarStereographicDataSet.h>
#include <VcstGeoAncInputItem.h>
#include <VcstMath.h>
//#include <Typedefs.h>
//#include <tereco.h>
#include <map>
#include <sstream>

class VcstCmnGeo;


const int MIN_LOS_PTS = 3;
const int MAX_LOS_PTS = 70;
const double MAX_CRCT = 10.e+0; // meters
const double MIN_EARTH_CTR_ANGLE = 1.0e-8;
const double NUM_BI_INTERP_PTS = 4.e+0;
const double BOTTOM_LOS = -600.0e+0; // meters
const double KM2METERS = 1.0e+3;
const double HORIZ_LIMIT = 20.e+0; // meters
const double FILL_AZM = -10.e+0;
const double PAD_VAL = 1.03e+0;

const int NORTH_KEY = 1000;
const int SOUTH_KEY = 3000;
const int MAX_TILE_ID = 5000;

const int TERRCORR_NOCRCT = 20;
const int TERRCORR_OCEAN_CRCT = 25;
const int TERRCORR_ERR = 15;

/****************************************************************************
number of rows and columns in one terrain DB box
 ***************************************/

const int TEBOX_ROWS = 1664;
const int TEBOX_COLS = 1664;

//***************************************************************************
// number of DB box rows and DB box columns in one hemisphere.
//*************************************

const int TEDB_ROWS_OF_BOXES = 32;
const int TEDB_COLS_OF_BOXES = 32;

//
// length of string containing box file name
const int BOX_FILE_STR = 32;

/******************************************************************************
NOTE: "No data by design" points are points that are more than 2 KM outside the
Equator.  They never get valid data in them.  Code which finds a "NOFILL"
flag as the data value has something wrong with it.

  TE_USGS_OCEAN - In GTOPO30 Terrain Data, ocean or ocean connected sea
  TE_BATHY_LAND - A land point in the bathymetric depth data
  TE_NOFILL_BIG - 16 bit number for point with no data by design
  TE_NOFILL_SML - 8 bit number for point with no data by design
    TE_INIT_BIG - point that should be filled by some other value
    TE_INIT_SML - point that should be filled by some other value
        TE_ZERO - a constant used to set variables to zero
 ******************************************************************************/

const int TE_USGS_OCEAN = -9999;
const int TE_BATHY_LAND = 1;
const int TE_NOFILL_BIG = -20000;
const int TE_NOFILL_SML = 120;
const int TE_INIT_BIG = -19999;
const int TE_INIT_SML = 119;
const int TE_ZERO = 0;

/****************************************************************************
These are the contents of one DB box
       mslhgt - Mean Sea Level surface height above NIMA GEOID, in meters
                TE_USGS_OCEAN, -9999, is deep ocean
                TE_NOFILL_BIG, -20000, is a point with no data
       maxhgt - Maximum terrain height within radius of MH_RADIUS kilometers
       minhgt - Minimum terrain height within radius of MH_RADIUS kilometers
        bathy - ETOPO2 Bathymetric depth of water, in meters
                TE_BATHY_LAND, +1, is a land or glacier point
                TE_NOFILL_BIG, -20000, is a point with no data
        sfcgp - Surface geopotential height
        egSep - Separation, NIMA 95 (Geoid - Ellipsoid) in meters
                TE_NOFILL_SML, +120, is a point with no data
          qst - QST-LWM (Quarterly Surface Type-Land Water Mask) data, 0 thru 20
 upIET_mslhgt - time SRTM30 height in the DB box was updated
  upIET_mmhgt - time max height within radius was determined
  upIET_bathy - time ETOPO2 depth in the DB box was updated
  upIET_sfcgp - time surface geopotential height was updated
    upIET_egs - time EG separation in the DB box was updated
    upIET_qst - time QST-LWM data updated.
 upIET_resolv - time that conflict resolution was done on this box
  updated_IET - time, in IET microseconds, that the DB box was Stored
       boxnum - DB box number
       boxrow - the row of boxes that this box is in
       boxcol - the column of boxes that this box is in
 ulrow_offset - DB grid row of the upper left corner pixel in the box
 ulcol_offset - DB grid column of the upper left corner pixel in the box
   num_nofill - number of points well outside Equator, never filled by data
         nOrS - character "N" for Northern hemisphere box, character "S" for
                 Southern hemisphere box
 ***************************************/

typedef struct {
    short mslhgt[TEBOX_ROWS][TEBOX_COLS];
    short maxhgt[TEBOX_ROWS][TEBOX_COLS];
    short minhgt[TEBOX_ROWS][TEBOX_COLS];
    short bathy[TEBOX_ROWS][TEBOX_COLS];
    short sfcgp[TEBOX_ROWS][TEBOX_COLS];
    char egSep[TEBOX_ROWS][TEBOX_COLS];
    unsigned char qst[TEBOX_ROWS][TEBOX_COLS];
    long long upIET_mslhgt;
    long long upIET_mmhgt;
    long long upIET_bathy;
    long long upIET_sfcgp;
    long long upIET_egs;
    long long upIET_qst;
    long long upIET_resolv;
    long long updated_IET;
    char boxfile[BOX_FILE_STR];
    int boxnum;
    int boxrow;
    int boxcol;
    int ulrow_offset;
    int ulcol_offset;
    int num_nofill;
    char nOrS;
} TerecoTileType;



int target_pt(double clat,
        double clon,
        double side_b,
        double azimuth,
        double *alat,
        double *alon);

//-------------------------------------------------------------------

typedef struct LOSPointType {
    double lat; // Latitude
    double lon; // Longitude
    double gRow; // Grid row
    double gCol; // Grid col
    double losHgt; // Line-Of-Sight(LOS) height
    double mslHgt; // Mean Sea Lvl Height
    double egsHgt; // Ellipsoid/Geoid separation height
    double eesHgt; // Ellipsoid to Earth Surface height
    double maxHgt; // Maximum height 32km around this point
    double minHgt; // Minimum height 32km around this point
    Hemisphere hem; // Hemisphere for this point
} LOS_POINT_TYPE;

typedef struct TerrainCorrectDataType {
    LOSPointType interpPts[2]; // Array w/ 2 points to interpolate
    double ellipLat; // Input ellipsoid latitude
    double ellipLon; // Input ellipsoid longitude
    double ellipSatAzm; // Input ellipsoid satellite azimuth
    double ellipSatZen; // Input ellipsoid satellite zenith
    double ellipEarthRad; // Geodetic Earth radius at ellipsoid
    double startDistAng; // Starting distance angle
    double endDistAng; // Ending distance angle
    double pntDistAng; // Distance angle between points
    int lowIdx; // Index of LOS point below terrain
    int upIdx; // Index of LOS point above terrain
    int numPts; // Number of LOS points calculated
} TERRAIN_CORRECT_DATA_TYPE;

typedef struct TerrainCorrectTerecoType {
    VcstGeoAncInputItem* ptr;
    bool tileErr;
} TERRAIN_CORRECT_TERECO_TYPE;

// class ProCmnDataItemEnv;

//------------------------------------------------------------------

class VcstTerrainCorrect {
public:

    /**
     * Constructor
     *
     */
    VcstTerrainCorrect();


    /**
     * Destructor
     */
    virtual ~VcstTerrainCorrect();


    /**
     * This is the main function for terrain correction.  This function
     * will return terrain corrected latitude/longitude and satellite
     * zenith angle.
     *
     * @param inLat          Input ellipsoid intersection latitude
     * @param inLon          Input ellipsoid intersection longitude
     * @param inSatAzm       Input satellite azimuth
     * @param inSatZen       Input satellite zenith
     * @param inRPos         Input ECR position
     * @param cmnGeoPtr      vis to a cmn geo instance
     * @param outEgs         Output ellipsoid-geoid separation
     * @param outLat         Output terrain corrected latitude
     * @param outLon         Output terrain corrected longitude
     * @param outHgt         Output height above the surface
     * @param outSatAzm      Output terrain corrected satellite azimuth
     * @param outSatZen      Output terrain corrected satellite zenith
     * @param outRange       Output terrain corrected satellite range
     *
     * @return               PRO_SUCCESS or PRO_FAIL
     */
    int terrainCorrection(
            const float inLat,
            const float inLon,
            const float inSatAzm,
            const float inSatZen,
            const double inRPos[VEC_SIZE],
            VcstCmnGeo *cmnGeoPtr,
            float& outEgs,
            float& outLat,
            float& outLon,
            float& outHgt,
            float& outSatAzm,
            float& outSatZen,
            float& outRange);


    /**
     * This function will retrieve the MSL height for the 
     * ellipsoid intersection specified
     *
     * @param inLat          Input latitude in radians
     * @param inLon          Input longitude in radians
     * @param cmnGeoPtr      visibility to a cmn geo instance
     * @param outPoint       Output structure containing terrain info
     *
     * @return PRO_SUCCESS or PRO_FAIL
     */
    int getTerrainInfo(const float inLat,
            const float inLon,
            VcstCmnGeo *cmnGeoPtr,
            LOSPointType& outPoint);


private:

    /**
     * Array containing ptrs to DEM data items in DMS
     */
    TERRAIN_CORRECT_TERECO_TYPE tileArr_[MAX_TILE_ID];

    /**
     * Pointer to an instance of ProCmnPolarTileDatabase
     */
    VcstPolarTileDatabase* tileDB_;

    /**
     * Array containing pointer to a DEM tile structures from DMS
     * No ownership.
     */
    TerecoTileType* tileTypeArr_[MAX_TILE_ID];


    /**
     * Data structure used for terrain correction
     */
    TerrainCorrectDataType tcData_;

    /**
     * Short name of the tile in DMS
     */
    std::ostringstream tileId_;

    /**
     * State holders for the DEBUG message throttle.  Protect against
     * flooding the logs with DEBUG messages during a maneuver.
     */
    int dbgMsgThrottle1_;
    int dbgMsgThrottle2_;
    int dbgMsgThrottle3_;

    /**
     * This function will create the Line-Of-Sight(LOS) points between
     * the satellite and the ellipsoid intersection point.
     *
     * @return     PRO_SUCCESS or PRO_FAIL
     */
    int setupLOSPoints(const LOSPointType startPt);

    /**
     * This function will calculate new lat/lon values for a LOS point
     * and then use the lat/lon values to get terrain DB grid coords.
     *
     * @param inAzm         Input azimuth angle
     * @param ioDistAng     Input/Output distance angle
     * @param point         Input/Output terrain LOS point calculating grid
     *                      coordinates for
     *
     * @return              PRO_SUCCESS or PRO_FAIL
     */
    int calcGridFromLatLon(double inAzm,
            double& ioDistAng,
            LOSPointType& point);

    /**
     * This function will retrieve a DEM tile from DMS if necessary
     * and fill an LOS point with the MSL, EEG, and EES height information
     *
     * @param point      Output point to place height information in
     *
     * return            PRO_SUCCESS or error code
     */
    int fillTerrainPoint(LOSPointType& point);

    /**
     * This function will retrieve height information of an LOS point
     * using bi-linear interpolation of the 4 points around it.
     *
     * @param point      Output point to place height information in
     * @param env        vis to my data environment
     *
     * return            PRO_SUCCESS or error code
     */
    int fillTerrainPointBiInterp(LOSPointType& point);


    /**
     * This function will calculate a terrain corrected latitude and
     * longitude and the MSL surface height.
     *
     * @param inRPos      Input ECR position
     * @param cmnGeoPtr     vis to a pro sdr cmn geo instance
     * @param outLat      Output terrain corrected latitude
     * @param outLon      Output terrain corrected longitude
     * @param outMslHgt   Output MSL surface height
     * @param outAzm      Output new satellite azimuth
     * @param outZen      Output new satellite zenith
     *
     * @return            PRO_SUCCESS or PRO_FAIL
     */
    int findTerrain(
            const double inRPos[VEC_SIZE],
            VcstCmnGeo *cmnGeoPtr,
            double& outLat,
            double& outLon,
            double& outMslHgt,
            double& outAzm,
            double& outZen);

    /**
     * This function will determine where the Line-Of-Sight(LOS)
     * intersects the terrain when the terrain is below the end point
     *
     * @param inStartPt        Input starting point.  At or below
     *                         ellipsoid intersection.
     * @param env              data environment we're operating under.
     *
     * @return                 PRO_SUCCESS or TERRCORR_ERR
     */
    int determineIntersectBelowEndPt(const LOSPointType& inStartPt);

    /**
     * This function will determine where the Line-Of-Sight(LOS)
     * intersects the terrain when the terrain is above the end point
     *
     * @param env              data environment we're operating under.
     * @return                 PRO_SUCCESS or TERRCORR_ERR
     */
    int determineIntersectAboveEndPt();

    /**
     * This function will calculate the latitude/longitude of a point
     * and then fill the point will data from the terrain DB.
     *
     * @param idx             Index of the LOS point to create
     * @param env             visibility to my data envirornment
     * @param outPoint        Output point that was created and filled
     *
     * @return                PRO_SUCCESS or TERRCORR_ERR
     */
    int calcAndFillPoint(const int idx,
            LOSPointType& outPoint);

    /**
     * This function will calculate the starting point grid coords
     * and fill the point with data from the terrain DB.
     *
     * @param outPoint        Output point to fill
     * @param env             data environment we're operating under
     *
     * @return                PRO_SUCCESS or TERRCORR_ERR
     */
    int calcStartPoint(LOSPointType& outPoint);

    /**
     * This function will calculate the LOS height for a point.'
     *
     * @param idx             Index of the point to calculate height for
     * @param outPoint        Output point to calculate LOS height for
     *
     * @return                Void
     */
    void calcLOSHgt(const int idx,
            LOSPointType& outPoint);


    /**
     * This function will calculate a satellite azimuth/zenith angle
     * based on a corrected lat/lon value.
     *
     * @param inLat         Input corrected latitude
     * @param inLon         Input corrected longitude
     * @param inMsl         Input MSL height of the pixel
     * @param inEgs         Input EGS height of the pixel
     * @param inRPos        Input ECR position of the satellite
     * @param cmnGeoPtr     vis to a pro sdr cmn geo instance
     * @param outAzm        Output satellite azimuth
     * @param outZen        Output satellite zenith
     *
     * @return Void
     */
    void calcTCSatAzmZen(
            const double inLat,
            const double inLon,
            const double inMsl,
            const double inEgs,
            const double inRPos[VEC_SIZE],
            VcstCmnGeo *cmnGeoPtr,
            double& outAzm,
            double& outZen);


    /**
     * This function will apply a correction to ocean pixels
     *
     * @param inPt          Input LOS point for the ellipsoid pixel
     * @param inRPos        Input ECR position of the satellite
     * @param cmnGeoPtr     vis to a pro cmn geo instance
     * @param outLat        Output corrected latitude
     * @param outLon        Output corrected longitude
     * @param outHt         Output MSL height
     * @param outAzm        Output corrected satellite azimuth
     * @param outZen        Output corrected satellite zenith
     *
     * @return Void
     */
    void applyOceanCorrection(
            const LOSPointType& inPt,
            const double inRPos[VEC_SIZE],
            VcstCmnGeo *cmnGeoPtr,
            float& outLat,
            float& outLon,
            float& outHt,
            float& outAzm,
            float& outZen);


    /**
     * This function will convert lat/lon to polarstereographic
     * grid coordinates
     *
     * @param point         LOS Point to get grid coords. for
     *
     * @return PRO_SUCCESS or PRO_FAIL
     */
    int convertLatLonToGrid(LOSPointType& point);

private:

};


#endif
