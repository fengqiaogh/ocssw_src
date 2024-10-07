
#ifndef L1C_MSI_PRIVATE_H
#define L1C_MSI_PRIVATE_H

#include <boost/geometry.hpp>                                                                               
#include <boost/geometry/geometries/point_xy.hpp>                                                           
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/linestring.hpp>
#include <openjpeg.h>

#include <proj.h>

#define maxBands 13 
#define numDetectors 12

typedef boost::geometry::model::d2::point_xy<double> Point_t;
typedef boost::geometry::model::polygon<Point_t> Polygon_t;

typedef struct msi_struct {
    double scene_start_time;
    char* orbit_direction;

    float esdist; // Earth-Sun distance correction factor
    float alt;

    char tileID[FILENAME_MAX]; // to hold the Sentinel TileID
    char imgDir[FILENAME_MAX]; // char array for directories
    char* granuleMetadataFile; //MTD_TL.xml for the granule - contains view angles
    char* datastripMetadataFile; // MTD_DS.xml file for the 'datastrip' - contains orbit info
    char* jp2[maxBands]; // char array for jp2 file name
    char* detectorFootprintFiles[maxBands]; // for simplicity, choose just one 20m band

    // members for projection
    int32_t CSCode; // EPSG Coordinate System Code
    char* UTMZone; // UTM Zone Number
    int32_t *ULCoord; // Upper left coordinates
    PJ *pj; // projection

    // orbit posvec
    int num_gps;
    double *position[3];
    double *gpstime;

    // detector footprint polygons
    // there are 12 detectors arrays - not all are filled for every "scene"
    Polygon_t detectorFootprints[numDetectors];
    // This delta time is the relative difference in time between the detectors
    double detectorDeltaTime[numDetectors][maxBands];
    // The time between each line (ne scan) use to propagate start time to scantime
    double lineTimeDelta;

    // tie-point view angles by band
    double **sensorZenith;
    double **sensorAzimuth;

    // members for decoding jp2 images
    opj_image_t* image[maxBands]; // opj_image_t pointer array for image for each band
    opj_stream_t* l_stream[maxBands]; // stream array
    opj_codec_t* l_codec[maxBands]; // array of handle to a decompressor
    opj_codestream_index_t* cstr_index[maxBands];
    opj_dparameters_t parameters[maxBands]; // decompression parameters
    opj_codestream_info_v2_t* cstr_info[maxBands]; //tile info

    uint32_t *buf; //buffer for data manipulations...meh.
} msi_t;

msi_t* createPrivateData_msi(int numBands);

#endif /* L1C_MSI_PRIVATE_H */

