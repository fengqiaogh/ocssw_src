/* Norman Kuring        22-Oct-1997 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h>
#include <time.h>
#include <math.h>

#include <GetStationInfo.h>

#include "swl0_proto.h"
#include <dfutils.h>
#include <timeutils.h>

/*
These global variables facilitate communication
between some of the functions defined in this file.
 */
static idDS ds_id;
static int32 numScans, numPixels;
static int16 tdi_global[8];
static char *hdfFile;
static int firstCallThisFile;
static char dataTypeString[5];
static int16 startYear, startDay, endDay;
static int32 startMillisec;
static int calibrationAppended = 0;

#define SENSOR          "Sea-viewing Wide Field-of-view Sensor (SeaWiFS)"
#define MISSIONCHAR     "Nominal orbit: inclination = 98.2 (Sun-synchronous); node = 12 noon local (descending); eccentricity = <0.002; altitude = 705 km; ground speed = 6.75 km/sec"
#define SENSORCHAR      "Number of bands = 8; number of active bands = 8; wavelengths per band (nm) = 412, 443, 490, 510, 555, 670, 765, 865; bits per pixel = 10; instantaneous field-of-view = 1.5835 mrad; pixels per scan = 1285; scan rate = 6/sec; sample rate = 7710/sec"

/****************************************************************************
Create a new HDF file and store some global attributes in it.
This function must be called before WriteScanData().
CloseL1aFile() should be called to finish up the file.
 *****************************************************************************/
int CreateL1aFile(
        char *path,
        swl0scene *scene,
        char *proccon,
        char *proclog,
        swl0ctl *l0ctl
        ) {

    unsigned char dataType;
    StationInfo stationInfo;
    int16 year;
    int32 millisec;
    int32 startpix, subsamp;

    /* for seadas version and meta */

    if (scene->type == HRPT)
        dataType = 16;
    else
        dataType = scene->mnftype;

    /* Get some ground station specific information. */
    PTB(GetStationInfo(l0ctl->stationInfoFile, &stationInfo));

    /* A few numbers are determined by the data type. */
    if (dataType == GACTYPE) {
        numPixels = NPIXGAC;
        startpix = SPIXGAC;
        subsamp = IPIXGAC;
    } else {
        numPixels = NPIXLAC;
        startpix = SPIXLAC;
        subsamp = IPIXLAC;
    }

    numScans = scene->nscan; /* Make it global. */

    /*
    Copy the output filename to a static global variable.
     */
    MALLOC(hdfFile, char, strlen(path) + 1);
    strcpy(hdfFile, path);

    /*
    Make sure the file does not already exist.
    if(fopen(hdfFile, "rb") != NULL){
      fprintf(stderr,
      "-E- %s line %d: File, %s , already exists.\n",
      __FILE__,__LINE__,hdfFile);
      return(FILE_ALREADY_EXISTS);
    }
     */

    /*
    Create the HDF file.
     */
    ds_id = startDS(hdfFile, DS_HDF, DS_WRITE, 0);
    if (ds_id.fid == FAIL) {
        fprintf(stderr,
                "-E- %s line %d: Could not create HDF file, %s .\n",
                __FILE__, __LINE__, hdfFile);
        return (HDF_FUNCTION_ERROR);
    }

    /*
    Set global variables that are also used
    in getting calibration data for the scene.
     */
    strncpy(dataTypeString, DTypeString(dataType), 4);
    DecomposeTime(scene->stime, &startYear, &startDay, &startMillisec);
    DecomposeTime(scene->etime, &year, &endDay, &millisec);
    firstCallThisFile = 0;

    /*
    Write out some global attributes.
     */
    PTB(SetChrGA(ds_id, "Product Name", basename(hdfFile)));
    PTB(SetChrGA(ds_id, "Title", "SeaWiFS Level-1A Data"));
    PTB(SetChrGA(ds_id, "Data Center", stationInfo.data_center));
    PTB(SetChrGA(ds_id, "Station Name", stationInfo.station_name));
    PTB(SetF32GA(ds_id, "Station Latitude", stationInfo.station_latitude));
    PTB(SetF32GA(ds_id, "Station Longitude", stationInfo.station_longitude));
    PTB(SetChrGA(ds_id, "Mission", "SeaStar SeaWiFS"));
    PTB(SetChrGA(ds_id, "Mission Characteristics", MISSIONCHAR));
    PTB(SetChrGA(ds_id, "Sensor", SENSOR));
    PTB(SetChrGA(ds_id, "Sensor Characteristics", SENSORCHAR));
    PTB(SetChrGA(ds_id, "Data Type", DTypeString(dataType)));
    PTB(SetChrGA(ds_id, "Replacement Flag", "ORIGINAL"));
    PTB(SetChrGA(ds_id, "Software ID", L01VERSION));
    PTB(SetChrGA(ds_id, "Processing Time", ydhmsf(now(), 'L')));
    PTB(SetChrGA(ds_id, "Input Files", scene->l0file));
    PTB(SetChrGA(ds_id, "Processing Control", proccon));
    PTB(SetChrGA(ds_id, "Processing Log", proclog));
    PTB(SetChrGA(ds_id, "Start Time", ydhmsf(scene->stime, 'G')));
    PTB(SetChrGA(ds_id, "End Time", ydhmsf(scene->etime, 'G')));
    PTB(SetChrGA(ds_id, "Scene Center Time", ydhmsf(scene->ctime, 'G')));
    PTB(SetChrGA(ds_id, "Node Crossing Time", ydhmsf(scene->node_time, 'G')));
    PTB(SetI16GA(ds_id, "Start Year", startYear));
    PTB(SetI16GA(ds_id, "Start Day", startDay));
    PTB(SetI32GA(ds_id, "Start Millisec", startMillisec));
    PTB(SetI16GA(ds_id, "End Year", year));
    PTB(SetI16GA(ds_id, "End Day", endDay));
    PTB(SetI32GA(ds_id, "End Millisec", millisec));
    PTB(SetI32GA(ds_id, "Orbit Number", scene->orbnum));
    PTB(SetChrGA(ds_id, "Start Node", scene->start_node));
    PTB(SetChrGA(ds_id, "End Node", scene->end_node));
    PTB(SetChrGA(ds_id, "NORAD Line 1", ""));
    PTB(SetChrGA(ds_id, "NORAD Line 2", ""));
    PTB(SetI32GA(ds_id, "Pixels per Scan Line", numPixels));
    PTB(SetI32GA(ds_id, "Number of Scan Lines", numScans));
    PTB(SetI32GA(ds_id, "LAC Pixel Start Number", startpix));
    PTB(SetI32GA(ds_id, "LAC Pixel Subsampling", subsamp));
    PTB(SetI32GA(ds_id, "Scene Center Scan Line", scene->center_scan_line));
    PTB(SetI32GA(ds_id, "Filled Scan Lines", 0));
    PTB(SetI32GA(ds_id, "FF Missing Frames", 0));
    PTB(SetI32GA(ds_id, "SDPS Missing Frames", 0));
    PTB(SetChrGA(ds_id, "Latitude Units", "degrees North"));
    PTB(SetChrGA(ds_id, "Longitude Units", "degrees East"));
    PTB(SetF32GA(ds_id, "Scene Center Latitude", scene->center_lat));
    PTB(SetF32GA(ds_id, "Scene Center Longitude", scene->center_lon));
    PTB(SetF32GA(ds_id, "Scene Center Solar Zenith", scene->center_solz));
    PTB(SetF32GA(ds_id, "Upper Left Latitude", scene->upper_left_lat));
    PTB(SetF32GA(ds_id, "Upper Left Longitude", scene->upper_left_lon));
    PTB(SetF32GA(ds_id, "Upper Right Latitude", scene->upper_right_lat));
    PTB(SetF32GA(ds_id, "Upper Right Longitude", scene->upper_right_lon));
    PTB(SetF32GA(ds_id, "Lower Left Latitude", scene->lower_left_lat));
    PTB(SetF32GA(ds_id, "Lower Left Longitude", scene->lower_left_lon));
    PTB(SetF32GA(ds_id, "Lower Right Latitude", scene->lower_right_lat));
    PTB(SetF32GA(ds_id, "Lower Right Longitude", scene->lower_right_lon));
    PTB(SetF32GA(ds_id, "Northernmost Latitude", scene->northern_lat));
    PTB(SetF32GA(ds_id, "Southernmost Latitude", scene->southern_lat));
    PTB(SetF32GA(ds_id, "Westernmost Longitude", scene->western_lon));
    PTB(SetF32GA(ds_id, "Easternmost Longitude", scene->eastern_lon));
    PTB(SetF32GA(ds_id, "Start Center Latitude", scene->start_center_lat));
    PTB(SetF32GA(ds_id, "Start Center Longitude", scene->start_center_lon));
    PTB(SetF32GA(ds_id, "End Center Latitude", scene->end_center_lat));
    PTB(SetF32GA(ds_id, "End Center Longitude", scene->end_center_lon));
    PTB(SetF32GA(ds_id, "Orbit Node Longitude", scene->node_lon));

    /* Add the tilt data. */
    PTB(AddTiltData(scene->ntilts,
            scene->tilt_flags,
            scene->tilt_ranges,
            scene->tilt_lats,
            scene->tilt_lons));

    /* Create the scan-line SDSes */
    PTB(CreateScanData(numScans, numPixels));

    return (LIFE_IS_GOOD);
}

/****************************************************************************
Create the level-1A SDSes that are to contain the scan-line data in the
currently open HDF file.
 *****************************************************************************/
int CreateScanData(int32 ns, int32 np) {

    PTB(CreateSDS(ds_id.fid,
            "msec", /* short name */
            "Scan-line time, milliseconds of day", /* long name */
            NULL, /* standard name */
            NULL, /* units */
            0, 86399999, /* valid range */
            1.0, 0.0, /* slope,offset */
            DFNT_INT32, /* HDF number type */
            1, /* rank */
            ns, 1, 1, /* dimension sizes */
            "Number of Scan Lines", NULL, NULL /* dimension names */
            ));

    PTB(CreateSDS(ds_id.fid,
            "eng_qual", /* short name */
            "Engineering data-out-of-range flags", /* long name */
            NULL, /* standard name */
            NULL, /* units */
            0, 0, /* valid range */
            1.0, 0.0, /* slope,offset */
            DFNT_UINT8, /* HDF number type */
            2, /* rank */
            ns, 4, 1, /* dimension sizes */
            "Number of Scan Lines", "4", NULL /* dimension names */
            ));

    PTB(CreateSDS(ds_id.fid,
            "s_flags", /* short name */
            "Scan-line quality flags", /* long name */
            NULL, /* standard name */
            NULL, /* units */
            0, 0, /* valid range */
            1.0, 0.0, /* slope,offset */
            DFNT_UINT8, /* HDF number type */
            2, /* rank */
            ns, 4, 1, /* dimension sizes */
            "Number of Scan Lines", "4", NULL /* dimension names */
            ));

    PTB(CreateSDS(ds_id.fid,
            "s_satp", /* short name */
            "Number of saturated pixels per band", /* long name */
            NULL, /* standard name */
            NULL, /* units */
            0, 0, /* valid range */
            1.0, 0.0, /* slope,offset */
            DFNT_INT16, /* HDF number type */
            2, /* rank */
            ns, 8, 1, /* dimension sizes */
            "Number of Scan Lines", "Number of Bands", NULL /* dimension names */
            ));

    PTB(CreateSDS(ds_id.fid,
            "s_zerop", /* short name */
            "Number of zero pixels per band", /* long name */
            NULL, /* standard name */
            NULL, /* units */
            0, 0, /* valid range */
            1.0, 0.0, /* slope,offset */
            DFNT_INT16, /* HDF number type */
            2, /* rank */
            ns, 8, 1, /* dimension sizes */
            "Number of Scan Lines", "Number of Bands", NULL /* dimension names */
            ));

    PTB(CreateSDS(ds_id.fid,
            "slat", /* short name */
            "Scan start-pixel latitude", /* long name */
            NULL, /* standard name */
            NULL, /* units */
            -90.0, 90.0, /* valid range */
            1.0, 0.0, /* slope,offset */
            DFNT_FLOAT32, /* HDF number type */
            1, /* rank */
            ns, 1, 1, /* dimension sizes */
            "Number of Scan Lines", NULL, NULL /* dimension names */
            ));

    PTB(CreateSDS(ds_id.fid,
            "slon", /* short name */
            "Scan start-pixel longitude", /* long name */
            NULL, /* standard name */
            NULL, /* units */
            -180.0, 180.0, /* valid range */
            1.0, 0.0, /* slope,offset */
            DFNT_FLOAT32, /* HDF number type */
            1, /* rank */
            ns, 1, 1, /* dimension sizes */
            "Number of Scan Lines", NULL, NULL /* dimension names */
            ));

    PTB(CreateSDS(ds_id.fid,
            "clat", /* short name */
            "Scan center-pixel latitude", /* long name */
            NULL, /* standard name */
            NULL, /* units */
            -90.0, 90.0, /* valid range */
            1.0, 0.0, /* slope,offset */
            DFNT_FLOAT32, /* HDF number type */
            1, /* rank */
            ns, 1, 1, /* dimension sizes */
            "Number of Scan Lines", NULL, NULL /* dimension names */
            ));

    PTB(CreateSDS(ds_id.fid,
            "clon", /* short name */
            "Scan center-pixel longitude", /* long name */
            NULL, /* standard name */
            NULL, /* units */
            -180.0, 180.0, /* valid range */
            1.0, 0.0, /* slope,offset */
            DFNT_FLOAT32, /* HDF number type */
            1, /* rank */
            ns, 1, 1, /* dimension sizes */
            "Number of Scan Lines", NULL, NULL /* dimension names */
            ));

    PTB(CreateSDS(ds_id.fid,
            "elat", /* short name */
            "Scan end-pixel latitude", /* long name */
            NULL, /* standard name */
            NULL, /* units */
            -90.0, 90.0, /* valid range */
            1.0, 0.0, /* slope,offset */
            DFNT_FLOAT32, /* HDF number type */
            1, /* rank */
            ns, 1, 1, /* dimension sizes */
            "Number of Scan Lines", NULL, NULL /* dimension names */
            ));

    PTB(CreateSDS(ds_id.fid,
            "elon", /* short name */
            "Scan end-pixel longitude", /* long name */
            NULL, /* standard name */
            NULL, /* units */
            -180.0, 180.0, /* valid range */
            1.0, 0.0, /* slope,offset */
            DFNT_FLOAT32, /* HDF number type */
            1, /* rank */
            ns, 1, 1, /* dimension sizes */
            "Number of Scan Lines", NULL, NULL /* dimension names */
            ));

    PTB(CreateSDS(ds_id.fid,
            "csol_z", /* short name */
            "Scan center-pixel solar zenith angle", /* long name */
            NULL, /* standard name */
            NULL, /* units */
            0.0, 180.0, /* valid range */
            1.0, 0.0, /* slope,offset */
            DFNT_FLOAT32, /* HDF number type */
            1, /* rank */
            ns, 1, 1, /* dimension sizes */
            "Number of Scan Lines", NULL, NULL /* dimension names */
            ));

    PTB(CreateSDS(ds_id.fid,
            "tilt", /* short name */
            "Tilt angle for scan line", /* long name */
            NULL, /* standard name */
            NULL, /* units */
            -20.1, 20.1, /* valid range */
            1.0, 0.0, /* slope,offset */
            DFNT_FLOAT32, /* HDF number type */
            1, /* rank */
            ns, 1, 1, /* dimension sizes */
            "Number of Scan Lines", NULL, NULL /* dimension names */
            ));

    PTB(CreateSDS(ds_id.fid,
            "sc_id", /* short name */
            "Spacecraft ID", /* long name */
            NULL, /* standard name */
            NULL, /* units */
            0, 0, /* valid range */
            1.0, 0.0, /* slope,offset */
            DFNT_INT16, /* HDF number type */
            2, /* rank */
            ns, 2, 1, /* dimension sizes */
            "Number of Scan Lines", "2", NULL /* dimension names */
            ));

    PTB(CreateSDS(ds_id.fid,
            "sc_ttag", /* short name */
            "Spacecraft time tag", /* long name */
            NULL, /* standard name */
            NULL, /* units */
            0, 0, /* valid range */
            1.0, 0.0, /* slope,offset */
            DFNT_INT16, /* HDF number type */
            2, /* rank */
            ns, 4, 1, /* dimension sizes */
            "Number of Scan Lines", "4", NULL /* dimension names */
            ));

    PTB(CreateSDS(ds_id.fid,
            "sc_soh", /* short name */
            "Spacecraft state-of-health data", /* long name */
            NULL, /* standard name */
            NULL, /* units */
            0, 0, /* valid range */
            1.0, 0.0, /* slope,offset */
            DFNT_UINT8, /* HDF number type */
            2, /* rank */
            ns, 775, 1, /* dimension sizes */
            "Number of Scan Lines", "775", NULL /* dimension names */
            ));

    PTB(CreateSDS(ds_id.fid,
            "inst_tlm", /* short name */
            "SeaWiFS instrument telemetry", /* long name */
            NULL, /* standard name */
            NULL, /* units */
            0, 0, /* valid range */
            1.0, 0.0, /* slope,offset */
            DFNT_INT16, /* HDF number type */
            2, /* rank */
            ns, 44, 1, /* dimension sizes */
            "Number of Scan Lines", "44", NULL /* dimension names */
            ));

    PTB(CreateSDS(ds_id.fid,
            "l1a_data", /* short name */
            "Level-1A data", /* long name */
            NULL, /* standard name */
            "radiance counts", /* units */
            0, 1023, /* valid range */
            1.0, 0.0, /* slope,offset */
            DFNT_INT16, /* HDF number type */
            3, /* rank */
            ns, np, 8, /* dimension sizes */
            "Number of Scan Lines", "Pixels per Scan Line", "Number of Bands" /* dimnames */
            ));

    PTB(CreateSDS(ds_id.fid,
            "start_syn", /* short name */
            "Start-synch pixel", /* long name */
            NULL, /* standard name */
            NULL, /* units */
            0, 0, /* valid range */
            1.0, 0.0, /* slope,offset */
            DFNT_INT16, /* HDF number type */
            2, /* rank */
            ns, 8, 1, /* dimension sizes */
            "Number of Scan Lines", "Number of Bands", NULL /* dimension names */
            ));

    PTB(CreateSDS(ds_id.fid,
            "stop_syn", /* short name */
            "Stop-synch pixel", /* long name */
            NULL, /* standard name */
            NULL, /* units */
            0, 0, /* valid range */
            1.0, 0.0, /* slope,offset */
            DFNT_INT16, /* HDF number type */
            2, /* rank */
            ns, 8, 1, /* dimension sizes */
            "Number of Scan Lines", "Number of Bands", NULL /* dimension names */
            ));

    PTB(CreateSDS(ds_id.fid,
            "dark_rest", /* short name */
            "Dark-restore pixel", /* long name */
            NULL, /* standard name */
            NULL, /* units */
            0, 0, /* valid range */
            1.0, 0.0, /* slope,offset */
            DFNT_INT16, /* HDF number type */
            2, /* rank */
            ns, 8, 1, /* dimension sizes */
            "Number of Scan Lines", "Number of Bands", NULL /* dimension names */
            ));

    PTB(CreateSDS(ds_id.fid,
            "gain", /* short name */
            "Band gain settings", /* long name */
            NULL, /* standard name */
            NULL, /* units */
            0, 3, /* valid range */
            1.0, 0.0, /* slope,offset */
            DFNT_INT16, /* HDF number type */
            2, /* rank */
            ns, 8, 1, /* dimension sizes */
            "Number of Scan Lines", "Number of Bands", NULL /* dimension names */
            ));

    PTB(CreateSDS(ds_id.fid,
            "tdi", /* short name */
            "Band time-delay and integration settings", /* long name */
            NULL, /* standard name */
            NULL, /* units */
            0, 255, /* valid range */
            1.0, 0.0, /* slope,offset */
            DFNT_INT16, /* HDF number type */
            2, /* rank */
            ns, 8, 1, /* dimension sizes */
            "Number of Scan Lines", "Number of Bands", NULL /* dimension names */
            ));

    PTB(CreateSDS(ds_id.fid,
            "inst_ana", /* short name */
            "Instrument analog telemetry", /* long name */
            NULL, /* standard name */
            NULL, /* units */
            0, 0, /* valid range */
            1.0, 0.0, /* slope,offset */
            DFNT_FLOAT32, /* HDF number type */
            2, /* rank */
            ns, 40, 1, /* dimension sizes */
            "Number of Scan Lines", "40", NULL /* dimension names */
            ));

    PTB(CreateSDS(ds_id.fid,
            "inst_dis", /* short name */
            "Instrument discrete telemetry", /* long name */
            NULL, /* standard name */
            NULL, /* units */
            0, 0, /* valid range */
            1.0, 0.0, /* slope,offset */
            DFNT_UINT8, /* HDF number type */
            2, /* rank */
            ns, 32, 1, /* dimension sizes */
            "Number of Scan Lines", "32", NULL /* dimension names */
            ));

    PTB(CreateSDS(ds_id.fid,
            "sc_ana", /* short name */
            "Spacecraft analog telemetry", /* long name */
            NULL, /* standard name */
            NULL, /* units */
            0, 0, /* valid range */
            1.0, 0.0, /* slope,offset */
            DFNT_FLOAT32, /* HDF number type */
            2, /* rank */
            ns, 40, 1, /* dimension sizes */
            "Number of Scan Lines", "40", NULL /* dimension names */
            ));

    PTB(CreateSDS(ds_id.fid,
            "sc_dis", /* short name */
            "Spacecraft discrete telemetry", /* long name */
            NULL, /* standard name */
            NULL, /* units */
            0, 0, /* valid range */
            1.0, 0.0, /* slope,offset */
            DFNT_UINT8, /* HDF number type */
            2, /* rank */
            ns, 40, 1, /* dimension sizes */
            "Number of Scan Lines", "40", NULL /* dimension names */
            ));

    PTB(CreateSDS(ds_id.fid,
            "scan_temp", /* short name */
            "Detector temperature counts", /* long name */
            NULL, /* standard name */
            NULL, /* units */
            0, 255, /* valid range */
            1.0, 0.0, /* slope,offset */
            DFNT_INT16, /* HDF number type */
            2, /* rank */
            ns, 8, 1, /* dimension sizes */
            "Number of Scan Lines", "Number of Bands", NULL /* dimension names */
            ));

    PTB(CreateSDS(ds_id.fid,
            "side", /* short name */
            "Mirror side for scan line", /* long name */
            NULL, /* standard name */
            NULL, /* units */
            0, 1, /* valid range */
            1.0, 0.0, /* slope,offset */
            DFNT_INT16, /* HDF number type */
            1, /* rank */
            ns, 1, 1, /* dimension sizes */
            "Number of Scan Lines", NULL, NULL /* dimension names */
            ));

    PTB(CreateSDS(ds_id.fid,
            "orb_vec", /* short name */
            "Orbit position vector at scan line time", /* long name */
            NULL, /* standard name */
            "kilometers", /* units */
            -7200.0, 7200.0, /* valid range */
            1.0, 0.0, /* slope,offset */
            DFNT_FLOAT32, /* HDF number type */
            2, /* rank */
            ns, 3, 1, /* dimension sizes */
            "Number of Scan Lines", "3", NULL /* dimension names */
            ));

    PTB(CreateSDS(ds_id.fid,
            "l_vert", /* short name */
            "Local vertical vector in ECEF frame", /* long name */
            NULL, /* standard name */
            NULL, /* units */
            -1.0, 1.0, /* valid range */
            1.0, 0.0, /* slope,offset */
            DFNT_FLOAT32, /* HDF number type */
            2, /* rank */
            ns, 3, 1, /* dimension sizes */
            "Number of Scan Lines", "3", NULL /* dimension names */
            ));

    PTB(CreateSDS(ds_id.fid,
            "sun_ref", /* short name */
            "Reference Sun vector in ECEF frame", /* long name */
            NULL, /* standard name */
            NULL, /* units */
            -1.0, 1.0, /* valid range */
            1.0, 0.0, /* slope,offset */
            DFNT_FLOAT32, /* HDF number type */
            2, /* rank */
            ns, 3, 1, /* dimension sizes */
            "Number of Scan Lines", "3", NULL /* dimension names */
            ));

    PTB(CreateSDS(ds_id.fid,
            "att_ang", /* short name */
            "Computed yaw, roll, pitch", /* long name */
            NULL, /* standard name */
            NULL, /* units */
            -180.0, 180.0, /* valid range */
            1.0, 0.0, /* slope,offset */
            DFNT_FLOAT32, /* HDF number type */
            2, /* rank */
            ns, 3, 1, /* dimension sizes */
            "Number of Scan Lines", "3", NULL /* dimension names */
            ));

    PTB(CreateSDS(ds_id.fid,
            "sen_mat", /* short name */
            "ECEF-to-sensor-frame matrix", /* long name */
            NULL, /* standard name */
            NULL, /* units */
            -1.0, 1.0, /* valid range */
            1.0, 0.0, /* slope,offset */
            DFNT_FLOAT32, /* HDF number type */
            3, /* rank */
            ns, 3, 3, /* dimension sizes */
            "Number of Scan Lines", "3", "3" /* dimension names */
            ));

    PTB(CreateSDS(ds_id.fid,
            "scan_ell", /* short name */
            "Scan-track ellipse coefficients", /* long name */
            NULL, /* standard name */
            NULL, /* units */
            0, 0, /* valid range */
            1.0, 0.0, /* slope,offset */
            DFNT_FLOAT32, /* HDF number type */
            2, /* rank */
            ns, 6, 1, /* dimension sizes */
            "Number of Scan Lines", "6", NULL /* dimension names */
            ));

    PTB(CreateSDS(ds_id.fid,
            "nflag", /* short name */
            "Navigation flags", /* long name */
            NULL, /* standard name */
            NULL, /* units */
            0, 0, /* valid range */
            1.0, 0.0, /* slope,offset */
            DFNT_INT32, /* HDF number type */
            2, /* rank */
            ns, 8, 1, /* dimension sizes */
            "Number of Scan Lines", "8", NULL /* dimension names */
            ));

    return (LIFE_IS_GOOD);
}

/****************************************************************************
Write one scan line's worth of each of the passed in data to the appropriate
SDS in the currently open HDF file.
 *****************************************************************************/
int WriteScanData(
        int32 scan,
        swl1rec *l1rec
        ) {

    int i;
    int32 msec = l1rec->msec;
    uint8 *eng_qual = l1rec->eng_qual;
    uint8 *s_flags = l1rec->s_flags;
    int16 *s_satp = l1rec->s_satp;
    int16 *s_zerop = l1rec->s_zerop;
    float32 slat = l1rec->slat;
    float32 slon = l1rec->slon;
    float32 clat = l1rec->clat;
    float32 clon = l1rec->clon;
    float32 elat = l1rec->elat;
    float32 elon = l1rec->elon;
    float32 csol_z = l1rec->csol_z;
    float32 tilt = l1rec->tilt;
    int16 *sc_id = l1rec->scid;
    int16 *sc_ttag = l1rec->ttag;
    uint8 *sc_soh = l1rec->soh;
    int16 *inst_tlm = l1rec->inst;
    int16 *l1a_data = &l1rec->data[0][0];
    int16 *start_syn = l1rec->startpix;
    int16 *stop_syn = l1rec->stoppix;
    int16 *dark_rest = l1rec->darkpix;
    int16 *gain = l1rec->gain;
    int16 *tdi = l1rec->tdi;
    float32 *inst_ana = l1rec->inst_ana;
    uint8 *inst_dis = l1rec->inst_dis;
    float32 *sc_ana = l1rec->sc_ana;
    uint8 *sc_dis = l1rec->sc_dis;
    int16 *scan_temp = l1rec->scan_temp;
    int16 side = l1rec->side;
    float32 *orb_vec = l1rec->orb_vec;
    float32 *l_vert = l1rec->l_vert;
    float32 *sun_ref = l1rec->sun_ref;
    float32 *att_ang = l1rec->att_ang;
    float32 *sen_mat = &l1rec->sen_mat[0][0];
    float32 *scan_ell = l1rec->scan_ell;
    int32 *nflag = (int32 *) l1rec->nflag;

    if (scan >= numScans) {
        fprintf(stderr, "-W- %s line %d: ", __FILE__, __LINE__);
        fprintf(stderr, "WriteScanData() called with scanline number (%d) ", scan);
        fprintf(stderr, "that is inappropriate to the number of scanlines (%d) ",
                numScans);
        fprintf(stderr, "in the current HDF file, %s .  Call ignored.\n", hdfFile);
        return (PROGRAMMER_BOOBOO);
    }

    PTB(sd_writedata(ds_id.fid, "msec", &msec, scan, 0, 0, 1, 1, 1));
    PTB(sd_writedata(ds_id.fid, "eng_qual", eng_qual, scan, 0, 0, 1, 4, 1));
    PTB(sd_writedata(ds_id.fid, "s_flags", s_flags, scan, 0, 0, 1, 4, 1));
    PTB(sd_writedata(ds_id.fid, "s_satp", s_satp, scan, 0, 0, 1, 8, 1));
    PTB(sd_writedata(ds_id.fid, "s_zerop", s_zerop, scan, 0, 0, 1, 8, 1));
    PTB(sd_writedata(ds_id.fid, "slat", &slat, scan, 0, 0, 1, 1, 1));
    PTB(sd_writedata(ds_id.fid, "slon", &slon, scan, 0, 0, 1, 1, 1));
    PTB(sd_writedata(ds_id.fid, "clat", &clat, scan, 0, 0, 1, 1, 1));
    PTB(sd_writedata(ds_id.fid, "clon", &clon, scan, 0, 0, 1, 1, 1));
    PTB(sd_writedata(ds_id.fid, "elat", &elat, scan, 0, 0, 1, 1, 1));
    PTB(sd_writedata(ds_id.fid, "elon", &elon, scan, 0, 0, 1, 1, 1));
    PTB(sd_writedata(ds_id.fid, "csol_z", &csol_z, scan, 0, 0, 1, 1, 1));
    PTB(sd_writedata(ds_id.fid, "tilt", &tilt, scan, 0, 0, 1, 1, 1));
    PTB(sd_writedata(ds_id.fid, "sc_id", sc_id, scan, 0, 0, 1, 2, 1));
    PTB(sd_writedata(ds_id.fid, "sc_ttag", sc_ttag, scan, 0, 0, 1, 4, 1));
    PTB(sd_writedata(ds_id.fid, "sc_soh", sc_soh, scan, 0, 0, 1, 775, 1));
    PTB(sd_writedata(ds_id.fid, "inst_tlm", inst_tlm, scan, 0, 0, 1, 44, 1));
    PTB(sd_writedata(ds_id.fid, "l1a_data", l1a_data, scan, 0, 0, 1, numPixels, 8));
    PTB(sd_writedata(ds_id.fid, "start_syn", start_syn, scan, 0, 0, 1, 8, 1));
    PTB(sd_writedata(ds_id.fid, "stop_syn", stop_syn, scan, 0, 0, 1, 8, 1));
    PTB(sd_writedata(ds_id.fid, "dark_rest", dark_rest, scan, 0, 0, 1, 8, 1));
    PTB(sd_writedata(ds_id.fid, "gain", gain, scan, 0, 0, 1, 8, 1));
    PTB(sd_writedata(ds_id.fid, "tdi", tdi, scan, 0, 0, 1, 8, 1));
    PTB(sd_writedata(ds_id.fid, "inst_ana", inst_ana, scan, 0, 0, 1, 40, 1));
    PTB(sd_writedata(ds_id.fid, "inst_dis", inst_dis, scan, 0, 0, 1, 32, 1));
    PTB(sd_writedata(ds_id.fid, "sc_ana", sc_ana, scan, 0, 0, 1, 40, 1));
    PTB(sd_writedata(ds_id.fid, "sc_dis", sc_dis, scan, 0, 0, 1, 40, 1));
    PTB(sd_writedata(ds_id.fid, "scan_temp", scan_temp, scan, 0, 0, 1, 8, 1));
    PTB(sd_writedata(ds_id.fid, "side", &side, scan, 0, 0, 1, 1, 1));
    PTB(sd_writedata(ds_id.fid, "orb_vec", orb_vec, scan, 0, 0, 1, 3, 1));
    PTB(sd_writedata(ds_id.fid, "l_vert", l_vert, scan, 0, 0, 1, 3, 1));
    PTB(sd_writedata(ds_id.fid, "sun_ref", sun_ref, scan, 0, 0, 1, 3, 1));
    PTB(sd_writedata(ds_id.fid, "att_ang", att_ang, scan, 0, 0, 1, 3, 1));
    PTB(sd_writedata(ds_id.fid, "sen_mat", sen_mat, scan, 0, 0, 1, 3, 3));
    PTB(sd_writedata(ds_id.fid, "scan_ell", scan_ell, scan, 0, 0, 1, 6, 1));
    PTB(sd_writedata(ds_id.fid, "nflag", nflag, scan, 0, 0, 1, 8, 1));

    /*
    The TDI values from the first call to this function are used
    in the AddCalData() function to retrieve calibration data
    that is stored in the level-1A file.
     */
    if (firstCallThisFile)
        for (i = 0; i < 8; i++)
            tdi_global[i] = tdi[i];

    return (LIFE_IS_GOOD);
}

/****************************************************************************
Finish accesses for the current file.  Each call to CreateL1aFile()
should have a corresponding call to this function.
 *****************************************************************************/
int CloseL1aFile(l1met *metrics) {

    /* Write out the file metrics as global attributes. */
    PTB(sd_setattr(ds_id.fid,
            "Gain 1 Saturated Pixels", DFNT_INT32, 8, metrics->gain1_satpix));
    PTB(sd_setattr(ds_id.fid,
            "Gain 2 Saturated Pixels", DFNT_INT32, 8, metrics->gain2_satpix));
    PTB(sd_setattr(ds_id.fid,
            "Gain 1 Non-Saturated Pixels", DFNT_INT32, 8, metrics->gain1_nonsatpix));
    PTB(sd_setattr(ds_id.fid,
            "Gain 2 Non-Saturated Pixels", DFNT_INT32, 8, metrics->gain2_nonsatpix));
    PTB(sd_setattr(ds_id.fid, "Zero Pixels", DFNT_INT32, 8, metrics->zeropix));
    PTB(sd_setattr(ds_id.fid,
            "Mean Gain 1 Radiance", DFNT_FLOAT32, 8, metrics->gain1_mean_rad));
    PTB(sd_setattr(ds_id.fid,
            "Mean Gain 2 Radiance", DFNT_FLOAT32, 8, metrics->gain2_mean_rad));

    PTB(AddCalData());

    PTB(MakeVgroups());

    if (SDend(ds_id.fid)) {
        fprintf(stderr, "-E- %s line %d: SDend(%d) failed for file, %s.\n",
                __FILE__, __LINE__, ds_id.fid, hdfFile);
        return (HDF_FUNCTION_ERROR);
    }
    free(hdfFile);
    return (LIFE_IS_GOOD);
}

/****************************************************************************
Store tilt information for the scene in the HDF file.
 *****************************************************************************/
int AddTiltData(
        int32 ntilts,
        int16 tilt_flags[20],
        int16 tilt_ranges[20][2],
        float32 tilt_lats[20][2][2],
        float32 tilt_lons[20][2][2]
        ) {

    /* Set unused array values to something eye-catching. :=) */
    /*
    for(i = ntilts; i < 20; i++){
      tilt_flags[i] = -9999;
      for(j = 0; j < 2; j++){
        tilt_ranges[i][j] = -9999;
        for(k = 0; k < 2; k++){
          tilt_lats[i][j][k] = -9999.9;
          tilt_lons[i][j][k] = -9999.9;
        }
      }
    }
     */


    PTB(CreateSDS(ds_id.fid,
            "ntilts", /* short name */
            "Number of scene tilt states", /* long name */
            NULL, /* standard name */
            NULL, /* units */
            0, 0, /* valid range */
            1.0, 0.0, /* slope,offset */
            DFNT_INT32, /* HDF number type */
            1, /* rank */
            1, 1, 1, /* dimension sizes */
            "1", NULL, NULL /* dimension names */
            ));

    PTB(CreateSDS(ds_id.fid,
            "tilt_flags", /* short name */
            "Tilt indicators", /* long name */
            NULL, /* standard name */
            NULL, /* units */
            -1, 3, /* valid range */
            1.0, 0.0, /* slope,offset */
            DFNT_INT16, /* HDF number type */
            1, /* rank */
            20, 1, 1, /* dimension sizes */
            "20", NULL, NULL /* dimension names */
            ));

    PTB(CreateSDS(ds_id.fid,
            "tilt_ranges", /* short name */
            "Scan-line number ranges of scene tilt states", /* long name */
            NULL, /* standard name */
            NULL, /* units */
            0, 0, /* valid range */
            1.0, 0.0, /* slope,offset */
            DFNT_INT16, /* HDF number type */
            2, /* rank */
            20, 2, 1, /* dimension sizes */
            "20", "2", NULL /* dimension names */
            ));

    PTB(CreateSDS(ds_id.fid,
            "tilt_lats", /* short name */
            "Latitudes of tilt-range scan line end points", /* long name */
            NULL, /* standard name */
            NULL, /* units */
            -90.0, 90.0, /* valid range */
            1.0, 0.0, /* slope,offset */
            DFNT_FLOAT32, /* HDF number type */
            3, /* rank */
            20, 2, 2, /* dimension sizes */
            "20", "2", "2" /* dimension names */
            ));

    PTB(CreateSDS(ds_id.fid,
            "tilt_lons", /* short name */
            "Longitudes of tilt-range scan line end points", /* long name */
            NULL, /* standard name */
            NULL, /* units */
            -180.0, 180.0, /* valid range */
            1.0, 0.0, /* slope,offset */
            DFNT_FLOAT32, /* HDF number type */
            3, /* rank */
            20, 2, 2, /* dimension sizes */
            "20", "2", "2" /* dimension names */
            ));

    PTB(sd_writedata(ds_id.fid, "ntilts", &ntilts, 0, 0, 0, 1, 1, 1));
    PTB(sd_writedata(ds_id.fid, "tilt_flags", tilt_flags, 0, 0, 0, 20, 1, 1));
    PTB(sd_writedata(ds_id.fid, "tilt_ranges", tilt_ranges, 0, 0, 0, 20, 2, 1));
    PTB(sd_writedata(ds_id.fid, "tilt_lats", tilt_lats, 0, 0, 0, 20, 2, 2));
    PTB(sd_writedata(ds_id.fid, "tilt_lons", tilt_lons, 0, 0, 0, 20, 2, 2));

    return (LIFE_IS_GOOD);
}

/****************************************************************************
Append calibration data to the HDF file.  The calibration data is read from
another HDF file pointed to by an environment variable.
 *****************************************************************************/
int AddCalData(void) {

    /*
      char          *envvar = "CAL_HDF_PATH";
      char          *calPath;
      int16         entry_year, entry_day, ref_year, ref_day, ref_minute;
      float32       temps[256][8], scan_mod[2][1285], mirror[2][8];
      float64       t_const[8], t_linear[8], t_quadratic[8];
      float32       cal_offs[8], counts[8][4][5], rads[8][4][5];

      calPath = getenv(envvar);
      if(calPath == NULL){
        fprintf(stderr,"-W- %s line %d: ",__FILE__,__LINE__);
        fprintf(stderr,"Environment variable, \"%s\", not set. ", envvar);
        fprintf(stderr,"Calibration data not appended to file, %s .\n", hdfFile);
        return(CALDATA_NOT_APPENDED);
      }

      if(
      get_cal(calPath, startYear, startDay, endDay,
              startMillisec, dataTypeString, tdi_global,
              &entry_year, &entry_day, &ref_year, &ref_day, &ref_minute,
              temps, scan_mod, mirror, t_const, t_linear, t_quadratic,
              cal_offs, counts, rads) < 0
      ){
        fprintf(stderr,"-W- %s line %d: ", __FILE__,__LINE__);
        fprintf(stderr,
        "get_cal(%s,%hd,%hd,%hd,%d,\"%s\",[%hd,%hd,%hd,%hd,%hd,%hd,%hd,%hd], ...) ",
        calPath,startYear,startDay,endDay,startMillisec,dataTypeString,
        tdi_global[0],tdi_global[1],tdi_global[2],tdi_global[3],
        tdi_global[4],tdi_global[5],tdi_global[6],tdi_global[7]);
        fprintf(stderr,"failed. ");
        fprintf(stderr,"Calibration data not appended to file, %s .\n", hdfFile);
        return(CALDATA_NOT_APPENDED);
      }
     */

    PTB(CreateSDS(ds_id.fid,
            "entry_year", /* short name */
            "Calibration entry year", /* long name */
            NULL, /* standard name */
            NULL, /* units */
            0, 0, /* valid range */
            1.0, 0.0, /* slope,offset */
            DFNT_INT16, /* HDF number type */
            1, /* rank */
            1, 1, 1, /* dimension sizes */
            "1", NULL, NULL /* dimension names */
            ));

    PTB(CreateSDS(ds_id.fid,
            "entry_day", /* short name */
            "Calibration entry day-of-year", /* long name */
            NULL, /* standard name */
            NULL, /* units */
            0, 0, /* valid range */
            1.0, 0.0, /* slope,offset */
            DFNT_INT16, /* HDF number type */
            1, /* rank */
            1, 1, 1, /* dimension sizes */
            "1", NULL, NULL /* dimension names */
            ));

    PTB(CreateSDS(ds_id.fid,
            "ref_year", /* short name */
            "Calibration reference year", /* long name */
            NULL, /* standard name */
            NULL, /* units */
            0, 0, /* valid range */
            1.0, 0.0, /* slope,offset */
            DFNT_INT16, /* HDF number type */
            1, /* rank */
            1, 1, 1, /* dimension sizes */
            "1", NULL, NULL /* dimension names */
            ));

    PTB(CreateSDS(ds_id.fid,
            "ref_day", /* short name */
            "Calibration reference day-of-year", /* long name */
            NULL, /* standard name */
            NULL, /* units */
            0, 0, /* valid range */
            1.0, 0.0, /* slope,offset */
            DFNT_INT16, /* HDF number type */
            1, /* rank */
            1, 1, 1, /* dimension sizes */
            "1", NULL, NULL /* dimension names */
            ));

    PTB(CreateSDS(ds_id.fid,
            "ref_minute", /* short name */
            "Calibration reference minute-of-day", /* long name */
            NULL, /* standard name */
            NULL, /* units */
            0, 0, /* valid range */
            1.0, 0.0, /* slope,offset */
            DFNT_INT16, /* HDF number type */
            1, /* rank */
            1, 1, 1, /* dimension sizes */
            "1", NULL, NULL /* dimension names */
            ));

    PTB(CreateSDS(ds_id.fid,
            "mirror", /* short name */
            "Mirror-side correction factors", /* long name */
            NULL, /* standard name */
            NULL, /* units */
            0, 0, /* valid range */
            1.0, 0.0, /* slope,offset */
            DFNT_FLOAT32, /* HDF number type */
            2, /* rank */
            2, 8, 1, /* dimension sizes */
            "Number of Sides", "Number of Bands", NULL /* dimension names */
            ));

    PTB(CreateSDS(ds_id.fid,
            "t_const", /* short name */
            "Time-dependent correction constant terms", /* long name */
            NULL, /* standard name */
            NULL, /* units */
            0, 0, /* valid range */
            1.0, 0.0, /* slope,offset */
            DFNT_FLOAT64, /* HDF number type */
            1, /* rank */
            8, 1, 1, /* dimension sizes */
            "Number of Bands", NULL, NULL /* dimension names */
            ));

    PTB(CreateSDS(ds_id.fid,
            "t_linear", /* short name */
            "Time-dependent correction linear coefficients", /* long name */
            NULL, /* standard name */
            NULL, /* units */
            0, 0, /* valid range */
            1.0, 0.0, /* slope,offset */
            DFNT_FLOAT64, /* HDF number type */
            1, /* rank */
            8, 1, 1, /* dimension sizes */
            "Number of Bands", NULL, NULL /* dimension names */
            ));

    PTB(CreateSDS(ds_id.fid,
            "t_quadratic", /* short name */
            "Time-dependent correction quadratic coefficients", /* long name */
            NULL, /* standard name */
            NULL, /* units */
            0, 0, /* valid range */
            1.0, 0.0, /* slope,offset */
            DFNT_FLOAT64, /* HDF number type */
            1, /* rank */
            8, 1, 1, /* dimension sizes */
            "Number of Bands", NULL, NULL /* dimension names */
            ));

    PTB(CreateSDS(ds_id.fid,
            "cal_offs", /* short name */
            "Calibration system offsets", /* long name */
            NULL, /* standard name */
            NULL, /* units */
            0, 0, /* valid range */
            1.0, 0.0, /* slope,offset */
            DFNT_FLOAT32, /* HDF number type */
            1, /* rank */
            8, 1, 1, /* dimension sizes */
            "Number of Bands", NULL, NULL /* dimension names */
            ));

    PTB(CreateSDS(ds_id.fid,
            "counts", /* short name */
            "Digital counts of calibration knees", /* long name */
            NULL, /* standard name */
            NULL, /* units */
            0, 1023, /* valid range */
            1.0, 0.0, /* slope,offset */
            DFNT_FLOAT32, /* HDF number type */
            3, /* rank */
            8, 4, 5, /* dimension sizes */
            "Number of Bands", "Number of Gains", "Number of Knees" /* dimension names */
            ));

    PTB(CreateSDS(ds_id.fid,
            "rads", /* short name */
            "Radiances of calibration knees", /* long name */
            NULL, /* standard name */
            NULL, /* units */
            0, 0, /* valid range */
            1.0, 0.0, /* slope,offset */
            DFNT_FLOAT32, /* HDF number type */
            3, /* rank */
            8, 4, 5, /* dimension sizes */
            "Number of Bands", "Number of Gains", "Number of Knees" /* dimension names */
            ));

    /*
      PTB( sd_writedata("entry_year"  , &entry_year , 0, 0, 0, 1, 1, 1)  );
      PTB( sd_writedata("entry_day"   , &entry_day  , 0, 0, 0, 1, 1, 1)  );
      PTB( sd_writedata("ref_year"    , &ref_year   , 0, 0, 0, 1, 1, 1)  );
      PTB( sd_writedata("ref_day"     , &ref_day    , 0, 0, 0, 1, 1, 1)  );
      PTB( sd_writedata("ref_minute"  , &ref_minute , 0, 0, 0, 1, 1, 1)  );
      PTB( sd_writedata("mirror"      , mirror      , 0, 0, 0, 2, 8, 1)  );
      PTB( sd_writedata("t_const"     , t_const     , 0, 0, 0, 8, 1, 1)  );
      PTB( sd_writedata("t_linear"    , t_linear    , 0, 0, 0, 8, 1, 1)  );
      PTB( sd_writedata("t_quadratic" , t_quadratic , 0, 0, 0, 8, 1, 1)  );
      PTB( sd_writedata("cal_offs"    , cal_offs    , 0, 0, 0, 8, 1, 1)  );
      PTB( sd_writedata("counts"      , counts      , 0, 0, 0, 8, 4, 5)  );
      PTB( sd_writedata("rads"        , rads        , 0, 0, 0, 8, 4, 5)  );
     */

    calibrationAppended = 1; /* signal to MakeVgroups() */

    return (LIFE_IS_GOOD);
}

/****************************************************************************
Associate various Scientific Data Sets into Vgroups.  I don't know what
useful purpose this serves.
 *****************************************************************************/
int MakeVgroups(void) {

    int32 h_id, v_id;

    h_id = Hopen(hdfFile, DFACC_RDWR, 0);
    if (h_id == FAIL) {
        fprintf(stderr, "-E- %s line %d: Hopen() failed for file, %s.\n",
                __FILE__, __LINE__, hdfFile);
        return (HDF_FUNCTION_ERROR);
    }
    Vstart(h_id);

    /* Scan-Line Attributes */
    PTB(v_attach(h_id, &v_id));
    Vsetclass(v_id, "Per Scan Data");
    Vsetname(v_id, "Scan-Line Attributes");
    PTB(AddSdsToVgroup(ds_id.fid, v_id, "msec"));
    PTB(AddSdsToVgroup(ds_id.fid, v_id, "eng_qual"));
    PTB(AddSdsToVgroup(ds_id.fid, v_id, "s_flags"));
    PTB(AddSdsToVgroup(ds_id.fid, v_id, "s_satp"));
    PTB(AddSdsToVgroup(ds_id.fid, v_id, "s_zerop"));
    PTB(AddSdsToVgroup(ds_id.fid, v_id, "slat"));
    PTB(AddSdsToVgroup(ds_id.fid, v_id, "slon"));
    PTB(AddSdsToVgroup(ds_id.fid, v_id, "clat"));
    PTB(AddSdsToVgroup(ds_id.fid, v_id, "clon"));
    PTB(AddSdsToVgroup(ds_id.fid, v_id, "elat"));
    PTB(AddSdsToVgroup(ds_id.fid, v_id, "elon"));
    PTB(AddSdsToVgroup(ds_id.fid, v_id, "csol_z"));
    PTB(AddSdsToVgroup(ds_id.fid, v_id, "tilt"));
    Vdetach(v_id);

    /* Raw SeaStar Data */
    PTB(v_attach(h_id, &v_id));
    Vsetclass(v_id, "Per Scan Data");
    Vsetname(v_id, "Raw SeaStar Data");
    PTB(AddSdsToVgroup(ds_id.fid, v_id, "sc_id"));
    PTB(AddSdsToVgroup(ds_id.fid, v_id, "sc_ttag"));
    PTB(AddSdsToVgroup(ds_id.fid, v_id, "sc_soh"));
    PTB(AddSdsToVgroup(ds_id.fid, v_id, "inst_tlm"));
    PTB(AddSdsToVgroup(ds_id.fid, v_id, "l1a_data"));
    PTB(AddSdsToVgroup(ds_id.fid, v_id, "start_syn"));
    PTB(AddSdsToVgroup(ds_id.fid, v_id, "stop_syn"));
    PTB(AddSdsToVgroup(ds_id.fid, v_id, "dark_rest"));
    PTB(AddSdsToVgroup(ds_id.fid, v_id, "gain"));
    PTB(AddSdsToVgroup(ds_id.fid, v_id, "tdi"));
    Vdetach(v_id);

    /* Converted Telemetry */
    PTB(v_attach(h_id, &v_id));
    Vsetclass(v_id, "Per Scan Data");
    Vsetname(v_id, "Converted Telemetry");
    PTB(AddSdsToVgroup(ds_id.fid, v_id, "inst_ana"));
    PTB(AddSdsToVgroup(ds_id.fid, v_id, "inst_dis"));
    PTB(AddSdsToVgroup(ds_id.fid, v_id, "sc_ana"));
    PTB(AddSdsToVgroup(ds_id.fid, v_id, "sc_dis"));
    PTB(AddSdsToVgroup(ds_id.fid, v_id, "scan_temp"));
    PTB(AddSdsToVgroup(ds_id.fid, v_id, "side"));
    Vdetach(v_id);

    /* Navigation */
    PTB(v_attach(h_id, &v_id));
    Vsetclass(v_id, "Per Scan Data");
    Vsetname(v_id, "Navigation");
    PTB(AddSdsToVgroup(ds_id.fid, v_id, "orb_vec"));
    PTB(AddSdsToVgroup(ds_id.fid, v_id, "l_vert"));
    PTB(AddSdsToVgroup(ds_id.fid, v_id, "sun_ref"));
    PTB(AddSdsToVgroup(ds_id.fid, v_id, "att_ang"));
    PTB(AddSdsToVgroup(ds_id.fid, v_id, "sen_mat"));
    PTB(AddSdsToVgroup(ds_id.fid, v_id, "scan_ell"));
    PTB(AddSdsToVgroup(ds_id.fid, v_id, "nflag"));
    Vdetach(v_id);

    /* Sensor Tilt */
    PTB(v_attach(h_id, &v_id));
    Vsetclass(v_id, "Per File Data");
    Vsetname(v_id, "Sensor Tilt");
    PTB(AddSdsToVgroup(ds_id.fid, v_id, "ntilts"));
    PTB(AddSdsToVgroup(ds_id.fid, v_id, "tilt_flags"));
    PTB(AddSdsToVgroup(ds_id.fid, v_id, "tilt_ranges"));
    PTB(AddSdsToVgroup(ds_id.fid, v_id, "tilt_lats"));
    PTB(AddSdsToVgroup(ds_id.fid, v_id, "tilt_lons"));
    Vdetach(v_id);

    /* Calibration */
    if (calibrationAppended) {
        PTB(v_attach(h_id, &v_id));
        Vsetclass(v_id, "Per File Data");
        Vsetname(v_id, "Calibration");
        PTB(AddSdsToVgroup(ds_id.fid, v_id, "entry_year"));
        PTB(AddSdsToVgroup(ds_id.fid, v_id, "entry_day"));
        PTB(AddSdsToVgroup(ds_id.fid, v_id, "ref_year"));
        PTB(AddSdsToVgroup(ds_id.fid, v_id, "ref_day"));
        PTB(AddSdsToVgroup(ds_id.fid, v_id, "ref_minute"));
        PTB(AddSdsToVgroup(ds_id.fid, v_id, "mirror"));
        PTB(AddSdsToVgroup(ds_id.fid, v_id, "t_const"));
        PTB(AddSdsToVgroup(ds_id.fid, v_id, "t_linear"));
        PTB(AddSdsToVgroup(ds_id.fid, v_id, "t_quadratic"));
        PTB(AddSdsToVgroup(ds_id.fid, v_id, "cal_offs"));
        PTB(AddSdsToVgroup(ds_id.fid, v_id, "counts"));
        PTB(AddSdsToVgroup(ds_id.fid, v_id, "rads"));
        Vdetach(v_id);
    }

    Vend(h_id);

    if (Hclose(h_id) != SUCCEED) {
        fprintf(stderr, "-E- %s line %d: Hclose(%d) failed for file, %s .\n",
                __FILE__, __LINE__, h_id, hdfFile);
        return (HDF_FUNCTION_ERROR);
    }
    return (LIFE_IS_GOOD);
}

/****************************************************************************
Construct a SeaWiFS level-1A filename from a time value and a data-type
value.  Return a pointer to the statically allocated filename string.
 *****************************************************************************/
char * L1aFilename(swl0ctl *l0ctl, double time, unsigned char dataType) {

    static char filename[24]; /* "Syyyydddhhmmss.L1A_tttt\0" */
    struct tm *t;
    time_t itime;
    double rint(double);

    itime = (time_t) rint(time); /* Round to nearest second. */
    t = gmtime(&itime);

    sprintf(
            filename,
            "S%4d%03d%02d%02d%02d.L1A_%.4s",
            t->tm_year + 1900,
            t->tm_yday + 1,
            t->tm_hour,
            t->tm_min,
            t->tm_sec,
            DataTypeString(l0ctl, dataType)
            );

    return (filename);
}

/****************************************************************************
Return a three- or four-character string to represent the input data type.
The input argument has the data-type value from the spacecraft ID in the
level-0 data or a value of 16 for HRPT data.  Unknown data-type values
cause an empty string, "", to be returned.
 *****************************************************************************/
char * DTypeString(unsigned char dataType) {
    switch (dataType) {
    case 0: return ("LAC");
    case 1: return ("LUN");
    case 2: return ("SOL");
    case 3: return ("IGC");
    case 4: return ("TDI");
    case 15: return ("GAC");
    case 16: return ("HRPT");
    default: return ("");
    }
}

/****************************************************************************
This function is essentially the same as DTypeString() except that
it returns "Hsta" when passed a 16, where "sta" is the 3-letter station
code as defined in the HRPT_STATION_IDENTIFICATION_FILE.  Unknown data-type
values are converted to an ASCII string and returned.
 *****************************************************************************/
char * DataTypeString(swl0ctl *l0ctl, unsigned char dataType) {

    static char type[5];
    char *t;
    int status;

    t = DTypeString(dataType);
    if (strcmp(t, "HRPT") == 0) {
        StationInfo stationInfo;
        status = GetStationInfo(l0ctl->stationInfoFile, &stationInfo);
        if (status != LIFE_IS_GOOD || stationInfo.code[0] == 0) {
            /* Code not found in StationInfo file. */
            sprintf(type, "Hxxx");
            fprintf(stderr, "-W- %s line %d: ", __FILE__, __LINE__);
            fprintf(stderr, "Station code not found; using \"xxx\".\n");
        } else {
            if (strlen(stationInfo.code) == 3)
                sprintf(type, "H%.3s", stationInfo.code);
            else
                sprintf(type, "%.4s", stationInfo.code);
        }
        return (type);
    } else if (*t == 0) {
        sprintf(type, "%03u", dataType);
        fprintf(stderr, "-W- %s line %d: ", __FILE__, __LINE__);
        fprintf(stderr, "Unknown data type.  ");
        fprintf(stderr, "Setting data-type string to \"%s\".\n", type);
        return (type);
    } else {
        return (t);
    }
}

/****************************************************************************
Return the year, day-of-year, and millisecond-of-day of the passed in
number of seconds since 1-Jan-1970 00:00:00.000 GMT .
 *****************************************************************************/
void DecomposeTime(
        double dtime,
        int16 *year,
        int16 *dayofyear,
        int32 *millisec
        ) {
    time_t itime = (time_t) dtime;
    struct tm *ts;
    double rint(double);

    ts = gmtime(&itime);
    *year = (int16) (ts->tm_year + 1900);
    *dayofyear = (int16) (ts->tm_yday + 1);
    *millisec = (int32) (floor(1000 * (ts->tm_hour * 3600 +
            ts->tm_min * 60 +
            ts->tm_sec +
            dtime - itime)));
}
