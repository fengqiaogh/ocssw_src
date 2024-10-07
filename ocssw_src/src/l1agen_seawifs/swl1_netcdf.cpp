/* Wei S. Jiang        Aug 15, 2023 */
#include <netcdf>
#include <iostream>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h>
#include <time.h>
#include <math.h>
#include <GetStationInfo.h>
#include "swl0_proto.h"
#include <timeutils.h>
#include "global_attrs.h" 

using namespace std;
using namespace netCDF;
using namespace netCDF::exceptions;


/*
These global variables facilitate communication
between some of the functions defined in this file.
 */
static int32 numScans, numPixels;
static int16 tdi_global[8];
static char *netcdfFile; // name of the netcdf file for output
static int firstCallThisFile;
static char dataTypeString[5];
static int16 startYear, startDay, endDay;
static int32 startMillisec;
static int calibrationAppended = 0;
static NcFile *dataFile; // global NetCDF file for calls outside the initial call

#define SENSOR          "Sea-viewing Wide Field-of-view Sensor (SeaWiFS)"
#define MISSIONCHAR     "Nominal orbit: inclination = 98.2 (Sun-synchronous); node = 12 noon local (descending); eccentricity = <0.002; altitude = 705 km; ground speed = 6.75 km/sec"
#define SENSORCHAR      "Number of bands = 8; number of active bands = 8; wavelengths per band (nm) = 412, 443, 490, 510, 555, 670, 765, 865; bits per pixel = 10; instantaneous field-of-view = 1.5835 mrad; pixels per scan = 1285; scan rate = 6/sec; sample rate = 7710/sec"


/****************************************************************************
 * Given a NcVar, add attributes passed into the parameters.
 * @param variable - reference to the NcVar that you will be adding attributes to
 * @param longName - long name description of the attribute 
 * @param unit - unit of the variable
 * @param hasValidRange - toggle for if min-max range is given. Add if true. 
 * @param validMin - min value the variable can take on if it has a valid range (default=0)
 * @param validMax - max value the variable can take on.
 * @param validMinMaxType - NetCDF type if min max is given
 * 
 * The valid range takes ints and can be converted to float, short, double, etc. 
 * as long as if the valid range's float/double does not contain percision data
 * 
 * ie. You can pass in 15 for 15.0 and NetCDF will convert it to 15.f via validMinMaxType
 * However, 15.1 will not work because the .1 is lost for int type.
*****************************************************************************/
void addAttribute_netcdf(
	NcVar &variable, 
	string longName, 		
	string unit,		
	bool hasValidRange, 	
	int validMin,			
	int validMax,			
	NcType validMinMaxType 
	) {
		// All variables have long names, no need to check if it's null
		variable.putAtt("long_name", longName);

		if (!unit.empty()) {
			variable.putAtt("unit", unit);
		}

		if (hasValidRange) {
			variable.putAtt("valid_min", validMinMaxType, validMin);
			variable.putAtt("valid_max", validMinMaxType, validMax);
		}
}

/****************************************************************************
 * Overloaded function to take in doubles as valid range. 
 * Original takes ints, but can be converted to shorts, floats or double
 * if the valid min/max does not have additional percision data.
 * 
 * This function takes care of instances that the original function cant 
 * get away with using int. 
 * ie. While 15.0 you can get away with passing in as 15, you cant for 15.1 or 
 * 15.12
*****************************************************************************/
void addAttribute_netcdf(
	NcVar &variable,  		
	string longName, 		
	string unit,			
	bool hasValidRange, 	
	double validMin,			
	double validMax,			
	NcType validMinMaxType 
	) {
		// All variables have long names, no need to check if it's null
		variable.putAtt("long name", longName);

		if (!unit.empty()) {
			variable.putAtt("unit", unit);
		}

		if (hasValidRange) {
			variable.putAtt("valid_min", validMinMaxType, validMin);
			variable.putAtt("valid_max", validMinMaxType, validMax);
		}
}


/****************************************************************************
Create a new NetCDF file and store some global attributes in it.
This function must be called before WriteScanData().
CloseL1aFile() should be called to finish up the file.
 *****************************************************************************/
extern "C" int CreateL1aFile_netcdf(
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
    MALLOC(netcdfFile, char, strlen(path) + 1);
    strcpy(netcdfFile, path);



/* 
        Create NetCDF File.
*/

	try {
		// Create the NetCDF file, 'replace' will ovrride the file if it exists already
		dataFile = new NcFile(netcdfFile, NcFile::replace);

		/***
			Set dimensions and keep its ref to set variable and attributes later
		***/ 
		NcDim dim_numTilts = dataFile->addDim("tilts", 20);
		NcDim dim_tiltRange = dataFile->addDim("tilt_ranges", 2);
		NcDim dim_numScans = dataFile->addDim("scans", numScans);
		NcDim dim_qualFlag = dataFile->addDim("quality_flags", 4); // 4 from old WriteScanData() 
		NcDim dim_numBands = dataFile->addDim("bands", 8); // 8 is from SENSORCHAR
		NcDim dim_scIds = dataFile->addDim("sc_ids", 2); 
		NcDim dim_scTimeTags = dataFile->addDim("sc_time_tags", 4); 
		NcDim dim_scSohs = dataFile->addDim("sc_sohs", 775); // todo, migh not be static
		NcDim dim_instTeleme = dataFile->addDim("inst_telemetries", 44); 
		NcDim dim_pixlePerLine = dataFile->addDim("pixels", numPixels); 
		NcDim dim_instAnaTeleme = dataFile->addDim("inst_analog_telemetries", 40); 
		NcDim dim_instDiscTeleme = dataFile->addDim("inst_discrete_telemetries", 32); 
		NcDim dim_scAnaTele = dataFile->addDim("sc_analog_telemetries", 40); 
		NcDim dim_scDiscTele = dataFile->addDim("sc_discrete_telemetries", 40);
		NcDim dim_vectorEle = dataFile->addDim("vector_elements", 3); 
		NcDim dim_scanTrackCoeff = dataFile->addDim("scan_track_coefficients", 6); 
		NcDim dim_navFlags = dataFile->addDim("navigation_flags", 8); 
		NcDim dim_numSides = dataFile->addDim("mirror_sides", 2);
		NcDim dim_numGains = dataFile->addDim("gains", 4); 
		NcDim dim_numKnees = dataFile->addDim("knees", 5); 
		
		/***
			Create variables in the NetCDF file and attach the diemtions that they use
			to it. Then, add attributes describing the variable and attach the data to 
			it if available. 

			If no data is attached right after adding the attributes, then it is 
			attached scan by scan in WriteScanData_netcdf() function
		***/ 

		//parameters that is reused by > 1 variables
        vector<NcDim> para_numTilt_tiltRange_tiltRange{dim_numTilts, dim_tiltRange, dim_tiltRange};
		vector<NcDim> para_numScan_numBands{dim_numScans, dim_numBands};
		vector<NcDim> para_numScan_qualFlag{dim_numScans, dim_qualFlag};
		vector<NcDim> para_numScans_vectorElements{dim_numScans, dim_vectorEle};
		vector<NcDim> para_numBands_numGains_numKnees{dim_numBands, dim_numGains, dim_numKnees};
		


		// tilt_flag()
		NcVar tiltFlag = dataFile->addVar("tilt_flags", NC_SHORT, dim_numTilts); 
		addAttribute_netcdf(
			tiltFlag,  				// reference to the variable
			"Tilt indicators", 		// long name
			"", 					// unit
			true, 					// has a min-max range					
			-1, 					// min value									
			3,						// max value														
			NC_SHORT				// min-max range type						
		);
		tiltFlag.putVar(scene->tilt_flags); // append the data to the variable


		// tilt ranges()
		vector<NcDim> tiltRangeParameters{dim_numTilts, dim_tiltRange};
		NcVar tiltRange = dataFile->addVar("tilt_ranges", NC_SHORT, tiltRangeParameters);
		addAttribute_netcdf(
			tiltRange,  									// NcVar reference
			"Scan-line number ranges of scene tilt states", // long name
			"", 											// unit
			false, 											// no min/xax, ignore everything  below	
			0, 											
			0,										
			NC_NAT											
		);
		tiltRange.putVar(scene->tilt_ranges); // append the data to the variable


		// tilt_lats()
		//vector<NcDim> tiltLatParameters{dim_numTilts, dim_tiltRange, dim_tiltRange};
		NcVar tiltLats = dataFile->addVar("tilt_lats", NC_FLOAT, para_numTilt_tiltRange_tiltRange);
		addAttribute_netcdf(
			tiltLats,  										
			"Latitudes of tilt-range scan line end points", 
			"", 											
			true, 											
			-90, 			
			90,											
			NC_FLOAT										
		);
		tiltLats.putVar(scene->tilt_lats);


		// tilt_lons
		//vector<NcDim> tiltLonParameters{dim_numTilts, dim_tiltRange, dim_tiltRange};
		NcVar tiltLons = dataFile->addVar("tilt_lons", NC_FLOAT, para_numTilt_tiltRange_tiltRange);
		addAttribute_netcdf(
			tiltLons, 
			"Longitudes of tilt-range scan line end points",
			"",
			true,
			-180,
			180,
			NC_FLOAT
		);
		tiltLons.putVar(scene->tilt_lons);


		// scan_time **not in the old data structure (maybe msec?)** 
		NcVar scanTime = dataFile->addVar("scan_time", NC_INT, dim_numScans);
		//scanTime.putAtt("_FillValue", NC_DOUBLE, -999.9);  // only one with "_FillValue" so doing it outside of func
		char* starttime = {unix2isodate(scene->stime, 'G')};
		string starttime_str(starttime);

		addAttribute_netcdf(
			scanTime,
			"Scan start time, milliseconds of day",
			"milliseconds since " + starttime_str.substr(0,10),
			true,
			0,
			86399999,
			NC_INT
		);
		free(starttime);


		// eng_qual()
		NcVar engQual = dataFile->addVar("eng_qual", NC_UBYTE, para_numScan_qualFlag);
		addAttribute_netcdf(
			engQual,
			"Engineering data-out-of-range flags",
			"",
			false,
			0,
			0,
			NC_NAT
		);

		// s_flag()
		NcVar sFlag = dataFile->addVar("s_flag", NC_UBYTE, para_numScan_qualFlag);
		addAttribute_netcdf(
			sFlag,
			"Scan-line quality flags",
			"",
			false,
			0,
			0,
			NC_NAT
		);

		// s_satp()
		NcVar s_satP= dataFile->addVar("s_satp", NC_SHORT, para_numScan_numBands);
		addAttribute_netcdf(
			s_satP,
			"Number of saturated pixels per band",
			"",
			false,
			0,
			0,
			NC_NAT
		);

		// s_zerop()
		NcVar s_zerop= dataFile->addVar("s_zerop", NC_SHORT, para_numScan_numBands);
		addAttribute_netcdf(
			s_zerop,
			"Number of zero pixels per band",
			"",
			false,
			0,
			0,
			NC_NAT
		);

		// slat()
		NcVar slat = dataFile->addVar("slat", NC_FLOAT, dim_numScans);
		addAttribute_netcdf(
			slat,
			"Scan start-pixel latitude",
			"",
			true,
			-90,
			90,
			NC_FLOAT
		);

		// slon()
		NcVar slon = dataFile->addVar("slon", NC_FLOAT, dim_numScans);
		addAttribute_netcdf(
			slon,
			"Scan start-pixel longitude",
			"",
			true,
			-180,
			180,
			NC_FLOAT
		);

		// clat()
		NcVar clat = dataFile->addVar("clat", NC_FLOAT, dim_numScans);
		addAttribute_netcdf(
			clat,
			"Scan center-pixel longitude",
			"",
			true,
			-90,
			90,
			NC_FLOAT
		);

		// clon()
		NcVar clon = dataFile->addVar("clon", NC_FLOAT, dim_numScans);
		addAttribute_netcdf(
			clon,
			"Scan center-pixel longitude",
			"",
			true,
			-180,
			180,
			NC_FLOAT
		);

		// elat()
		NcVar elat = dataFile->addVar("elat", NC_FLOAT, dim_numScans);
		addAttribute_netcdf(
			elat,
			"Scan end-pixel longitude",
			"",
			true,
			-90,
			90,
			NC_FLOAT
		);

		// elon()
		NcVar elon = dataFile->addVar("elon", NC_FLOAT, dim_numScans);
		addAttribute_netcdf(
			elon,
			"Scan end-pixel longitude",
			"",
			true,
			-180,
			180,
			NC_FLOAT
		);

		// csol()
		NcVar csol_z = dataFile->addVar("csol_z", NC_FLOAT, dim_numScans);
		addAttribute_netcdf(
			csol_z,
			"Scan center-pixel solar zenith angle",
			"",
			true,
			0,
			180,
			NC_FLOAT
		);

		// tilt()
		NcVar tilt = dataFile->addVar("tilt", NC_FLOAT, dim_numScans);
		addAttribute_netcdf(
			tilt,
			"Tilt angle for scan line",
			"",
			true,
			-20.1,
			20.1,
			NC_FLOAT
		);

		// sc_id() *Spacecraft Id*
		vector<NcDim> scIdParameter{dim_numScans, dim_scIds};
		NcVar sc_id = dataFile->addVar("sc_id", NC_SHORT, scIdParameter);
		addAttribute_netcdf(
			sc_id,
			"Spacecraft ID",
			"",
			false,
			0,
			0,
			NC_NAT
		);

		// sc_ttag()
		vector<NcDim> scTtagParameter{dim_numScans, dim_scTimeTags};
		NcVar sc_ttag = dataFile->addVar("sc_ttag", NC_SHORT, scTtagParameter);
		addAttribute_netcdf(
			sc_ttag,
			"Spacecraft time tag",
			"",
			false,
			0,
			0,
			NC_NAT
		);

		// sc_soh()
		vector<NcDim> scSohParameter{dim_numScans, dim_scSohs};
		NcVar sc_soh = dataFile->addVar("sc_soh", NC_UBYTE, scSohParameter);
		addAttribute_netcdf(
			sc_soh,
			"Spacecraft state-of-health data",
			"",
			false,
			0,
			0,
			NC_NAT
		);

		// inst_tlm()
		vector<NcDim> instTlmParameter{dim_numScans, dim_instTeleme};
		NcVar sc_tlm = dataFile->addVar("inst_tlm", NC_SHORT, instTlmParameter);
		addAttribute_netcdf(
			sc_tlm,
			"SeaWiFS instrument telemetry",
			"",
			false,
			0,
			0,
			NC_NAT
		);

		// l1a_data()
		vector<NcDim> l1aDataParameter{dim_numScans, dim_pixlePerLine, dim_numBands};
		NcVar l1a_data = dataFile->addVar("l1a_data", NC_SHORT, l1aDataParameter);
		addAttribute_netcdf(
			l1a_data,
			"Level-1A data",
			"radiance counts",
			true,
			0,
			1023,
			NC_SHORT
		);

		// start_syn()
		NcVar startSyn = dataFile->addVar("start_syn", NC_SHORT, para_numScan_numBands);
		addAttribute_netcdf(
			startSyn,
			"Start-synch pixel",
			"",
			false,
			0,
			0,
			NC_NAT
		);

		// start_syn() note: uses same parameters as start_sync
		NcVar stopSyn = dataFile->addVar("stop_syn", NC_SHORT, para_numScan_numBands);
		addAttribute_netcdf(
			stopSyn,
			"Stop-synch pixel",
			"",
			false,
			0,
			0,
			NC_NAT
		);

		// dark_rest()
		NcVar darkRest = dataFile->addVar("dark_rest", NC_SHORT, para_numScan_numBands);
		addAttribute_netcdf(
			darkRest,
			"Dark-restore pixel",
			"",
			false,
			0,
			0,
			NC_NAT
		);

		// gain()
		NcVar gain = dataFile->addVar("gain", NC_SHORT, para_numScan_numBands);
		addAttribute_netcdf(
			gain,
			"Band gain settings",
			"",
			true,
			0,
			3,
			NC_SHORT
		);

		// tdi()
		NcVar tdi = dataFile->addVar("tdi", NC_SHORT, para_numScan_numBands);
		addAttribute_netcdf(
			tdi,
			"Band time-delay and integration settings",
			"",
			true,
			0,
			255,
			NC_SHORT
		);

		// inst_ana()
		vector<NcDim> instAnaParameter{dim_numScans, dim_instAnaTeleme};
		NcVar instAna = dataFile->addVar("inst_ana", NC_FLOAT, instAnaParameter);
		addAttribute_netcdf(
			instAna,
			"Instrument analog telemetry",
			"",
			false,
			0,
			0,
			NC_FLOAT
		);

		// inst_dis()
		vector<NcDim> instDisParameter{dim_numScans, dim_instDiscTeleme};
		NcVar instDis = dataFile->addVar("inst_dis", NC_UBYTE, instDisParameter);
		addAttribute_netcdf(
			instDis,
			"Instrument discrete telemetry",
			"",
			false,
			0,
			0,
			NC_NAT
		);

		// sc_ana()
		vector<NcDim> scAnaParameter{dim_numScans, dim_scAnaTele};
		NcVar scAna = dataFile->addVar("sc_ana", NC_FLOAT, scAnaParameter);
		addAttribute_netcdf(
			scAna,
			"Spacecraft analog telemetry",
			"",
			false,
			0,
			0,
			NC_NAT
		);

		// sc_dis()
		vector<NcDim> scDisParameter{dim_numScans, dim_scDiscTele};
		NcVar scDis = dataFile->addVar("sc_dis", NC_UBYTE, scDisParameter);
		addAttribute_netcdf(
			scDis,
			"Spacecraft discrete telemetry",
			"",
			false,
			0,
			0,
			NC_NAT
		);

		// scan_temp()
		NcVar scanTemp = dataFile->addVar("scan_temp", NC_SHORT, para_numScan_numBands);
		addAttribute_netcdf(
			scanTemp,
			"Detector temperature counts",
			"",
			true,
			0,
			255,
			NC_SHORT
		);

		// side()
		NcVar side = dataFile->addVar("side", NC_SHORT, dim_numScans);
		addAttribute_netcdf(
			side,
			"Mirror side for scan line",
			"",
			true,
			0,
			1,
			NC_SHORT
		);

		// orb_vec()
		NcVar orbVec = dataFile->addVar("orb_vec", NC_FLOAT, para_numScans_vectorElements);
		addAttribute_netcdf(
			orbVec,
			"Orbit position vector at scan line time",
			"kilometers",
			true,
			-7200,
			7200,
			NC_FLOAT
		);

		// l_vert()
		NcVar l_vert = dataFile->addVar("l_vert", NC_FLOAT,para_numScans_vectorElements);
		addAttribute_netcdf(
			l_vert,
			"Local vertical vector in ECEF frame",
			"kilometers",
			true,
			-1,
			1,
			NC_FLOAT
		);

		// sun_ref()
		NcVar sun_ref = dataFile->addVar("sun_ref", NC_FLOAT, para_numScans_vectorElements);
		addAttribute_netcdf(
			sun_ref,
			"Reference Sun vector in ECEF frame",
			"",
			true,
			-1,
			1,
			NC_FLOAT
		);

		// att_ang()
		NcVar att_ang = dataFile->addVar("att_ang", NC_FLOAT, para_numScans_vectorElements);
		addAttribute_netcdf(
			att_ang,
			"Computed yaw, roll, pitch",
			"",
			true,
			-180,
			180,
			NC_FLOAT
		);

		// sen_mat()
		vector<NcDim> sen_mat_Parameter{dim_numScans, dim_vectorEle, dim_vectorEle};
		NcVar sen_mat = dataFile->addVar("sen_mat", NC_FLOAT, sen_mat_Parameter);
		addAttribute_netcdf(
			sen_mat,
			"ECEF-to-sensor-frame matrix",
			"",
			true,
			-1,
			1,
			NC_FLOAT
		);

		// scan_ell()
		vector<NcDim> scan_ell_Parameter{dim_numScans, dim_scanTrackCoeff};
		NcVar scan_ell = dataFile->addVar("scan_ell", NC_FLOAT, scan_ell_Parameter);
		addAttribute_netcdf(
			scan_ell,
			"Scan-track ellipse coefficients",
			"",
			false,
			0,
			0,
			NC_NAT
		);

		// nflag()
		vector<NcDim> nFlagParameter{dim_numScans, dim_navFlags};
		NcVar nFlag = dataFile->addVar("nflag", NC_INT, nFlagParameter);
		addAttribute_netcdf(
			nFlag,
			"Navigation flags",
			"",
			false,
			0,
			0,
			NC_NAT
		);

		// mirror()
		vector<NcDim> mirrorParameter{dim_numSides, dim_numBands};
		NcVar mirror = dataFile->addVar("mirror", NC_FLOAT, mirrorParameter);
		addAttribute_netcdf(
			mirror,
			"Mirror-side correction factors",
			"",
			false,
			0,
			0,
			NC_NAT
		);

		// t_const()
		NcVar t_const = dataFile->addVar("t_const", NC_DOUBLE, dim_numBands);
		addAttribute_netcdf(
			t_const,
			"Time-dependent correction constant terms",
			"",
			false,
			0,
			0,
			NC_NAT
		);

		// t_linear()
		NcVar t_linear = dataFile->addVar("t_linear", NC_DOUBLE, dim_numBands);
		addAttribute_netcdf(
			t_linear,
			"Time-dependent correction linear coefficients",
			"",
			false,
			0,
			0,
			NC_NAT
		);

		// t_quadratic()
		NcVar t_quadratic = dataFile->addVar("t_quadratic", NC_DOUBLE, dim_numBands);
		addAttribute_netcdf(
			t_quadratic,
			"Time-dependent correction quadratic coefficients",
			"",
			false,
			0,
			0,
			NC_NAT
		);

		// cal_offs()
		NcVar cal_offs = dataFile->addVar("cal_offs", NC_FLOAT, dim_numBands);
		addAttribute_netcdf(
			cal_offs,
			"Calibration system offsets",
			"",
			false,
			0,
			0,
			NC_NAT
		);

		// counts()
		NcVar counts = dataFile->addVar("counts", NC_FLOAT, para_numBands_numGains_numKnees);
		addAttribute_netcdf(
			counts,
			"Digital counts of calibration knees",
			"",
			true,
			0,
			1023,
			NC_FLOAT
		);


		// rads()
		NcVar rads = dataFile->addVar("rads", NC_FLOAT, para_numBands_numGains_numKnees);
		addAttribute_netcdf(
			rads,
			"Radiances of calibration knees",
			"",
			false,
			0,
			0,
			NC_NAT
		);

		
		// Set global seawif attributes

		int pixelOffset = 0;  // required as putting 0 into putAtt throws an error
		dataFile->putAtt("title", "SeaWiFS Level-1A Data");
		dataFile->putAtt("instrument", "SeaWiFS");
		dataFile->putAtt("platform", "Orbview-2");
		dataFile->putAtt("product_name", basename(netcdfFile));
		dataFile->putAtt("processing_version", "V2");
		dataFile->putAtt("processing_level", "L1A");
		dataFile->putAtt("cdm_data_type", "swath");
		dataFile->putAtt("orbit_number", NC_INT, scene->orbnum);  // NC_INT is 4 bytes == 32 bits
		// dataFile->putAtt("history", " "); // need history .nc file to use the util. 
		dataFile->putAtt("time_coverage_start", unix2isodate(scene->stime, 'G'));
		dataFile->putAtt("time_coverage_end", unix2isodate(scene->etime, 'G'));
		dataFile->putAtt("startDirection", scene->start_node);
		dataFile->putAtt("endDirection", scene->end_node);
		// dataFile->putAtt("day_night_flag", "Day");
		dataFile->putAtt("data_type", DTypeString(dataType));
		// dataFile->putAtt("replacement_flag", basename(netcdfFile)); // same as product_name
		dataFile->putAtt("LAC_pixel_start_number", NC_INT, startpix); 
		dataFile->putAtt("LAC_pixel_subsampling", NC_INT, subsamp);
		dataFile->putAtt("pixel_offset", NC_INT, pixelOffset);
		dataFile->putAtt("data_center", stationInfo.data_center);
		dataFile->putAtt("station_name", stationInfo.station_name);
		dataFile->putAtt("station_latitude",NC_FLOAT, stationInfo.station_latitude);
		dataFile->putAtt("station_longitude", NC_FLOAT, stationInfo.station_longitude);
		dataFile->putAtt("mission", "SeaStar SeaWiFS");
		dataFile->putAtt("mission_characteristics", MISSIONCHAR);
		dataFile->putAtt("sensor", SENSOR);
		dataFile->putAtt("sensor_characteristics", SENSORCHAR);
		dataFile->putAtt("input_files", scene->l0file);
		dataFile->putAtt("scene_center_time", unix2isodate(scene->ctime, 'G'));
		dataFile->putAtt("node_crossing_time", unix2isodate(scene->node_time, 'G'));
		// dataFile->putAtt("calibration_entry_time", "...");
		// dataFile->putAtt("calibration_reference_time", "...");
		dataFile->putAtt("orbit_node_longitude", NC_FLOAT, scene->node_lon);

		// Set standard global attributes
		string history = get_history(dataFile);
		string doi = "";
		set_global_attrs(dataFile, history, doi);

	} catch(NcException& e) {
		cerr << "-E- Error with creating NetCDF File: " << e.what() << endl;
		return 1; // main file checks if it is not 0 for error
	}

    /*
    Set global variables that are also used
    in getting calibration data for the scene.
     */
    strncpy(dataTypeString, DTypeString(dataType), 4);
    DecomposeTime(scene->stime, &startYear, &startDay, &startMillisec);
    DecomposeTime(scene->etime, &year, &endDay, &millisec);
    firstCallThisFile = 0;

    return (LIFE_IS_GOOD);
}

/****************************************************************************
 * Given a NcVar's name, find it in the current opened NetCDF file.
 * Put the single value data into the file.
 * 
 * Function is using a void pointer to data because you will be passing it a 
 * reference to a data that's already initialized. 
 * 
 * @param ncVarName             name of the NcVar to write the subset data to
 * @param rowIndex              row location of where to the put the data
 * @param colIndex              col location of where to put the data. If 1D arr, this is 0.
 * @param data                  memory address of the data being added
 * @return                      1 on success and 0 on failure 
 
 *****************************************************************************/
int writeSingleSubsetVarData_netcdf(
	string ncVarName,
	size_t rowIndex,		
	size_t colIndex,				
	const void* data		

) {
	try {
		NcVar variable = dataFile->getVar(ncVarName);
		vector<size_t> writeToLocation = {rowIndex, colIndex};
		variable.putVar(writeToLocation, data);

	} catch (NcException& e) {
		cout << "-E- Error writing data to the NcVar " << ncVarName << ". Err: " << e.what() << e.errorCode() << endl;
		return 0;
	}
	return 1;

}

/****************************************************************************
 * Same function as 'writeSingleSubsetVarData' except this function writes
 * a 2D array into the NetCDF.
 * 
 * @param ncVarName             name of the NcVar to write the subset data to
 * @param rowIndex              row location of where to the put the data
 * @param colIndex              col location of where to put the data.
 * @param subsetArrRowSize      row size of the array being added
 * @param subsetArrColSize      col size of the array being added
 * @param data                  memory address of the data being added
 * @return                      0 on success and 0 on failure 
 *****************************************************************************/
int write2DSubsetVarData_netcdf(
	string ncVarName,
	size_t rowIndex,		
	size_t colIndex,				
	size_t subsetArrRowSize,				
	size_t subsetArrColSize,			
	const void* data		

) {

	try {
		NcVar variable = dataFile->getVar(ncVarName);
		vector<size_t> writeToLocation = {rowIndex, colIndex};
                vector<size_t> newDataSize = {subsetArrRowSize, subsetArrColSize};
		variable.putVar(writeToLocation, newDataSize, data);
		
	} catch (NcException& e) {
		cout << "-E- Error writing data to the NcVar " << ncVarName << ". Err: " << e.what() << e.errorCode() << endl;
		return 1;
	}
	return 0;
}

/****************************************************************************
 * Same function as 'writeSingleSubsetVarData' except this function writes
 * a 3D array into the NetCDF.
 * 
 * @param ncVarName             name of the NcVar to write the subset data to
 * @param blockIndex            which block to put the 2D data into
 * @param rowIndex              row location of where to the put the data
 * @param colIndex              col location of where to put the data.
 * @param subsetArrBlockSize    how many blocks are in the arr being added
 * @param subsetArrRowSize      row size of the array being added
 * @param subsetArrColSize      col size of the array being added
 * @param data                  memory address of the data being added
 * @return                      0 on success and 1 on failure 
 *****************************************************************************/
int write3DSubsetVarData_netcdf(
	string ncVarName,
        size_t blockIndex,
	size_t rowIndex,		
	size_t colIndex,	
        size_t subsetArrBlockSize,			
	size_t subsetArrRowSize,				
	size_t subsetArrColSize,			
	const void* data		

) {

	try {
		NcVar variable = dataFile->getVar(ncVarName);
		vector<size_t> writeToLocation = {blockIndex, rowIndex, colIndex};
	        variable.putVar(writeToLocation, vector<size_t>{subsetArrBlockSize, subsetArrRowSize, subsetArrColSize}, data);

	} catch (NcException& e) {
		cout << "-E- Error writing data to the NcVar " << ncVarName << ". Err: " << e.what() << e.errorCode() << endl;
		return 1;
	}

	return 0;

}

/****************************************************************************
Write one scan line's worth of each of the passed in data to the appropriate
SDS in the currently open HDF file.
 *****************************************************************************/
extern "C" int WriteScanData_netcdf(
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
        fprintf(stderr, "in the current HDF file, %s .  Call ignored.\n", netcdfFile);
        return (PROGRAMMER_BOOBOO);
    }

	try {
		
		writeSingleSubsetVarData_netcdf("scan_time", scan, 0, &msec);
		write2DSubsetVarData_netcdf("eng_qual", scan, 0, 1, 4, eng_qual);
		write2DSubsetVarData_netcdf("s_flag", scan, 0, 1, 4, s_flags);
		write2DSubsetVarData_netcdf("s_satp", scan, 0, 1, 8, s_satp);
		write2DSubsetVarData_netcdf("s_zerop", scan, 0, 1, 8, s_zerop);
		writeSingleSubsetVarData_netcdf("slat", scan, 0, &slat);
		writeSingleSubsetVarData_netcdf("slon", scan, 0, &slon);
		writeSingleSubsetVarData_netcdf("clat", scan, 0, &clat);
		writeSingleSubsetVarData_netcdf("clon", scan, 0, &clon);
		writeSingleSubsetVarData_netcdf("elat", scan, 0, &elat);
		writeSingleSubsetVarData_netcdf("elon", scan, 0, &elon);
		writeSingleSubsetVarData_netcdf("csol_z", scan, 0, &csol_z);
		writeSingleSubsetVarData_netcdf("tilt", scan, 0, &tilt);
		write2DSubsetVarData_netcdf("sc_id", scan, 0, 1, 2, sc_id);
		write2DSubsetVarData_netcdf("sc_ttag", scan, 0, 1, 4, sc_ttag);
		write2DSubsetVarData_netcdf("sc_soh", scan, 0, 1, 775, sc_soh);
		write2DSubsetVarData_netcdf("inst_tlm", scan, 0, 1, 44, inst_tlm);
		write3DSubsetVarData_netcdf("l1a_data", scan, 0, 0, 1, numPixels, 8, l1a_data);
		write2DSubsetVarData_netcdf("start_syn", scan, 0, 1, 8, start_syn);
		write2DSubsetVarData_netcdf("stop_syn", scan, 0, 1, 8, stop_syn);
		write2DSubsetVarData_netcdf("dark_rest", scan, 0, 1, 8, dark_rest);
		write2DSubsetVarData_netcdf("gain", scan, 0, 1, 8, gain);
		write2DSubsetVarData_netcdf("tdi", scan, 0, 1, 8, tdi);
		write2DSubsetVarData_netcdf("inst_ana", scan, 0, 1, 40, inst_ana);
		write2DSubsetVarData_netcdf("inst_dis", scan, 0, 1, 32, inst_dis);
		write2DSubsetVarData_netcdf("sc_ana", scan, 0, 1, 40, sc_ana);
		write2DSubsetVarData_netcdf("sc_dis", scan, 0, 1, 40, sc_dis);
		write2DSubsetVarData_netcdf("scan_temp", scan, 0, 1, 8, scan_temp);
		writeSingleSubsetVarData_netcdf("side", scan, 0, &side);
		write2DSubsetVarData_netcdf("orb_vec", scan, 0, 1, 3, orb_vec);
		write2DSubsetVarData_netcdf("l_vert", scan, 0, 1, 3, l_vert);
		write2DSubsetVarData_netcdf("sun_ref", scan, 0, 1, 3, sun_ref);
		write2DSubsetVarData_netcdf("att_ang", scan, 0, 1, 3, att_ang);
		write3DSubsetVarData_netcdf("sen_mat", scan, 0, 0, 1, 3, 3, sen_mat);
		write2DSubsetVarData_netcdf("scan_ell", scan, 0, 1, 6, scan_ell);
		write2DSubsetVarData_netcdf("nflag", scan, 0, 1, 8, nflag);
                	
	} catch(NcException& e) {
		cout << "Error in WriteScanData_netcdf. Error: " << e.what() << endl;
		return 0;
	}

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
extern "C" int CloseL1aFile_netcdf(l1met *metrics) {

    // /* Write out the file metrics as global attributes. */
    // PTB(sd_setattr(ds_id.fid,
    //         "Gain 1 Saturated Pixels", DFNT_INT32, 8, metrics->gain1_satpix));
    // PTB(sd_setattr(ds_id.fid,
    //         "Gain 2 Saturated Pixels", DFNT_INT32, 8, metrics->gain2_satpix));
    // PTB(sd_setattr(ds_id.fid,
    //         "Gain 1 Non-Saturated Pixels", DFNT_INT32, 8, metrics->gain1_nonsatpix));
    // PTB(sd_setattr(ds_id.fid,
    //         "Gain 2 Non-Saturated Pixels", DFNT_INT32, 8, metrics->gain2_nonsatpix));
    // PTB(sd_setattr(ds_id.fid, "Zero Pixels", DFNT_INT32, 8, metrics->zeropix));
    // PTB(sd_setattr(ds_id.fid,
    //         "Mean Gain 1 Radiance", DFNT_FLOAT32, 8, metrics->gain1_mean_rad));
    // PTB(sd_setattr(ds_id.fid,
    //         "Mean Gain 2 Radiance", DFNT_FLOAT32, 8, metrics->gain2_mean_rad));

    AddCalData();

    // PTB(MakeVgroups());

    // if (SDend(ds_id.fid)) {
    //     fprintf(stderr, "-E- %s line %d: SDend(%d) failed for file, %s.\n",
    //             __FILE__, __LINE__, ds_id.fid, hdfFile);
    //     return (HDF_FUNCTION_ERROR);
    // }
    free(netcdfFile);
	dataFile->close();
    return (LIFE_IS_GOOD);
}


/****************************************************************************
Append calibration data to the HDF file.  The calibration data is read from
another HDF file pointed to by an environment variable.
 *****************************************************************************/
int AddCalData(void) {

    //   char          *envvar = "CAL_HDF_PATH";
    //   char          *calPath;
    //   int16         entry_year, entry_day, ref_year, ref_day, ref_minute;
    //   float32       temps[256][8], scan_mod[2][1285], mirror[2][8];
    //   float64       t_const[8], t_linear[8], t_quadratic[8];
    //   float32       cal_offs[8], counts[8][4][5], rads[8][4][5];

    //   calPath = getenv(envvar);
    //   if(calPath == NULL){
    //     fprintf(stderr,"-W- %s line %d: ",__FILE__,__LINE__);
    //     fprintf(stderr,"Environment variable, \"%s\", not set. ", envvar);
    //     fprintf(stderr,"Calibration data not appended to file, %s .\n", hdfFile);
    //     return(CALDATA_NOT_APPENDED);
    //   }

    //   if(
    //   get_cal(calPath, startYear, startDay, endDay,
    //           startMillisec, dataTypeString, tdi_global,
    //           &entry_year, &entry_day, &ref_year, &ref_day, &ref_minute,
    //           temps, scan_mod, mirror, t_const, t_linear, t_quadratic,
    //           cal_offs, counts, rads) < 0
    //   ){
    //     fprintf(stderr,"-W- %s line %d: ", __FILE__,__LINE__);
    //     fprintf(stderr,
    //     "get_cal(%s,%hd,%hd,%hd,%d,\"%s\",[%hd,%hd,%hd,%hd,%hd,%hd,%hd,%hd], ...) ",
    //     calPath,startYear,startDay,endDay,startMillisec,dataTypeString,
    //     tdi_global[0],tdi_global[1],tdi_global[2],tdi_global[3],
    //     tdi_global[4],tdi_global[5],tdi_global[6],tdi_global[7]);
    //     fprintf(stderr,"failed. ");
    //     fprintf(stderr,"Calibration data not appended to file, %s .\n", hdfFile);
    //     return(CALDATA_NOT_APPENDED);
    //   }


    // PTB(CreateSDS(ds_id.fid,
    //         "entry_year", /* short name */
    //         "Calibration entry year", /* long name */
    //         NULL, /* standard name */
    //         NULL, /* units */
    //         0, 0, /* valid range */
    //         1.0, 0.0, /* slope,offset */
    //         DFNT_INT16, /* HDF number type */
    //         1, /* rank */
    //         1, 1, 1, /* dimension sizes */
    //         "1", NULL, NULL /* dimension names */
    //         ));

    // PTB(CreateSDS(ds_id.fid,
    //         "entry_day", /* short name */
    //         "Calibration entry day-of-year", /* long name */
    //         NULL, /* standard name */
    //         NULL, /* units */
    //         0, 0, /* valid range */
    //         1.0, 0.0, /* slope,offset */
    //         DFNT_INT16, /* HDF number type */
    //         1, /* rank */
    //         1, 1, 1, /* dimension sizes */
    //         "1", NULL, NULL /* dimension names */
    //         ));

    // PTB(CreateSDS(ds_id.fid,
    //         "ref_year", /* short name */
    //         "Calibration reference year", /* long name */
    //         NULL, /* standard name */
    //         NULL, /* units */
    //         0, 0, /* valid range */
    //         1.0, 0.0, /* slope,offset */
    //         DFNT_INT16, /* HDF number type */
    //         1, /* rank */
    //         1, 1, 1, /* dimension sizes */
    //         "1", NULL, NULL /* dimension names */
    //         ));

    // PTB(CreateSDS(ds_id.fid,
    //         "ref_day", /* short name */
    //         "Calibration reference day-of-year", /* long name */
    //         NULL, /* standard name */
    //         NULL, /* units */
    //         0, 0, /* valid range */
    //         1.0, 0.0, /* slope,offset */
    //         DFNT_INT16, /* HDF number type */
    //         1, /* rank */
    //         1, 1, 1, /* dimension sizes */
    //         "1", NULL, NULL /* dimension names */
    //         ));

    // PTB(CreateSDS(ds_id.fid,
    //         "ref_minute", /* short name */
    //         "Calibration reference minute-of-day", /* long name */
    //         NULL, /* standard name */
    //         NULL, /* units */
    //         0, 0, /* valid range */
    //         1.0, 0.0, /* slope,offset */
    //         DFNT_INT16, /* HDF number type */
    //         1, /* rank */
    //         1, 1, 1, /* dimension sizes */
    //         "1", NULL, NULL /* dimension names */
    //         ));

    // PTB(CreateSDS(ds_id.fid,
    //         "mirror", /* short name */
    //         "Mirror-side correction factors", /* long name */
    //         NULL, /* standard name */
    //         NULL, /* units */
    //         0, 0, /* valid range */
    //         1.0, 0.0, /* slope,offset */
    //         DFNT_FLOAT32, /* HDF number type */
    //         2, /* rank */
    //         2, 8, 1, /* dimension sizes */
    //         "Number of Sides", "Number of Bands", NULL /* dimension names */
    //         ));

    // PTB(CreateSDS(ds_id.fid,
    //         "t_const", /* short name */
    //         "Time-dependent correction constant terms", /* long name */
    //         NULL, /* standard name */
    //         NULL, /* units */
    //         0, 0, /* valid range */
    //         1.0, 0.0, /* slope,offset */
    //         DFNT_FLOAT64, /* HDF number type */
    //         1, /* rank */
    //         8, 1, 1, /* dimension sizes */
    //         "Number of Bands", NULL, NULL /* dimension names */
    //         ));

    // PTB(CreateSDS(ds_id.fid,
    //         "t_linear", /* short name */
    //         "Time-dependent correction linear coefficients", /* long name */
    //         NULL, /* standard name */
    //         NULL, /* units */
    //         0, 0, /* valid range */
    //         1.0, 0.0, /* slope,offset */
    //         DFNT_FLOAT64, /* HDF number type */
    //         1, /* rank */
    //         8, 1, 1, /* dimension sizes */
    //         "Number of Bands", NULL, NULL /* dimension names */
    //         ));

    // PTB(CreateSDS(ds_id.fid,
    //         "t_quadratic", /* short name */
    //         "Time-dependent correction quadratic coefficients", /* long name */
    //         NULL, /* standard name */
    //         NULL, /* units */
    //         0, 0, /* valid range */
    //         1.0, 0.0, /* slope,offset */
    //         DFNT_FLOAT64, /* HDF number type */
    //         1, /* rank */
    //         8, 1, 1, /* dimension sizes */
    //         "Number of Bands", NULL, NULL /* dimension names */
    //         ));

    // PTB(CreateSDS(ds_id.fid,
    //         "cal_offs", /* short name */
    //         "Calibration system offsets", /* long name */
    //         NULL, /* standard name */
    //         NULL, /* units */
    //         0, 0, /* valid range */
    //         1.0, 0.0, /* slope,offset */
    //         DFNT_FLOAT32, /* HDF number type */
    //         1, /* rank */
    //         8, 1, 1, /* dimension sizes */
    //         "Number of Bands", NULL, NULL /* dimension names */
    //         ));

    // PTB(CreateSDS(ds_id.fid,
    //         "counts", /* short name */
    //         "Digital counts of calibration knees", /* long name */
    //         NULL, /* standard name */
    //         NULL, /* units */
    //         0, 1023, /* valid range */
    //         1.0, 0.0, /* slope,offset */
    //         DFNT_FLOAT32, /* HDF number type */
    //         3, /* rank */
    //         8, 4, 5, /* dimension sizes */
    //         "Number of Bands", "Number of Gains", "Number of Knees" /* dimension names */
    //         ));

    // PTB(CreateSDS(ds_id.fid,
    //         "rads", /* short name */
    //         "Radiances of calibration knees", /* long name */
    //         NULL, /* standard name */
    //         NULL, /* units */
    //         0, 0, /* valid range */
    //         1.0, 0.0, /* slope,offset */
    //         DFNT_FLOAT32, /* HDF number type */
    //         3, /* rank */
    //         8, 4, 5, /* dimension sizes */
    //         "Number of Bands", "Number of Gains", "Number of Knees" /* dimension names */
    //         ));

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
Construct a SeaWiFS level-1A filename from a time value and a data-type
value.  Return a pointer to the statically allocated filename string.
 *****************************************************************************/
extern "C" char * L1aFilename_netcdf(swl0ctl *l0ctl, double time, unsigned char dataType) {

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