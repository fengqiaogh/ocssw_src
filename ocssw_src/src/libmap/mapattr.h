#ifndef MAPATTR_H /* avoid re-inclusion */
#define MAPATTR_H

/*  Global (file-level) Attributes  */
/*-----------------------------------------------------------------------------
     Mission and Documentation
-----------------------------------------------------------------------------*/
#define L3M_PNAME     "Product Name"
#define L3M_TITLE      "Title"
#define L3M_TITLE_VAL    "SeaWiFS Level-3 Standard Mapped Image"
#define L3M_DCENTER  "Data Center"
#define L3M_DCENTER_VAL    "NASA/GSFC SeaWiFS Data Processing Center"
#define L3M_STATION  "Station Name"
#define L3M_STATION_VAL  "Wallops Flight Facility"
#define L3M_STLAT  "Station Latitude"
#define L3M_STLAT_VAL   37.9272
#define L3M_STLON  "Station Longitude"
#define L3M_STLON_VAL  -75.4753f 
#define L3M_MISSION  "Mission"
#define L3M_MISSION_VAL    "SeaStar SeaWiFS"
#define L3M_MSNCHAR    "Mission Characteristics"
#define L3M_MSNCHAR_VAL  "Nominal orbit: inclination = 98.2 (Sun-synchronous); node = 12 noon local (descending); eccentricity = <0.002; altitude = 705 km; ground speed = 6.75 km/sec"
#define L3M_SENSOR  "Sensor"
#define L3M_SENSOR_NAME  "Sensor Name"
#define L3M_SENSOR_VAL   "Sea-viewing Wide Field-of-view Sensor (SeaWiFS)"
#define L3M_SNSCHAR   "Sensor Characteristics"

#define L3M_SNSCHAR_VAL        "Number of bands = 8; number of active bands = 8; wavelengths per band (nm) = 412, 443, 490, 510, 555, 670, 765, 865; bits per pixel = 10; instantaneous field-of-view = 1.5835 mrad; pixels per scan = 1285; scan rate = 6/sec; sample rate = 7710/sec"

#define L3M_PRODTYPE  "Product Type"
#define L3M_PVERSION  "Processing Version"
#define L3M_REPLACE  "Replacement Flag"
#define L3M_SOFTID  "Software ID"
#define L3M_SOFTNM  "Software Name"
#define L3M_SOFTVER  "Software Version"
#define L3M_PTIME  "Processing Time"
#define L3M_INFILES  "Input Files"
#define L3M_PROCCON  "Processing Control"
#define L3M_PROCLOG  "Processing Log"
#define L3M_INPARMS  "Input Parameters"
#define L3M_FLAG_NAMES  "L2 Flag Names"
#define L3M_FLAG_USE  "L2 Flag Usage"
#define L3M_ENG_Q_USE  "L2 Engineering Quality Usage"
/*-----------------------------------------------------------------------------
     Data Time
-----------------------------------------------------------------------------*/
#define L3M_PSYEAR  "Period Start Year"
#define L3M_PSDAY  "Period Start Day"
#define L3M_PEYEAR  "Period End Year"
#define L3M_PEDAY  "Period End Day"
#define L3M_STIME  "Start Time"
#define L3M_ETIME  "End Time"
#define L3M_SYEAR  "Start Year"
#define L3M_SDAY  "Start Day"
#define L3M_SMSEC  "Start Millisec"
#define L3M_EYEAR  "End Year"
#define L3M_EDAY  "End Day"
#define L3M_EMSEC  "End Millisec"
#define L3M_ORBIT  "Orbit"
#define L3M_SORBIT  "Start Orbit"
#define L3M_EORBIT  "End Orbit"
/*-----------------------------------------------------------------------------
     Scene Coordinates
-----------------------------------------------------------------------------*/
#define L3M_MAPPROJ  "Map Projection"
#define L3M_MAPPROJ_VAL  "Equidistant Cylindrical"
#define L3M_LATUNITS  "Latitude Units"
#define L3M_LATUNITS_VAL "degrees North"
#define L3M_LONUNITS  "Longitude Units"
#define L3M_LONUNITS_VAL "degrees East"
#define L3M_NLAT  "Northernmost Latitude"
#define L3M_NLAT_VAL  90.0
#define L3M_SLAT  "Southernmost Latitude"
#define L3M_SLAT_VAL  -90.0
#define L3M_WLON  "Westernmost Longitude"
#define L3M_WLON_VAL  -180.0
#define L3M_ELON  "Easternmost Longitude"
#define L3M_ELON_VAL  180.0
#define L3M_LAT_STEP  "Latitude Step"
#define L3M_LON_STEP  "Longitude Step"
#define L3M_SWLAT  "SW Point Latitude"
#define L3M_SWLON  "SW Point Longitude"
/*-----------------------------------------------------------------------------
     Data Description 
-----------------------------------------------------------------------------*/
#define SDS_NAME  "l3m_data"
#define LONGNAME  "long_name"
#define LONGNAME_VAL  "Level 3 Standard Mapped Image Data"

#define L3M_DATABINS  "Data Bins"
#define L3M_NROWS  "Number of Lines"
#define L3M_NCOLS  "Number of Columns"
#define L3M_MEASURE  "Measure"
#define L3M_PARAMETER  "Parameter"
#define L3M_UNITS  "Units"
#define L3M_SCALING  "Scaling"
#define L3M_LOG_SCALE  "logarithmic"
#define L3M_LINEAR_SCALE  "linear"
#define L3M_SC_EQN  "Scaling Equation"
#define L3M_LOG_EQN   "Base**((Slope*l3m_data) + Intercept) = Parameter value"
#define L3M_LINEAR_EQN  "(Slope*l3m_data) + Intercept = Parameter value"
#define L3M_BASE  "Base"
#define L3M_SLOPE  "Slope"
#define L3M_INTERCEPT  "Intercept"
#define L3M_MIN   "Data Minimum"
#define L3M_MAX   "Data Maximum"
/*---------------------------------------------------------------------------*/

#endif /* MAPATTR_H */

