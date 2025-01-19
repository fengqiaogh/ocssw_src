/*
 * sensorDefs.h
 *
 *  Created on: Oct 29, 2013
 *      Author: dshea
 */

#ifndef SENSOR_DEFS_H_
#define SENSOR_DEFS_H_

/* unique sensor IDs */
#define SEAWIFS        0

#define OCTS           2
#define AVHRR          3

#define CZCS           5
#define MODIST         6
#define MODISA         7
#define OCM1           8
#define OCM2           9
#define MERIS         10
#define VIIRSN        11 /* VIIRS NPP */
#define OCRVC         12
#define HICO          13
#define GOCI          14
#define OLIL8         15
#define AQUARIUS      16

#define AVIRIS        18
#define PRISM         19
#define OLCIS3A       20
#define SGLI          21
#define MSIS2A        22 /* MSI Sentinel-2A */
#define L5TM          23
#define L7ETMP        24
#define VIIRSJ1       25 /* VIIRS J1 */
#define MSIS2B        26 /* MSI Sentinel-2B */
#define HAWKEYE       27
#define MISR          28
#define OLCIS3B       29
#define OCI           30 

#define VIIRSJ2       32 /* VIIRS J2 */
#define OLIL9         33 /* OLI LANDSAT9 */
#define SPEXONE       34
#define HARP2         35
#define HARP          36
#define SENSOR_NUM    37

#define MAXWAVE_VIS  720
#define MINWAVE_IR  3000

/* subsensor IDs */
#define SEAWIFS_GAC    0
#define SEAWIFS_LAC    1
#define MODIS_TERRA    2
#define MODIS_AQUA     3
#define VIIRS_NPP      4
#define VIIRS_J1       5
#define MSI_S2A        6
#define MSI_S2B        7
#define OLCI_S3A       8
#define OLCI_S3B       9
#define VIIRS_J2       10
#define OLI_L8         11
#define OLI_L9         12
#define SUBSENSOR_NUM  13

/* instrument IDs */
#define INSTRUMENT_SEAWIFS    0

#define INSTRUMENT_OCTS       2
#define INSTRUMENT_AVHRR      3

#define INSTRUMENT_CZCS       5
#define INSTRUMENT_MODIS      6
#define INSTRUMENT_OCM        7
#define INSTRUMENT_OCM2       8
#define INSTRUMENT_MERIS      9
#define INSTRUMENT_VIIRS      10
#define INSTRUMENT_OCRVC      11
#define INSTRUMENT_HICO       12
#define INSTRUMENT_GOCI       13
#define INSTRUMENT_OLI        14
#define INSTRUMENT_AQUARIUS   15

#define INSTRUMENT_AVIRIS     17
#define INSTRUMENT_PRISM      18
#define INSTRUMENT_OLCI       19
#define INSTRUMENT_SGLI       20
#define INSTRUMENT_MSI        21
#define INSTRUMENT_L5TM       22
#define INSTRUMENT_L7ETMP     23
#define INSTRUMENT_HAWKEYE    24
#define INSTRUMENT_MISR       25
#define INSTRUMENT_OCI        26

#define INSTRUMENT_SPEXONE    28
#define INSTRUMENT_HARP2      29
#define INSTRUMENT_HARP       30
#define INSTRUMENT_NUM        31


#endif /* SENSOR_DEFS_H_ */
