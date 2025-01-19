#include <sensorInfo.h>

#include <stddef.h>
#include <string.h>
#include <stdlib.h>

#include <genutils.h>

// sensor name indexed by sensorId
// NOTE - some of the sensors have been deleted, but we have to keep
//        the placeholder to make the index match
static const char* sensorName[] = {
    "SeaWiFS",  // 0
    "MOS(Delete)",      // 1 deleted
    "OCTS",     // 2
    "AVHRR",    // 3
    "OSMI(Delete)",     // 4 deleted
    "CZCS",     // 5
    "MODIST",   // 6
    "MODISA",   // 7
    "OCM1",     // 8
    "OCM2",     // 9
    "MERIS",    // 10
    "VIIRSN",   // 11
    "OCRVC",    // 12
    "HICO",     // 13
    "GOCI",     // 14
    "OLIL8",    // 15
    "Aquarius", // 16
    "OCIA(Delete)",     // 17 deleted
    "AVIRIS",   // 18
    "PRISM",    // 19
    "OLCIS3A",  // 20
    "SGLI",     // 21
    "MSIS2A",   // 22
    "L5TM",     // 23
    "L7ETMP",   // 24
    "VIIRSJ1",  // 25
    "MSIS2B",   // 26
    "HAWKEYE",  // 27
    "MISR",     // 28
    "OLCIS3B",  // 29
    "OCI",      // 30
    "OCIS(Delete)",     // 31 deleted
    "VIIRSJ2",  // 32
    "OLIL9",    // 33
    "SPEXONE",  // 34
    "HARP2",    // 35
    "HARP"      // 36
};

// instrument name indexed by sensorId
static const char *instrumentName[] = {
    "SeaWiFS",  // 0
    "MOS(Delete)",      // 1 deleted
    "OCTS",     // 2
    "AVHRR",    // 3
    "OSMI(Delete)",     // 4 deleted
    "CZCS",     // 5
    "MODIS",    // 6
    "MODIS",    // 7
    "OCM",      // 8
    "OCM-2",    // 9
    "MERIS",    // 10
    "VIIRS",    // 11
    "OCRVC",    // 12
    "HICO",     // 13
    "GOCI",     // 14
    "OLI",      // 15
    "Aquarius", // 16
    "OCIA(Delete)",     // 17 deleted
    "AVIRIS",   // 18
    "PRISM",    // 19
    "OLCI",     // 20
    "SGLI",     // 21
    "MSI",      // 22
    "L5TM",     // 23
    "L7ETMP",   // 24
    "VIIRS",    // 25
    "MSI",      // 26
    "HAWKEYE",  // 27
    "MISR",     // 28
    "OLCI",     // 29
    "OCI",      // 30
    "OCIS(Delete)",     // 31 deleted
    "VIIRS",    // 32
    "OLI",      // 33
    "SPEXONE",  // 34
    "HARP2",    // 35
    "HARP"      // 36
};

// platform name indexed by sensorId
static const char *platformName[] = {
    "Orbview-2",    // 0
    "IRS-P3(Delete)",       // 1
    "ADEOS",        // 2
    "AVHRR",        // 3
    "KOMPSAT(Delete)",      // 4
    "Nimbus-7",     // 5
    "Terra",        // 6
    "Aqua",         // 7
    "IRS-P4",       // 8
    "Oceansat-2",   // 9
    "Envisat",      // 10
    "Suomi-NPP",    // 11
    "OCRVC",        // 12
    "ISS",          // 13
    "COMS",         // 14
    "Landsat-8",    // 15
    "SAC-D",        // 16
    "PACE(Delete)",         // 17
    "AVIRIS",       // 18
    "PRISM",        // 19
    "Sentinel-3A",  // 20
    "GCOM_C",       // 21
    "Sentinel-2A",  // 22
    "L5TM",         // 23
    "L7ETMP",       // 24
    "JPSS-1",       // 25
    "Sentinel-2B",  // 26
    "Seahawk1",     // 27
    "Terra",        // 28
    "Sentinel-3B",  // 29
    "PACE",         // 30
    "PACE(Delete)",         // 31
    "JPSS-2",       // 32
    "Landsat-9",    // 33
    "PACE",         // 34
    "PACE",         // 35
    "Air-HARP"      // 36
};

// sensor directory indexed by sensorId
static const char *sensorDir[] = {
    "seawifs",  // 0
    "mos(Delete)",      // 1 deleted
    "octs",     // 2
    "avhrr",    // 3
    "osmi(Delete)",     // 4 deleted
    "czcs",     // 5
    "modis",    // 6
    "modis",    // 7
    "ocm1",     // 8
    "ocm2",     // 9
    "meris",    // 10
    "viirs",    // 11
    "ocrvc",    // 12
    "hico",     // 13
    "goci",     // 14
    "oli",      // 15
    "aquarius", // 16
    "ocia(Delete)",     // 17 deleted
    "aviris",   // 18
    "prism",    // 19
    "olci",     // 20
    "sgli",     // 21
    "msi",      // 22
    "l5tm",     // 23
    "l7etmp",   // 24
    "viirs",    // 25
    "msi",      // 26
    "hawkeye",  // 27
    "misr",     // 28
    "olci",     // 29
    "oci",      // 30
    "ocis(Delete)",     // 31 deleted
    "viirs",    // 32
    "oli",      // 33
    "spexone",  // 34
    "harp2",    // 35
    "harp"      // 36
};

// subsensor directory indexed by subsensorId
static const char *subsensorDir[] = {
    "gac",      // 0
    "lac",      // 1
    "terra",    // 2
    "aqua",     // 3
    "npp",      // 4
    "j1",       // 5
    "s2a",      // 6
    "s2b",      // 7
    "s3a",      // 8
    "s3b",      // 9
    "j2",       // 10
    "l8",       // 11
    "l9"        // 12
};

// instrument ID indexed by sensorId
static const int instrumentId[] = {
    INSTRUMENT_SEAWIFS, // 0
    1,
    INSTRUMENT_OCTS,    // 2
    INSTRUMENT_AVHRR,   // 3
    1,
    INSTRUMENT_CZCS,    // 5
    INSTRUMENT_MODIS,   // 6
    INSTRUMENT_MODIS,   // 7
    INSTRUMENT_OCM,     // 8
    INSTRUMENT_OCM2,    // 9
    INSTRUMENT_MERIS,   // 10
    INSTRUMENT_VIIRS,   // 11
    INSTRUMENT_OCRVC,   // 12
    INSTRUMENT_HICO,    // 13
    INSTRUMENT_GOCI,    // 14
    INSTRUMENT_OLI,     // 15
    INSTRUMENT_AQUARIUS,// 16
    1,
    INSTRUMENT_AVIRIS,  // 18
    INSTRUMENT_PRISM,   // 19
    INSTRUMENT_OLCI,    // 20
    INSTRUMENT_SGLI,    // 21
    INSTRUMENT_MSI,     // 22
    INSTRUMENT_L5TM,    // 23
    INSTRUMENT_L7ETMP,  // 24
    INSTRUMENT_VIIRS,   // 25
    INSTRUMENT_MSI,     // 26
    INSTRUMENT_HAWKEYE, // 27
    INSTRUMENT_MISR,    // 28
    INSTRUMENT_OLCI,    // 29
    INSTRUMENT_OCI,     // 30
    1,
    INSTRUMENT_VIIRS,   // 32
    INSTRUMENT_OLI,     // 33
    INSTRUMENT_SPEXONE, // 34
    INSTRUMENT_HARP2,   // 35
    INSTRUMENT_HARP     // 36
};

// instrument name indexed by instrumentId
static const char *instrumentNameByInstrumentId[] = {
    "SeaWiFS",  // 0
    "MOS(Delete)",      // 1 deleted
    "OCTS",     // 2
    "AVHRR",    // 3
    "OSMI(Delete)",     // 4 deleted
    "CZCS",     // 5
    "MODIS",    // 6
    "OCM",      // 7
    "OCM-2",    // 8
    "MERIS",    // 9
    "VIIRS",    // 10
    "OCRVC",    // 11
    "HICO",     // 12
    "GOCI",     // 13
    "OLI",      // 14
    "Aquarius", // 15
    "OCIA(Delete)",     // 16 deleted
    "AVIRIS",   // 17
    "PRISM",    // 18
    "OLCI",     // 19
    "SGLI",     // 20
    "MSI",      // 21
    "L5TM",     // 22
    "L7ETMP",   // 23
    "HAWKEYE",  // 24
    "MISR",     // 25
    "OCI",      // 26
    "OCIS(Delete)",     // 27 deleted
    "SPEXONE",  // 28
    "HARP2",    // 29
    "HARP"      // 30
};


/**
 * Get the name of the sensor for this sensorId.
 * 
 * @param sensorId sensor identifier to look up.
 * @return name of sensor or NULL if not found.
 */
const char* sensorId2SensorName(int sensorId) {
    if(sensorId < 0)
        return NULL;
    if(sensorId >= SENSOR_NUM)
        return NULL;
    return sensorName[sensorId];
}

/**
 * Get the name of the instrument for this sensorId.
 * 
 * @param sensorId sensor identifier to look up.
 * @return name of the instrument or NULL if not found.
 */
const char* sensorId2InstrumentName(int sensorId) {
    if(sensorId < 0)
        return NULL;
    if(sensorId >= SENSOR_NUM)
        return NULL;
    return instrumentName[sensorId];
}

/**
 * Get the name of the platform for this sensorId.
 * 
 * @param sensorId sensor identifier to look up.
 * @return name of the platform or NULL if not found.
 */
const char* sensorId2PlatformName(int sensorId) {
    if(sensorId < 0)
        return NULL;
    if(sensorId >= SENSOR_NUM)
        return NULL;
    return platformName[sensorId];
}

/**
 * Get the name of the sensor directory for this sensorId.
 * 
 * @param sensorId sensor identifier to look up.
 * @return sensor directory or NULL if not found.
 */
const char* sensorId2SensorDir(int sensorId) {
    if(sensorId < 0)
        return NULL;
    if(sensorId >= SENSOR_NUM)
        return NULL;
    return sensorDir[sensorId];
}

/**
 * Get the name of the subsensor directory for this subsensorId.
 * 
 * @param subsensorId subsensor identifier to look up.
 * @return subsensor directory or NULL if not found.
 */
const char* subsensorId2SubsensorDir(int subsensorId) {
    if(subsensorId < 0)
        return NULL;
    if(subsensorId >= SENSOR_NUM)
        return NULL;
    return subsensorDir[subsensorId];
}

/**
 * Get the instrument ID for this sensor ID
 * 
 * @param sensorId sensor identifier to look up.
 * @return int instrument ID for this sensor.  -1 if invalid sensorId.
 */
int sensorId2InstrumentId(int sensorId) {
    if(sensorId < 0)
        return -1;
    if(sensorId >= SENSOR_NUM)
        return -1;
    return instrumentId[sensorId];
}

/**
 * get the name of the instrument for this instrument ID
 * 
 * @param instrumentId instrument identifier to look up
 * @return const char* instrument name or NULL if not found
 */
const char* instrumentId2InstrumentName(int instrumentId) {
    if(instrumentId < 0)
        return NULL;
    if(instrumentId >= INSTRUMENT_NUM)
        return NULL;
    return instrumentNameByInstrumentId[instrumentId];
}

/**
 * lookup the ID for a sensor name
 *
 * @param sensorName name of the sensor to lookup
 * @return sensor ID for the given sensor, -1 if not found.
 */
int sensorName2SensorId(const char* name) {
    int i;

    if(name == NULL)
        return -1;
    
    // convert the number if an int was sent in
    if (isValidInt(name)) {
        i = atoi(name);
        if (i >= 0 && i < SENSOR_NUM) {
            return i;
        }
    }

    if (strcasecmp(name, "hmodisa") == 0)
        return MODISA;
    if (strcasecmp(name, "hmodist") == 0)
        return MODIST;

    for (i = 0; i < SENSOR_NUM; i++) {
        if (strcasecmp(sensorName[i], name) == 0)
            return i;
    }

    return -1;
}

/**
 * lookup the sensorID for a given instrument and platform
 *
 * @param instrument instrument to use for sensorID look up
 * @param platform platform to use for sensorID lookup
 * @return the matching sensorID, or -1 if not found
 */
int instrumentPlatform2SensorId(const char* instrument, const char* platform) {
    int i;

    if(instrument == NULL  || platform == NULL)
        return -1;
    
    // if multiple platforms, look up the first one
    char *tmpInstrument = strdup(instrument);
    char* ptr = strchr(tmpInstrument, ',');
    if(ptr) *ptr = '\0';

    char *tmpPlatform = strdup(platform);
    ptr = strchr(tmpPlatform, ',');
    if(ptr) *ptr = '\0';

    for (i = 0; i < SENSOR_NUM; i++) {
        if (strcasecmp(instrumentName[i], tmpInstrument) == 0)
            if (strcasecmp(platformName[i], tmpPlatform) == 0) {
                free(tmpInstrument);
                free(tmpPlatform);
                return i;
            }
    }
    free(tmpInstrument);
    free(tmpPlatform);
    return -1;
}

/**
 * return the subsensorId for the given sensorId
 * @param sensorId sensorId to lookup
 * @return the subsensorId, or -1 if not found
 */
int sensorId2SubsensorId(int sensorId) {
    switch(sensorId) {
    case MODIST:
        return MODIS_TERRA;
    case MODISA:
        return MODIS_AQUA;
    case VIIRSN:
        return VIIRS_NPP;
    case VIIRSJ1:
        return VIIRS_J1;
    case MSIS2A:
        return MSI_S2A;
    case MSIS2B:
        return MSI_S2B;
    case OLCIS3A:
        return OLCI_S3A;
    case OLCIS3B:
        return OLCI_S3B;
    case VIIRSJ2:
        return VIIRS_J2;
    case OLIL8:
        return OLI_L8;
    case OLIL9:
        return OLI_L9;
    default:
        return -1;
    }
}

/**
 * fine the sensor name given the instrument and platform
 * @param instrument name of instrument
 * @param platform name of platform
 * @return sensor name
 */
const char* instrumentPlatform2SensorName(const char* instrument, const char* platform) {
    return sensorId2SensorName(instrumentPlatform2SensorId(instrument, platform));
}
