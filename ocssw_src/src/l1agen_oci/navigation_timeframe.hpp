#ifndef _NAVIGATION_TIMEFRAME_
#define _NAVIGATION_TIMEFRAME_


// hold time for scan time and attitude time to determine valid navigation coverage
struct NavigationTimeFrame {

    // navigation time is the attitude time extracted from HKT file.
    double navigationStartTime = BAD_FLT;
    double navigationEndTime = BAD_FLT;

    // start and end times are updated as you extract the times from ancillary packets.
    double scanStartTime = BAD_FLT;
    double scanEndTime = BAD_FLT;


    bool isStartScanTimeSet() {
        return scanStartTime != BAD_FLT;
    }

    void reset() {
        navigationStartTime = BAD_FLT;
        navigationEndTime = BAD_FLT;
        scanStartTime = BAD_FLT;
        scanEndTime = BAD_FLT;
    }

    /**
     * @brief Checks first and last scan and attitude time for good nav coverage
     * @return true if first scan time > first attitude time or 
     *         or last scan time > last attitude time
     */
    bool isGoodNavigationCoverage() {    
        // for a valid navigation coverage, att time has to cover scan time
        // ------|----------------|------------------|------------------|---------
        // (1st att time)   (1st scan time)   (last scan time)   (last att time)
        return (navigationStartTime <= scanStartTime && 
            scanEndTime <= navigationEndTime);
    }
};

#endif