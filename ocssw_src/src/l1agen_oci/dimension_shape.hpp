#ifndef _DIMENSION_SHAPE_HPP_
#define _DIMENSION_SHAPE_HPP_

#include <cstddef> // for size_t


// tracks UNLIMITED NetCDF dimensions for an L1A file. When appending multiple times for the same
// variable, tracking the current shape of the dimension so you know the index to start appending
// the next set of data to.


class DimensionShape {

    public:
        DimensionShape();
        ~DimensionShape();
        size_t numScans = 0;
        size_t numMceScans = 0;
        size_t numScaScans = 0;
        size_t attRecords = 0;      // attitude
        size_t orbRecords = 0;      // orbital
        size_t tlmPackets = 0;      // telemetry packets
        size_t tiltSamples = 0;

        // science_data is written line by line and not in batches (multiple lines).
        // this variable tracks the current line so that if there is a datatype change
        // when the data changes back, it knows where to start writing for the next
        // science data
        size_t scienceScanNum = 0; 

        void incrementNumScansShape(size_t incrementBy);

        void incrementNumMceScanShape(size_t incrementBy);

        void incrementNumScaScanShape(size_t incrementBy);

        void incrementAttRecordsShape(size_t incrementBy);

        void incrementOrbRecordsShape(size_t incrementBy);

        void incrementTlmPacketsShape(size_t incrementBy);

        void incrementTiltSampleShape(size_t incrementBy);

        // always by 1 because science data is written line by line
        void incrementScienceScanNumShape();

        // overloaded version to take in 2 increment numbers, but take the larger one
        void incrementNumScansShape(size_t val1, size_t val2);

        void incrementNumMceScanShape(size_t val1, size_t val2);

        void incrementNumScaScanShape(size_t val1, size_t val2);
    
};

#endif