#include "dimension_shape.hpp"


// Constructor and Destructor takes no arguments
DimensionShape::DimensionShape() {}

DimensionShape::~DimensionShape() {}

void DimensionShape::incrementNumScansShape(size_t incrementBy) {
    numScans += incrementBy;
}

void DimensionShape::incrementNumMceScanShape(size_t incrementBy) {
    numMceScans += incrementBy;
}

void DimensionShape::incrementNumScaScanShape(size_t incrementBy) {
    numScaScans += incrementBy;
}

void DimensionShape::incrementAttRecordsShape(size_t incrementBy) {
    attRecords += incrementBy;
}

void DimensionShape::incrementOrbRecordsShape(size_t incrementBy) {
    orbRecords += incrementBy;
}

void DimensionShape::incrementTlmPacketsShape(size_t incrementBy) {
    tlmPackets += incrementBy;
}

void DimensionShape::incrementTiltSampleShape(size_t incrementBy) {
    tiltSamples += incrementBy;
}

void DimensionShape::incrementScienceScanNumShape() {
    scienceScanNum += 1;
}

// overloaded version to increment by the larger of the 2 increment values
void DimensionShape::incrementNumScansShape(size_t val1, size_t val2) {
    if (val1 < val2) {
        numScans += val2;
        return;
    }
    numScans += val1;
}

void DimensionShape::incrementNumMceScanShape(size_t val1, size_t val2) {
    if (val1 < val2) {
        numMceScans += val2;
        return;
    }
    numMceScans += val1;
}

void DimensionShape::incrementNumScaScanShape(size_t val1, size_t val2) {
    if (val1 < val2) {
        numScaScans += val2;
        return;
    }
    numScaScans += val1;
}
