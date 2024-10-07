
#include <cassert>
#include <stdexcept>

#include "scaledNcVar.hpp"

using namespace netCDF;
using namespace std;

const double TEST_VALUE = 20.1;  // To be put and expected after getting
const double TEST_BADVALUE = -1;
const double TEST_SCALE = 2.0;
const double TEST_OFFSET = 1.0;
const int TEST_COMPRESS_VALUE = (TEST_VALUE - TEST_OFFSET) / TEST_SCALE;  // expected raw data value
const double TEST_UNCOMPRESSED_VALUE = (TEST_COMPRESS_VALUE * TEST_SCALE) + TEST_OFFSET;
const int TEST_SIZE = 50;

void testGetVar(const vector<double> &data) {
    for (double datum : data) {
        if (datum != TEST_VALUE) {
            cerr << "-E- " << datum << " != " << TEST_VALUE << endl;
        }
    }
}

int testReading() {
    NcFile inFile("test.nc", NcFile::read);
    NcGroup group = inFile.getGroup("testGroup");
    vector<double> rawData(TEST_SIZE);

    // Scale and offset should have been ignored
    vector<double> dataScaledDouble(TEST_SIZE);
    ScaledNcVar scaledDouble = group.getVar("scaledDouble");
    scaledDouble.getVar(dataScaledDouble.data());
    for (double datum : dataScaledDouble) {
        if (datum != TEST_VALUE) {
            cerr << "-E- Scaled double: " << datum << " != " << TEST_VALUE << endl;
        }
    }

    // Scale and offset should NOT have been ignored, but ScaledNcVar::getVar() should undo it
    vector<double> dataScaledShort(TEST_SIZE);
    ScaledNcVar scaledShort = group.getVar("scaledShort");
    scaledShort.getVar(dataScaledShort.data());  // Undo scale and offset
    scaledShort.NcVar::getVar(rawData.data());   // Get what was actually written
    for (double datum : dataScaledShort) {
        if (datum != TEST_UNCOMPRESSED_VALUE) {
            cerr << "-E- Scaled short " << datum << " != " << TEST_UNCOMPRESSED_VALUE << endl;
        }
    }
    for (double datum : rawData)
        if (datum != TEST_COMPRESS_VALUE) {
            cerr << "-E- Scaled short raw data " << datum << " != " << TEST_COMPRESS_VALUE << endl;
        }

    // If ScaledNcVar has a user-defined "bad value", then, when it finds a fill value in a file, it should
    // interpret that fill value as that user-defined bad value
    vector<double> dataBadValue(TEST_SIZE);
    vector<double> rawDataBadValue(TEST_SIZE);
    ScaledNcVar badValue = group.getVar("badValue");
    badValue.assignBadValue(TEST_BADVALUE);
    badValue.NcVar::getVar(rawDataBadValue.data());  // Raw fill value
    badValue.getVar(dataBadValue.data());            // Fill value translated to badValue
    for (size_t i = 0; i < dataBadValue.size(); i++) {
        if (i % 2 == 0) {
            if (dataBadValue[i] != TEST_BADVALUE) {
                cerr << "-E- Bad value: " << dataBadValue[i] << " != " << TEST_BADVALUE << endl;
            }
        }
    }
    for (size_t i = 0; i < rawData.size(); i++) {
        if (i % 2 == 0) {
            if (rawDataBadValue[i] != BAD_FLT) {
                cerr << "-E- Bad value: " << rawDataBadValue[i] << " != " << BAD_FLT << endl;
            }
        }
    }

    // Test getting a range of data from a variable instead of the whole thing
    vector<double> dataVectored(10);
    ScaledNcVar vectored = group.getVar("vectored");
    vectored.getVar({2, 0}, {1, 10}, dataVectored.data());
    for (double value : dataVectored)
        if (value != TEST_VALUE)
            cerr << "-E- Vectored: " << value << " != " << TEST_VALUE << endl;

    return EXIT_SUCCESS;
}

int testWriting() {
    int status = EXIT_SUCCESS;
    vector<double> pretendDataDouble(10);

    for (size_t i = 0; i < 10; i++) {
        pretendDataDouble[i] = i + 1;
    }

    NcFile outFile("test.nc", NcFile::replace);
    auto dimX = outFile.addDim("x", 10);
    auto dimY = outFile.addDim("y", 5);
    NcGroup testGroup = outFile.addGroup("testGroup");
    vector<double> doubleData(TEST_SIZE, TEST_VALUE);
    vector<double> smokeArray(TEST_SIZE, .012);
    smokeArray[1] = BAD_FLT;
    vector<double> compressable(TEST_SIZE, TEST_COMPRESS_VALUE);
    vector<double> dataWithBadValue(TEST_SIZE, TEST_VALUE);
    for (size_t i = 0; i < dataWithBadValue.size(); i++) {
        if (i % 2 == 0)
            dataWithBadValue[i] = TEST_BADVALUE;
    }

    // Should ignore scale/offset
    ScaledNcVar scaledDouble = testGroup.addVar("scaledDouble", ncDouble, {dimY, dimX});
    try {
        scaledDouble.setScaleFactors(TEST_SCALE, TEST_OFFSET);
    } catch (const invalid_argument &e) {
        cout << "Caught exception - " << e.what() << endl;
    }
    scaledDouble.putVar(doubleData.data());

    // Should NOT ignore scale/offset
    ScaledNcVar scaledShort = testGroup.addVar("scaledShort", ncShort, {dimY, dimX});
    scaledShort.setScaleFactors(TEST_SCALE, TEST_OFFSET, 1);
    scaledShort.putVar(doubleData.data());

    // Populating from product.xml
    ScaledNcVar ndvi = newScaledNcVar(testGroup, "ndvi", {dimY, dimX});
    ndvi.putVar(doubleData.data());

    // Also populating from product.xml, this product has a scale and offset. Should NOT ignore scale/offset
    ScaledNcVar smoke = newScaledNcVar(testGroup, "smoke", {dimY, dimX});
    smoke.putVar(smokeArray.data());

    // Putting a variable in specified spots (like a line out of a whole file)
    ScaledNcVar vectored = testGroup.addVar("vectored", ncDouble, {dimY, dimX});
    vectored.putVar({2, 0}, {1, 10}, doubleData.data());

    // When a ScaledNcVar finds a user-defined "bad value", it should write the assigned (or default) fill
    // value instead
    ScaledNcVar badValue = testGroup.addVar("badValue", ncDouble, {dimY, dimX});
    badValue.assignFillValue(BAD_FLT);
    badValue.assignBadValue(TEST_BADVALUE);
    badValue.putVar(dataWithBadValue.data());

    // Expect an exception if a fill value that isn't representable is passed to assignFillValue
    ScaledNcVar outOfRange = testGroup.addVar("outOfRange", ncByte, {dimY, dimX});
    try {
        outOfRange.assignFillValue(257);
    } catch (const out_of_range &e) {
        cout << "Caught exception - " << e.what() << endl;
    }

    outFile.close();

    return status;
}

int main() {
    return testWriting() + testReading();  // Nonzero if any fail
}