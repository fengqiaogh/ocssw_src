#include <memory>
#include <cmath>
#include <functional>
#include <iomanip>

#include "scaled_nc_var.hpp"
#include "scaled_nc_group.hpp"
#include "scaled_nc_file.hpp"

using namespace netCDF;
using namespace std;

typedef pair<bool, string> TestsResult;

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

TestsResult testReading() {
    TestsResult result = {true, ""};
    NcFile inFile("test.nc", NcFile::read);
    NcGroup group = inFile.getGroup("testGroup");
    vector<double> rawData(TEST_SIZE);

    // Scale and offset should have been ignored
    vector<double> dataScaledDouble(TEST_SIZE);
    ScaledNcVar scaledDouble = group.getVar("scaledDouble");
    scaledDouble.getVar(dataScaledDouble.data());
    for (double datum : dataScaledDouble) {
        if (datum != TEST_VALUE) {
            result.first = false;
            result.second += "Scaled double: " + to_string(datum) + " != " + to_string(TEST_VALUE) + "\n";
        }
    }

    // Scale and offset should NOT have been ignored, but ScaledNcVar::getVar() should undo it
    vector<double> dataScaledShort(TEST_SIZE);
    ScaledNcVar scaledShort = group.getVar("scaledShort");
    scaledShort.getVar(dataScaledShort.data());  // Undo scale and offset
    scaledShort.NcVar::getVar(rawData.data());   // Get what was actually written
    for (double datum : dataScaledShort) {
        if (datum != TEST_UNCOMPRESSED_VALUE) {
            result.first = false;
            result.second +=
                "Scaled short " + to_string(datum) + " != " + to_string(TEST_UNCOMPRESSED_VALUE) + "\n";
        }
    }
    for (double datum : rawData)
        if (datum != TEST_COMPRESS_VALUE) {
            result.first = false;
            result.second +=
                "Scaled short raw data " + to_string(datum) + " != " + to_string(TEST_COMPRESS_VALUE) + "\n";
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
                result.first = false;
                result.second +=
                    "Bad value: " + to_string(dataBadValue[i]) + " != " + to_string(TEST_BADVALUE) + "\n";
            }
        }
    }
    for (size_t i = 0; i < rawData.size(); i++) {
        if (i % 2 == 0) {
            if (rawDataBadValue[i] != BAD_FLT) {
                result.first = false;
                result.second +=
                    "Bad value: " + to_string(rawDataBadValue[i]) + " != " + to_string(BAD_FLT) + "\n";
            }
        }
    }

    // Test getting a range of data from a variable instead of the whole thing
    vector<double> dataVectored(10);
    ScaledNcVar vectored = group.getVar("vectored");
    vectored.getVar({2, 0}, {1, 10}, dataVectored.data());
    for (double value : dataVectored)
        if (value != TEST_VALUE) {
            result.first = false;
            result.second += "Vectored: " + to_string(value) + " != " + to_string(TEST_VALUE) + "\n";
        }
    ScaledNcVar nanvar = group.getVar("nanvar");
    vector<double> nanVec(nanvar.getDimsSize());
    nanvar.getVar(nanVec.data());
    for (double value : nanVec)
        if (isnan(value)) {
            result.first = false;
            result.second += "NaN value found in nanvar\n";
        }

    return result;
}

TestsResult testWriting() {
    TestsResult result = {true, ""};
    vector<double> pretendDataDouble(10);

    for (size_t i = 0; i < 10; i++) {
        pretendDataDouble[i] = i + 1;
    }

    try {
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
            result.second += "Test 1: Caught expected exception - " + string(e.what()) + "\n";
        }
        scaledDouble.putVar(doubleData.data());

        // Should NOT ignore scale/offset
        ScaledNcVar scaledShort = testGroup.addVar("scaledShort", ncShort, {dimY, dimX});
        scaledShort.setScaleFactors(TEST_SCALE, TEST_OFFSET, 1);
        scaledShort.putVar(doubleData.data());

        // Populating from product.xml
        ScaledNcVar ndvi = newScaledNcVar(testGroup, "ndvi", {dimY, dimX});
        ndvi.putVar(doubleData.data());

        // Also populating from product.xml, this product has a scale and offset. Should NOT ignore
        // scale/offset
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

        ScaledNcVar defaultSNcVar;
        netCDF::NcVar defaultNcVar = testGroup.addVar("defaultVar", ncByte, {dimY, dimX});
        ScaledNcVar copy;
        copy = defaultNcVar;

        netCDF::NcVar fromGroup = testGroup.getVar("badValue");  // Could be any var
        ScaledNcVar copiedFromGroup;
        copiedFromGroup = fromGroup;

        // Expect an exception if a fill value that isn't representable is passed to assignFillValue
        ScaledNcVar outOfRange = testGroup.addVar("outOfRange", ncByte, {dimY, dimX});
        try {
            outOfRange.assignFillValue(257);
            result.first = false;
            result.second += "Test 2: Expected out_of_range exception, but none was thrown\n";
        } catch (const out_of_range &e) {
            result.second += "Test 2: Caught expected exception - " + string(e.what()) + "\n";
        }

        // What does NaN do?
        ScaledNcVar nanVar = testGroup.addVar("nanvar", ncDouble, {dimY, dimX});
        double zero = 1 - 1;  // Quieting a warning
        doubleData[0] = 1 / zero;
        nanVar.putVar(doubleData.data());

        outFile.close();
    } catch (const exception &e) {
        result.first = false;
        result.second += "Unexpected exception: " + string(e.what()) + "\n";
    }

    return result;
}
TestsResult testScaledNcGroupMultimapGetVars() {
    TestsResult result = {true, ""};

    string fileName = "scaledNcTest.nc";
    ScaledNcFile file(fileName, NcFile::replace);
    ScaledNcGroup group = file.addGroup("groupTest");
    group.addVar("firstVar", NC_FLOAT);
    group.addVar("secondVar", NC_SHORT);
    multimap<string, ScaledNcVar> varMap = group.getVars();
    // Test 1: Check if the correct number of variables were retrieved
    if (varMap.size() != 2) {
        result.first = false;
        result.second += "   Test 1: Expected 2 variables, but got " + to_string(varMap.size()) + "\n";
    }

    // Test 2: Verify that the retrieved variable names match the expected ones
    set<string> expectedVarNames = {"firstVar", "secondVar"};
    set<string> actualVarNames;
    for (const auto &pair : varMap) {
        actualVarNames.insert(pair.first);
    }
    if (expectedVarNames != actualVarNames) {
        result.first = false;
        result.second += "   Test 2: Variable names do not match expected names\n";
    }

    // Test 3 and 4: Ensure each variable has the correct data type
    for (const auto &pair : varMap) {
        const string &name = pair.first;
        const ScaledNcVar &var = pair.second;
        if (name == "firstVar" && var.getType() != NC_FLOAT) {
            result.first = false;
            result.second += "   Test 3: firstVar is not of type NC_FLOAT\n";
        } else if (name == "secondVar" && var.getType() != NC_SHORT) {
            result.first = false;
            result.second += "   Test 4: secondVar is not of type NC_SHORT\n";
        }
    }

    if (remove(fileName.c_str()) != 0) {
        cerr << "-E- Couldn't delete " << fileName << endl;
    }

    return result;
}

TestsResult testScaledNcGroupSetGetVars() {
    TestsResult result = {true, ""};

    string fileName = "scaledNcTest.nc";
    ScaledNcFile file(fileName, NcFile::replace);
    ScaledNcGroup group = file.addGroup("groupTest");
    group.addVar("firstVar", NC_FLOAT);
    group.addVar("secondVar", NC_SHORT);

    set<ScaledNcVar> varSet = group.getVars("firstVar");

    // Test 1: Check if the set contains exactly one variable
    if (varSet.size() != 1) {
        result.first = false;
        result.second += "   Test 1: Expected 1 variable, but got " + to_string(varSet.size()) + "\n";
    }

    if (!varSet.empty()) {
        const ScaledNcVar &var = *varSet.begin();

        // Test 2: Verify the name of the retrieved variable
        if (var.getName() != "firstVar") {
            result.first = false;
            result.second +=
                "   Test 2: Expected variable name 'firstVar', but got '" + var.getName() + "'\n";
        }

        // Test 3: Check if the variable type is correct
        if (var.getType() != NC_FLOAT) {
            result.first = false;
            result.second += "   Test 3: Expected variable type NC_FLOAT, but got a different type\n";
        }
    } else {
        // Test 4: Ensure the set is not empty
        result.first = false;
        result.second += "   Test 4: varSet is empty\n";
    }

    // Test 5: Verify that getVars returns an empty set for a non-existent variable name
    set<ScaledNcVar> emptySet = group.getVars("nonExistentVar");
    if (!emptySet.empty()) {
        result.first = false;
        result.second += "   Test 5: getVars with non-existent name should return empty set\n";
    }

    if (remove(fileName.c_str()) != 0) {
        cerr << "-E- Couldn't delete " << fileName << endl;
    }

    return result;
}

TestsResult testScaledNcGroupGetVar() {
    TestsResult result = {true, ""};

    string fileName = "scaledNcTest.nc";
    ScaledNcFile file(fileName, NcFile::replace);
    ScaledNcGroup group = file.addGroup("groupTest");
    group.addVar("firstVar", NC_FLOAT);

    ScaledNcVar firstVar = group.getVar("firstVar");

    // Test 1: Ensure getVar returns a non-null
    if (firstVar.isNull()) {
        result.first = false;
        result.second += "   Test 1: getVar returned a null ScaledNcVar\n";
    }

    // Test 2: Type should stay consistent
    if (firstVar.getType().getTypeClass() != NC_FLOAT) {
        result.first = false;
        result.second += "   Test 2: Variable type is not NC_FLOAT as expected\n";
    }

    // Test 3: Getting a non-existent variable should throw an exception
    try {
        ScaledNcVar nonExistentVar = group.getVar("nonExistentVar");
        if(!nonExistentVar.isNull()) {
            result.first = false;
            result.second += "   Test 3: getVar for a non-existent variable should be NULL Var\n";
        }
    } catch (const netCDF::exceptions::NcException &e) {
        // Exception caught
        result.first = false;
        result.second += "   Test 3: getVar should not throw an NcException for a non-existent variable\n";
    } catch (const std::invalid_argument &e) {
        // Exception caught
        result.first = false;
        result.second += "   Test 3: getVar should not throw an exception for a non-existent variable\n";
    }

    // Test 4: Getting a variable with an empty name should throw an exception
    try {
        ScaledNcVar emptyNameVar = group.getVar("");
        if(!emptyNameVar.isNull()) {
            result.first = false;
            result.second += "   Test 4: getVar should be a NULL Var for an empty variable name\n";
        }
    } catch (const netCDF::exceptions::NcException &e) {
        // Exception caught
        result.first = false;
        result.second += "   Test 4: getVar should not throw an NcException for an empty variable name\n";
    } catch (const std::invalid_argument &e) {
        // Exception caught
        result.first = false;
        result.second += "   Test 4: getVar should not throw an exception for an empty variable name\n";
    }

    // Test 5: Getting a variable from a different group should throw an exception
    ScaledNcGroup anotherGroup = file.addGroup("anotherGroup");
    anotherGroup.addVar("secondVar", NC_DOUBLE);
    try {
        ScaledNcVar secondVar = group.getVar("secondVar");
        if(!secondVar.isNull()) {
            result.first = false;
            result.second +=
                "   Test 5: getVar should return NULL Var when accessing a variable from a different group\n";
        }
    } catch (const netCDF::exceptions::NcException &e) {
        // Exception caught as expected
        result.first = false;
        result.second +=
            "   Test 5: getVar should not throw an NcException when accessing a variable from a different group\n";
    } catch (const std::invalid_argument &e) {
        // Exception caught as expected
        result.first = false;
        result.second +=
            "   Test 5: getVar should not throw an exception when accessing a variable from a different group\n";
    }

    if (remove(fileName.c_str()) != 0) {
        cerr << "-E- Couldn't delete " << fileName << endl;
    }

    return result;
}

TestsResult testScaledNcGroupAddVar() {
    TestsResult result = {true, ""};

    string fileName = "scaledNcTest.nc";
    ScaledNcFile file(fileName, NcFile::replace);
    ScaledNcGroup group = file.addGroup("groupTest");

    // Test 1: Add a variable with name and type
    ScaledNcVar var1 = group.addVar("var1", NC_FLOAT);
    if (var1.isNull() || var1.getName() != "var1" || var1.getType() != NC_FLOAT) {
        result.first = false;
        result.second += "   Test 1: Failed to add variable with name and type\n";
    }

    // Test 2: Add a variable with name, type name, and dimension name
    NcDim dim1 = file.addDim("dim1", 10);
    ScaledNcVar var2 = group.addVar("var2", "double", "dim1");
    if (var2.isNull() || var2.getName() != "var2" || var2.getType() != NC_DOUBLE || var2.getDimCount() != 1) {
        result.first = false;
        result.second += "   Test 2: Failed to add variable with name, type name, and dimension name\n";
    }

    // Test 3: Add a variable with name, type, and dimension
    ScaledNcVar var3 = group.addVar("var3", NC_INT, dim1);
    if (var3.isNull() || var3.getName() != "var3" || var3.getType() != NC_INT || var3.getDimCount() != 1) {
        result.first = false;
        result.second += "   Test 3: Failed to add variable with name, type, and dimension\n";
    }

    // Test 4: Add a variable with name, type name, and vector of dimension names
    NcDim dim2 = file.addDim("dim2", 5);
    vector<string> dimNames = {"dim1", "dim2"};
    ScaledNcVar var4 = group.addVar("var4", "float", dimNames);
    if (var4.isNull() || var4.getName() != "var4" || var4.getType() != NC_FLOAT || var4.getDimCount() != 2) {
        result.first = false;
        result.second +=
            "   Test 4: Failed to add variable with name, type name, and vector of dimension names\n";
    }

    // Test 5: Add a variable with name, type, and vector of dimensions
    vector<NcDim> dims = {dim1, dim2};
    ScaledNcVar var5 = group.addVar("var5", NC_DOUBLE, dims);
    if (var5.isNull() || var5.getName() != "var5" || var5.getType() != NC_DOUBLE || var5.getDimCount() != 2) {
        result.first = false;
        result.second += "   Test 5: Failed to add variable with name, type, and vector of dimensions\n";
    }

    // Test 6: Add a variable with an existing name (should fail)
    try {
        group.addVar("var1", NC_INT);
        result.first = false;
        result.second += "   Test 6: Adding variable with existing name should throw an exception\n";
    } catch (const netCDF::exceptions::NcException &e) {
        // Exception caught as expected
    }

    // Test 7: Add a variable with an empty name (should fail)
    try {
        group.addVar("", NC_CHAR);
        result.first = false;
        result.second += "   Test 7: Adding variable with empty name should throw an exception\n";
    } catch (const netCDF::exceptions::NcException &e) {
        // Exception caught as expected
    }

    // Test 8: Add a variable with invalid type (should fail)
    try {
        group.addVar("invalidVar", static_cast<NcType::ncType>(999));
        result.first = false;
        result.second += "   Test 8: Adding variable with invalid type should throw an exception\n";
    } catch (const netCDF::exceptions::NcException &e) {
        // Exception caught as expected
    }

    if (remove(fileName.c_str()) != 0) {
        cerr << "-E- Couldn't delete " << fileName << endl;
    }

    return result;
}

TestsResult testScaledNcGroupGetGroups() {
    TestsResult result = {true, ""};

    string fileName = "scaledNcGroupGetGroupsTest.nc";
    ScaledNcFile file(fileName, NcFile::replace);

    // Test 1: Add and retrieve groups
    ScaledNcGroup group1 = file.addGroup("group1");
    ScaledNcGroup group2 = file.addGroup("group2");
    ScaledNcGroup nestedGroup = group1.addGroup("nestedGroup");

    multimap<string, ScaledNcGroup> groups = file.getGroups();
    if (groups.size() != 2) {
        result.first = false;
        result.second += "   Test 1: Expected 2 groups, but got " + to_string(groups.size()) + "\n";
    }

    // Test 2: Verify group names
    set<string> expectedNames = {"group1", "group2"};
    set<string> actualNames;
    for (const auto &pair : groups) {
        actualNames.insert(pair.first);
    }
    if (expectedNames != actualNames) {
        result.first = false;
        result.second += "   Test 2: Group names do not match expected names\n";
    }

    // Test 3: Get a specific group
    set<ScaledNcGroup> group1Set = file.getGroups("group1");
    if (group1Set.size() != 1) {
        result.first = false;
        result.second += "   Test 3: Expected 1 group, but got " + to_string(group1Set.size()) + "\n";
    }

    // Test 4: Try to get a non-existent group
    set<ScaledNcGroup> nonExistentSet = file.getGroups("nonexistent");
    if (!nonExistentSet.empty()) {
        result.first = false;
        result.second += "   Test 4: Expected empty set for non-existent group\n";
    }

    // Test 5: Get nested group
    multimap<string, ScaledNcGroup> nestedGroups = group1.getGroups();
    if (nestedGroups.size() != 1 || nestedGroups.begin()->first != "nestedGroup") {
        result.first = false;
        result.second += "   Test 5: Failed to retrieve nested group\n";
    }

    // Test 6: Get groups from an empty group
    ScaledNcGroup emptyGroup = file.addGroup("emptyGroup");
    multimap<string, ScaledNcGroup> emptyGroupChildren = emptyGroup.getGroups();
    if (!emptyGroupChildren.empty()) {
        result.first = false;
        result.second += "   Test 6: Expected empty group to have no children\n";
    }

    // Test 7: Test getParentGroup() for root group and nested groups
    ScaledNcGroup parentOfRoot = file.getParentGroup();
    if (!parentOfRoot.isNull()) {
        result.first = false;
        result.second += "   Test 7: Root group should have null parent\n";
    }

    ScaledNcGroup nestedGroup1 = file.addGroup("nestedGroup1");
    ScaledNcGroup nestedGroup2 = nestedGroup1.addGroup("nestedGroup2");

    ScaledNcGroup parentOfNested1 = nestedGroup1.getParentGroup();
    if (parentOfNested1.isNull() || parentOfNested1.getName() != file.getName()) {
        result.first = false;
        result.second += "   Test 7: Incorrect parent for nestedGroup1\n";
    }

    ScaledNcGroup parentOfNested2 = nestedGroup2.getParentGroup();
    if (parentOfNested2.isNull() || parentOfNested2.getName() != nestedGroup1.getName()) {
        result.first = false;
        result.second += "   Test 7: Incorrect parent for nestedGroup2\n";
    }

    if (remove(fileName.c_str()) != 0) {
        cerr << "-E- Couldn't delete " << fileName << endl;
    }

    return result;
}

TestsResult testScaledNcGroupGetGroup() {
    TestsResult result = {true, ""};

    string fileName = "scaledNcGroupGetGroupTest.nc";
    ScaledNcFile file(fileName, NcFile::replace);

    // Test 1: Add and retrieve a group
    ScaledNcGroup group1 = file.addGroup("group1");
    ScaledNcGroup retrievedGroup = file.getGroup("group1");
    if (retrievedGroup.isNull() || retrievedGroup.getName() != "group1") {
        result.first = false;
        result.second += "   Test 1: Failed to retrieve added group\n";
    }

    // Test 2: Try to get a non-existent group (should return a null ScaledNcGroup)
    ScaledNcGroup nonExistentGroup = file.getGroup("nonexistent");
    if (!nonExistentGroup.isNull()) {
        result.first = false;
        result.second += "   Test 2: Getting non-existent group should return a null ScaledNcGroup\n";
    }

    // Test 3: Add a nested group and retrieve it
    ScaledNcGroup nestedGroup = group1.addGroup("nestedGroup");
    ScaledNcGroup retrievedNestedGroup = group1.getGroup("nestedGroup");
    if (retrievedNestedGroup.isNull() || retrievedNestedGroup.getName() != "nestedGroup") {
        result.first = false;
        result.second += "   Test 3: Failed to retrieve nested group\n";
    }

    // Test 4: Try to get a group with an empty name (should return a null ScaledNcGroup)
    ScaledNcGroup emptyNameGroup = file.getGroup("");
    if (!emptyNameGroup.isNull()) {
        result.first = false;
        result.second += "   Test 4: Getting empty-name group should return a null ScaledNcGroup\n";
    }

    // Test 5: Verify that getGroup returns the same group object for multiple calls
    ScaledNcGroup group1FirstCall = file.getGroup("group1");
    ScaledNcGroup group1SecondCall = file.getGroup("group1");
    if (group1FirstCall.getName() != group1SecondCall.getName() ||
        group1FirstCall.getId() != group1SecondCall.getId()) {
        result.first = false;
        result.second += "   Test 5: getGroup should return the same group object for multiple calls\n";
    }

    if (remove(fileName.c_str()) != 0) {
        cerr << "-E- Couldn't delete " << fileName << endl;
    }

    return result;
}

TestsResult testScaledNcFileCreate() {
    TestsResult result = {true, ""};

    // Test 1: Create a new file
    string fileName = "scaledNcFileTest.nc";
    try {
        ScaledNcFile file(fileName, NcFile::replace);
        if (file.isNull()) {
            result.first = false;
            result.second += "   Test 1: Failed to create a new file\n";
        }
    } catch (const netCDF::exceptions::NcException &e) {
        result.first = false;
        result.second += "   Test 1: Exception while creating a new file: " + string(e.what()) + "\n";
    }

    // Test 2: Open an existing file
    try {
        ScaledNcFile file(fileName, NcFile::read);
        if (file.isNull()) {
            result.first = false;
            result.second += "   Test 2: Failed to open an existing file\n";
        }
    } catch (const netCDF::exceptions::NcException &e) {
        result.first = false;
        result.second += "   Test 2: Exception while opening an existing file: " + string(e.what()) + "\n";
    }

    // Test 3: Try to open a non-existent file (should fail)
    try {
        ScaledNcFile file("nonexistent.nc", NcFile::read);
        result.first = false;
        result.second += "   Test 3: Opening a non-existent file should throw an exception\n";
    } catch (const netCDF::exceptions::NcException &e) {
        // Exception caught as expected
    }

    if (remove(fileName.c_str()) != 0) {
        cerr << "-E- Couldn't delete " << fileName << endl;
    }

    return result;
}

TestsResult testScaledNcFileAddGroup() {
    TestsResult result = {true, ""};

    string fileName = "scaledNcFileTest.nc";
    ScaledNcFile file(fileName, NcFile::replace);

    // Test 1: Add a new group
    ScaledNcGroup group1 = file.addGroup("group1");
    if (group1.isNull() || group1.getName() != "group1") {
        result.first = false;
        result.second += "   Test 1: Failed to add a new group\n";
    }

    // Test 2: Add a nested group
    ScaledNcGroup nestedGroup = group1.addGroup("nestedGroup");
    if (nestedGroup.isNull() || nestedGroup.getName() != "nestedGroup") {
        result.first = false;
        result.second += "   Test 2: Failed to add a nested group\n";
    }

    // Test 3: Try to add a group with an existing name (should fail)
    try {
        file.addGroup("group1");
        result.first = false;
        result.second += "   Test 3: Adding a group with an existing name should throw an exception\n";
    } catch (const netCDF::exceptions::NcException &e) {
        // Exception caught as expected
    }

    // Test 4: Try to add a group with an empty name (should fail)
    try {
        file.addGroup("");
        result.first = false;
        result.second += "   Test 4: Adding a group with an empty name should throw an exception\n";
    } catch (const netCDF::exceptions::NcException &e) {
        // Exception caught as expected
    }

    if (remove(fileName.c_str()) != 0) {
        cerr << "-E- Couldn't delete " << fileName << endl;
    }

    return result;
}

TestsResult testScaledNcFileGetGroups() {
    TestsResult result = {true, ""};

    string fileName = "scaledNcFileTest.nc";
    ScaledNcFile file(fileName, NcFile::replace);

    file.addGroup("group1");
    file.addGroup("group2");

    // Test 1: Get all groups
    multimap<string, ScaledNcGroup> groups = file.getGroups();
    if (groups.size() != 2) {
        result.first = false;
        result.second += "   Test 1: Expected 2 groups, but got " + to_string(groups.size()) + "\n";
    }

    // Test 2: Verify group names
    set<string> expectedNames = {"group1", "group2"};
    set<string> actualNames;
    for (const auto &pair : groups) {
        actualNames.insert(pair.first);
    }
    if (expectedNames != actualNames) {
        result.first = false;
        result.second += "   Test 2: Group names do not match expected names\n";
    }

    // Test 3: Get a specific group
    set<ScaledNcGroup> group1Set = file.getGroups("group1");
    if (group1Set.size() != 1) {
        result.first = false;
        result.second += "   Test 3: Expected 1 group, but got " + to_string(group1Set.size()) + "\n";
    }

    // Test 4: Try to get a non-existent group
    set<ScaledNcGroup> nonExistentSet = file.getGroups("nonexistent");
    if (!nonExistentSet.empty()) {
        result.first = false;
        result.second += "   Test 4: Expected empty set for non-existent group\n";
    }

    // Test 5: Get parent group
    ScaledNcGroup parentGroup = file.getParentGroup();
    if (!parentGroup.isNull()) {
        result.first = false;
        result.second += "   Test 5: Expected null parent group for root group\n";
    }

    // Test 6: Get parent group of a non-root group
    ScaledNcGroup childGroup = file.addGroup("childGroup");
    ScaledNcGroup parentOfChild = childGroup.getParentGroup();
    if (parentOfChild.isNull() || parentOfChild.getName() != file.getName()) {
        result.first = false;
        result.second += "   Test 6: Incorrect parent group for child group\n";
    }

    if (remove(fileName.c_str()) != 0) {
        cerr << "-E- Couldn't delete " << fileName << endl;
    }

    return result;
}

TestsResult testScaledNcFileAddVar() {
    TestsResult result = {true, ""};

    string fileName = "scaledNcFileTest.nc";
    ScaledNcFile file(fileName, NcFile::replace);

    // Test 1: Add a variable with name and type
    ScaledNcVar var1 = file.addVar("var1", NC_FLOAT);
    if (var1.isNull() || var1.getName() != "var1" || var1.getType() != NC_FLOAT) {
        result.first = false;
        result.second += "   Test 1: Failed to add variable with name and type\n";
    }

    // Test 2: Add a variable with name, type name, and dimension name
    NcDim dim1 = file.addDim("dim1", 10);
    ScaledNcVar var2 = file.addVar("var2", "double", "dim1");
    if (var2.isNull() || var2.getName() != "var2" || var2.getType() != NC_DOUBLE || var2.getDimCount() != 1) {
        result.first = false;
        result.second += "   Test 2: Failed to add variable with name, type name, and dimension name\n";
    }

    // Test 3: Add a variable with name, type, and dimension
    ScaledNcVar var3 = file.addVar("var3", NC_INT, dim1);
    if (var3.isNull() || var3.getName() != "var3" || var3.getType() != NC_INT || var3.getDimCount() != 1) {
        result.first = false;
        result.second += "   Test 3: Failed to add variable with name, type, and dimension\n";
    }

    // Test 4: Add a variable with name, type name, and vector of dimension names
    NcDim dim2 = file.addDim("dim2", 5);
    vector<string> dimNames = {"dim1", "dim2"};
    ScaledNcVar var4 = file.addVar("var4", "float", dimNames);
    if (var4.isNull() || var4.getName() != "var4" || var4.getType() != NC_FLOAT || var4.getDimCount() != 2) {
        result.first = false;
        result.second +=
            "   Test 4: Failed to add variable with name, type name, and vector of dimension names\n";
    }

    // Test 5: Add a variable with name, type, and vector of dimensions
    vector<NcDim> dims = {dim1, dim2};
    ScaledNcVar var5 = file.addVar("var5", NC_DOUBLE, dims);
    if (var5.isNull() || var5.getName() != "var5" || var5.getType() != NC_DOUBLE || var5.getDimCount() != 2) {
        result.first = false;
        result.second += "   Test 5: Failed to add variable with name, type, and vector of dimensions\n";
    }

    // Test 6: Add a variable with an existing name (should fail)
    try {
        file.addVar("var1", NC_INT);
        result.first = false;
        result.second += "   Test 6: Adding variable with existing name should throw an exception\n";
    } catch (const netCDF::exceptions::NcException &e) {
        // Exception caught as expected
    }

    // Test 7: Add a variable with an empty name (should fail)
    try {
        file.addVar("", NC_CHAR);
        result.first = false;
        result.second += "   Test 7: Adding variable with empty name should throw an exception\n";
    } catch (const netCDF::exceptions::NcException &e) {
        // Exception caught as expected
    }

    // Test 8: Add a variable with invalid type (should fail)
    try {
        file.addVar("invalidVar", static_cast<NcType::ncType>(999));
        result.first = false;
        result.second += "   Test 8: Adding variable with invalid type should throw an exception\n";
    } catch (const netCDF::exceptions::NcException &e) {
        // Exception caught as expected
    }

    if (remove(fileName.c_str()) != 0) {
        cerr << "-E- Couldn't delete " << fileName << endl;
    }

    return result;
}

TestsResult testScaledNcFileGetVars() {
    TestsResult result = {true, ""};

    string fileName = "scaledNcFileTest.nc";
    ScaledNcFile file(fileName, NcFile::replace);

    // Add dimensions and variables
    NcDim dim1 = file.addDim("dim1", 10);
    NcDim dim2 = file.addDim("dim2", 5);
    file.addVar("var1", NC_FLOAT);
    file.addVar("var2", NC_DOUBLE, dim1);
    file.addVar("var3", NC_INT, {dim1, dim2});

    // Test 1: Get all variables
    multimap<string, ScaledNcVar> allVars = file.getVars();
    if (allVars.size() != 3) {
        result.first = false;
        result.second += "   Test 1: Expected 3 variables, but got " + to_string(allVars.size()) + "\n";
    }

    // Test 2: Get variables by name
    set<ScaledNcVar> var1Set = file.getVars("var1");
    if (var1Set.size() != 1) {
        result.first = false;
        result.second += "   Test 2: Expected 1 variable, but got " + to_string(var1Set.size()) + "\n";
    }

    // Test 3: Get non-existent variable
    set<ScaledNcVar> nonExistentSet = file.getVars("nonexistent");
    if (!nonExistentSet.empty()) {
        result.first = false;
        result.second += "   Test 3: Expected empty set for non-existent variable\n";
    }

    if (remove(fileName.c_str()) != 0) {
        cerr << "-E- Couldn't delete " << fileName << endl;
    }

    return result;
}

TestsResult testScaledNcFileGetVar() {
    TestsResult result = {true, ""};

    string fileName = "scaledNcFileGetVarTest.nc";
    ScaledNcFile file(fileName, NcFile::replace);

    // Add a variable to the file
    NcDim dim = file.addDim("dim", 10);
    ScaledNcVar originalVar = file.addVar("testVar", NC_FLOAT, dim);

    // Test 1: Get an existing variable
    ScaledNcVar retrievedVar = file.getVar("testVar");
    if (retrievedVar.isNull() || retrievedVar.getName() != "testVar") {
        result.first = false;
        result.second += "   Test 1: Failed to retrieve existing variable\n";
    }

    // Test 2: Try to get a non-existent variable (should throw an exception)
    try {
        ScaledNcVar nonExistentVar = file.getVar("nonExistentVar");
        if(!nonExistentVar.isNull()) {
            result.first = false;
            result.second += "   Test 2: Getting non-existent variable should return a NULL Var\n";
        }
    } catch (const netCDF::exceptions::NcException &e) {
        // Exception caught
        result.first = false;
        result.second += "   Test 2: Getting variable with empty name should not throw NcException\n";
    } catch (const std::invalid_argument &e) {
        // Exception caught
        result.first = false;
        result.second += "   Test 2: Getting variable with empty name should not throw exception\n";
    }

    // Test 3: Try to get a variable with an empty name (should return NULL Var)
    try {
        ScaledNcVar emptyNameVar = file.getVar("");
        if(!emptyNameVar.isNull()) {
            // Exception caught as expected
            result.first = false;
            result.second += "   Test 3: Getting variable with empty name should return a NULL Var\n";
        }
    } catch (const netCDF::exceptions::NcException &e) {
        // Exception caught
        result.first = false;
        result.second += "   Test 3: Getting variable with empty name should not throw NcException\n";
    } catch (const std::invalid_argument &e) {
        // Exception caught
        result.first = false;
        result.second += "   Test 3: Getting variable with empty name should not throw exception\n";
    }

    // Test 4: Check if the retrieved variable has the correct type
    if (retrievedVar.getType() != NC_FLOAT) {
        result.first = false;
        result.second += "   Test 4: Retrieved variable has incorrect type\n";
    }

    // Test 5: Check if the retrieved variable has the correct dimension
    if (retrievedVar.getDimCount() != 1 || retrievedVar.getDim(0).getName() != "dim") {
        result.first = false;
        result.second += "   Test 5: Retrieved variable has incorrect dimension\n";
    }

    if (remove(fileName.c_str()) != 0) {
        cerr << "-E- Couldn't delete " << fileName << endl;
    }

    return result;
}

TestsResult testScaledNcFileDestructor() {
    TestsResult result = {true, ""};

    string fileName = "scaledNcFileDestructorTest.nc";
    {
        ScaledNcFile file(fileName, NcFile::replace);
        file.addGroup("testGroup");
        file.addVar("testVar", NC_FLOAT);
    }  // file goes out of scope here, destructor should be called

    // Try to open the file again to ensure it was properly closed
    try {
        ScaledNcFile checkFile(fileName, NcFile::read);
        // If we can open the file, it was properly closed
    } catch (const netCDF::exceptions::NcException &e) {
        result.first = false;
        result.second += "   Test: File was not properly closed by destructor\n";
    }

    if (remove(fileName.c_str()) != 0) {
        cerr << "-E- Couldn't delete " << fileName << endl;
    }

    return result;
}

int main() {
    bool allTestsPassed = true;
    TestsResult testsResult = {true, ""};  // Result, which tests if any failed
    vector<pair<string, function<TestsResult()>>> tests = {
        {"ScaledNcGroup getVars() test", testScaledNcGroupMultimapGetVars},
        {"ScaledNcGroup getVars(string) test", testScaledNcGroupSetGetVars},
        {"ScaledNcGroup getVar(string) test", testScaledNcGroupGetVar},
        {"ScaledNcGroup addVar test", testScaledNcGroupAddVar},
        {"ScaledNcGroup getGroups test", testScaledNcGroupGetGroups},
        {"ScaledNcGroup getGroup test", testScaledNcGroupGetGroup},
        {"ScaledNcFile creation test", testScaledNcFileCreate},
        {"ScaledNcFile addGroup test", testScaledNcFileAddGroup},
        {"ScaledNcFile getGroups test", testScaledNcFileGetGroups},
        {"ScaledNcFile addVar test", testScaledNcFileAddVar},
        {"ScaledNcFile getVars test", testScaledNcFileGetVars},
        {"ScaledNcFile getVar test", testScaledNcFileGetVar},
        {"ScaledNcFile destructor test", testScaledNcFileDestructor},
        {"ScaledNcVar writing test", testWriting},
        {"ScaledNcVar reading test", testReading}};

    int totalTests = tests.size();
    int numTestsPassed = 0;

    for (const pair<string, function<TestsResult()>> &test : tests) {
        testsResult = test.second();
        string result = testsResult.first ? "PASSED" : "FAILED";
        cout << left << setw(40) << test.first << setfill('.') << setw(50) << right << result << endl;
        if (testsResult.first) {
            numTestsPassed++;
        } else {
            cout << testsResult.second << endl;
            allTestsPassed = false;
        }
    }

    double percentPassed = (static_cast<double>(numTestsPassed) / totalTests) * 100.0;
    cout << endl
         << "Percentage of tests passed: " << fixed << setprecision(2) << percentPassed << "%" << endl;

    return allTestsPassed ? 0 : 1;  // Return 0 if all tests pass, 1 otherwise
}
