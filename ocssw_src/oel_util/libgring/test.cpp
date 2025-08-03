#include "gring.h"
#include <utility>  // std::pair
#include <string>
#include <vector>
#include <functional>
#include <iostream>
#include <iomanip>

using namespace std;

typedef pair<bool, string> TestsResult;

constexpr double COMPARISON_PRECISION = 0.00001;

TestsResult testMinsAndMaxes(float expectedMin, float expectedMax, float actualMin, float actualMax) {
    float expectedVsActualMinDiff = expectedMin - actualMin;
    float expectedVsActualMaxDiff = expectedMax - actualMax;

    bool minDiffWithinBounds = expectedVsActualMinDiff <= COMPARISON_PRECISION;
    bool maxDiffWithinBounds = expectedVsActualMaxDiff <= COMPARISON_PRECISION;

    if (minDiffWithinBounds && maxDiffWithinBounds) {
        return {true, "Test passed"};
    } else {
        return {false, "    Expected min: " + to_string(expectedMin) + ", max: " + to_string(expectedMax) +
                       ". \n    Got min:       " + to_string(actualMin) + ", max: " + to_string(actualMax)};
    }
}

TestsResult testGringLonExtremesWhenBothArePositive() {
    vector<float> lats = {1.0, 2.0, 3.0};
    Gring gring;

    vector<float> lons = {1.0, 2.0, 3.0};
    gring.processScan(lats.data(), lons.data(), lats.size(), 1, 0);

    float expectedMin = 1.0;
    float expectedMax = 3.0;
    float actualMin = gring.getGeospatialLongitudeMin();
    float actualMax = gring.getGeospatialLongitudeMax();

    return testMinsAndMaxes(expectedMin, expectedMax, actualMin, actualMax);
}

TestsResult testGringLonExtremesWhenBothAreNegative() {
    vector<float> lats = {1.0, 2.0, 3.0};
    Gring gring;

    vector<float> lons = {-1.0, -2.0, -3.0};
    gring.processScan(lats.data(), lons.data(), lats.size(), 1, 0);

    float expectedMin = -3.0;
    float expectedMax = -1.0;
    float actualMin = gring.getGeospatialLongitudeMin();
    float actualMax = gring.getGeospatialLongitudeMax();

    return testMinsAndMaxes(expectedMin, expectedMax, actualMin, actualMax);
}

TestsResult testGringLonExtremesWhenDataSpansPrimeMeridian() {
    vector<float> lats = {1.0, 2.0, 3.0, 4.0};
    Gring gring;

    vector<float> lons = {-10.0, -5.0, 5.0, 10.0};
    gring.processScan(lats.data(), lons.data(), lats.size(), 1, 0);

    float expectedMin = -10.0;
    float expectedMax = 10.0;
    float actualMin = gring.getGeospatialLongitudeMin();
    float actualMax = gring.getGeospatialLongitudeMax();

    return testMinsAndMaxes(expectedMin, expectedMax, actualMin, actualMax);
}

TestsResult testGringLonExtremesWhenDataSpansAntimeridian() {
    vector<float> lats = {1.0, 2.0, 3.0, 4.0};
    Gring gring;

    vector<float> lons = {175.0, 179.0, -179.0, -175.0};
    gring.processScan(lats.data(), lons.data(), lats.size(), 1, 0);

    float expectedMin = 175.0;
    float expectedMax = -175.0;
    float actualMin = gring.getGeospatialLongitudeMin();
    float actualMax = gring.getGeospatialLongitudeMax();

    return testMinsAndMaxes(expectedMin, expectedMax, actualMin, actualMax);
}

TestsResult testGringLonExtremesFullGlobe() {
    vector<float> lats = {90, 45, 0, -45, 90};
    Gring gring;

    vector<float> lons = {180, 90, 0, -90, -180};
    gring.processScan(lats.data(), lons.data(), lats.size(), 1, 0);

    float expectedMin = 180;
    float expectedMax = -180;
    float actualMin = gring.getGeospatialLongitudeMin();
    float actualMax = gring.getGeospatialLongitudeMax();

    return testMinsAndMaxes(expectedMin, expectedMax, actualMin, actualMax);
}

TestsResult testGetSameResultPACE_OCI20240404T141256L1Bephnc() {

    vector<float> latsFirstScan = {-3.19418, 1.83129};
    vector<float> lonsFirstScan = {-29.85341, -6.20588};
    vector<float> latsLastScan = {5.63751, 10.69493};
    vector<float> lonsLastScan = {-31.78984, -7.91932};

    Gring gring;
    gring.processScan(latsFirstScan.data(), lonsFirstScan.data(), latsFirstScan.size(), 0, 1);
    gring.processScan(latsLastScan.data(), lonsLastScan.data(), latsLastScan.size(), 1, 1);

    float expectedLonMin = -31.789984;
    float expectedLonMax = -6.20588;
    float actualLonMin = gring.getGeospatialLongitudeMin();
    float actualLonMax = gring.getGeospatialLongitudeMax();

    return testMinsAndMaxes(expectedLonMin, expectedLonMax, actualLonMin, actualLonMax);

}

TestsResult testGringLatExtremesWhenBothAreNegative() {
    vector<float> lons = {1.0, 2.0, 3.0};
    vector<float> lats = {-1.0, -2.0, -3.0};

    Gring gring;

    gring.processScan(lats.data(), lons.data(), lats.size(), 0, 0);
    float expectedMin = -3.0;
    float expectedMax = -1.0;
    float actualMin = gring.getGeospatialLatitudeMin();
    float actualMax = gring.getGeospatialLatitudeMax();

    return testMinsAndMaxes(expectedMin, expectedMax, actualMin, actualMax);
}

TestsResult testGringLatExtremesWhenBothArePositive() {
    vector<float> lats = {1.0, 2.0, 3.0};
    vector<float> lons = {-1.0, -2.0, -3.0};

    Gring gring;

    gring.processScan(lats.data(), lons.data(), lats.size(), 0, 0);
    float expectedMin = 1.0;
    float expectedMax = 3.0;
    float actualMin = gring.getGeospatialLatitudeMin();
    float actualMax = gring.getGeospatialLatitudeMax();

    return testMinsAndMaxes(expectedMin, expectedMax, actualMin, actualMax);
}

TestsResult testGringLatExtremesWhenDataSpansEquator() {
    vector<float> lats = {-1.0, 0, 1.0};
    vector<float> lons = {-1.0, -2.0, -3.0};

    Gring gring;

    gring.processScan(lats.data(), lons.data(), lats.size(), 0, 0);
    float expectedMin = -1.0;
    float expectedMax = 1.0;
    float actualMin = gring.getGeospatialLatitudeMin();
    float actualMax = gring.getGeospatialLatitudeMax();

    return testMinsAndMaxes(expectedMin, expectedMax, actualMin, actualMax);
}

void printTestResult(const pair<string, function<TestsResult()>> &test, TestsResult testResult) {

    const string& testName = test.first;
    const string result = testResult.first ? "PASSED" : "FAILED";
    const int totalWidth = 65;
    int dotCount = totalWidth - testName.length() - result.length();
    if (dotCount < 0) 
        dotCount = 0;

    cout << testName
        << string(dotCount, '.')
        << result
        << endl;

    if (!testResult.first) {
        cout << testResult.second << endl;
    }
}

TestsResult testLonMaxOnFirstLine() {
    vector<vector<float>> lats = { // Not a factor, but still needed
        {0.0, 0.0},
        {0.0, 0.0}
    };

    vector<vector<float>> lons = {
        {0.0, 2.0},
        {0.0, 1.0}
    };

    Gring gring;

    for (size_t i = 0; i < lons.size(); i++) {
        gring.processScan(lats[i].data(), lons[i].data(), lats.size(), i, 1);
    }
    
    float expectedMin = lons[0][0];
    float expectedMax = lons[0][1];
    float actualMin = gring.getGeospatialLongitudeMin();
    float actualMax = gring.getGeospatialLongitudeMax();

    return testMinsAndMaxes(expectedMin, expectedMax, actualMin, actualMax);
}

TestsResult testLonMaxOnFirstLineAllNegative() {
    vector<vector<float>> lats = { // Not a factor, but still needed
        {0.0, 0.0},
        {0.0, 0.0}
    };

    vector<vector<float>> lons = {
        {-61.45958, -60.00514},
        {-63.752, -62.60792}
    };

    Gring gring;

    for (size_t i = 0; i < lons.size(); i++) {
        gring.processScan(lats[i].data(), lons[i].data(), lats.size(), i, 1);
    }
    
    float expectedMin = lons[1][0];
    float expectedMax = lons[0][1];
    float actualMin = gring.getGeospatialLongitudeMin();
    float actualMax = gring.getGeospatialLongitudeMax();

    return testMinsAndMaxes(expectedMin, expectedMax, actualMin, actualMax);
}

vector<pair<string, function<TestsResult()>>> getTests() {
    return {
        {"Test GringLonExtremes when both are positive", testGringLonExtremesWhenBothArePositive},
        {"Test GringLonExtremes when both are negative", testGringLonExtremesWhenBothAreNegative},
        {"Test GringLonExtremes when data spans prime meridian",
         testGringLonExtremesWhenDataSpansPrimeMeridian},
        {"Test GringLonExtremes when data spans antimeridian", testGringLonExtremesWhenDataSpansAntimeridian},
        {"Test GringLonExtremes full globe", testGringLonExtremesFullGlobe},
        {"Test PACE OCI definitive ephemeris ctest", testGetSameResultPACE_OCI20240404T141256L1Bephnc},
        {"Test lat extremes when both are positive", testGringLatExtremesWhenBothArePositive},
        {"Test lat extremes when both are negative", testGringLatExtremesWhenBothAreNegative},
        {"Test lat extremes when data spans the equator", testGringLatExtremesWhenDataSpansEquator},
        {"Test lon extremes when max value is on first line", testLonMaxOnFirstLine},
        {"Test lon extremes when all negative and max on first line", testLonMaxOnFirstLineAllNegative},
    };
}

int main() {
    vector<pair<string, function<TestsResult()>>> tests = getTests();

    size_t testsPassed = 0;

    for (const pair<string, function<TestsResult()>> &test : tests) {
        TestsResult testResult = test.second();

        printTestResult(test, testResult);

        if (testResult.first) {
            testsPassed++;
        }
    }

    double percentPassed = (static_cast<double>(testsPassed) / tests.size()) * 100.0;
    cout << endl
         << testsPassed << " of " << tests.size() << " tests passed (" << fixed << setprecision(2)
         << percentPassed << "\%)" << endl;

    return testsPassed != tests.size();
}
