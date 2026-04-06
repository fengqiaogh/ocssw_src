#include "gring.h"
#include <utility>  // std::pair
#include <string>
#include <vector>
#include <functional>
#include <iostream>
#include <iomanip>

#include <sstream>  // For CoutRedirector

class OutputRedirector {
   public:
    OutputRedirector(std::ostream& stream = std::cerr) : oldBuf_(stream.rdbuf()), stream_(stream) {
        stream_.rdbuf(capture_.rdbuf());
    }

    ~OutputRedirector() {
        stream_.rdbuf(oldBuf_);  // restore original buffer
    }

    std::string str() const {
        return capture_.str();
    }

   private:
    std::ostringstream capture_;
    std::streambuf* oldBuf_;
    std::ostream& stream_;
};

using namespace std;

typedef pair<bool, string> TestsResult;

OutputRedirector _; // Don't want cerr output while testing

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
                           ". \n    Got min:       " + to_string(actualMin) +
                           ", max: " + to_string(actualMax)};
    }
}

TestsResult testGringLonExtremesWhenBothArePositive() {
    vector<float> lats = {1.0, 2.0, 3.0};
    Gring gring;

    vector<float> lons = {1.0, 2.0, 3.0};
    gring.tryIncludeScan(lats.data(), lons.data(), lats.size(), 0);

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
    gring.tryIncludeScan(lats.data(), lons.data(), lats.size(), 1, 0);

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
    gring.tryIncludeScan(lats.data(), lons.data(), lats.size(), 1, 0);

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
    gring.tryIncludeScan(lats.data(), lons.data(), lats.size(), 1, 0);

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
    gring.tryIncludeScan(lats.data(), lons.data(), lats.size(), 1, 0);

    float expectedMin = -180;
    float expectedMax = 180;
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
    gring.tryIncludeScan(latsFirstScan.data(), lonsFirstScan.data(), latsFirstScan.size(), 1);
    gring.tryIncludeScan(latsLastScan.data(), lonsLastScan.data(), latsLastScan.size(), 1);

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

    gring.tryIncludeScan(lats.data(), lons.data(), lats.size(), 0, 0);
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

    gring.tryIncludeScan(lats.data(), lons.data(), lats.size(), 0, 0);
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

    gring.tryIncludeScan(lats.data(), lons.data(), lats.size(), 0, 0);
    float expectedMin = -1.0;
    float expectedMax = 1.0;
    float actualMin = gring.getGeospatialLatitudeMin();
    float actualMax = gring.getGeospatialLatitudeMax();

    return testMinsAndMaxes(expectedMin, expectedMax, actualMin, actualMax);
}

void printTestResult(const pair<string, function<TestsResult()>>& test, TestsResult testResult) {
    const string& testName = test.first;
    const string result = testResult.first ? "PASSED" : "FAILED";
    const int totalWidth = 65;
    int dotCount = totalWidth - testName.length() - result.length();
    if (dotCount < 0)
        dotCount = 0;

    cout << testName << string(dotCount, '.') << result << endl;

    if (!testResult.first) {
        cout << testResult.second << endl;
    }
}

TestsResult testLonMaxOnFirstLine() {
    vector<vector<float>> lats = {// Not a factor, but still needed
                                  {0.0, 0.0},
                                  {0.0, 0.0}};

    vector<vector<float>> lons = {{0.0, 2.0}, {0.0, 1.0}};

    Gring gring;

    for (size_t i = 0; i < lons.size(); i++) {
        gring.tryIncludeScan(lats[i].data(), lons[i].data(), lats.size(), 1);
    }

    float expectedMin = lons[0][0];
    float expectedMax = lons[0][1];
    float actualMin = gring.getGeospatialLongitudeMin();
    float actualMax = gring.getGeospatialLongitudeMax();

    return testMinsAndMaxes(expectedMin, expectedMax, actualMin, actualMax);
}

TestsResult testLonMaxOnFirstLineAllNegative() {
    vector<vector<float>> lats = {// Not a factor, but still needed
                                  {0.0, 0.0},
                                  {0.0, 0.0}};

    vector<vector<float>> lons = {{-61.45958, -60.00514}, {-63.752, -62.60792}};

    Gring gring;

    for (size_t i = 0; i < lons.size(); i++) {
        gring.tryIncludeScan(lats[i].data(), lons[i].data(), lats.size(), 1);
    }

    float expectedMin = lons[1][0];
    float expectedMax = lons[0][1];
    float actualMin = gring.getGeospatialLongitudeMin();
    float actualMax = gring.getGeospatialLongitudeMax();

    return testMinsAndMaxes(expectedMin, expectedMax, actualMin, actualMax);
}

TestsResult testWktLonsPositive() {
    vector<double> lats = {1.0, 1.0};
    vector<double> lons = {1.0, 2.0};

    Gring gring;
    gring.tryIncludeScan(lats.data(), lons.data(), lons.size(), 0);
    lats = {2.0, 2.0};
    lons = {2.0, 3.0};
    gring.tryIncludeScan(lats.data(), lons.data(), lons.size(), 0);

    string expected =
        "POLYGON((1.00000 1.00000,2.00000 1.00000,3.00000 2.00000,2.00000 2.00000,1.00000 1.00000))";
    string returned = gring.getGeospatialBounds();
    string expectedVsReturned = "     Expected: " + expected + "\n     Got:      " + returned;

    return {expected == returned, expectedVsReturned};
}

TestsResult testWktLonsNegative() {
    vector<double> lats = {1.0, 1.0};
    vector<double> lons = {-1.0, -2.0};

    Gring gring;
    gring.tryIncludeScan(lats.data(), lons.data(), lons.size(), 0);
    lats = {2.0, 2.0};
    lons = {-2.0, -3.0};
    gring.tryIncludeScan(lats.data(), lons.data(), lons.size(), 0);

    string expected =
        "POLYGON((-1.00000 1.00000,-2.00000 2.00000,-3.00000 2.00000,-2.00000 1.00000,-1.00000 1.00000))";
    string returned = gring.getGeospatialBounds();
    string expectedVsReturned = "     Expected: " + expected + "\n     Got:      " + returned;

    return {expected == returned, expectedVsReturned};
}

TestsResult testWktLonsMixed() {
    vector<double> lats = {0.0, 0.0};
    vector<double> lons = {-2.0, -1.0};

    Gring gring;
    gring.tryIncludeScan(lats.data(), lons.data(), lons.size(), 0);
    lats = {2.0, 2.0};
    lons = {0.0, 1.0};
    gring.tryIncludeScan(lats.data(), lons.data(), lons.size(), 0);

    string expected =
        "POLYGON((-2.00000 0.00000,-1.00000 0.00000,1.00000 2.00000,0.00000 2.00000,-2.00000 0.00000))";
    string returned = gring.getGeospatialBounds();
    string expectedVsReturned = "     Expected: " + expected + "\n     Got:      " + returned;

    return {expected == returned, expectedVsReturned};
}

TestsResult testWktLatsPositive() {
    vector<double> lats = {1.0, 2.0};
    vector<double> lons = {1.0, 1.0};

    Gring gring;
    gring.tryIncludeScan(lats.data(), lons.data(), lons.size(), 0);
    lats = {2.0, 3.0};
    lons = {2.0, 2.0};
    gring.tryIncludeScan(lats.data(), lons.data(), lons.size(), 0);

    string expected =
        "POLYGON((1.00000 1.00000,2.00000 2.00000,2.00000 3.00000,1.00000 2.00000,1.00000 1.00000))";
    string returned = gring.getGeospatialBounds();
    string expectedVsReturned = "     Expected: " + expected + "\n     Got:      " + returned;

    return {expected == returned, expectedVsReturned};
}

TestsResult testWktLatsNegative() {
    vector<double> lats = {-1.0, -2.0};
    vector<double> lons = {1.0, 1.0};

    Gring gring;
    gring.tryIncludeScan(lats.data(), lons.data(), lons.size(), 0);
    lats = {-2.0, -3.0};
    lons = {2.0, 2.0};
    gring.tryIncludeScan(lats.data(), lons.data(), lons.size(), 0);

    string expected =
        "POLYGON((1.00000 -1.00000,1.00000 -2.00000,2.00000 -3.00000,2.00000 -2.00000,1.00000 -1.00000))";
    string returned = gring.getGeospatialBounds();
    string expectedVsReturned = "     Expected: " + expected + "\n     Got:      " + returned;

    return {expected == returned, expectedVsReturned};
}

TestsResult testWktLatsMixed() {
    vector<double> lats = {-2.0, -1.0};
    vector<double> lons = {0.0, 0.0};

    Gring gring;
    gring.tryIncludeScan(lats.data(), lons.data(), lons.size(), 0);
    lats = {0.0, 1.0};
    lons = {2.0, 2.0};
    gring.tryIncludeScan(lats.data(), lons.data(), lons.size(), 0);

    string expected =
        "POLYGON((0.00000 -2.00000,2.00000 0.00000,2.00000 1.00000,0.00000 -1.00000,0.00000 -2.00000))";
    string returned = gring.getGeospatialBounds();
    string expectedVsReturned = "     Expected: " + expected + "\n     Got:      " + returned;

    return {expected == returned, expectedVsReturned};
}

TestsResult testGringLatsMixed() {
    vector<double> lats = {-2.0, -1.0};
    vector<double> lons = {0.0, 0.0};

    Gring gring;
    gring.tryIncludeScan(lats.data(), lons.data(), lons.size(), 0);
    lats = {0.0, 1.0};
    lons = {2.0, 2.0};
    gring.tryIncludeScan(lats.data(), lons.data(), lons.size(), 0);

    string expectedLons = "0.00000,2.00000,2.00000,0.00000";
    string expectedLats = "-1.00000,1.00000,0.00000,-2.00000";
    string expectedSequence = "1,2,3,4";
    string gotLons, gotLats, gotSequence;
    gring.getGringStrings(gotLons, gotLats, gotSequence);

    string expected = expectedLons + "\n        " + expectedLats + "\n        " + expectedSequence;
    string returned = gotLons + "\n        " + gotLats + "\n        " + gotSequence;
    string expectedVsReturned = "    Expected:\n        " + expected + "\n    Got:\n        " + returned;

    return {expected == returned, expectedVsReturned};
}

TestsResult testGringLatsNegative() {
    vector<double> lats = {-1.0, -2.0};
    vector<double> lons = {1.0, 1.0};

    Gring gring;
    gring.tryIncludeScan(lats.data(), lons.data(), lons.size(), 0);
    lats = {-2.0, -3.0};
    lons = {2.0, 2.0};
    gring.tryIncludeScan(lats.data(), lons.data(), lons.size(), 0);

    string expectedLons = "1.00000,2.00000,2.00000,1.00000";
    string expectedLats = "-1.00000,-2.00000,-3.00000,-2.00000";
    string expectedSequence = "1,2,3,4";
    string gotLons, gotLats, gotSequence;
    gring.getGringStrings(gotLons, gotLats, gotSequence);

    string expected = expectedLons + "\n        " + expectedLats + "\n        " + expectedSequence;
    string returned = gotLons + "\n        " + gotLats + "\n        " + gotSequence;
    string expectedVsReturned = "    Expected:\n        " + expected + "\n    Got:\n        " + returned;

    return {expected == returned, expectedVsReturned};
}

TestsResult testGringLatsPositive() {
    vector<double> lats = {1.0, 2.0};
    vector<double> lons = {1.0, 1.0};

    Gring gring;
    gring.tryIncludeScan(lats.data(), lons.data(), lons.size(), 0);
    lats = {2.0, 3.0};
    lons = {2.0, 2.0};
    gring.tryIncludeScan(lats.data(), lons.data(), lons.size(), 0);

    string expectedLons = "1.00000,2.00000,2.00000,1.00000";
    string expectedLats = "2.00000,3.00000,2.00000,1.00000";
    string expectedSequence = "1,2,3,4";
    string gotLons, gotLats, gotSequence;
    gring.getGringStrings(gotLons, gotLats, gotSequence);

    string expected = expectedLons + "\n        " + expectedLats + "\n        " + expectedSequence;
    string returned = gotLons + "\n        " + gotLats + "\n        " + gotSequence;
    string expectedVsReturned = "    Expected:\n        " + expected + "\n    Got:\n        " + returned;

    return {expected == returned, expectedVsReturned};
}

TestsResult testGringLonsPositive() {
    vector<double> lats = {1.0, 1.0};
    vector<double> lons = {1.0, 2.0};

    Gring gring;
    gring.tryIncludeScan(lats.data(), lons.data(), lons.size(), 0);
    lats = {2.0, 2.0};
    lons = {2.0, 3.0};
    gring.tryIncludeScan(lats.data(), lons.data(), lons.size(), 0);

    string expectedLons = "1.00000,2.00000,3.00000,2.00000";
    string expectedLats = "1.00000,2.00000,2.00000,1.00000";
    string expectedSequence = "1,2,3,4";
    string gotLons, gotLats, gotSequence;
    gring.getGringStrings(gotLons, gotLats, gotSequence);

    string expected = expectedLons + "\n        " + expectedLats + "\n        " + expectedSequence;
    string returned = gotLons + "\n        " + gotLats + "\n        " + gotSequence;
    string expectedVsReturned = "    Expected:\n        " + expected + "\n    Got:\n        " + returned;

    return {expected == returned, expectedVsReturned};
}

TestsResult testGringLonsNegative() {
    vector<double> lats = {1.0, 1.0};
    vector<double> lons = {-1.0, -2.0};

    Gring gring;
    gring.tryIncludeScan(lats.data(), lons.data(), lons.size(), 0);
    lats = {2.0, 2.0};
    lons = {-2.0, -3.0};
    gring.tryIncludeScan(lats.data(), lons.data(), lons.size(), 0);

    string expectedLons = "-2.00000,-3.00000,-2.00000,-1.00000";
    string expectedLats = "1.00000,2.00000,2.00000,1.00000";
    string expectedSequence = "1,2,3,4";
    string gotLons, gotLats, gotSequence;
    gring.getGringStrings(gotLons, gotLats, gotSequence);

    string expected = expectedLons + "\n        " + expectedLats + "\n        " + expectedSequence;
    string returned = gotLons + "\n        " + gotLats + "\n        " + gotSequence;
    string expectedVsReturned = "    Expected:\n        " + expected + "\n    Got:\n        " + returned;

    return {expected == returned, expectedVsReturned};
}

TestsResult testGringLonsMixed() {
    vector<double> lats = {0.0, 0.0};
    vector<double> lons = {-2.0, -1.0};

    Gring gring;
    gring.tryIncludeScan(lats.data(), lons.data(), lons.size(), 0);
    lats = {2.0, 2.0};
    lons = {0.0, 1.0};
    gring.tryIncludeScan(lats.data(), lons.data(), lons.size(), 0);

    string expectedLons = "-2.00000,0.00000,1.00000,-1.00000";
    string expectedLats = "0.00000,2.00000,2.00000,0.00000";
    string expectedSequence = "1,2,3,4";
    string gotLons, gotLats, gotSequence;
    gring.getGringStrings(gotLons, gotLats, gotSequence);

    string expected = expectedLons + "\n        " + expectedLats + "\n        " + expectedSequence;
    string returned = gotLons + "\n        " + gotLats + "\n        " + gotSequence;
    string expectedVsReturned = "    Expected:\n        " + expected + "\n    Got:\n        " + returned;

    return {expected == returned, expectedVsReturned};
}

TestsResult testWktStartingNegativeCrossing180CounterClockwiseStaysCounterClockwise() {
    vector<double> firstLons = {-179.0, 179.0};
    vector<double> firstLats = {1.0, 1.0};
    vector<double> secondLons = {-179.0, 179.0};
    vector<double> secondLats = {0.0, 0.0};

    Gring gring;
    gring.tryIncludeScan(firstLats.data(), firstLons.data(), firstLons.size(), 0);
    gring.tryIncludeScan(secondLats.data(), secondLons.data(), secondLons.size(), 0);

    // Expect BCDAB
    string expected =
        "POLYGON((-179.00000 1.00000,179.00000 1.00000,179.00000 0.00000,-179.00000 0.00000,-179.00000 "
        "1.00000))";
    string got = gring.getGeospatialBounds();
    return {expected == got, "    Expected " + expected + " \n    Got      " + got};
}

TestsResult testWktStartingPositiveCrossing180CounterClockwiseStaysCounterClockwise() {
    vector<double> firstLons = {179.0, -179.0};   // A, B
    vector<double> firstLats = {0.0, 0.0};        // A, B
    vector<double> secondLons = {179.0, -179.0};  // D, C
    vector<double> secondLats = {1.0, 1.0};       // D, C

    Gring gring;
    gring.tryIncludeScan(firstLats.data(), firstLons.data(), firstLons.size(), 0);
    gring.tryIncludeScan(secondLats.data(), secondLons.data(), secondLons.size(), 0);

    // Expect ABCDA
    string expected =
        "POLYGON((179.00000 0.00000,-179.00000 0.00000,-179.00000 1.00000,179.00000 1.00000,179.00000 "
        "0.00000))";
    string got = gring.getGeospatialBounds();
    return {expected == got, "    Expected " + expected + " \n    Got      " + got};
}

TestsResult testWktStartingNegativeCrossing180ClockwiseGoesCounterClockwise() {
    vector<double> firstLons = {-179.0, 179.0};   // B, A
    vector<double> firstLats = {1.0, 1.0};        // B, A
    vector<double> secondLons = {-179.0, 179.0};  // D, C
    vector<double> secondLats = {0.0, 0.0};       // D, C

    Gring gring;
    gring.tryIncludeScan(firstLats.data(), firstLons.data(), firstLons.size(), 0);
    gring.tryIncludeScan(secondLats.data(), secondLons.data(), secondLons.size(), 0);

    // Expect BADCB
    string expected =
        "POLYGON((-179.00000 1.00000,179.00000 1.00000,179.00000 0.00000,-179.00000 0.00000,-179.00000 "
        "1.00000))";
    string got = gring.getGeospatialBounds();
    return {expected == got, "    Expected " + expected + " \n    Got      " + got};
}

TestsResult testWktStartingPositiveCrossing180ClockwiseGoesCounterClockwise() {
    vector<double> firstLons = {179.0, -179.0};   // D, C
    vector<double> firstLats = {0.0, 0.0};        // D, C
    vector<double> secondLons = {179.0, -179.0};  // A, B
    vector<double> secondLats = {1.0, 1.0};       // A, B

    Gring gring;
    gring.tryIncludeScan(firstLats.data(), firstLons.data(), firstLons.size(), 0);
    gring.tryIncludeScan(secondLats.data(), secondLons.data(), secondLons.size(), 0);

    // Expect ABCDA
    string expected =
        "POLYGON((179.00000 0.00000,-179.00000 0.00000,-179.00000 1.00000,179.00000 1.00000,179.00000 "
        "0.00000))";
    string got = gring.getGeospatialBounds();
    return {expected == got, "    Expected " + expected + " \n    Got      " + got};
}

TestsResult testWktIntersectingPolygonProducesEmptyString() {
    /*
       Cross centered at (5, 5). This creates an intersecting polygon, and Gring should fail to create a WKT
       string
    */

    vector<double> verticalLons = {5.0, 5.0};
    vector<double> verticalLats = {0.0, 10.0};
    vector<double> horizontalLons = {5.0, 10.0};
    vector<double> horizontalLats = {5.0, 5.0};

    Gring gring;

    gring.tryIncludeScan(verticalLats.data(), verticalLons.data(), verticalLons.size(), 0);
    gring.tryIncludeScan(horizontalLats.data(), horizontalLons.data(), horizontalLons.size(), 1);

    string expected = "";  // Gring shouldn't produce a WKT string when the data intersects
    string got = gring.getGeospatialBounds();

    return {expected == got, "   Expected \"\", got " + got};
}

TestsResult testWktWithDataSpanningPoleIsCcw() {
    vector<double> latsBeforeCrossing = {76.0, 70.0};
    vector<double> latsWhileCrossing = {86.0, 80.0};
    vector<double> latsAfterCrossing = {66.0, 60.0};

    vector<double> lonsBeforeCrossing = {60.0, 90.0};
    vector<double> lonsWhileCrossing = {0.0, 180.0};
    vector<double> lonsAfterCrossing = {-60.0, -90.0};

    Gring gring;
    gring.tryIncludeScan(latsBeforeCrossing.data(), lonsBeforeCrossing.data(), lonsBeforeCrossing.size(), 0);
    gring.tryIncludeScan(latsWhileCrossing.data(), lonsWhileCrossing.data(), lonsWhileCrossing.size(), 1);
    gring.tryIncludeScan(latsAfterCrossing.data(), lonsAfterCrossing.data(), lonsWhileCrossing.size(), 2);

    string expected =
        "POLYGON((60.00000 76.00000,-60.00000 66.00000,-90.00000 60.00000,90.00000 70.00000,60.00000 "
        "76.00000))";
    string got = gring.getGeospatialBounds();

    return {expected == got, "    Expected " + expected + " \n    Got      " + got};
}

TestsResult testGringCrossingPoleIsClockwise() {
    vector<double> latsBeforeCrossing = {76.0, 70.0};
    vector<double> latsWhileCrossing = {86.0, 80.0};
    vector<double> latsAfterCrossing = {66.0, 60.0};

    vector<double> lonsBeforeCrossing = {60.0, 90.0};
    vector<double> lonsWhileCrossing = {0.0, 180.0};
    vector<double> lonsAfterCrossing = {-60.0, -90.0};

    Gring gring;
    gring.tryIncludeScan(latsBeforeCrossing.data(), lonsBeforeCrossing.data(), lonsBeforeCrossing.size(), 0);
    gring.tryIncludeScan(latsWhileCrossing.data(), lonsWhileCrossing.data(), lonsWhileCrossing.size(), 1);
    gring.tryIncludeScan(latsAfterCrossing.data(), lonsAfterCrossing.data(), lonsWhileCrossing.size(), 2);

    string expectedLons = "60.00000,-60.00000,-90.00000,90.00000";
    string expectedLats = "76.00000,66.00000,60.00000,70.00000";
    string expectedSequence = "1,2,3,4";
    string gotLons, gotLats, gotSequence;
    gring.getGringStrings(gotLons, gotLats, gotSequence);

    string expected = expectedLons + "\n        " + expectedLats + "\n        " + expectedSequence;
    string returned = gotLons + "\n        " + gotLats + "\n        " + gotSequence;
    string expectedVsReturned = "    Expected:\n        " + expected + "\n    Got:\n        " + returned;

    return {expected == returned, expectedVsReturned};
}

TestsResult testDirectionIsPositiveWhenAllScansAscend() {
    vector<double> lons = {0.0, 1.0};
    vector<double> lats1 = {0.0, 0.0};
    vector<double> lats2 = {1.0, 1.0};
    vector<double> lats3 = {2.0, 2.0};


    Gring gring(20.0, 0, 0); // Flow will be checked every scan. Other params are default
    gring.tryIncludeScan(lats1.data(), lons.data(), lons.size(), 0);
    gring.tryIncludeScan(lats2.data(), lons.data(), lons.size(), 1);
    gring.tryIncludeScan(lats3.data(), lons.data(), lons.size(), 2);

    int satelliteDirection = gring.getSatelliteDirection();

    return {satelliteDirection > 0,
            "    Satellite direction was not ascending: " + std::to_string(satelliteDirection)};
}

TestsResult testDirectionIsNegativeWhenAllScansDescend() {
    vector<double> lons = {0.0, 1.0};
    vector<double> lats1 = {2.0, 2.0};
    vector<double> lats2 = {1.0, 1.0};
    vector<double> lats3 = {0.0, 0.0};


    Gring gring(20.0, 0, 0); // Flow will be checked every scan. Other params are default
    gring.tryIncludeScan(lats1.data(), lons.data(), lons.size(), 0);
    gring.tryIncludeScan(lats2.data(), lons.data(), lons.size(), 1);
    gring.tryIncludeScan(lats3.data(), lons.data(), lons.size(), 2);

    int satelliteDirection = gring.getSatelliteDirection();

    return {satelliteDirection < 0,
            "    Satellite direction was not descending: " + std::to_string(satelliteDirection)};
}

TestsResult testDirectionIs0WhenDataIsCenteredOnAPole() {
    // This would pretty much never happen with real data, but is provided for completeness of testing
    vector<double> lons1 = {0.0, 1.0};
    vector<double> lons2 = {1.0, 1.0}; // Prevents lat/lon point duplication
    vector<double> lats1 = {89.0, 89.0};
    vector<double> lats2 = {90.0, 90.0};
    vector<double> lats3 = {89.0, 89.0};


    Gring gring(20.0, 0, 0); // Flow will be checked every scan. Other params are default
    gring.tryIncludeScan(lats1.data(), lons1.data(), lons1.size(), 0);
    gring.tryIncludeScan(lats2.data(), lons2.data(), lons2.size(), 1);
    gring.tryIncludeScan(lats3.data(), lons2.data(), lons2.size(), 2);

    int satelliteDirection = gring.getSatelliteDirection();

    return {satelliteDirection == 0,
            "    Satellite direction was not undetermined: " + std::to_string(satelliteDirection)};
}

TestsResult testWktIsEmptyBeforeProcessing() {
    Gring gring;

    string wkt = gring.getGeospatialBounds();
    return {wkt == "", "    WKT was not empty: " + wkt};
}

TestsResult testGringStringsEmptyBeforeProcessing() {
    Gring gring;

    string lats, lons, sequence;
    gring.getGringStrings(lons, lats, sequence);

    bool allStringsEmpty = lons == "" && lats == "" && sequence == "";
    return {allStringsEmpty, "    " + lons + " " + lats + " " + sequence};
}

TestsResult testExtremaAreLimitsBeforeProcessing() {
    Gring gring;

    array<double, 4> extrema = gring.getGeospatialExtremes();
    bool extremaAreLimits = 
        extrema[0] == numeric_limits<float>::max() &&
        extrema[1] == numeric_limits<float>::lowest() &&
        extrema[2] == numeric_limits<float>::max() &&
        extrema[3] == numeric_limits<float>::lowest();


    string extremaCsv = to_string(extrema[0]) + ',' + to_string(extrema[1]) + ',' + to_string(extrema[2]) + ',' + to_string(extrema[3]);
    return {extremaAreLimits, "    " + extremaCsv};
}

TestsResult testDirectionIs0BeforeProcessing() {
    Gring gring;

    int direction = gring.getSatelliteDirection();

    return {direction == 0, "    Direction was not 0: " + std::to_string(direction)};
}

TestsResult testPointIsNotIncludedIfFlagged() {
    vector<int32_t> flagMasks = {0xfffff, 0xfffff}; // This would flag all pixels
    vector<double> lons = {0.0, 1.0};
    vector<double> lats1 = {0.0, 0.0};
    vector<double> lats2 = {1.0, 1.0};
    vector<double> lats3 = {2.0, 2.0};

    Gring gring;
    gring.setInvalidPixelMask(0x02);
    gring.tryIncludeScan(lats1.data(), lons.data(), lons.size(), 0, flagMasks.data()); 
    gring.tryIncludeScan(lats2.data(), lons.data(), lons.size(), 1);
    gring.tryIncludeScan(lats3.data(), lons.data(), lons.size(), 2);

    double latMin = gring.getGeospatialLatitudeMin(); 

    return {latMin == 1.0, "     Flagged scan was included"};
}

TestsResult testPointIsNotIncludedInWktIfFlagged() {
    vector<int32_t> flagMasks = {0xfffff, 0xfffff}; // This would flag all pixels
    vector<double> lons = {0.0, 1.0};
    vector<double> lats1 = {0.0, 0.0};
    vector<double> lats2 = {1.0, 1.0};
    vector<double> lats3 = {2.0, 2.0};

    Gring gring;
    gring.setInvalidPixelMask(0x02);
    gring.tryIncludeScan(lats1.data(), lons.data(), lons.size(), 0, flagMasks.data()); 
    gring.tryIncludeScan(lats2.data(), lons.data(), lons.size(), 1);
    gring.tryIncludeScan(lats3.data(), lons.data(), lons.size(), 2);

    string got = gring.getGeospatialBounds();
    string expected =
        "POLYGON((0.00000 1.00000,1.00000 1.00000,1.00000 2.00000,0.00000 2.00000,0.00000 1.00000))";

    return {expected == got, "     Expected: " + expected + "\n     Got:     " + got};
}

TestsResult testPointIsNotIncludedInGringStringsIfFlagged() {
    vector<int32_t> flagMasks = {0xfffff, 0xfffff}; // This would flag all pixels
    vector<double> lons = {0.0, 1.0};
    vector<double> lats1 = {0.0, 0.0};
    vector<double> lats2 = {1.0, 1.0};
    vector<double> lats3 = {2.0, 2.0};

    Gring gring;
    gring.setInvalidPixelMask(0x02);
    gring.tryIncludeScan(lats1.data(), lons.data(), lons.size(), 0, flagMasks.data()); 
    gring.tryIncludeScan(lats2.data(), lons.data(), lons.size(), 1);
    gring.tryIncludeScan(lats3.data(), lons.data(), lons.size(), 2);

    string gotLats, gotLons, gotSequence;
    gring.getGringStrings(gotLons, gotLats, gotSequence);

    string expectedLats = "1.00000,2.00000,2.00000,1.00000";
    string expectedLons = "0.00000,0.00000,1.00000,1.00000";
    string expectedSequence = "1,2,3,4";

    string got = gotLats + ' ' + gotLons + ' ' + gotSequence;
    string expected = expectedLats + ' ' + expectedLons + ' ' + expectedSequence;
    string errorMessage = "    Expected: " + expected + "\n    Got:      " + got;

    return {expected == got, errorMessage};
}

TestsResult testScansWithFillGetIncludedInGring() {
    vector<double> lons =  {-32767.0, 0.0, 1.0, -32767.0};
    vector<double> lats1 = {-32767.0, 0.0, 0.0, -32767.0};
    vector<double> lats2 = {-32767.0, 1.0, 1.0, -32767.0};
    vector<double> lats3 = {-32767.0, 2.0, 2.0, -32767.0};

    Gring gring;
    gring.tryIncludeScan(lats1.data(), lons.data(), lons.size(), 0); 
    gring.tryIncludeScan(lats2.data(), lons.data(), lons.size(), 1);
    gring.tryIncludeScan(lats3.data(), lons.data(), lons.size(), 2);

    string gotLats, gotLons, gotSequence;
    gring.getGringStrings(gotLons, gotLats, gotSequence);

    string expectedLats = "0.00000,2.00000,2.00000,0.00000";
    string expectedLons = "0.00000,0.00000,1.00000,1.00000";
    string expectedSequence = "1,2,3,4";

    string got = gotLats + ' ' + gotLons + ' ' + gotSequence;
    string expected = expectedLats + ' ' + expectedLons + ' ' + expectedSequence;
    string errorMessage = "    Expected: " + expected + "\n    Got:      " + got;

    return {expected == got, errorMessage};
}

TestsResult testScansWithFillGetIncludedInWkt() {
    vector<double> lons =  {-32767.0, 0.0, 1.0, -32767.0};
    vector<double> lats1 = {-32767.0, 0.0, 0.0, -32767.0};
    vector<double> lats2 = {-32767.0, 1.0, 1.0, -32767.0};
    vector<double> lats3 = {-32767.0, 2.0, 2.0, -32767.0};

    Gring gring;
    gring.tryIncludeScan(lats1.data(), lons.data(), lons.size(), 0); 
    gring.tryIncludeScan(lats2.data(), lons.data(), lons.size(), 1);
    gring.tryIncludeScan(lats3.data(), lons.data(), lons.size(), 2);

    string expected =
        "POLYGON((0.00000 0.00000,1.00000 0.00000,1.00000 2.00000,0.00000 2.00000,0.00000 0.00000))";
    string got = gring.getGeospatialBounds();

    return {expected == got, "     Expected: " + expected + "\n     Got:      " + got};
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
        {"Test WKT string lons all positive", testWktLonsPositive},
        {"Test WKT string lons all negative", testWktLonsNegative},
        {"Test WKT string lons mixed", testWktLonsMixed},
        {"Test WKT string lats all positive", testWktLatsPositive},
        {"Test WKT string lats all negative", testWktLatsNegative},
        {"Test WKT string lats mixed", testWktLatsMixed},
        {"Test Gring strings with positive lats", testGringLatsPositive},
        {"Test Gring strings with negative lats", testGringLatsNegative},
        {"Test Gring strings with mixed-sign lats", testGringLatsMixed},
        {"Test Gring strings with positive lons", testGringLonsPositive},
        {"Test Gring strings with negative lons", testGringLonsNegative},
        {"Test Gring strings with mixed-sign lons", testGringLonsMixed},
        {"Test WKT stays CCW crossing 180 starting negative",
         testWktStartingNegativeCrossing180CounterClockwiseStaysCounterClockwise},
        {"Test WKT stays CCW crossing 180 starting positive",
         testWktStartingPositiveCrossing180CounterClockwiseStaysCounterClockwise},
        {"Test CW WKT goes CCW crossing 180 starting negative",
         testWktStartingNegativeCrossing180ClockwiseGoesCounterClockwise},
        {"Test CW WKT goes CCW crossing 180 starting positive",
         testWktStartingPositiveCrossing180ClockwiseGoesCounterClockwise},
        {"Test WKT is empty when scans intersect", testWktIntersectingPolygonProducesEmptyString},
        {"Test WKT is CCW when data spans pole", testWktWithDataSpanningPoleIsCcw},
        {"Test Gring is clockwise when data spans pole", testGringCrossingPoleIsClockwise},
        {"Test direction is positive when satellite only ascends", testDirectionIsPositiveWhenAllScansAscend},
        {"Test direction is negative when satellite only descends",
         testDirectionIsNegativeWhenAllScansDescend},
        {"Test direction is zero when lats center on a pole", testDirectionIs0WhenDataIsCenteredOnAPole},
        {"Test WKT is empty before processing", testWktIsEmptyBeforeProcessing},
        {"Test Gring strings are empty before processing", testGringStringsEmptyBeforeProcessing},
        {"Test extrema are numeric limits before processing", testExtremaAreLimitsBeforeProcessing},
        {"Test satellite direction is 0 before processing", testDirectionIs0BeforeProcessing},
        {"Test point is not included if the scan is flagged", testPointIsNotIncludedIfFlagged},
        {"Test point is not included in WKT if the scan is flagged", testPointIsNotIncludedInWktIfFlagged},
        {"Test flagged point is not included in Gring strings",
         testPointIsNotIncludedInGringStringsIfFlagged},
        {"Test Scans bracketed by -32767.0 are included in Gring", testScansWithFillGetIncludedInGring},
        {"Test Scans bracketed by -32767.0 are included in WKT", testScansWithFillGetIncludedInWkt},

    };
}

int main() {
    vector<pair<string, function<TestsResult()>>> tests = getTests();

    size_t testsPassed = 0;

    for (const pair<string, function<TestsResult()>>& test : tests) {
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
