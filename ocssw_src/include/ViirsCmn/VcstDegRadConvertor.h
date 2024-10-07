/**************************************************************************
 *
 * NAME: VcstDegRadConvertor.h
 *
 * DESCRIPTION: A template class meant to convert degrees to radians and
 *              vice versa for different data types.
 *              Due to its template nature any defined types are allowed.
 *              User types must have the
 *              =, *, +, / , and == operators overloaded in order to
 *              utilize this class.
 *
 *              A separate class, VcstDegRadConvertorData, is also defined
 *              to hold static arrays used for bounds-checking performed by
 *              template class functions in VcstDegRadConvertor.
 *
 * Adapted directly from ProCmnDegRadConvertor.h published by
 * Raytheon Company.
 *
 **************************************************************************/

#ifndef VCSTDEGRADCONVERTOR_H
#define VCSTDEGRADCONVERTOR_H

#include <cmath>
#include <VcstLogger.h>
#include <iostream>
#include <sstream>

//
//--------------------------------------------------------------------------
//

//   *** W A R N I N G ***
// Do not rearrange order of enum values below without care - code
// dependencies exist.
//
// The following enum, struct, and array are for range-checking conversion results
// if selected (optional) by caller.  Default case is "NOCHECK" where no range
// checks are performed on converted data.

enum RangeCheckTypeEnum {
    NOCHECK = 0,
    CKDEG90TO90,
    CKDEG0TO90,
    CKDEG180TO180,
    CKDEG0TO180,
    CKDEG360TO360,
    CKDEG0TO360,
    CKRADPIO2TOPIO2,
    CKRAD0TOPIO2,
    CKRADPITOPI,
    CKRAD0TOPI,
    CKRAD2PITO2PI,
    CKRAD0TO2PI
};

// rangeCheckBounds: type for float upper/lower limits for bounds-checks.

typedef struct {
    float fLower; // float Lower Bound.
    float fUpper; // float Upper Bound.
} rangeCheckBounds;

// rangeCheckBounds64: type for double upper/lower limits for bounds-checks.

typedef struct {
    double bLower; // double Lower Bound.
    double bUpper; // double Upper Bound.
} rangeCheckBounds64;

//
//--------------------------------------------------------------------------
//

class VcstDegRadConvertorData {
public:

    /**
     * rangeCheck[]: array used for float specializations for bounds-checking results
     * of DEG <--> RAD conversions.
     */
    static rangeCheckBounds rangeCheck[13];

    /**
     * rangeCheck64[]: array used for double specializations for bounds-checking results
     * of DEG <--> RAD conversions.
     */
    static rangeCheckBounds64 rangeCheck64[13];

protected:

private:
    /**
     *  Make private so a object can not be instantiated.
     */
    VcstDegRadConvertorData();

    /**
     *  Destructor
     */
    ~VcstDegRadConvertorData();

};

//
//--------------------------------------------------------------------------
//
//   *** W A R N I N G ***
// Do not rearrange order of enum values below without care - code
// dependencies exist.

const int NUMCHECKS = 13; // Number of bounds checks 

// Array "rangeCheck[]" which follows is used for bounds-checking the
// results of float degree-radian conversions.

//static//
rangeCheckBounds VcstDegRadConvertorData::rangeCheck[NUMCHECKS] = {
    { 0.0, 0.0}, // NOCHECK
    { -90.0, 90.0}, // CKDEG90TO90
    { 0.0, 90.0}, // CKDEG0TO90
    { -180.0, 180.0}, // CKDEG180TO180
    { 0.0, 180.0}, // CKDEG0TO180
    { -360.0, 360.0}, // CKDEG360TO360
    { 0.0, 360.0}, // CKDEG0TO360
    { -FLOAT32_PIO2, FLOAT32_PIO2}, // CKRADPIO2TOPIO2
    { 0.0, FLOAT32_PIO2}, // CKRAD0TOPIO2
    { -FLOAT32_PI, FLOAT32_PI}, // CKRADPITOPI
    { 0.0, FLOAT32_PI}, // CKRAD0TOPI
    { -(FLOAT32_PI * 2), (FLOAT32_PI * 2)}, // CKRAD2PITO2PI
    { 0.0, (FLOAT32_PI * 2)} // CKRAD0TO2PI
};

//
//--------------------------------------------------------------------------
//
//   *** W A R N I N G ***
// Do not rearrange order of enum values below without care - code
// dependencies exist.

// Array "rangeCheck64[]" which follows is used for bounds-checking the
// results of double degree-radian conversions.

//static//
rangeCheckBounds64 VcstDegRadConvertorData::rangeCheck64[NUMCHECKS] = {
    { 0.0,
        0.0}, // NOCHECK
    { -90.0, 90.0}, // CKDEG90TO90
    { 0.0, 90.0}, // CKDEG0TO90
    { -180.0, 180.0}, // CKDEG180TO180
    { 0.0, 180.0}, // CKDEG0TO180
    { -360.0, 360.0}, // CKDEG360TO360
    { 0.0, 360.0}, // CKDEG0TO360
    { -PIO2, PIO2}, // CKRADPIO2TOPIO2
    { 0.0, PIO2}, // CKRAD0TOPIO2
    { -PI, PI}, // CKRADPITOPI
    { 0.0, PI}, // CKRAD0TOPI
    { -TWOPI, TWOPI}, // CKRAD2PITO2PI
    { 0.0, TWOPI} // CKRAD0TO2PI
};

//
//--------------------------------------------------------------------------
//

template<class _T>
class VcstDegRadConvertor {
public:

    /**
     * degToRadConversion - converts an array from degrees to radian values,
     * could be an array of 1
     */
    static int degToRadConversion(_T* in, _T* out, int count, int checkType =
            NOCHECK);

    /**
     * radToDegConversion - converts an array from radians to degree values,
     * could be an array of 1
     */
    static int radToDegConversion(_T* in, _T* out, int count, int checkType =
            NOCHECK);

protected:

private:
    /**
     *  Make private so a object can not be instantiated.
     */
    VcstDegRadConvertor();

    /**
     *  Destructor
     */
    ~VcstDegRadConvertor();

};

//------------------------------------------------------------------------------
//  Template Class Implementation
//------------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//-Constructor-

template<class _T>
VcstDegRadConvertor<_T>::VcstDegRadConvertor() {
    // Empty constructor
}

//-Destructor-

template<class _T>
VcstDegRadConvertor<_T>::~VcstDegRadConvertor() {
    // Empty destructor
}

//-----------------------------------------------------------------------------
// degToRadConversion - Takes an array of values in the degree format, 
//                      populates an output array of values in radian format
// Main template.
//-----------------------------------------------------------------------------

template<class _T>
int VcstDegRadConvertor<_T>::degToRadConversion(_T* inPtr, _T* outPtr,
        int count, int checkType) {
    // Main template
    // NOT IMPLEMENTED - DO NOT USE

    return (VCST_FAIL);
}

//-----------------------------------------------------------------------------
// degToRadConversion - Takes an array of values in the degree format, 
//                      populates an output array of values in radian format
// float specialization.
//-----------------------------------------------------------------------------

template<>
int VcstDegRadConvertor<float>::degToRadConversion(float* inPtr, float* outPtr,
        int count, int checkType) {
    int status = VCST_SUCCESS;
    bool rangeCheck = false;
    float fLower = 0.0;
    float fUpper = 0.0;

    // Validate range-checking option
    if (checkType != NOCHECK) {
        if ((checkType >= CKRADPIO2TOPIO2) && (checkType <= CKRAD0TO2PI)) {
            rangeCheck = true;
            fLower = VcstDegRadConvertorData::rangeCheck[checkType].fLower;
            fUpper = VcstDegRadConvertorData::rangeCheck[checkType].fUpper;
        } else {
            DEBUG_HIGH(VcstLogger::getLogger(),
                    "degToRadConversion failed checkType");
            return VCST_FAIL;
        }
    }

    // loop through the elements to be converted
    int nFillValues = 0;
    for (int i = 0; i < count; ++i) {
        // Error if NaN encountered.
        if (std::isnan(*inPtr)) {
            std::ostringstream oss;
            oss << "degToRadConversion detected NaN at " << (i);
            DEBUG_HIGH(VcstLogger::getLogger(), oss.str());
            return (VCST_FAIL);
        }

        // Bypass convert for fill values <= FLOAT32_FILL_TEST (-999.0)
        if (*inPtr > FLOAT32_FILL_TEST) {
            *outPtr = *inPtr * DEG2RAD;
            if (rangeCheck == true) {
                if (*outPtr < fLower) {
                    *outPtr = fLower;
                } else if (*outPtr > fUpper) {
                    *outPtr = fUpper;
                }
            }

        } else {
            *outPtr = *inPtr;
            nFillValues++;
        }
        // increment the pointers
        outPtr++;
        inPtr++;
    } // end for

    return (status);
}

//-----------------------------------------------------------------------------
// degToRadConversion - Takes an array of values in the degree format, 
//                      populates an output array of values in radian format
// double specialization.
//-----------------------------------------------------------------------------

template<>
int VcstDegRadConvertor<double>::degToRadConversion(double* inPtr,
        double* outPtr, int count, int checkType) {
    int status = VCST_SUCCESS;
    bool rangeCheck = false;
    double bLower = 0.0;
    double bUpper = 0.0;

    // Validate range-checking option
    if (checkType != NOCHECK) {
        if ((checkType >= CKRADPIO2TOPIO2) && (checkType <= CKRAD0TO2PI)) {
            rangeCheck = true;
            bLower = VcstDegRadConvertorData::rangeCheck64[checkType].bLower;
            bUpper = VcstDegRadConvertorData::rangeCheck64[checkType].bUpper;
        } else {
            DEBUG_HIGH(VcstLogger::getLogger(),
                    "degToRadConversion failed checkType");
            return VCST_FAIL;
        }
    }

    // loop through the elements to be converted
    int nFillValues = 0;
    for (int i = 0; i < count; ++i) {
        // Error if NaN encountered.
        if (std::isnan(*inPtr)) {
            std::ostringstream oss;
            oss << "degToRadConversion detected NaN at " << (i);
            DEBUG_HIGH(VcstLogger::getLogger(), oss.str());
            return (VCST_FAIL);
        }

        // Bypass convert for fill values <= FLOAT64_FILL_TEST (-999.0)
        if (*inPtr > FLOAT64_FILL_TEST) {
            *outPtr = *inPtr * DEG2RAD;
            if (rangeCheck == true) {
                if (*outPtr < bLower) {
                    *outPtr = bLower;
                } else if (*outPtr > bUpper) {
                    *outPtr = bUpper;
                }
            }
        } else {
            *outPtr = *inPtr;
            nFillValues++;
        }
        // increment the pointers
        outPtr++;
        inPtr++;
    } // end for

    return (status);
}

//-----------------------------------------------------------------------------
// radToDegConversion - Takes an array of values in the radian format, 
//                      populates an output array of values in degree format
// Main template.
//-----------------------------------------------------------------------------

template<class _T>
int VcstDegRadConvertor<_T>::radToDegConversion(_T* inPtr, _T* outPtr,
        int count, int checkType) {
    // Main template
    // NOT IMPLEMENTED - DO NOT USE

    return (VCST_FAIL);
}

//-----------------------------------------------------------------------------
// radToDegConversion - Takes an array of values in the radian format, 
//                      populates an output array of values in degree format
// float specialization.
//-----------------------------------------------------------------------------

template<>
int VcstDegRadConvertor<float>::radToDegConversion(float* inPtr, float* outPtr,
        int count, int checkType) {
    int status = VCST_SUCCESS;
    bool rangeCheck = false;
    float fLower = 0.0;
    float fUpper = 0.0;

    if (checkType != NOCHECK) {
        if ((checkType > NOCHECK) && (checkType <= CKDEG0TO360)) {
            rangeCheck = true;
            fLower = VcstDegRadConvertorData::rangeCheck[checkType].fLower;
            fUpper = VcstDegRadConvertorData::rangeCheck[checkType].fUpper;
        } else {
            DEBUG_HIGH(VcstLogger::getLogger(),
                    "radToDegConversion failed checkType");
            return VCST_FAIL;
        }
    }

    // loop through the elements to be converted
    int nFillValues = 0;
    for (int i = 0; i < count; ++i) {
        // Error if NaN encountered.
        if (std::isnan(*inPtr)) {
            std::ostringstream oss;
            oss << "radToDegConversion detected NaN at " << (i);
            DEBUG_HIGH(VcstLogger::getLogger(), oss.str());

            return (VCST_FAIL);
        }

        // Bypass convert for fill values <= FLOAT32_FILL_TEST (-999.0)
        if (*inPtr > FLOAT32_FILL_TEST) {
            *outPtr = *inPtr * RAD2DEG;
            if (rangeCheck == true) {
                if (*outPtr < fLower) {
                    *outPtr = fLower;
                } else if (*outPtr > fUpper) {
                    *outPtr = fUpper;
                }
            }
        } else {
            *outPtr = *inPtr;
            nFillValues++;
        }
        // increment the pointers
        outPtr++;
        inPtr++;
    } // end for

    return (status);
}

//-----------------------------------------------------------------------------
// radToDegConversion - Takes an array of values in the radian format, 
//                      populates an output array of values in degree format
// double specialization.
//-----------------------------------------------------------------------------

template<>
int VcstDegRadConvertor<double>::radToDegConversion(double* inPtr,
        double* outPtr, int count, int checkType) {
    int status = VCST_SUCCESS;
    bool rangeCheck = false;
    double bLower = 0.0;
    double bUpper = 0.0;

    // Validate range-checking option
    if (checkType != NOCHECK) {
        if ((checkType > NOCHECK) && (checkType <= CKDEG0TO360)) {
            rangeCheck = true;
            bLower = VcstDegRadConvertorData::rangeCheck64[checkType].bLower;
            bUpper = VcstDegRadConvertorData::rangeCheck64[checkType].bUpper;
        } else {
            DEBUG_HIGH(VcstLogger::getLogger(),
                    "radToDegConversion failed checkType");
            return VCST_FAIL;
        }
    }

    // loop through the elements to be converted
    int nFillValues = 0;
    for (int i = 0; i < count; ++i) {
        // Error if NaN encountered.
        if (std::isnan(*inPtr)) {
            std::ostringstream oss;
            oss << "radToDegConversion detected NaN at " << (i);
            DEBUG_HIGH(VcstLogger::getLogger(), oss.str());
            return (VCST_FAIL);
        }

        // Bypass convert for fill values <= FLOAT64_FILL_TEST (-999.0)
        if (*inPtr > FLOAT64_FILL_TEST) {
            *outPtr = *inPtr * RAD2DEG;
            if (rangeCheck == true) {
                if (*outPtr < bLower) {
                    *outPtr = bLower;
                } else if (*outPtr > bUpper) {
                    *outPtr = bUpper;
                }
            }
        } else {
            *outPtr = *inPtr;
            nFillValues++;
        }
        // increment the pointers
        outPtr++;
        inPtr++;
    } // end for

    return (status);
}

#endif
