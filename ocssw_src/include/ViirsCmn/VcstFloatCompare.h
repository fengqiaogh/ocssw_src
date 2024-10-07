/**************************************************************************
 * 
 * NAME: VcstFloatCompare
 *
 * DESCRIPTION: Compare two floating point values for equality.  Two param
 *              method uses system defined epsilon.  Three param method uses
 *              a user provided error delta value.  This class addresses
 *              SSPM E.4.1.3.5 Never compare 2 float or double values using
 *              the == operator.  Class is templated on type T, being either
 *              Float32 (float) or Float64 (double).  Typedefs are provided
 *              for convenience.
 *
 * Adapted directly from the IDPS file FloatCompare.h published
 * by Raytheon Company
 *
 **************************************************************************/

#ifndef FloatCompare_h
#define FloatCompare_h

#include <limits>
#include <cmath>

template<typename T>
class FloatCompare {
public:

    static bool equal(T firstValue, T secondValue);
    static bool equal(T firstValue, T secondValue, T errorDelta);

    static bool notequal(T firstValue, T secondValue);
    static bool notequal(T firstValue, T secondValue, T errorDelta);
};

typedef FloatCompare<float> Float32Compare;
typedef FloatCompare<double> Float64Compare;

template<typename T>
inline bool FloatCompare<T>::equal(T firstValue, T secondValue) {
    const T EPSILON = std::numeric_limits < T > ::epsilon();

    if ((fabs(firstValue - secondValue)) <= EPSILON)
        return true;
    else
        return false;
}

template<typename T>
inline bool FloatCompare<T>::equal(T firstValue, T secondValue, T errorDelta) {
    if ((fabs(firstValue - secondValue)) <= errorDelta)
        return true;
    else
        return false;
}

template<typename T>
inline bool FloatCompare<T>::notequal(T firstValue, T secondValue) {
    return !equal(firstValue, secondValue);
}

template<typename T>
inline bool FloatCompare<T>::notequal(T firstValue, T secondValue,
        T errorDelta) {
    return !equal(firstValue, secondValue, errorDelta);
}

#endif
