/*******************************************************************************
 *
 * NAME:  VcstMath.h
 *
 * DESCRIPTION: This is the header file for the VcstMath object class.
 *
 * Adapted directly from the IDPS file ProCmnMath.h published
 * by Raytheon Company.
 *
 *******************************************************************************/

#ifndef VcstMath_h
#define VcstMath_h

#include <math.h>
#include <VcstCmnConsts.h>

/**
 * The VcstMath class contains general math functions.
 */

class VcstMath {
public:

    /**
     * This function computes the magnitude of a vector.
     *
     * @param   aVector   Input vector
     *
     * @return  double
     */
    static double calculateMag(const double aVector[VEC_SIZE]);

    /**
     * This function computes the magnitude of a quaternion.
     *
     * @param   aQuat     Input quaternion
     *
     * @return  double
     */
    static double calcQuatMag(const double aQuat[QUAT_SIZE]);

    /**
     * This function produces the cross product of 2 vectors.
     *
     * aUvector      Input vector
     * aVvector      Input vector
     * anOutVector   Output cross product
     *
     * @return PRO_SUCCESS or an error code
     */
    static void calculateCross(const double aUvector[VEC_SIZE],
            const double aVvector[VEC_SIZE], double anOutVector[VEC_SIZE]);

    /**
     * This function calculates the dot product of 2 vectors.
     *
     * @param  aUvector   Input vector
     * @param  aVvector   Input vector
     *
     * @return double
     */
    static double calculateDot(const double aUvector[VEC_SIZE],
            const double aVvector[VEC_SIZE]);

    /**
     * This function calculates the dot product of 2 quaternions.
     *
     * @param  aUquat     Input quaternion
     * @param  aVquat     Input quaternion
     *
     * @return double
     */
    static double calculateDotQuat(const double aUquat[QUAT_SIZE],
            const double aVquat[QUAT_SIZE]);

    /**
     * This function will apply a scalar value to each component
     * of a vector.
     *
     * @param aScalar  Input scalar value
     * @param aVector  Input/Output vector with scalar applied
     *
     * @return void
     */
    static void applyScalar(const double aScalar, double aVector[VEC_SIZE]);

    /**
     * This function will apply a scalar value to each component
     * of a quaternion.
     *
     * @param aScalar  Input scalar value
     * @param aQuat    Input/Output quaternion with scalar applied
     *
     * @return void
     */
    static void applyScalarQuat(const double aScalar, double aQuat[QUAT_SIZE]);

    /**
     * This function calculates the conjugate of a quaternion.
     *
     * The conjugate of a quaternion Q, is Q* which has the
     * property that Q X Q* = 1 using quaternion multiplication.
     * Where 1 is the identity quaternion, or the real value 1.0.
     * In quaternion notation (q1,q2,q3,q4) this is (0,0,0,1).
     * The conjugate of a quaternion Q, is formed by reversing
     * the signs of the vector part and not the scalar, so:
     *    Q* = (-q1, -q2, -q3, q4)
     *
     * @param aQuat      Input quaternion Q
     * @param outQuat    Output quaternion conjugate Q*
     *
     * @return  Void
     */
    static void conjugateQuat(const double aQuat[QUAT_SIZE],
            double outQuat[QUAT_SIZE]);

    /**
     * This function will multiply two quaternions together.
     *
     * q1        Input quaternion
     * q2        Input quaternion
     * outQuat   Output quaternion product
     *
     * @return PRO_SUCCESS or an error code
     */
    static void quatMultiply(const double q1[QUAT_SIZE],
            const double q2[QUAT_SIZE], double outQuat[QUAT_SIZE]);

    /**
     * This alternate function will multiply two quaternions together
     * using common matrix notation.
     *
     * q1        Input quaternion
     * q2        Input quaternion
     * outQuat   Output quaternion product
     *
     * @return PRO_SUCCESS or an error code
     */
    static void quatMultiplyMatrix(const double q1[QUAT_SIZE],
            const double q2[QUAT_SIZE], double outQuat[QUAT_SIZE]);

    /**
     * This function will multiply two matrices together.
     *
     * @param aM1        Input matrix
     * @param aM2        Input matrix
     * @param outMatrix  Output product of the 2 input matrices
     *
     * @return void
     */
    static void matrixMultiply(const double aM1[VEC_SIZE][VEC_SIZE],
            const double aM2[VEC_SIZE][VEC_SIZE],
            double outMatrix[VEC_SIZE][VEC_SIZE]);

    /**
     * This function will multiply two matrices (MByN NByP) together.
     *
     * @param dimM       Input Dimension of M.
     * @param dimN       Input Dimension of N.
     * @param dimP       Input Dimension of P.
     * @param aMByN      Input matrix
     * @param aNByP      Input matrix
     * @param outMatrix  Output product of the 2 input matrices
     *
     * @return void
     */
    static void matrixMultiply(const int dimM, const int dimN, const int dimP,
            const double *aMByN, const double *aNByP, double *outMatrix);

    /**
     * This function will transpose a matrix.
     *
     * @param inMatrix  Input matrix
     * @param outMatrix Output transpose of input matrix
     *
     * @return void
     */
    static void transposeMatrix(const double inMatrix[VEC_SIZE][VEC_SIZE],
            double outMatrix[VEC_SIZE][VEC_SIZE]);

    /**
     * This function will transpose a quaternion.
     * NOTE: for quaternions, the transpose is the inverse.
     * For the quaternion Q = q1 q2 q3 q4,
     * the inverse is -Q = -q1 -q2 -q3 -q4
     *
     * @param inQuat   Input quaternion
     * @param outQuat  Output transpose of input quaternion
     *
     * @return void
     */
    static void transposeQuat(const double inQuat[QUAT_SIZE],
            double outQuat[QUAT_SIZE]);

    /**
     * This function will print the contents of a matrix.
     * This method is for debugging and unit test only.
     *
     * @param aMatrix    Input matrix
     *
     * @return  Void
     */
    static void printMatrixCOut(const double aMatrix[VEC_SIZE][VEC_SIZE]);

    /**
     * This function will print the contents of a vector.
     * This method is for debugging and unit test only.
     *
     * @param aVector    Input vector
     *
     * @return  Void
     */
    static void printVectorCOut(const double aVector[VEC_SIZE]);

    /**
     * This function will print the contents of a quaternion.
     * This method is for debugging and unit test only.
     *
     * @param aQuat      Input quaternion
     *
     * @return  Void
     */
    static void printQuaternionCOut(const double aQuat[QUAT_SIZE]);

    /**
     * This function multiplies a matrix times a vector.
     *
     * @param aMatrix     Input matrix
     * @param aVector     Input vector
     *
     * @return     Void
     */
    static void matrixVectorProduct(const double aMatrix[VEC_SIZE][VEC_SIZE],
            const double aVector[VEC_SIZE], double outVec[VEC_SIZE]);

    /**
     * This function creates a unit vector.
     *
     * @param aVector      Input vector
     * @param unitVector   Output unit vector
     *
     * @return    Void
     */
    static void calculateUnitVector(const double aVector[VEC_SIZE],
            double unitVector[VEC_SIZE]);

    /**
     * This function normalizes a quaternion.
     * This is the analog to calculateUnitVector() for quaternions.
     * By definition, |Q| = sqrt(q1*q1 + q2*q2 + q3*q3 + q4*q4) = 1.
     *
     * @param aQuat        Input/Output quaternion
     *
     * @return    Void
     */
    static void normalizeQuat(double aQuat[QUAT_SIZE]);

    /**
     * This function creates the inverse of a 3x3 matrix.
     *
     * @param matA      Input matrix
     * @param matR      Output inverse of input matrix
     *
     * @return    Void
     */
    static void inverseMatrix(const double matA[VEC_SIZE][VEC_SIZE],
            double matR[VEC_SIZE][VEC_SIZE]);

    /**
     * This function creates a rotation matrix through the
     * specified angle about the specified axis.
     *
     * @param angle     Input angle of positive rotation
     * @param axis      Input axis of positive rotation
     * @param matR      Output rotation matrix
     *
     * @return    Void
     */
    static void rotationMatrix(const double angle, const short axis,
            double outMatrix[VEC_SIZE][VEC_SIZE]);

    /**
     * This function calculates the determinant of a matrix.
     *
     * @param matA     Input matrix
     *
     * @return         double
     */
    static double calculateDeterminant(const double matA[VEC_SIZE][VEC_SIZE]);

    /**
     * This function calculates the cofactor of a matrix.
     *
     * @param matA      Input matrix
     * @param cof       Output cofactor matrix
     *
     * @return  Void
     */
    static void cofactorMatrix(const double matA[VEC_SIZE][VEC_SIZE],
            double cof[VEC_SIZE][VEC_SIZE]);

};

#endif
