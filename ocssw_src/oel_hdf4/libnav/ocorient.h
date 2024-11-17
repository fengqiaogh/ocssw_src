#ifndef OCSSW_OCORIENT_H
#define OCSSW_OCORIENT_H

/**
 * c  Computes matrix product of 3x3 matrices xm1 and xm2
 * @param xm1 R*4      I      Input Matrix
 * @param xm2 R*4      I      Input Matrix
 * @param xm3 R*4      O      Output Matrix
 */
void matmpy(const float xm1[3][3], const float xm2[3][3], float xm3[3][3]);


/**
 * @brief   c  Computes coordinate transformation matrix corresponding to Euler
            c  sequence.  The order of angles in the input array is yaw, roll,
            c  pitch; according to OSC, the order of the rotations is the reverse
            c  of this; the roll and pitch angles are about the negative Y and Z
            c  axes, respectively, while the yaw angle is about the positive X axis.
            c  Reference:  Wertz, Appendix E; OSC, personal communication
            c  Modification history:
            c
            c
            c  Modified to change order of rotations to pitch, roll, yaw (-Z, -Y, -X)
            c  F. S. Patt, GSC, September 29, 1996.
            c
            c  Removed negative signs on Y and Z rotations for OCTS; order of rotations
            c  is yaw, pitch, roll (Z, Y, X)
 * @param a R*4      I      Input Array of Euler Angles (degrees)
 * @param xm  R*4      O      Output Transformation Matrix
 */
void oceuler(float a[3], float xm[3][3]);

#endif  // OCSSW_OCORIENT_H
