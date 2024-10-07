#include "viirs_sim_sdr.h"
#include <math.h>
#define H_CON 6.6262e-34  /* Plank's const J s^-1 */
#define C_CON 2.99792458e8  /* speed of light m s^-1 */
#define K_CON 1.3807e-23  /* Boltzman's constant J K^-1 */
#define C1 2. * H_CON * pow( C_CON, 2. ) * 1.e24
#define C2 H_CON * 1.e6 * C_CON / K_CON

float bbt_2_rad(float bbt, float lam)
/*-----------------------------------------------------------------------------
    Program:   bbt_2_rad.c

    Description: convert a BBT into a radiance for VIIRS use

    Returns, radiance in W m^-2 sec^-1 sr^-1

    Arguments:
        Type      Name         I/O   Description
        ----      ----         ---   -----------
        float     bbt           I    brightness temperature in degrees K
        float     lam           I    wavelength in microns = nm / 10^3

    Modification history:

    W. Robinson, SAIC   8 Jul 2010  Original development

----------------------------------------------------------------------------*/ {
    float rad;
    /*
     *  perform computation
     */
    rad = C1 / (pow(lam, 5.) * (exp(C2 / (lam * bbt)) - 1.));

    return rad;
}
