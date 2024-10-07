#include <math.h>

/* Returns SeaWiFS water-vapor transmittance */
float water_vapor(int iband, float uH2O, float airmass) {
    static float a[4] = {-8.62884, -6.94310, -5.81033, -5.51066};
    static float b[4] = {0.766159, 0.813607, 0.617758, 0.678041};

    int i = iband - 4;

    if (i < 0 || i > 3)
        return (1.0);
    else
        return (exp(-exp(a[i] + b[i] * log(airmass * uH2O))));
}

void water_vapor_(int *iband, float *uH2O, float *airmass, float *t_h2o) {
    *t_h2o = water_vapor(*iband, *uH2O, *airmass);
    return;
}


/*
From: Jacques Descloitres <jack@modland.nascom.nasa.gov>

I used Anne Vermeulen's computations to produce a fit of the water vapor 
transmittance. She once used 6S to compute the gaseous transmittance in many 
different conditions, using the exact SeaWiFS bands.
These computations use the temperature and pressure vertical profiles of a 
US62 atmospheric profile.
These computations give a comprehensive set of transmittance values in all 
SeaWiFS bands for a whole range of viewing conditions and water vapor content.

Based on these simulations, I came up with the following expression to 
approximate the transmittance:
tH2O=exp(-exp(a+b*log(m*uH2O))),
where m is the air mass: 1/mus+1/muv
      uH2O is the total water vapor content in g/cm2

Band 5: a=-8.62884  b=0.766159  RMSE=9.21441e-06
Band 6: a=-6.9431   b=0.813607  RMSE=7.44725e-05
Band 7: a=-5.81033  b=0.617758  RMSE=9.91917e-05
Band 8: a=-5.51066  b=0.678041  RMSE=0.000112495

Note that this approximation is close to what is used for the MODIS land 
atmospheric correction.

As you can see the root mean square error is lower than 1e-4, which is much 
more accurate than necessary.
 */
