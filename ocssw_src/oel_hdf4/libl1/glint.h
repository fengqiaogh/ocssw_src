#include "math.h"
/**
 * @brief Get the glint iqu object
 * 
 * @param senz sensor zenith angle
 * @param solz solar zenith angle
 * @param raz relative azimuth angle, Gordon's definition, (sena - sola - PI)
 * @param ws wind speed
 * @param chi wind direction
 * @param glint_coef glitter coefficient
 * @param glint_coef_q Q stokes component
 * @param glint_coef_u U stokes compoment
 */
void getglint_iqu(float senz, float solz, float raz, float ws, float zero, float *glint_coef, float *glint_coef_q,
             float *glint_coef_u,double nw);