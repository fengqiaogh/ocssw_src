#include "l12_proto.h"

void nlw_outband(int32_t evalmask, int32_t sensorID, float wave[], int32_t nwave, float Lw[], float nLw[], float outband_correction[],float derv_f_nlw[], float *dratio) {
    static float *a0;
    static float *a1;
    static float *a2;
    static int32_t ib1;
    static int32_t ib2;
    static int firstCall = 1;
    float derv_ratio_nlw1, derv_ratio_nlw2;

    int32_t ib;
    float ratio;
    float f;

    if (firstCall) {
        firstCall = 0;
        rdsensorinfo(sensorID, evalmask, "ooblw01", (void **) &a0);
        rdsensorinfo(sensorID, evalmask, "ooblw02", (void **) &a1);
        rdsensorinfo(sensorID, evalmask, "ooblw03", (void **) &a2);
        if (sensorID == CZCS)
            ib1 = windex(443., wave, nwave);
        else
            ib1 = windex(490., wave, nwave);
        ib2 = windex(550., wave, nwave); /* not 555 for HMODIS */
    }

    if (nLw[ib1] > 0.0 && nLw[ib2] > 0.0) {

        ratio = nLw[ib1] / nLw[ib2];

        if(derv_f_nlw){
            derv_ratio_nlw1=1/nLw[ib2];
            derv_ratio_nlw2=-ratio/nLw[ib2];
            f=pow(derv_ratio_nlw1,2.)*derv_f_nlw[ib1*nwave+ib1];
            f+=pow(derv_ratio_nlw2,2.)*derv_f_nlw[ib2*nwave+ib2];
            f+=2*derv_ratio_nlw1*derv_ratio_nlw2*derv_f_nlw[ib1*nwave+ib2];
            *dratio=f;
        }

        /* limit needs to be wavelength-specific
         if ((evalmask & NEWOOB) > 0) 
             ratio = MIN(nLw[ib1]/nLw[ib2],4.0);
         else 
             ratio = nLw[ib1]/nLw[ib2];
         */

        // should be ratio < 0 = 0, ratio > 10 = 10 ???

        if (ratio > 0.0 && ratio <= 10.0) {
            for (ib = 0; ib < nwave; ib++) {
                f = (a2[ib] * ratio + a1[ib]) * ratio + a0[ib];
                nLw[ib] = nLw[ib] * f;
                Lw [ib] = Lw [ib] * f;
                outband_correction[ib] = f;

                if(derv_f_nlw){
                    f=2*a2[ib]*ratio+a1[ib];
                    derv_f_nlw[ib]=f;
                    //derv_f_nlw[ib*nwave+ib1]=f*derv_ratio_nlw1;
                   // derv_f_nlw[ib*nwave+ib2]=f*derv_ratio_nlw2;
                }
            }
        } else {
            if (derv_f_nlw)
                *dratio = 0.;
        }
    }

}
