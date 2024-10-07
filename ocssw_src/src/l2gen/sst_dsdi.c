#include <netcdf.h>
#include "l1.h"
#include "l2_struc.h"
#include "math.h"
#include "sst_dsdi.h"
#include "l12_proto.h"


float dust_correction(float dustExtinction, float csenz, float BT39, float BT85, float BT11, float BT12, int32_t sensorID) {

    size_t ncoefficients = 13;
    int32_t i;
    float dsdi;
    float correction = 0.0;
    static float* coef;
    static float dsdiThreshold = 0.8;
    static float dustAODThreshold = 0.05;
    if (!coef) {
        coef = (float *) calloc(ncoefficients, sizeof (float));
        if (strlen(input->dsdicoeffile)) {
            int32_t coeffDimID, coefficientID;
            int32_t ncid;
            if (nc_open(input->dsdicoeffile, NC_NOWRITE, &ncid) == NC_NOERR) {
                printf("Loading SST DSDI coefficients from %s:\n", input->dsdicoeffile);

                if (nc_inq_dimid(ncid, "coefficient", &coeffDimID) != NC_NOERR) {
                    printf("Whoops! something is wrong reading the SST DSDI coefficient file: %s\n", input->dsdicoeffile);
                    exit(EXIT_FAILURE);
                }
                nc_inq_dimlen(ncid, coeffDimID, &ncoefficients);

                if (nc_inq_varid(ncid, "coefficients", &coefficientID) == NC_NOERR) {
                    nc_get_var_float(ncid, coefficientID, coef);
                }
                // Get DSDI threshold
                if (nc_get_att_float(ncid, NC_GLOBAL, "dsdi_threshold", &dsdiThreshold) != NC_NOERR) {
                    printf("Whoops! something is wrong reading the SST DSDI threshold attribute from file: %s\n", input->dsdicoeffile);
                    exit(EXIT_FAILURE);
                }
                // Get DSDI threshold
                if (nc_get_att_float(ncid, NC_GLOBAL, "dust_threshold", &dustAODThreshold) != NC_NOERR) {
                    printf("Did not find dust AOD threshold attribute from file: %s\n", input->dsdicoeffile);
                    printf("using default of 0.025\n");
                }
                nc_close(ncid);
                // print coefficients to be used
                printf("DSDI coefficients:");
                for (i = 0; i < ncoefficients; i++) {
                    printf("%9.6f ", coef[i]);
                }
                printf("\nDSDI threshold: %f\n",dsdiThreshold);
                printf("\nDSDI Dust AOD threshold: %f\n",dsdiThreshold);
            } else {
                printf("\nCannot load DSDI coefficients file: %s\n", input->dsdicoeffile);
                exit(EXIT_FAILURE);
            }
        }
    }

    float secantTheta = 1. / csenz;
    float S0 = secantTheta - 1.;
    dustExtinction *= secantTheta;

    dsdi = coef[0] + (coef[1] + coef[2] * S0) * (BT39 - BT12) +
            (coef[3] + coef[4] * S0) * (BT39 - BT85) +
            (coef[5] + coef[6] * S0) * (BT11 - BT12) +
            (coef[7] + coef[8] * S0) * pow(BT11 - BT12, 2) +
            (coef[9] * sqrtf(dustExtinction) + coef[10]);

    if (dustExtinction > dustAODThreshold && dsdi > dsdiThreshold) {
//        printf("dsdi: %f dust: %f S: %f csenz: %f\n",dsdi, dustExtinction,secantTheta,csenz);
        correction = coef[11] * dsdi + coef[12];
    }

    return correction;
}
