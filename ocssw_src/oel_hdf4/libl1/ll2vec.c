#include <math.h>
#include "l1.h"

static float raddeg = 180. / PI;

int ll2vec(float *ll, float *vec)
/*******************************************************************

   ll2vec

   purpose: convert lat, lon (or az, ele) into vectors in x, y, z

   Returns type: int - 0 if all is OK

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      float *           ll               I      latitude longitude pair
      float *           vec              O      x, y, z vector with x
                                                along the 0 lat, lon,
                                                y at 0 lat, 90 lon and z
                                                at 90 lat
                  (assumed to exist - array passed in)

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       10 Sep 2010     Original development

 *******************************************************************/ {
    int32_t err;
    float latr, lonr;

    if ((*ll > 90.) || (*ll < -90.)) {
        *vec = -999.;
        *(vec + 1) = -999.;
        *(vec + 2) = -999.;
        err = -1;
    } else {
        latr = *ll / raddeg;
        lonr = *(ll + 1) / raddeg;
        /*
         *  get the components
         */
        *vec = cos(lonr) * cos(latr);
        *(vec + 1) = sin(lonr) * cos(latr);
        *(vec + 2) = sin(latr);
        err = 0;
    }
    return err;
}

int vec2ll(float *vec, float *ll)
/*******************************************************************

   vec2ll

   purpose: convert a vectors in x, y, z into lat, lon (or az, ele)

   Returns type: int - 0 if all is OK

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      float *           vec              I      x, y, z vector with x
                                                along the 0 lat, lon,
                                                y at 0 lat, 90 lon and z
                                                at 90 lat
      float *           ll               O      latitude longitude pair
                  (assumed to exist - array passed in)

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       10 Sep 2010     Original development

 *******************************************************************/ {
    int32_t err;
    float lonr, latr, vls, vlen;

    vls = pow(*vec, 2.) + pow(*(vec + 1), 2.) + pow(*(vec + 2), 2.);
    if (vls <= 0) {
        *ll = -999.;
        *(ll + 1) = -999.;
        err = -1;
    } else {
        /*
         *  get the vector lengths for normalization and derive lat, lon
         */
        vlen = sqrt(vls);
        latr = asin(*(vec + 2) / vlen);
        lonr = atan2(*(vec + 1) / vlen, *vec / vlen);

        *ll = latr * raddeg;
        *(ll + 1) = lonr * raddeg;
        err = 0;
    }
    return err;
}
