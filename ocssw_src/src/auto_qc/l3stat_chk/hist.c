#include <math.h>
#include <libgen.h>
#include <stdio.h>
#include <string.h>
#include "hist_proto.h"

#define NPARM 12
#define NBIN 100

static float hmin[NPARM] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.8, 0.};
static float hmax[NPARM] = {5., 4., 4., 3., 2., 2., 1., 1., 0.5, 20.,
    1.3, 0.5};
static float step[NPARM];
static int nhi[NPARM], nlo[NPARM], nbin[NPARM][NBIN];
char bfile[412];

void h_init(char *file)
/*******************************************************************

   h_init

   purpose: initialize the histogram accumulation

   Returns type: void none

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      char *            file             I      level-3 file name
                                                we'll add _bhist to 
                                                make binary output file

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       26-Nov-1997     Original development

 *******************************************************************/ {
    int i, j;
    /*
     *  initialize the bin arrays
     */
    for (i = 0; i < NPARM; i++) {
        nhi[i] = 0;
        nlo[i] = 0;
        step[i] = (hmax[i] - hmin[i]) / NBIN;
        for (j = 0; j < NBIN; j++) {
            nbin[i][j] = 0.;
        }
    }
    /*
     * create the name of the binary histogram file
     */
    /*  for tests
     strcpy( bfile, basename( file ) );
     */
    strcpy(bfile, file);
    strcat(bfile, "_bhist");
}

void h_accum(float value, int parm_num)
/*******************************************************************

   h_accum

   purpose: accumulate the histogram parameters
   Returns type: void none

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      float             value            I      value to accumulate
      int               parm_num         I      number of the parameter that 
                                                is being accumulated
   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       26-Nov-1997     Original development

 *******************************************************************/ {
    int bin;

    if (value < hmin[parm_num]) {
        nlo[parm_num]++;
    } else if (value > hmax[parm_num]) {
        nhi[parm_num]++;
    } else {
        /*
         * find the bin that this value is added to
         */
        bin = (value - hmin[parm_num]) / step[parm_num];
        if (bin > (NBIN - 1)) bin = NBIN - 1;

        nbin[parm_num][bin]++;
    }
}

void h_out(int nbins, double *sum, double *sumsq, double *loparm, double *hiparm)
/*******************************************************************

   h_out

   purpose: output the results 
   Returns type: void none

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      int               nbins            I      # observations
      double *          sum              I      sum for 12 parameters
      double *          sumsq            I      sum of squares for 12 params
      double *          loparm           I      low parameter for 12 params
      double *          hiparm           I      high parameter for 12 params
                                                is being accumulated
   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       26-Nov-1997     Original development

 *******************************************************************/ {
    FILE *stream;
    int bins, i;
    float mean[NPARM], sd[NPARM], lo[NPARM], hi[NPARM];
    /*
     *  open the file
     */
    if ((stream = fopen(bfile, "w")) == NULL) {
        printf("\n\nNOTE, failure to open binary histogram file\n'%s'\n",
                bfile);
        printf("No action will be taken but this is a problem\n");
    } else {
        bins = NBIN;
        /*
         *  write the # bins and the file name in first 416 bytes
         */
        fwrite(&bins, sizeof ( int), 1, stream);
        fwrite(bfile, sizeof ( char), 412, stream);

        /*
         *  then, write the parameter info
         */
        fwrite(hmin, sizeof ( float), NPARM, stream);
        fwrite(hmax, sizeof ( float), NPARM, stream);
        fwrite(nlo, sizeof ( int), NPARM, stream);
        fwrite(nhi, sizeof ( int), NPARM, stream);
        fwrite(nbin, sizeof ( int), (NPARM * NBIN), stream);
        /*
         *  compute the mean, sd
         */
        for (i = 0; i < NPARM; i++) {
            if (nbins > 0) {
                mean[i] = sum[i] / nbins;
                sd[i] = sumsq[i] / nbins - mean[i] * mean[i];
                if (sd[i] > 0.) {
                    sd[i] = sqrt(sd[i]);
                } else
                    sd[i] = 0;
            } else {
                mean[i] = 0;
                sd[i] = 0;
            }

            lo[i] = loparm[i];
            hi[i] = hiparm[i];
        }
        fwrite(mean, sizeof ( float), NPARM, stream);
        fwrite(sd, sizeof ( float), NPARM, stream);
        fwrite(lo, sizeof ( float), NPARM, stream);
        fwrite(hi, sizeof ( float), NPARM, stream);
        /*
         *  and close the file
         */
        fclose(stream);
    }
}
