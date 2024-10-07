#include "l1czcs.h"

int cz_dat_alloc(int nscans, int n_ctl_pt, int r_mode,
        l1_data_struc *l1_data)
/*******************************************************************

   cz_dat_alloc

   purpose: allocate the data storage in the czcs data structure

   Returns type: int - 0 - no problems

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      int               nscans           I      # scan lines of data
      int               n_ctl_pt         I      # pixel control points
      int               r_mode           I      read mode for how much of the
                                                data arrays to allocate:
                                                0 - allocate all data arrays
                                                1 - allocate time and qual only
      l1_data_struc *   l1_data         I/O     data arrays

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson, SAIC 10 Aug 2004      Original development
      W. Robinson, SAIC 19 Dec 2005      add the pos_err position error array
                                         allocation

 *******************************************************************/
 {
    int j;

    l1_data->msec = (int *) malloc(nscans * sizeof ( int));
    l1_data->cal_sum =
            (unsigned char *) malloc(nscans * 5 * sizeof ( unsigned char));
    l1_data->cal_scan =
            (unsigned char *) malloc(nscans * 6 * sizeof ( unsigned char));
    /*
     *  for the rest, only allocate if mode 0
     */
    if (r_mode == 0) {
        for (j = 0; j < 6; j++)
            l1_data->counts[j] = (unsigned char *)
            malloc(nscans * NCZCS_PIX * sizeof ( unsigned char));
        l1_data->tilt = (float *) malloc(nscans * sizeof ( float));
        l1_data->slat = (float *) malloc(nscans * sizeof ( float));
        l1_data->slon = (float *) malloc(nscans * sizeof ( float));
        l1_data->clat = (float *) malloc(nscans * sizeof ( float));
        l1_data->clon = (float *) malloc(nscans * sizeof ( float));
        l1_data->elat = (float *) malloc(nscans * sizeof ( float));
        l1_data->elon = (float *) malloc(nscans * sizeof ( float));
        l1_data->orb_vec =
                (float *) malloc(nscans * 3 * sizeof ( float));
        l1_data->att_ang =
                (float *) malloc(nscans * 3 * sizeof ( float));
        l1_data->pos_err = (float *) malloc(nscans * sizeof ( float));
        l1_data->slope =
                (float *) malloc(nscans * 6 * sizeof ( float));
        l1_data->intercept =
                (float *) malloc(nscans * 6 * sizeof ( float));
        l1_data->gain = (short *) malloc(nscans * sizeof ( short));
        l1_data->ctl_pt_rows = (int *) malloc(nscans * sizeof ( int));
        l1_data->ctl_pt_cols = (int *) malloc(n_ctl_pt * sizeof (int));
        l1_data->ctl_pt_lat =
                (float *) malloc(nscans * n_ctl_pt * sizeof ( float));
        l1_data->ctl_pt_lon =
                (float *) malloc(nscans * n_ctl_pt * sizeof ( float));
        /*
         *  for extra outputs of geometry and calibrated radiance data
         */
#ifdef GEOM_CAL
        l1_data->sen_zen =
                (float *) malloc(nscans * NCZCS_PIX * sizeof ( float));
        l1_data->sen_az =
                (float *) malloc(nscans * NCZCS_PIX * sizeof ( float));
        l1_data->sol_zen =
                (float *) malloc(nscans * NCZCS_PIX * sizeof ( float));
        l1_data->sol_az =
                (float *) malloc(nscans * NCZCS_PIX * sizeof ( float));
        l1_data->all_lat =
                (float *) malloc(nscans * NCZCS_PIX * sizeof ( float));
        l1_data->all_lon =
                (float *) malloc(nscans * NCZCS_PIX * sizeof ( float));
        l1_data->Lt_443 =
                (float *) malloc(nscans * NCZCS_PIX * sizeof ( float));
        l1_data->Lt_520 =
                (float *) malloc(nscans * NCZCS_PIX * sizeof ( float));
        l1_data->Lt_550 =
                (float *) malloc(nscans * NCZCS_PIX * sizeof ( float));
        l1_data->Lt_670 =
                (float *) malloc(nscans * NCZCS_PIX * sizeof ( float));
        l1_data->Lt_750 =
                (float *) malloc(nscans * NCZCS_PIX * sizeof ( float));
        l1_data->Lt_11500 =
                (float *) malloc(nscans * NCZCS_PIX * sizeof ( float));
#endif
    }
    /*
     *  and return
     */
    return 0;
}

void cz_dat_free(l1_data_struc *l1_data, int r_mode)
/*******************************************************************

   cz_dat_free

  purpose: easy way to fre up all the array space for the SDSes

  Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      l1_data_struc *   l1_data          I      structure for czcs data
      int               r_mode           I      read mode, as above

  Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson, SAIC 6 Aug 2004      Original development

 *******************************************************************/
 {
    int i;
    /*
     *  just free the allocated arrays
     */
    free(l1_data->msec);
    free(l1_data->cal_sum);
    free(l1_data->cal_scan);
    if (r_mode == 0) {
        for (i = 0; i < 6; i++)
            free(l1_data->counts[i]);
        free(l1_data->tilt);
        free(l1_data->slat);
        free(l1_data->slon);
        free(l1_data->clat);
        free(l1_data->clon);
        free(l1_data->elat);
        free(l1_data->elon);
        free(l1_data->orb_vec);
        free(l1_data->att_ang);
        free(l1_data->pos_err);
        free(l1_data->slope);
        free(l1_data->intercept);
        free(l1_data->gain);
        free(l1_data->ctl_pt_rows);
        free(l1_data->ctl_pt_cols);
        free(l1_data->ctl_pt_lat);
        free(l1_data->ctl_pt_lon);
        /*
         *  The extra values of satellite view, full navigation and calibrated 
         *  radiances need not be accounted for as they are for test only
         */
    }
    return;
}
