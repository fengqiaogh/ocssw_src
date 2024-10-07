#include "l1stat.h"
#include "l1stat_proto.h"

extern int32 stat_status;

int rd_nav_t(int32 sdfid, int32 nscans, nav_t_str *nav_t)
/*******************************************************************

   rd_nav_t

   purpose: read a set of navigation and time values from the 
            L1 dataset

   Returns type: none

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      int32             sdfid           I       SD interface ID
      int32             nscans          I       # of scan lines
      nav_t_str *       nav_t           O       structure containing the 
                                                pointers to the nav sds data

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson, GSC  29-Dec-1999     Original development

 *******************************************************************/ {
    int32 start[3], edge[3];
    /*
     *  Allocate space for the navigation buffers
     */

    if ((nav_t->orb_vec = (float32 *) calloc(nscans * 3, sizeof (float32))) == NULL) {
        sprintf(err_msg, "chk_nav: cannot allocate memory for reading orb_vec");
        return FAIL;
    }

    if ((nav_t->sun_ref = (float32 *) calloc(nscans * 3, sizeof (float32))) == NULL) {
        sprintf(err_msg, "chk_nav: cannot allocate memory for reading sun_ref");
        return FAIL;
    }

    if ((nav_t->sen_mat = (float32 *) calloc(nscans * 3 * 3, sizeof (float32))) == NULL) {
        sprintf(err_msg, "chk_nav: cannot allocate memory for reading sen_mat");
        return FAIL;
    }

    if ((nav_t->scan_ell = (float32 *) calloc(nscans * 6, sizeof (float32))) == NULL) {
        sprintf(err_msg, "chk_nav: cannot allocate memory for reading scan_ell");
        return FAIL;
    }

    if ((nav_t->msec = (int32_t *) calloc(nscans, sizeof (int32_t))) == NULL) {
        sprintf(err_msg, "chk_nav: cannot allocate memory for reading msec");
        return FAIL;
    }
    if ((nav_t->nflag = (int32 *) calloc(nscans * 8, sizeof (int32))) == NULL) {
        sprintf(err_msg, "chk_nav: cannot allocate memory for reading nflag");
        return FAIL;
    }

    /*
     *  Set start and end dims for reading orb_vec sds
     */
    start[0] = start[1] = start[2] = 0;
    edge[0] = nscans;
    edge[1] = 3;
    edge[2] = 0;
    if ((rdslice(sdfid, ORBVEC, start, edge, (VOIDP) nav_t->orb_vec)) < 0) {
        stat_status = stat_status | 1;
        return FAIL;
    }

    /*
     *  Set start and end dims and read sun_ref sds
     */
    start[0] = start[1] = start[2] = 0;
    edge[0] = nscans;
    edge[1] = 3;
    edge[2] = 0;
    if ((rdslice(sdfid, SUNREF, start, edge, (VOIDP) nav_t->sun_ref)) < 0) {
        stat_status = stat_status | 1;
        return FAIL;
    }

    /*
     *  Set start and end dims and read sen_mat sds
     */
    start[0] = start[1] = start[2] = 0;
    edge[0] = nscans;
    edge[1] = 3;
    edge[2] = 3;
    if ((rdslice(sdfid, SENMAT, start, edge, (VOIDP) nav_t->sen_mat)) < 0) {
        stat_status = stat_status | 1;
        return FAIL;
    }

    /*
     *  Set start and end dims and read scan_ell sds
     */
    start[0] = start[1] = start[2] = 0;
    edge[0] = nscans;
    edge[1] = 6;
    edge[2] = 0;
    if ((rdslice(sdfid, SCANELL, start, edge, (VOIDP) nav_t->scan_ell)) < 0) {
        stat_status = stat_status | 1;
        return FAIL;
    }

    /*
     *  Set start and end dims and read msec sds
     */
    start[0] = start[1] = start[2] = 0;
    edge[0] = nscans;
    edge[1] = 0;
    edge[2] = 0;
    if ((rdslice(sdfid, "msec", start, edge, (VOIDP) nav_t->msec)) < 0) {
        stat_status = stat_status | 1;
        return FAIL;
    }

    /*
     *  Set start and end dims and read nflag sds
     */
    start[0] = start[1] = start[2] = 0;
    edge[0] = nscans;
    edge[1] = 8;
    edge[2] = 0;
    if ((rdslice(sdfid, NFLAG, start, edge, (VOIDP) nav_t->nflag)) < 0) {
        stat_status = stat_status | 1;
        return FAIL;
    }
    return SUCCEED;
}
