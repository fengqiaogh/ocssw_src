/*-----------------------------------------------------------------------------
    Function: stray_light_gac

    Returns: int32 (status)
        Returns as 0, "not done", to indicate that the routine did not 
        return any useful information for the last call, but was filling
        its buffers with scan line data; and 1, "done", to indicate that 
        the l1b_data and sl_flag arrays for scan line sl_scan have been 
        set.

    Description:
        The function stray_light_gac processes GAC scan lines by detecting
        bright targets (BTs) and labelling pixels as being within the BTs 
        or adjacent to them where they are likely to be contaminated by
        stray light (ST) from the BTs.  The Level-1B data for the ST pixels
        in the along-scan direction are corrected for ST contamination 
        according to the method in Barnes et al. (1995).  The routine 
        processes one scan line at a time using a rolling window for along-
        track processing of three scan lines for GAC data.  If along-track
        processing is indicated (nscan < 3), each scan line is processed
        individually without the rolling window. 

    Parameters: (in calling order)
        Type      Name         I/O      Description
        ----      ----         ---      -----------
        int32   *initial        I       Flag, should be set to "true" (a value
                                        of 1) if it is the fst call for a scene
        float32  Ltyp_frac      I       fraction of Ltypical for band 8
        float32  Styp_frac      I       fraction value that will be applied
                                        to knee value to calculate straylight
                                        threshold
        int32    nscans         I       the number of scan lines in the scene
                                        (to be processed).  Must be set to <5
                                        to process only in the along scan
        int32    nsamples       I       number of pixel data values in 
                                        l1b_data and sl_flag
        int32    scan_no        I       scan line number of l1b_data input
        int16	 gn		I 	Band 8 gain of scan line scan_no 
        float32	*rads		I	Radiances of calibration knees[8,4,5]  	
        float32 *l1b_data      I/O      array of nsamples long containing 
                                        l1b_data of scan line scan_no
        int32   *sl_scan        O       scan line number for which the 
                                        returned l1b_data apply
        int32   *sl_flag        O       array of nsamples long containing
                                        codes identifying BT pixels and the 
                                        proximity of others to BTs.
    Notes:

    Modification history:
          Programmer     Organization      Date      Description of change
        --------------   ------------    --------    ---------------------
        W. Robinson      SAIC            01/28/02    enable correction for
                                                     2nd pixel from BT and
                                                     fix problem in flagging
                                                     and correcting
        Joel Gales       Futuretech      10/16/00    1) Set AT in flag_buf
                                                     when initially setting
                                                     BT.
                                                     2) Set left edge found
                                                     if first pixel in scan
                                                     is above knee_val
        Lakshmi Kumar    Hughes STX      02/12/97    Added Styp_frac to input
                                                     parameter list
                                                     Removed PROTOTYPE defns.
                                                     corrected algorithm that
                                                     finds edges
        Lakshmi Kumar	 Hughes STX	 01/26/96    
        Lakshmi Kumar    Hughes STX      10/12/95    Original development
------------------------------------------------------------------------------*/

#include "st_lt.h"
#include <mfhdf.h>

int32_t stray_light_gac(int32_t *initial, float Ltyp_frac, float Styp_frac,
        int32_t nscans, int32_t nsamples, int32_t scan_no, int16_t gn,
        float *rads, float *l1b_data, int32_t *sl_scan, int32_t *sl_flag) {

    /*-----------------------------------------------------------------------------
     **  i		-- dummy index
     **  n		-- pixel pointer on the scan line being processed
     **  p		-- dummy index
     **  band	-- dummy index
     **  status	-- returns as 0 (NOTDONE) to indicate that the routine did not
                            return any useful information for the last call,
                            but was filling its buffers with scan line data; 
                            and 1 (DONE) to indicate that the l1b_data and sl_flag
                            arrays for scan line sl_scan have been set
     **  edge 	-- flag that indicates whether an edge was found or not
     **  scans_done  -- number of scans processed (if nscans < 5, scans_done = 2
     **  scan	-- pointer to the scan line being processed in the buffers--
                            the central scan line in a rolling window of five
                            consecutive scan lines
     **  scan_p1	-- pointer to the scan line following scan in the buffers
                            (scan plus 1)
     **  prev_scans	-- the number of scans preceding the one represented by scan
     **  next_scans  -- the number of scans following the one represented by scan
     **  line_no	-- a 5-integer array of the scan line numbers of the scan 
                            lines in the buffers
     **  flag_buf	-- a 2-dimensional buffer (5 x nsamples) for marking BT and
                            SL codes for five consecutive scan lines
     **  scan_m1	-- pointer to the scan line preceding scan lin the buffers
                            (scan minus 1)
     **  Ltypical_8	-- Ltypical value for band 8; = 1.09
     **  Ltyp_thresh -- value used to test if the band-8 radiance difference from
                            one pixel to the next indicates the edge
     **  l1b_buf	-- a 3-dimensional buffer (8 x 5 x nsamples) for storing five
                            consecutive scan lines of Level-1B data for 8 bands
     **  GAC_RRANGE 	-- Maximum number of pixels to the right of a BT that are
                            affected by ST from that BT; = 3
     **  GAC_LRANGE	-- Maximum number of pixels to the left (along scan) of a 
                            BT that are affected by ST from that BT; = 3
     **  K		-- BT response constants for along-scan stray light correction
                            array of two dimensions: 8 bands by GAC_LRANGE+GAC_RRANGE+1;
                            values are from Table 8 of Barnes et al. (1995),
                            except for value of central pixel which equals 1 minus
                            table value
    ----------------------------------------------------------------------------*/

    int32 i, n, p, farthest, status, pix;
    float32 delta, Ltyp_thresh, Ltypical_8 = 1.09;
    float32 knee_val = 0, stray_thresh = 0;
    int32 edge = 0, right_edge = -1, left_edge = -1;
    static int32 scans_done, scan, scan_p1, prev_scans, next_scans;
    static int32 scan_m1, rt_edge_found = FALSE;
    static int32 line_no[GACLINES], flag_buf[GACLINES][GACSAMPS];
    static float32 l1b_buf[MAXBANDS][GACLINES][GACSAMPS];
    float32 Ctyp_frac = 1.25, Ctyp_thresh;

    /*  If first call, initialize the counters and pointers and fill buffer
     **  with first scan line.  If along-track processing is requested (nscans > 2)
     **  then return for more scans to fill the buffers; otherwise, set 
     **  scans_done = 2. 
     */

    if (*initial == TRUE) {
        scans_done = 0;
        scan = 1;
        scan_p1 = 0;
        prev_scans = 0;
        next_scans = 0;
        line_no[scan - 1] = scan_no;
        for (i = 0; i < MAXBANDS; i++)
            memcpy(l1b_buf[i][scan - 1], &l1b_data[i * nsamples],
                (sizeof (float)*nsamples));
        for (i = 0; i < nsamples; i++)
            flag_buf[scan - 1][i] = BLANK;
        status = NOTDONE;
        if (nscans > GACLINES - 1) {
            *initial = FALSE;
            return status;
        }
    }

    if (*initial == FALSE) {
        status = NOTDONE;
        /* If only 1 scan has been obtained, get the 2nd and return for more */
        if (scan_p1 == 0) {
            scan_p1 = 2;
            next_scans = 1;
            line_no[scan_p1 - 1] = scan_no;
            for (i = 0; i < MAXBANDS; i++)
                memcpy((float *) l1b_buf[i][scan_p1 - 1],
                    &l1b_data[i * nsamples], sizeof (float)*nsamples);
            for (i = 0; i < nsamples; i++)
                flag_buf[scan_p1 - 1][i] = BLANK;
        }

        /* If only 2 scans have been obtained, rotate the array pointers
         ** (since the 1st scan line pointed to by scan has been processed),
         ** get the 3rd scan line, and process scan = 2.
         */
        if (scans_done == 1) {
            scan_m1 = 1;
            scan = 2;
            scan_p1 = 3;
            prev_scans = 1;
            next_scans = 1;
            line_no[scan_p1 - 1] = scan_no;
            for (i = 0; i < MAXBANDS; i++)
                memcpy((float *) l1b_buf[i][scan_p1 - 1],
                    &l1b_data[i * nsamples], sizeof (float)*nsamples);
            for (i = 0; i < nsamples; i++)
                flag_buf[scan_p1 - 1][i] = BLANK;
        }

        /*  If the buffers are full (as happens when scans_done > 1), rotate
         **  the array pointers */
        if (scans_done > 1) {
            scan_m1 = (scan_m1) % 3 + 1;
            scan = (scan) % 3 + 1;
            scan_p1 = (scan_p1) % 3 + 1;

            /*  If not all scan lines have been obtained, get the next one into
             **  scan_p1 buffers */
            if (scans_done < nscans - 1) {
                line_no[scan_p1 - 1] = scan_no;
                for (i = 0; i < MAXBANDS; i++)
                    memcpy((float *) l1b_buf[i][scan_p1 - 1],
                        &l1b_data[i * nsamples], sizeof (float)*nsamples);
                for (i = 0; i < nsamples; i++)
                    flag_buf[scan_p1 - 1][i] = BLANK;
            }
                /*  Otherwise, finish ones that haven't been processed.  If all have
                 **  been processed (next_scans < 0), flush out those remaining */
            else {
                next_scans = next_scans - 1;
                if (next_scans < 0) {
                    *sl_scan = line_no[scan_m1 - 1];
                    for (i = 0; i < MAXBANDS; i++)
                        memcpy((float *) &l1b_data[i * nsamples],
                            l1b_buf[i][scan_m1 - 1], sizeof (float)*nsamples);
                    for (i = 0; i < nsamples; i++)
                        sl_flag[i] = flag_buf[scan_m1 - 1][i];
                    status = DONE;
                    return status;
                } /* end if */
            } /* end if-else */
        } /* end if */
    } /* end if */

    if (*initial == TRUE) {
        scans_done = 2;
        scan_m1 = scan;
    }


    /* setup clear water threshold value Ctyp_thresh */
    Ctyp_thresh = Ctyp_frac * Ltypical_8;

    /*  if input radiance is greater than its knee, mark that pix as BT */
    knee_val = rads[7 * GAINS * KNEES + gn * KNEES + 1];
    stray_thresh = Styp_frac * knee_val;

    for (pix = 0; pix < nsamples; pix++) {
        if (l1b_buf[7][scan - 1][pix] > knee_val) {
            flag_buf[scan - 1][pix] = BT;
            if ((prev_scans >= 1) && (flag_buf[scan_m1 - 1][pix]) != BT)
                flag_buf[scan_m1 - 1][pix] = AT;
            if ((next_scans >= 1) && (flag_buf[scan_p1 - 1][pix] != BT))
                flag_buf[scan_p1 - 1][pix] = AT;
        }
    } /* end for */



    /*  Start loop to process the scan line pointed to by scan */
    rt_edge_found = FALSE;
    n = 0;
    while (n < nsamples - 1) {
        edge = NO;

        if (n == 0 && flag_buf[scan - 1][n] == BT) {
            edge = LEFT;
            rt_edge_found = FALSE;
        }


        while (edge == NO && (n < nsamples - 1)) {

            if (flag_buf[scan - 1][n] != BT || flag_buf[scan - 1][n + 1] != BT) {
                delta = l1b_buf[7][scan - 1][n + 1] - l1b_buf[7][scan - 1][n];
                if (delta > 0)
                    Ltyp_thresh =
                        max((Ltyp_frac * (l1b_buf[7][scan - 1][n + 1] - Ltypical_8)),
                        Ltyp_frac * Ltypical_8);
                else
                    Ltyp_thresh =
                        max((Ltyp_frac * (l1b_buf[7][scan - 1][n] - Ltypical_8)),
                        Ltyp_frac * Ltypical_8);
                if (delta > Ltyp_thresh && l1b_buf[7][scan - 1][n + 1] > stray_thresh) {
                    edge = LEFT;
                    rt_edge_found = FALSE;
                }
                if (-1 * delta > Ltyp_thresh && l1b_buf[7][scan - 1][n] > stray_thresh) {
                    if (!(rt_edge_found)) {
                        edge = RIGHT;
                        rt_edge_found = TRUE;
                    }
                }

            }
            n++;
        }

        if (edge != NO) {
            n--;

            /*  If right edge is found first, set the scan-line start as the left
             **  edge, and correct the ST pixels which are 1 and two pixels away 
             **  from the right edge. 
             **  CF_right = along scan stray-light GAC correctionn factors for
             **  pixels to the right of a BT; array of two dimensions: 8 bands by
             **  GAC_RRANGE; values are from Table 15 of Barnes et al. (1995).
             **  The above table "CF_right" is defined in cf.h file */

            if (edge == RIGHT) {
                right_edge = n;
                for (i = n; i >= 0; i--) if (l1b_buf[7][scan - 1][i] < Ctyp_thresh) break;
                if (i < 0) i = 0;
                left_edge = i;
                p = right_edge + 1;
                farthest = min((right_edge + GAC_RRANGE), (nsamples - 1));
                while (p <= farthest) {
                    if (flag_buf[scan - 1][p] == BT)
                        p = farthest + 1;
                    else {
                        if (flag_buf[scan - 1][p] == BLANK) {
                            flag_buf[scan - 1][p] = p - right_edge; /*flag*/
                            if (p == (right_edge + GAC_RRANGE)) {
                                for (i = 0; i < MAXBANDS; i++)
                                    l1b_buf[i][scan - 1][p] += /*correction*/
                                        CF_right[i][2] * l1b_buf[i][scan - 1][right_edge];
                            } else if (p == (right_edge + GAC_RRANGE - 1)) {
                                for (i = 0; i < MAXBANDS; i++)
                                    l1b_buf[i][scan - 1][p] += /*correction*/
                                        CF_right[i][1] * l1b_buf[i][scan - 1][right_edge];
                            }
                        }
                    } /* end else */
                    p += 1;
                } /* end while */
            }/* end if */
            else {
                /* if a left edge is found, correct the 2nd and third pixels 
                 ** prior to the left edge. 
                 ** CF_left = along-scan stray-light GAC correction factors for 
                 ** pixels to the left of a BT; array of two dimensions: 8 bands by 
                 ** Lrange; values are from Table 15 of Barnes et al. (1995).
                 ** The table "CF_left" is defined in cf.h file */

                left_edge = n + 1;
                p = n;
                farthest = max(left_edge - GAC_LRANGE, 0);
                while (p >= farthest) {
                    if (flag_buf[scan - 1][p] == BT)
                        p = farthest - 1;
                    else {
                        if (flag_buf[scan - 1][p] == BLANK) {
                            flag_buf[scan - 1][p] = left_edge - p; /* flag */
                            if (p == (left_edge - GAC_LRANGE)) {
                                for (i = 0; i < MAXBANDS; i++)
                                    l1b_buf[i][scan - 1][p] += /* correction*/
                                        (CF_left[i][2] * l1b_buf[i][scan - 1][left_edge]);
                            } else if (p == (left_edge - GAC_LRANGE + 1)) {
                                for (i = 0; i < MAXBANDS; i++)
                                    l1b_buf[i][scan - 1][p] += /* correction*/
                                        (CF_left[i][1] * l1b_buf[i][scan - 1][left_edge]);
                            }
                        } else if (flag_buf[scan - 1][p] > 0) {
                            flag_buf[scan - 1][p] += 1000 * (left_edge - p); /* flag */
                            if (p == (left_edge - GAC_LRANGE)) {
                                for (i = 0; i < MAXBANDS; i++)
                                    l1b_buf[i][scan - 1][p] += /* correction*/
                                        (CF_left[i][2] * l1b_buf[i][scan - 1][left_edge]);
                            } else if (p == (left_edge - GAC_LRANGE + 1)) {
                                for (i = 0; i < MAXBANDS; i++)
                                    l1b_buf[i][scan - 1][p] += /* correction*/
                                        (CF_left[i][1] * l1b_buf[i][scan - 1][left_edge]);
                            }
                        }
                    } /* end else */
                    p = p - 1;
                } /* end while */

                /* search for the associated right edge */
                right_edge = nsamples - 1;
                n = n + 1;
                while ((n < nsamples - 1) && (edge == LEFT)) {
                    if (flag_buf[scan - 1][n] != BT || flag_buf[scan - 1][n + 1] != BT) {
                        if (l1b_buf[7][scan - 1][n + 1] < Ctyp_thresh) goto R_found;
                        delta = l1b_buf[7][scan - 1][n + 1] - l1b_buf[7][scan - 1][n];
                        if (delta > 0)
                            Ltyp_thresh =
                                max((Ltyp_frac * (l1b_buf[7][scan - 1][n + 1] - Ltypical_8)),
                                Ltyp_frac * Ltypical_8);
                        else
                            Ltyp_thresh =
                                max((Ltyp_frac * (l1b_buf[7][scan - 1][n] - Ltypical_8)),
                                Ltyp_frac * Ltypical_8);
                        if (-1 * delta > Ltyp_thresh &&
                                l1b_buf[7][scan - 1][n] > stray_thresh) {
                        R_found:
                            edge = RIGHT;
                            rt_edge_found = TRUE;
                        }
                    } /*if BT */
                    n++;
                } /*while ((n < nsamples - 1) && (edge == LEFT) */

                if (edge == LEFT) /* Cound not find right_edge */
                    rt_edge_found = TRUE; /* setting the last pix as rt edge*/
                /* When the right edge is found, apply the correction to the
                 ** 2nd, third pixels down from the right edge */
                /* WDR fix flag_buf check to account for corrections from 
                   previous right edges */
                if (edge == RIGHT) {
                    right_edge = n - 1;
                    p = right_edge + 1;
                    farthest = min((right_edge + GAC_RRANGE), (nsamples - 1));
                    while (p <= farthest) {
                        if (flag_buf[scan - 1][p] == BT)
                            p = farthest + 1;
                        else {
                            if ((flag_buf[scan - 1][p] == BLANK) ||
                                    (flag_buf[scan - 1][p] > 0)) {
                                flag_buf[scan - 1][p] = p - right_edge; /* flag */
                                if (p == (right_edge + GAC_RRANGE)) {
                                    for (i = 0; i < MAXBANDS; i++) /* correction */
                                        l1b_buf[i][scan - 1][p] +=
                                            (CF_right[i][2] * l1b_buf[i][scan - 1][right_edge]);
                                } else if (p == (right_edge + GAC_RRANGE - 1)) {
                                    for (i = 0; i < MAXBANDS; i++) /* correction */
                                        l1b_buf[i][scan - 1][p] +=
                                            (CF_right[i][1] * l1b_buf[i][scan - 1][right_edge]);
                                }
                            }
                        } /* end else */
                        p += 1;
                    } /* end while */
                } /* end if */
            } /* end else */
            /* Label all pixels between both edges as BT pixels.  Some of the
             ** pixels might already be BT, will get set again to avoid another 
             ** if check   */

            for (p = left_edge; p <= right_edge; p++) {
                flag_buf[scan - 1][p] = BT;
                if ((prev_scans >= 1) && (flag_buf[scan_m1 - 1][p]) != BT)
                    flag_buf[scan_m1 - 1][p] = AT;
                if ((next_scans >= 1) && (flag_buf[scan_p1 - 1][p] != BT))
                    flag_buf[scan_p1 - 1][p] = AT;
            } /* end for */
        } /* end if */

        if (rt_edge_found && n < nsamples - 1)
            n = right_edge + 1;

    } /* end while */


    /* If three scan lines have been obtained (and one has been processed), 
     ** return the earliest scan line (scan_m1) in the window.  */

    scans_done = scans_done + 1;
    if (scans_done > 1) {
        *sl_scan = line_no[scan_m1 - 1];

        for (i = 0; i < MAXBANDS; i++)
            memcpy(&l1b_data[i * nsamples],
                l1b_buf[i][scan_m1 - 1], (sizeof (float)*nsamples));
        for (i = 0; i < nsamples; i++)
            sl_flag[i] = flag_buf[scan_m1 - 1][i];

        status = DONE;
    } /* end if */
    return status;
}

