/*
 *******************************************************************************

 !C-INC

 !Description:
     calibration information include file

 !Input Parameters: N/A
 !Output Parameters: N/A

 !Revision History:

        $Id: rawcal.h,v 1.1 2010/08/07 18:44:48 sue Exp $

        $Log: rawcal.h,v $
        Revision 1.1  2010/08/07 18:44:48  sue
        seadas 6.1 l2gen with modis dust and avhrr.
	
        Revision 1.3  2002/08/26 11:46:26  sue
        Update copyright notices for year 2002.
	
        Revision 1.2  1997/11/14 14:12:54  sue
        Fix prologs.

        Revision 1.1.1.1  1996/09/17 16:40:07  angel

 * Revision 1.3  1993/12/18  17:02:51  angel
 * ANSIfy the headers.
 *
 * Revision 1.2  1992/07/13  13:24:15  angel
 * Change from int32 to INT32, etc.  This is the start of support for the
 * Alpha chip.
 *
 * Revision 1.1  1992/03/12  21:16:46  angel
 * Initial revision
 *
 */

/*
  !Team-unique Header:

  Copyright 1988-2002 by Rosenstiel School of Marine and Atmospheric Science,
  University of Miami, Miami, Florida.

                        All Rights Reserved

  Permission to use, copy, modify, and distribute this software and its
  documentation for non-commercial purposes and without fee is hereby granted,
  provided that the above copyright notice appear in all copies and that both
  that copyright notice and this permission notice appear in supporting
  documentation, and that the names of University of Miami and/or RSMAS not be
  used in advertising or publicity pertaining to distribution of the software
  without specific, written prior permission.

  UNIVERSITY OF MIAMI DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
  INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS, IN NO EVENT
  SHALL UNIVERSITY OF MIAMI BE LIABLE FOR ANY SPECIAL, INDIRECT OR CONSEQUENTIAL
  DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS,
  WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING
  OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.

  !References and credits:
     Written by:
     University of Miami
     Rosenstiel School for Marine and Atmospheric Science
     Division of Meteorology and Physical Oceanography
     4600 Rickenbacker Cswy
     Miami,Fl
     Contact: SWalsh@rsmas.miami.edu

  !Design Notes

  !END*************************************************************
 */

#ifndef _RAWCAL_H_
#define _RAWCAL_H_

#include "stdint.h"

#define MAX_RUNCAL 75

typedef struct {
    int32_t numraw;
    int32_t intrvl;

    struct _rawCalEnt {
        float scan;
        float numsmp;
        double telem[1 * 5];
        double teletg[3 * 4 * 5];
        double bckscn[5][30];
        double space[5][50];
    } rawcal[MAX_RUNCAL];
} RAW_CAL;

/*
extern int RunCalXdr _ANSI_ARGS_((RUN_CAL *, int));
 */

#endif /* _RAWCAL_H_ */
