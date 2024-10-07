/*
 *******************************************************************************

 !C-INC

 !Description:
     calibration structure include file

 !Input Parameters: N/A
 !Output Parameters: N/A

 !Revision History:

        $Id: runcal.h,v 1.2 2012/04/26 18:55:32 sue Exp $

        $Log: runcal.h,v $
        Revision 1.2  2012/04/26 18:55:32  sue
            Changes made by seadas group for seadas6.3, newer fortran, and/or 64 bit mode.
	
        Revision 1.3  2002/08/26 11:46:26  sue
        Update copyright notices for year 2002.
	
        Revision 1.2  1997/11/14 14:12:59  sue
        Fix prologs.

        Revision 1.1.1.1  1996/09/17 16:40:07  angel

 * Revision 1.2  1992/07/13  13:24:18  angel
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

#ifndef _RUNCAL_H_
#define _RUNCAL_H_

#include "stdint.h"

#define MAX_RUNCAL 75

typedef struct {
    int32_t numcal;
    int32_t intrvl;

    struct _runCalEnt {
        float cenlin;
        float prtemp;
        float slope[3];
        float intcp[3];
    } runcal[MAX_RUNCAL];
} RUN_CAL;

#endif /* _RUNCAL_H_ */
