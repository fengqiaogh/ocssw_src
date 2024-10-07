/***************************************************************

 !C

    function ftrim(str, len)
    !Description:	
        Remove trailing blanks from a FORTRAN string.
    !Input Parameters:
        str(char *) - string
        len(int *) - length of string
    !Output Parameters:
        str(char *) - updated string
    !Returns:
        ftrim(int) - length of updated string

 !Revision History:
    $Log: ftrim.c,v $
    Revision 1.1  2010/08/07 18:44:32  sue
    seadas 6.1 l2gen with modis dust and avhrr.

    Revision 1.1  2008/08/15 20:30:25  sue
    NODC versions of the pathfinder processing programs.

    Revision 1.5  2002/08/23 20:50:22  sue
    Update copyright notices to year 2002.

    Revision 1.4  1997/11/20 19:37:40  jim
    Fix prologs.

    Revision 1.3  1997/09/22 16:39:14  jim
    Change included function test program to not use a prohibited function.

    Revision 1.2  1996/09/16 22:14:14  sue
    Fix prologues and get rid of compiler warnings.

    Revision 1.1  1996/09/16 15:07:43  sue
    Modis versions of rtlib routines.

 * Revision 1.3  1994/07/25  18:16:19  angel
 * Include "rtlib.h".
 *
 * Revision 1.2  1993/11/10  17:45:36  jim
 * Only add _ on certain unix systems, never on vms.
 *
 * Revision 1.1  1992/02/26  19:57:51  angel
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

 !References and Credits:
        Written by:
        University of Miami
        Rosensteil School for Marine and Atmospheric Science
        Division of Meteorology and Oceanography Oceanography
        4600 Rickenbacker Cswy.
        Miami, Fl
        33149

        contact: SWalsh@rsmas.miami.edu

 !Design Notes:

 !END****************************************************************
 */
// #include "rtlib.h"

/* To use this test, add a line at the end of main to print out buf,
 * or change the sprintf to printf
 */
#ifdef TEST

main() {
    char s[80];
    int i;
    int l = 80;
    char buf[80];

    for (i = 0; i < 80; i++)
        s[i] = ' ';
    /*
    s[0] = 'a';
    s[1] = 'b';
     */
    l = ftrim(s, &l);
    sprintf(buf, "%d %s\n", l, s);
}
#endif

int ftrim_(str, len)
char *str;
int *len;
{
    char *s;
    int l = *len;
    int i;

    i = l - 1;
    for (s = str + i; i >= 0; i--, s--) {
        if (*s != ' ') {
            if (i != l - 1)
                *(s + 1) = 0;
            return i + 1;
        }
    }
    str[0] = 0;
    return 0;
}
