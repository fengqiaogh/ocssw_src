/*
        subroutine jddate(jd,i,j,k)
C
C       This routine computes the calendar date corresponding to
C       a given Julian day.  This code was brazenly copied from
C       a routine written by Myron Shear for CSC on Julian Day 1.
C       
C       ARGUMENT        TYPE    I/O     DESCRIPTION     
C       __________________________________________________________
C        JD             I*4      I      Julian Day (reference Jan 1, 4713 BC)
C        I              I*4      O      Year 
C        J              I*4      O      Month   
C        K              I*4      0      Day of Month
C
        l = jd + 68569
        n = 4*l/146097
        l = l - (146097*n + 3)/4
        i = 4000*(l+1)/1461001
        l = l - 1461*i/4 + 31
        j = 80*l/2447
        k = l - 2447*j/80
        l = j/11
        j = j + 2 - 12*l
        i = 100*(n-49) + i + l
        return
        end
 */

/*
        subroutine jdate(jd,i,k)
C
C       This routine computes the year and day-of-year corresponding 
C       to a given Julian day.  This algorithm is designed for the 
C       period 1900 - 2100. 
C       
C       ARGUMENT        TYPE    I/O     DESCRIPTION     
C       __________________________________________________________
C        JD             I*4      I      Julian Day (reference Jan 1, 4713 BC)
C        I              I*4      O      Year 
C        K              I*4      0      Day of Year
C
c       Program written by:     Frederick S. Patt
c                               General Sciences Corporation
c                               May 12, 1993

c       Compute days since January 0, 1900
        l = jd - 2415020

c       Compute years since 1900
        i = 4*l/1461

c       Compute day-of-year
        k = l - 1461*(i-1)/4 - 365 

c       Add first two digits of year
        i = i + 1900
        return
        end
 */

/*
      function jd(i,j,k)
c
c
c    This function converts a calendar date to the corresponding Julian
c    day starting at noon on the calendar date.  The algorithm used is
c    from Van Flandern and Pulkkinen, Ap. J. Supplement Series 41, 
c    November 1979, p. 400.
c
c
c       Arguments
c     
c       Name    Type    I/O     Description
c       ----    ----    ---     -----------
c       i       I*4      I      Year - e.g. 1970
c       j       I*4      I      Month - (1-12)
c       k       I*4      I      Day  - (1-31)
c       jd      I*4      O      Julian day
c
c     external references
c     -------------------
c      none
c
c
c     Written by Frederick S. Patt, GSC, November 4, 1992
c
c
        jd = 367*i - 7*(i+(j+9)/12)/4 + 275*j/9 + k + 1721014

c  This additional calculation is needed only for dates outside of the 
c   period March 1, 1900 to February 28, 2100
c       jd = jd + 15 - 3*((i+(j-9)/7)/100+1)/4
        return
        end
 */

int jd_c(int y, int m, int d) {
    return ( 367 * y -
            7 * (y + (m + 9) / 12) / 4 +
            275 * m / 9 +
            d + 1721014);
}
