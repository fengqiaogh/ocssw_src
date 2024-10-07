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
