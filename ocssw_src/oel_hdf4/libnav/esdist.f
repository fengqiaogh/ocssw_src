
        real*8 function esdist(iyr,iday,msec)
c
c  This subroutine computes the earth-sun distance in AU. It uses the model 
c  referenced in The Astronomical Almanac for 1984, Section S (Supplement).
c
c       Arguments:
c
c       Name    Type    I/O     Description
c       --------------------------------------------------------
c       IYR     I*4      I      Year, four digits (i.e, 1993)
c       IDAY    I*4      I      Day of year (1-366)
c       MSEC    I*4      I      milliseconds of day 
c       RS      R*4      O      Magnitude of the Sun vector (AU)
c
c       Subprograms referenced:
c
c       JD              Computes Julian day from calendar date
c
c       Coded by:  Frederick S. Patt, GSC, November 2, 1992
c       Adapted from sun2000 to esdist by B. Franz, Oct 2003.
c
        implicit none
        integer*4 iyr, iday, msec, jd
        real*8    t, gs
c
        real*8    radeg
        integer*4 imon
c
        data radeg /57.29577951/
        data imon  /1/


c   Compute floating point days since Jan 1.5, 2000 
c   Note that the Julian day starts at noon on the specified date

        t = jd(iyr,imon,iday) 
     *    - 2451545.0d0 
     *    + (msec/1000.d0 - 43200.d0)/86400.d0

c  Compute mean anomaly

        gs = 357.52772d0 + 0.9856002831d0*t

c  Compute solar distance (AU)
        esdist = 1.00014D0 - 0.01671D0*cos(gs/radeg) 
     *         - 0.00014D0*cos(2.0D0*gs/radeg)

        return
        end
