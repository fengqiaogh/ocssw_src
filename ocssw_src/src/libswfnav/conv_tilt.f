      subroutine conv_tilt(ang_fm,ang_mm,tilt)

c $Header$
c $Log$
c
c
c  conv_tilt(ang_fm,ang_mm,tilt)
c
c  Purpose: Convert tilt motor telemetry angles to instrument tilt
c
c  Calling Arguments:
c
c  Name         Type    I/O     Description
c  --------     ----    ---     -----------
c  ang_fm       R*4      I      Fixed motor angle
c  ang_mm       R*4      I      Moving motor angle
c  tilt         R*4      I      Instrument Tilt angle
c
c  By: Frederick S. Patt, GSC, 25 April 1996
c
c  Notes:       
c
c  Algorithm described in "SeaWiFS Tilt Telemetry Analysis", F. Patt, 
c  informal memorandum, July 21, 1994
c
c  Modification History:
c
c  Adjusted tilt telemetry angle calibrations to force agreement between
c  tilt angles computed from telemetry and static values from SBRS.
c  F. S. Patt, GSC, November 10, 1997
        
        implicit none
        real*4 ang_fm,ang_mm,tilt
        real*4 alen,blen,rlen,llen
        real*4 fscal,foff,fref,mscal,moff,mref,toff
        real*4 loc_f,loc_m,apri,bpri,tpri,dtf,dtm
        real*8 pi,radeg,re,rem,f,omf2,omegae
        common /gconst/pi,radeg,re,rem,f,omf2,omegae

c  The values below are given in the referenced memorandum
c    Dimensions for tilt components 
        data alen/7.088/, blen/6.000/, rlen/1.025/, llen/4.070/

c    Reference angles for telemetry measurements and tilt angle
        data mref/235.0/, fref/80.0/, toff/-35.068/

c  These are the corrections for the tilt motor angles from the 
c    nominal conversions provided by SBRC
        data fscal/1.0/, foff/-1.5/, mscal/.964/, moff/0.0/


c  Apply corrections to motor angles
        loc_f = fscal*ang_fm + foff
        loc_m = mscal*ang_mm + moff

c  Convert tilt as described in memo

c  Compute pivot to link distances (equations 1 and 2)
        apri = sqrt(alen**2 + rlen**2 -
     *    2.0*alen*rlen*cos((loc_f+fref)/radeg))
        bpri = sqrt(blen**2 + rlen**2 -
     *    2.0*blen*rlen*cos((loc_m+mref)/radeg))

c  Compute angle theta-prime (equation 3)
        tpri = acos((apri**2 + bpri**2 - llen**2)/(2.0*apri*bpri))

c  Compute corrections to theta-prime for fixed and moving motors
c    (equations 4 and 5)
        dtf = asin(rlen*sin((loc_f+fref)/radeg)/apri)
        dtm = asin(rlen*sin((loc_m+mref)/radeg)/bpri)

c  Compute tilt angle (equations 6 and 7)
        tilt = (tpri - dtf - dtm)*radeg + toff

        return
        end
