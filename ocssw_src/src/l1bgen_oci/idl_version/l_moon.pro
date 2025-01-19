	pro l_moon,iyr,iday,sec,xmnr,rm
;  Computes unit Moon vector in geocentric rotating coodinates, using 
;  subprograms to compute inertial Moon vector and Greenwich hour angle
;
;
;	Arguments:
;
;	Name	Type	I/O	Description
;	--------------------------------------------------------
;	IYR	I*4	 I	Year, four digits (i.e, 1993)
;	IDAY	I*4	 I	Day of year (1-366)
;	SEC	R*8	 I	Seconds of day 
;	XMNR(3)	R*4	 O	Unit Moon vector in geocentric rotating 
;				 coordinates
;	RM	R*4	 O	Earth-to-Moon distance (km)
;
;	Subprograms referenced:
;
;	MOON2000	Computes inertial Moon vector
;	GHA2000		Computes Greenwich sidereal angle
;
;	Coded by:  Frederick S. Patt, GSC, October 4, 1992
;
;	Modification History:
;
;       IDL version created, F.S. Patt, SAIC, September 5, 2003.

;       implicit real*8 (a-h,o-z)
;	real*4 xmnr(3),xm(3),rm
;       common /gconst/pi,radeg,re,rem,f,omf2,omegae
	
;  Get unit Moon vector in geocentric inertial coordinates
        moon2000,iyr,iday,sec,xm,rm,beta

;  Get Greenwich mean sideral angle
	day = iday + sec/864.d2 
	gha2000,iyr,day,gha
	ghar = gha/!radeg

;  Transform Moon vector into geocentric rotating frame
        xmnr = xm
	xmnr(0,*) = xm(0,*)*cos(ghar) + xm(1,*)*sin(ghar)
	xmnr(1,*) = xm(1,*)*cos(ghar) - xm(0,*)*sin(ghar)
;	xmnr(2,*) = xm(2,*)

	return
	end

	pro moon2000,iyr,iday,sec,xmn,rm,xlma,beta
;
;  This subroutine computes the Moon vector in geocentric inertial 
;  (equatorial) coodinates.  It uses the model described in "Low-precision 
;  Formula for Planetary Positions", T.C. Van Flandern and K.F. Pulkkinen,
;  Ap. J. Supplement Series, 41:391-411, 1979.  The accuracy of the Moon 
;  vector is approximately 1 arcminute.
;
;	Arguments:
;
;	Name	Type	I/O	Description
;	--------------------------------------------------------
;	IYR	I*4	 I	Year, four digits (i.e, 1993)
;	IDAY	I*4	 I	Day of year (1-366)
;	SEC	R*8	 I	Seconds of day 
;	XMN(3)	R*4	 O	Unit Moon vector in geocentric inertial 
;				 coordinates of date
;	RM	R*4	 O	Magnitude of the Moon vector (km)
;
;	Subprograms referenced:
;
;	JD		Computes Julian day from calendar date
;	EPHPARMS	Computes mean solar longitude and anomaly and
;			 mean lunar lontitude and ascending node
;	NUTATE		Compute nutation corrections to lontitude and 
;			 obliquity
;
;	Coded by:  Frederick S. Patt, GSC, September 29, 1993


;	implicit real*8 (a-h,o-z)
;	real*4 xmn(3),xme(3),rm
;	integer*4 iyr,iday,imon/1/
        imon = 1
;	common /gconst/pi,radeg,re,rem,f,omf2,omegae
        re = 6378.137d0
;	common /nutcm/dpsi,eps,nutime

;   Compute days since Jan 2, 1980 and add seconds as fraction of day
	jd,iyr,imon,iday,jday
	t = jday - 2451545.0d0 + (sec-43200.d0)/86400.d0

;  Compute solar ephemeris parameters
	ephparms,t,xls,gs,xlm,omega

;  Check if need to compute nutation corrections for this day
;	nt = t
;	if (nt.ne.nutime) then
;	  nutime = nt
	  nutate,t,xls,gs,xlm,omega,dpsi,eps
;	end if

;  Compute Moon mean anomaly, argument of latitude and mean elongation
	gm = 134.96292d0 + 13.06499295d0*t 	
	gm = gm mod 360.d0

	fm = xlm - omega
	fm = fm mod 360.d0

	dm = xlm - xls
	dm = dm mod 360.d0


;  Compute lunar distance (Re)

	rm = 60.36298d0 - 3.27746d0*cos(gm/!radeg) $
     		- 0.57994d0*cos((gm - 2.d0*dm)/!radeg) $
     		- 0.46357d0*cos(2.d0*dm/!radeg) $
     		- 0.08904d0*cos(2.d0*gm/!radeg) $
     		+ 0.03865d0*cos((2.d0*gm - 2.d0*dm)/!radeg) $
     		- 0.03237d0*cos((2.d0*dm - gs)/!radeg) $
     		- 0.02688d0*cos((gm + 2.d0*dm)/!radeg) $
     		- 0.02358d0*cos((gm - 2.d0*dm + gs)/!radeg) $
     		- 0.02030d0*cos((gm - gs)/!radeg) $
     		+ 0.01719d0*cos(dm/!radeg) $
     		+ 0.01671d0*cos((gm + gs)/!radeg) $
     		+ 0.01247d0*cos((gm - 2.d0*fm)/!radeg)


;  Compute lunar longitude
	dlm = 22640.d0*sin(gm/!radeg) $
     		- 4586.d0*sin((gm - 2.d0*dm)/!radeg) $
     		+ 2370.d0*sin(2.d0*dm/!radeg) $
     		+  769.d0*sin(2.d0*gm/!radeg) $
     		-  668.d0*sin(gs/!radeg) $
     		-  412.d0*sin(2.d0*fm/!radeg) $
     		-  212.d0*sin((2.d0*gm - 2.d0*dm)/!radeg) $
     		-  206.d0*sin((gm - 2.d0*dm + gs)/!radeg) $
     		+  192.d0*sin((gm + 2.d0*dm)/!radeg) $
     		+  165.d0*sin((2.d0*dm - gs)/!radeg) $
     		+  148.d0*sin((gm - gs)/!radeg) $
     		-  125.d0*sin(dm/!radeg) $
     		-  110.d0*sin((gm + gs)/!radeg) $
     		-   55.d0*sin((2.d0*fm - 2.d0*dm)/!radeg) $
     		-   45.d0*sin((gm + 2.d0*fm)/!radeg) $
     		+   40.d0*sin((gm - 2.d0*fm)/!radeg) $
     		-   38.d0*sin((gm - 4.d0*dm)/!radeg) $
     		+   36.d0*sin(3.d0*gm/!radeg) $
     		-   31.d0*sin((2.d0*gm - 4.d0*dm)/!radeg) $
     		+   28.d0*sin((gm - 2.d0*dm - gs)/!radeg) $
     		-   24.d0*sin((2.d0*dm + gs)/!radeg) $
     		+   19.d0*sin((gm - dm)/!radeg) $
     		+   18.d0*sin((dm + gs)/!radeg) $
     		+   15.d0*sin((gm + 2.d0*dm - gs)/!radeg) $
     		+   14.d0*sin((2.d0*gm + 2.d0*dm)/!radeg) $
     		+   14.d0*sin((4.d0*dm)/!radeg) $
     		-   13.d0*sin((3.d0*gm - 2.d0*dm)/!radeg)

	xlma = xlm + dlm/3600.d0 + dpsi
; 	if (xlma lt 0.d0) then xlma = xlma + 360.d0

;   Compute lunar ecliptic latitude 
	beta = 18461.d0*sin(fm/!radeg) $
     		+ 1010.d0*sin((gm + fm)/!radeg) $
     		+ 1000.d0*sin((gm - fm)/!radeg) $
     		-  624.d0*sin((fm - 2.d0*dm)/!radeg) $
     		-  199.d0*sin((gm - fm - 2.d0*dm)/!radeg) $
     		-  167.d0*sin((gm + fm - 2.d0*dm)/!radeg) $
     		+  117.d0*sin((fm + 2.d0*dm)/!radeg) $
     		+   62.d0*sin((2.d0*gm + fm)/!radeg) $
     		+   33.d0*sin((gm - fm + 2.d0*dm)/!radeg) $
     		+   32.d0*sin((2.d0*gm - fm)/!radeg) $
     		-   30.d0*sin((fm - 2.d0*dm + gs)/!radeg) $
     		-   16.d0*sin((2.d0*gm + fm - 2.d0*dm)/!radeg) $
     		+   15.d0*sin((gm + fm + 2.d0*dm)/!radeg) $
     		+   12.d0*sin((fm - 2.d0*dm - gs)/!radeg)

	beta = beta/3600.d0

;	type *, xlm, fm, dlm/3600., beta

;   Compute unit Moon vector in ecliptic plane
        xme = fltarr(3,n_elements(t))
	xme(0,*) = cos(xlma/!radeg)*cos(beta/!radeg)
	xme(1,*) = sin(xlma/!radeg)*cos(beta/!radeg)
	xme(2,*) = sin(beta/!radeg)

;   Rotate vector to equatorial plane
        xmn = xme
;	xmn(0) = xme(0)
	xmn(1,*) = xme(1,*)*cos(eps/!radeg) - xme(2,*)*sin(eps/!radeg)
	xmn(2,*) = xme(1,*)*sin(eps/!radeg) + xme(2,*)*cos(eps/!radeg)

	rm = rm*re
;	type *,' lon,lat,dist=',xlma,beta,rm

	return
	end

      pro jd,i,j,k,jday
;
;
;    This function converts a calendar date to the corresponding Julian
;    day starting at noon on the calendar date.  The algorithm used is
;    from Van Flandern and Pulkkinen, Ap. J. Supplement Series 41, 
;    November 1979, p. 400.
;
;
;	Arguments
;     
;     	Name    Type 	I/O 	Description
;     	----	---- 	--- 	-----------
;     	i	I*4  	 I 	Year - e.g. 1970
;     	j       I*4  	 I  	Month - (1-12)
;     	k       I*4  	 I  	Day  - (1-31)
;     	jday    I*4  	 O  	Julian day
;
;     external references
;     -------------------
;      none
;
;
;     Written by Frederick S. Patt, GSC, November 4, 1992
;
;
	i = long(i)
	j = long(j)
	k = long(k)
      	jday = 367*i - 7*(i+(j+9)/12)/4 + 275*j/9 + k + 1721014

;  This additional calculation is needed only for dates outside of the 
;   period March 1, 1900 to February 28, 2100
;     	jday = jday + 15 - 3*((i+(j-9)/7)/100+1)/4
      	return
      	end


