	pro gha2000,iyr,day,gha,xlm,omega,eps

;  This subroutine computes the Greenwich hour angle in degrees for the
;  input time.  It uses the model referenced in The Astronomical Almanac
;  for 1984, Section S (Supplement) and documented for the SeaWiFS
;  Project in "Constants and Parameters for SeaWiFS Mission Operations",
;  in TBD.  It includes the correction to mean sideral time for nutation
;  as well as precession.
;
;	implicit real*8 (a-h,o-z)
	imon=1
;	common nutcm,dpsi,eps,nutime
;	radeg = 180.d0/3.14159265359d0

;  Compute days since J2000
	iday = long(day)
	fday = day - iday
	jd,iyr,1,iday,jday
	t = jday - 2451545.5d0 + fday
	
;  Compute Greenwich Mean Sidereal Time	(degrees)
	gmst = 100.4606184d0 + 0.9856473663d0*t + 2.908d-13*t*t

;  Check if need to compute nutation correction for this day
;	nt = t
;	if (nt ne nutime) then begin
;	  nutime = nt
	  ephparms,t,xls,gs,xlm,omega
	  nutate,t,xls,gs,xlm,omega,dpsi,eps
;	endif

;  Include apparent time correction and time-of-day
	gha = gmst + dpsi*cos(eps/!radeg) + fday*360.d0
	gha = gha mod 360.d0
	neg = where(gha lt 0.d0)
        if (neg(0) gt -1) then gha(neg) = gha(neg) + 360.d0
;
	return
	end

	pro ephparms,t,xls,gs,xlm,omega

;  This subroutine computes ephemeris parameters used by other Mission
;  Operations routines:  the solar mean longitude and mean anomaly, and
;  the lunar mean longitude and mean ascending node.  It uses the model
;  referenced in The Astronomical Almanac for 1984, Section S 
;  (Supplement) and documented for the SeaWiFS Project in "Constants
;  and Parameters for SeaWiFS Mission Operations", in TBD.  These
;  parameters are used to compute the solar longitude and the nutation
;  in longitude and obliquity.
;
;	implicit real*8 (a-h,o-z)
;	radeg = 180.d0/3.14159265359d0

;  Sun Mean Longitude 		
	xls = 280.46592d0 + 0.9856473516d0*t
	xls = xls mod 360.d0
 
;  Sun Mean Anomaly		
	gs = 357.52772d0 + 0.9856002831d0*t 
	gs = gs mod 360.d0

;  Moon Mean Longitude		
	xlm = 218.31643d0 + 13.17639648d0*t 
	xlm = xlm mod 360.d0

;  Ascending Node of Moon's Mean Orbit 	
	omega = 125.04452d0 - 0.0529537648d0*t 
	omega = omega mod 360.d0
	
	return
	end

	pro nutate,t,xls,gs,xlm,omega,dpsi,eps,epsm

;  This subroutine computes the nutation in longitude and the obliquity
;  of the ecliptic corrected for nutation.  It uses the model referenced
;  in The Astronomical Almanac for 1984, Section S (Supplement) and 
;  documented for the SeaWiFS Project in "Constants and Parameters for 
;  SeaWiFS Mission Operations", in TBD.  These parameters are used to 
;  compute the apparent time correction to the Greenwich Hour Angle and 
;  for the calculation of the geocentric Sun vector.  The input 
;  ephemeris parameters are computed using subroutine ephparms.  Terms 
;  are included to 0.1 arcsecond.
;
;	implicit real*8 (a-h,o-z)
;	radeg = 180.d0/3.14159265359d0

	
;  Nutation in Longitude
	dpsi = - 17.1996*sin(omega/!radeg) 	$
     	 	+ 0.2062*sin(2.*omega/!radeg)	$
     	     	- 1.3187*sin(2.*xls/!radeg) 	$
     		+ 0.1426*sin(gs/!radeg) 	$
     		- 0.2274*sin(2.*xlm/!radeg)
	
;  Mean Obliquity of the Ecliptic	
	epsm = 23.439291d0 - 3.560d-7*t 

;  Nutation in Obliquity 
	deps = 9.2025*cos(omega/!radeg) + 0.5736*cos(2.*xls/!radeg)

;  True Obliquity of the Ecliptic 
	eps = epsm + deps/3600.d0

	dpsi = dpsi/3600.d0
        return
	end

