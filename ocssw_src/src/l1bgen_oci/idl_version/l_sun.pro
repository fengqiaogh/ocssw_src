	pro l_sun,iyr,iday,sec,sunr,rs
;
;  Computes unit Sun vector in geocentric rotating coodinates, using 
;  subprograms to compute inertial Sun vector and Greenwich hour angle

;	implicit real*8 (a-h,o-z)
;	real*4 sunr(3),su(3),rs
	
;	radeg = 180.d0/3.14159265359d0

;  Get unit Sun vector in geocentric inertial coordinates
	sun2000,iyr,iday,sec,su,rs

;  Get Greenwich mean sideral angle
	day = iday + sec/864.d2 
	gha2000,iyr,day,gha
	ghar = gha/!radeg

;  Transform Sun vector into geocentric rotating frame
	n = n_elements(day)
	sunr = fltarr(3,n)
	sunr(0,*) = su(0,*)*cos(ghar) + su(1,*)*sin(ghar)
	sunr(1,*) = su(1,*)*cos(ghar) - su(0,*)*sin(ghar)
	sunr(2,*) = su(2,*)

	return
	end

	pro sun2000,iyr,iday,sec,sun,rs
;
;  This subroutine computes the Sun vector in geocentric inertial 
;  (equatorial) coodinates.  It uses the model referenced in The 
;  Astronomical Almanac for 1984, Section S (Supplement) and documented
;  for the SeaWiFS Project in "Constants and Parameters for SeaWiFS
;  Mission Operations", in TBD.  The accuracy of the Sun vector is
;  approximately 0.1 arcminute.
;
;	implicit real*8 (a-h,o-z)
;	real*4 sun(3),rs
	xk=0.0056932		;Constant of aberration 
	imon=1
;	common nutcm,dpsi,eps,nutime
;	radeg = 180.d0/3.14159265359d0

;   Compute floating point days since Jan 1.5, 2000 
;    Note that the Julian day starts at noon on the specified date
	jd,iyr,imon,iday,jday
	t = jday - 2451545.0d0 + (sec-43200.d0)/86400.d0
	n = n_elements(t)
	sun=fltarr(3,n)

;  Compute solar ephemeris parameters
	ephparms,t,xls,gs,xlm,omega

;  Check if need to compute nutation corrections for this day
;	nt = t
;	if (nt ne nutime) then begin
;	  nutime = nt
	  nutate,t,xls,gs,xlm,omega,dpsi,eps
;	endif

;  Compute planet mean anomalies
;   Venus Mean Anomaly 	
	g2 = 50.40828 + 1.60213022*t 	
	g2 = g2 mod 360.d0

;   Mars Mean Anomaly 		
	g4 = 19.38816 + 0.52402078*t 	
	g4 = g4 mod 360.d0

;  Jupiter Mean Anomaly 
	g5 = 20.35116 + 0.08309121*t 	
	g5 = g5 mod 360.d0

;  Compute solar distance (AU)
	rs = 1.00014-0.01671*cos(gs/!radeg)-0.00014*cos(2.*gs/!radeg)

;  Compute Geometric Solar Longitude 
	dls = 	(6893. - 4.6543463D-4*t)*sin(gs/!radeg) $
     		+ 72.*sin(2.*gs/!radeg) 		$
     		- 7.*cos((gs - g5)/!radeg) 		$
     		+ 6.*sin((xlm - xls)/!radeg) 		$ 
     		+ 5.*sin((4.*gs - 8.*g4 + 3.*g5)/!radeg) $
     		- 5.*cos((2.*gs - 2.*g2)/!radeg) 	$
     		- 4.*sin((gs - g2)/!radeg) 		$
     		+ 4.*cos((4.*gs - 8.*g4 + 3.*g5)/!radeg) $
     		+ 3.*sin((2.*gs - 2.*g2)/!radeg) 	$
     		- 3.*sin(g5/!radeg) 			$
     		- 3.*sin((2.*gs - 2.*g5)/!radeg)

	xlsg = xls + dls/3600.d0

;  Compute Apparent Solar Longitude; includes corrections for nutation 
;   in longitude and velocity aberration
	xlsa = xlsg + dpsi - xk/rs

;   Compute unit Sun vector 
	sun(0,*) = cos(xlsa/!radeg)
	sun(1,*) = sin(xlsa/!radeg)*cos(eps/!radeg)
	sun(2,*) = sin(xlsa/!radeg)*sin(eps/!radeg)

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

