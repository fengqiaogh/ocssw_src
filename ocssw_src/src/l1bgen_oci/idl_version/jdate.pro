	pro jdate,jd,i,k
;
;	This routine computes the year and day-of-year corresponding 
;	to a given Julian day.  This algorithm is designed for the 
;	period 1900 - 2100. 
;	
;	ARGUMENT	TYPE	I/O	DESCRIPTION	
;	__________________________________________________________
;	 JD		I*4	 I	Julian Day (reference Jan 1, 4713 BC)
;	 I		I*4	 O	Year 
;	 K		I*4	 0	Day of Year
;
;	Program written by:	Frederick S. Patt
;				General Sciences Corporation
;				May 12, 1993

;  	Compute days since January 0, 1900
	l = jd - 2415020

; 	Compute years since 1900
	i = 4*l/1461

;	Compute day-of-year
	k = l - 1461*(i-1)/4 - 365 

; 	Add first two digits of year
	i = i + 1900
	return
	end	
