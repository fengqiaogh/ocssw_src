	pro jddate,jd,i,j,k
;
;	This routine computes the calendar date corresponding to
;	a given Julian day.  This code was brazenly copied from
;	a routine written by Myron Shear for CSC on Julian Day 1.
;	
;	ARGUMENT	TYPE	I/O	DESCRIPTION	
;	__________________________________________________________
;	 JD		I*4	 I	Julian Day (reference Jan 1, 4713 BC)
;	 I		I*4	 O	Year 
;	 J		I*4	 O	Month	
;	 K		I*4	 0	Day of Month

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
