 pro j2000_to_mod,iy,idy,sec,j2mod

; Get J2000 to MOD (precession) transformation

; Arguments:
;
; Name 		Type 	I/O 	Description
; --------------------------------------------------------
; iy      	 I   	 I  	Year, four digits
; idy     	 I   	 I  	Day of year
; sec     	R*8  	 I  	Seconds of day
; j2mod(3,3)	 R 	 O  	J2000 to MOD matrix

 t = jday(iy,1,idy) - 2451545.5d0 + sec/864.d2
 t = t/36525

 j2mod = dblarr(3,3)
 zeta0 = ( 2306.2181*t + 0.302*t^2 + 0.018*t^3 )/!radeg/36.d2
 thetap = ( 2004.3109*t - 0.4266*t^2 - 0.04160*t^3 )/!radeg/36.d2
 xip = ( 2306.2181*t + 1.095*t^2 + 0.018*t^3 )/!radeg/36.d2

 j2mod(0,0) = -sin(zeta0)*sin(xip) + cos(zeta0)*cos(xip)*cos(thetap)
 j2mod(0,1) = -cos(zeta0)*sin(xip) - sin(zeta0)*cos(xip)*cos(thetap)
 j2mod(0,2) = -cos(xip)*sin(thetap)
 j2mod(1,0) = sin(zeta0)*cos(xip) + cos(zeta0)*sin(xip)*cos(thetap)
 j2mod(1,1) = cos(zeta0)*cos(xip) - sin(zeta0) * sin(xip) * cos(thetap)
 j2mod(1,2) = -sin(xip)*sin(thetap)
 j2mod(2,0) = cos(zeta0)*sin(thetap)
 j2mod(2,1) = -sin(zeta0)*sin(thetap)
 j2mod(2,2) = cos(thetap)

 return
 end


