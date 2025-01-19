pro j2000_to_ecr,iy,idy,sec,ecmat

; Get J2000 to ECEF transformation matrix

; Arguments:
;
; Name 		Type 	I/O 	Description
; --------------------------------------------------------
; iy      	 I   	 I  	Year, four digits
; idy     	 I   	 I  	Day of year
; sec(*)     	R*8  	 I  	Seconds of day
; ecmat(3,3,*)	 R 	 O  	J2000 to ECEF matrix

 daysec = 86400.d0 ; Second per day
 ns = n_elements(sec)

; Get nutation and UT1-UTC (once per run)
 get_nut,iy,idy,xnut
 get_ut1,iy,idy,ut1utc
 
; Compute Greenwich hour angle for time of day

 day = idy + (sec+ut1utc)/daysec
 gha2000,iy,day,gha
 gham = dblarr(3,3,ns)
 gham(0,0,*) = cos(gha/!radeg)
 gham(1,1,*) = cos(gha/!radeg)
 gham(2,2,*) = 1.d0
 gham(0,1,*) = sin(gha/!radeg)
 gham(1,0,*) = -sin(gha/!radeg)

 ecmat = dblarr(3,3,ns)

 for i=0,ns-1 do begin

; Get transformation from J2000 to mean-of-date inertial
   j2000_to_mod,iy,idy,sec(i),j2mod

; Combine all transformations
   ecmat(*,*,i) = gham(*,*,i)#transpose(xnut)#j2mod

 endfor

 return
 end

