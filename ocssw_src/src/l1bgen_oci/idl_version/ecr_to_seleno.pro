 pro ecr_to_seleno,iy,idy,sec,selmat

; Get ECR to selenographic transformation matrix

; Arguments:
;
; Name 		Type 	I/O 	Description
; --------------------------------------------------------
; iy      	 I   	 I  	Year, four digits
; idy     	 I   	 I  	Day of year
; sec(*)     	R*8  	 I  	Seconds of day
; selmat(3,3,*)  R 	 O  	ECR to selenographic matrix

 ximeq = -1.535 ; Inclination of the Moon's equator to the ecliptic 
 daysec = 86400.d0 ; Second per day
 ns = n_elements(sec)

; Compute Greenwich hour angle rotation for time of day
 get_ut1,iy,idy,ut1utc
 day = idy + (sec+ut1utc)/daysec
 gha2000,iy,day,gha,xlm,omega,eps
 gham = dblarr(3,3,ns)
 gham[0,0,*] = cos(gha/!radeg)
 gham[1,1,*] = cos(gha/!radeg)
 gham[2,2,*] = 1.d0
 gham[0,1,*] = sin(gha/!radeg)
 gham[1,0,*] = -sin(gha/!radeg)

; Compute TOD to ecliptic transformation
 tod2ecl = dblarr(3,3)
 tod2ecl[0,0] = 1.d0
 tod2ecl[1,1] = cos(eps[0]/!radeg)
 tod2ecl[2,2] = tod2ecl[1,1]
 tod2ecl[1,2] = sin(eps[0]/!radeg)
 tod2ecl[2,1] = -tod2ecl[1,2]
 
 selmat = dblarr(3,3,ns)
; ecl2sel = dblarr(3,3)
 fm = xlm - omega - 180

 for i=0,ns-1 do begin
; Compute ECL to selenographic transformation
   euler313, [omega[i], ximeq, fm[i]], ecl2sel
;   stop

; Combine all transformations
   selmat(*,*,i) = ecl2sel#tod2ecl#transpose(gham[*,*,i])

 endfor

; stop
 return
 end

