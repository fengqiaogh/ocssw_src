pro orb_interp,torb,p,v,time,posi,veli,orbfl,ecr=ec
;  
;  pro orb_interp(to,p,v,time,pi,vi)
;
;  Purpose: Interpolate orbit position and velocity vectors to scan line times
;  
;
;  Calling Arguments:
;
;  Name         Type    I/O     Description
;  --------     ----    ---     -----------
;  torb(*)	double 	 I      Array of orbit vector times in seconds of day
;  p(3,*)       float	 I	Array of orbit position vectors for
;                                each time torb
;  v(3,*)       float	 I	Array of orbit velocity vectors for
;                                each time torb
;  time(*)	double	 I	Array of time in seconds of day 
;				 for every scan line
;  pi(3,*)	float	 O	Array of interpolated positions
;  vi(3,*)	float	 O	Array of interpolated velocities
;  orbfl(*)	int	 O	Interpolated orbit flags (0=good)
;
;
;  By: Frederick S. Patt, GSC, August 10, 1993
;
;  Notes:  Method uses cubic polynomial to match positions and velocities
;   at input data points.
;
;  Modification History:
;
;  Created IDL version from FORTRAN code.  F.S. Patt, SAIC, November 29, 2006
;
  
  norb = n_elements(torb)
  nlines = n_elements(time)
  posi = dblarr(3,nlines)
  veli = dblarr(3,nlines)
  orbfl = intarr(nlines)
  a0 = dblarr(3)&a1 = dblarr(3)&a2 = dblarr(3)&a3 = dblarr(3)&
  
  i1 = 0

;  Make sure that first orbit vector precedes first scan line
  k = where (time lt torb(0))
  if (k(0) ne -1) then begin          
     posi(*,k) = 0.0
     veli(*,k) = 0.0   
     orbfl(k) = 1
     print, 'Scan line times before available orbit data'
     i1 = max(k) + 1
  endif
  
  while (i1 lt nlines) do begin
;  Find input orbit vectors bracketing scan

     ind = where(torb(0:norb-2) le time(i1) and torb(1:norb-1) ge time(i1))
     ind = ind(0)
     if (ind(0) ne -1) then begin
;  Set up cubic interpolation
        dt = torb(ind+1) - torb(ind)
        dt = dt(0)
        a0(*) = p(*,ind)
        a1(*) = v(*,ind)*dt
        if (dt ge 3) then begin
          a2(*) = 3.d0*p(*,ind+1) - 3.d0*p(*,ind) - 2.d0*v(*,ind)*dt - v(*,ind+1)*dt
          a3(*) = 2.d0*p(*,ind) - 2.d0*p(*,ind+1) + v(*,ind)*dt + v(*,ind+1)*dt
        endif else begin
          a2(*) = (v(*,ind+1)-v(*,ind))*dt/2
          a3(*) = 0.d0
        endelse

;  Interpolate orbit position and velocity components to scan line time
        ii = where (time ge torb(ind) and time le torb(ind+1))
        x = (time(ii) - torb(ind))/dt
        x2 = x*x
        x3 = x2*x
        for j=0,2 do begin
           posi(j,ii) = a0(j) + a1(j)*x + a2(j)*x2 + a3(j)*x3
           veli(j,ii) = (a1(j) + 2.*a2(j)*x + 3.*a3(j)*x2)/dt
        endfor
        orbfl(ii) = 0
        i1 = max(ii) + 1
     endif else begin
        posi(*,i1:nlines-1) = 0.0
        veli(*,i1:nlines-1) = 0.0
        orbfl(ii) = 1
        i1 = nlines
        print, 'Scan line times after available orbit data'
     endelse
     
  endwhile
  
  return
  
end	
