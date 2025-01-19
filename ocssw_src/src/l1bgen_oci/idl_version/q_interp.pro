pro q_interp,tq,q,time,qi

;  pro q_interp
;
;  Purpose: Interpolate quaternions to scan line times
;  
;
;  Calling Arguments:
;
;  	Name         	Type  	I/O 	Description
;  	--------     	---- 	--- 	-----------
;  	tq(*)		double 	 I  	Array of quaternion times in seconds of day
;  	q(4,*)       	float	 I	Array of quaternions each time tq
;  	time(*)		double	 I	Array of time in seconds of day for every scan line
;  	qi(4,*)		float	 I	Array of interpolated quaternions


  nq = n_elements(tq)
  nlines = n_elements(time)
  qi = fltarr(4,nlines)
  e = fltarr(3)
  qri = fltarr(4)
  qri(3) = 1.0
  
  i1 = 0
  inq = -1

;  Make sure that first quaternion precedes first scan line
  k = where (time lt tq(0))
  if (k(0) ne -1) then begin          
     qi(*,k) = 0.0
     print, 'Scan line times before available quaternions'
     i1 = max(k) + 1
  endif
  
  for i=i1,nlines-1 do begin
;  Find input quaternions bracketing scan

     ind = where(tq(0:nq-2) le time(i) and tq(1:nq-1) ge time(i))
     if (ind(0) ne -1) then begin

;  Set up quaternion interpolation
	if (inq ne ind(0)) then begin
	  inq = ind(0)
	  dt = tq(inq+1) - tq(inq)
	  qinv,q(*,inq),qin
	  qprod,qin,q(*,inq+1),qr
	  e(*) = qr(0:2)
	  sto2 = sqrt(total(e*e))
	  e = e/sto2
	endif

; Interpolate quaternion to scan times
        x = (time(i) - tq(inq))/dt
	qri(0:2) = e*sto2*x
	qprod,q(*,inq),qri,qp
	qi(*,i) = qp
      endif else begin
	print, 'Scan line times after available quaternions'
      endelse

  endfor

  return
  end

