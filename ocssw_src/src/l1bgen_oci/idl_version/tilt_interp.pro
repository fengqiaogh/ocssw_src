  pro tilt_interp,tt,tin,time,tilt

;  pro tilt_interp
;
;  Purpose: Interpolate quaternions to scan line times
;  
;
;  Calling Arguments:
;
;  	Name         	Type  	I/O 	Description
;  	--------     	---- 	--- 	-----------
;  	tt(*)		double 	 I  	Array of tilt times in seconds of day
;  	tin(*)       	float	 I	Array of tilt angles
;  	time(*)		double	 I	Array of time in seconds of day for every scan line
;  	tilt(*)		float	 I	Array of interpolated tilt angles

  nt = n_elements(tt)
  nlines = n_elements(time)
  tilt = fltarr(nlines)
  
  i1 = 0
  inq = -1

;  Make sure that first tilt precedes first scan line
  k = where (time lt tt[0])
  if (k(0) ne -1) then begin          
     tilt[k] = tin[0]
     print, 'Scan line times before available tilt angles'
     i1 = max(k) + 1
  endif
  
  i = i1
  while (i lt nlines) do begin
  ;  Find input tilt times bracketing scan
    ind = where(tt[0:nt-2] le time[i] and tt[1:nt-1] ge time[i])
    if (ind[0] ne -1) then begin
      inq = ind[0]
      j = where(time ge tt[inq] AND time le tt[inq+1])
      if (tin[inq] eq tin[inq+1]) then begin
        tilt[j] = tin[inq]
      endif else begin
; Interpolate quaternion to scan times
        x = (time[j] - tt[inq])/(tt[inq+1]-tt[inq])
        tilt[j] = (1-x)*tin[inq] + x*tin[inq+1]
      endelse
    endif else begin
      print, 'Scan line times after available tilt angles'
      tilt[i:nlines-1] = tin[nt-1]
      return
    endelse
    i = i + n_elements(j)
  endwhile

  return
  end
