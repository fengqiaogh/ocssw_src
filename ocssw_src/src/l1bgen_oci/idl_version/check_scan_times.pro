  pro check_scan_times, sstime, spin, sfl

; Routine to check for and fill in missing scan start times.  
;  The spin IDs are used to ensure correct interpolation.
;  Missing times are flagged at the scan level

;  Calling Arguments
;
;  Name		Type 	I/O	Description
;
;  sstime(*)	R*8	 I	scan start times
;  spin(*)	I*4	 I	spin IDs
;  sfl(*)	I*2	 O	scan quality flag (1 = missing time)

  ns = n_elements(sstime)
  sfl = intarr(ns)

; Check for fill time values
  kf = where(sstime lt 0)
  if (kf(0) eq -1) then return

; Interpolate valid times to fill missing values
  nf = n_elements(kf)
  spnf = spin(kf)
  kv = where(sstime ge 0)
  nv = n_elements(kv)
  spnv = spin(kv)
  for i=0,nf-1 do begin
    sfl[kf[i]] = 1
;  Check for missing time before valid time
    if (spnf[i] lt spnv[0]) then $
	sstime[spnf[i]] = sstime[kv[0]] - (sstime[kv[1]]-sstime[kv[0]])*(spnv[0]-spnf[i])/(spnv[1]-spnv[0]) $
    else if (spnf(i) gt spnv(nv-1)) then $
	sstime[kf[i]] = sstime[kv[nv-1]] + (sstime[kv[nv-1]]-sstime[kv[nv-2]])*(spnf[i]-spnv[nv-1])/(spnv(nv-1)-spnv(nv-2)) $
    else begin
      iv = max(where(spnv lt spnf(i)))
      sstime(kf(i)) = sstime(kv(iv)) + (sstime(kv(iv+1))-sstime(kv(iv)))*(spnf(i)-spnv(iv))/(spnv(iv+1)-spnv(iv))
    endelse
  endfor

  return
  end

