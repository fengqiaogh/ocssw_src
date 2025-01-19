  pro check_scan_times, sstime, sfl

; Routine to check for and fill in missing scan start times
;  Missing times are flagged at the scan level

;  Calling Arguments
;
;  Name		Type 	I/O	Description
;
;  sstime(*)	R*8	 I	scan start times
;  sfl(*)	I*2	 O	scan quality flag (1 = missing time)

  ns = n_elements(sstime)
  sfl = intarr(ns)

; Check for fill time values
  kf = where(sstime lt 0)
  if (kf(0) eq -1) then return

; Interpolate valid times to fill missing values
  nf = n_elements(kf)
  kv = where(sstime ge 0)
  nv = n_elements(kv)
  for i=0,nf-1 do begin
    sfl(kf(i)) = 1
;  Check for missing time before valid time
    if (kf(i) lt kv(0)) then $
	sstime(kf(i)) = sstime(kv(0)) - (sstime(kv(1))-sstime(kv(0)))*(kv(0)-kf(i))/(kv(1)-kv(0)) $
    else if (kf(i) gt kv(nv-1)) then $
	sstime(kf(i)) = sstime(kv(nv-1)) + (sstime(kv(nv-1))-sstime(kv(nv-2)))*(kf(i)-kv(nv-1))/(kv(nv-1)-kv(nv-2)) $
    else begin
      iv = max(where(kv lt kf(i)))
      sstime(kf(i)) = sstime(kv(iv)) + (sstime(kv(iv+1))-sstime(kv(iv)))*(kf(i)-kv(iv))/(kv(iv+1)-kv(iv))
    endelse
  endfor

  return
  end

