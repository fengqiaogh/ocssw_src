  pro get_hyst_corr, hysttime, hystamp, sdn, theta, hyst

;  Program to compute OCI SWIR band hysteresis corrections

;  Name		Type 	I/O	Description
;
;  hysttime(*)	R*4	 I	Exponential time constants (pixels)
;  hystamp(*)	R*4	 I	Exponential amplitudes
;  sdn(*)	R*4	 I	Dark-corrected counts
;  theta(*)	R*4	 I	Scan angles
;  hyst(*)	R*4	 O	Hysteresis correction

  nexp = n_elements(hysttime)
  npix = n_elements(sdn)
  hc = fltarr(npix,nexp)
  hyst = fltarr(npix)
  if (npix lt 2) then return
  dtheta = theta[1:*]-theta

; Compute exponential decay constants
  e = exp(-1.0/hysttime)
  
; Loop through pixels
  for i=1,npix-1 do begin
;  Loop through exponentials
    if (dtheta[i-1] lt 0.1) then begin
      for j=0,nexp-1 do begin
      	hc[i,j] = hc[i-1,j]*e[j] + sdn[i-1]*hystamp[j]
      	hyst[i] = hyst[i] + hc[i,j]
      endfor
    endif else begin
      hc[i,*] = 0.0
      hyst[i] = 0.0
    endelse
  endfor
;  stop

  return
  end
