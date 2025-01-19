  pro get_oci_dark, iscn, hside, ndsc, nskp, iags, iagd, jagg, dfill, dark, dc, iret

; Program to generate dark corrections for OCI data by averaging the dark collect data
;  and correcting for bit shift/truncation if necessary

;  Name		Type 	I/O	Description
;
;  iscn		 I	 I	Scan number for which to compute dark count
;  hside(*)	 I	 I	HAM side for all scans
;  ndsc		 I	 I	Number of dark collect scans to average (should be an odd number)
;  nskp		 I	 I	Number of dark pixels to skip at the start of the collect
;  iags		 I	 I	Spatial aggregation for science data
;  iagd		 I	 I	Spatial aggregation for dark collect
;  jagg(*)	 I	 I	Spectral aggregation for each tap
;  dfill	 I	 I	Fill value for data
;  dark(*,*,*)	 I	 I	Dark collect data for granule
;  dc(*)	R*4	 O	Dark correction for each band
;  iret		 I	 O	Return code: 	0 = OK
;						1 = used adjacent scan(s)
;						-1 = no valid dark data

; Determine number of bands per tap for hyperspectral data
  ntaps = n_elements(jagg)
  nbndt = intarr(ntaps)
  if (ntaps eq 16) then begin; hyperspectral bands
    kt = where(jagg gt 0)
    nbndt(kt) = 32/jagg(kt)
  endif else nbndt = 9
  nbnds = fix(total(nbndt))

; Select data for HAM side and determine scan indices
  kh = where(hside eq hside(iscn))
  nkh = n_elements(kh)
  js = where(kh eq iscn)
  is1 = js
  is2 = js

; Check for valid dark collect data within specified range
  ndscl = ndsc
  kv = -1
  iret = -1
  while (kv(0) eq -1 AND ndscl le nkh) do begin  
    if (ndscl gt 1) then begin
      is1 = js - ndscl/2
      is2 = js + ndscl/2
;  Check for start or end of granule
      if (is1 lt 0) then begin
      	is1 = 0
      	is2 = ndscl - 1
      endif
      if (is2 ge nkh) then begin
      	is1 = nkh - ndscl
      	is2 = nkh - 1
      endif
    endif
    kv = where(dark(nskp:*,0,kh(is1:is2)) ne dfill)
;  If no valid dark data, expand scan range
    if (kv(0) eq -1) then ndscl = ndscl + 2 else iret = 0
  endwhile
  if (iret eq -1) then return
  if (ndscl gt ndsc) then iret = 1

  darks = dark(*,*,kh(is1:is2))

; Loop through taps and compute dark correction
  ibnd = 0
  dc = fltarr(nbnds)
  for i=0,ntaps-1 do begin
    if (jagg(i) gt 0) then begin
      ddiv = 1. ; division factor for bit shift
      doff = 0. ; offset for bit shift
      if (iags*jagg(i) gt 4) then begin
	ddiv = iagd*jagg(i)/4.0
	doff = (ddiv-1)/(2*ddiv)
      endif
      for j=0,nbndt(i)-1 do begin
	dc(ibnd+j) = 0.
	nv = 0
	for k=0,ndscl-1 do begin
	  kv = where(darks(nskp:*,ibnd+j,k) ne dfill)
	  if (kv(0) ne -1) then begin
  	    dc(ibnd+j) = dc(ibnd+j) + total(darks(nskp+kv,ibnd+j,k))
	    nv = nv + n_elements(kv)
	  endif
	endfor
	dc(ibnd+j) = dc(ibnd+j)/(nv*ddiv) - doff
      endfor
      ibnd = ibnd + nbndt(i)
    endif
  endfor
  
  return
  end
