  pro cross_corr_oci, nib, npix, ncpix, ccoef, banddata, cross
  
; Program to compute along-scan stray light and crosstalk for OCI data from one FPA

;  Name				Type 	I/O	Description
;
;  nib				Int	 I	Number of L1A bands in focal plane
;  npix			Int	 I 	Number of pixels in focal plane
;  ncpix			Int	 I	Number of influence pixels around each pixel
;  ccoef(nib,nib,nscpix)	R*4	 I	Influence coefficients
;  banddata(npix,nib)		R*4	 I	Band data from focal plane
;  cross(npix,nib)		R*4	 O	Crosstalk correction

; Create correction array
  cross = fltarr(npix,nib)
  
; Create band array with padding
  nco2 = ncpix/2
  bandpad = fltarr(npix + ncpix-1,nib)
  ipndx = indgen(npix)
  bandpad[nco2:nco2+npix-1,*] = banddata
  
; Loop through influence pixels
  for j=0,ncpix-1 do cross = cross + transpose(ccoef[*,*,j]#transpose(bandpad[ipndx+j,*]))

  return
  end

