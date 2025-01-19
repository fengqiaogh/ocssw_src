  pro get_oci_lin_corr, nib, pdim, gains, dn, k5

;  Program to compute OCI linearity correction from coefficients and dn.  The correction
;  is a 4th-order polynomial of dn


;  Name		Type 	I/O	Description
;
;  nib		 I	 I	Number of bands in L1A file
;  pdim	 I	 I	Number of science pixels per scan
;  gains	struct	 I  	Structure of calibration coefficeints
;  dn(pdim,nib) R*4	 I	Array of dark-corrected DN
;  k5(pdim,nib) R*4	 I	Array of linearity correction factors

  k5 = fltarr(pdim,nib)

;  Loop through bands
  for i=0,nib-1 do k5[*,i] = poly(dn[*,i],gains.k5_coef[*,i])

  return
  end




