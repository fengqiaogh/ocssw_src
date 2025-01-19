  pro get_oci_rvs_corr, nib, pdim, hside, gains, theta, k4

; Program to compute RVS correction from coefficients and scan angle

;  Name		Type 	I/O	Description
;
;  nib		 I	 I	Number of bands in L1A file
;  pdim		 I	 I	Number of science pixels per scan
;  hside	 I	 I	HAM side
;  gains	struct	 I  	Structure of calibration coefficeints
;  theta(*)	R*4	 I	Array of scan angles in degrees
;  k4(pdim,nib) R*4	 O	Array of RVS correction factors

  k4 = fltarr(pdim,nib)
  k4(*,*) = 1.0
  if (theta[0] gt -60 AND theta[pdim-1] lt 60) then begin
    rvsdim = gains.ldims(3)
    for i=0,rvsdim-1 do for j=0,nib-1 do k4(*,j) = k4(*,j) + gains.k4_coef(i,hside,j)*theta^(i+1)
    return
  endif else begin
    if (theta[0] le -60) then begin
      k = where(theta le -60)
      for j=0,nib-1 do k4(k,j) = gains.k4_sca(hside,j)
    endif
  endelse
   
  return
  end

