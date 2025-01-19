  pro get_oci_temp_corr, nib, gains, K3T, caltemps, k3

;  Program to compute OCI temperature corrections using selected temperatures 

;  Name		Type 	I/O	Description
;
;  nib		 I	 I	Number of bands
;  gains	struct	 I	Structure containing calibration coefficients
;  K3T(*)	R*4	 I	Array of reference temperatures
;  caltemps(*)	R*4	 I	Array of instrument temperatures
;  k3(*)	R*4	 O	Array of temperature correction factors

  tempdim = gains.ldims(1)
  tcdim = gains.ldims(2)

; 
  td = caltemps - K3T
  k3 = fltarr(nib)
  k3(*) = 1.0
  for i=0,tempdim-1 do for j=0,tcdim-1 do k3(*) = k3(*) - gains.k3_coef(j,i,*)*td(i)^(j+1)

  return
  end

  
  
