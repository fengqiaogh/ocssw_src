  pro make_oci_gains, nib, iyr, iday, stime, k2t, board_id, iagg, jagg, cal_lut, gmat, gains

; Program to generate the OCI hyperspectral gains for a granule and focal plane

;  Name		Type 	I/O	Description
;
;  nib		 I	 I	Number of bands in L1A file
;  iyear	 I	 I	Granule year
;  iday		 I	 I	Granule day of year
;  stime	 I	 I	Granule start time
;  k2t(*)	 I	 I	Dates of K2 epochs (days since 1 January 2000)
;  board_id	 I	 I	MCE board ID
;  iagg		 I 	 I	Spatial aggregation factor
;  jagg(16)	 I	 I	Spectral aggregation factors
;  cal_lut	struct	 I	Input calibration LUT
;  gmat(nib,512) R*4	 I	Matrix to convert gains from unaggregated values to L1A spectral aggregation
;  gains	struct	 O	Structure of calibration coefficients for granule

; Generate gain structure

  timedim = cal_lut.ldims[0]
  tempdim = cal_lut.ldims[1]
  tcdim = cal_lut.ldims[2]
  rvsdim = cal_lut.ldims[3]
  nldim = cal_lut.ldims[4]
  msdim = cal_lut.ldims[5]
  gains = { 	ldims:intarr(7), $
		k1k2:fltarr(msdim, nib), $ ; absolute and temporal gain
  		k3_coef:fltarr(tcdim, tempdim, nib), $ ; temperature correction
  		k4_coef:fltarr(rvsdim, msdim, nib), $ ; RVS correction
  		k4_sca:fltarr(msdim, nib), $ ; RVS correction for SCA
  		k5_coef:dblarr(nldim, nib), $ ; linearity correction
  		sat_thres:lonarr(nib) } ; Saturation thresholds
  gains.ldims = cal_lut.ldims
  bd_id = board_id mod 2

  if (n_elements(jagg) eq 16) then begin; Hyperspectral bands
    hyper = 1
    iaf = intarr(nib)
    ib = 0
    for i=0,15 do begin
      if (jagg[i] gt 0) then begin
        nb = 32/jagg[i]
        iaf[ib:ib+nb-1] = 4/min([iagg*jagg[i],4])
        ib = ib + nb
      endif
    endfor
  endif else hyper = 0
  
; Mirror-side dependent gains
  for ms = 0,1 do begin

; Get temporal gain and combine with absolute gain
    d2 = jday(iyr, 1, iday) - 2451545 + stime/864.d2
    kd = max(where(d2 gt k2t))
    if (kd lt timedim-1) then begin
      ff = (d2 - k2t[kd])/(k2t[kd+1] - k2t[kd])
      k2 = cal_lut.k2[kd,ms,*]*(1.0-ff) + cal_lut.K2[kd+1,ms,*]*ff
    endif else k2 = cal_lut.K2(kd,ms,*)
    gains.K1K2[ms,*] = gmat#transpose((cal_lut.K1[ms,*]*k2))
    if (hyper) then gains.K1K2[ms,*] = gains.K1K2[ms,*]*iaf

; Generate RVS coefficents
    if (hyper) then begin
      for i=0,rvsdim-1 do gains.k4_coef[i,ms,*] = gmat#transpose(cal_lut.K4_coef[i,0,ms,*])
      gains.k4_sca[ms,*] = gmat#transpose(cal_lut.K4_sca[0,ms,*])
    endif else begin ; Select based on MCE side
      for i=0,rvsdim-1 do gains.k4_coef[i,ms,*] = gmat#transpose(cal_lut.K4_coef[i,bd_id,ms,*])
      gains.k4_sca[ms,*] = gmat#transpose(cal_lut.K4_sca[bd_id,ms,*])
    endelse

  endfor

; Generate temperature coefficients
  for i=0,tcdim-1 do for j=0,tempdim-1 do gains.k3_coef[i,j,*] = gmat#transpose(cal_lut.K3_coef[i,j,*])

; Generate linearity coefficients
  for i=0,nldim-1 do begin
    gains.k5_coef[i,*] = gmat#transpose(cal_lut.K5_coef[i,*])
    if (hyper) then gains.k5_coef[i,*] = gains.k5_coef[i,*]*iaf^i
  endfor 
  
;  Generate saturation thresholds
  gains.sat_thres = (gmat#cal_lut.sat_thres)
  if (hyper) then gains.sat_thres = gains.sat_thres/iaf

  return
  end

