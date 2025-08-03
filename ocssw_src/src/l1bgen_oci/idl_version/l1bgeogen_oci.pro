  pro l1bgeogen_oci, l1afile, lutfile, glutfile, l1bfile, CLUTFILE=clutfile, $
    RAD=rad, NOAGG=noagg, EPHFILE=ephfile, NOGEO=nogeo, SELENO = seleno

; IDL prototype of OCI radiometric calibration program

  print,'Generating file ', l1bfile

  keystring = ''  
  if (KEYWORD_SET(nogeo)) then begin
    geoloc = 0
    keystring = strjoin([keystring, ' /NOGEO']) 
  endif else geoloc = 1

  if (KEYWORD_SET(rad)) then begin 
    ociref = 0 
    keystring = strjoin([keystring, ' /RAD']) 
  endif else ociref = 1
  if (ociref AND (NOT geoloc)) then begin
    print,'Reflectances requested without geolocation - program exiting'
    return
  end
  if (KEYWORD_SET(NOAGG)) then begin
    noagg = 1
    keystring = strjoin([keystring, ' /NOAGG'])  
  endif else noagg = 0
  if (KEYWORD_SET(CLUTFILE)) then begin
    cross = 1
    keystring = strjoin([keystring, ' CLUTFILE=',clutfile])
  endif else cross = 0
  if (KEYWORD_SET(SELENO)) then begin
    seleno = 1
    ociref = 0
    keystring = strjoin([keystring, ' /RAD /SELENO'])  
  endif else seleno = 0
  if (KEYWORD_SET(ephfile)) then begin
    def = 1
    keystring = strjoin([keystring, ' EPHFILE=',ephfile])
  endif else def = 0

; Read calibration LUT
  lutid = ncdf_open(lutfile)
;  LUT structures for each band set
  read_oci_cal_lut, lutid, "blue", blue_lut
  read_oci_cal_lut, lutid, "red", red_lut
  read_oci_cal_lut, lutid, "SWIR", swir_lut
;  Common group parameters
  cgid = ncdf_ncidinq(lutid,'common')
  ncdf_varget, cgid, 'blue_wavelength',bwave
  ncdf_varget, cgid, 'blue_F0',bf0
  ncdf_varget, cgid, 'red_wavelength',rwave
  ncdf_varget, cgid, 'red_F0',rf0
  ncdf_varget, cgid, 'SWIR_wavelength',swave
  ncdf_varget, cgid, 'SWIR_F0',sf0
  ncdf_varget, cgid, 'K2t', K2t
  ncdf_varget, cgid, 'K3T', K3T
;  SWIR hysteresis parameters
  sgid = ncdf_ncidinq(lutid,'SWIR')
  ncdf_varget, sgid, 'hyst_time_const',hysttime
  ncdf_varget, sgid, 'hyst_amplitude',hystamp

  ncdf_close,lutid

; Open and read data from L1A file
  l1aid = ncdf_open(l1afile)

; Get date (change this when year and day are added to time field)
  ncdf_attget, l1aid, 'time_coverage_start', tstart, /GLOBAL
  tstart = string(tstart)
  tsplit = strsplit(tstart,'T',/EXTRACT)
  tdate = strsplit(tsplit[0],'-',/EXTRACT)
  jdate, jday(fix(tdate[0]), fix(tdate[1]), fix(tdate[2])), iyr, iday

; Get numbers of blue and red bands
  bdimid = ncdf_dimid(l1aid, 'blue_bands')
  ncdf_diminq, l1aid, bdimid, dimnm, bbands
  rdimid = ncdf_dimid(l1aid, 'red_bands')
  ncdf_diminq, l1aid, rdimid, dimnm, rbands

; Scan time, spin ID and HAM side
  sgid = ncdf_ncidinq(l1aid,'scan_line_attributes')
  ncdf_varget, sgid, 'scan_start_time', sstime
  ncdf_varget, sgid, 'spin_ID', spin
  k = where(spin gt 0)
  nscan = n_elements(k)
  sstime = sstime[k]
  ns = n_elements(k)
  ncdf_varget, sgid, 'HAM_side', hside
  hside = hside[k]

; Get spatial and spectral aggregation
  mgid = ncdf_ncidinq(l1aid,'spatial_spectral_modes')
  ncdf_varget, mgid, 'spatial_zone_data_type', dtype
  ncdf_varget, mgid, 'spatial_zone_lines', lines
  ncdf_varget, mgid, 'spatial_aggregation', iagg
  ncdf_varget, mgid, 'blue_spectral_mode', bagg
  ncdf_varget, mgid, 'red_spectral_mode', ragg
  ka = min(where(dtype ne 0 AND dtype ne 2 AND dtype ne 10))
  solcal = 0
  if (dtype(ka) eq 3) then begin
    print,'Solar cal granule'
    solcal = 1
    geoloc = 0
  endif

  if (geoloc) then begin
    if (seleno) then begin
      if (def) then $
    	selenolocate_oci, l1afile, glutfile, l1bfile, pcdim, psdim, board_id, distcorr, evtime, scana, scans, solz, qflg, iret, EPHFILE=ephfile $
      else selenolocate_oci, l1afile, glutfile, l1bfile, pcdim, psdim, board_id, distcorr, evtime, scana, scans, solz, qflg, iret
    endif else begin
      if (def) then $
        geolocate_oci, l1afile, glutfile, l1bfile, pcdim, psdim, board_id, distcorr, evtime, scana, scans, solz, qflg, iret, EPHFILE=ephfile $          
      else geolocate_oci, l1afile, glutfile, l1bfile, pcdim, psdim, board_id, distcorr, evtime, scana, scans, solz, qflg, iret
    endelse
    
    csolz = cos(solz/(!radeg*100))
    kfl = where(qflg eq 0)
    if ((kfl[0] eq -1 OR iret ne 0) AND ociref) then begin
      print,'No valid geolocation can be performed; reflectance cannot be computed'
      return
    end
  
  endif else begin
; Create L1B file and data arrays
    spawn, 'ncgen OCI_Level-1B_Data_Structure.cdl -o Otemp.nc'
    file_move, 'Otemp.nc', l1bfile, /OVERWRITE
; If solar cal granule
    if (solcal) then begin
; Compute EV time and angles for SCA
      get_oci_sca_angles, l1afile, glutfile, l1bfile, pcdim, distcorr, board_id, $
        evtime, rtain, rtaaz, solin, solaz, face, scpos, scana, iret
      if (iret ne 0) then begin
        print,'Unable to compute angles for solar cal file'
        return
      endif
      psdim = pcdim
      scans = scana
      csolz = cos(solin/!radeg)
    endif else begin
    
; Get number of EV lines and offset from scan start time to EV mid-time
; Read geo LUT
      read_oci_geo_lut,glutfile,geo_lut
      read_mce_tlm, l1aid, geo_lut, revpsec, ppr_off, secpline, board_id, mspin, ot_10us, enc_count, hamenc, rtaenc, iret
      if (iret eq 1) then begin
        print,'No valid MCE telemetry in file'
        return
      endif
      if (iret eq 2) then print, 'OCI static mode'
      if (iret eq 3) then print, 'OCI reverse scan'

      get_ev, secpline, dtype, lines, iagg, pcdim, psdim, ev_toff, clines, slines, deltc, delts, iret
      if (iret lt 0) then begin
        print,'No science collect in file ',l1afile
        ncdf_close,l1aid
        return
      endif
;  Check for and fill in missing scan times
      check_scan_times, sstime, sfl
      evtime = sstime + ev_toff
    endelse
    qflg = bytarr(pcdim,nscan)
  endelse

; Generate matrices for spectral and gain aggregation
  get_agg_mat, iagg[ka], bagg, noagg, bib, bbb, bamat, bgmat
; Crosstalk correction for blue bands only
  if (cross) then begin
    make_cross_coef_mat, 'blue', clutfile, iagg[ka], bib, bgmat, ncpix, bcmat, iret
    if (iret ne 0) then begin
      print,'Error in generating crosstalk coefficient matrix'
      return
    endif    
  endif
  get_agg_mat, iagg[ka], ragg, noagg, rib, rbb, ramat, rgmat
  if (bib ne bbands) then begin
    print,'Number of blue bands in file ',l1afile,' not consistent with spectral aggregation'
    ncdf_close,l1aid
    return
  endif else if (bib lt 4) then print,'No blue bands in file ',l1afile
  if (rib ne rbands) then begin
    print,'Number of red bands in file ',l1afile,' not consistent with spectral aggregation'
    ncdf_close,l1aid
    return
  endif else if (rib lt 4) then print,'No red bands in file ',l1afile
  swb = 9

; Get dark collect location
  kd = min(where(dtype eq 2))
  if (kd eq -1) then begin
    print,'No dark collect in file',l1afile
    ncdf_close,l1aid
    return
  end
  ldc = 0
  lds = 0
  for i=0,kd-1 do begin
    if (dtype[i] ne 0 AND dtype[i] ne 10) then begin
      ldc = ldc + lines[i]/iagg[i]
      lds = lds + lines[i]/8
    endif
  endfor
  ndc = lines[kd]/iagg[kd]
  nds = lines[kd]/8

; Generate band gain structures from LUTs, date/time and gain matrices
  if (bib ge 4) then make_oci_gains, bib, iyr, iday, evtime[0], k2t, board_id, iagg[ka], bagg, blue_lut, bgmat, bgains
  if (rib ge 4) then make_oci_gains, rib, iyr, iday, evtime[0], k2t, board_id, iagg[ka], ragg, red_lut, rgmat, rgains
  sgmat = fltarr(swb,swb)
  for i=0,swb-1 do sgmat[i,i] = 1.0
  make_oci_gains, swb, iyr, iday, evtime[0], k2t, 0, 1, 1, swir_lut, sgmat, sgains
;  stop

; Read selected temperature fields and interpolate to scan times
  ntemps = bgains.ldims[1] + sgains.ldims[1]
  get_oci_cal_temps, l1aid, ntemps, evtime, caltemps
  nctemps = bgains.ldims[1]

; Get L1A data object IDs and fill values
  egid = ncdf_ncidinq(l1aid,'science_data')
  bid = ncdf_varid( egid, 'sci_blue')
  ncdf_attget, egid, bid, '_FillValue', bfill
  rid = ncdf_varid( egid, 'sci_red')
  ncdf_attget, egid, rid, '_FillValue', rfill
  sid = ncdf_varid( egid, 'sci_SWIR')
  ncdf_attget, egid, sid, '_FillValue', sfill

; Read dark collects from science data arrays and check/flag noise spikes
  start = [ldc, 0, 0]
  if (bib ge 4) then begin
    dims = [ndc, bib, nscan]
    ncdf_varget, egid, bid, bdark, OFFSET = start, COUNT = dims
    filter_dark_noise, bib, bfill, bdark 
  endif
  if (rib ge 4) then begin
    dims = [ndc, rib, nscan]
    ncdf_varget, egid, rid, rdark, OFFSET = start, COUNT = dims
    filter_dark_noise, rib, rfill, rdark 
  endif
  start = [lds, 0, 0]
  dims = [nds, swb, nscan]
  ncdf_varget, egid, sid, sdark, OFFSET = start, COUNT = dims
  filter_dark_noise, swb, sfill, sdark 
;  stop

; Create data arrays
  l1bid = ncdf_open(l1bfile,/WRITE)
  if (ociref) then make_l1b_data_objects, l1bid, nscan, pcdim, bbb, rbb, psdim, swb, gid, nid, dids $
    else make_l1b_data_objects, l1bid, nscan, pcdim, bbb, rbb, psdim, swb, gid, nid, dids, /RAD

;  Calculate band centers and solar irradiances for aggregated hyperspectral bands
;   Irradiances at L1A spectral aggregation for reflectance calculation
  if (bib ge 4) then begin
    b1bwave = bamat#bgmat#bwave 
    b1bf0 = bamat#bgmat#bf0 
  endif else b1bwave = 0.
  if (rib ge 4) then begin
    r1bwave = ramat#rgmat#rwave 
    r1bf0 = ramat#rgmat#rf0 
  endif else r1bwave = 0.

  ndsc = 1 ; number of scans of dark data to average; will make this an input parameter
  nskp = 0 ; number of dark pixels to skip; will make this an input parameter
  nsskp = nskp*iagg[kd]/8 ; adjust SWIR skip factor in case of CCD aggregation LT 8
  bicount = [pcdim, bib, 1]
  ricount = [pcdim, rib, 1]
  sicount = [psdim, swb, 1]
  bcount = [pcdim, 1, bbb]
  rcount = [pcdim, 1, rbb]
  scount = [psdim, 1, swb]

; Calibrated data variables
  bdn = fltarr(pcdim, bib)
  bdnc = fltarr(pcdim, bib)
  rdn = fltarr(pcdim, rib)
  sdn = fltarr(psdim, swb)
  bcal = fltarr(pcdim, bib)
  rcal = fltarr(pcdim, rib)
  scal = fltarr(psdim, swb)
  bcalb = fltarr(pcdim, bbb)
  rcalb = fltarr(pcdim, rbb)
  scalb = fltarr(psdim, swb)
  bsat = bytarr(bib)
  bqual = bytarr(pcdim, bbb)
  rsat = bytarr(rib)
  rqual = bytarr(pcdim, rbb)
  ssat = bytarr(swb)
  squal = bytarr(psdim, swb)
  
; lim on solar zenith for reflectance conversion
  csolz_min = cos(88/!radeg)

; Read, calibrate and write science data
  for iscn=0,nscan-1 do begin

    if ((iscn mod 100) eq 0) then print,'Calibrating scan',iscn

;  Check for valid mirror side
    if (hside[iscn] eq 0 OR hside[iscn] eq 1) then begin

;  Get scan angle
    if (geoloc OR solcal) then begin
      thetap = scana[*,iscn]*!radeg
      thetas = scans[*,iscn]*!radeg
    endif else begin
      get_oci_vecs, pcdim, geo_lut, ev_toff, spin[iscn], hside[iscn], clines, deltc, revpsec, $
	ppr_off, board_id, mspin, ot_10us, enc_count, hamenc, rtaenc, pview, thetap, iret
      thetap = thetap*!radeg
      get_oci_vecs, psdim, geo_lut, ev_toff, spin[iscn], hside[iscn], slines, delts, revpsec, $
	ppr_off, board_id, mspin, ot_10us, enc_count, hamenc, rtaenc, sview, thetas, iret
      thetas = thetas*!radeg
    endelse
    istart = [0, 0, iscn]
    start = [0, iscn, 0]
    ncdf_varput, nid, dids[6], thetap, OFFSET = start[0:1], COUNT = bcount[0:1]
    ncdf_varput, nid, dids[7], thetas, OFFSET = start[0:1], COUNT = scount[0:1]

;  Blue bands
    if (bib ge 4) then begin

;  Read L1A science data
      ncdf_varget, egid, bid, bsci, OFFSET = istart, COUNT = bicount
            
;  Compute dark offset, correct data, and apply absolute and temporal gain and temperature correction
      get_oci_dark, iscn, hside, ndsc, nskp, iagg[ka], iagg[kd], bagg, bfill, bdark, bdc, iret
      if (ociref) then bval = where(bsci[*,0] ne bfill AND qflg[*,iscn] eq 0 AND csolz[*,iscn] ge csolz_min) $
        else bval = where(bsci[*,0] ne bfill)
      bcalb[*,*] = -32767.
      if (iret ne -1 AND bval[0] ne -1) then begin
      	for j=0,bib-1 do bdn[*,j] = bsci[*,j] - bdc[j] ; Need to save dn for linearity correction
      	
;  Compute and apply crosstalk correction for blue band only
	if (cross) then begin
	  cross_corr_oci, bib, pcdim, ncpix, bcmat, bdn, bcross
	  bdnc = bdn - bcross
	endif else bdnc = bdn

;  Compute and apply temperature correction      	
      	get_oci_temp_corr, bib, bgains, K3T[0:nctemps-1], caltemps[0:nctemps-1,iscn], k3
	for j=0,bib-1 do bcal[*,j] = bgains.k1k2[hside[iscn],j]*k3[j]*bdnc[*,j]

;  Compute and apply RVS and linearity
      	get_oci_rvs_corr, bib, pcdim, hside[iscn], bgains, thetap, k4
      	get_oci_lin_corr, bib, pcdim, bgains, bdn, k5
      	bcal = k5*k4*bcal
      	
;  Aggregate to L1B bands and compute reflectance
      	bcalb[bval,*] = transpose(bamat#transpose(bcal[bval,*]))
	if (ociref) then for j=0,bbb-1 do bcalb[bval,j] = bcalb[bval,j]*!pi*distcorr/(b1bf0[j]*csolz[bval,iscn])

; Check for saturation 
	bqual[*,*] = 0
	for j=0,pcdim-1 do begin
	  bsat[*] = (bdn[j,*] ge bgains.sat_thres); Check on dark-corrected dn
	  ls = where((bamat#bsat) gt 0); Aggregate results to L1B bands
	  if (ls[0] ne -1) then bqual[j,ls] = bqual[j,ls] + 1
	endfor
      endif
      
;  Output to L1B file
      ncdf_varput, gid, dids[0], bcalb, OFFSET = start, COUNT = bcount
      ncdf_varput, gid, dids[3], bqual, OFFSET = start, COUNT = bcount
    endif

;  Red bands
    if (rib ge 4) then begin

;  Read L1A science data
      ncdf_varget, egid, rid, rsci, OFFSET = istart, COUNT = ricount
            
;  Compute dark offset, correct data, and apply absolute and temporal gain and temperature correction
      get_oci_dark, iscn, hside, ndsc, nskp, iagg[ka], iagg[kd], ragg, rfill, rdark, rdc, iret
      if (ociref) then rval = where(rsci[*,0] ne rfill AND qflg[*,iscn] eq 0 AND csolz[*,iscn] ge csolz_min) $
        else rval = where(rsci[*,0] ne rfill)
      rcalb[*,*] = -32767.
      if (iret ne -1 AND rval[0] ne -1) then begin
      	get_oci_temp_corr, rib, rgains, K3T[0:nctemps-1], caltemps[0:nctemps-1,iscn], k3
      	for j=0,rib-1 do begin
	  rdn[*,j] = rsci[*,j] - rdc[j] ; Need to save dn for linearity correction
	  rcal[*,j] = k3[j]*rgains.k1k2[hside[iscn],j]*rdn[*,j]
      	endfor

;  Compute and apply RVS and linearity
      	get_oci_rvs_corr, rib, pcdim, hside[iscn], rgains, thetap, k4
      	get_oci_lin_corr, rib, pcdim, rgains, rdn, k5
      	rcal = k5*k4*rcal

;  Aggregate to L1B bands and compute reflectance
      	rcalb[rval,*] = transpose(ramat#transpose(rcal[rval,*]))
	if (ociref) then for j=0,rbb-1 do rcalb[rval,j] = rcalb[rval,j]*!pi*distcorr/(r1bf0[j]*csolz[rval,iscn])
      	
; Check for saturation
	rqual[*,*] = 0
	for j=0,pcdim-1 do begin
	  rsat[*] = (rdn[j,*] ge rgains.sat_thres); Check on dark-corrected dn
	  ls = where((ramat#rsat) gt 0); Aggregate results to L1B bands
	  if (ls[0] ne -1) then rqual[j,ls] = rqual[j,ls] + 1
	endfor 
      endif

;  Output to L1B file
      ncdf_varput, gid, dids[1], rcalb, OFFSET = start, COUNT = rcount
      ncdf_varput, gid, dids[4], rqual, OFFSET = start, COUNT = rcount
;      stop
    endif

;  SWIR bands
   
;  Read L1A science data
    dims = [psdim, swb, 1]
    ncdf_varget, egid, sid, ssci, OFFSET = istart, COUNT = sicount
            
;  Compute dark offset, correct data, and apply absolute and temporal gain and corrections
      get_oci_dark, iscn, hside, ndsc, nsskp, 1, 1, 1, sfill, sdark, sdc, iret
      scalb[*,*] = -32767.
      if (iret ne -1) then begin
;  Compute temperature, RVS and linearity
      	get_oci_temp_corr, swb, sgains, K3T[nctemps:*], caltemps[nctemps:*,iscn], k3
      	get_oci_rvs_corr, swb, psdim, hside[iscn], sgains, thetas, k4
      	for j=0,swb-1 do sdn[*,j] = ssci[*,j] - sdc[j] 
      	get_oci_lin_corr, swb, psdim, sgains, sdn, k5
;  Interpolate solar zenith if needed
	if (ociref) then begin
	  csolzs = fltarr(psdim)
	  if (pcdim gt psdim) then begin
	    for i=0,psdim-1 do begin
	      j=where((thetas[i] ge thetap) AND (thetas[i] lt thetap[1:*]))
	      csolzs[i] = ((thetas[i] - thetap[j])*csolz[j+1,iscn] + (thetap[j+1] - thetas[i])*csolz[j,iscn])/(thetap[j+1] - thetap[j])
	    endfor
	  endif else csolzs[*] = csolz[*,iscn]
	endif
      	for j=0,swb-1 do begin
      	
          if (ociref) then sval = where(ssci[*,j] ne sfill AND qflg[*,iscn] eq 0 AND csolz[*,iscn] ge csolz_min) $
       		else sval = where(ssci[*,j] ne sfill); May be different for every SWIR band
          if (sval[0] ne -1) then begin
;	    sdn[sval,j] = ssci[sval,j] - sdc[j] ; Need to save dn for linearity correction
	    get_hyst_corr, hysttime[*,j], hystamp[*,j], sdn[sval,j], thetas[sval], hyst
	    scal[sval,j] = k3[j]*sgains.k1k2[hside[iscn],j]*(sdn[sval,j]-hyst)
      	    scal[sval,j] = k5[sval,j]*k4[sval,j]*scal[sval,j]
	    if (ociref) then scalb[sval,j] = scal[sval,j]*!pi*distcorr/(sf0[j]*csolzs[sval]) else scalb[sval,j] = scal[sval,j]
	  endif
      	endfor
;	stop

; Check for saturation
	squal[*,*] = 0
	for j=0,psdim-1 do begin
	  ssat[*] = (ssci[j,*] ge sgains.sat_thres); Check on uncorrected dn
	  squal[j,*] = squal[j,*] + ssat
	endfor
      endif

;  Output to L1B file
      ncdf_varput, gid, dids[2], scalb, OFFSET = start, COUNT = scount
      ncdf_varput, gid, dids[5], squal, OFFSET = start, COUNT = scount

    endif else print,'No mirror side index for scan',iscn
;    stop

  endfor

  ncdf_close,l1aid

; Write spectral band information
;  b1bf0 = bamat#bibf0 
;  r1bf0 = ramat#ribf0
  pdim = blue_lut.ldims[6]
  b1bm12 = fltarr(pdim,2,bbb)
  b1bm13 = fltarr(pdim,2,bbb)
  r1bm12 = fltarr(pdim,2,rbb)
  r1bm13 = fltarr(pdim,2,rbb)
  for im = 0,1 do begin
    for ip = 0,pdim-1 do begin
      b1bm12[ip,im,*] = bamat#bgmat#transpose(blue_lut.m12_coef[ip,im,*])
      b1bm13[ip,im,*] = bamat#bgmat#transpose(blue_lut.m13_coef[ip,im,*])
      r1bm12[ip,im,*] = ramat#rgmat#transpose(red_lut.m12_coef[ip,im,*])
      r1bm13[ip,im,*] = ramat#rgmat#transpose(red_lut.m13_coef[ip,im,*])
    endfor
  endfor
  write_band_parms, l1bid, b1bwave, b1bf0, r1bwave, r1bf0, swave, sf0, b1bm12, $
    b1bm13, r1bm12, r1bm13, swir_lut.m12_coef, swir_lut.m13_coef 

; Write scan-level metadata
  mgid = ncdf_ncidinq(l1bid,'scan_line_attributes')
;  Scan time
  stid = ncdf_varid(mgid,'time')
  ncdf_varput, mgid, stid, evtime
;  Store date as attribute
  yds2timstr,iyr,iday,0.d0,tref
  refstr = strjoin(['seconds since ',tref])
  ncdf_attput, mgid, stid, 'units', refstr, /CHAR
;  HAM side
  msid = ncdf_varid(mgid,'HAM_side')
  ncdf_varput, mgid, msid, hside
  
; Write granule metadata
  write_oci_l1b_granule_metadata,l1bid,iyr,iday,evtime[0],evtime[nscan-1],l1bfile,l1afile,lutfile,glutfile,keystring

; Close L1B file
  print,'Closing file ', l1bfile
  ncdf_close,l1bid

  return
  end


