  pro selenolocate_oci, l1afile, lutfile, l1bfile, pcdim, psdim, board_id, distcorr, evtime, scana, scans, solz, qflg, iret, EPHFILE=ephfile

; IDL prototype of OCI selenolocation program

  print,'Selenolocating file ', l1bfile

; Read geolocation LUT
  read_oci_geo_lut, lutfile, geo_lut

; Open and read data from L1A file

  l1aid = ncdf_open(l1afile)

; Scan time, spin ID and HAM side
  sgid = ncdf_ncidinq(l1aid,'scan_line_attributes')
  ncdf_varget, sgid, 'scan_start_time', sstime
  ncdf_varget, sgid, 'spin_ID', spin
  k = where(spin gt 0)
  sstime = sstime(k)
  ns = n_elements(k)
  ncdf_varget, sgid, 'HAM_side', hside
  hside = hside(k)

; Navigation fields
  if (KEYWORD_SET(ephfile)) then def = 1 else def = 0

  ngid = ncdf_ncidinq(l1aid,'navigation_data')
  ncdf_varget, ngid, 'att_time', atime
  ncdf_varget, ngid, 'att_quat', quat
  ka = where(atime gt 0)
  atime = atime(ka)
  quat = quat(*,ka)

  if (NOT def) then begin
    ncdf_varget, ngid, 'orb_time', otime
    ncdf_varget, ngid, 'orb_pos', pos
    ncdf_varget, ngid, 'orb_vel', vel
    ko = where(otime gt 0)
    otime = otime(ko)
    pos = pos(*,ko)
    vel = vel(*,ko)
  endif else read_def_eph, ephfile, iyd, idyd, otime, pos, vel
  
  ncdf_varget, ngid, 'tilt', tiltin
  ncdf_varget, ngid, 'tilt_time', ttime

;  MCE telemetry
;  read_mce_tlm, l1aid, revpsec, ppr_off, secpline, board_id, mspin, ot_10us, enc_count, hamenc, rtaenc
  read_mce_tlm, l1aid, geo_lut, revpsec, ppr_off, secpline, board_id, mspin, ot_10us, enc_count, hamenc, rtaenc, iret

  if (iret eq 1) then return
  if (iret eq 2) then print, 'OCI static mode'
  if (iret eq 3) then print, 'OCI reverse scan'

; Need date from L1A file
  ncdf_attget, l1aid, 'time_coverage_start', tstart, /GLOBAL
  tstart = string(tstart)
  tsplit = strsplit(tstart,'T',/EXTRACT)
  tdate = strsplit(tsplit(0),'-',/EXTRACT)
  jd0 = jday(fix(tdate(0)), fix(tdate(1)), fix(tdate(2)))
  jdate, jd0, iyr, iday
; Check definitive ephemeris time range if used
  if (def) then begin
    otime = otime - (jd0 - jday(iyd,1,idyd))*86400
    norb = n_elements(otime)
    if ((sstime[0] lt otime[0]) OR (sstime[ns-1] gt otime[norb-1])) then begin
      print, 'Definitive ephemeris time range does not cover data time range"
      iret = 1
      return
    endif
  endif

; Get scan dimension for geolocation data arrays
  sdimid = ncdf_dimid(l1aid, 'number_of_scans')
  ncdf_diminq, l1aid, sdimid, dimnm, sdim
  sdim = min([sdim,ns])

; Check for missing times
  check_scan_times, sstime, tfl
  sfl = bytarr(sdim)
  sfl[*] = 2*tfl 

; Transform orbit (if needed) and attitude from J2000 to ECR
  if (def) then begin
    j1 = max(where(otime le sstime[0]))
    j2 = min(where(otime ge sstime[ns-1]))
    omegae = 7.29211585494d-5
    norb = j2 - j1 + 1
    otime = otime[j1:j2]
    posr = fltarr(3,norb)
    velr = fltarr(3,norb)
    j2000_to_ecr, iyr, iday, otime, ecmat
    for i=0,norb-1 do begin
      posr(*,i) = ecmat[*,*,i]#pos[*,i+j1]
      velr(*,i) = ecmat[*,*,i]#vel[*,i+j1]
      velr(0,i) = velr[0,i] + posr[1,i]*omegae
      velr(1,i) = velr[1,i] - posr[0,i]*omegae
    endfor
  endif else begin
    posr = pos/1000
    velr = vel/1000
  endelse

  natt = n_elements(atime)
  quatr = fltarr(4,natt)
  j2000_to_ecr, iyr, iday, atime, ecmat
  for i=0,natt-1 do begin
    mtoq, transpose(ecmat[*,*,i]), ecq
    qprod, ecq, quat(*,i), qt2
    quatr(*,i) = qt2
  endfor

; Get number of EV lines and offset from scan start time to EV mid-time
  mgid = ncdf_ncidinq(l1aid,'spatial_spectral_modes')
  ncdf_varget, mgid, 'spatial_zone_data_type', dtype
  ncdf_varget, mgid, 'spatial_zone_lines', lines
  ncdf_varget, mgid, 'spatial_aggregation', iagg
  get_ev, secpline, dtype, lines, iagg, pcdim, psdim, ev_toff, clines, slines, deltc, delts, iret
  if (iret lt 0) then begin
    print,'No science view in file ',l1afile
    return
  end
  evtime = sstime + ev_toff

; Interpolate orbit, attitude and tilt to scan times
  orb_interp, otime, posr, velr, evtime, posi, veli
  q_interp, atime, quatr, evtime, quati
  tilt_interp, ttime, tiltin, evtime, tilt
  tilt = tilt + geo_lut.tilt_home ; Add tilt home position to angles
  jn = where(tilt lt geo_lut.tilt_angles[0])
  if (jn[0] ne -1) then tilt[jn] = geo_lut.tilt_angles[0]
  jp = where(tilt gt geo_lut.tilt_angles[1])
  if (jp[0] ne -1) then tilt[jp] = geo_lut.tilt_angles[1]
; Set flag for tilt change 
  kt = where(tilt[1:*] ne tilt)
  if (kt[0] gt -1) then sfl[kt[1:*]] = sfl[kt[1:*]] + 1
  

; Create data arrays

  xlon = fltarr(pcdim,sdim)
  xlat = fltarr(pcdim,sdim)
  solz = intarr(pcdim,sdim)
  sola = intarr(pcdim,sdim)
  senz = intarr(pcdim,sdim)
  sena = intarr(pcdim,sdim)
  hgt = intarr(pcdim,sdim)
  xlon[*,*] = -32767.
  xlat[*,*] = -32767.
  solz[*,*] = -32767
  sola[*,*] = -32767
  senz[*,*] = -32767
  sena[*,*] = -32767
  hgt[*,*] = -32767
;  range = intarr(pcdim,sdim)
  qflg = bytarr(pcdim,sdim)
  ang = fltarr(3,sdim)
  scana = fltarr(pcdim,sdim)
  scans = fltarr(psdim,sdim)

; Get Sun vectors
  l_sun, iyr, iday, evtime, sunr, rs
  distcorr = rs[sdim/2]^2
  au = 149597870.7d0
  
; Get Moon vectors and compute selenographic transformation
  l_moon, iyr, iday, evtime, xmoon, rmoon
  ecr_to_seleno, iyr, iday, evtime, selmat

; Geolocate each scan line

  for iscn=0,sdim-1 do begin

    if ((iscn mod 100) eq 0) then print,'Geolocating scan',iscn

; Get S/C-to-sensor matrix
    get_sc_to_oci, geo_lut, tilt(iscn), sc_to_oci

; Convert quaternion to matrix
    qtom, quati(*,iscn), qmat
    smat = sc_to_oci#qmat#transpose(selmat[*,*,iscn])

; Compute attitude angles (informational only)
    mat2rpy, posi(*,iscn), veli(*,iscn), qmat, rpy
    ang(*,iscn) = rpy

; Get scan ellipse coefficients
    posmr = posi[*,iscn] - xmoon[*,iscn]*rmoon[iscn]
    posm = selmat[*,*,iscn]#posmr
    velm = selmat[*,*,iscn]*veli[*,iscn]
    scan_ell_moon,posm,smat,coef

; Generate pointing vector and relative time arrays in instrument frame
    get_oci_vecs, pcdim, geo_lut, ev_toff, spin(iscn), hside[iscn], clines, deltc, revpsec, $ 
	ppr_off, board_id, mspin, ot_10us, enc_count, hamenc, rtaenc, pview, theta, iretp
    sfl(iscn) = sfl(iscn) + 4*iretp
    scana[*,iscn] = theta
; Also need scan angles for SWIR bands in case they are different
    get_oci_vecs, psdim, geo_lut, ev_toff, spin(iscn), hside[iscn], slines, delts, revpsec, $ 
	ppr_off, board_id, mspin, ot_10us, enc_count, hamenc, rtaenc, pviews, thetas, irets
    scans[*,iscn] = thetas
;    stop

; Get Moon-to-Sun vector in selenographic frame
    sunm = sunr[*,iscn]*rs[iscn]*au - xmoon[*,iscn]*rmoon[iscn]
    suns = selmat[*,*,iscn]#sunm/sqrt(total(sunm*sunm))

; Geolocate pixels
    oci_selenonav, posm, velm, smat, coef, suns, $
	pview, pcdim, deltc, xlt, xln, slz, sla, snz, sna, rng, qfl

    kq = where(qfl eq 0)
    if (kq[0] ne -1) then begin
; PLACEHOLDER FOR TERRAIN CORRECTION
      hgt(kq,iscn) = 0

      xlat(kq,iscn) = xlt[kq]
      xlon(kq,iscn) = xln[kq]
      solz(kq,iscn) = slz[kq]*100
      sola(kq,iscn) = sla[kq]*100
      senz(kq,iscn) = snz[kq]*100
      sena(kq,iscn) = sna[kq]*100
;    range(*,iscn) = (rng-400)*10
    endif
    qflg(*,iscn) = qfl
;    if (iscn eq 239 OR iscn eq 822) then stop

  endfor

; Generate geolocation file template
  spawn,'ncgen OCI_Level-1B_Data_Structure.cdl -o Otemp.nc'
  file_move,'Otemp.nc',l1bfile,/OVERWRITE

; Write data to file

  l1bid = ncdf_open(l1bfile,/WRITE)
  make_geo_data_objects, l1bid, sdim, pcdim, /SELENO

; Scan times, HAM sides and quality flag
  sgid = ncdf_ncidinq(l1bid,'scan_line_attributes')
  write_ncdf_data_object, sgid, 'time', evtime
  write_ncdf_data_object, sgid, 'HAM_side', hside
  write_ncdf_data_object, sgid, 'scan_quality_flags', sfl

; Navigation data
  print,'Writing navigation data'
  ngid = ncdf_ncidinq(l1bid,'navigation_data')
  write_ncdf_data_object, ngid, 'att_quat', quati
  write_ncdf_data_object, ngid, 'att_ang', ang
  write_ncdf_data_object, ngid, 'orb_pos', posi
  write_ncdf_data_object, ngid, 'orb_vel', veli
  write_ncdf_data_object, ngid, 'tilt_angle', tilt
;  write_ncdf_data_object, ngid, 'sun_ref', sunr

; Geolocation data
  print,'Writing geolocation data'
  ggid = ncdf_ncidinq(l1bid,'selenolocation_data')
  write_ncdf_data_object, ggid, 'latitude', xlat
  write_ncdf_data_object, ggid, 'longitude', xlon
  write_ncdf_data_object, ggid, 'solar_zenith', solz
  write_ncdf_data_object, ggid, 'solar_azimuth', sola
  write_ncdf_data_object, ggid, 'sensor_zenith', senz
  write_ncdf_data_object, ggid, 'sensor_azimuth', sena
  write_ncdf_data_object, ggid, 'height', hgt
  write_ncdf_data_object, ggid, 'quality_flag', qflg

; PLACEHOLDER FOR GEOSPATIAL METADATA CALCULATION

; Granule_level metadata
  write_geo_granule_metadata,l1aid,l1bid,l1bfile
  write_ncdf_geo_metadata,l1bid,pcdim,sdim,xlat,xlon
  ncdf_attput, l1bid, 'earth_sun_distance_correction', distcorr, /GLOBAL

  ncdf_close,l1aid
  ncdf_close,l1bid

  return
  end

