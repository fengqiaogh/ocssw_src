  pro get_oci_sca_angles, l1afile, lutfile, l1bfile, pcdim, distcorr, board_id, $
    evtime, rtain, rtaaz, solin, solaz, face, scpos, scan_angles, iret
  
; IDL program to determine incidence and azimuth angles of the RTA and Sun on the SCA 
;  diffuser face during an OCI solar cal

;       Arguments
;
;       Name    	Type    I/O     Description
;       ----    	----    ---     -----------
;	l1afile		string	 I	Name and path of OCI L1A file
;	lutfile		string	 I	Name and path of OCI geolocation LUT
;	evtime		R*8	 O	Center time of Earth view 
;	rtain		R*4	 O	Incidence angle of RTA view on the diffuser (0 - 90 degrees)
;	rtaaz		R*4	 O	Azimuth angle of RTA view on the diffuser (0 - 180 degrees)
;	solin		R*4	 O	Incidence angle of Sun vector on the diffuser (0 - 90 degrees)
;	solaz		R*4	 O	Incidence angle of RTA view on the diffuser (0 - 180 degrees)
;	face		 I	 O	Diffuser face: 1 = daily bright, 2 = dim, 3 = monthly bright
;	scpos		R*4	 O	SCA position angle relative to fixed position for face
;	scan_angles	R*4	 O	Scan angles for RTA view of the diffuser
;	iret		 I	 O	Return code: 1 = failure


  print,'Processing file ', l1afile

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
  ngid = ncdf_ncidinq(l1aid,'navigation_data')
  ncdf_varget, ngid, 'att_time', atime
  ncdf_varget, ngid, 'att_quat', quat
  ka = where(atime gt 0)
  atime = atime(ka)
  quat = quat(*,ka)
  ncdf_varget, ngid, 'orb_time', otime
  ncdf_varget, ngid, 'orb_pos', pos
  ncdf_varget, ngid, 'orb_vel', vel
  ko = where(otime gt 0)
  otime = otime(ko)
  pos = pos(*,ko)
  vel = vel(*,ko)
  ncdf_varget, ngid, 'tilt', tiltin
  ncdf_varget, ngid, 'tilt_time', ttime

;  MCE telemetry
;  read_mce_tlm, l1aid, revpsec, ppr_off, secpline, board_id, mspin, ot_10us, enc_count, hamenc, rtaenc
  read_mce_tlm, l1aid, geo_lut, revpsec, ppr_off, secpline, board_id, mspin, ot_10us, enc_count, hamenc, rtaenc, iret

  if (iret eq 1) then return

; Need date from L1A file
  ncdf_attget, l1aid, 'time_coverage_start', tstart, /GLOBAL
  tstart = string(tstart)
  tsplit = strsplit(tstart,'T',/EXTRACT)
  tdate = strsplit(tsplit(0),'-',/EXTRACT)
  jdate, jday(fix(tdate(0)), fix(tdate(1)), fix(tdate(2))), iyr, iday

; Get scan dimension for geolocation data arrays
  sdimid = ncdf_dimid(l1aid, 'number_of_scans')
  ncdf_diminq, l1aid, sdimid, dimnm, sdim
  sdim = min([sdim,ns])

; Check for missing times
  check_scan_times, sstime, sfl

; Transform attitude from J2000 to ECR
  quatr = quat
  posr = pos/1000
  velr = vel/1000
  omegae = 7.29211585494d-5

  natt = n_elements(atime)
  for i=0,natt-1 do begin
    j2000_to_ecr, iyr, iday, atime(i), ecmat
    mtoq, transpose(ecmat), ecq
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
  endif
  evtime = sstime + ev_toff

; Interpolate orbit, attitude and tilt to scan times
  orb_interp, otime, posr, velr, evtime, posi, veli
  q_interp, atime, quatr, evtime, quati
  tilt_interp, ttime, tiltin, evtime, tilt
  tilt = tilt + geo_lut.tilt_home ; Add tilt home position to angles
; Set flag for tilt change 
  kt = where(tilt[1:*] ne tilt)
  if (kt[0] gt -1) then sfl[kt[1:*]] = sfl[kt[1:*]] + 1
  

; Create data arrays

  rtain = fltarr(pcdim,sdim)
  rtaaz = fltarr(pcdim,sdim)
  solin = fltarr(pcdim,sdim)
  solaz = fltarr(pcdim,sdim)
  face = intarr(sdim)
  scpos = fltarr(sdim)
  scan_angles = fltarr(pcdim,sdim)
;  ang = fltarr(3,sdim)

; Get Sun vectors
  l_sun, iyr, iday, evtime, sunr, rs
  distcorr = rs[sdim/2]^2
  
; Get SCA positions  
  get_oci_sca_position, l1aid, scatime, scaspin, scapos

; Process each scan line

  for iscn=0,sdim-1 do begin

    if ((iscn mod 100) eq 0) then print,'Processing scan',iscn

; Get S/C-to-sensor matrix
    get_sc_to_oci, geo_lut, tilt(iscn), sc_to_oci
    
; Get OCI-to-SCA matrix
    get_oci_to_sca, spin(iscn), scaspin, scapos, oci_to_sca, fce, delta
    face[iscn] = fce
    scpos[iscn] = delta

; Convert quaternion to matrix
    qtom, quati(*,iscn), qmat
    smat = sc_to_oci#qmat

; Compute attitude angles (informational only)
;    mat2rpy, posi(*,iscn), veli(*,iscn), qmat, rpy
;    ang(*,iscn) = rpy

; Generate pointing vector and relative time arrays in instrument frame
    get_oci_vecs, pcdim, geo_lut, ev_toff, spin(iscn), hside[iscn], clines, deltc, revpsec, $ 
	ppr_off, board_id, mspin, ot_10us, enc_count, hamenc, rtaenc, pview, theta, iretp
    scan_angles[*,iscn] = theta;*!radeg
	
; Compute incidence and azimuth angles
    pvsca = -oci_to_sca#pview 
    sunsca = oci_to_sca#smat#sunr[*,iscn]
    rtain[*,iscn] = !radeg*acos(pvsca[2,*])
    rtaaz[*,iscn] = !radeg*abs(atan(-pvsca[1,*],-pvsca[0,*]))
    solin[*,iscn] = !radeg*acos(sunsca[2])
    solaz[*,iscn] = !radeg*abs(atan(-sunsca[1],-sunsca[0]))
;    stop

  endfor

  ncdf_close,l1aid

; Write data to L1B file

  if (l1bfile ne '') then begin
  l1bid = ncdf_open(l1bfile,/WRITE)
  make_geo_data_objects, l1bid, sdim, pcdim

; Scan times, HAM sides and quality flag
  sgid = ncdf_ncidinq(l1bid,'scan_line_attributes')
  write_ncdf_data_object, sgid, 'scan_quality_flags', sfl

; Navigation data
  print,'Writing navigation data'
  ngid = ncdf_ncidinq(l1bid,'navigation_data')
  write_ncdf_data_object, ngid, 'att_quat', quati
  write_ncdf_data_object, ngid, 'orb_pos', posi
  write_ncdf_data_object, ngid, 'orb_vel', veli
  write_ncdf_data_object, ngid, 'tilt_angle', tilt

; Geolocation data
  print,'Writing geolocation data'
  ggid = ncdf_ncidinq(l1bid,'geolocation_data')
  solz = fix(solin*100)
  write_ncdf_data_object, ggid, 'solar_zenith', solz
  sola = fix(solaz*100)
  write_ncdf_data_object, ggid, 'solar_azimuth', sola
  senz = fix(rtain*100)
  write_ncdf_data_object, ggid, 'sensor_zenith', senz
  sena = fix(rtaaz*100)
  write_ncdf_data_object, ggid, 'sensor_azimuth', sena

; PLACEHOLDER FOR GEOSPATIAL METADATA CALCULATION

; Granule_level metadata
  ncdf_attput, l1bid, 'earth_sun_distance_correction', distcorr, /GLOBAL
  ncdf_close, l1bid
  endif
;  stop

  return
  end

