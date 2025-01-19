  pro get_oci_vecs, pdim, geo_lut, ev_toff, spin, hside, slines, delt, revpsec, $
	ppr_off, board_id, mspin, ot_10us, enc_count, hamenc, rtaenc, pview, theta, iret

; This program generates the OCI Earth view vectors for one spin.  
;  It uses MCE telemetry and encoder data.  Further refinements will be made
;  as the instrument optics model and test results become available.  
; Reference: "Use of OCI Telemetry to Determine Pixel Line-of-Sight", F. Patt, 2020-05-18

;  Calling Arguments
;
;  Name		Type 	I/O	Description
;
;  pdim		  I 	 I	Number of pixels in Earth view
;  geo_lut	struct	 I	Structure containing geolocation LUT parameters
;  ev_toff	 R*8	 I	Time offset from PPR to center of Earth view
;  spin		 I*4	 I	Spin number for science data
;  hside	  I	 I	HAM side
;  slines(pdim)	 I*2	 I	OCI line numbers corresponding to science pixel centers
;  delt(pdim)	 R*8	 I	Offset in seconds from EV mid-time to pixel times
;  revpsec	 R*8	 I	Ideal rotation rate of RTA (RPS)
;  ppr_off	 I*4	 I	PPR offset in encoder counts
;  board_id	 I*2	 I	MCE board ID
;  mspin(*)	 I*4	 I	Spin numbers for MCE data
;  ot_10us(*)	 I*4	 I	10 us time of first encoder value (not currently used)
;  enc_count(*)	 I*2	 I	Number of encoder values for each spin
;  hamenc(*,*)	 R*4	 I	HAM encoder data (offsets from ideal spin in arcseconds)
;  rtaenc(*,*)	 R*4	 I	RTA encoder data (offsets from ideal spin in arcseconds)
;  pview(3,pdim) R*4	 O	View vectors in the instrument frame
;  theta(pdim)	 R*4	 O	Scan angles for science pixels
;  iret		  I	 O	Return code (0 = OK, 1 = no MCE data for spin)

  max_enc_cts = 2L^17
  dtenc = 1.d-3
  rad2asec = !radeg*36.d2
  pview = fltarr(3,pdim)
  bd_id = board_id mod 2

; Compute scan angle corresponding to PPR
  pprang = 2*!pi*(ppr_off - geo_lut.RTA_nadir[bd_id])/max_enc_cts
  if (pprang gt !pi) then pprang = pprang - 2*!pi
  
; Compute ideal scan angles for science pixels
  toff = delt + ev_toff
  theta = pprang + 2*!pi*revpsec*toff

; Interpolate encoder data to pixel times and add to scan angles
  thetacor = fltarr(pdim)
  iret = 0
  isp = where(mspin eq spin)
  if (isp[0] eq -1) then begin
    print,' No MCE encoder data for spin',spin
    iret = 1
  endif else begin
    tenc = dindgen(enc_count[isp[0]])*dtenc ; Encoder sample times at 1 KHz
    ip = 0
    while (ip lt pdim) do begin
      ke = where(tenc le toff[ip] AND tenc[1:*] gt toff[ip])
;  There will be multiple science pixels per encoder sample, so do all at once
      jp = where(toff ge tenc[ke[0]] AND toff lt tenc[ke[0]+1])
;      if (jp[0] eq -1) then stop
      ft = (toff[jp] - tenc[ke[0]])/dtenc
      thetacor[jp] = (1-ft)*rtaenc[ke[0],isp[0]] + ft*rtaenc[ke[0]+1,isp[0]] $
	- ((1-ft)*hamenc[ke[0],isp[0]] + ft*hamenc[ke[0]+1,isp[0]] + geo_lut.HAM_CT_angles[hside])*0.236
      njp = n_elements(jp)
      ip = ip + njp
    endwhile
  endelse

; Calculate planarity deviations and view vectors
  ascan = poly(theta,geo_lut.as_planarity)
  atrack = poly(theta,geo_lut.at_planarity)
  theta = theta - thetacor/rad2asec
  pview[0,*] = -sin(atrack/rad2asec)
  thetac = theta - ascan/rad2asec
  pview(1,*) = sin(thetac)
  pview(2,*) = cos(thetac)
;  stop

  return
  end
