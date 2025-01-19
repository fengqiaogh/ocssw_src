  pro read_mce_tlm,l1aid,geo_lut,revpsec,ppr_off,secpline,board_id,mspin,ot_10us,enc_count,hamenc,rtaenc, iret

;  Calling Arguments
;
;  Name		Type 	I/O	Description
;
;  l1aid	Int	 I	L1A file ID from NetCDF open

;  revpsec	R*8	 O	Telescope rotation speed in revs-per-second
;  geo_lut	structure O	Structure containing LUT parameters
;  ppr_off	I*4	 O	Offset of PPR pulse (start of scan)
;  secpline	R*8	 O	Seconds per line 	
;  board_id	I*2	 O	MCE board ID
;  mspin(*)	I*2	 O	Spin IDs from MCE packets
;  ot_10us(*)	I*4	 O 	on-time_10us (determines first encoder time in spin)
;  enc_count(*) I*2	 O 	Number of vOCI_GEO_LUT_V1.ncalid encoder counts per spin
;  hamenc(200,*) R*4	 O	HAM deviations (from encoder) (arcseconds)	
;  rtaenc(200,*) R*4	 O	RTA deviations (from encoder) (arcseconds)

  egid = ncdf_ncidinq(l1aid,'engineering_data')
  ncdf_varget, egid, 'MCE_telemetry', mtlm
  ncdf_varget, egid, 'MCE_encoder_data', enc
  ncdf_varget, egid, 'MCE_spin_ID',mspin
  ncdf_varget, egid, 'DDC_telemetry',ddctlm

; Check for missing MCE telemetry
  iret = 0
  ks = where(mspin ge 0)
  if (ks[0] eq -1) then begin
    print,'No MCE telemetry in L1A file'
    iret = 1
    return
  endif

  mtsize = size(mtlm)
  ntlm = mtsize[2]
  max_enc_cts = 2L^17
;  clock = [1.36d8,1.d8]
  clock = [geo_lut.master_clock, geo_lut.mce_clock]

; Get ref_pulse_divider and compute commanded rotation rate
  ref_pulse_div = swap_endian(long(mtlm[0:7,0],0,2))
  ref_pulse_div = ref_pulse_div mod 2L^24
  ref_pulse_sel = mtlm[9,0]/128
  revpsec = clock[ref_pulse_sel]/2/max_enc_cts/(ref_pulse_div[ref_pulse_sel]/256.d0 + 1)

; Check for static or reverse scan  
  avg_step_spd = swap_endian(fix(mtlm(428:431,0),0,2))
  if (abs(avg_step_spd[0]) lt 1000) then begin ; Static mode
    iret = 2
    revpsec = 0.d0
  endif else if (avg_step_spd[0] lt 0) then begin ; Reverse spin
    revpsec = -revpsec
    iret = 3
  endif

; Get PPR offset and on-time_10us
  ppr_off = swap_endian(long(mtlm[8:11,*],0)) mod max_enc_cts
  ot_10us = swap_endian(long(mtlm[368:371,*],0,ntlm)) mod 2L^2
  
; Get MCE board ID  
  board_id = fix(mtlm[322,0]/16)

; Get TDI time and compute time increment per line
  tditime = swap_endian(fix(ddctlm[346:347],0))
  secpline = (tditime+1)/clock(0) 

; Get valid encoder count, HAM and RTA encoder data
  enc_count = mtlm[475,*]
;  enc_s = 81./2560
  hamenc = fltarr(200,ntlm)
  rtaenc = fltarr(200,ntlm)
;  hamenc[*,*] = enc[0,*,*]*enc_s
;  rtaenc[*,*] = enc[1,*,*]*enc_s
  hamenc[*,*] = enc[0,*,*]*geo_lut.HAM_enc_scale
  rtaenc[*,*] = enc[1,*,*]*geo_lut.RTA_enc_scale

  return
  end
 
