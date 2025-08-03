  pro get_oci_sca_position, fid, scatime, scaspin, scapos
  
; IDL procedure to read and convert the SCA position from an OCI L1A file

;       Arguments
;
;       Name    	Type    I/O     Description
;       ----    	----    ---     -----------
;	fid		string	 I	L1A file NetCDF ID
;	scatime(*)	double	 O	Array of scan times for SCA spin IDs
;	scaspin(*)	long	 O	Array of SCA telemetry spin IDs
;	scapos		double	 O	Array of SCA positions (degrees

; This routine is based on the following input from Tom Capon:

; The SCA position based on the encoder measurement can be found in the following pseudo mnemonic, and its conversion from the packet:
; SCA_AbsEncoderDegrees_conv  = 330 - (oci.mce.solcal.GtSCAAbsoluteEncoderCount / 11930464.7 )

; From the ITOS DB entries for the SCA telemetry packet:
; mnemonicName					type	byteOffset	bitOffset	fieldLengthInBits	valueLengthInBits
; oci.mce.solcal.GtSCAAbsoluteEncoderCount	U1234	348		0		32			32

; The SCA telemetry array in the OCI L1A excludes the first 16 bytes of the packet (CCSDS primary and secondary headers and spin ID),
;  so the byte offset is reduced by this amount.


; Read SCA telemetry and spin ID from L1A file
; Open file and get group ID
;  fid = ncdf_open(l1afile)
  group = 'engineering_data'
  ngid = ncdf_ncidinq(fid,group)
; Check for SCA telemetry in file
  r = ncdf_varid(ngid, 'SCA_telemetry')
  if (r eq -1) then begin
    print,'No SCA telemetry in file ',l1afile
    return
  endif
; Read data arrays
  ncdf_varget, ngid, 'SCA_telemetry', scatlm
  ncdf_varget, ngid, 'SCA_spin_ID', scaspin
  group = 'scan_line_attributes'
  ngid = ncdf_ncidinq(fid,group)
  ncdf_varget, ngid, 'scan_start_time', stime
  ncdf_varget, ngid, 'spin_ID', spin
;  ncdf_close,fid

; Extract SCA Absolute Encoder Count and convert to angle
  ksca = where(scaspin gt 0)
  nsca = n_elements(ksca)
;  ntlm = n_elements(scatlm[0,*]) ; Number of telemetry samples in file
  sca_abs_enc_count = swap_endian(ulong(scatlm[332:335,ksca],0,nsca)) ; Convert byte array to unsigned long and big- to little-endian
  sca_enc_scal = 360.d0/2.d0^32 ; Equivalent to scale factor shown above
  scapos = 330.d0 - sca_abs_enc_count*sca_enc_scal; Perform conversion from Capon
  
  scatime = dblarr(nsca)
  scatime[*] = -999.d0
  for i=0,nsca-1 do begin
    j = where(spin eq scaspin[i])
    if (j ne -1) then scatime[i] = stime[j]
  endfor
  
  k = where(scatime gt 0)
  scatime = scatime[k]
  scaspin = scaspin[k]
  scapos = scapos[k]
  
  return
  end
  
  
