  pro get_oci_cal_temps, l1aid, ntemps, sstime, caltemps

; Program to read temperatures used in calibration from L1A file and interpolate to scan times
; The order of temperatures is:
;  Lens housings: blue CCD side, blue grating side, red CCD side, red grating side
;  CCDs: red right, red left, blue right, blue left
;  SDA detector 1 - 16
;  AOB 1, 2, 3, 4, 7, 8
;  MOSB near MLA

;  Calling Arguments
;
;  Name			Type 	I/O	Description
;
;  l1aid		 I	 I	NetCDF ID for L1A file
;  ntemps		 I	 I	Size of temperature array
;  sstime(*)		R*8	 I	Time of each scan (seconds of day)
;  caltemps(32,*) 	R*4	 O	Temperatures at each scan time

; Create temperature array
  nsc = n_elements(sstime)
  caltemps = fltarr(ntemps,nsc)

; Read temperature fields from L1A file
  egid = ncdf_ncidinq(l1aid,'engineering_data')
  ncdf_varget, egid, 'DAUC_temperatures',dauctemp
  ncdf_varget, egid, 'DAUC_temp_time',dauctime
  ncdf_varget, egid, 'ICDU_thermisters',icdutherm
  ncdf_varget, egid, 'TC_tlm_time',icdutime
  
; Indices of required temperatures
; MLA and lens housings (blue CCD, blue grating, red CCD, red grating)
  iicdu = [23,24,25,26,11]
; Red and blue CCDs, SDA detectors, AOB 1, 2, 3, 4, 7, 8
  idauc = [5,6,12,13,indgen(16)+14,30,31,32,33,36,37]
  
; Loop through temperatures
  dstime = sstime - sstime[0]
; ICDU thermistors
  k = where(icdutime gt 0)
  ditime = icdutime[k] - sstime[0]
  for i=0,3 do begin
    p = poly_fit(ditime,icdutherm[iicdu[i],k],1)
    caltemps[i,*] = p[0] + p[1]*dstime
  endfor
  p = poly_fit(ditime,icdutherm[iicdu[4],k],1)
  caltemps[30,*] = p[0] + p[1]*dstime
  
; DAUC temperatures
  k = where(dauctime gt 0)
  ddtime = dauctime[k] - sstime[0]
  for i=0,25 do begin 
    p = poly_fit(ddtime,dauctemp[idauc[i],k],1)
    caltemps[i+4,*] = p[0] + p[1]*dstime
  endfor

  return
  end

