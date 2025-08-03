  pro filter_dark_noise, nbnds, dfill, dark

; Program to check OCI dark collect data for noise spikes and set to fill values

;  Name	Type 	I/O	Description
;
;  nbnds	 I	 I	Number of bands
;  dfill	 I	 I	Fill value for data
;  dark(*,*,*)	 I	I/O	Dark collect data for granule

  rejfac = 4.5 ; applied to IQR, equivalent to 6 sigma
; Loop through bands
  for i=0,nbnds-1 do begin

;  Get median and IAR
    dtmp = dark[*,i,*]
    dm = median(dtmp)
    di = iqr(dtmp)
    kn = where(abs(dtmp - dm) gt rejfac*di)
    if (kn[0] ne -1) then dtmp[kn] = dfill
    dark[*,i,*] = dtmp
    
  endfor
  
  return
  end

