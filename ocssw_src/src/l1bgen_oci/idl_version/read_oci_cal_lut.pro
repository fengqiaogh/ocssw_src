  pro read_oci_cal_lut, lutid, bands, cal_lut, ldims, iret

; Program to read the OCI calibration LUT for one set of bands into a structure

;  Name		Type 	I/O	Description
;
;  lutid	string	 I	NetCDF ID of calibration LUT file
;  bands	string	 I	OCI band set ("red", "blue", "SWIR")
;  geo_lut	structure O	Structure containing LUT parameters
;  iret		Int	 O	Return code (-1 = invalid band set)

  iret = 0
; Check for valid band set
  if (bands ne 'red' AND bands ne 'blue' AND bands ne 'SWIR') then begin
    print,'Invalid band set ',bands
    iret = -1
    return
  endif

; Get dimensions
  dimnm = strjoin([bands,'_bands'])
  mdimid = ncdf_dimid(lutid,dimnm)
  ncdf_diminq, lutid, mdimid, dimnm, banddim
  mdimid = ncdf_dimid(lutid,'number_of_times')
  ncdf_diminq, lutid, mdimid, dimnm, timedim
  if (bands eq 'red' or bands eq 'blue') then begin
    mdimid = ncdf_dimid(lutid,'number_of_CCD_temperatures')
    ncdf_diminq, lutid, mdimid, dimnm,tempdim
    mcedim = 1
  endif else begin
    mdimid = ncdf_dimid(lutid,'number_of_SWIR_temperatures')
    ncdf_diminq, lutid, mdimid, dimnm,tempdim
    mcedim = 2
  endelse
  mdimid = ncdf_dimid(lutid,'number_of_T_coefficients')
  ncdf_diminq, lutid, mdimid, dimnm, tcdim
  mdimid = ncdf_dimid(lutid,'number_of_RVS_coefficients')
  ncdf_diminq, lutid, mdimid, dimnm, rvsdim
  mdimid = ncdf_dimid(lutid,'number_of_nonlinearity_coefficients')
  ncdf_diminq, lutid, mdimid, dimnm, nldim
  mdimid = ncdf_dimid(lutid,'number_of_polarization_coefficients')
  ncdf_diminq, lutid, mdimid, dimnm, poldim
  msdim = 2
;  if (bands eq 'SWIR') then mcedim = 2 else mcedim = 1 ; safe to assume this

; Create LUT structure
  cal_lut = { 	ldims:intarr(7), $
		K1:fltarr(msdim, banddim), $
		K2:fltarr(timedim, msdim, banddim), $
  		K3_coef:fltarr(tcdim, tempdim, banddim), $
  		K4_coef:fltarr(rvsdim, mcedim, msdim, banddim), $
  		K4_sca:fltarr(mcedim, msdim, banddim), $
  		K5_coef:dblarr(nldim, banddim), $
  		sat_thres:lonarr(banddim), $
  		m12_coef:fltarr(poldim, msdim, banddim), $
  		m13_coef:fltarr(poldim, msdim, banddim) }
  cal_lut.ldims = [timedim, tempdim, tcdim, rvsdim, nldim, msdim, poldim]

; Read file and populate structure
  cgid = ncdf_ncidinq(lutid,bands)
  ncdf_varget, cgid, 'K1', K1
  cal_lut.K1 = K1
  ncdf_varget, cgid, 'K2', K2
  cal_lut.K2 = K2
  ncdf_varget, cgid, 'K3_coef', K3
  cal_lut.K3_coef = K3
  ncdf_varget, cgid, 'K4_coef', K4
  cal_lut.K4_coef = K4
  ncdf_varget, cgid, 'K4_sca', K4s
  cal_lut.K4_sca = K4s
  ncdf_varget, cgid, 'K5_coef', K5
  cal_lut.K5_coef = K5
  ncdf_varget, cgid, 'sat_thres', sat
  cal_lut.sat_thres = sat
  ncdf_varget, cgid, 'm12_coef', m12
  cal_lut.m12_coef = m12
  ncdf_varget, cgid, 'm13_coef', m13
  cal_lut.m13_coef = m13

  return
  end

