  pro write_band_parms,l1bid,b1bwave,b1bf0,r1bwave,r1bf0,bm12,bm13,rm12,rm13,sm12,sm13


; Program to write spectral band parameters to an OCI L1B file

; 	Name		Type 	I/O	Description
;
;    	l1bid   	int   	 I  	File ID for L1B file
;	b1bwave(*)	R*4	 I	Array of blue band center wavelengths
;	b1bf0(*)	R*4	 I	Array of blue band solar irradiances
;	r1bwave(*)	R*4	 I	Array of red band center wavelengths
;	r1bf0(*)	R*4	 I	Array of red band solar irradiances
;	bm12(*,*,*)	R*4	 I	Array of blue band m12 coefficients
;	bm13(*,*,*)	R*4	 I	Array of blue band m13 coefficients
;	rm12(*,*,*)	R*4	 I	Array of red band m12 coefficients
;	rm13(*,*,*)	R*4	 I	Array of red band m13 coefficients
;	sm12(*,*,*)	R*4	 I	Array of SWIR band m12 coefficients
;	sm13(*,*,*)	R*4	 I	Array of SWIR band m13 coefficients


; Open the group containing the band data
  gname = 'sensor_band_parameters'
  gid = ncdf_ncidinq(l1bid,gname)

; Get blue and red band and polarization dimensions
  bdimid = ncdf_dimid(l1bid, 'blue_bands')
;  ncdf_diminq, l1aid, bdimid, dimnm, bbands
  rdimid = ncdf_dimid(l1bid, 'red_bands')
;  ncdf_diminq, l1aid, rdimid, dimnm, rbands
  pdimid = ncdf_dimid(l1bid, 'polarization_coefficients')
  hdimid = ncdf_dimid(l1bid, 'HAM_sides')

; Define data objects
  ncdf_control,l1bid, /VERBOSE
  ncdf_control,l1bid, /REDEF
  bdid = ncdf_vardef(gid, 'blue_wavelength', bdimid, /FLOAT) 
  bfid = ncdf_vardef(gid, 'blue_solar_irradiance', bdimid, /FLOAT) 
  rdid = ncdf_vardef(gid, 'red_wavelength', rdimid, /FLOAT)
  rfid = ncdf_vardef(gid, 'red_solar_irradiance', rdimid, /FLOAT)
  dim = [pdimid, hdimid, bdimid]
  bm12id = ncdf_vardef(gid, 'blue_m12_coef', dim, /FLOAT)
  bm13id = ncdf_vardef(gid, 'blue_m13_coef', dim, /FLOAT)
  dim = [pdimid, hdimid, rdimid]
  rm12id = ncdf_vardef(gid, 'red_m12_coef', dim, /FLOAT)
  rm13id = ncdf_vardef(gid, 'red_m13_coef', dim, /FLOAT)
  ncdf_control,l1bid, /ENDEF

; Set attributes
  ncdf_attput, gid, bdid, "_FillValue", -32767, /FLOAT
  ncdf_attput, gid, bdid, "long_name", "Band center wavelengths for bands from blue CCD", /CHAR
  ncdf_attput, gid, bdid, "valid_min", 305., /FLOAT
  ncdf_attput, gid, bdid, "valid_max", 610., /FLOAT
  ncdf_attput, gid, bdid, "units", "nm", /CHAR
  ncdf_attput, gid, bfid, "_FillValue", -32767, /FLOAT
  ncdf_attput, gid, bfid, "long_name", "Mean extraterrestrial solar irradiance at 1 astronomical unit for the wavelengths of the blue CCD", /CHAR
  ncdf_attput, gid, bfid, "valid_min", 90., /FLOAT
  ncdf_attput, gid, bfid, "valid_max", 215., /FLOAT
  ncdf_attput, gid, bfid, "units", "W m^-2 um^-1", /CHAR
  ncdf_attput, gid, bm12id, "_FillValue", -32767, /FLOAT
  ncdf_attput, gid, bm12id, "long_name", "Blue band M12/M11 polynomial coefficients", /CHAR
  ncdf_attput, gid, bm12id, "units", "dimensionless", /CHAR
  ncdf_attput, gid, bm13id, "_FillValue", -32767, /FLOAT
  ncdf_attput, gid, bm13id, "long_name", "Blue band M13/M11 polynomial coefficients", /CHAR
  ncdf_attput, gid, bm13id, "units", "dimensionless", /CHAR
  ncdf_attput, gid, rdid, "_FillValue", -32767, /FLOAT
  ncdf_attput, gid, rdid, "long_name", "Band center wavelengths for bands from red CCD", /CHAR
  ncdf_attput, gid, rdid, "valid_min", 595., /FLOAT
  ncdf_attput, gid, rdid, "valid_max", 900., /FLOAT
  ncdf_attput, gid, rdid, "units", "nm", /CHAR
  ncdf_attput, gid, rfid, "_FillValue", -32767, /FLOAT
  ncdf_attput, gid, rfid, "long_name", "Mean extraterrestrial solar irradiance at 1 astronomical unit for the wavelengths of the red CCD", /CHAR
  ncdf_attput, gid, rfid, "valid_min", 90., /FLOAT
  ncdf_attput, gid, rfid, "valid_max", 180., /FLOAT
  ncdf_attput, gid, rfid, "units", "W m^-2 um^-1", /CHAR
  ncdf_attput, gid, rm12id, "_FillValue", -32767, /FLOAT
  ncdf_attput, gid, rm12id, "long_name", "Red band M12/M11 polynomial coefficients", /CHAR
  ncdf_attput, gid, rm12id, "units", "dimensionless", /CHAR
  ncdf_attput, gid, rm13id, "_FillValue", -32767, /FLOAT
  ncdf_attput, gid, rm13id, "long_name", "Red band M13/M11 polynomial coefficients", /CHAR
  ncdf_attput, gid, rm13id, "units", "dimensionless", /CHAR

; Write data
  ncdf_varput, gid, bdid, b1bwave
  ncdf_varput, gid, rdid, r1bwave
  ncdf_varput, gid, bm12id, bm12
  ncdf_varput, gid, bm13id, bm13
  ncdf_varput, gid, bfid, b1bf0
  ncdf_varput, gid, rfid, r1bf0
  ncdf_varput, gid, rm12id, rm12
  ncdf_varput, gid, rm13id, rm13
  swave = [939.716, 1038.317, 1250.375, 1248.550, 1378.169, 1619.626, 1618.036, 2130.593, 2258.432]
  swdid = ncdf_varid( gid, 'SWIR_wavelength')
  ncdf_varput, gid, swdid, swave
  spass = [45, 80, 30, 30, 15, 75, 75, 50, 75]
  sbdid = ncdf_varid( gid, 'SWIR_bandpass')
  ncdf_varput, gid, sbdid, spass
  sf0 = [81.959, 67.021, 44.350, 44.490, 35.551, 23.438, 23.476, 9.122, 7.427]
  sfdid = ncdf_varid( gid, 'SWIR_solar_irradiance')
  ncdf_varput, gid, sfdid, sf0
  sm12id = ncdf_varid(gid, 'SWIR_m12_coef')
  ncdf_varput, gid, sm12id, sm12
  sm13id = ncdf_varid(gid, 'SWIR_m13_coef')
  ncdf_varput, gid, sm13id, sm13

  return
  end
