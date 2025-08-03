  pro make_geo_data_objects, fid, nscn, ncps, SELENO=seleno

; IDL procedure to create the dimensions and NetCDF data objects for the geolocation
;  data for an OCI L1B file

;	Arguments
;     
;     	Name    Type 	I/O 	Description
;     	----	---- 	--- 	-----------
;       fid     int      I      File ID for the open L1B file
;	nscn	int	 I	Number of scans in L1A file
;	ncps	int	 I	Number of CCD band pixels 

  dnames = ['latitude','longitude','height','sensor_azimuth','sensor_zenith','solar_azimuth','solar_zenith','quality_flag']
  dimnames = ['number_of_scans','ccd_pixels']
  dimsizes = [nscn, ncps]

  dimids = intarr(2)
  dimids(0) = ncdf_dimid(fid, dimnames(0))
;  dimids(5) = ncdf_dimid(fid, dimnames(5))

; Create pixel dimension
  ncdf_control,fid,/REDEF
  dimids(1) = ncdf_dimdef(fid, dimnames(1), dimsizes(1))

; Open the group containing the EV data
  gname = 'geolocation_data'
  gid = ncdf_groupdef(fid,gname)

; Define data objects
  dids = intarr(8)
  dim = [dimids(1), dimids(0)]
  for i=0,1 do dids(i) = ncdf_vardef(gid, dnames(i), dim, /FLOAT) 
  for i=2,6 do dids(i) = ncdf_vardef(gid, dnames(i), dim, /SHORT) 
  dids(7) = ncdf_vardef(gid, dnames(7), dim, /BYTE) 

  ncdf_control,fid,/ENDEF

; Set attributes
;  Variable specific
  if (KEYWORD_SET(seleno)) then begin
    ncdf_attput, gid, dids(0), "long_name", "Selenographic latitudes of pixel locations", /CHAR  
    ncdf_attput, gid, dids(1), "long_name", "Selenographic longitudes of pixel locations", /CHAR
  endif else begin
    ncdf_attput, gid, dids(0), "long_name", "Latitudes of pixel locations", /CHAR
    ncdf_attput, gid, dids(1), "long_name", "Longitudes of pixel locations", /CHAR
  endelse
  ncdf_attput, gid, dids(0), "units", "degrees north", /CHAR
  ncdf_attput, gid, dids(0), "valid_min", -90., /FLOAT
  ncdf_attput, gid, dids(0), "valid_max", 90., /FLOAT
  ncdf_attput, gid, dids(1), "units", "degrees east", /CHAR
  ncdf_attput, gid, dids(1), "valid_min", -90., /FLOAT
  ncdf_attput, gid, dids(1), "valid_max", 90., /FLOAT
  ncdf_attput, gid, dids(2), "long_name", "Terrain height at pixel locations", /CHAR
  ncdf_attput, gid, dids(2), "units", "meters", /CHAR
  ncdf_attput, gid, dids(2), "valid_min", -1000, /SHORT
  ncdf_attput, gid, dids(2), "valid_max", 10000, /SHORT
  ncdf_attput, gid, dids(3), "long_name", "Sensor azimuth angle at pixel locations", /CHAR
  ncdf_attput, gid, dids(3), "valid_min", -18000, /SHORT
  ncdf_attput, gid, dids(4), "long_name", "Sensor zenith angle at pixel locations", /CHAR
  ncdf_attput, gid, dids(4), "valid_min", 0, /SHORT
  ncdf_attput, gid, dids(5), "long_name", "Solar azimuth angle at pixel locations", /CHAR
  ncdf_attput, gid, dids(5), "valid_min", -18000, /SHORT
  ncdf_attput, gid, dids(6), "long_name", "Solar zenith angle at pixel locations", /CHAR
  ncdf_attput, gid, dids(6), "valid_min", 0, /SHORT
  ncdf_attput, gid, dids(7), "long_name", "Geolocation pixel quality flags", /CHAR
  ncdf_attput, gid, dids(7), "flag_masks", [1, 2, 4], /BYTE
  ncdf_attput, gid, dids(7), "flag_meanings", "Off_Earth Input_invalid Terrain_bad", /CHAR

;  Common attributes
  for i=0,1 do ncdf_attput, gid, dids(i), "_FillValue", -32767., /FLOAT
  for i=3,6 do begin
    ncdf_attput, gid, dids(i), "units", "degrees", /CHAR
    ncdf_attput, gid, dids(i), "_FillValue", -32767, /SHORT
    ncdf_attput, gid, dids(i), "valid_max", 18000, /SHORT
    ncdf_attput, gid, dids(i), "scale_factor", 0.01, /FLOAT
    ncdf_attput, gid, dids(i), "add_offset", 0.0, /FLOAT
  endfor

return
end


