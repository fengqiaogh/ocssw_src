  pro make_l1b_data_objects, fid, nscn, ncps, nbbs, nrbs, nsps, nswb, gid, nid, dids, RAD = rad

; IDL procedure to create the dimensions and NetCDF data objects for the science data
;  and quality flags for an OCI L1B file

;	Arguments
;     
;     	Name    Type 	I/O 	Description
;     	----	---- 	--- 	-----------
;       fid     int      I      File ID for the open L1B file
;	nscn	int	 I	Number of scans in L1A file
;	ncps	int	 I	Number of CCD band pixels 
;	nbbs	int	 I	Number of blue bands
;	nrbs	int	 I	Number of red bands
;	nsps	int	 I	Number of SWIR band pixels
;	nswb	int	 I	Number of SWIR bands
;    	gid     int      O      Group ID for calibrated data
;    	nid     int      O      Group ID for navigation data
;    	dids    long(8)  O      Array of dataset IDs (blue, red, SWIR, blue flags, red flags, SWIR flags, CCD scan angles, SWIR scan angles)

  if (KEYWORD_SET(rad)) then begin
    dnames = ['Lt_blue','Lt_red','Lt_SWIR','qual_blue','qual_red','qual_SWIR','CCD_scan_angles','SWIR_scan_angles']
    runits = "W m-2 sr-1"
    lparm = "Radiance"
    vmax = [800., 700., 300.]
  endif else begin
    dnames = ['rhot_blue','rhot_red','rhot_SWIR','qual_blue','qual_red','qual_SWIR','CCD_scan_angles','SWIR_scan_angles']
    runits = "dimensionless"
    lparm = "Reflectance"
    vmax = [1.3, 1.3, 1.3]
  endelse
  dimnames = ['number_of_scans','ccd_pixels','SWIR_pixels','blue_bands','red_bands','SWIR_bands']
  dimsizes = [nscn, ncps, nsps, nbbs, nrbs, nswb]

  dimids = intarr(6)
  dimids[0] = ncdf_dimid(fid, dimnames[0])
  dimids[5] = ncdf_dimid(fid, dimnames[5])

; Create dimensions
  ncdf_control,fid,/REDEF
  
; Check if CCD pixel dimension already exists
  dimids[1] = ncdf_dimid(fid,dimnames[1])
  if (dimids[1] eq -1) then dimids[1] = ncdf_dimdef(fid, dimnames[1], dimsizes[1])
;  for i=2,4 do dimids(i) = ncdf_dimdef(fid, dimnames(i), dimsizes(i))
  dimids[3] = ncdf_dimdef(fid, dimnames[3], dimsizes[3])
  dimids[4] = ncdf_dimdef(fid, dimnames[4], dimsizes[4])

; Open the groups
  gname = 'observation_data'
  gid = ncdf_ncidinq(fid,gname)
  nname = 'navigation_data'
  nid = ncdf_ncidinq(fid,nname)

; Define data objects
  dids = intarr(8)
  dim = [dimids[1], dimids[0], dimids[3]]
  dids[0] = ncdf_vardef(gid, dnames[0], dim, /FLOAT) 
  dids[3] = ncdf_vardef(gid, dnames[3], dim, /BYTE)
;  dim = [dimids[1], dimids[0], dimids[4]]
  dim[2] = dimids[4]
  dids[1] = ncdf_vardef(gid, dnames[1], dim, /FLOAT)
  dids[4] = ncdf_vardef(gid, dnames[4], dim, /BYTE)
  dim2 = dim[0:1]
  dids[6] = ncdf_vardef(nid, dnames[6], dim2, /FLOAT)
;  dim = [dimids[2], dimids[0], dimids[5]]
  dim[2] = dimids[5]
  if (nsps ne ncps) then begin
    dimids[2] = ncdf_dimdef(fid, dimnames[2], dimsizes[2])
    dim[0] = dimids[2]
  endif
  dids[2] = ncdf_vardef(gid, dnames[2], dim, /FLOAT)
  dids[5] = ncdf_vardef(gid, dnames[5], dim, /BYTE)
  dim2 = dim[0:1]
  dids[7] = ncdf_vardef(nid, dnames[7], dim2, /FLOAT)

  ncdf_control,fid,/ENDEF

; Set attributes
;  Band-set specific
  ncdf_attput, gid, dids[0], "long_name", strjoin(["Top of Atmosphere Blue Band ",lparm]), /CHAR
  ncdf_attput, gid, dids[1], "long_name", strjoin(["Top of Atmosphere Red Band ",lparm]), /CHAR
  ncdf_attput, gid, dids[2], "long_name", strjoin(["Top of Atmosphere SWIR Band ",lparm]), /CHAR
  ncdf_attput, gid, dids[3], "long_name", "Blue Band Quality Flag", /CHAR
  ncdf_attput, gid, dids[4], "long_name", "Red Band Quality Flag", /CHAR
  ncdf_attput, gid, dids[5], "long_name", "SWIR Band Quality Flag", /CHAR
  ncdf_attput, nid, dids[6], "long_name", "Scan angles for blue and red band science pixels", /CHAR
  ncdf_attput, nid, dids[7], "long_name", "Scan angles for SWIR band science pixels", /CHAR

;  Common attributes
  for i=0,2 do begin
    ncdf_attput, gid, dids(i), "units", runits, /CHAR
    ncdf_attput, gid, dids(i), "_FillValue", -32767., /FLOAT
    ncdf_attput, gid, dids(i), "valid_min", 0., /FLOAT
    ncdf_attput, gid, dids(i), "valid_max", vmax[i], /FLOAT
    ncdf_attput, gid, dids(i+3), "_FillValue", 255, /BYTE
    ncdf_attput, gid, dids(i+3), "flag_masks", [1, 2, 4], /BYTE
    ncdf_attput, gid, dids(i+3), "flag_meanings", "saturation", /CHAR
  endfor
  for i=6,7 do begin
    ncdf_attput, nid, dids(i), "units", "degrees", /CHAR
    ncdf_attput, nid, dids(i), "_FillValue", -32767., /FLOAT
    ncdf_attput, nid, dids(i), "valid_min", -110., /FLOAT
    ncdf_attput, nid, dids(i), "valid_max", 250., /FLOAT
  endfor

return
end

