pro write_ncdf_geo_metadata,fid,pdim,sdim,xlat,xlon

; IDL procedure to generate and write granule-level geolocation attributes
;  to a NetCDF file

;	Arguments
;     
;     	Name    Type	I/O 	Description
;     	----	----	--- 	-----------
;       fid   	int   	 I  	File ID
;	pdim	int	 I	Pixel dimension (number of columns)
;	sdim	int	 I	Scan dimension (number of rows)
;	xlat	float	 I	Array of geolocation latitudes
; 	xlon	float	 I	Array of geolocation longitudes	

; Get lat/lon min/max and write to file
  k = where(xlat ge -90 AND xlon gt -180)
  xlamin = min(xlat[k])
  xlamax = max(xlat[k])
  xlomin = min(xlon[k])
  xlomax = max(xlon[k])
; Check for dateline crossing, non-polar granule
  if (xlomin lt -179.9 AND xlomax gt 179.9 AND xlamin gt -89.9 AND xlamax lt 89.9) then begin
    kp = where(xlon[k] ge 0)
    xlomin = min(xlon[k[kp]])
    kn = where(xlon[k] lt 0)
    xlomax = max(xlon[k[kn]])
  endif
    
  ncdf_attput,fid,'geospatial_lat_min',xlamin,/GLOBAL,/FLOAT  
  ncdf_attput,fid,'geospatial_lat_max',xlamax,/GLOBAL,/FLOAT  
  ncdf_attput,fid,'geospatial_lon_min',xlomin,/GLOBAL,/FLOAT  
  ncdf_attput,fid,'geospatial_lon_max',xlomax,/GLOBAL,/FLOAT  

; Get g-ring and write to file (need to implement a more robust search)
  grlat = [xlat[pdim-1,0],xlat[0,0],xlat[0,sdim-1],xlat[pdim-1,sdim-1]]  
  grlon = [xlon[pdim-1,0],xlon[0,0],xlon[0,sdim-1],xlon[pdim-1,sdim-1]]  
  ncdf_attput,fid,'gringpointlongitude',grlon,/GLOBAL,/FLOAT
  ncdf_attput,fid,'gringpointlatitude',grlat,/GLOBAL,/FLOAT

return

end







