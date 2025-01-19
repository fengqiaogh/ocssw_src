pro write_geo_granule_metadata,l1aid,geoid,geofile

; IDL procedure to generate and write granule-level attributes
;  to a Hawkeye geolocation file

;	Arguments
;     
;     	Name    Type	I/O 	Description
;     	----	----	--- 	-----------
;       l1aid   int   	 I  	File ID for L1A file
;       geoid   int   	 I     	File ID for geolocation file
;       geofile string	 I	Geolocation file name

 print,'Writing granule metadata'

; Get start and end times from L1A file and write to geo file
 ncdf_attget, l1aid, 'time_coverage_start', tstart, /GLOBAL
 ncdf_attget, l1aid, 'time_coverage_end', tend, /GLOBAL
 ncdf_attput, geoid, 'time_coverage_start', tstart, /GLOBAL, /CHAR
 ncdf_attput, geoid, 'time_coverage_end', tend, /GLOBAL, /CHAR

; Get creation time and write to file
 ct = bin_date(systime(0))
 ctims = strarr(6)
 ctims(0) = string(ct(0),format='(i4)')
 ctims(1:5) = string(ct(1:5),format='(i02)')
 cdate = strjoin(ctims(0:2),'-')
 ctime = strjoin(ctims(3:5),':')
 tcreate = strjoin([cdate,'T',ctime,'Z'])
 ncdf_attput, geoid, 'date_created', tcreate, /GLOBAL, /CHAR

; Write product file name
 ncdf_attput, geoid, 'product_name', geofile, /GLOBAL, /CHAR

return

end







