pro write_oci_l1b_granule_metadata,fid,iyr,iday,stime,etime,l1bfile,l1afile,lutfile,glutfile,keystring

; IDL procedure to generate and write granule-level attributes
;  to a NetCDF file

;	Arguments
;     
;     	Name    Type	I/O 	Description
;     	----	----	--- 	-----------
;       fid     int   	 I  	File ID
;       iyr     int   	 I     	Year
;	iday	int	 I	Day of year
;	stime		 I	Start time (seconds of day)
;	etime		 I	End time
;       l1bfile string	 I	L1B file path and name
;       l1afile string	 I	L1A file path and name
;       lutfile string	 I	Calibration LUT file path and name
;       glutfile string	 I	Geolocation LUT file path and name
;       keystring string  I	Keyword string

 print,'Writing granule metadata'

; Get start and end times and write to file
 yds2timstr,iyr,iday,stime,tstart
 ncdf_attput, fid, 'time_coverage_start', tstart, /GLOBAL, /CHAR
 yds2timstr,iyr,iday,etime,tend
 ncdf_attput, fid, 'time_coverage_end', tend, /GLOBAL, /CHAR

; Get creation time and write to file
 ct = bin_date(systime(0))
 ctims = strarr(6)
 ctims(0) = string(ct(0),format='(i4)')
 ctims(1:5) = string(ct(1:5),format='(i02)')
 cdate = strjoin(ctims(0:2),'-')
 ctime = strjoin(ctims(3:5),':')
 tcreate = strjoin([cdate,'T',ctime,'Z'])
 ncdf_attput, fid, 'date_created', tcreate, /GLOBAL, /CHAR

; Create history string
 l1a = parse_file_name(l1afile)
 l1b = parse_file_name(l1bfile)
 lut = parse_file_name(lutfile)
 glut = parse_file_name(glutfile)
 history = strjoin(['l1bgeogen_oci', l1a, lut, glut, l1b, keystring],' ')
 ncdf_attput, fid, 'history', history, /GLOBAL, /CHAR

; Write product file name
 ncdf_attput, fid, 'product_name', l1b, /GLOBAL, /CHAR
 
 return

 end







