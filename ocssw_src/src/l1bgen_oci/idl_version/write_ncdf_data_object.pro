pro write_ncdf_data_object,gid,dname,data,start

; IDL procedure to write data to a NetCDF data object

;	Arguments
;     
;     	Name    Type 	  I/O 	 Description
;     	----	---- 	  --- 	 -----------
;       gid     int        I     Group ID for data object
;       dname   string     I     Data object name
;       data               I     Data array to write
; 	start   int(*)     I     Optional array of dimension offsets

; Open data object
did = ncdf_varid(gid,dname)

; Get dimensions from data array
isize = size(data)
ndimd = isize(0)
ldimd = isize(1:ndimd)

; Check for start array
if (n_elements(start) eq 0) then begin
  ncdf_varput, gid, did, data
endif else begin
  ncdf_varput, gid, did, data, OFFSET=start
endelse

return
end
