pro yds2timstr,iyr,idy,secd,timstr

; IDL procedure to convert time from year, day, seconds to an
; ASCII time string

;       Arguments
;
;       Name    Type    I/O     Description
;       ----    ----    ---     -----------
;       iyr     int      I      Year
;       idy     int      I      Day of year
;       secd             I      Seconds of day (long or float)
;       timstr  sting    O      ASCII time string

; Convert year and day to year, month, day
jddate,jday(iyr,1,idy),iyr,mn,idm
cyr = string(iyr,format='(i4)')
cmn = string(mn,format='(i02)')
cdy = string(idm,format='(i02)')
cdat2 = strjoin([cyr,cmn,cdy],'-')

; Convert seconds to hours, minutes, seconds
sod2hms,secd,ih,mn,sec
chr = string(ih,format='(i02)')
cmi = string(mn,format='(i02)')
cis = string(sec,format='(f09.6)')
ctim2 = strjoin([chr,cmi,cis],':')

; Concatenate date and time strings
timstr = strjoin([cdat2,'T',ctim2,'Z'])

return
end
