pro timstr2yds,timstr,iyr,idy,secd

; IDL procedure to convert time from an ASCII time string to 
;    year, day, seconds

;       Arguments
;
;       Name    Type    I/O     Description
;       ----    ----    ---     -----------
;       timstr  sting    I      ASCII time string (yyyy-mm-ddThh:mm:ss.sss
;       iyr     int      O      Year
;       idy     int      O      Day of year
;       secd    double   O      Seconds of day

; Separate date and time strings
datetime = strsplit(timstr,'T',/extract)

; Convert year, month, day string to year and day integers
ymd = strsplit(datetime(0),'-',/extract)
iymd = long(ymd)
if (n_elements(iymd) eq 1) then iymd = [iymd/10000, (iymd mod 10000)/100, iymd mod 100]
jdate,jday(iymd(0),iymd(1),iymd(2)),iyr,idy

; Convert hours, minutes, seconds string to seconds of day
hms = strsplit(datetime(1),':',/extract)
hms = double(hms)
if (n_elements(hms) eq 1) then hms = [hms/10000, (hms mod 10000)/100, hms mod 100]
secd = fix(hms(0))*36.d2 + fix(hms(1))*6.d1 + hms(2)

return
end
