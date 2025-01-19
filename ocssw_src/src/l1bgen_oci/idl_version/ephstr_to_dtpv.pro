  pro ephstr_to_dtpv,ephstr,iyr,iday,sec,p,v
  
; Program to parse time, position and velocity from one line of a PACE ephemeris file
  
; Parse date/time and convert to YDS
  tpv = strsplit(ephstr,' ',/extract)
  timstr2yds, tpv[0], iyr, iday, sec
;  datetpv = strsplit(ephstr,'T',/extract)
;  date = strsplit(datetpv[0],'-',/extract)
;  iyr = uint(date[0])
;  mm = uint(date[1])
;  id = uint(date[2])
; Parse and convert time
;  tpv = strsplit(datetpv[1],' ',/extract)
;  time = strsplit(tpv[0],':',/extract)
;  ih = uint(time[0])
;  mn = uint(time[1])
;  sc = double(time[2])
; Extract position and velocity
  p = double(tpv[1:3])
  v = double(tpv[4:6])
  
  return
  end
  
