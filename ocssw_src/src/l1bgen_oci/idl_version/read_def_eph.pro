  pro read_def_eph, ephfile, iyr, iday, otime, posj, velj

; Program to read a PACE predicted ephemeris file and convert vectors to ECR

  get_lun,elun
  openr,elun,ephfile
; Skip header
  head=strarr(1)
  while (head ne 'META_STOP') do begin
    readf,elun,head
  endwhile
  readf,elun,head  

; Set up arrays
;  32 hours or 1921 records per file
  numrec = 1921    
  orb = dblarr(7,numrec)
  posj = dblarr(3,numrec)
  velj = dblarr(3,numrec)
  otime = dblarr(numrec)
  p = dblarr(3)
  v = dblarr(3)
  sc = 0.d0
  daysec = 86400.d0
  ephstr = strarr(1)

; Loop through records in file (32 hours at 1-minute intervals)
  for i=0,numrec-1 do begin
    readf,elun,ephstr
    ephstr_to_dtpv,ephstr,iy,idy,sc,p,v
    jd = jday(iy,1,idy)
; Save year and day-of-year from first record
    if (i eq 0) then begin
      iyr = iy
;      jd0 = jday(iyr,1,0)
;      iday = jd - jd0
      iday = idy
    endif
;    idy = jd - jd0
; Convert hours, minutes, seconds to seconds from start of first day
    sec = sc + (idy - iday)*daysec
    if (i eq 0 OR sc eq 0.d0) then begin
      print,i,idy,sc
    endif
    posj[*,i] = p
    velj[*,i] = v
    otime[i] = sec
  endfor

  close,elun
  free_lun,elun
  end
  
