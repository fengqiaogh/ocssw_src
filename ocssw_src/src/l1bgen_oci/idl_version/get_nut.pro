pro get_nut,iyr,day,xnut

;  Program to get nutation matrix for TOD-to-MOD conversion

;  Calling Arguments:

;  Name         Type    I/O     Description
;  --------     ----    ---     -----------
;  iyr          I*4      I      Year
;  iday         I*4      I      Day of year
;  xnut(3,3)    R*8      O      Nutation matrix

xnut = fltarr(3,3)
imon = 1
iday = fix(day)
t = jday(iyr,imon,iday) - 2451545.5d0 + day - iday
ephparms,t,xls,gs,xlm,omega
nutate,t,xls,gs,xlm,omega,dpsi,eps,epsm

xnut(0,0) = cos(dpsi/!radeg)
xnut(1,0) = -sin(dpsi/!radeg)*cos(epsm/!radeg)
xnut(2,0) = -sin(dpsi/!radeg)*sin(epsm/!radeg)
xnut(0,1) = sin(dpsi/!radeg)*cos(eps/!radeg)
xnut(1,1) = cos(dpsi/!radeg)*cos(eps/!radeg)*cos(epsm/!radeg) $
            + sin(eps/!radeg)*sin(epsm/!radeg)
xnut(2,1) = cos(dpsi/!radeg)*cos(eps/!radeg)*sin(epsm/!radeg) $
            - sin(eps/!radeg)*cos(epsm/!radeg)
xnut(0,2) = sin(dpsi/!radeg)*sin(eps/!radeg)
xnut(1,2) = cos(dpsi/!radeg)*sin(eps/!radeg)*cos(epsm/!radeg) $
            - cos(eps/!radeg)*sin(epsm/!radeg)
xnut(2,2) = cos(dpsi/!radeg)*sin(eps/!radeg)*sin(epsm/!radeg) $
            + cos(eps/!radeg)*cos(epsm/!radeg)

return
end
     
