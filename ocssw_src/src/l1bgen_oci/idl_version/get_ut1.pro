pro get_ut1,iyr,idy,ut1utc

; Program to get UT1-UTC from utcpole.dat file

mjd = jday(iyr,1,idy) - 2400000
tjd = string(mjd,format='(i5)')

spawn,strjoin(['grep',tjd,'$OCVARROOT/modis/utcpole.dat > utcpole.tmp'],' ')

get_lun,ulun
openr,ulun,'utcpole.tmp'
ijd = 0L
while (ijd ne mjd and NOT EOF(ulun)) do readf,ulun,ijd,xn,xe,yn,ye,ut1utc,ute

if (ijd ne mjd) then ut1utc = 0.0
close,ulun
free_lun,ulun

return
end
