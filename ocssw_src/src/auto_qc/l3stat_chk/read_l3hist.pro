PRO read_l3hist, file, tot, hists, xaxes, hmin, nlo, hmax, nhi, $
     mean, sd, loparm, hiparm
;
;  read_l3hist will read the binary histogram file and
;  get the histograms of the 12 parameters 
;
;  file  input file to read
;    Outputs
;  tot    total # points found
;  hists  12 size 100 histogram arrays
;  xaxes  12 x axes for the plot
;  hmin   12 histogram minimums
;  nlo    # below histogram low for 12 params
;  hmax   12 histogram maximums
;  nhi    # above histogram high for 12 params
;  mean   size 12 mean of data
;  sd     size 12 sd of data
;  loparm size 12 lowest param value found
;  hiparm size 12 highest param value found
;
;
;  open the binary file
openr, unit, file, /get_lun

;  declare all the variables needed to read into
bins=0L
ifile=bytarr(412)
hmin=fltarr(12)
hmax=fltarr(12)
nlo=lonarr(12)
nhi=lonarr(12)
mean = fltarr(12)
sd = fltarr(12)
loparm = fltarr(12)
hiparm = fltarr(12)
hists = lonarr(100,12)
;
;  read it in
readu, unit, bins, ifile, hmin, hmax, nlo, nhi, hists, mean, sd, loparm, hiparm
print, "getting data from file:"
print, string( ifile )
;
;  get the x axes
xaxes = fltarr(100,12)
sep = ( hmax - hmin ) / 100
for  i = 0, bins - 1  do $
  BEGIN
  xaxes(i, * ) = hmin(*) + sep * ( i + .5 )
  END

tot = total( hists(*,0) ) + nlo(0) + nhi(0)
;
;  close the file
close, unit
;
;and end
return
end
