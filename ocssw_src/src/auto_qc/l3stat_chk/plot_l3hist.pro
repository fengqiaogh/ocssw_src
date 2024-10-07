PRO plot_l3hist, file, mode, hists, xaxes
;
;  plot_l3hist will read the binary histogram file and
;  get the histograms of the 12 parameters and plot them
;
;  file  input file to read
;  mode  0 to just return the histogram arrays 1 to plot also
;  hists  12 size 100 histogram arrays
;  xaxes  12 x axes for the plot
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
hists = lonarr(100,12)
;
;  read it in
readu, unit, bins, ifile, hmin, hmax, nlo, nhi, hists
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

;
;  close the file
close, unit
;
;and end
return
end
