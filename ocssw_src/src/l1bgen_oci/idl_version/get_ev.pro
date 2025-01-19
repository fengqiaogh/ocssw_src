  pro get_ev, secpline, dtype, lines, iagg, pcdim, psdim, ev_toff, clines, slines, deltc, delts, iret, DARK=dark

; Program to determine Earth view pixels and time offsets from PPR
;  Uses scan mode table fields to determine number of pixels offset

;  Calling Arguments
;
;  Name		Type 	I/O	Description
;
;  dtype(10)	I*2	 I	Data type index for each data zone (0=no data, 1=Earth view)
;  lines(10)    I*2 	 I 	Number of lines (unaggregated) in each zone
;  iagg(10)	I*2	 I	Spatial aggregation in each zone
;  pcdim		Int	 O	Number of hyperspectral science pixels
;  psdim		Int	 O	Number of SWIR science pixels
;  ev_toff	R*8	 O	Offset in seconds from PPR to EV mid-time
;  clines	R*4	 O	Line numbers corresponding to science pixel centers
;  slines	R*4	 O	Line numbers corresponding to SWIR pixel centers
;  deltc	R*8	 O	Time offset from EV mid-time for each science pixel
;  delts	R*8	 O	Time offset from EV mid-time for each SWIR pixel
;  iret		Int	 O	Return code (0 = normal, -1 = no Earth view)
;  DARK		Int	 I	Keyword to include dark collect

  ltype = [0,1,0,1,1,1,1,1,1,1,0,0,0]
  if (KEYWORD_SET(dark)) then ltype[2] = 1

; Find end of no-data zone
  iz = 0
  line0 = 0
  iret = -1
  while (dtype[iz] eq 0) do begin
    line0 = line0 + lines[iz]
    iz = iz + 1
  endwhile
  if (iz eq 10) then return

; Find number of pixels in Earth views
  clines = fltarr(32400)
  slines = fltarr(4050)
  pcdim = 0
  psdim = 0
  linen = line0
  for i=iz,9 do begin
; Check for not dark or no-data
;    if (dtype(i) eq 1 OR dtype(i) eq 3 OR dtype(i) eq 4 OR dtype(i) eq 6 OR dtype(i) eq 7 OR dtype(i) eq 9) then begin
;    if (dtype(i) ne 0 AND dtype(i) ne 2 AND dtype(i) lt 10) then begin
    if (ltype[dtype[i]]) then begin
      np = lines[i]/iagg[i]
      clines[pcdim:(pcdim+np-1)] = linen + findgen(np)*iagg[i] + (iagg[i])/2. - 64
      pcdim = pcdim + np
      ns = lines[i]/8
      slines(psdim:(psdim+ns-1)) = linen + findgen(ns)*8 + 4 - 64
      psdim = psdim + ns
      iret = 0
    endif
    linen = linen + lines[i]
  endfor
  if (iret eq -1) then return

; Calculate times
  clines = clines[0:pcdim-1]
  deltc = secpline*clines
  slines = slines[0:psdim-1]
  delts = secpline*slines
  ev_toff = (deltc[0] + deltc[pcdim-1])/2
  deltc = deltc - ev_toff
  delts = delts - ev_toff

  return
  end



