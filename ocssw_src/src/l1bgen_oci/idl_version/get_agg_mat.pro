  pro get_agg_mat, iagg, jagg, noagg, nib, nbb, amat, gmat

; Program to generate matrices to aggregate gain corrections for the input band aggregations and
;  to aggregate instrument bands to L1B bands

;  Name	Type 	I/O	Description
;
;  iagg	 I	 I	Spatial aggregation factor
;  jagg(16)	 I	 I	Aggregation factor/enable flag for each tap
;  noagg	 L	 I	L1A/L1B aggregation flag (1 = no aggregation)
;  nib		 I	 O	Number of instrument bands based on aggregation factors
;  nbb		 I	 O	Number of L1B bands after software aggregation
;  amat(nbb,nib) R*4	 O	Matrix to aggregate instrument to L1B bands
;  gmat(nib,512) R*4	 O	Matrix to aggregate gains from CCD to instrument output

; Determine number of instrument bands
;  Check for active taps
  ia = where(jagg gt 0)
  if (ia[0] eq -1) then begin
    print,'All taps disabled'
    nib = 1
    nbb = 1
    return
  endif
  ntb = intarr(16)
  ntb[ia] = 32/jagg[ia]
  nib = fix(total(ntb))
  gmat = fltarr(nib, 512)

; Populate gain aggregation matrix
  ii = 0
  itt = intarr(2,16) ; tap start/end locations for instrument bands
  for i=0,15 do begin
    if (jagg[i] gt 0) then begin
;      iaf = min([iagg*jagg[i],4]); Gain correction for bit-shifting 
      for k=0, ntb[i]-1 do begin
        ic = 32*i
	kj = k*jagg[i]
	gmat[ii+k,ic+kj:ic+kj+jagg[i]-1] = 1.0/jagg[i];*4/iaf; assumes LUT gains are 8x1
      endfor
      itt[0,i] = ii
      ii = ii + ntb[i]
      itt[1,i] = ii-1
    endif
  endfor
  
; If no spectral aggregation, identity matrix
  if (noagg) then begin
    nbb = nib
    amat = fltarr(nib,nib)
    for j=0,nib-1 do amat[j,j] = 1.0
    return
  endif else begin 

; Compute number of bands for 8x aggregation with overlapping bands
  nia = n_elements(ia)
;  nbb = ntb[0]*3/4 + 1
  nbb = ntb[ia[0]]*3/4 + 1
;  for i=1,15 do if (jagg[i] ge jagg[i-1]) then nbb = nbb + ntb[i] else $
;	nbb = nbb + ntb[i]*3/4 + ntb[i-1]/4
  for i=1,nia-1 do if (jagg[ia[i]] ge jagg[ia[i]-1]) then nbb = nbb + ntb[ia[i]] else $
	nbb = nbb + ntb[ia[i]]*3/4 + ntb[ia[i]-1]/4
  amat = fltarr(nbb,nib)
  
; First tap

  for k=0,ntb[ia[0]]*3/4 do amat(k,k:k+ntb[ia[0]]/4-1) = jagg[ia[0]]/8.0
  ib = k
;  ic = k

; Remaining taps
  for i=ia[1],15 do begin
    if (ntb[i] gt 0) then begin
      if (ntb[i] ge ntb[i-1]) then begin ; Transition resolution determined by preceding tap
      	nr = ntb[i-1]/4 - 1; Remaining bands using preceding tap
      	if (nr gt 0) then begin
          for k=0,nr-1 do begin
	    k1 = nr - k - 1
	    k2 = (k+1)*ntb[i]/ntb[i-1] - 1
	    amat[ib+k,itt[1,i-1]-k1:itt[1,i-1]] = jagg[i-1]/8.0
	    amat[ib+k,itt[0,i]:itt[0,i]+k2] = jagg[i]/8.0
      	  endfor
      	  ib = ib + nr
;	ic = ic + nr
      	endif
      endif else begin ; Transition resolution determined using current tap
      	nr = ntb[i]/4 - 1; Remaining bands using previous tap
      	if (nr gt 0) then begin
      	  for k=0,nr-1 do begin
	    k1 = (nr-k)*ntb[i-1]/ntb[i] - 1
	    k2 = k
	    amat[ib+k,itt[1,i-1]-k1:itt[1,i-1]] = jagg[i-1]/8.0
	    amat[ib+k,itt[0,i]:itt[0,i]+k2] = jagg[i]/8.0
      	  endfor
   	  ib = ib + nr
      	endif
      endelse

; Remaining bands using this tap
      for k=0,ntb[i]*3/4 do amat[ib+k,itt[0,i]+k:itt[0,i]+k+ntb[i]/4-1] = jagg[i]/8.0
      ib = ib + k
    endif
  endfor
  endelse

  return
  end    	

