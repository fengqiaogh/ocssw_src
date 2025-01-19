  pro make_cross_coef_mat, band, clutfile, iagg, nib, gmat, ncpix, cmat, iret
  
; Program to read the OCI cross-correlation matrix for a band
; and aggregate the coefficients for the instrument configuration

;  Name		Type 	I/O	Description
;
;  Band	string	 I	Band ID (blue, red or SWIR)
;  clutfile	string	 I	Crosstalk influence coeffient LUT
;  iagg	int	 I	Spatial aggregation factor
;  nib		int	 I	Number of bands in data
;  gmat[*,*]	float	 I	Gain aggregation matrix
;  ncpix	int	 O	Number of influence pixels
;  cmat[*,*,*]	float	 O	Crosstalk coefficient matrix	
;  iret	int	 O	Return code (-1 = failure)

; Check for valid band 
  if (band ne 'blue' and band ne 'red' and band ne 'SWIR') then begin
    print,'Invalid band specification'
    iret = -1
    return
  endif  

; Open file and read coefficients and dimensions
  bandstr = strjoin([band, '_cross_coef'])
  fid = ncdf_open(clutfile)
  if (fid eq -1) then begin
    print,'Invalid crosstalk LUT file name'
    iret = -1
    return
  endif
 
  ncdf_varget, fid, bandstr, cmat_in
  pdimid = ncdf_dimid(fid, 'number_of_pixels')
  ncdf_diminq, fid, pdimid, dimnm, ncpix
  pdimid = ncdf_dimid(fid, 'number_of_bands')
  ncdf_diminq, fid, pdimid, dimnm, nbands
  ncdf_close,fid
  
; If SWIR bands, use matrix as is
  if (band eq 'SWIR') then begin
    cmat = cmat_in
    return
  endif
  
; Generate correction matrix for aggregated bands
;  This requires a unitized version of the gain aggregation matrix
;  gmat is sized for 512 bands, the crosstalk coefficients are for 480 bands.
  
  If (band eq 'blue') then gmat_tmp = gmat[*,(512-nbands):*]
  If (band eq 'red') then gmat_tmp = gmat[*,0:nbands-1]
  gmatu = gmat_tmp
  k = where(gmatu ne 0)
  gmatu[k] = 1.0
  cmat_tmp = fltarr(nib, nib, ncpix)
  for i=0,ncpix-1 do cmat_tmp[*,*,i] = gmat_tmp#cmat_in[*,*,i]#transpose(gmatu)
  
; If spatial aggregation is not 8
  if (iagg ne 8) then begin
    iagf = 8/iagg
    ncpixc = ncpix*iagf
    cmat = fltarr(nib, nib, ncpixc)
    indf = indgen(ncpix)*iagf
    for i=0,iagf-1 do cmat[*,*,indf+i] = cmat_tmp/iagf
    ncpix = ncpixc
  
  endif else cmat = cmat_tmp
;  stop  
  return
  end
 
  
