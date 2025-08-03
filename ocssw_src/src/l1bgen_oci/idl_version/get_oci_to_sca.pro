  pro get_oci_to_sca, spin, scaspin, scapos, oci_to_sca, face, delta

; IDL routine to compute OCI-to-SCA diffuser transformation matrix
;  SCA frame is +Z normal to face and +Y downward

;       Arguments
;
;       Name    	Type    I/O     Description
;       ----    	----    ---     -----------
;	spin		 I	 I	OCI spin ID
;	scaspin(*)	 I	 I	Array of SCA spin IDs
;	scapos(*)	 I	 I	Array of SCA position angles
;	oci_to_sca(3,3) R*4	 O	OCI-to-SCA transformation matrix
;	face		 I	 O	Diffuser face: 1 = daily bright, 2 = dim, 
;					 3 = monthly bright, -1 = no SCA data for spin
;	delta		R*4	 O	Diffuser position delta from nominal value

; Based on figure provided by Daniel Senai-Alemou

  posnom = [83.425, 203.425, 323.425]
  
; Determine SCA position to use; interpolate if necessary
  nsca = n_elements(scaspin)
  ins = where(scaspin eq spin)
  if (ins eq -1) then begin
    face = -1
    if (spin gt scaspin[0] AND spin lt scaspin[nsca-1]) then ins = max(where(scaspin lt spin))
    if (spin eq scaspin[0]-1) then ins = 0
    if (spin eq scaspin[nsca-1]+1) then ins=nsca-2
    if (ins eq -1) then return
    spos = ((spin - scaspin[ins])*scapos[ins+1] + (scaspin[ins+1] - spin)*scapos[ins]) / $
      (scaspin[ins+1] - scaspin[ins])
  endif else spos = scapos[ins]

; Get difuser face index and compute view angle
  face = fix(spos - 23.425)/120 + 1
  delta = spos - posnom[face-1]
  v_ang = 36.0 - delta
  
; Compute transformation matrix
  ox = [-cos(v_ang/!radeg), -sin(v_ang/!radeg), 0.]
  oy = [0., 0., 1]
  oz = [-sin(v_ang/!radeg), cos(v_ang/!radeg), 0.]
  oci_to_sca = fltarr(3,3)
  
  oci_to_sca[0,*] = ox
  oci_to_sca[1,*] = oy
  oci_to_sca[2,*] = oz
  
  return
  end
  
  
  
