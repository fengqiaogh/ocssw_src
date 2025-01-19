  pro get_sc_to_oci, geo_lut, tilt, sc_to_oci

; Program to compute the spacecraft-to-OCI transformation matrix 
;  from the tilt angle, axis and OCI alignment matrix
;  Eventually need to add the capability to read this from files provided by MOC

;  Calling Arguments
;
;  Name		Type 	I/O	Description
;
;  geo_lut	struct   I	Geolocation LUT structure
;  tilt		R*4	 I	Tilt angle from S/C telemetry
;  sc_to_oci(3,3) R*4	 O	S/C to OCI transformation matrix

;  t_axis = [0.0, 1.0, 0.0] ; Assume perfect tilt axis alignment for now
;  t2oci = fltarr(3,3) ; Tilt mechanism to OCI alignment
;  for i=0,2 do t2oci(i,i) = 1.0 ; Assume perfect OCI alignment for now

; Model tilt rotation using a quaternion
  qt = fltarr(4)
  qt(0:2) = geo_lut.tilt_axis*sin(tilt/2/!radeg)
  qt(3) = cos(tilt/2/!radeg)
  qtom,qt,tiltm

; Combine tilt and alignments
;  sc2tiltp = tiltm#transpose(geo_lut.sc_to_tilt)
  sc2tiltp = transpose(geo_lut.sc_to_tilt)#tiltm
  sc2ocim = geo_lut.tilt_to_oci_mech#sc2tiltp
  sc_to_oci = geo_lut.oci_mech_to_oci_opt#sc2ocim

  return
  end

