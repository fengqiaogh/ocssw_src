	pro mat2rpy,pos,vel,smat,rpy

;  This program calculates the attitude angles from the ECEF orbit vector and 
;  attitude matrix.  The rotation order is (1,2,3).
;
;  The reference ellipsoid uses an equatorial radius of 6378.137 km and
;  a flattening factor of 1/298.257 (WGS 1984). 

;  Calling Arguments
;
;  Name		Type 	I/O	Description
; 
;  pos(3)	R*4	 I	Orbit position vector (ECEF)
;  vel(3)	R*4	 I 	Orbit velocity vector (ECEF)
;  smat(3,3)	R*4	 I	Sensor attitude matrix (ECEF to sensor)
;  rpy(3)	R*4	 O	Attitude angles (roll, pitch, yaw)

  re = 6378.137d0&rem = 6371.d0
  f = 1./298.257d0
  omegae = 7.29211585494d-5
  omf2 = (1.d0-f)*(1.d0-f)

;  Compute constants for navigation model using Earth radius values
  rd=1.d0/omf2

;  Determine local vertical reference axes
  v = vel
  p = pos
  v(0) = v(0) - p(1)*omegae
  v(1) = v(1) + p(0)*omegae

;  Compute Z axis as local nadir vector
  pm = sqrt(p(0)*p(0)+p(1)*p(1)+p(2)*p(2))
  omf2p = (omf2*rem + pm - rem)/pm
  pxy = p(0)*p(0)+p(1)*p(1)
  temp = sqrt(p(2)*p(2) + omf2p*omf2p*pxy)
  z = fltarr(3)
  z(0) = -omf2p*p(0)/temp
  z(1) = -omf2p*p(1)/temp
  z(2) = -p(2)/temp

;  Compute Y axis along negative orbit normal
  on=crossp(v,z)
  onm = sqrt(on(0)*on(0)+on(1)*on(1)+on(2)*on(2))
  y = -on/onm

;  Compute X axis to complete orthonormal triad (velocity direction)
  x = crossp(y,z)

;  Store local vertical reference vectors in matrix 
  om = fltarr(3,3)
  om(0,*)=x
  om(1,*)=y
  om(2,*)=z

;  Compute orbital-to-spacecraft matrix
  rm = smat#invert(om)

;  Compute attitude angles
  rpy = fltarr(3)
  rpy(0) = !radeg*atan(-rm(2,1),rm(2,2))
;  cosp = sqrt(rm(2,1)^2+rm(2,2)^2)
;  if (rm(2,2) lt 0) then cosp = -cosp
;  rpy(1) = !radeg*atan(rm(2,0),cosp)
  rpy(1) = !radeg*asin(rm(2,0))
  rpy(2) = !radeg*atan(-rm(1,0),rm(0,0))
  
  return
  end

