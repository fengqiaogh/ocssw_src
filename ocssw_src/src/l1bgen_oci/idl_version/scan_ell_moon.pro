	pro scan_ell_moon,pos,smat,coef

;  This program calcualtes the coefficients which
;  represent the Moon scan track in the sensor frame.
;
;  The reference ellipsoid uses an equatorial radius of 1738.1 km and
;  a flattening factor of 1/827.667. 

;  Calling Arguments
;
;  Name		Type 	I/O	Description
;
;  pos(3)	R*4	 I	ECR Orbit Position Vector (km)
;  smat(3,3)	R*4	 I	Sensor Orientation Matrix
;  coef(10)	R*4	 O	Scan path coefficients

 	coef=dblarr(10)
	r = 1738.1d0
	f = 1./827.667d0
	omf2 = (1.d0-f)*(1.d0-f)

;  Compute constants for navigation model using Earth radius values
	rd=1.d0/omf2
	
;  Compute coefficients of intersection ellipse in scan plane
	sm = smat
	p = pos
	coef(0) = 1.d0 + (rd-1.d0)*sm(0,2)*sm(0,2)
	coef(1) = 1.d0 + (rd-1.d0)*sm(1,2)*sm(1,2)
	coef(2) = 1.d0 + (rd-1.d0)*sm(2,2)*sm(2,2)
	coef(3) = (rd - 1.d0)*sm(0,2)*sm(1,2)*2.d0
	coef(4) = (rd - 1.d0)*sm(0,2)*sm(2,2)*2.d0
	coef(5) = (rd - 1.d0)*sm(1,2)*sm(2,2)*2.d0
	coef(6) = (sm(0,0)*p(0) + sm(0,1)*p(1) + sm(0,2)*p(2)*rd)*2.d0
	coef(7) = (sm(1,0)*p(0) + sm(1,1)*p(1) + sm(1,2)*p(2)*rd)*2.d0
	coef(8) = (sm(2,0)*p(0) + sm(2,1)*p(1) + sm(2,2)*p(2)*rd)*2.d0
	coef(9) = p(0)*p(0) + p(1)*p(1) + p(2)*p(2)*rd - r*r

	return
	end

