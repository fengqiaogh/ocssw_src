  pro oci_geonav,pos,vel,smat,coef,sunr,pview,npix,delt,xlat,xlon,solz,sola, $
	senz,sena,range,qflag,eclipse,rs,rmoon
;
;  This subroutine performs navigation of a scanning sensor on the 
;  surface of an ellipsoid based on an input orbit position vector and 
;  spacecraft orientation matrix.  It uses a closed-form algorithm for 
;  determining points on the ellipsoidal surface which involves 
;  determining the intersection of the scan plan with the ellipsoid.  
;  The sensor view vectors in the sensor frame are passed in as a 3xN array.
;
;  The reference ellipsoid is set according to the scan 
;  intersection coefficients in the calling sequence; an equatorial 
;  radius of 6378.137 km. and a flattening factor of 1/298.257 are 
;  used by both the Geodetic Reference System (GRS) 1980 and the 
;  World Geodetic System (WGS) 1984.
;
;  It then computes geometric parameters using the pixel locations on
;  the Earth, the spaecraft position vector and the unit Sun vector in
;  the geocentric rotating reference frame.  The outputs are arrays of
;  geodetic latitude and longitude, solar zenith and azimuth and sensor
;  zenith and azimuth.  The azimuth angles are measured from local
;  North toward East.  Flag values of 999. are returned for any pixels
;  whose scan angle is past the Earth's horizon.
;
;  Reference: "Exact closed-form geolocation geolocation algorithm for
;  Earth survey sensors", F. S. Patt and W. W. Gregg, IJRS, Vol. 15
;  No. 18, 1994.

;  Calling Arguments
;
;  Name		Type 	I/O	Description
;
;  pos(3)	R*4	 I	ECR Orbit Position Vector (km) at scan
;                                 mid-time
;  vel(3)       R*4      I      ECR Orbit Velocity Vector (km/sec)
;  smat(3,3)	R*4	 I	Sensor Orientation Matrix
;  coef(10)	R*4	 I	Scan path coefficients
;  sunr(3)	R*4	 I	Sun unit vector in ECR frame
;  pview(3,*)   R*4      I      Array of sensor view vectors
;  npix         R*4      I      Number of pixels to geolocate
;  delt(*)	R*8	 I	Pixel time offsets from scan mid-time		
;  xlat(*)	R*4	 O	Pixel geodetic latitudes
;  xlon(*)	R*4	 O	Pixel geodetic longitudes
;  solz(*)	R*4	 O	Pixel solar zenith angles
;  sola(*)	R*4	 O	Pixel solar azimuth angles
;  senz(*)	R*4	 O	Pixel sensor zenith angles
;  sena(*)	R*4	 O	Pixel sensor azimuth angles
;  range(*)     R*4      O      Sensor-to-pixel ranges
;  qflag(*)	byte	 O	Quality flag (1 = off Earth, 2 = solar eclipse)
;  eclipse	I	 I	Check for solar eclipse (1 = yes)
;  rs		R*4	 I	Earth-Sun distance (AU)		
;  rmoon(3)	R*4	 I	Earth-Moon vector (km) in ECR frame 	
;
;	Subprograms Called:
;
;	CROSSP		Compute cross product of two vectors
;
;
;	Program written by:	Frederick S. Patt
;				General Sciences Corporation
;				October 20, 1992
;
;	Modification History:
;
;       Created universal version of the SeaWiFS geolocation algorithm
;       by specifying view vectors as an input.  F. S. Patt, SAIC, 11/17/09

    
    xlat=fltarr(npix)& xlon=fltarr(npix)
    solz=fltarr(npix)& sola=fltarr(npix)
    senz=fltarr(npix)& sena=fltarr(npix)
    range=fltarr(npix)& qflag = bytarr(npix)
    geovec=fltarr(3,npix)
    vv=fltarr(3,npix)
    x1=fltarr(3,npix)
    ea=fltarr(3,npix)& no=fltarr(3,npix)& up=fltarr(3,npix)
    rh=fltarr(3,npix)& sl=fltarr(3,npix)& rl=fltarr(3,npix)

; Earth ellipsoid parameters 
    f = 1/298.257
    omf2 = (1.d0-f)*(1.d0-f)

; Move scan view vectors to local array
    vv = pview[*,0:npix-1]

; Compute sensor-to-surface vectors for all scan angles
;  Compute terms for quadradic equation
    o = coef[0]*vv[0,*]*vv[0,*] + coef[1]*vv[1,*]*vv[1,*] $
      + coef[2]*vv[2,*]*vv[2,*] + coef[3]*vv[0,*]*vv[1,*] $
      + coef[4]*vv[0,*]*vv[2,*] + coef[5]*vv[1,*]*vv[2,*]
    
    p = coef[6]*vv[0,*] + coef[7]*vv[1,*] + coef[8]*vv[2,*]
    q = coef[9]
    r = p*p-4.d0*q*o

;  Check for scan past edge of Earth

    xlat[*] = -32767.
    xlon[*] = -32767.
    range[*] = -32767.
    qflag[*] = 1

    n1 = where(r ge 0.)
    if (n1[0] eq -1) then return
;  Solve for magnitude of sensor-to-pixel vector and compute components
    d = (-p[n1]-sqrt(r[n1]))/(2.d0*o[n1])
;  Check for negative distance
    n = where(d gt 0)
    if (n[0] eq -1) then return    
    n1 = n1[n]
    d = d[n]
    x1[0,n1]=d*vv[0,n1]
    x1[1,n1]=d*vv[1,n1]
    x1[2,n1]=d*vv[2,n1]

;  Convert velocity vector to ground speed
    re = 6378.137
    omegae = 7.29211585494d-5
    pm = sqrt(pos[0]^2+pos[1]^2 + pos[2]^2)
    clatg = sqrt(pos[0]^2+pos[1]^2)/pm
    rg = re*(1.-f)/sqrt(1.-(2.-f)*f*clatg^2)
    v = vel
;    v[0] = v[0] - pos[1]*omegae
;    v[1] = v[1] + pos[0]*omegae
    v = v*rg/pm

;  Transform vector from sensor to geocentric frame
    rh = transpose(smat)#x1
    for k=0,2 do geovec[k,n1] =  pos[k] + rh[k,n1] + v[k]*delt[n1]
    
;    Compute the local vertical, East and North unit vectors  
    uxy = geovec[0,*]*geovec[0,*]+geovec[1,*]*geovec[1,*]
    temp = sqrt(geovec[2,*]*geovec[2,*] + omf2*omf2*uxy)
    up[0,*] = omf2*geovec[0,*]/temp
    up[1,*] = omf2*geovec[1,*]/temp
    up[2,*] = geovec[2,*]/temp
    upxy = sqrt(up[0,*]*up[0,*]+up[1,*]*up[1,*])
    ea[0,*] = -up[1,*]/upxy
    ea[1,*] = up[0,*]/upxy
    ea[2,*] = 0.0
;	    no=crossp[up,ea]
    no[0,*] = -up[2,*]*ea[1,*]
    no[1,*] = up[2,*]*ea[0,*]
    no[2,*] = up[0,*]*ea[1,*] - up[1,*]*ea[0,*]

;    Compute geodetic latitude and longitude
    xlat[n1] = !radeg*asin(up[2,n1])
    xlon[n1] = !radeg*atan(up[1,n1],up[0,n1])
    range[n1] = d
    qflag[n1] = 0

;     Transform the pixel-to-spacecraft and Sun vectors into the local
;	frame
    rl[0,*] = -ea[0,*]*rh[0,*] - ea[1,*]*rh[1,*] - ea[2,*]*rh[2,*]
    rl[1,*] = -no[0,*]*rh[0,*] - no[1,*]*rh[1,*] - no[2,*]*rh[2,*]
    rl[2,*] = -up[0,*]*rh[0,*] - up[1,*]*rh[1,*] - up[2,*]*rh[2,*]
    sl[0,*] = sunr#ea
    sl[1,*] = sunr#no
    sl[2,*] = sunr#up

;    Compute the solar zenith and azimuth
    solz[n1]=!radeg*atan(sqrt(sl[0,n1]*sl[0,n1]+sl[1,n1]*sl[1,n1]),sl[2,n1])
;      Check for zenith close to zero
    n2 = where(solz gt 0.01)
    sola[n2] = !radeg*atan(sl[0,n2],sl[1,n2])

;    Compute the sensor zenith and azimuth
    senz[n1]=!radeg*atan(sqrt(rl[0,n1]*rl[0,n1]+rl[1,n1]*rl[1,n1]),rl[2,n1])
;      Check for zenith close to zero
    n2 = where(senz gt 0.01)
    sena[n2] = !radeg*atan(rl[0,n2],rl[1,n2])
    
; If eclipse check required
    if (eclipse) then begin

; Get Sun angular radius
      radsun = 0.53313/rs/!radeg/2 ; Reference Wertz, App. l, Table L-4
      
; Get vectors from Earth locations to Moon
      geo2moon = fltarr(3,npix)
      for k=0,2 do geo2moon[k,n1] = rmoon[k] - geovec[k,n1]
      dmoon = sqrt(geo2moon[0,n1]^2 + geo2moon[1,n1]^2 + geo2moon[2,n1]^2)
; Compute angle between location-Moon and Sun vectors
      smdiff = fltarr(3,npix)
      for k=0,2 do smdiff[k,n1] = geo2moon[k,n1]/dmoon - sunr[k]
      smdiffm = sqrt(smdiff[0,n1]^2 + smdiff[1,n1]^2 + smdiff[2,n1]^2)
      
; Get Moon angular radius
      radmoon = 1738.2/dmoon
; Set eclipse flag where angle is less than combined Sun and Moon radii
      ecl = where(smdiffm lt (radsun+radmoon))
      if (ecl[0] ne -1) then qflag[n1[ecl]] = qflag[n1[ecl]] + 2

    endif    
      
  return
  end
