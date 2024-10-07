C Incident ray
C       i = ( sin(PHI)       , 0              , cos(PHI) )
C
C Reflected ray
C       o = ( sin(MU)cos(NU) , sin(MU)sin(NU) , cos(MU)  )
C
C Reflection angle
C       i . o = cos(2 OMEGA) = sin(PHI)sin(MU)cos(NU)+cos(PHI)cos(MU)
C
C Normal to plane of reflection
C       n = ( -sin(BETA)cos(ALPHA) , -sin(BETA)sin(ALPHA) , cos(BETA) )
C
C Reflection plane orientation
C       n = i+o
C
C       n = ( sin(PHI)+sin(MU)cos(NU) , sin(MU)sin(NU) , cos(PHI)+cos(MU) )/
C            (2 cos(OMEGA))
C
C  Tilt of plane
C       cos(BETA) = (cos(PHI)+cos(MU))/(2 cos(OMEGA))
C
C  Rotation of plane
C       (i+n)^2 = 4 cos(OMEGA/2)
C
C       cos(OMEGA) = cos(BETA)cos(PHI)-sin(BETA)sin(PHI)cos(ALPHA)
C
C   See papers by Cox and Munk (1954a,b; 1956) for derivation of equations.


        subroutine getglint(X1,X2,X3,X4,X5,X6)
C
C  Calculate sun gliter coefficient
C
C   X1  Angle MU       (sensor zenith angle)
C   X2  Angle PHI      (solar zenith angle)
C   X3  Angle NU       (sensor-sun azimuth)
C   X4  Wind speed     (m/s)
C   X5  Wind direction (radians)
C   X6  Radiance       (ignoring the atmosphere)
C
C #include "sensor_cmn.fin"
        save
        real pi
        parameter (pi=3.1415926)
        real X2, X3, X4, X5, acoss, RHO
        real X1, Y1, RAD, X, Y2, OMEGA, BETA, ALPHA, SIGC, SIGU, Y3
        real CHI, ALPHAP, SWIG, ETA, PROB, X6, EXPON
        real Y4

        RAD(X) = X*pi/180.

        Y4 = max(X4,0.001)

        Y1 = RAD(X1)
        if (Y1.eq.0.) Y1 = 1.e-7
        Y2 = RAD(X2)
        if (Y2.eq.0.) Y2 = 1.e-7
        Y3 = RAD(X3)
        OMEGA  = ACOSS(COS(Y1)*COS(Y2)-SIN(Y1)*SIN(Y2)*COS(Y3))/2.
        if (OMEGA.eq.0.) OMEGA = 1.e-7
        BETA   = ACOSS((COS(Y1)+COS(Y2))/(2.*COS(OMEGA)))
        if (BETA.eq.0.) BETA = 1.e-7
        ALPHA  = ACOSS((COS(BETA)*COS(Y2)-COS(OMEGA))/(SIN(BETA)*SIN(Y2)))
        if (SIN(Y3).lt.0.) ALPHA = -ALPHA

c Isotropic wind
c        if (glintOn .ge. 2) then
c           ! from Ebuchi & Kizu
c            SIGC = SQRT(0.0048 + 0.00152*Y4)
c            SIGU = SQRT(0.0053 + 0.000671*Y4)
c        else
c           ! from Cox & Munk
            SIGC   = .04964*SQRT(Y4)
            SIGU   = .04964*SQRT(Y4)
c        endif
        CHI    = X5
        ALPHAP = ALPHA+CHI
        SWIG   = SIN(ALPHAP)*TAN(BETA)/SIGC
        ETA    = COS(ALPHAP)*TAN(BETA)/SIGU
        EXPON  = -(SWIG**2+ETA**2)/2. 
        if (EXPON.lt.-30.) EXPON = -30.         ! trap underflow
        if (EXPON.gt.+30.) EXPON = +30.         ! trap overflow
        PROB   = EXP(EXPON)/(2.*pi*SIGU*SIGC)
        call REFLEC(OMEGA,RHO)

c Normal distribution
        X6     = RHO*PROB/(4.*COS(Y1)*COS(BETA)**4)

        return
        end


        subroutine getglint_iqu(X1,X2,X3,X4,X5,X6,X7,X8)

C  Calculate sun gliter coefficient, including polarization
C
C       X1(real) - Angle MU       (sensor zenith angle)
C       X2(real) - Angle PHI      (solar zenith angle)
C       X3(real) - Angle NU       (sensor-sun azimuth)
C       X4(real) - Wind speed     (m/s)
C       X5(real) - Wind direction
C       X6(real) - Sun glitter coefficient
C       X7(real) - Q/I for glitter
C       X8(real) - U/I for glitter
C
C #include "sensor_cmn.fin"
        save
        real pi
        parameter (pi=3.1415926)
        real X2, X3, X4, X5, acoss, RHO_PLUS
        real X1, Y1, RAD, X, Y2, OMEGA, BETA, ALPHA, SIGC, SIGU
        real CHI, ALPHAP, SWIG, ETA, PROB, X6, EXPON
        real X7, X8, RHO_MINUS, SR, CR, ROT_ANG, C2R, S2R, ASINN
        real Y4

        RAD(X) = X*pi/180.

        Y4 = max(X4,0.001)

        Y1 = RAD(X1)
        if (Y1.eq.0.) Y1 = 1.e-7
        Y2 = RAD(X2)
        if (Y2.eq.0.) Y2 = 1.e-7
        Y3 = RAD(X3)
        OMEGA  = ACOSS(COS(Y1)*COS(Y2)-SIN(Y1)*SIN(Y2)*COS(Y3))/2.
        if (OMEGA.eq.0.) OMEGA = 1.e-7
        BETA   = ACOSS((COS(Y1)+COS(Y2))/(2.*COS(OMEGA)))
        if (BETA.eq.0.) BETA = 1.e-7
        ALPHA  = ACOSS((COS(BETA)*COS(Y2)-COS(OMEGA))/(SIN(BETA)*SIN(Y2)))
        if (SIN(Y3).lt.0.) ALPHA = -ALPHA

c Isotropic wind
c        if (glintOn .ge. 2) then
c           ! from Ebuchi & Kizu
c            SIGC = SQRT(0.0048 + 0.00152*Y4)
c            SIGU = SQRT(0.0053 + 0.000671*Y4)
c        else
c           ! from Cox & Munk
            SIGC   = .04964*SQRT(Y4)
            SIGU   = .04964*SQRT(Y4)
c        endif
        CHI    = X5
        ALPHAP = ALPHA+CHI
        SWIG   = SIN(ALPHAP)*TAN(BETA)/SIGC
        ETA    = COS(ALPHAP)*TAN(BETA)/SIGU
        EXPON  = -(SWIG**2+ETA**2)/2.
        if (EXPON.lt.-30.) EXPON = -30.         ! trap underflow
        if (EXPON.gt.+30.) EXPON = +30.         ! trap overflow
        PROB   = EXP(EXPON)/(2.*pi*SIGU*SIGC)
        call REFLEC_BOTH(OMEGA,RHO_PLUS,RHO_MINUS)

c Normal distribution
        X6     = RHO_PLUS*PROB/(4.*COS(Y1)*COS(BETA)**4)

c Polarization components
        if (OMEGA .gt. .0001) then
          CR = (cos(Y2) - cos(2.*OMEGA)*cos(Y1))/(sin(2.*OMEGA)*sin(Y2))
          SR = sin(Y2)*sin(pi-Y3) / sin(2.*OMEGA)
          ROT_ANG = sign(1., CR)*asinn(SR)
        else
          ROT_ANG = pi/2.
        endif

        C2R = cos(2.*ROT_ANG)
        S2R = sin(2.*ROT_ANG)

        X7 =  C2R * RHO_MINUS / RHO_PLUS ! q_ov_i
        X8 = -S2R * RHO_MINUS / RHO_PLUS ! u_ov_i

        return
        end


        subroutine REFLEC (X1,X3)
C
C  X1  Incident angle (radians)
C  X3  Reflectance
C
C       n1 sin(x1) = n2 sin(x2)
C
C                    tan(x1-x2)**2
C       Refl(par ) = -------------
C                    tan(x1+x2)**2
C
C                    sin(x1-x2)**2
C       Refl(perp) = -------------
C                    sin(x1+x2)**2
C
C  Where:
C       x1  Incident angle
C       n1  Index refraction of Air
C       x2  Refracted angle
C       n2  Index refraction of Water
C
      real X1, X2, X3, ref
*                                       ! Index refraction of sea water
        REF = 4./3.
        if (X1.lt..00001) then
          X3 = .0204078
        else
          X2 = ASIN(SIN(X1)/REF)
          X3 = (SIN(X1-X2)/SIN(X1+X2))**2+(TAN(X1-X2)/TAN(X1+X2))**2
          X3 = X3/2.
        end if
        return
        end


        subroutine REFLEC_BOTH (X1,X3,X4)
C
C  X1  Incident angle (radians)
C  X3  Reflectance sum
C  X4  Reflectance difference
C
C       n1 sin(x1) = n2 sin(x2)
C
C                    tan(x1-x2)**2
C       Refl(par ) = -------------
C                    tan(x1+x2)**2
C
C                    sin(x1-x2)**2
C       Refl(perp) = -------------
C                    sin(x1+x2)**2
C
C  Where:
C       x1  Incident angle
C       n1  Index refraction of Air
C       x2  Refracted angle
C       n2  Index refraction of Water
C
      real X1, X2, X3, X4, REF
      real PERP, PAR
c                                       ! Index refraction of sea water
        REF = 4./3.
        if (X1.lt..00001) then
          X3 = .0204078
          X4 = 0.
        else
          X2 = ASIN(SIN(X1)/REF)
          PERP = (SIN(X1-X2)/SIN(X1+X2))**2
          PAR  = (TAN(X1-X2)/TAN(X1+X2))**2
          X3 =  PERP+PAR
          X3 = X3/2.
          X4 = -PERP+PAR
          X4 = X4/2.
        end if
        return
        end


        function ACOSS (X1)
        real    ACOSS, X1
C
C   This calculates ACOS(X) when X is near + or - 1
C
        if (X1.ge.1.) then
          ACOSS = 0.
        else if (X1.le.-1.) then
          ACOSS = 3.141592654
        else
          ACOSS = acos(X1)
        end if
        return
        end


        function ASINN (X1)
        real    ASINN, X1
C
C   This calculates ASIN(X) when X is near + or - 1
C
        if (X1.ge.1.) then
          ASINN = 3.141592654/2.0
        else if (X1.le.-1.) then
          ASINN =-3.141592654/2.0
        else
          ASINN = asin(X1)
        end if
        return
        end
