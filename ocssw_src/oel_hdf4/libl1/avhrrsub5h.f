C
C@***s* PROJECT_PFSST/l2gen_both/avhrrsub5.f
C
C This header contains documentation required by the NOAA Climate Data Record
C Program (CDRP), which is managed at the NOAA National Climatic Data Center (NCDC).
C Only the code that applies to AVHRR SST data is documented in this header. 
C
C The AVHRR Pathfinder Sea Surface Temperature (PFSST) processing code was originally 
C developed at the University of Miami.  In 2010, the code was integrated into
C the multi-sensor SeaWiFS Data Analysis System (SeaDAS) obtained from NASA GSFC.
C SeaDAS was used for processing the PFSST beginning with the Pathfinder Version
C 5.2 (PFV5.2) dataset, produced jointly by the University of Miami and the NOAA
C National Oceanographic Data Center (NODC).  These data are provided to the 
C public and archived at NODC, and have been transitioned along with the production 
C code and documentation to the CDRP at NCDC.
C 
C This NOAA required header is specifically written for Pathfinder SST and may
C not be relevant to other sensors or products processed by SeaDAS. Please 
C review the SEADAS software distribution policy for public domain software 
C located at http://seadas.gsfc.nasa.gov/copying.html for more information and 
C documentation
C
C NAME
C    avhrrsub5.f
C
C LOCATION
C $OCSSWROOT 
C 
C PURPOSE
C    Collection of subroutines needed to compute various atmospheric parameters
C    for AVHRR visible channels 1 and 2.
C
C DESCRIPTION
C    Functions and subroutines to Compute Rayleigh/atmospheric parameters for satellite pass.
C     This is the implementation for avhrr channels 1 and 2.
C
C    subroutine AVLOOPH (satz, solz, delphi, RAYLY, AERSOL, AGLINT) - The subroutine AVLOOPH calculates
C      and collects coefficients and parameters need for the atmospheric correction.
C
C    subroutine AVCONSH, is the constant initializer for the subroutine AVLOOPH.
C
C    function AVHRRSUB5H -wrapper for the subroutines AVCONSH, AVINITH, and AVLOOPH,for persistent
C      storage.  This routine is never called directly.
C
C
C NOAA PFSST-SEADAS BUILD VERSION
C Pathfinder SST V5.2 code built with SEADAS version 6.3 64 bit L2gen_both for  
C CDR processed at University of Miami/RSMAS
C
C PRIMARY SEADAS CODE DOCUMENTATION NASA
C For complete documentation of multi sensor SEADAS code see 
C http://seadas.gsfc.nasa.gov/seadas/doc/toplevel/sds_program.html
C
C AUTHOR
C                                                   
C                              
C   PFSST project embedded code 
C     Susan Walsh
C     University of Miami/RSMAS
C
C CREATION DATE 
C  2010  
C
C COPYRIGHT
C  THIS SOFTWARE AND ITS DOCUMENTATION ARE CONSIDERED TO BE IN THE PUBLIC DOMAIN AND
C  THUS ARE AVAILABLE FOR UNRESTRICTED PUBLIC USE.  THEY ARE FURNISHED "AS IS." THE 
C  AUTHORS, THE UNITED STATES GOVERNMENT, ITS INSTRUMENTALITIES, OFFICERS, EMPLOYEES, 
C  AND AGENTS MAKE NO WARRANTY, EXPRESS OR IMPLIED, AS TO THE USEFULNESS OF THE 
C  SOFTWARE AND DOCUMENTATION FOR ANY PURPOSE. THEY ASSUME NO RESPONSIBILITY (1) FOR
C  THE USE OF THE SOFTWARE AND DOCUMENTATION; OR (2) TO PROVIDE TECHNICAL SUPPORT 
C  TO USERS.
C
C MODIFICATION HISTORY
C 
C  See CVS revision control history embedded in actual file.
C               
C LANGUAGE
C Fortran   
C
C*****
C#############################################################
C
C !F90
C
C !Description:
C    Function to Compute Rayleigh/atmospheric parameters for satellite pass.
C    This is the implementation for avhrr channels 1 and 2.
C
C Subroutines and Functions:
C    integer function AVHRRSUB5H
C    !Description:
C          Wrapper for the subroutines AVCONSH, AVINITH, and AVLOOPH,for persistent
C       storage.  This routine is never called directly.
C    !Input Parameters:
C      none
C    !Output Parameters:
C       None
C
C    entry AVCONSH (LUNDK, inpixsiz, JDAY)
C    !Decsription:
C      The subroutine AVCONSH, is the constant initializer for the subroutine
C       AVLOOPH.  The subroutine AVCON begins by obtaining ozone optical
C      thickness with a call to getozone.  Clear water radiance values and
C      the solar constant are then corrected for time of year. The subroutine
C      rayget is then called to obtain the Rayleigh optical thickness, tau.
C      The subroutine AVCONSH then converts epsilon values to an internal form.
C    !Input Parameters:
C      LUNDK(integer) - Logical unit for file access
C      inpixsiz(integer) - number of pixels per line
C      JDAY(integer) - julian day of start of data
C    !Return value:
C      (integer) - status error messages from subroutines
C
C    entry AVLOOPH (satz, solz, delphi, RAYLY, AERSOL, AGLINT)
C    !Decsription:
C      The subroutine AVLOOPH calculates and collects coefficients and
C       parameters need for the atmospheric correction.
C       The subroutine begins by calculating solar angles and the "airmass"
C       parameters.  The routine then calculates the Rayleigh radiance,
C       <rayly>, given the present geometry.
C       The sun glitter coefficent, <aglint>, is obtained by a call to the
C       subroutine gliter.  The routine then calls the subroutine hmf8 to obtain
C       the aerosol radiance coefficient and calculates the aerosol structure
C       function, <aerosol>.
C    !Input Parameters:
C      satz(real array, size ANCHOR_1) - satellite zenith angle in degrees
C      solz(real array, size ANCHOR_1) - solar zenith angle in degrees
C      delphi(real array, size ANCHOR_1) - delta azimuth angle in degrees
C    !Output Parameters:
C      RAYLY(real array, size ANCHOR_1,2) - Rayleigh radiance by channel
C      AERSOL(real array, size ANCHOR_1) - Aerosol structure factor
C      AGLINT(real array, size ANCHOR_1) - Sun glitter coefficient
C
C !Revision History:
C

C !Team Unique Header:
C
C Copyright 1988-1998 by Rosenstiel School of Marine and Atmospheric Science,
C University of Miami, Miami, Florida.
C
C                       All Rights Reserved
C
C Permission to use, copy, modify, and distribute this software and its
C documentation for non-commercial purposes and without fee is hereby granted,
C provided that the above copyright notice appear in all copies and that both
C that copyright notice and this permission notice appear in supporting
C documentation, and that the names of University of Miami and/or RSMAS not be
C used in advertising or publicity pertaining to distribution of the software
C without specific, written prior permission.
C
C UNIVERSITY OF MIAMI DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
C INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS, IN NO EVENT
C SHALL UNIVERSITY OF MIAMI BE LIABLE FOR ANY SPECIAL, INDIRECT OR CONSEQUENTIAL
C DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS,
C WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING
C OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
C
C !References and Credits:
C    Written by:
C    University of Miami
C    Rosensteil School for Marine and Atmospheric Sciences
C    Division of Meteorology and Oceanography
C    4600 Rickenbacker Cswy
C    Miami,Fl
C    33149
C    Contact: SWalsh@rsmas.miami.edu
C
C !Design Notes:
C
C !End###############################################################

       function avhrrsub5h (jdummy)
       implicit none

C $COMSIZ - Maximum message length
C $NAMSIZ - Maximum filename length
C $ARRSIZ - Maximum pixels per line (lines/screen)

       integer, parameter :: C_OMSIZ = 512
       integer, parameter :: N_AMSIZ = 384
       integer, parameter :: A_RRSIZ = 8704

       integer, parameter :: ANCHOR_C = 1
       integer, parameter :: ANCHOR_ = 2048

       integer, parameter :: MULT_ = 1

       integer, parameter :: MX_BOXSIZ = 9

C      Local definitions.

       integer, parameter :: NUMGEO = 5

       real*8, parameter :: P_I = 3.1415926535898d0

C Read Rayleigh multiple scattering polarization data

       integer, parameter :: M_M = 41
       integer, parameter :: N_N = 45
       integer, parameter :: C_H = 2

       integer*4      jdummy;
       integer*4      lenstr;
       integer*4      isflag(c_h);
       real*4      theta (m_m);
       real*4      theta0(n_n);
       real*4      tau   (n_n,c_h);
       real*4      phi0  (m_m,n_n,c_h);
       real*4      phi1  (m_m,n_n,c_h);
       real*4      phi2  (m_m,n_n,c_h);
       real*4      phicof(c_h);

       integer*4      avhrrsub5h, avconsh, avlooph;

       logical      needray;
       integer*4      lundk;
       integer*4      inpixsiz;
       integer*4      ii;
       integer*4      nk, ianchr;
       integer*4      ni, nj, jnpixsiz;
       integer*4      jday;
       real*8      thesat;
       real*4      thesatd, thesund, fsun, fsat;
       real*4      phi00, phi01, phi10, phi11, phi21, phi23;
       real*4      thesun, sunmui, aercof;
       real*4      deg;
       real*4      rad;
       real*4      angle, c00;
       real*4      to(c_h)
       real*4      f1ln(c_h)
       real*4      f0var, x, fday, gday;
       real*8  x8;
       real*4      b00;

C      ANLY      CALEPS
C
C
C

       real*4      rayly (c_h);      ! Rayleigh radiance by channel
       real*4      aglint;
       real*4      aersol(c_h);
       
       real*4      satz
       real*4      solz
       real*4      delphi

       real*8      deg8;

C      System common blocks.

C      Variables to save across calls

       save      f1ln;
       save      to;
       save      needray;
       save      theta, theta0, phi0, phi1, phi2, tau, isflag;
       save      ianchr;
       save      jnpixsiz;

C      TO   - Ozone absorbtion

       data      to(1),  to(2)  / .090, .000 /;

       data      needray / .true. /;      ! Load tables first time

C
C ASF's
C
C  NOTE: The earth is closer to the sun in the winter.
C
C f0var: 7 % yearly variation
       f0var(x,fday)=x*(
     .      (1.0+.0167*sngl(dcos(2.0d0*p_i*(dble(fday)-3.0d0)/365.0d0)))**2 );
C       deg(x) = x*180./sngl(p_i)
       deg8(x8) = x8*180.0d0/p_i
       rad(x) = x*sngl(p_i)/180.

       stop 'AVHRRSUB5H is a package'

C
C Constant initializer
C

       entry avconsh (lundk, inpixsiz, jday)

C calibration information must be in cal1... for calaryin and
C                                    cal7... for calaryout
       jnpixsiz = inpixsiz;      ! Local copy
       ianchr = jnpixsiz;      ! do it for every pixel now

C  Get ozone optical thickness.

C  Correct solar constant and clear water radiance for time of year.

      gday = amod(float(jday),1000.);    ! Pick out julian day of year
      do ii=1,c_h
          f1ln(ii) = f0var(sngl(p_i)*100.,gday); ! Percent albedo
      end do

      if (needray) then

          !  Just load tables first time through.

          call rayget(lundk,theta,theta0,phi0,phi1,phi2,tau,isflag)
          do nk=1,c_h

              if (tau(1,nk) .eq. 0.) then
                  tau(1,nk) = 1.        ! Fake up missing channels
              end if

          end do

          needray = .false.
      end if

      avconsh=(1);
      return

C

      entry avlooph (satz, solz, delphi, rayly, aersol, aglint)

C      THESAT   Zenith angle to satellite from ground element
C      B00      Air-mass between satellite and viewed ground element

C      use value from class file instead of calculated value
      thesat = rad(satz)
      b00 = 1./sngl(dcos(thesat));      ! Air-mass to sat from ground

C      THESUN   Zenith angle to sun from ground element
C      SUNMUI   Air-mass between sun and viewed ground element

C      Geocentric latitude for SUNANG
      thesun = rad(solz)     ! airmass_avhrr needs solz as radians, not degrees
      call airmass_avhrr(thesun, sunmui)

C      C00      Total air-mass from sun to ground to satellite

      c00 = b00+sunmui;      ! Sun path + satellite path

C      ANGLE    Azimuth between view from satellite at ground and sun
C      RAYLY    Rayleigh radiance by channel

      angle = rad(delphi)

      if (sunmui .eq. 0.) then
          do nk=1,c_h
              phicof(nk) = 0.;        ! Out of range
          end do

      else
          thesatd = sngl(deg8(thesat));
          if (thesatd > theta(m_m)) then
              go to 901               ! Value out of range!
          end if

          ni=INT((((thesatd-theta(1))/(theta(m_m)-theta(1)))*float(m_m-1))+1.);
          do while (ni >= 1 .and. ni < m_m)
              if (theta(ni) <= thesatd .and.
     .                  theta(ni+1) >= thesatd) then
                  exit;

              else if (theta(ni+1) < thesatd) then
                  ni = ni+1;

              else  !if (THETA(NI) > THESATD)
                  ni = ni-1;
              end if

          end do

C          thesund = deg(thesun);
          thesund = solz;
          if (thesund > theta0(n_n)) then
              go to 901        ! Value out of range!
          end if

          nj=INT((((thesund-theta0(1))/(theta0(n_n)-theta0(1)))*float(n_n-1))+1.);
          do while (nj >= 1 .and. nj < n_n)
              if (theta0(nj) <= thesund .and.
     .            theta0(nj+1) >= thesund) then
                  exit;

              else if (theta0(nj+1) < thesund) then
                  nj = nj+1;

              else  !if (THETA0(NJ) > THESUND)
                  nj = nj-1;
              end if

          end do

          if (ni >= 1 .and. ni < m_m .and.
     .        nj >= 1 .and. nj < n_n) then
              fsat   = (thesatd-theta (ni))/(theta (ni+1)-theta (ni));
              fsun   = (thesund-theta0(nj))/(theta0(nj+1)-theta0(nj));
              do nk=1,c_h

                if (isflag(nk) .eq. 0) then
                  phi00  = phi0(ni  ,nj  ,nk)+
     .                          phi1(ni  ,nj  ,nk)*cos(angle)+
     .                          phi2(ni  ,nj  ,nk)*cos(2*angle);
                  phi10  = phi0(ni+1,nj  ,nk)+
     .                          phi1(ni+1,nj  ,nk)*cos(angle)+
     .                          phi2(ni+1,nj  ,nk)*cos(2*angle);
                  phi01  = phi0(ni  ,nj+1,nk)+
     .                          phi1(ni  ,nj+1,nk)*cos(angle)+
     .                          phi2(ni  ,nj+1,nk)*cos(2*angle);
                  phi11  = phi0(ni+1,nj+1,nk)+
     .                          phi1(ni+1,nj+1,nk)*cos(angle)+
     .                          phi2(ni+1,nj+1,nk)*cos(2*angle);
                  phi21  = phi00*(1-fsun)+phi01*fsun;
                  phi23  = phi10*(1-fsun)+phi11*fsun;
                  phicof(nk) = phi21*(1.-fsat)+phi23*fsat;

                else
                  phicof(nk) = 0.        ! No data for this channel
                end if

              end do

          else
              go to 901        ! Value out of range!
          end if

      end if

      go to 902

901      continue
      !  Handle value out of range
      do nk=1,c_h
          phicof(nk) = 0.;    ! Out of range
      end do

902      continue
      rayly(1) = phicof(1)*f1ln(1)*exp(-to(1)*c00);
      rayly(2) = phicof(2)*f1ln(2)*exp(-to(2)*c00);

      !  Compute an aerosol radiance -> tau conversion coefficient

      call hmf8(sngl(thesat),thesun,angle,1.,1.,aercof);
      aersol(1) = aercof*f1ln(1)*exp(-to(1)*c00);
      aersol(2) = aercof*f1ln(2)*exp(-to(2)*c00);

      !  Compute gliter coefficient
      !    Not sure why use getglint here instead of letting atmocor1 call getglint_iqu as it does for modis
      !    Only difference between the two is that getglint_iqu returns more variables.

C      call gliter(sngl(thesat),thesun,angle,6.,0.,aglint);
      call getglint(satz, solz, delphi, 6.0, 0.0, aglint);

      avlooph=(1)
      return

      end
