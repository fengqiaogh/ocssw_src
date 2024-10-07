C
C@***s* PROJECT_PFSST/l2gen_both/raygetpol.f
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
C      raygetpol.f
C
C LOCATION
C  $OCSSWROOT 
C
C 
C PURPOSE
C      get coefficients needed to compute Rayleigh radiance
C
C DESCRIPTION
C  
C      The subroutine RAYGET called by l2gen_both obtains the azimuthal coefficients for Rayleigh
C      radiance calculations. The coefficients are also obtained for various
C      solar and satellite zenith angles.  The coefficients obtained by this
C      routine are used in computing calibrated radiance for AVHRR channels 1 and 2.
C
C NOAA PFSST-SEADAS BUILD VERSION
C      Pathfinder SST V5.2 code built with SEADAS version 6.3 64 bit l2gen_both for  
C      CDR processed at University of Miami/RSMAS
C
C PRIMARY SEADAS CODE DOCUMENTATION NASA
C      For complete documentation of multi sensor SEADAS code see 
C      http://seadas.gsfc.nasa.gov/seadas/doc/toplevel/sds_program.html
C      http://seadas.gsfc.nasa.gov/seadas/doc/l2gen/l2gen.html
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
C  Fortran   
C
C*****
C##############################################################
C
C !F90
C
C !Description:
C    The subroutine RAYGET obtains the azimuthal coefficients for Rayleigh
C    radiance calculations. The coefficients are also obtained for various
C    solar and satellite zenith angles.  The coefficients obtained by this
C    routine are for avhrr channels 1 and 2.
C
C !Subroutines and Functions:
C   subroutine RAYGET(LUN,THETA,THETA0,PHI0,PHI1,PHI2,TAU,ISFLAG)
C   !Input Parameters:
C       LUN(integer) - logical unit for file access
C   !Output Parameters:
C       THETA(real array, size M_M) - Viewing angles
C       THETA0(real array, size N_N) - Sun angles
C       PHI0(real array, size M_M,N_N,C_H) - Azimuthal function
C       PHI1(real array, size M_M,N_N,C_H) - Azimuthal function
C       PHI2(real array, size M_M,N_N,C_H) - Azimuthal function
C       TAU(real array, size N_N,C_H) - Rayleigh optical depth
C       ISFLAG(integer array, size C_H) - channel present flag
C
C !Revision History:
C !Team_Unique Header:
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
C !End##############################################################

       subroutine rayget(lun,theta,theta0,phi0,phi1,phi2,tau,isflag)
       implicit none


C
       integer, parameter :: N_AMSIZ = 384

C Read Rayleigh multiple scattering polarization data

       integer, parameter :: M_M = 41
       integer, parameter :: N_N = 45
       integer, parameter :: C_H = 2

       integer*4 lenstr;
       integer*4 lun
       real*4    theta (m_m)                ! Viewing angles
       real*4    theta0(n_n)                ! Sun angles
       real*4    phi0  (m_m,n_n,c_h)        ! Azimuthal function
       real*4    phi1  (m_m,n_n,c_h)        !
       real*4    phi2  (m_m,n_n,c_h)        !
       real*4    tau   (n_n,c_h)            ! Rayleigh optical depth
       integer*4    isflag(c_h)             ! Channel present flag

       character*(80) msgbuf
       character caldir*255
       character*(n_amsiz) filenm
       character*7 files(c_h)
       integer*4    i, j, k

       integer*4    ii, jj
       integer*4    flen

10     format(//(1x,5f13.0))
20     format(//(1x,2f13.0))

C       call mergec (files(1), "/avhch1"//char(0))
C       call mergec (files(2), "/avhch2"//char(0))
       files(1) = '/avhch1';
       files(2) = '/avhch2';

       call getenv('OCDATAROOT', caldir)
       flen = lenstr(caldir)
       if (flen .eq. 0) then
           write(*,*)
     &       '-E- RAYGETPOL: Environment variable OCDATAROOT undefined'
           call exit(1)
       end if

       do k=1,c_h

C           call mergec (filenm, caldir)
C           call append (filenm, files(k))

C           call append (filenm, "pol.new"//char(0))    ! zero interpolated data

           filenm = caldir(1:flen)//'/avhrr/cal'//files(k)(1:7)//'pol.new'
           write(*,*) 'Loading '//filenm(1:lenstr(filenm))
           open(unit=lun,file=filenm,status='old',form='formatted',
     &         iostat=isflag(k))
           if (isflag(k) .eq. 0) then
           read(lun,10) theta
           do i=1,n_n

               read(lun,20) theta0(i),tau(i,k)
               read(lun,10) (phi0(j,i,k),j=1,m_m)

C# note: .new files have 0's for theta0 == 0
               read(lun,10) (phi1(j,i,k),j=1,m_m)
               read(lun,10) (phi2(j,i,k),j=1,m_m)
           end do

           close(unit=lun)
           cycle            ! found data for this channel
           end if

           write(*,*)
     &       '-E- RAYGETPOL: No data found for ', filenm
           call exit(1)
       end do

       end
