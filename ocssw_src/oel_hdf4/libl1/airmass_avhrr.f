C
C@***s* PROJECT_PFSST/l2gen_both/airmass_avhrr.f
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
C         airmass_avhrr.f 
C
C LOCATION
C $OCSSWROOT 
C
C PURPOSE
C          Calculate atmospheric pathlength from viewing geometry.
C
C DESCRIPTION
C         The subroutine airmass_avhrr calculates the airmass from sun to ground
C        using Kasten's formula.
C       NOTE: All angles must be in radians.
C
C NOAA PFSST-SEADAS BUILD VERSION
C Pathfinder SST V5.2 code built with SEADAS version 6.3 64 bit l2gen_both for  
C CDR processed at University of Miami/RSMAS.
C
C PRIMARY SEADAS CODE DOCUMENTATION NASA
C For complete documentation of multi sensor SEADAS code see 
C http://seadas.gsfc.nasa.gov/doc/l2gen/l2gen.html 
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
C
C INPUTS
C          X5(real) - Angle PHI   (solar zenith angle)
C        
C
C OUTPUTS
C         X7(real) - Air mass 
C      
C  
C
C LANGUAGE
C Fortran   
C
C@*****
c !F90
c
c
c   subroutine AIRMASS_AVHRR(X5,X7)
c   !Description:
c        The subroutine airmass_avhrr calculates the airmass from sun to ground
c        using Kasten's formula.
c           NOTE: All angles must be in radians.
c   !Input Parameters:
c        X5(real) - Angle PHI      (solar zenith angle)
c   !Output Parameters:
c        X7(real) - Air mass
c
c $Id: airmass_avhrr.f,v 1.4 2012/05/07 20:10:09 sue Exp $
c
c $Log: airmass_avhrr.f,v $
c Revision 1.4  2012/05/07 20:10:09  sue
c Add or modify the PFSST headers.
c
c Revision 1.3  2012/04/26 18:55:31  sue
c     Changes made by seadas group for seadas6.3, newer fortran, and/or 64 bit mode.
c
c Revision 1.2  2011/10/17 15:52:36  kay
c Added NOAA headers
c
c Revision 1.1  2010/08/07 18:44:27  sue
c seadas 6.1 l2gen with modis dust and avhrr.
c
c Revision 1.1  2008/08/15 20:30:20  sue
c NODC versions of the pathfinder processing programs.
c
c Revision 1.1  2008/01/22 20:17:47  sue
c Initial version of pathnlch which reads the aci hdf data and geolocation files
c and creates hdf files.  This version was made just to test the aci hdf files.
c
c Revision 1.3  2002/08/23 17:42:08  sue
c Update copyright notices to year 2002.
c
c Revision 1.2  2000/06/14 16:08:46  jim
c Rearrange code to eliminate divide by zero.  Simplify logic.  Add comments.
c
c Revision 1.1  2000/03/07 22:58:53  jim
c Routine to compute airmass for zenith angles near 90 degrees.
c
c

c !Team-Unique Header:
C
C Copyright 1988-2002 by Rosenstiel School of Marine and Atmospheric Science,
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
c
c !References and Credits:
c        Written by:
c        University of Miami
c        Rosensteil School for Marine and Atmospheric Science
c        Division of Meteorology and Physical Oceanography
c        4600 Rickenbacker Cswy
c        Miami,Fl
c        33149
C        contact:SWalsh@rsmas.miami.edu
c
c !Design Notes:
c
c !END########################################################

        SUBROUTINE AIRMASS_AVHRR(X5,X7)
C
C  X5 IS THE SOLAR ZENITH ANGLE
C  X7 IS THE AIR MASS CORRECTED FOR ATMOSPHERIC REFRACTION
C      USING KASTEN'S FORMULA
C
C   ALL ANGLES MUST BE IN RADIANS
C
        real        X5, X7
        real        CT, ST, S, Z, R4, R5, P1, P2
        real        ANG, RAD, xx
        integer        IU

        real        A1, A2, A3
        data    A1, A2, A3/-.001867,-.002875,-.0008083/
        save        A1, A2, A3

        ANG(xx) = xx*180./3.14159265
        RAD(xx) = xx*3.14159265/180.

C
C  IS APPARENT SUN BELOW THE HORIZON?
C
C  IF IT IS, SET AIRMASS TO ZERO AS A FLAG.
C
        IF (ANG(X5).GT.90.57) THEN
          X7 = 0.
          GO TO 1
        END IF
C
        CT = COS(X5)
C       ST = SIN(X5)
C       IF (ANG(X5).GT.90.57) GO TO 1
        IF (CT.LT.1./3.86) GO TO 2
C
C  AWAY FROM HORIZON
C
        S = 1./CT
C       IF (S.GT.3.86 .OR. S.LT.0.) GO TO 2
        ST = S-1.
        X7 = S+A1*ST+A2*ST**2+A3*ST**3
        GO TO 1
C
C  VERY CLOSE TO HORIZON
C
2       CONTINUE
        R4 = 0.
        Z = ANG(X5)
        DO 3 IU=1,4
        P1 = Z-R4
        P2 = -116.94+4.41925*P1-.056623*P1*P1+.00024364*P1*P1*P1
        R4 = 10.**P2
3       CONTINUE
        R4 = Z-R4
        R5 = COS(RAD(R4))+.15*(90-R4+3.885)**(-1.253)
        X7 = 1./R5
1       CONTINUE
        RETURN
        END
