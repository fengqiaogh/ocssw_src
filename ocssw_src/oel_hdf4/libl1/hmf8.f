C
C@***s* PROJECT_PFSST/l2gen_both/hmf8.f
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
C       hmf8.f
C
C LOCATION
C $OCSSWROOT
C
C PURPOSE
C       Calculate the normalized aerosol radiance (lacking tau_a).
C
C DESCRIPTION
C       The subroutine hmf8 is the algorithm for a non-absorbing standard marine
C       aerosol with 80% humidity.
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
C       X1(real) - Satellite zenith angle
C       X2(real) - Sun zenith angle
C       X3(real) - Delta azimuth
C       x4(real) - Optical thickness (unused)
C       X5(real) - Solar constant
C
C
C OUTPUTS
C       X6(real) - Normalized aerosol radiance (lacking tau_a)
C
C
C LANGUAGE
C Fortran
C
C@*****
C
C !F90
C
C
C   subroutine HMF8 (X1,X2,X3,X4,X5,X6)
C   !Description:
C       The subroutine hmf8 is the algorithm for a non-absorbing standard marine
C       aerosol with 80% humidity.
C       NOTE: All angles must be in radians.
C   !Input Parameters:
C       X1(real) - Satellite zenith angle
C       X2(real) - Sun zenith angle
C       X3(real) - Delta azimuth
C       x4(real) - Optical thickness (unused)
C       X5(real) - Solar constant
C   !Output Parameters:
C       X6(real) - Normalized aerosol radiance (lacking tau_a)
C
C !Revision History:
C
C $Id: hmf8.f,v 1.3 2012/05/07 20:10:13 sue Exp $
C
C $Log: hmf8.f,v $
C Revision 1.3  2012/05/07 20:10:13  sue
C Add or modify the PFSST headers.
C
C Revision 1.2  2012/04/26 18:55:32  sue
C     Changes made by seadas group for seadas6.3, newer fortran, and/or 64 bit mode.
C
C Revision 1.1  2010/08/07 18:44:36  sue
C seadas 6.1 l2gen with modis dust and avhrr.
C
C Revision 1.1  2008/08/15 20:30:27  sue
C NODC versions of the pathfinder processing programs.
C
C Revision 1.1  2008/01/22 20:17:49  sue
C Initial version of pathnlch which reads the aci hdf data and geolocation files
C and creates hdf files.  This version was made just to test the aci hdf files.
C
C Revision 1.4  2002/08/23 17:42:09  sue
C Update copyright notices to year 2002.
C
C Revision 1.3  1999/03/26 16:02:04  sue
C Change copyright date to 1999.
C
C Revision 1.2  1996/09/06 13:32:58  kay
C update prologs
C
C Revision 1.1  1996/09/05  12:48:01  sue
C New library which contains all of the routines and include files which
C are common to both modcol and modsst (except the routines which have
C to do with binning which are in binshr).
C
C Revision 1.5  1996/05/08 15:44:39  kay
C  update prolog, add bang
C
C Revision 1.4  1996/03/11  14:20:10  kay
C update prologs
C
C Revision 1.3  1995/12/04  14:31:56  kay
C dd prolog
C
C Revision 1.2  1992/11/09  20:17:22  sue
C Make all constants double or single where appropriate, and use functions
C to specifically convert types.
C
C Revision 1.1  1992/02/26  19:47:11  angel
C Initial revision
C

C !Team-Unique Header:
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
C
C !References and credits:
C        Written by:
C         University of Miami
C        Rosensteil School for Marine and Atmospheric Science
C        Division of Meteorology and physical Oceanography
C        4600 Rickenbacker Cswy
C        Miami,Fl
C        33149
C        Contact: SWalsh@rsmas.miami.edu
C
C !Design Notes:
C
C !END####################################################################

        subroutine HMF8 (X1,X2,X3,X4,X5,X6)
C
C        HMF8: Standard marine aerosol, 80% humidity, non-absorbing.
C
C  X1  Satellite zenith angle
C  X2  Sun zenith angle
C  X3  Delta azimuth
C  X4  Optical thickness (unused)
C  X5  Solar constant
C  X6  Normalized aerosol radiance (lacking tau_a)
C
C  ALL ANGLES MUST BE IN RADIANS
C
        real PI
        parameter (PI = 3.14159265358979)

        real        X1, X2, X3, X4, X5, X6
        real*8      alpha
        real*8      g1
        real*8      g2

        real*8      Pa,        ctheta
        real*8      f,        g
        real*8      A,        B
        real        R1,        R2

        data        alpha        / 0.983d0/
        data        g1        / 0.82d0 /
        data        g2        /-0.55d0 /

C        ASF's

        f(ctheta,g) = (1d0-g*g)/(1d0+g*g-2d0*g*ctheta)**1.5d0
        Pa(ctheta) = alpha*f(ctheta,g1)+(1d0-alpha)*f(ctheta,g2)

C        Start of code

C Outgoing reflectivity
        call REFLEC(X1,R1)
C Incoming reflectivity
        call REFLEC(X2,R2)

C cos(theta-)
        A  = dble(-COS(X1)*COS(X2)+SIN(X1)*SIN(X2)*COS(X3))
C cos(theta+)
        B  = dble(+COS(X1)*COS(X2)+SIN(X1)*SIN(X2)*COS(X3))

        X6 = X5/(4*PI*cos(X1)) * (sngl(Pa(A)) + (R1+R2)*sngl(Pa(B)))
        return
        end
