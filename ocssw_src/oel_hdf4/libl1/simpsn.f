C
C@***s* PROJECT_PFSST/l2gen_both/simpsn.f
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
C         simpsn.f 
C
C LOCATION
C $OCSSWROOT 
C
C PURPOSE
C          Calculate atmospheric pathlength from viewing geometry.
C
C DESCRIPTION
C         The function simpsn is used to integrate over a given set of ordinates for
C         equally spaced abscissas using simpson(s) method or using entry point simpne to integrate
C         analytically a 3-point lagrangian interpolation polynomial fitting the ordinates.
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
C          X(real) - FOR SIMPSN, THIS IS A SINGLE VARIABLE SPECIFYING THE INCREMENT
C                      FOR THE EQUALLY SPACED ABSCISSAS.
C                  - FOR SIMPNE, THIS IS A VECTOR OF LENGTH NUM
C                      SPECIFYING THE UNEQUALLY SPACED ABSCISSAS.
C
C          Y(real) - VECTOR OF ORDINATES.
C
C          NUM(integer) - LENGTH OF VECTOR Y (SHOULD BE .GT. 2).
C        
C
C OUTPUTS
C          IER(integer) - ERROR FLAG.
C                          = 32  NUM IS .LT. 3, NO EXECUTION.
C                          =  0  NO ERROR.
C      
C  
C
C LANGUAGE
C Fortran   
C
C@*****
c*************************************************************
c 
c !F77
c
c FUNCTION SIMPSN (X,Y,NUM,IER)
c !Description:
c    TO INTEGRATE OVER A GIVEN SET OF ORDINATES FOR
C    EQUALLY SPACED ABSCISSAS USING SIMPSON(S)
c    METHOD OR USING ENTRY POINT SIMPNE TO INTEGRATE
C    ANALYTICALLY A 3-POINT LAGRANGIAN INTERPOLATION
C    POLYNOMIAL FITTING THE ORDINATES.
c !Input Parameters:
c           REAL  X
C                        . FOR SIMPSN, THIS IS A SINGLE VARIABLE
C                          SPECIFYING THE INCREMENT FOR THE EQUALLY
C                          SPACED ABSCISSAS.
C                        . FOR SIMPNE, THIS IS A VECTOR OF LENGTH NUM
C                          SPECIFYING THE UNEQUALLY SPACED ABSCISSAS.
C
C            REAL Y
C                          VECTOR OF ORDINATES.
C
C            INTEGER   NUM
C                          LENGTH OF VECTOR Y (SHOULD BE .GT. 2).
C
c
c !Output Parameters:
c            INTEGER   IER
C                          ERROR FLAG.
C                          = 32  NUM IS .LT. 3, NO EXECUTION.
C                          =  0  NO ERROR.
c
c !Revision History:
c
c
c $Id: simpsn.f,v 1.4 2012/06/14 21:54:49 sue Exp $
c
c $Log: simpsn.f,v $
c Revision 1.4  2012/06/14 21:54:49  sue
c Fix extra comment starts in the .c files, and remove tabs in the .f
c files to get rid of compiler warnings.
c
c Revision 1.3  2012/05/07 20:10:15  sue
c Add or modify the PFSST headers.
c
c Revision 1.2  2012/04/26 18:55:32  sue
c     Changes made by seadas group for seadas6.3, newer fortran, and/or 64 bit mode.
c
c Revision 1.1  2010/08/07 18:44:49  sue
c seadas 6.1 l2gen with modis dust and avhrr.
c
c Revision 1.1  2008/08/15 20:30:36  sue
c NODC versions of the pathfinder processing programs.
c
c Revision 1.3  2002/08/23 20:53:07  sue
c Update copyright notices to year 2002.
c
c Revision 1.2  1997/11/20 19:39:20  jim
c Fix prologs.
c
c Revision 1.1  1996/09/12 17:44:20  sue
c Modis version of ulib.
c
c Revision 1.1  1992/02/26 19:57:42  angel
c Initial revision
c
c
c !Team-unique Header: 
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
c
c       Written by:
c       University of Miami
c       Rosensteil School for Marine and Atmospheric Science
c       Division of Meteorology and Physical Oceanography
c       4600 Rickenbacker Cswy
c       Miami,Fl
c       33149
c
c       contact: SWalsh@rsmas.miami.edu
c
c !Design Notes:
c
c !ENDcccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
C     NCAR SCIENTIFIC SUBROUTINE LIBRARY ROUTINE@ SIMPSN
      FUNCTION SIMPSN (X,Y,NUM,IER)
C
C
C DIMENSION OF           X(NUM),Y(NUM)
C ARGUMENTS
C
C LATEST REVISION        MAY, 1974
C
C PURPOSE                TO INTEGRATE OVER A GIVEN SET OF ORDINATES FOR
C                        EQUALLY SPACED ABSCISSAS USING SIMPSON(S)
C                        METHOD OR USING ENTRY POINT SIMPNE TO INTEGRATE
C                        ANALYTICALLY A 3-POINT LAGRANGIAN INTERPOLATION
C                        POLYNOMIAL FITTING THE ORDINATES.
C
C ACCESS CARDS           *FORTRAN,S=ULIB,N=SIMPSN
C                        *COSY
C
C USAGE                  FOR EQUALLY SPACED ABSCISSAS, USE
C                            YINT = SIMPSN(X,Y,NUM,IER)  .
C                        FOR UNEQUALLY SPACED ABSCISSAS, USE
C                            YINT = SIMPNE(X,Y,NUM,IER)  .
C
C ARGUMENTS
C
C ON INPUT               X
C                        . FOR SIMPSN, THIS IS A SINGLE VARIABLE
C                          SPECIFYING THE INCREMENT FOR THE EQUALLY
C                          SPACED ABSCISSAS.
C                        . FOR SIMPNE, THIS IS A VECTOR OF LENGTH NUM
C                          SPECIFYING THE UNEQUALLY SPACED ABSCISSAS.
C
C                        Y
C                          VECTOR OF ORDINATES.
C
C                        NUM
C                          LENGTH OF VECTOR Y (SHOULD BE .GT. 2).
C
C ON OUTPUT              IER
C                          ERROR FLAG.
C                          = 32  NUM IS .LT. 3, NO EXECUTION.
C                          =  0  NO ERROR.
C
C ENTRY POINTS           SIMPSN, SIMPNE
C
C SPECIAL CONDITIONS     IF NUM IS AN EVEN NUMBER, SIMPSON(S) RULE IS
C                        APPLIED TO INTERVAL X(1) THROUGH X(NUM-1) AND
C                        THE RESULT IS ADDED TO THE INTEGRAL OVER
C                        X(NUM-1) TO X(NUM) OF THE THREE-POINT
C                        LAGRANGIAN INTERPOLATING POLYNOMIAL THROUGH
C                        X(NUM-2), X(NUM-1) AND X(NUM).
C
C COMMON BLOCKS          NONE
C
C I/O                    NONE
C
C PRECISION              SINGLE
C
C REQUIRED ULIB          NONE
C ROUTINES
C
C SPECIALIST             WILLIAM B. FRYE, NCAR, BOULDER, COLORADO  80303
C
C LANGUAGE               FORTRAN
C
C HISTORY                STANDARDIZED MAY 1974.
C
C
C
C
C     IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      INTEGER     NUM
      INTEGER     IER, LIM, I
      REAL*8      SSUM
      REAL        X(NUM)     ,Y(NUM)     ,SIMPSN, SIMPNE
      REAL*8      DEL(3)     ,PI(3)      ,G(3)
      REAL*8      THIRD, HALF, TWO3R, E, F, FEINTS, DELPRD
      REAL*8      X1, X2, X3
      INTEGER     N
      DATA THIRD,HALF,TWO3R/0.333333333333333,0.5,0.666666666666666/

      IF (NUM .GT. 2) GO TO 101
      IER = 32
      SIMPSN = 0.
      RETURN
  101 CONTINUE
      IER = 0
      N = NUM
      IF (MOD(NUM,2) .EQ. 0) N = NUM-1
      SIMPSN = HALF*(Y(1)+Y(N))+2.*Y(N-1)
      IF (N .EQ. 3) GO TO 103
      LIM = N-2
      DO 102 I=2,LIM,2
         SIMPSN = SIMPSN+2.*Y(I)+Y(I+1)
  102 CONTINUE
  103 SIMPSN = TWO3R*X(1)*SIMPSN
      IF (MOD(NUM,2) .EQ. 1) RETURN
      E = X(1)*X(1)
      F = X(1)*E
      FEINTS = X(1)
      DEL(1) = X(1)
      DEL(2) = -2.*X(1)
      DEL(3) = X(1)
      G(1) = X(1)
      G(2) = 0.0
      G(3) = -X(1)
      PI(1) = 0.0
      PI(2) = -E
      PI(3) = 0.0
      DELPRD = -2.*F
      N = N-1
      GO TO 107
      ENTRY SIMPNE (X,Y,NUM,IER)
C
C X AND Y ARE COORDINATES OF THE POINTS (IN ARRAY FORM)
C
      IF (NUM .GT. 2) GO TO 104
      IER = 32
      SIMPSN = 0.
      RETURN
  104 CONTINUE
      IER = 0
      SIMPSN = 0.
      N = 1
  105 CONTINUE
      X1 = X(N)
      X2 = X(N+1)
      X3 = X(N+2)
      E = X3*X3-X1*X1
      F = X3*X3*X3-X1*X1*X1
      FEINTS = X3-X1
      DEL(1) = X3-X2
      DEL(2) = -FEINTS
      DEL(3) = X2-X1
      G(1) = X2+X3
      G(2) = X1+X3
      G(3) = X1+X2
      PI(1) = X2*X3
      PI(2) = X1*X3
      PI(3) = X1*X2
      DELPRD = DEL(1)*DEL(2)*DEL(3)
      SSUM = 0.0d0
      DO 106 I=1,3
         SSUM = SSUM+Y(N-1+I)*DEL(I)*(THIRD*F-G(I)*HALF*E+PI(I)*FEINTS)
  106 CONTINUE
      SIMPSN = SIMPSN-SSUM/DELPRD
      N = N+2
      IF (NUM .GT. (N+1)) GO TO 105
      IF (MOD(NUM,2) .NE. 0) RETURN
      N = NUM-2
      X3 = X(NUM)
      X2 = X(NUM-1)
      X1 = X(NUM-2)
      E = X3*X3-X2*X2
      F = X3*X3*X3-X2*X2*X2
      FEINTS = X3-X2
      DEL(1) = FEINTS
      DEL(2) = X1-X3
      DEL(3) = X2-X1
      G(1) = X2+X3
      G(2) = X1+X3
      G(3) = X1+X2
      PI(1) = X2*X3
      PI(2) = X1*X3
      PI(3) = X1*X2
      DELPRD = DEL(1)*DEL(2)*DEL(3)
  107 SSUM = 0.0d0
      DO 108 I=1,3
         SSUM = SSUM+Y(N-1+I)*DEL(I)*(THIRD*F-G(I)*HALF*E+PI(I)*FEINTS)
  108 CONTINUE
      SIMPSN = SIMPSN-SSUM/DELPRD
      RETURN
      END
