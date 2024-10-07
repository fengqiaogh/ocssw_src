      SUBROUTINE IDENTY                                                         
     I          (EPTIME, IMATCH, DANGTL, DMAGTL, PANGTL, TANGTL,                
     I           TMINCO, SMAGLM, IQLIMT, MAXCAT, IDNCAT, DATCAT,                
     I           IFXCLM, NUMCLM, TIMCLM, GCICLM, BRICLM, VRMCLM,                
     I           VRPCLM, NOBCLM, MRKCLM, NUMSTR, KLMSTR,                        
     O           MRKSTR, LBLCLM, IDFCLM, NRFCLM, MAPCLM, SKYCLM, IDFHST,        
     O           IRCODE)                                                        
C-----------------------------------------------------------------------        
C   MODULE NAME: STIDENTY                                                       
C                                                                               
C   PURPOSE: TO IDENTIFY CLUMPS WITH STARS IN THE REFERENCE CATALOG             
C            AND TO CORRECT MATCHES FOR EARTH AND S/C VELOCITY.                 
C                                                                               
C   ARGUMENT LIST:                                                              
C   ARGUMENT I/O TYPE DIMENS DESCRIPTION                                        
C   -------- --- ---- ------ -----------                                        
C   EPTIME   I   R*8         EPOCH TIME FOR SOLUTION (SECS SINCE 1/1/72)        
C   IMATCH   I   I*4         MATCH TYPE (1=DIRECT,2=DOUBLET,3=TRIPLET)          
C   DANGTL   I   R*8         DIRECT MATCH MAX ANGULAR DIF TOLERANCE             
C   DMAGTL   I   R*8         DIRECT MATCH MAX MAGNITUDE DIF TOLERANCE           
C   PANGTL   I   R*8         PAIRWISE MATCH MAX ANGULAR DIF TOLERANCE           
C   TANGTL   I   R*8         TRIPLET MATCH MAX ANGULAR TOLERANCE                
C   TMINCO   I   R*8         TRIPLET MATCH MINIMUM COLINEARITY ANGLE            
C   SMAGLM   I   R*4         STAR MAGNITUDE LIMIT (MINIMUM MAGNITUDE)           
C   IQLIMT   I   I*4    6    CATALOG STAR QUALITY LIMITS                        
C                                (1=VARIABLILITY,2=COLOR,3=MULTIPLICITY,        
C                                 4=NEAR NEIGHBORS, 5=POSITION ERROR,           
C                                 6=MAGNITUDE ERROR)                            
C   MAXCAT   I   I*4         MAX NUMBER OF REFERENCE STARS IN CATALOG           
C   IDNCAT   I   I*4    *    REFERENCE STAR SKYMAP ID NUMBERS                   
C   DATCAT   I   R*4   7,*   CATALOG DATA (1=X AXIS,2=Y AXIS,3=Z AXIS,          
C                                          4=MAGNITUDE,5=MOTION,                
C                                          6=QUALITY,7=COLOR)                   
C   IFXCLM   I   I*4         MAX NUMBER OF REF STAR MATCHES PER CLUMP           
C                            (I.E. 1ST DIM OF MAPCLM & SKYCLM ARRAYS)           
C   NUMCLM   I   I*4         NUMBER OF CLUMPS                                   
C   TIMCLM   I   R*8    *    AVERAGE CLUMP TIMES                                
C   GCICLM   I   R*4   3,*   AVERAGE CLUMP POSITIONS GCI                        
C   BRICLM   I   R*4    *    AVERAGE CLUMP MAGNITUDES                           
C   VRMCLM   I   R*4    *    CLUMP MAGNITUDE VARIANCE                           
C   VRPCLM   I   R*4    *    CLUMP POSITION VARIANCE                            
C   NOBCLM   I   I*4    *    NUMBER OF OBSERVATIONS PER CLUMP                   
C   MRKCLM   I   I*4    *    CLUMP STATUS FLAGS                                 
C   NUMSTR   I   I*4         NUMBER OF OBSERVATIONS                             
C   KLMSTR   I   I*4    *    CLUMP NUMBER FOR EACH OBSERVATION                  
C   MRKSTR   I O I*4    *    STATUS FLAG FOR EACH OBSERVATION                   
C   LBLCLM   I O I*4    *    CLUMP LABELS                                       
C   IDFCLM     O I*4    *    CLUMP IDENTIFICATION FLAGS                         
C   NRFCLM     O I*4    *    NUMBER OF REFERENCE STARS MATCHES PER CLUMP        
C   MAPCLM     O I*4  10,*   LIST OF REFERENCE STAR #'S MATCHED TO CLUMP        
C   SKYCLM     O R*4 10,3,*  LIST OF CORRECTED REFERENCE STAR POSITIONS         
C                                    (MATCH#,AXIS,CLUMP#)                       
C   IDFHST    IO I*4    *    FHST NUMBER FOR EACH CLUMP                         
C   IRCODE     O I*4         ERROR FLAG (0=NO ERROR, 1= ERROR)                  
C                                                                               
C                                                                               
C   COMMON BLOCK VARIABLES USED:                                                
C   COMMON   VAR    I/O   VAR    I/O   VAR    I/O   VAR    I/O                  
C   ------   ---    ---   ---    ---   ---    ---   ---    ---                  
C   CMLUNS   LUCAT  I                                                           
C   CMDEBG   LEVDBG I     LUDBUG I                                              
C   CMSMSG   IMSGNM I     IVARLN I     IDSTFG I     IRC    I                    
C            C$VDAT I     C$SBID I                                              
      IMPLICIT NONE                                                             
C       ++INCLUDE STCMDEBG                                                       
C       ++INCLUDE STCMLUNS                                                       
C       ++INCLUDE STCMSMSG                                                       
      INTEGER*4 LEVDBG(8),LUDBUG
      COMMON /CMDEBG/LEVDBG,LUDBUG
C                                                                               
C                                                                               
C   EXTERNAL FILES REFERENCED:                                                  
C   FILENAME      OPERATION   FORTRAN UNIT ID                                   
C   --------      ---------   ---------------                                   
C   NONE                                                                        
C                                                                               
C   EXTERNAL REFERENCES:                                                        
C   --------------------------------------------------------------------        
C   CORECT - CORRECT REFERENCE STAR MATCHES FOR SPACECRAFT AND                  
C              EARTH VELOCITY ABBERATION                                        
C   DMATCH - DIRECT MATCH ALGORITHM                                             
C   DOUBLT - DOUBLET MATCH ALGORITHM                                            
C   GETCAT - GET STAR CATALOG DATA                                              
C   SORTCL - SORT CLUMPS BY RESULTS OF DIRECT MATCH                             
C   TRIPLT - TRIPLET MATCH ALGORITHM                                            
C   UTMSG  - MESSAGE OUTPUT UTILITY                                             
C                                                                               
C   SUBROUTINE CALLED FROM:                                                     
C   --------------------------------------------------------------------        
C     STIDSTAR - STAR IDENTIFICATION PROCESSING DRIVER                          
C                                                                               
C   CONSTRAINTS, RESTRICTIONS, MESSAGES, NOTES:                                 
C   --------------------------------------------------------------------        
C   NONE                                                                        
C                                                                               
C   REQUIREMENTS REFERENCES:                                                    
C   --------------------------------------------------------------------        
C   UARS FDSS SPECS PAGES 3.1.1.5-17 TO 3.1.1.5-19 (F3-1A TO F3-4)              
C                                                                               
C   DEVELOPMENT HISTORY:                                                        
C   DATE      AUTHOR              DESCRIPTION                                   
C   --------  ------              -----------                                   
C   8 / 1/88  R.J. BURLEY         DESIGN                                        
C   5 /17/89  R.J. BURLEY         CODED                                         
C   10/10/89  R.J. BURLEY         DIRECT MATCH ALGORITHM NAME CHANGED           
C                                 FROM DIRECT TO DMATCH, BECAUSE DIRECT         
C                                 IS A GESS RESERVED WORD.                      
C   11/ 2/89  R.J. BURLEY         USE UTMSG TO ADVISE USER WHEN MATCHING        
C                                 ALGORITHM DEFAULTS FROM IMATCH.               
C   11/ 6/89  R.J. BURLEY         ADD LBLCLM GESS ARRAY                         
C   02/ 3/92  C. C. YEH           ADD IDFHST (MTASS-11)                         
C    9/30/92  D. MUCCI            MTASS 151  REMOVED EPTIME FROM THE            
C                                 CALL TO STCORECT AND ADDED MRKCLM.            
C-----------------------------------------------------------------------        
C   METHOD:                                                                     
C     CALL GETCAT TO GET STAR CATALOG                                           
C     IF (AN ERROR OCCURED IN GETCAT) WRITE MESSAGE AND ABORT                   
C                                                                               
C     CALL DMATCH TO PERFORM DIRECT MATCH ALGORITHM                             
C     IF (AN ERROR OCCURED IN DMATCH) WRITE MESSAGE AND ABORT                   
C                                                                               
C     CALL CORECT TO PERFORM VELOCITY ABBERATION CORRECTION                     
C     IF (AN ERROR OCCURED IN CORECT) WRITE MESSAGE AND ABORT                   
C                                                                               
C     IF (NUMCLM > 1) THEN                                                      
C       CALL SORTCL TO SORT CLUMPS IN ASCEND ORDER OF REF STAR MATCHES          
C     ENDIF                                                                     
C                                                                               
C     IF ((IMATCH=3).AND.(NUMCLM >= 3)) THEN                                    
C        CALL TRIPLT TO PERFORM TRIPLET MATCH ALGORITHM                         
C     ENDIF                                                                     
C                                                                               
C     IF ((IMATCH=2).AND.(NUMCLM >= 2)) THEN                                    
C        CALL DOUBLT TO PERFORM DOUBLET MATCH ALGORITHM                         
C        IF (# VALID DOUBLETS = 0) THEN                                         
C          CALL UTMSG TO ADVISE USER THAT OUTPUT WILL BE FROM DMATCH            
C        ENDIF                                                                  
C     ELSE IF ((IMATCH=3).AND.( # VALID TRIPLETS = 0).AND.                      
C              (NUMCLM >= 2)) THEN                                              
C        CALL UTMSG TO ADVISE USER THAT STARID IS DEFAULTING TO DOUBLT          
C        CALL DOUBLT TO PERFORM DOUBLET MATCH ALGORITHM                         
C        IF (# VALID DOUBLETS = 0) THEN                                         
C          CALL UTMSG TO ADVISE USER THAT OUTPUT WILL BE FROM DMATCH            
C        ENDIF                                                                  
C     ENDIF                                                                     
C     RETURN                                                                    
C-----------------------------------------------------------------------        
C                                                                               
C     * DEFINE PARAMETER VARIABLES                                              
      REAL*8    EPTIME      , DANGTL        , DMAGTL   ,  PANGTL                
      REAL*8    TANGTL      , TMINCO        , TIMCLM(*)                         
C                                                                               
      REAL*4    SMAGLM      , DATCAT(7,*)                                       
      REAL*4    GCICLM(3,*) , BRICLM(*)     , VRMCLM(*)                         
      REAL*4    VRPCLM(*)   , SKYCLM(10,3,*)                                    
C                                                                               
      INTEGER*4 IQLIMT(6)   , IMATCH        , MAXCAT   , IDNCAT(*)              
      INTEGER*4 NUMCLM      , NOBCLM(*)     , MRKCLM(*), NUMSTR                 
      INTEGER*4 KLMSTR(*)   , MRKSTR(*)     , IDFCLM(*), NRFCLM(*)              
      INTEGER*4 LBLCLM(*)   , MAPCLM(10,*)  , IFXCLM   , IRCODE                 
      INTEGER*4 IDFHST(*)                                                       
C                                                                               
C     * DECLARE LOCAL VARIABLES                                                 
      INTEGER*4 IERR        , NUMTRP        , NUMDUB   , NUMCAT  ,LUCAT
C                                                                               
C                                      INITIALIZE ROUTINE                       
      IF (LEVDBG(7) .NE. 0) WRITE (LUDBUG,1000)                                 
      IRCODE = 0                                                                
C                                                                               
CC   -GET STAR DATA FROM CATALOG                                                
C                                                                               
      CALL GETCAT(LUCAT , EPTIME, MAXCAT, SMAGLM, IQLIMT,                       
     1            NUMCAT, IDNCAT, DATCAT, IERR)                                 
      IF (IERR .NE. 0) THEN                                                     
        IF (LEVDBG(7) .NE. 0) WRITE (LUDBUG,6000) IERR                          
        IRCODE = 1                                                              
        GO TO 9999                                                              
      ENDIF                                                                     
C                                                                               
CC   -PERFORM DIRECT MATCH                                                      
C                                                                               
      CALL DMATCH (DANGTL, DMAGTL, IFXCLM, NUMCAT,                              
     1             IDNCAT, DATCAT, NUMSTR, KLMSTR,                              
     2             NUMCLM, LBLCLM, GCICLM, BRICLM, NOBCLM,                      
     3             MRKCLM, MRKSTR, IDFCLM, NRFCLM,                              
     4             MAPCLM, SKYCLM, IERR)                                        
      IF (IERR .NE. 0) THEN                                                     
        IF (LEVDBG(7) .NE. 0) WRITE (LUDBUG,6001) IERR                          
        IRCODE = 1                                                              
        GO TO 9999                                                              
      ENDIF                                                                     
C                                                                               
CC   -PERFORM VELOCITY ABBERATION CORRECTION ON MATCHED STARS                   
C                                                                               
      CALL CORECT(IDNCAT, DATCAT, NUMCLM, LBLCLM, TIMCLM,                      
     *            NRFCLM, MAPCLM, MRKCLM, SKYCLM, IERR)                        
      IF (IERR .NE. 0) THEN                                                    
        IF (LEVDBG(7) .NE. 0) WRITE (LUDBUG,6002) IERR                         
        IRCODE = 1                                                             
        GO TO 9999                                                             
      ENDIF                                                                    
C                                                                               
CC   -SORT CLUMPS BY # OF MATCHES                                               
C                                                                               
      IF (NUMCLM .GT. 1) THEN                                                   
        CALL SORTCL (NUMCLM, NUMSTR, IFXCLM,                                    
     1           LBLCLM, TIMCLM, BRICLM, GCICLM, VRMCLM, VRPCLM, NOBCLM,        
     2           MRKCLM, IDFCLM, NRFCLM, MAPCLM, SKYCLM, KLMSTR, IDFHST)        
      ENDIF                                                                     
C                                                                               
CC   -PERFORM TRIPLET MATCHING                                                  
C                                                                               
      IF ((IMATCH .EQ. 3).AND.(NUMCLM .GE. 3)) THEN                             
        CALL TRIPLT (TANGTL, TMINCO, NUMCLM, LBLCLM, GCICLM, NOBCLM,            
     1               MRKCLM, NUMSTR, KLMSTR, MRKSTR, SKYCLM,                    
     2               NRFCLM, MAPCLM, IDFCLM, NUMTRP)                            
        IF (LEVDBG(7) .GE. 4) WRITE (LUDBUG,4000) NUMTRP                        
      ENDIF                                                                     
C                                                                               
CC   -PERFORM PAIRWISE MATCHING                                                 
C                                                                               
      IF ((IMATCH .EQ. 2).AND.(NUMCLM .GE. 2)) THEN                             
        CALL DOUBLT (PANGTL, NUMCLM, NUMSTR, KLMSTR, LBLCLM, GCICLM,            
     1               MRKCLM, NRFCLM, MAPCLM, SKYCLM, MRKSTR, IDFCLM,            
     2               NUMDUB)                                                    
C        IF (NUMDUB .EQ. 0) THEN                                                
C          IMSGNM = 63                                                          
C          IVARLN = 0                                                           
C          CALL UTMSG (C$SBID, IMSGNM, IVARLN, C$VDAT, IDSTFG, IRC)             
C        ENDIF                                                                  
        IF (LEVDBG(7) .GE. 4) WRITE (LUDBUG,4001) NUMDUB                        
      ELSE IF ((IMATCH .EQ. 3).AND.(NUMTRP .EQ. 0).AND.                         
     1         (NUMCLM .GE. 2)) THEN                                            
C        IMSGNM = 62                                                            
C        IVARLN = 0                                                             
C        CALL UTMSG (C$SBID, IMSGNM, IVARLN, C$VDAT, IDSTFG, IRC)               
        WRITE(*,*)'NO TRIPLETS FOUND -- ATTEMPTING PAIR MATCHING'         
        CALL DOUBLT (PANGTL, NUMCLM, NUMSTR, KLMSTR, LBLCLM, GCICLM,            
     1               MRKCLM, NRFCLM, MAPCLM, SKYCLM, MRKSTR, IDFCLM,            
     2               NUMDUB)                                                    
C        IF (NUMDUB .EQ. 0) THEN                                                
C          IMSGNM = 63                                                          
C          IVARLN = 0                                                           
C          CALL UTMSG (C$SBID, IMSGNM, IVARLN, C$VDAT, IDSTFG, IRC)             
C        ENDIF                                                                  
        IF (LEVDBG(7) .GE. 4) WRITE (LUDBUG,4001) NUMDUB                        
      ENDIF                                                                     
9999  CONTINUE                                                                  
      IF (LEVDBG(7) .NE. 0) WRITE (LUDBUG,2000)                                 
C                                                                               
C    -FORMAT SECTION                                                            
C                                                                               
1000  FORMAT(' *** ENTER IDENTY 92/09/30 ***')                                  
2000  FORMAT(' *** EXIT  IDENTY ***')                                           
4000  FORMAT(' TRIPLT MATCH: NUMBER OF VALID TRIPLETS = ',I8)                   
4001  FORMAT(' DOUBLT MATCH: NUMBER OF VALID DOUBLETS = ',I8)                   
6000  FORMAT(' ABEND IDENTY: ERROR RETURN FROM GETCAT = ',I8)                   
6001  FORMAT(' ABEND IDENTY: ERROR RETURN FROM DMATCH = ',I8)                   
6002  FORMAT(' ABEND IDENTY: ERROR RETURN FROM CORECT = ',I8)                   
      RETURN                                                                    
      END                                                                       
