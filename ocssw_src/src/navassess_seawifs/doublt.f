      SUBROUTINE DOUBLT                                                        
     I          (PANGTL, NUMCLM, NUMSTR, KLMSTR, LBLCLM, GCICLM,               
     O           MRKCLM, NRFCLM, MAPCLM, SKYCLM, MRKSTR, IDFCLM, NUMDUB)       
C-----------------------------------------------------------------------       
C   MODULE NAME:  STDOUBLT                                                     
C                                                                              
C                                                                              
C   PURPOSE: TO IDENTIFY STAR CLUMPS USING PAIRWISE STAR MATCHING.             
C                                                                              
C                                                                              
C   ARGUMENT LIST:                                                             
C   ARGUMENT I/O TYPE DIMENS DESCRIPTION                                       
C   -------- --- ---- ------ -----------                                       
C   PANGTL   I   R*8         MAX ANGULAR SEPARATION FOR MATCH (DEG)            
C   NUMCLM   I   I*4         NUMBER OF CLUMPS                                  
C   NUMSTR   I   I*4         NUMBER OF OBSERVATIONS                            
C   KLMSTR   I   I*4    *    CLUMP NUMBER FOR EACH OBSERVATION                 
C   LBLCLM   I   I*4    *    CLUMP LABELS                                      
C   GCICLM   I   R*4   3,*   AVG POSITION VECTOR (GCI) FOR EACH CLUMP          
C   MRKCLM   I O I*4    *    STATUS FLAG FOR EACH CLUMP                        
C   NRFCLM   I O I*4    *    NUMBER OF REFERENCE STARS MATCHED TO CLUMP        
C   MAPCLM   I O I*4  10,*   LIST OF REFERENCE STAR #'S MATCHED TO CLUMP       
C   SKYCLM   I O R*4 10,3,*  LIST OF CORRECTED REFERENCE STAR POSITIONS        
C   MRKSTR   I O I*4    *    REJECTION FLAG FOR EVERY OBSERVATION              
C   IDFCLM     O I*4    *    IDENTIFICATION FLAG FOR EACH CLUMP                
C   NUMDUB     O I*4         NUMBER OF DOUBLETS CHECKED                        
C                                                                              
C                                                                              
C   COMMON BLOCK VARIABLES USED:                                               
C   COMMON   VAR    I/O   VAR    I/O   VAR    I/O   VAR    I/O                 
C   ------   ---    ---   ---    ---   ---    ---   ---    ---                 
C   CMDEBG   LEVDBG I     LUDBUG I                                             
C   CMCONV   DTR    I                                                          
C                                                                              
C      ++INCLUDE STCMDEBG                                                      
C      ++INCLUDE AECMCONV                                                      
      INTEGER*4 LEVDBG(8),LUDBUG
      COMMON /CMDEBG/LEVDBG,LUDBUG
      REAL*8 PI,RADEG,RE,REM,F,OMF2,OMEGAE
      COMMON /GCONST/PI,RADEG,RE,REM,F,OMF2,OMEGAE
C                                                                              
C   EXTERNAL FILES REFERENCED:                                                 
C   FILENAME      OPERATION   FORTRAN UNIT ID                                  
C   --------      ---------   ---------------                                  
C   NONE                                                                       
C                                                                              
C   EXTERNAL REFERENCES:                                                       
C   --------------------------------------------------------------------       
C   UTDANGLE - FUNCTION TO COMPUTE ANGULAR SEPARATION BETWEEN 2 VECTORS        
C                                                                              
C   SUBROUTINE CALLED FROM:                                                    
C   --------------------------------------------------------------------       
C     STIDENTY - STAR MATCHING DRIVER                                          
C                                                                              
C                                                                              
C   CONSTRAINTS, RESTRICTIONS, MESSAGES, NOTES:                                
C   --------------------------------------------------------------------       
C   NONE                                                                       
C                                                                              
C   REQUIREMENTS REFERENCES:                                                   
C   --------------------------------------------------------------------       
C   UARS FDSS SPECS 3.1.1.5 FUNCTION F3                                        
C                                                                              
C                                                                              
C   DEVELOPMENT HISTORY:                                                       
C   DATE     AUTHOR              DESCRIPTION                                   
C   -------- ------              -----------                                   
C   8/ 2/88  R.J. BURLEY         DESIGN                                        
C   5/16/89  R.J. BURLEY         CODED                                         
C  10/17/89  R.J. BURLEY         ADD FLAG UNIDENTIFIED CLUMPS, AND             
C                                CORRECT FLAGGING OF OBSERVATIONS.             
C  10/19/89  R.J. BURLEY         CORRECT ALGORITHM INITIALIZATION TO           
C                                RESET MRKCLM AND MRKSTR VALUES FLAGGED        
C                                AS UNIDENTIFIED BY TRIPLT ALGORITHM,          
C                                IN THE EVENT THAT DOUBLT IS BEING USED        
C                                AS A DEFAULT BACKUP FOR TRIPLT.               
C                                MAPCLM VALUES FOR MATCHED CLUMPS SHOULD       
C                                BE SET TO STAR INDEX OF BEST MATCH.           
C  11/ 6/89  R.J. BURLEY         ADD LBLCLM GESS ARRAY                         
C-----------------------------------------------------------------------       
C   METHOD:                                                                    
C   CONVERT PAIRWISE ANGULAR SEPARATION TOLERANCE TO RADIANS                   
C   RESET ALL CLUMP ID FLAGS TO 0 TO INDICATE UNIDENTIFIED AND                 
C    RESET ALL CLUMP AND OBSERVATION STATUS FLAGS TO ZERO IF THE               
C    FLAG INDICATES IT WAS DROPPED DUE TO UNIDENTIFICATION BY TRIPLT.          
C                                                                              
C                                                                              
C!  TEST ALL COMBINATIONS OF CLUMP PAIRS                                       
C   DO FOR ICLMA = 1 TO (NUMCLM-1)                                             
C     DO FOR ICLMB = (ICLMA+1),NUMCLM                                          
C                                                                              
C       IF ((BOTH CLUMPS HAVE 1 OR MORE REFERENCE STAR MATCHES).AND.           
C           (AT LEAST 1 OF THE CLUMPS IS NOT YET IDENTIFIED)) THEN             
C         INCREMENT COUNT OF VALID PAIRS                                       
C         USE DANGLE TO COMPUTE THE CLUMP SEPARATION ANGLE                     
C                                                                              
C                                                                              
C         SET THE NUMBER OF MATCHES FOR THIS CLUMP PAIR TO ZERO                
C         CURRENT SMALLEST DIFFERENCE = SOME MAXIMUM HIGH VALUE                
C                                                                              
C!        TEST ALL COMBINATIONS OF REFERENCE STAR PAIRS                        
C         DO FOR ALL REF STARS MATCHED TO ICLMA BY DIRECT MATCH                
C           DO FOR ALL REF STARS MATCHED TO ICLMB BY DIRECT MATCH              
C             USE DANGLE TO COMPUTE THE STAR SEPARATION ANGLE                  
C             ANGDIF = ABSOLUTE VALUE OF THE DIFFERENCE BETWEEN THE            
C                      CLUMP SEPARATION ANGLE AND STAR SEPARATION ANGLE        
C             IF ((ANGDIF <= CURRENT SMALLEST DIFFERENCE).AND.                 
C                 (ANGDIF < PAIRWISE SEPARATION TOLERANCE )) THEN              
C               IF (ANGDIF = CURRENT SMALLEST DIFFERENCE) THEN                 
C                 INCREMENT NUMBER OF MATCHES                                  
C               ELSE                                                           
C                 SET NUMBER OF MATCHES TO 1                                   
C               ENDIF                                                          
C               SAVE ANGDIF AS CURRENT SMALLEST AND SAVE REF STAR #'S          
C             ENDIF                                                            
C           ENDDO FOR                                                          
C         ENDDO FOR                                                            
C                                                                              
C         IF (MATCHES = 0) THEN                                                
C           MARK CLUMPS AS UN+IDENTIFIED IF NOT YET IDENTIFIED BY              
C            SETTING IDFCLM FLAG TO 0                                          
C         ELSE IF (MATCHES > 1) THEN                                           
C           MARK CLUMPS AS QUESTIONABLE IF NOT YET IDENTIFIED BY               
C            SETTING IDFCLM FLAG TO 1                                          
C         ELSE IF (MATCHES = 1) THEN                                           
C           MARK CLUMPS AS IDENTIFIED WITH REFERENCE STARS BY                  
C            SETTING IDFCLM FLAG TO 2                                          
C           DELETE ANY OTHER REFERENCE STARS MATCHED TO CLUMPS                 
C         ENDIF                                                                
C       ENDIF                                                                  
C     ENDDO FOR                                                                
C   ENDDO FOR                                                                  
C                                                                              
C   DO FOR ALL CLUMPS                                                          
C     IF (CLUMP HAS NOT BEEN IDENTIFIED BY THE DOUBLT ALGORITHM) THEN          
C        SET ITS MRKCLM FLAG TO 15 TO INDICATE SO.                             
C   ENDDO FOR                                                                  
C   DO FOR ALL OBSERVATIONS                                                    
C     IF (THE CLUMP THAT OBSERVATION IS IN IS UNIDENTIFIED) THEN               
C       SET ITS MRKSTR STATUS FLAG TO 17 TO INDICATE DROPPED DUE               
C        TO UNIDENTIFICATION IN DOUBLT MATCHING ALGORITHM.                     
C     ENDIF                                                                    
C   ENDDO FOR                                                                  
C                                                                              
C   RETURN                                                                     
C-----------------------------------------------------------------------       
C                                                                              
C     * DEFINE PARAMETER VARIABLES                                             
      REAL*8    PANGTL                                                         
C                                                                              
      REAL*4    GCICLM(3,*) , SKYCLM(10,3,*)                                   
C                                                                              
      INTEGER*4 NUMCLM      , NUMSTR      , LBLCLM(*), KLMSTR(*)               
      INTEGER*4 MRKCLM(*)   , IDFCLM(*)   , NRFCLM(*)                          
      INTEGER*4 MAPCLM(10,*), MRKSTR(*)   , NUMDUB                             
C                                                                              
C     * DECLARE LOCAL VARIABLES                                                
      REAL*8    RANGTL   , SEPANG   , REFANG   , ANGDIF   , TOLER              
      REAL*8    E2CLMA(3), E2CLMB(3), E2REF1(3), E2REF2(3)                     
      REAL*8    BESTA    , RMAX     , DANGLE                                   
      INTEGER*4 ICLMA    , ICLMB    , IMATCH   , IREF1    , IREF2              
      INTEGER*4 IBEST1   , IBEST2   , ISTR     , I                             
      LOGICAL*4 L_MTCH   , L_NOID
      DATA      TOLER /0.99D0/                                                 
      DATA      RMAX  /999999.0D0/                                             
C                                                                              
C                                                                              
C                                      INITIALIZE ROUTINE                      
      IF (LEVDBG(7) .NE. 0) WRITE (LUDBUG,1000)                                
C                                                                              
C                                      INITIALIZE ALGORITHM                    
      NUMDUB = 0                                                               
C      RANGTL = PANGTL * DTR                                                    
      RANGTL = PANGTL / RADEG                                                  
      DO 110 ICLMA = 1,NUMCLM                                                  
        IDFCLM(ICLMA) = 0                                                      
        IF (MRKCLM(ICLMA) .EQ. 15) THEN                                        
          MRKCLM(ICLMA) = 0                                                    
          DO 100 ISTR = 1,NUMSTR                                               
            IF ((KLMSTR(ISTR) .EQ. ICLMA).AND.                                 
     1          (MRKSTR(ISTR) .EQ. 17)) MRKSTR(ISTR) = 0                       
100       CONTINUE                                                             
        ENDIF                                                                  
110   CONTINUE                                                                 
C                                                                              
C                                      PERFORM PAIRWISE CLUMP MATCHING         
      DO 500 ICLMA = 1,(NUMCLM-1)                                              
        DO 400 ICLMB = (ICLMA+1),NUMCLM                                        
C                                                                              
C                                      TEST IF CLUMPS MAKE VALID PAIR          
          IF ((NRFCLM(ICLMA) .GE. 1).AND.(NRFCLM(ICLMB) .GE. 1)) THEN          
            L_MTCH = .TRUE.                                                    
          ELSE                                                                 
            L_MTCH = .FALSE.                                                   
          ENDIF                                                                
          IF ((IDFCLM(ICLMA) .EQ. 0).OR.(IDFCLM(ICLMB) .EQ. 0)) THEN           
            L_NOID = .TRUE.                                                    
          ELSE                                                                 
            L_NOID = .FALSE.                                                   
          ENDIF                                                                
          IF (L_MTCH .AND. L_NOID) THEN                                        
C                                                                              
C                                      VALID CLUMP PAIR                        
            NUMDUB = NUMDUB + 1                                                
C                                                                              
C                                      COMPUTE ANGLE BETWEEN CLUMPS            
            E2CLMA(1) = DBLE(GCICLM(1,ICLMA))                                  
            E2CLMA(2) = DBLE(GCICLM(2,ICLMA))                                  
            E2CLMA(3) = DBLE(GCICLM(3,ICLMA))                                  
            E2CLMB(1) = DBLE(GCICLM(1,ICLMB))                                  
            E2CLMB(2) = DBLE(GCICLM(2,ICLMB))                                  
            E2CLMB(3) = DBLE(GCICLM(3,ICLMB))                                  
            SEPANG    = DANGLE(E2CLMA, E2CLMB, TOLER, ANGDIF)                  
C                                                                              
C                                      PERFORM PAIRWISE STAR MATCHING          
            IMATCH = 0                                                         
            BESTA  = RMAX                                                      
            DO 300 IREF1 = 1,NRFCLM(ICLMA)                                     
              DO 200 IREF2 = 1,NRFCLM(ICLMB)                                   
C                                                                              
C                                      COMPUTE ANGLE BETWEEN STARS             
                E2REF1(1) = DBLE(SKYCLM(IREF1,1,ICLMA))                        
                E2REF1(2) = DBLE(SKYCLM(IREF1,2,ICLMA))                        
                E2REF1(3) = DBLE(SKYCLM(IREF1,3,ICLMA))                        
                E2REF2(1) = DBLE(SKYCLM(IREF2,1,ICLMB))                        
                E2REF2(2) = DBLE(SKYCLM(IREF2,2,ICLMB))                        
                E2REF2(3) = DBLE(SKYCLM(IREF2,3,ICLMB))                        
                REFANG    = DANGLE(E2REF1, E2REF2, TOLER, ANGDIF)              
C                                                                              
C                                      COMPUTE DIFFERENCE BETWEEN SUMS         
                ANGDIF    = DABS(REFANG - SEPANG)                              
C                                                                              
C                                      SAVE BEST FIT WITHIN RANGTL LIMIT       
                IF ((ANGDIF .LE. BESTA).AND.(ANGDIF .LT. RANGTL)) THEN         
                  IF (ANGDIF .EQ. BESTA) THEN                                  
                    IMATCH = IMATCH + 1                                        
                  ELSE                                                         
                    IMATCH = 1                                                 
                  ENDIF                                                        
                  BESTA  = ANGDIF                                              
                  IBEST1 = IREF1                                               
                  IBEST2 = IREF2                                               
                ENDIF                                                          
C                                                                              
C                                      INTERMEDIATE DEBUG                      
                IF (LEVDBG(7) .GE. 4) THEN                                     
                  WRITE (LUDBUG,4000) LBLCLM(ICLMA), LBLCLM(ICLMB),            
     1                                SEPANG                                   
                  WRITE (LUDBUG,4010) IREF1, IREF2, MAPCLM(IREF1,ICLMA),       
     1                                MAPCLM(IREF2,ICLMB), REFANG              
                  WRITE (LUDBUG,4020) ANGDIF, BESTA, IMATCH                    
                ENDIF                                                          
C                                                                              
C                                                                              
200           CONTINUE                                                         
300         CONTINUE                                                           
C                                                                              
C                                                                              
            IF (IMATCH .EQ. 0) THEN                                            
C                                      MARK CLUMPS AS UNIDENTIFIED             
              IF (IDFCLM(ICLMA) .NE. 2) IDFCLM(ICLMA) = 0                      
              IF (IDFCLM(ICLMB) .NE. 2) IDFCLM(ICLMB) = 0                      
            ELSE IF (IMATCH .GT. 1) THEN                                       
C                                      MARK CLUMPS AS QUESTIONABLE             
              IF (IDFCLM(ICLMA) .NE. 2) IDFCLM(ICLMA) = 1                      
              IF (IDFCLM(ICLMB) .NE. 2) IDFCLM(ICLMB) = 1                      
            ELSE                                                               
C                                                                              
C                                      MARK CLUMPS AS IDENTIFIED               
              IDFCLM(ICLMA) = 2                                                
              IDFCLM(ICLMB) = 2                                                
C                                                                              
C                                      CLUMPS NOW HAVE ONLY 1 MATCH            
              NRFCLM(ICLMA) = 1                                                
              NRFCLM(ICLMB) = 1                                                
C                                                                              
C                                      SAVE BEST MATCH AT FRONT OF LIST        
              MAPCLM(1,ICLMA)   = MAPCLM(IBEST1,ICLMA)                         
              MAPCLM(1,ICLMB)   = MAPCLM(IBEST2,ICLMB)                         
C                                                                              
C                                      SAVE BEST MATCH STAR POSITIONS          
              SKYCLM(1,1,ICLMA) = SKYCLM(IBEST1,1,ICLMA)                       
              SKYCLM(1,2,ICLMA) = SKYCLM(IBEST1,2,ICLMA)                       
              SKYCLM(1,3,ICLMA) = SKYCLM(IBEST1,3,ICLMA)                       
              SKYCLM(1,1,ICLMB) = SKYCLM(IBEST2,1,ICLMB)                       
              SKYCLM(1,2,ICLMB) = SKYCLM(IBEST2,2,ICLMB)                       
              SKYCLM(1,3,ICLMB) = SKYCLM(IBEST2,3,ICLMB)                       
            ENDIF                                                              
          ENDIF                                                                
400     CONTINUE                                                               
500   CONTINUE                                                                 
C                                                                              
C                                      FLAG UNIDENTIFIED CLUMPS                
      DO 600 ICLMA=1,NUMCLM                                                    
        IF ((MRKCLM(ICLMA) .EQ. 0).AND.                                        
     1      (IDFCLM(ICLMA) .EQ. 0)) MRKCLM(ICLMA) = 15                         
600   CONTINUE                                                                 
C                                                                              
C                                      FLAG UNIDENTIFIED OBSERVATIONS          
      DO 700 ISTR = 1,NUMSTR                                                   
        IF (MRKCLM(KLMSTR(ISTR)) .EQ. 15) MRKSTR(ISTR) = 17                    
700   CONTINUE                                                                 
C                                                                              
C                                                                              
C                                      OUTGOING DEBUG                          
      IF (LEVDBG(7) .GE. 3) THEN                                               
        WRITE (LUDBUG,5000)                                                    
        DO 950 ICLMA=1,NUMCLM                                                  
          WRITE (LUDBUG,5010) LBLCLM(ICLMA),MRKCLM(ICLMA),IDFCLM(ICLMA)        
          DO 900 IMATCH=1,NRFCLM(ICLMA)                                        
            WRITE (LUDBUG,5020) IMATCH, MAPCLM(IMATCH,ICLMA),                  
     1                          (SKYCLM(IMATCH,I,ICLMA),I=1,3)                 
900       CONTINUE                                                             
950     CONTINUE                                                               
      ENDIF                                                                    
C                                                                              
C                                      NORMAL TERMINATION                      
      IF (LEVDBG(7) .NE. 0) WRITE (LUDBUG,2000)                                
      RETURN                                                                   
C                                                                              
C                                      FORMAT SECTION                          
1000  FORMAT(' *** ENTER DOUBLT ***')                                          
2000  FORMAT(' *** EXIT  DOUBLT ***')                                          
4000  FORMAT(' INTERMEDIATE DEBUG: '/,                                         
     1       4X,'CLUMP PAIR=',I8,2X,I8,' SEPARATION=',D14.8)                   
4010  FORMAT(4X,'STAR PAIR =',I8,2X,I8,' CATALOG INDEX=',I8,2X,I8,             
     1       ' SEPARATION=',D14.8)                                             
4020  FORMAT(4X,'ANGDIF=',D14.8,' BEST YET=',D14.8,' MATCH #=',I8)             
5000  FORMAT(' STAR CLUMP TABLE AFTER DOUBLT ALGORITHM')                       
5010  FORMAT(' CLUMP=',I6,' STATUS=',I6,' IDFLAG=',I6)                         
5020  FORMAT(10X,'MATCH#',I4,' CATINDEX=',I8,' POSITION=',3(F13.6,1X))         
      END                                                                      
