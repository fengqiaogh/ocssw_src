      SUBROUTINE TRIPLT                                                         
     I          (TANGTL, TMINCO, NUMCLM, LBLCLM, GCICLM, NOBCLM,                
     I           MRKCLM, NUMSTR, KLMSTR,                                        
     O           MRKSTR, SKYCLM, NRFCLM, MAPCLM,                                
     O           IDFCLM, NUMTRP)                                                
C-----------------------------------------------------------------------        
C   MODULE NAME:  STTRIPLT                                                      
C                                                                               
C                                                                               
C   PURPOSE: TO IDENTIFY STAR CLUMPS USING TRIPLET STAR MATCHING.               
C                                                                               
C                                                                               
C   ARGUMENT LIST:                                                              
C   ARGUMENT I/O TYPE DIMENS DESCRIPTION                                        
C   -------- --- ---- ------ -----------                                        
C   TANGTL   I   R*8         MAX ANGULAR SEPARATION FOR TRIPLET MATCH           
C   TMINCO   I   R*8         MINIMUM COLINEARITY ANGLE FOR VALID TRIPLET        
C   NUMCLM   I   I*4         NUMBER OF CLUMPS                                   
C   LBLCLM   I   I*4         CLUMP LABELS                                       
C   GCICLM   I   R*4   3,*   AVERAGE POSITION VEC (GCI) FOR EACH CLUMP          
C   NOBCLM   I   I*4    *    NUMBER OF OBSERVATIONS PER CLUMP                   
C   MRKCLM   I O I*4    *    STATUS FLAG FOR EACH CLUMP                         
C   NUMSTR   I   I*4         NUMBER OF OBSERVATIONS                             
C   KLMSTR   I   I*4    *    CLUMP NUMBER FOR EACH OBSERVATION                  
C   MRKSTR   I O I*4    *    STATUS FLAG FOR EACH OBSERVATION                   
C   SKYCLM   I O R*4 10,3,*  LIST OF CORRECTED REFERENCE STAR POSITIONS         
C   NRFCLM   I O I*4    *    NUMBER OF REFERENCE STARS MATCHED TO CLUMP         
C   MAPCLM   I O I*4  10,*   LIST OF REFERENCE STAR #'S MATCHED TO CLUMP        
C   IDFCLM     O I*4    *    IDENTIFICATION FLAG FOR EACH CLUMP                 
C   NUMTRP     O I*4         NUMBER OF VALID CLUMP TRIPLETS CHECKED             
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
C   STCOLINE - TEST COLINEAR ACCEPTABILITY OF CLUMP TRIPLET                     
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
C   UARS FDSS SPECS 3.1.1.5 (FUNCTION 3)                                        
C                                                                               
C                                                                               
C   DEVELOPMENT HISTORY:                                                        
C   DATE     AUTHOR              DESCRIPTION                                    
C   -------- ------              -----------                                    
C   8/ 3/88  R.J. BURLEY         DESIGN                                         
C   5/18/89  R.J. BURLEY         CODED                                          
C  10/17/89  R.J. BURLEY         ADD FLAGGING OF UNIDENTIFIED CLUMPS AND        
C                                CORRECT FLAGGING OF OBSERVATIONS.              
C  10/19/89  R.J. BURLEY         PUT CATALOG STAR INDEX NUMBER IN MAPCLM        
C                                ARRAY INSTEAD OF BEST MATCH NUMBER.            
C  10/24/89  R.J. BURLEY         CONVERT MINIMUM COLINEARITY ANGLE TO           
C                                TO RADIANS.                                    
C  11/ 2/89  R.J. BURLEY         MODIFY DEBUG                                   
C  11/ 6/89  R.J. BURLEY         ADD LBLCLM GESS ARRAY                          
C-----------------------------------------------------------------------        
C   METHOD:                                                                     
C   CONVERT TRIPLET ANGULAR SEPARATION TOLERANCE TO RADIANS                     
C   CONVERT MINIMUM COLINEARITY ANGLE TO RADIANS                                
C   SET NUMTRP TO ZERO                                                          
C   RESET ALL CLUMP ID FLAGS AS UNIDENTIFIED IN IDFCLM ARRAY                    
C                                                                               
C   DO FOR ICLMPA = 1 TO (NUMCLM-2)                                             
C     DO FOR ICLMPB = (ICLMPA+1) TO (NUMCLM-1)                                  
C       DO FOR ICLMPC = (ICLMPB+1) TO NUMCLM                                    
C                                                                               
C!        * TEST VALIDITY OF CLUMP TRIPLET                                      
C         CALL COLINE TO TEST COLINEAR ACCEPTABILITY OF CLUMP TRIPLET           
C         IF (ALL HAVE >=1 REF STAR).AND.(ALL HAVE MRKCLM FLAG = 0).AND.        
C           (1 OR MORE IS UNIDENTIFIED).AND.(COLINEARITY > TIMINCO) THEN        
C             INCREMENT NUMBER OF VALID TRIPLETS COUNT                          
C             COMPUTE SUM OF ANGULAR SEPARATIONS OF CLUMP TRIPLET               
C             SET NUMBER OF TRIPLET MATCHES TO ZERO                             
C             CURRENT SMALLEST DIFFERENCE = SOME MAXIMUM VALUE                  
C                                                                               
C!            * TEST ALL REFERENCE STAR COMBINATIONS                            
C             DO FOR ALL REFERENCE STARS MATCHED TO CLUMPA                      
C               DO FOR ALL REFERENCE STARS MATCHED TO CLUMPB                    
C                 DO FOR ALL REFERENCE STARS MATCHED TO CLUMPC                  
C                   COMPUTE SUM OF ANG SEPS OF REF STAR TRIPLET                 
C                   COMPUTE DIFFERENCE BETWEEN SUM OF CLUMP TRIPLET             
C                    ANGLES AND STAR TRIPLET ANGLES                             
C                   IF ((DIFFERENCE <= CURRENT SMALLEST DIFFERENCE).AND.        
C                       (DIFFERENCE < TOLERANCE)) THEN                          
C                     IF (DIFFERENCE = CURRENT SMALLEST DIFFERENCE) THEN        
C                       INCREMENT NUMBER OF MATCHES                             
C                     ELSE                                                      
C                       NUMBER OF MATCHES = 1                                   
C                     ENDIF                                                     
C                     SET CURRENT SMALLEST DIFFERENCE = DIFFERENCE              
C                     SAVE THE REF STAR #'S WHICH MADE THIS TRIPLET             
C                   ENDIF                                                       
C                 ENDDO FOR                                                     
C               ENDDO FOR                                                       
C             ENDDO FOR                                                         
C                                                                               
C             IF (NUMBER OF MATCHES = 0) THEN                                   
C               MARK CLUMPS AS UNIDENTIFIED IF NOT YET IDENTIFIED               
C             ELSE IF (NUMBER OF MATCHES > 1) THEN                              
C               MARK CLUMPS AS QUESTIONABLE IF NOT YET IDENTIFIED               
C             ELSE                                                              
C               MARK CLUMPS AS IDENTIFIED, SAVE REFERENCE STAR NUMBERS          
C               WIPE OUT ANY OTHER REFERENCE STARS MATCHED TO CLUMP             
C             ENDIF                                                             
C           ENDIF                                                               
C         ENDIF                                                                 
C                                                                               
C       ENDDO FOR                                                               
C     ENDDO FOR                                                                 
C   ENDDO FOR                                                                   
C                                                                               
C   DO FOR ALL CLUMPS                                                           
C     IF (CLUMP HAS NOT BEEN IDENTIFIED BY TRIPLT ALGORITHM) THEN               
C       SET ITS MRKCLM FLAG TO 15.                                              
C   ENDDO FOR                                                                   
C   DO FOR ALL OBSERVATIONS                                                     
C     SET MRKSTR FLAG TO 17 IF IT IS IN A CLUMP MARKED AS UNIDENTIFIED          
C   ENDDO FOR                                                                   
C   RETURN                                                                      
C-----------------------------------------------------------------------        
C                                                                               
C     * DEFINE PARAMETER VARIABLES                                              
      REAL*8    TANGTL      , TMINCO                                            
C                                                                               
      REAL*4    GCICLM(3,*) , SKYCLM(10,3,*)                                    
C                                                                               
      INTEGER*4 NUMCLM      , LBLCLM(*)   , NOBCLM(*)                           
      INTEGER*4 MRKCLM(*)   , IDFCLM(*)   , NUMSTR       , KLMSTR(*)            
      INTEGER*4 MRKSTR(*)   , NRFCLM(*)   , MAPCLM(10,*) , NUMTRP               
C                                                                               
C                                                                               
C     * DECLARE LOCAL VARIABLES                                                 
      REAL*8    RANGTL   , ANGDIF   , TOLER    , BESTA , RMAX  , COANGL         
      REAL*8    CANGAB   , CANGAC   , CANGBC   , CANGLE, DANGLE, RMINCO         
      REAL*8    RANG12   , RANG13   , RANG23   , RANGLE                         
      REAL*8    E2CLMA(3), E2CLMB(3), E2CLMC(3)                                 
      REAL*8    E2REF1(3), E2REF2(3), E2REF3(3)                                 
      INTEGER*4 ICLMA    , ICLMB    , ICLMC    , IREF1  , IREF2 , IREF3         
      INTEGER*4 IMATCH   , IBEST1   , IBEST2   , IBEST3 , ISTR                  
      LOGICAL*4 L_GOOD   , L_MTCH   , L_NOID                                    
      DATA      TOLER /0.99D0/                                                  
      DATA      RMAX  /999999.0D0/                                              
C                                                                               
C                                                                               
C                                      INITIALIZE ROUTINE                       
      IF (LEVDBG(7) .NE. 0) WRITE (LUDBUG,1000)                                 
      IF (LEVDBG(7) .GE. 1) WRITE (LUDBUG,3999) TANGTL,TMINCO                   
C                                                                               
C                                      INITIALIZE ALGORITHM                     
      NUMTRP = 0                                                                
C      RANGTL = TANGTL * DTR                                                    
C      RMINCO = TMINCO * DTR                                                    
      RANGTL = TANGTL / RADEG                                                   
      RMINCO = TMINCO / RADEG                                                   
      IF (LEVDBG(7) .GE. 1) WRITE (LUDBUG,3999) RANGTL,RMINCO                   
      DO 100 ICLMA = 1,NUMCLM                                                   
        IDFCLM(ICLMA) = 0                                                       
100   CONTINUE                                                                  
C
      ISTR = 1
      DOWHILE(NRFCLM(ISTR).LT.1)
        ISTR = ISTR + 1
      END DO                                                                               
C                                      PERFORM TRIPLET CLUMP MATCHING           
      DO 700 ICLMA = ISTR,(NUMCLM-2)                                               
        DO 600 ICLMB = (ICLMA+1),(NUMCLM-1)                                     
          DO 500 ICLMC = (ICLMB+1),NUMCLM                                       
C                                                                               
C                                      TEST IF CLUMPS MAKE VALID TRIPLET        
            IF ((MRKCLM(ICLMA) .EQ. 0).AND.(MRKCLM(ICLMB) .EQ. 0).AND.          
     1          (MRKCLM(ICLMC) .EQ. 0)) THEN                                    
              L_GOOD = .TRUE.                                                   
            ELSE                                                                
              L_GOOD = .FALSE.                                                  
            ENDIF                                                               
            IF ((NRFCLM(ICLMA) .GE. 1).AND.(NRFCLM(ICLMB) .GE. 1).AND.          
     1          (NRFCLM(ICLMC) .GE. 1)) THEN                                    
              L_MTCH = .TRUE.                                                   
            ELSE                                                                
              L_MTCH = .FALSE.                                                  
            ENDIF                                                               
            IF ((IDFCLM(ICLMA) .EQ. 0).OR.(IDFCLM(ICLMB) .EQ. 0).OR.            
     1          (IDFCLM(ICLMC) .EQ. 0)) THEN                                    
              L_NOID = .TRUE.                                                   
            ELSE                                                                
              L_NOID   = .FALSE.                                                
            ENDIF                                                               
            IF (L_GOOD .AND. L_MTCH .AND. L_NOID) THEN                       
            CALL COLINE (ICLMA, ICLMB, ICLMC, LBLCLM, GCICLM, COANGL)           
            IF (COANGL .GT. RMINCO) THEN                                      
C                                                                               
C                                      VALID CLUMP TRIPLET                      
              NUMTRP = NUMTRP + 1                                               
C                                                                               
C                                      FIND SUM OF ANGULAR SEPARATIONS          
              E2CLMA(1) = DBLE(GCICLM(1,ICLMA))                                 
              E2CLMA(2) = DBLE(GCICLM(2,ICLMA))                                 
              E2CLMA(3) = DBLE(GCICLM(3,ICLMA))                                 
              E2CLMB(1) = DBLE(GCICLM(1,ICLMB))                                 
              E2CLMB(2) = DBLE(GCICLM(2,ICLMB))                                 
              E2CLMB(3) = DBLE(GCICLM(3,ICLMB))                                 
              E2CLMC(1) = DBLE(GCICLM(1,ICLMC))                                 
              E2CLMC(2) = DBLE(GCICLM(2,ICLMC))                                 
              E2CLMC(3) = DBLE(GCICLM(3,ICLMC))                                 
              CANGAB    = DANGLE(E2CLMA, E2CLMB, TOLER, ANGDIF)                 
              CANGBC    = DANGLE(E2CLMB, E2CLMC, TOLER, ANGDIF)                 
              CANGAC    = DANGLE(E2CLMA, E2CLMC, TOLER, ANGDIF)                 
              CANGLE    = CANGAB + CANGBC + CANGAC                              
C                                                                               
C                                      PERFORM TRIPLET STAR MATCHING            
              IMATCH = 0                                                        
              BESTA  = RMAX                                                     
              DO 400 IREF1 = 1,NRFCLM(ICLMA)                                    
                DO 300 IREF2 = 1,NRFCLM(ICLMB)                                  
                  DO 200 IREF3 = 1,NRFCLM(ICLMC)                                
C                                                                               
C                                      COMPUTE SUM OF ANGLES BETWEEN            
C                                      REFERENCE STARS                          
                    E2REF1(1) = DBLE(SKYCLM(IREF1,1,ICLMA))                     
                    E2REF1(2) = DBLE(SKYCLM(IREF1,2,ICLMA))                     
                    E2REF1(3) = DBLE(SKYCLM(IREF1,3,ICLMA))                     
                    E2REF2(1) = DBLE(SKYCLM(IREF2,1,ICLMB))                     
                    E2REF2(2) = DBLE(SKYCLM(IREF2,2,ICLMB))                     
                    E2REF2(3) = DBLE(SKYCLM(IREF2,3,ICLMB))                     
                    E2REF3(1) = DBLE(SKYCLM(IREF3,1,ICLMC))                     
                    E2REF3(2) = DBLE(SKYCLM(IREF3,2,ICLMC))                     
                    E2REF3(3) = DBLE(SKYCLM(IREF3,3,ICLMC))                     
                    RANG12    = DANGLE(E2REF1, E2REF2, TOLER, ANGDIF)           
                    RANG13    = DANGLE(E2REF1, E2REF3, TOLER, ANGDIF)           
                    RANG23    = DANGLE(E2REF2, E2REF3, TOLER, ANGDIF)           
                    RANGLE    = RANG12 + RANG13 + RANG23                        
C                                                                               
C                                      COMPUTE DIFFERENCE BETWEEN SUMS          
C                    ANGDIF = DABS(CANGLE - RANGLE)                           
                    ANGDIF = DABS(CANGAB - RANG12)                              
     *                     + DABS(CANGAC - RANG13)                              
     *                     + DABS(CANGBC - RANG23)                              
C                                                                               
C                                      SAVE BEST FIT WITHIN RANGTL LIMIT        
                    IF ((ANGDIF.LE.BESTA).AND.(ANGDIF.LT.RANGTL)) THEN          
                       IF (ANGDIF .EQ. BESTA) THEN                              
                         IMATCH = IMATCH + 1                                    
                       ELSE                                                     
                         IMATCH = 1                                             
                       ENDIF                                                    
                       BESTA  = ANGDIF                                          
                       IBEST1 = IREF1                                           
                       IBEST2 = IREF2                                           
                       IBEST3 = IREF3                                           
                    ENDIF                                                       
C                                                                               
C                                      INTERMEDIATE DEBUG                       
                    IF (LEVDBG(7) .GE. 4) THEN                                  
                      WRITE (LUDBUG,4000) LBLCLM(ICLMA),                        
     1                     LBLCLM(ICLMB), LBLCLM(ICLMC), CANGLE                 
                      WRITE (LUDBUG,4010) IREF1, IREF2, IREF3,                  
     1                     MAPCLM(IREF1,ICLMA), MAPCLM(IREF2,ICLMB),            
     2                     MAPCLM(IREF3,ICLMC), RANGLE                          
                      WRITE (LUDBUG,4020) ANGDIF, BESTA, IMATCH                 
                    ENDIF                                                       
C                                                                               
C                                                                               
200               CONTINUE                                                      
300             CONTINUE                                                        
400           CONTINUE                                                          
C                                                                               
C                                                                               
              IF (IMATCH .EQ. 0) THEN                                           
C                                      MARK CLUMPS AS UNIDENTIFIED              
                IF (IDFCLM(ICLMA) .NE. 2) IDFCLM(ICLMA) = 0                     
                IF (IDFCLM(ICLMB) .NE. 2) IDFCLM(ICLMB) = 0                     
                IF (IDFCLM(ICLMC) .NE. 2) IDFCLM(ICLMC) = 0                     
              ELSE IF (IMATCH .GT. 1) THEN                                      
C                                      MARK CLUMPS AS QUESTIONABLE              
                IF (IDFCLM(ICLMA) .NE. 2) IDFCLM(ICLMA) = 1                     
                IF (IDFCLM(ICLMB) .NE. 2) IDFCLM(ICLMB) = 1                     
                IF (IDFCLM(ICLMC) .NE. 2) IDFCLM(ICLMB) = 1                     
              ELSE                                                              
C                                                                               
C                                      MARK CLUMPS AS IDENTIFIED                
                IDFCLM(ICLMA) = 2                                               
                IDFCLM(ICLMB) = 2                                               
                IDFCLM(ICLMC) = 2                                               
C                                                                               
C                                      CLUMPS NOW HAVE ONLY 1 MATCH             
                NRFCLM(ICLMA) = 1                                               
                NRFCLM(ICLMB) = 1                                               
                NRFCLM(ICLMC) = 1                                               
C                                                                               
C                                      SAVE BEST MATCH AT FRONT OF LIST         
                MAPCLM(1,ICLMA)   = MAPCLM(IBEST1,ICLMA)                        
                MAPCLM(1,ICLMB)   = MAPCLM(IBEST2,ICLMB)                        
                MAPCLM(1,ICLMC)   = MAPCLM(IBEST3,ICLMC)                        
C                                                                               
C                                      SAVE BEST MATCH STAR POSITIONS           
                SKYCLM(1,1,ICLMA) = SKYCLM(IBEST1,1,ICLMA)                      
                SKYCLM(1,2,ICLMA) = SKYCLM(IBEST1,2,ICLMA)                      
                SKYCLM(1,3,ICLMA) = SKYCLM(IBEST1,3,ICLMA)                      
                SKYCLM(1,1,ICLMB) = SKYCLM(IBEST2,1,ICLMB)                      
                SKYCLM(1,2,ICLMB) = SKYCLM(IBEST2,2,ICLMB)                      
                SKYCLM(1,3,ICLMB) = SKYCLM(IBEST2,3,ICLMB)                      
                SKYCLM(1,1,ICLMC) = SKYCLM(IBEST3,1,ICLMC)                      
                SKYCLM(1,2,ICLMC) = SKYCLM(IBEST3,2,ICLMC)                      
                SKYCLM(1,3,ICLMC) = SKYCLM(IBEST3,3,ICLMC)                      
              ENDIF                                                             
            ENDIF                                                               
            ENDIF                                                               
C                                                                               
C                                                                               
500       CONTINUE                                                              
600     CONTINUE                                                                
700   CONTINUE                                                                  
C                                                                               
C                                      FLAG UNIDENTIFIED CLUMPS                 
      DO 800 ICLMA=1,NUMCLM                                                     
        IF ((MRKCLM(ICLMA) .EQ. 0).AND.                                         
     1      (IDFCLM(ICLMA) .EQ. 0)) MRKCLM(ICLMA) = 15                          
800   CONTINUE                                                                  
C                                                                               
C                                      FLAG UNIDENTIFIED OBSERVATIONS           
      DO 900 ISTR = 1,NUMSTR                                                    
        IF (MRKCLM(KLMSTR(ISTR)) .EQ. 15) MRKSTR(ISTR) = 17                     
900   CONTINUE                                                                  
C                                                                               
C                                                                               
C                                      OUTGOING DEBUG                           
      IF (LEVDBG(7) .GE. 3) THEN                                                
        WRITE (LUDBUG,5000)                                                     
        DO 950 ICLMA=1,NUMCLM                                                   
          WRITE (LUDBUG,5010) LBLCLM(ICLMA),MRKCLM(ICLMA),IDFCLM(ICLMA)         
          DO 930 IMATCH=1,NRFCLM(ICLMA)                                         
            WRITE (LUDBUG,5020) IMATCH, MAPCLM(IMATCH,ICLMA),                   
     1                          (SKYCLM(IMATCH,I,ICLMA),I=1,3)                  
930       CONTINUE                                                              
950     CONTINUE                                                                
      ENDIF                                                                     
C                                                                               
C                                      NORMAL TERMINATION                       
      IF (LEVDBG(7) .NE. 0) WRITE (LUDBUG,2000)                                 
      RETURN                                                                    
C                                                                               
C                                      FORMAT SECTION                           
1000  FORMAT(' *** ENTER TRIPLT ***')                                           
2000  FORMAT(' *** EXIT  TRIPLT ***')                                           
3999  FORMAT(' TRIPLET PARAMETERS: TANGTL=',D12.6,' TMINCO=',D12.6)             
4000  FORMAT(' INTERMEDIATE DEBUG: '/,                                          
     1       4X,'CLUMP TRIPLT=',I8,2X,I8,2X,I8,' SEPARATION=',D12.6)            
4010  FORMAT(4X,'STAR TRIPLT =',I8,2X,I8,2X,I8,' CATALOG #S=',                  
     1       I8,2X,I8,2X,I8,' SEPARATION=',D12.6)                               
4020  FORMAT(4X,'ANGDIF=',D12.6,' BEST YET=',D12.6,' MATCH #=',I8)              
5000  FORMAT(' STAR CLUMP TABLE AFTER TRIPLT ALGORITHM')                        
5010  FORMAT(' CLUMP=',I6,' STATUS=',I6,' IDFLAG=',I6)                          
5020  FORMAT(10X,'MATCH#',I4,' CATINDEX=',I8,' POSITION=',3(F13.6,1X))          
      END
