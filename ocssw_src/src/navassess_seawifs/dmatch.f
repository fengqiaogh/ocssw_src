      SUBROUTINE DMATCH (DANGTL, DMAGTL, IFXCLM, NUMCAT,                       
     I                   IDNCAT, DATCAT, NUMSTR, KLMSTR,                       
     I                   NUMCLM, LBLCLM, GCICLM, BRICLM, NOBCLM,               
     O                   MRKCLM, MRKSTR, IDFCLM, NRFCLM,                       
     O                   MAPCLM, SKYCLM, IRCODE)                               
C-----------------------------------------------------------------------       
C   MODULE NAME: STDMATCH                                                      
C                                                                              
C   PURPOSE: TO PERFORM DIRECT MATCH ALGORITHM ON CLUMPS.                      
C                                                                              
C                                                                              
C   ARGUMENT LIST:                                                             
C   ARGUMENT I/O TYPE DIMENS DESCRIPTION                                       
C   -------- --- ---- ------ -----------                                       
C   DANGTL   I   R*8         MAX ANGULAR SEPARATION FOR MATCH (DEGREES)        
C   DMAGTL   I   R*8         MAX MAGNITUDE DIFFERENCE FOR MATCH                
C   IFXCLM   I   I*4         MAX NUMBER OF REF STAR MATCHES :ARRAY DIM         
C   NUMCAT   I   I*4         NUMBER OF REFERENCE STARS IN CATALOG              
C   IDNCAT   I   I*4    *    SKYMAP REFERENCE STAR NUMBERS IN CATALOG          
C   DATCAT   I   R*4   7,*   CATALOG DATA (1=X AXIS,2=Y AXIS,3=Z AXIS,         
C                                          4=MAGNITUDE,5=MOTION,               
C                                          6=QUALITY,7=COLOR)                  
C   NUMSTR   I   I*4         NUMBER OF OBSERVATIONS                            
C   KLMSTR   I   I*4    *    CLUMP NUMBER FOR EACH OBSERVATION                 
C   NUMCLM   I   I*4         NUMBER OF CLUMPS PER TRACKER                      
C   LBLCLM   I   I*4    *    CLUMP LABELS                                      
C   GCICLM   I   R*4   3,*   AVERAGE POSITION VEC (GCI) FOR EACH CLUMP         
C   BRICLM   I   R*4    *    AVERAGE MAGNITUDE FOR EACH CLUMP                  
C   NOBCLM   I   I*4    *    NUMBER OF OBSERVATIONS IN EACH CLUMP              
C   MRKCLM   I O I*4    *    STATUS FLAG FOR EACH CLUMP                        
C   MRKSTR   I O I*4    *    STATUS FLAG FOR EACH OBSERVATION                  
C   IDFCLM     O I*4    *    CLUMP IDENTIFICATION FLAG                         
C   NRFCLM     O I*4    *    NUMBER OF REFERENCE STARS MATCHED TO CLUMP        
C   MAPCLM     O I*4  10,*   LIST OF REF STAR NUMBERS MATCHED TO CLUMP         
C   SKYCLM     O R*4 10,3,*  LIST OF CORRECTED REFERENCE STAR POSITIONS        
C   IRCODE     O I*4         ERROR FLAG (0=NO ERROR, 1=ERROR)                  
C                                                                              
C                                                                              
C   COMMON BLOCK VARIABLES USED:                                               
C   COMMON   VAR    I/O   VAR    I/O   VAR    I/O   VAR    I/O                 
C   ------   ---    ---   ---    ---   ---    ---   ---    ---                 
C   CMDEBG   LEVDBG I     LUDBUG I                                             
C   CMSMSG   IMSGNM I     IVARLN I     IDSTFG I     IRC    I                   
C            C$VDAT I     C$SBID I                                             
C   CMCONV   DTR    I                                                          
C                                                                              
C       ++INCLUDE STCMDEBG                                                      
C       ++INCLUDE STCMSMSG                                                      
C       ++INCLUDE AECMCONV                                                      
      INTEGER*4 LEVDBG(8),LUDBUG
      COMMON /CMDEBG/LEVDBG,LUDBUG                                             
      REAL*8 PI,RADEG,RE,REM,F,OMF2,OMEGAE
      COMMON /GCONST/PI,RADEG,RE,REM,F,OMF2,OMEGAE
C                                                                              
C   EXTERNAL FILES REFERENCED:                                                 l
C   FILENAME      OPERATION   FORTRAN UNIT ID                                  
C   --------      ---------   ---------------                                  
C   NONE                                                                       
C                                                                              
C   EXTERNAL REFERENCES:                                                       
C   --------------------------------------------------------------------       
C   STDIRINI - INITIALIZE DIRECT MATCH ALGORITHM AND DATA                      
C   UTDANGLE - FUNCTION TO COMPUTE ANGULAR SEPARATION BETWEEN 2 VECTORS        
C   UTMSG    - OUTPUT AN ERROR MESSAGE                                         
C                                                                              
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
C   UARS FDSS SPECS PAGE 3.1.1.5-17 F3-1A TO F3-1B                             
C                                                                              
C   DEVELOPMENT HISTORY:                                                       
C   DATE      AUTHOR              DESCRIPTION                                  
C   --------  ------              -----------                                  
C   8/2/88    R.J. BURLEY         DESIGN                                       
C   5/8/89    R.J. BURLEY         CODED                                        
C   10/10/89  R.J. BURLEY         ROUTINE NAME CHANGED FROM DIRECT TO          
C                                 DMATCH BECAUSE DIRECT IS A RESERVED          
C                                 NAME IN GESS.                                
C   10/12/89  R.J. BURLEY         IMPROVE DEBUG                                
C   11/ 2/89  R.J. BURLEY         CORRECT IDFCLM VALUES                        
C   11/ 6/89  R.J. BURLEY         ADD LBLCLM GESS ARRAY                        
C-----------------------------------------------------------------------       
C   METHOD:                                                                    
C     CALL DIRINI TO INITIALIZE CLUMP AND OBSERVATION ID FIELDS                
C                                                                              
C!    * MATCH CLUMPS WITH REFERENCE STARS                                      
C     DO FOR ALL CLUMPS                                                        
C       IF (MRKCLM FLAG INDICATES THAT CLUMP IS GOOD) THEN                     
C         DO FOR ALL REFERENCE STARS IN CATALOG  (1 TO NUMCAT)                 
C           USE DANGLE TO COMPUTE ANGULAR SEPARATION BETWEEN AVG CLUMP         
C             POSITION AND REFERENCE STAR POSITION                             
C           COMPUTE ABSOLUTE VALUE OF MAGNITUDE DIFFERENCE BETWEEN CLUMP       
C             AND REFERENCE STAR                                               
C           IF ((ANGULAR SEPARATION < DANGTL).AND.                             
C               (MAGNITUDE DIFFERENCE < DMAGTL)) THEN                          
C!            WITHIN TOLERANCE, SAVE REF STAR NUMBER IN LIST OF MATCHES        
C!             FOR CURRENT CLUMP AS FOLLOWS:                                   
C             NRFCLM(CLUMP#) = NRFCLM(CLUMP#)+1                                
C             IF (NRFCLM(CLUMP#) > IFXCLM) THEN                                
C               IRCODE = 1   OVERFILL OF MATCH ARRAY SPACE                     
C               ABORT TO ERROR_HANDLE                                          
C             ELSE                                                             
C               MAPCLM(NRFCLM(CLUMP#),CLUMP#)=IDNCAT(REF STAR#)                
C             ENDIF                                                            
C           ENDIF                                                              
C         ENDDO FOR                                                            
C       ENDIF                                                                  
C     ENDDO FOR                                                                
C                                                                              
C     DO FOR ALL CLUMPS (1 TO NUMCLM)                                          
C       IF (# OF REFERENCE STAR MATCHES FOR CLUMP = 0) THEN                    
C         IDFCLM(CLUMP#) = 0    MARK CLUMP AS UNIDENTIFIED                     
C         SET MRKCLM FLAG TO 14 TO INDICATE UNIDENTIFIED                       
C         DO FOR ALL OBSERVATIONS                                              
C           IF (OBSERVATION IS IN CURRENT CLUMP AND THE OBSERVATION            
C               HAS NOT YET BEEN DROPPED) THEN                                 
C             SET FLAG IN MRKSTR ARRAY TO 16 INDICATING THE OBSERVATION        
C               HAS BEEN DROPPED DUE TO NON-IDENTIFICATION                     
C           ENDIF                                                              
C         ENDDO FOR                                                            
C       ELSE IF (# OF REFERENCE STAR MATCHES FOR CLUMP = 1) THEN               
C         IDFCLM(CLUMP#) = 2    MARK CLUMP AS IDENTIFIED                       
C       ELSE MORE THAN 1 REFERENCE STAR MATCH FOR CLUMP                        
C         IDFCLM(CLUMP#) = 1    MARK CLUMP AS QUESTIONABLE MATCH               
C         COMPUTE BEST POSITION MATCH FROM CANDIDATE REF STARS AND             
C           PLACE THIS REF STAR AT BEGINNING OF REF STAR LIST FOR              
C           THIS CLUMP.                                                        
C       ENDIF                                                                  
C     ENDDO FOR                                                                
C                                                                              
C                                                                              
C     RETURN                                                                   
C     -------------                                                            
C     ERROR_HANDLE:                                                            
C     DO CASE OF ERROR CONDITION                                               
C       IRCODE = 1:  CALL UTMSG TO OUTPUT ERROR MESSAGE                        
C     ENDCASE                                                                  
C     RETURN                                                                   
C-----------------------------------------------------------------------       
C                                                                              
C     * DEFINE PARAMETER VARIABLES                                             
      REAL*8    DANGTL      , DMAGTL                                           
C                                                                              
      REAL*4    DATCAT(7,*) , GCICLM(3,*), BRICLM(*), SKYCLM(10,3,*)           
C                                                                              
      INTEGER*4 IFXCLM      , NUMCAT     , IDNCAT(*)                           
      INTEGER*4 NUMSTR      , KLMSTR(*)  , MRKSTR(*)                           
      INTEGER*4 NUMCLM      , LBLCLM(*)  , NOBCLM(*)                           
      INTEGER*4 MRKCLM(*)   , IDFCLM(*)  , NRFCLM(*)                           
      INTEGER*4 MAPCLM(10,*), IRCODE                                           
C                                                                              
C     * DECLARE LOCAL VARIABLES                                                
      REAL*8    ANGSEP  , ANGDIF   , TOLER    , DANGLE, BESTA                  
      REAL*8    E2CLM(3), E2REF(3) , UE2REF(3), RANGTL, RMAGDF, RMAX           
      INTEGER*4 ICLM    , ISTAR    , IOBS     , IBEST , ITEMP , IMATCH         
      INTEGER*4 I                                                              
      DATA TOLER /0.99D0/                                                      
      DATA RMAX  /999999.0D0/                                                  
C                                                                              
C                                      INITIALIZE ROUTINE                      
      IRCODE = 0                                                               
      IF (LEVDBG(7) .NE. 0) WRITE (LUDBUG,1000)                                
C                                                                              
C                                      INITIALIZE ALGORITHM                    
      CALL DIRINI (NUMSTR, NUMCLM, IFXCLM,                                     
     1             MRKCLM, IDFCLM, NRFCLM, MAPCLM, SKYCLM, MRKSTR)             
C                                                                              
C                                                                              
C                                      CONVERT ANGULAR SEPARATION              
C                                      TOLERANCE TO RADIANS                    
C     RANGTL = DANGTL * DTR                                                    
      RANGTL = DANGTL / RADEG                                                 
      IF (LEVDBG(7) .GE. 1) WRITE (LUDBUG,4000) RANGTL, DMAGTL                 
C                                                                              
C                                                                              
C                                      MATCH CLUMPS WITH REFERENCE STARS       
      DO 200 ICLM = 1,NUMCLM                                                   
        IF (MRKCLM(ICLM) .EQ. 0) THEN                                          
C                                                                              
C                                      INTERMEDIATE DEBUG                      
          IF (LEVDBG(7) .GE. 4) THEN                                           
            WRITE (LUDBUG,4010) LBLCLM(ICLM),(GCICLM(I,ICLM),I=1,3),           
     1                          BRICLM(ICLM)                                   
          ENDIF                                                                
C                                                                              
C                                      COMPARE CLUMP WITH ALL STARS            
          DO 100 ISTAR = 1,NUMCAT                                              
C                                                                              
C                                      COMPUTE ANGLE BETWEEN VECTORS           
            E2CLM(1) = DBLE(GCICLM(1,ICLM))                                    
            E2CLM(2) = DBLE(GCICLM(2,ICLM))                                    
            E2CLM(3) = DBLE(GCICLM(3,ICLM))                                    
            E2REF(1) = DBLE(DATCAT(1,ISTAR))                                   
            E2REF(2) = DBLE(DATCAT(2,ISTAR))                                   
            E2REF(3) = DBLE(DATCAT(3,ISTAR))                                   
            ANGSEP   = DANGLE(E2CLM, E2REF, TOLER, ANGDIF)                     
C                                                                              
C                                      COMPUTE MAGNITUDE DIFFERENCE            
            RMAGDF = DBLE(ABS(BRICLM(ICLM)-DATCAT(4,ISTAR)))                   
C                                                                              
C                                      TEST FOR MATCH                          
            IF ((ANGSEP .LT. RANGTL).AND.(RMAGDF .LT. DMAGTL)) THEN            
              NRFCLM(ICLM) = NRFCLM(ICLM) + 1                                  
C                                                                              
C                                      INTERMEDIATE DEBUG                      
              IF (LEVDBG(7) .GE. 4) THEN                                       
                WRITE (LUDBUG,4020) NRFCLM(ICLM), ISTAR, IDNCAT(ISTAR),        
     1                 ANGSEP, RMAGDF                                          
              ENDIF                                                            
              IF (NRFCLM(ICLM) .GT. IFXCLM) THEN                               
                IRCODE = 1                                                     
                GO TO 9999                                                     
              ELSE                                                             
                MAPCLM(NRFCLM(ICLM),ICLM) = ISTAR                              
              ENDIF                                                            
            ENDIF                                                              
100       CONTINUE                                                             
        ENDIF                                                                  
200   CONTINUE                                                                 
C                                                                              
C                                                                              
C                                      MATCH OBSERVATIONS WITH THE             
C                                      REFERENCE STARS FROM ABOVE              
      DO 500 ICLM = 1,NUMCLM                                                   
        IF (NRFCLM(ICLM) .EQ. 0) THEN                                          
C                                      FLAG CLUMP AND ALL OBSERVATIONS         
C                                      IN CLUMP AS UNIDENTIFIED                
          IDFCLM(ICLM) = 0                                                     
          MRKCLM(ICLM) = 14                                                    
          DO 300 IOBS  = 1,NUMSTR                                              
            IF ((KLMSTR(IOBS).EQ. ICLM).AND.(MRKSTR(IOBS).EQ. 0)) THEN         
              MRKSTR(IOBS) = 16                                                
            ENDIF                                                              
300       CONTINUE                                                             
        ELSE IF (NRFCLM(ICLM) .EQ. 1) THEN                                     
C                                      FLAG CLUMP AS IDENTIFIED                
          IDFCLM(ICLM) = 2                                                     
        ELSE                                                                   
C                                      FLAG CLUMP AS QUESTIONABLE MATCH        
C                                      AND INITIALIZE BEST MATCH LOOP          
          IDFCLM(ICLM) = 1                                                     
          E2CLM(1) = DBLE(GCICLM(1,ICLM))                                      
          E2CLM(2) = DBLE(GCICLM(2,ICLM))                                      
          E2CLM(3) = DBLE(GCICLM(3,ICLM))                                      
          BESTA    = RMAX                                                      
          IBEST    = 0                                                         
C                                                                              
C                                      COMPUTE BEST MATCH, BY POSITION         
C                                      OF CANDIDATE REFERENCE STARS            
          DO 400 IMATCH = 1,NRFCLM(ICLM)                                       
            ISTAR    = MAPCLM(IMATCH,ICLM)                                     
            E2REF(1) = DBLE(DATCAT(1,ISTAR))                                   
            E2REF(2) = DBLE(DATCAT(2,ISTAR))                                   
            E2REF(3) = DBLE(DATCAT(3,ISTAR))                                   
            ANGSEP = DANGLE(E2CLM, E2REF, TOLER, ANGDIF)                       
            IF (ANGSEP .LT. BESTA) THEN                                        
              BESTA = ANGSEP                                                   
              IBEST = IMATCH                                                   
            ENDIF                                                              
400       CONTINUE                                                             
C                                                                              
C                                      MOVE BEST MATCH TO FRONT OF LIST        
          ITEMP = MAPCLM(1,ICLM)                                               
          MAPCLM(1,ICLM) = MAPCLM(IBEST,ICLM)                                  
          MAPCLM(IBEST,ICLM) = ITEMP                                           
        ENDIF                                                                  
500   CONTINUE                                                                 
C                                                                              
C                                                                              
C                                      NORMAL TERMINATION                      
      IF (LEVDBG(7) .NE. 0) WRITE (LUDBUG,2000)                                
      RETURN                                                                   
C                                                                              
C                                      ERROR HANDLING                          
9999  CONTINUE                                                                 
      IF (LEVDBG(7) .NE. 0) WRITE (LUDBUG,3000)                                
      IF (IRCODE .EQ. 1) THEN                                                  
c        IMSGNM = 60                                                            
c        IVARLN = 0                                                             
c        CALL UTMSG (C$SBID, IMSGNM, IVARLN, C$VDAT, IDSTFG, IRC)               
        IF (LEVDBG(7) .GE. 4) WRITE (LUDBUG,6000) ICLM                         
      ENDIF                                                                    
      RETURN                                                                   
C                                                                              
C                                      FORMAT SECTION                          
1000  FORMAT(' *** ENTER DMATCH ***')                                          
2000  FORMAT(' *** EXIT  DMATCH ***')                                          
3000  FORMAT(' *** ABEND DMATCH ***')                                          
4000  FORMAT(' DIRECT MATCH MAX ANGSEP=',D14.6,' RADIANS AND',                 
     1       ' MAGDIF MAX=',D14.6)                                             
4010  FORMAT(' WORKING ON CLUMP ',I6,' GCI=',3(F14.6,1X),' MAG=',F14.6)        
4020  FORMAT('         MATCH ',I4,' STARINDEX=',I6,' SKYMAP#=',I10,            
     1               ' ANGSEP=',D14.6,' MAGDIF=',D14.6)                        
6000  FORMAT(' DIRECT: MATCH TABLE OVERFILL WHILE MATCHING CLUMP ',I6)         
      END                                                                      
