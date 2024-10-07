      SUBROUTINE SORTCL                                                         
     I          (NUMCLM, NUMSTR, IFXCLM,                                        
     O           LBLCLM, TIMCLM, BRICLM, GCICLM, VRMCLM, VRPCLM, NOBCLM,        
     O           MRKCLM, IDFCLM, NRFCLM, MAPCLM, SKYCLM, KLMSTR, IDFHST)        
C-----------------------------------------------------------------------        
C   MODULE NAME: STSORTCL                                                       
C                                                                               
C   PURPOSE:TO BUBBLE SORT THE CLUMP TABLE IN INCREASING ORDER BY NUMBER        
C           OF DIRECT MATCH REFERENCE STAR MATCHES, AND ALTER THE STAR          
C           OBSERVATION TABLE ACCORDINGLY.                                      
C                                                                               
C                                                                               
C   ARGUMENT LIST:                                                              
C   ARGUMENT I/O TYPE DIMENS DESCRIPTION                                        
C   -------- --- ---- ------ -----------                                        
C   NUMCLM   I   I*4         NUMBER OF CLUMPS                                   
C   NUMSTR   I   I*4         NUMBER OF OBSERVATIONS                             
C   IFXCLM   I   I*4         MAX NUMBER OF REF STAR MATCHES PER CLUMP           
C   LBLCLM   I O I*4    *    CLUMP LABLES                                       
C   TIMCLM   I O R*8    *    AVERAGE CLUMP TIMES                                
C   BRICLM   I O R*4    *    AVERAGE CLUMP MAGNITUDES                           
C   GCICLM   I O R*4   3,*   AVERAGE CLUMP POSITIONS (GCI)                      
C   VRMCLM   I O R*4    *    CLUMP MAGNITUDE VARIANCES                          
C   VRPCLM   I O R*4    *    CLUMP POSITION VARIANCES                           
C   NOBCLM   I O I*4    *    NUMBER OF OBSERVATIONS IN EACH CLUMP               
C   MRKCLM   I O I*4    *    STATUS FLAG FOR EACH CLUMP                         
C   IDFCLM   I O I*4    *    IDENTIFICATION FLAG FOR EACH CLUMP                 
C   NRFCLM   I O I*4    *    # OF REFERENCE STAR MATCHES FOR EACH CLUMP         
C   MAPCLM   I O I*4  10,*   SKYMAP ID NUMBERS OF REFERENCE STAR MATCHES        
C   SKYCLM   I O R*4 10,3,*  LIST OF CORRECTED REFERENCE STAR POSITIONS         
C                                (MATCH#,AXIS,CLUMP#)                           
C   KLMSTR   I O I*4    *    CLUMP NUMBER FOR EACH OBSERVATION                  
C   IDFHST   I O I*4    *    FHST NUMBER FOR EACH CLUMP                         
C                                                                               
C   COMMON BLOCK VARIABLES USED:                                                
C   COMMON   VAR    I/O   VAR    I/O   VAR    I/O   VAR    I/O                  
C   ------   ---    ---   ---    ---   ---    ---   ---    ---                  
C   CMDEBG   LEVDBG I     LUDBUG I                                              
C                                                                               
      IMPLICIT NONE                                                             
C       ++INCLUDE STCMDEBG
      INTEGER*4 LEVDBG(8),LUDBUG
      COMMON /CMDEBG/LEVDBG,LUDBUG                                    
C                                                                               
C   EXTERNAL FILES REFERENCED:                                                  
C   FILENAME      OPERATION   FORTRAN UNIT ID                                   
C   --------      ---------   ---------------                                   
C   NONE                                                                        
C                                                                               
C   EXTERNAL REFERENCES:                                                        
C   --------------------------------------------------------------------        
C   STSWAPCL - SWAP ALL CLUMP ELEMENTS                                          
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
C   NONE                                                                        
C                                                                               
C   DEVELOPMENT HISTORY:                                                        
C   DATE     AUTHOR              DESCRIPTION                                    
C   -------- ------              -----------                                    
C   8/ 3/88  R.J. BURLEY         DESIGN                                         
C   5/16/89  R.J. BURLEY         CODED                                          
C  11/06/89  R.J. BURLEY         ADD LBLCLM GESS ARRAYS                         
C  02/04/92  C.C. YEH            ADD IDFHST (MTASS-11)                          
C-----------------------------------------------------------------------        
C   METHOD:                                                                     
C     SORTED = .FALSE.                                                          
C     DO WHILE NOT YET SORTED                                                   
C       SORTED = .TRUE.                                                         
C       DO FOR CLUMP# = 2, NUMCLM                                               
C         IF (#DIRECT MATCHES FOR CURRENT CLUMP IS LESS THAN                    
C             #DIRECT MATCHES FOR PREVIOUS CLUMP)                               
C           SORTED = .FALSE.                                                    
C           CALL SWAPCL TO SWAP ALL CLUMP ELEMENTS FROM PREVIOUS                
C            CLUMP TO CURRENT CLUMP                                             
C           DO FOR ALL OBSERVATIONS                                             
C             IF (OBSERVATION IS IN PREVIOUS CLUMP) THEN                        
C               MARK IT AS IN CLUMP                                             
C             ELSE IF (OBSERVATION IS IN CLUMP) THEN                            
C               MARK IT AS IN PREVIOUS CLUMP                                    
C             ENDIF                                                             
C           ENDDO FOR                                                           
C         ENDIF                                                                 
C       ENDDO FOR                                                               
C     ENDDO WHILE                                                               
C     RETURN                                                                    
C-----------------------------------------------------------------------        
C                                                                               
C     * DEFINE PARAMETER VARIABLES                                              
      REAL*8    TIMCLM(*)                                                       
C                                                                               
      REAL*4    BRICLM(*), GCICLM(3,*), VRMCLM(*)   , VRPCLM(*)                 
      REAL*4    SKYCLM(10,3,*)                                                  
C                                                                               
      INTEGER*4 NUMCLM      , NUMSTR   , IFXCLM     , LBLCLM(*)                 
      INTEGER*4 NOBCLM(*)   , MRKCLM(*), IDFCLM(*)  , NRFCLM(*)                 
      INTEGER*4 MAPCLM(10,*), KLMSTR(*), IDFHST(*)                              
C                                                                               
C     * DECLARE LOCAL VARIABLES                                                 
      INTEGER*4 ICLM  , IPREV, IOBS, ITEMP                                      
      LOGICAL*4 L_DONE                                                          
C                                                                               
C                                      INITIALIZE ROUTINE                       
      WRITE(*,*) 'ENTERING SORTCL'
      IF (LEVDBG(7) .NE. 0) WRITE (LUDBUG,1000)                                 
C                                                                               
C                                                                               
C                                      SORT CLUMPS BY # OF MATCHES              
      L_DONE = .FALSE.                                                          
100   CONTINUE                                                                  
      IF (.NOT. L_DONE) THEN                                                    
        L_DONE = .TRUE.                                                         
        DO 300 ICLM = 2,NUMCLM                                                  
          IPREV = ICLM - 1                                                      
          IF (NRFCLM(ICLM) .LE. NRFCLM(IPREV)) THEN                             
           IF ((NRFCLM(ICLM) .LT. NRFCLM(IPREV)).OR.
     *         (VRPCLM(ICLM) .LT. VRPCLM(IPREV))) THEN 
C                                                                               
C                                      NOT IN ORDER, SWAP CLUMPS                
            ITEMP = IDFHST(ICLM)                                                
            IDFHST(ICLM) = IDFHST(IPREV)                                        
            IDFHST(IPREV) = ITEMP                                               
            L_DONE = .FALSE.                                                    
            CALL SWAPCL (ICLM  , IPREV , IFXCLM,                                
     1                   LBLCLM, TIMCLM, BRICLM, GCICLM, VRMCLM, VRPCLM,        
     2                   NOBCLM, MRKCLM, IDFCLM, NRFCLM, MAPCLM, SKYCLM)        
C                                                                               
C                                      UPDATE OBSERVATION CLUMP #'S             
            DO 200 IOBS = 1,NUMSTR                                              
              IF (KLMSTR(IOBS) .EQ. IPREV) THEN                                 
                KLMSTR(IOBS) = ICLM                                             
              ELSE IF (KLMSTR(IOBS) .EQ. ICLM) THEN                             
                KLMSTR(IOBS) = IPREV                                            
              ENDIF                                                             
200         CONTINUE                                                            
           ENDIF
          ENDIF                                                                 
300     CONTINUE                                                                
        GO TO 100                                                               
      ENDIF                                                                     
C                                                                               
C                                      NORMAL TERMINATION                       
      IF (LEVDBG(7) .NE. 0) WRITE (LUDBUG,2000)                                 
      WRITE(*,*) 'EXIT SORTCL'
      RETURN                                                                    
C                                                                               
C                                      FORMAT SECTION                           
1000  FORMAT(' *** ENTER SORTCL *** ')                                          
2000  FORMAT(' *** EXIT  SORTCL *** ')                                          
      END                                                                      

