      SUBROUTINE DIRINI                                                        
     I          (NUMSTR, NUMCLM, IFXCLM,                                       
     O           MRKCLM, IDFCLM, NRFCLM, MAPCLM, SKYCLM, MRKSTR)               
C-----------------------------------------------------------------------       
C   MODULE NAME: STDIRINI                                                      
C                                                                              
C                                                                              
C   PURPOSE:  TO INITIALIZE DIRECT MATCH ALGORITHM BY INITIALIZING             
C             SEVERAL FIELDS IN THE CLUMP AND OBSERVATION DATA.                
C                                                                              
C                                                                              
C   ARGUMENT LIST:                                                             
C   ARGUMENT I/O TYPE DIMENS DESCRIPTION                                       
C   -------- --- ---- ------ -----------                                       
C   NUMSTR   I   I*4         NUMBER OF OBSERVATIONS                            
C   NUMCLM   I   I*4         NUMBER OF CLUMPS                                  
C   IFXCLM   I   I*4         MAX # OF REF STAR MATCHES                         
C                              (FIRST DIM. FOR MAPCLM & SKYCLM ARRAYS)         
C   MRKCLM   I O I*4    *    CLUMP STATUS FLAGS                                
C   IDFCLM   I O I*4    *    CLUMP IDENTIFICATION FLAGS                        
C   NRFCLM   I O I*4    *    NUMBER OF REFERENCE STAR MATCHES PER CLUMP        
C   MAPCLM   I O I*4  10,*   LIST OF REFERENCE STAR MATCH CATALOG #'S          
C   SKYCLM   I O R*4 10,3,*  LIST OF CORRECTED REFERENCE STAR POSITIONS        
C   MRKSTR   I O I*4    *    OBSERVATION STATUS FLAGS                          
C                                                                              
C                                                                              
C   COMMON BLOCK VARIABLES USED:                                               
C   COMMON   VAR    I/O   VAR    I/O   VAR    I/O   VAR    I/O                 
C   ------   ---    ---   ---    ---   ---    ---   ---    ---                 
C   CMDEBG   LEVDBG I     LUDBUG I                                             
C                                                                              
C      ++INCLUDE STCMDEBG                                                      
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
C   NONE                                                                       
C                                                                              
C   SUBROUTINE CALLED FROM:                                                    
C   --------------------------------------------------------------------       
C     STDIRECT - DIRECT MATCH ALGORITHM                                        
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
C   8/ 8/88  R.J. BURLEY         DESIGN                                        
C   5/11/89  R.J. BURLEY         CODED                                         
C-----------------------------------------------------------------------       
C   METHOD:                                                                    
C     DO FOR ALL OBSERVATIONS                                                  
C       IF (OBS HAS BEEN DROPPED BECAUSE OF NON-IDENTIFICATION) THEN           
C         SET STATUS FLAG FOR THIS OBSERVATION TO 0                            
C       ENDIF                                                                  
C     ENDDO FOR                                                                
C     DO FOR ALL CLUMPS                                                        
C       IF ((MRKCLM(CLUMP#) = 0).OR.(MRKCLM(CLUMP#) INDICATES DROPPED          
C         DUE TO NON-IDENTIFICATION)) THEN                                     
C         SET CLUMP STATUS FLAG TO 0                                           
C         SET CLUMP IDENTIFICATION FLAG TO 0                                   
C         SET # OF REFERENCE STAR MATCHES TO 0                                 
C         SET ALL ELEMENTS OF LIST OF MATCHED SKYMAP #'S TO 0                  
C         SET ALL ELEMENTS OF LIST OF CORRECTED REF STAR POSITIONS TO 0        
C       ENDIF                                                                  
C     ENDDO FOR                                                                
C     RETURN                                                                   
C-----------------------------------------------------------------------       
C                                                                              
C     * DEFINE PARAMETER VARIABLES                                             
      REAL*4  SKYCLM(10,3,*)                                                   
C                                                                              
      INTEGER*4 NUMSTR   , NUMCLM   , IFXCLM                                   
      INTEGER*4 MRKCLM(*), IDFCLM(*), NRFCLM(*), MAPCLM(10,*), MRKSTR(*)       
C                                                                              
C     * DECLARE LOCAL VARIABLES                                                
      INTEGER*4 IOBS, ICLM, IMATCH                                             
C                                                                              
C                                      INITIALIZE ROUTINE                      
      IF (LEVDBG(7) .NE. 0) WRITE (LUDBUG,1000)                                
C                                                                              
C                                                                              
C                                      UNDROP UNIDENTIFIED OBSERVATIONS        
      DO 100 IOBS = 1,NUMSTR                                                   
        IF ((MRKSTR(IOBS) .EQ. 16).OR.(MRKSTR(IOBS) .EQ. 17)) THEN             
          MRKSTR(IOBS) = 0                                                     
        ENDIF                                                                  
100   CONTINUE                                                                 
C                                                                              
C                                      UNDROP UNIDENTIFIED CLUMPS AND          
C                                      RESET IDENTIFICATION ARRAYS             
      DO 300 ICLM = 1,NUMCLM                                                   
        IF ((MRKCLM(ICLM) .EQ. 0).OR.(MRKCLM(ICLM) .EQ. 14).OR.                
     1      (MRKCLM(ICLM) .EQ. 15)) THEN                                       
          MRKCLM(ICLM)  = 0                                                    
          IDFCLM(ICLM)  = 0                                                    
          NRFCLM(ICLM)  = 0                                                    
          DO 200 IMATCH = 1,IFXCLM                                             
            MAPCLM(IMATCH,ICLM) = 0                                            
            SKYCLM(IMATCH,1,ICLM) = 0.0                                        
            SKYCLM(IMATCH,2,ICLM) = 0.0                                        
            SKYCLM(IMATCH,3,ICLM) = 0.0                                        
200       CONTINUE                                                             
        ENDIF                                                                  
300   CONTINUE                                                                 
C                                                                              
C                                      NORMAL TERMINATION                      
      IF (LEVDBG(7) .NE. 0) WRITE (LUDBUG,2000)                                
      RETURN                                                                   
C                                                                              
C                                      FORMAT SECTION                          
1000  FORMAT(' *** ENTER DIRINI ***')                                          
2000  FORMAT(' *** EXIT  DIRINI ***')                                          
      END                                                                      
