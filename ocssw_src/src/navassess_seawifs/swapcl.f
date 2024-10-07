      SUBROUTINE SWAPCL                                                         
     I    (ICLMP1, ICLMP2, IFXCLM,                                              
     O     LBLCLM, TIMCLM, BRICLM, GCICLM, VRMCLM, VRPCLM, NOBCLM,              
     O     MRKCLM, IDFCLM, NRFCLM, MAPCLM, SKYCLM)                              
C-----------------------------------------------------------------------        
C   MODULE NAME: STSWAPCL                                                       
C                                                                               
C                                                                               
C   PURPOSE: TO SWAP ALL ELEMENTS OF 2 DIFFERENT CLUMPS DURING SORT.            
C                                                                               
C                                                                               
C   ARGUMENT LIST:                                                              
C   ARGUMENT I/O TYPE DIMENS DESCRIPTION                                        
C   -------- --- ---- ------ -----------                                        
C   ICLMP1   I   I*4         FIRST OF 2 CLUMPS TO SWAP                          
C   ICLMP2   I   I*4         SECOND CLUMP TO SWAP                               
C   IFXCLM   I   I*4         MAX NUMBER OF REF STAR MATCHES PER CLUMP           
C   LBLCLM   I O I*4    *    CLUMP LABELS                                       
C   TIMCLM   I O R*8    *    AVERAGE CLUMP TIMES                                
C   BRICLM   I O R*4    *    AVERAGE CLUMP MAGNITUDES                           
C   GCICLM   I O R*4   3,*   AVERAGE CLUMP POSITIONS (GCI)                      
C   VRMCLM   I O R*4    *    CLUMP MAGNITUDE VARIANCES                          
C   VRPCLM   I O R*4    *    CLUMP POSITION VARIANCES                           
C   NOBCLM   I O I*4    *    NUMBER OF OBSERVATIONS IN EACH CLUMP               
C   MRKCLM   I O I*4    *    STATUS FLAG FOR EACH CLUMP                         
C   IDFCLM   I O I*4    *    IDENTIFICATION FLAG FOR EACH CLUMP                 
C   NRFCLM   I O I*4    *    # OF REFERENCE STAR MATCHES FOR EACH CLUMP         
C   MAPCLM   I O I*4  10,*   SKYMAP ID NUMBER OF REFERENCE STAR MATCHES         
C   SKYCLM   I O R*4 10,3,*  LIST OF CORRECTED REFERENCE STAR POSITIONS         
C                               (MATCH#,AXIS,CLUMP#)                            
C                                                                               
C   COMMON BLOCK VARIABLES USED:                                                
C   COMMON   VAR    I/O   VAR    I/O   VAR    I/O   VAR    I/O                  
C   ------   ---    ---   ---    ---   ---    ---   ---    ---                  
C   CMDEBG   LEVDBG I     LUDBUG I                                              
C                                                                               
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
C     NONE                                                                      
C                                                                               
C   SUBROUTINE CALLED FROM:                                                     
C   --------------------------------------------------------------------        
C     STSORTCL - STAR MATCHING SORT ROUTINE                                     
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
C   8/ 3/88   R.J. BURLEY         DESIGN                                        
C   5/16/89   R.J. BURLEY         CODED                                         
C  11/ 6/89   R.J. BURLEY         ADD LBLCLM GESS ARRAY                         
C-----------------------------------------------------------------------        
C   METHOD:                                                                     
C     PUT ALL ELEMENTS OF ICLMP1 INTO TEMP VARIABLES                            
C     PUT ALL ELEMENTS OF ICLMP2 INTO ICLMP1                                    
C     PUT ALL ELEMENTS OF TEMP VARIABLES INTO ICLMP2                            
C     RETURN                                                                    
C-----------------------------------------------------------------------        
C                                                                               
C     * DEFINE PARAMETER VARIABLES                                              
      REAL*8    TIMCLM(*)                                                       
C                                                                               
      REAL*4    BRICLM(*), GCICLM(3,*), VRMCLM(*)   , VRPCLM(*)                 
      REAL*4    SKYCLM(10,3,*)                                                  
C                                                                               
      INTEGER*4 ICLMP1   , ICLMP2     , IFXCLM      , LBLCLM(*)                 
      INTEGER*4 NOBCLM(*), MRKCLM(*)  , IDFCLM(*)   , NRFCLM(*)                 
      INTEGER*4 MAPCLM(10,*)                                                    
C                                                                               
C     * DECLARE LOCAL VARIABLES                                                 
      REAL*8    TEMP8                                                           
      REAL*4    TEMP4                                                           
      INTEGER*4 ITEMP, IAXIS, IMATCH                                            
C                                                                               
C                                      INITIALIZE ROUTINE                       
C      IF (LEVDBG(7) .NE. 0) WRITE (LUDBUG,1000)                                 
C                                                                               
C                                      SWAP CLUMP LABELS                        
      ITEMP = LBLCLM(ICLMP1)                                                    
      LBLCLM(ICLMP1) = LBLCLM(ICLMP2)                                           
      LBLCLM(ICLMP2) = ITEMP                                                    
C                                                                               
C                                      SWAP CLUMP TIME                          
      TEMP8 = TIMCLM(ICLMP1)                                                    
      TIMCLM(ICLMP1) = TIMCLM(ICLMP2)                                           
      TIMCLM(ICLMP2) = TEMP8                                                    
C                                                                               
C                                      SWAP CLUMP POSITIONS                     
      DO 100 IAXIS=1,3                                                          
        TEMP4 = GCICLM(IAXIS,ICLMP1)                                            
        GCICLM(IAXIS,ICLMP1) = GCICLM(IAXIS,ICLMP2)                             
        GCICLM(IAXIS,ICLMP2) = TEMP4                                            
100   CONTINUE                                                                  
C                                                                               
C                                      SWAP CLUMP MAGNITUDE                     
      TEMP4 = BRICLM(ICLMP1)                                                    
      BRICLM(ICLMP1) = BRICLM(ICLMP2)                                           
      BRICLM(ICLMP2) = TEMP4                                                    
C                                                                               
C                                      SWAP CLUMP MAGNITUDE VARIANCE            
      TEMP4 = VRMCLM(ICLMP1)                                                    
      VRMCLM(ICLMP1) = VRMCLM(ICLMP2)                                           
      VRMCLM(ICLMP2) = TEMP4                                                    
C                                                                               
C                                      SWAP CLUMP POSITION VARIANCE             
      TEMP4 = VRPCLM(ICLMP1)                                                    
      VRPCLM(ICLMP1) = VRPCLM(ICLMP2)                                           
      VRPCLM(ICLMP2) = TEMP4                                                    
C                                                                               
C                                      SWAP NUMBER OF OBSERVATIONS              
      ITEMP = NOBCLM(ICLMP1)                                                    
      NOBCLM(ICLMP1) = NOBCLM(ICLMP2)                                           
      NOBCLM(ICLMP2) = ITEMP                                                    
C                                                                               
C                                      SWAP CLUMP STATUS FLAGS                  
      ITEMP = MRKCLM(ICLMP1)                                                    
      MRKCLM(ICLMP1) = MRKCLM(ICLMP2)                                           
      MRKCLM(ICLMP2) = ITEMP                                                    
C                                                                               
C                                      SWAP IDENTIFICATION FLAGS                
      ITEMP = IDFCLM(ICLMP1)                                                    
      IDFCLM(ICLMP1) = IDFCLM(ICLMP2)                                           
      IDFCLM(ICLMP2) = ITEMP                                                    
C                                                                               
C                                      SWAP NUMBER OF MATCHES                   
      ITEMP = NRFCLM(ICLMP1)                                                    
      NRFCLM(ICLMP1) = NRFCLM(ICLMP2)                                           
      NRFCLM(ICLMP2) = ITEMP                                                    
C                                                                               
C                                      SWAP MATCH CATALOG #'S                   
      DO 200 IMATCH=1,10                                                        
        ITEMP = MAPCLM(IMATCH,ICLMP1)                                           
        MAPCLM(IMATCH,ICLMP1) = MAPCLM(IMATCH,ICLMP2)                           
        MAPCLM(IMATCH,ICLMP2) = ITEMP                                           
200   CONTINUE                                                                  
C                                                                               
C                                      SWAP CORRECTED REFERENCE                 
C                                      STAR POSITIONS                           
      DO 400 IMATCH=1,10                                                        
        DO 300 IAXIS=1,3                                                        
          TEMP4 = SKYCLM(IMATCH,IAXIS,ICLMP1)                                   
          SKYCLM(IMATCH,IAXIS,ICLMP1) = SKYCLM(IMATCH,IAXIS,ICLMP2)             
          SKYCLM(IMATCH,IAXIS,ICLMP2) = TEMP4                                   
300     CONTINUE                                                                
400   CONTINUE                                                                  
C                                                                               
C                                      NORMAL TERMINATION                       
C      IF (LEVDBG(7) .NE. 0) WRITE (LUDBUG,2000)                                 
      RETURN                                                                    
C                                                                               
C                                      FORMAT SECTION                           
1000  FORMAT(' *** ENTER SWAPCL *** ')                                          
2000  FORMAT(' *** EXIT  SWAPCL *** ')                                          
      END                                                                       
