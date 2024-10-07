      SUBROUTINE COLINE                                                        
     I          (ICLMA, ICLMB, ICLMC, LBLCLM, GCICLM,                          
     O           COANGL)                                                       
C-----------------------------------------------------------------------       
C   MODULE NAME:  STCOLINE                                                     
C                                                                              
C                                                                              
C   PURPOSE:  TO COMPUTE COLINEAR FUNCTION FOR STAR CLUMP TRIPLET              
C                                                                              
C                                                                              
C   ARGUMENT LIST:                                                             
C   ARGUMENT I/O TYPE DIMENS DESCRIPTION                                       
C   -------- --- ---- ------ -----------                                       
C   ICLMA    I   I*4         CLUMP NUMBER OF FIRST CLUMP OF TRIPLET            
C   ICLMB    I   I*4         CLUMP NUMBER OF SECOND CLUMP OF TRIPLET           
C   ICLMC    I   I*4         CLUMP NUMBER OF THIRD CLUMP OF TRIPLET            
C   LBLCLM   I   I*4         CLUMP LABLES                                      
C   GCICLM   I   R*4   3,*   AVERAGE CLUMP POSITIONS GCI                       
C   COANGL     O R*8         COLINEARITY ANGLE                                 
C                                                                              
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
C   UTDANGLE - FUNCTION TO COMPUTE ANGULAR SEPARATIONS BETWEEN VECTORS         
C                                                                              
C                                                                              
C   SUBROUTINE CALLED FROM:                                                    
C   --------------------------------------------------------------------       
C     STTRIPLT - STAR MATCHING TRIPLET ALGORITHM                               
C                                                                              
C                                                                              
C   CONSTRAINTS, RESTRICTIONS, MESSAGES, NOTES:                                
C   --------------------------------------------------------------------       
C   NONE                                                                       
C                                                                              
C   REQUIREMENTS REFERENCES:                                                   
C   --------------------------------------------------------------------       
C   UARS FDSS SPECS 3.1.1.5  (F5.3)                                            
C                                                                              
C   DEVELOPMENT HISTORY:                                                       
C   DATE     AUTHOR              DESCRIPTION                                   
C   -------- ------              -----------                                   
C   8/  8/88  R.J. BURLEY         DESIGN                                       
C   5/ 18/89  R.J. BURLEY         CODED                                        
C   10/24/89  R.J. BURLEY         SUBTRACT LARGEST ANGULAR SEPARATION          
C                                 FROM THE SUM OF THE OTHER 2, NOT             
C                                 VICE VERSA.                                  
C   11/ 6/89  R.J. BURLEY         ADD CLUMP LABLES                             
C-----------------------------------------------------------------------       
C   METHOD:                                                                    
C     COMPUTE ROTATION ANGLES BETWEEN 3 AVG CLUMP POSITIONS                  
C     RETURN THE SMALLEST ANGLE OF THE 3                             
C     REFERENCE:  WERTZ, APPENDIX A, EQ. A-2
C-----------------------------------------------------------------------       
C                                                                              
C     * DEFINE PARAMETER VARIABLES                                             
      REAL*8    COANGL                                                         
C                                                                              
      REAL*4    GCICLM(3,*)                                                    
C                                                                              
      INTEGER*4 LBLCLM(*),  ICLMA   , ICLMB    , ICLMC                         
C                                                                              
C     * DECLARE LOCAL VARIABLES                                                
      REAL*8    DANGLE   , ANGDIF   , TOLER                                    
      REAL*8    E2CLMA(3), E2CLMB(3), E2CLMC(3)                                
      REAL*8    ENCLMA(3), ENCLMB(3), ENCLMC(3), EMAG                     
      REAL*8    CANG12   , CANG13   , CANG23   , PI , PIO2
      REAL*8    CANG1    , CANG2    , CANG3
      REAL*4    CLMA(3)  , CLMB(3)  , CLMC(3)
      INTEGER*4 IAXIS                                                          
      DATA      TOLER /0.99D0/                                                 
      DATA PI/3.14159265359D0/,PIO2/1.570796327/
                                                                               
C                                                                              
C                                      INITIALIZE ROUTINE                      
C      IF (LEVDBG(7) .NE. 0) WRITE (LUDBUG,1000)                                
C                                                                              
C                                      COMPUTE ROTATION ANGLES             
      DO 100 IAXIS = 1,3                                                       
        E2CLMA(IAXIS) = GCICLM(IAXIS,ICLMA)                                    
        E2CLMB(IAXIS) = GCICLM(IAXIS,ICLMB)                                    
        E2CLMC(IAXIS) = GCICLM(IAXIS,ICLMC)                                    
100   CONTINUE                                                                 
      CANG12 = DANGLE(E2CLMA,E2CLMB, TOLER, ANGDIF)                            
      CANG13 = DANGLE(E2CLMA,E2CLMC, TOLER, ANGDIF)                            
      CANG23 = DANGLE(E2CLMB,E2CLMC, TOLER, ANGDIF)                            

      CANG1 = ACOS((COS(CANG23) - COS(CANG12)*COS(CANG13))/
     *  (SIN(CANG12)*SIN(CANG13)))
      CANG2 = ACOS((COS(CANG13) - COS(CANG12)*COS(CANG23))/
     *  (SIN(CANG12)*SIN(CANG23)))
      CANG3 = ACOS((COS(CANG12) - COS(CANG13)*COS(CANG23))/
     *  (SIN(CANG13)*SIN(CANG23)))
C                                                                              
C                                      FIND THE SMALLEST ANGLE            
      COANGL = DMIN1(CANG1,CANG2,CANG3)
C                                                                              
C                                      DEBUG                                   
      IF (LEVDBG(7) .GE. 4) THEN                                               
        WRITE (LUDBUG,4000) LBLCLM(ICLMA), LBLCLM(ICLMB),                      
     1                      LBLCLM(ICLMC), COANGL                              
      ENDIF                                                                    
C                                                                              
C                                      NORMAL TERMINATION                      
C      IF (LEVDBG(7) .NE. 0) WRITE (LUDBUG,2000)                                
      RETURN                                                                   
C                                                                              
C                                      FORMAT SECTION                          
1000  FORMAT(' *** ENTER COLINE ***')                                          
2000  FORMAT(' *** EXIT  COLINE ***')                                          
C4000  FORMAT(' STCOLINE DEBUG: ',                                             
4000  FORMAT(' COLINEAR ANGLE BETWEEN CLUMPS',3(2X,I6),' = ',D12.6)          
      END                                                                      
