CCCC                                                                            
C                                                                               
      SUBROUTINE DUNVEC (A, AU, AMAG)                                           
C                                                                               
    1 FORMAT(' *** SUBROUTINE DUNVEC  MARCH 11, 1987 ***')                      
C                                                                               
C-----LANGUAGE - VS FORTRAN                                                     
C                                                                               
C-----FUNCTION -  UNITIZES A VECTOR AND COMPUTES ITS MAGNITUDE                  
C                                                                               
C-----MATHEMATICAL METHOD -                                                     
C                                                                               
C-----ARGUMENTS -                                                               
C     ARGUMENT    TYPE  IO  DESCRIPTION                                         
C     --------    ----  --  -----------                                         
C     A(3)        R*8   I   INPUT VECTOR                                        
C     AU(3)       R*8   O   OUTPUT UNIT VECTOR OF A                             
C     AMAG        R*8   O   MAGNITUDE OF A                                      
C                                                                               
C-----EXTERNAL REFERENCES - NONE                                                
C                                                                               
C-----CALLED BY - ANY                                                           
C                                                                               
C-----COMMONS REFERENCED - NONE                                                 
C                                                                               
C-----ERROR HANDLING - NONE                                                     
C                                                                               
C-----FILES REFERENCED - NONE                                                   
C                                                                               
C-----DESIGNER - R. COON  CSC  03/01/87                                         
C                                                                               
C-----PROGRAMMER - R. COON  CSC  MARCH 11, 1987                                 
C                                                                               
C-----VERIFIED BY - K. P. WASCZKIEWICZ  CSC  MARCH 11, 1987                     
C                                                                               
C-----MODIFICATIONS -                                                           
C     NAME DATE DESCRIPTION                                                     
C M. WOOLSEY  JAN.89  INSTALLED UNDERFLOW ERROR TRAP                            
CCCC                                                                            
C     DECLARE ARGUMENTS                                                         
      REAL*8 A(3), AU(3), AMAG                                                  
C                                                                               
      AMAG = DSQRT(A(1)**2 + A(2)**2 + A(3)**2)                                 
      IF (AMAG .NE. 0.0D+00 .AND. DABS (AMAG) .LT. 1.0D32) THEN                 
        IF (DABS (A (1)) .GT. 1.0D-32) THEN                                     
          AU(1)=A(1)/AMAG                                                       
        ELSE                                                                    
          AU (1) = 0.D0                                                         
        END IF                                                                  
        IF (DABS (A (2)) .GT. 1.0D-32) THEN                                     
          AU(2)=A(2)/AMAG                                                       
        ELSE                                                                    
          AU (2) = 0.D0                                                         
        END IF                                                                  
        IF (DABS (A (3)) .GT. 1.0D-32) THEN                                     
          AU(3)=A(3)/AMAG                                                       
        ELSE                                                                    
          AU (3) = 0.D0                                                         
        END IF                                                                  
      ELSE                                                                      
        AU(1)=0.0D+00                                                           
        AU(2)=0.0D+00                                                           
        AU(3)=0.0D+00                                                           
      END IF                                                                    
C                                                                               
      RETURN                                                                    
      END
