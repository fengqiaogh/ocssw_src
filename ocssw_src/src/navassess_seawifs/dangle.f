**                                                                              
** INVOCATION NAME:  DANGLE                                                     
**                                                                              
** PURPOSE:  FUNCTION TO COMPUTE THE ANGLE BETWEEN TWO UNIT VECTORS             
**                                                                              
** INVOCATION METHOD:  X = DANGLE(A, B, TOLER, ANG)                             
**                                                                              
** ARGUMENT LIST:                                                               
** NAME    TYPE  USE PURPOSE                                                    
** A(3)    R*8    I  FIRST UNIT VECTOR                                          
** B(3)    R*8    I  SECOND UNIT VECTOR                                         
** TOLER   R*8    I  MAXIMUM MAGNITUDE OF DOT-PRODUCT BEFORE CROSS-             
**                   PRODUCT IS USED TO COMPUTE ANGLE.  IF DOT-PRODUCT          
**                   IS GREATER THAN TOLER, THEN THE ANGLE IS ARCSIN            
**                   OF THE MAGNITUDE OF THE CROSS-PRODUCT, OTHERWISE,          
**                   THE ARCCOS OF THE DOT-PRODUCT IS USED.                     
**                   (RECOMMENDED VALUE IS 0.99D0)                              
** ANG     R*8    O  ANGLE BETWEEN UNIT VECTORS A AND B IN RADIANS              
** DANGLE  R*8    O  ANGLE BETWEEN UNIT VECTORS A AND B IN RADIANS              
**                                                                              
** FILE/RECORD REFERENCES: NONE                                                 
**                                                                              
** GLOBAL REFERENCES: NONE                                                      
**                                                                              
** EXTERNAL REFERENCES: NONE                                                    
**                                                                              
** INTERNAL VARIABLES:                                                          
** VAR       TYPE  PURPOSE                                                      
** ACROSB(3) R*8   CROSS-PRODUCT OF A AND B                                     
** ADOTB     R*8   DOT-PRODUCT OF A AND B                                       
** CROSMG    R*8   MAGNITUDE OF CROSS-PRODUCT                                   
**                                                                              
** NOTES:                                                                       
**                                                                              
** CHANGE HISTORY:                                                              
** AUTHOR             CHG-ID   MMM.YY    CHG-SUMMARY                            
** B.GROVEMAN                  JUL.88    PDL                                    
** M.WOOLSEY                   AUG.88    CERTIFIED PDL                          
** R.COON                      NOV.88    ORIGINAL CODE                          
** B.GROVEMAN                  NOV.88    CERTIFIED CODE                         
** M. WOOLSEY                  AUG.89    CORRECTED ARCSIN METHOD                
**                                        RETURNING ANSWERS > PI / 2            
**                                                                              
** BEGIN PDL                                                                    
**      COMPUTE DOT-PRODUCT OF A AND B                                          
**      IF (DOT-PRODUCT LESS THAN TOLER) THEN                                   
**         ANGLE = ARCCOS(DOT-PRODUCT)                                          
**      ELSE                                                                    
**         COMPUTE CROSS-PRODUCT OF A AND B                                     
**         COMPUTE THE MAGNITUDE OF THE CROSS-PRODUCT                           
**         ANGLE = ARCSIN(MAGNITUDE OF CROSS-PRODUCT)                           
**         IF (DOT-PRODUCT IS LESS THAN ZERO) ANGLE = PI - ANGLE                
**      END IF                                                                  
**                                                                              
**      RETURN                                                                  
** END PDL                                                                      
**                                                                              
** MAIN ENTRY POINT TO DANGLE                                                   
************************************************************************        
C                                                                               
      FUNCTION DANGLE (A, B, TOLER, ANG)                                        
C                                                                               
C-----DECLARE ARGUMENTS                                                         
      REAL*8 A(3), B(3), TOLER, ANG, DANGLE                                     
C                                                                               
C-----DECLARE LOCAL VARIABLES                                                   
      REAL*8 ACROSB(3), ADOTB, CROSMG                                           
      REAL*8 PI /3.141592653589793D+00/                                         
C                                                                               
      ADOTB = A(1)*B(1) + A(2)*B(2) + A(3)*B(3)                                 
      IF (ABS(ADOTB) .LE. ABS(TOLER)) THEN                                      
C                          USE DOT-PRODUCT                                      
        ANG = ACOS(ADOTB)                                                       
        DANGLE = ANG                                                            
      ELSE                                                                      
C                          USE CROSS PRODUCT FOR GREATER ACCURACY               
        ACROSB(1) = A(2)*B(3) - A(3)*B(2)                                       
        ACROSB(2) = A(3)*B(1) - A(1)*B(3)                                       
        ACROSB(3) = A(1)*B(2) - A(2)*B(1)                                       
        CROSMG = SQRT(ACROSB(1)**2 + ACROSB(2)**2 + ACROSB(3)**2)               
        ANG = ASIN(CROSMG)                                                      
        IF (ADOTB .LT. 0.0D0)  ANG = PI - ANG                                   
        DANGLE = ANG                                                            
      END IF                                                                    
C                                                                               
      RETURN                                                                    
      END
