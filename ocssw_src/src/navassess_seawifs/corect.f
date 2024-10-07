      SUBROUTINE CORECT(IDNCAT, DATCAT, NUMCLM, LBLCLM, TIMCLM,                      
     *            NRFCLM, MAPCLM, MRKCLM, SKYCLM, IERR)                        

      REAL*8    TIMCLM(*)                         
C                                                                               
      REAL*4    DATCAT(7,*)                                       
      REAL*4    SKYCLM(10,3,*)                                    
C                                                                               
      INTEGER*4 IDNCAT(*)              
      INTEGER*4 NUMCLM      , MRKCLM(*), NUMSTR                 
      INTEGER*4 NRFCLM(*)              
      INTEGER*4 LBLCLM(*)   , MAPCLM(10,*)  
C                                                                               
C     * DECLARE LOCAL VARIABLES                                                 
      INTEGER*4 IERR        , NUMTRP        , NUMDUB   , NUMCAT  ,LUCAT
C                                                                               

      DO I=1,NUMCLM
        DO J=1,NRFCLM(I)
          DO K=1,3
            SKYCLM(J,K,I) = DATCAT(K,MAPCLM(J,I))
          END DO
        END DO
      END DO

      IERR = 0
      RETURN
      END
