module cubeio

      INTEGER, PARAMETER :: INT16 = SELECTED_INT_KIND(4)

contains
!***********************************************************************
!--- Subroutine to open an input imaging spectrometer data cube ---
      SUBROUTINE OPENINFILE(LUN_IN, I_RET)

      INTEGER, intent(in)  :: LUN_IN
      INTEGER, intent(out) :: I_RET
      CHARACTER (LEN = 1000) :: FINAV
      COMMON /INCUBE/  FINAV

      OPEN(LUN_IN,file=FINAV,access='stream',status='old',             &
     & form='unformatted',action='read', iostat=ios)
      I_RET = ios

      END SUBROUTINE OPENINFILE
!---------------------------------------------------------------------

!--- Subroutine to open an output reflectance data cube ---
      SUBROUTINE OPENOUTFILE(LUN_OUT, I_RET)

      INTEGER, intent(in) :: LUN_OUT
      INTEGER, intent(out) :: I_RET
      CHARACTER (LEN = 1000) :: FOCUB
      COMMON /OUTCUBE/  FOCUB

      OPEN(LUN_OUT,file=FOCUB,access='stream', status='unknown',        &
     & form='unformatted',action='write', iostat=ios)
      I_RET = ios

      END SUBROUTINE OPENOUTFILE
!---------------------------------------------------------------------
 
!--- Subroutine to open an output water vapor image file ---
      SUBROUTINE OPENVAPFILE(LUN_VAP, I_RET)

      INTEGER, intent(in) :: LUN_VAP
      INTEGER, intent(out) :: I_RET
      CHARACTER (LEN = 1000) :: FOH2O
      COMMON /OUTH2OVAP/ FOH2O

      OPEN(LUN_VAP,file=FOH2O,access='stream',status='unknown',         &   
     & form='unformatted',action='write', iostat=ios)
      I_RET = ios

      END SUBROUTINE OPENVAPFILE
!---------------------------------------------------------------------


!***********************************************************************
!--- Subroutine to read a slice of spatial & spectral data ----
      SUBROUTINE RD_SLICE(LUN_IN,NSAMPS,NBANDS,SORDER,BUFFER)
      
      INTEGER, intent(in) :: LUN_IN,NSAMPS,NBANDS,SORDER
      INTEGER(INT16), DIMENSION(NBANDS*NSAMPS), intent(out) :: BUFFER

      INTEGER(INT16), DIMENSION(NSAMPS,NBANDS) :: BUFFER_BIL
      INTEGER(INT16), DIMENSION(NBANDS,NSAMPS) :: BUFFER_BIP

      READ(LUN_IN) BUFFER

      IF(SORDER == 2) then
        BUFFER_BIL = reshape(BUFFER, (/ NSAMPS, NBANDS /) )
        BUFFER_BIP = TRANSPOSE(BUFFER_BIL)
        BUFFER     = reshape(BUFFER_BIP, (/NBANDS*NSAMPS/))
      END IF

      END SUBROUTINE RD_SLICE
!---------------------------------------------------------------------


!--- Subroutine to write a slice of spatial & spectral data ----
      SUBROUTINE WT_SLICE(LUN_OUT,NSAMPS,NBANDS,SORDER,BUFFER)
    
      INTEGER, intent(in) :: LUN_OUT,NSAMPS,NBANDS,SORDER

      INTEGER(INT16), DIMENSION(NBANDS*NSAMPS), intent(inout) :: BUFFER

      INTEGER(INT16), DIMENSION(NBANDS,NSAMPS) :: BUFFER_BIP
      INTEGER(INT16), DIMENSION(NSAMPS,NBANDS) :: BUFFER_BIL

      IF(SORDER == 1) WRITE(LUN_OUT) BUFFER

      IF(SORDER == 2) then
        BUFFER_BIP = reshape(BUFFER, (/NBANDS,NSAMPS/))
        BUFFER_BIL = transpose(BUFFER_BIP)
        BUFFER     = reshape(BUFFER_BIL, (/NBANDS*NSAMPS/))
        WRITE(LUN_OUT) BUFFER
      END IF

      END SUBROUTINE WT_SLICE
!---------------------------------------------------------------------


!--- Subroutine to write a 'line' of water vapor data ----
      SUBROUTINE WT_LINE(LUN_VAP,NSAMPS,H2OBUF)

      INTEGER, intent(in) :: LUN_VAP,NSAMPS

      INTEGER(INT16), DIMENSION(NSAMPS), intent(inout) :: H2OBUF

      WRITE(LUN_VAP) H2OBUF

      END SUBROUTINE WT_LINE

!***********************************************************************


!--- Subroutine to close the input imaging spectrometer data cube ---
      SUBROUTINE CLOSEINFILE(LUN_IN)
         INTEGER :: LUN_IN
         CLOSE(LUN_IN, status='keep')
      END SUBROUTINE CLOSEINFILE
!---------------------------------------------------------------------


!--- Subroutine to close the output reflectance data cube ---
      SUBROUTINE CLOSEOUTFILE(LUN_OUT)
         INTEGER :: LUN_OUT
         CLOSE(LUN_OUT, status='keep')
      END SUBROUTINE CLOSEOUTFILE


!--- Subroutine to close the output water vapor image file ---
      SUBROUTINE CLOSEVAPFILE(LUN_VAP)
         INTEGER :: LUN_VAP
         CLOSE(LUN_VAP, status='keep')
      END SUBROUTINE CLOSEVAPFILE

!***********************************************************************
end module cubeio
