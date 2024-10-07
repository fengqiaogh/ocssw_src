        real*4 function vmag(vec)
c
c  real*4 function vmag(vec)
c
c  Purpose: create the magnitude of a vector of size 3
c
c  Calling Arguments:
c
c  Name         Type    I/O     Description
c  --------     ----    ---     -----------
c  vec          R*4      I      size 3 vector
c  vmag         R*4      O      magnitude of the vector
c
c  By: W. Robinson, GSC, 13 Apr 93
c
c  Notes:  
c
c  Modification History:
c
      implicit none
      real*4 vec(3)
c
c
      vmag = sqrt( vec(1) * vec(1) + vec(2) * vec(2) + vec(3) * vec(3) )
c
      return
      end
