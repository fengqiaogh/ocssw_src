        integer*4 function lenstr(string)
c
c  Function lenstr(string)
c
c  Purpose: find the string length to the first blank or 0
c
c  Calling Arguments:
c
c  Name         Type    I/O     Description
c  --------     ----    ---     -----------
c  string       C*(*)    I      string to check
c  lenstr       I*4      O      length of string to first blank
c
c  By: W. Robinson, GSC, 26 Mar 93
c
c  Notes:  Any longer descriptions needed
c
c  Modification History:
c
      implicit none
c
      character string*(*)
      integer*4 i
c
c
c     loop through and find the first blank or 0
c
      lenstr = len( string )
      do i = 1, len( string )
         if( string(i:i) .eq. ' '  .or.  string(i:i) .eq. char(0) ) then
            lenstr = i - 1
            go to 990
         end if
      end do
  990 return
      end
