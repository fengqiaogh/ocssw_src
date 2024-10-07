      subroutine filenv(infil,outfil)
c
c  filenv(infil,outfil)
c
c  Purpose: detect any environmental names in file name and convert
c           it to the proper name
c
c  Calling Arguments:
c
c  Name         Type    I/O     Description
c  --------     ----    ---     -----------
c  infil        C*(*)    I      input file name to check
c  outfil       C*(*)    O      output file name
c
c  By: W. Robinson, GSC, 26 Mar 93
c
c  Notes:  this only looks at the start of the string for
c       a '$' to indicate an environment variable and 
c       then translates it to the output and adds the 
c       rest of the file description.
c
c  Modification History:
c
c  Enlarged envar and string to length 256 for seadas.  B. A. Franz,
c  November 14, 1997.
c
      implicit none
c
      character infil*(*), outfil*(*)
c
      character envar*255, string*256
      integer*4 iptr, lenstr
c
c
c       start, see if an environment variable is involved
c
      iptr = 1
      string = infil
      outfil = ''
      if( string(1:1) .eq. '$' ) then
c
c          extract the environmental variable
c
         iptr = index( string, '/' ) - 1
         envar = string(2:iptr)
         iptr = iptr + 1
c
c          translate environmental variable to output file
c
         call getenv( envar, outfil )
c
      end if
c
c       add the rest of the file name to the output file
c
      outfil = outfil(1:lenstr(outfil)) // string(iptr:lenstr(string))
c
c       and exit
c
      return
      end
