      subroutine read_double(data,loc,con,raw)

c $Header$
c $Log$
c
c  This subroutine retrieves a 8-byte floating point value from a specified 
c  location in the raw data array and converts it to a 4-byte real variable.
c
c  April 3, 1995 by Frederick S. Patt, GSC
c
c Modification History
c
c  Changed integer*1 to byte for Sun OS compatibility, B. A. Franz,
c  GAC, November 14, 1997.


      real*8 dtmp8
      real*4 con(2)
      integer*4 loc(2)
      byte raw(*),temp(8)
      equivalence (temp,dtmp8)

      dtmp8 = 0.0
#ifdef LINUX
      do i=1,8
        temp(i) = raw(loc(1)+8-i)
      end do
#else
      do i=1,8
        temp(i) = raw(i+loc(1)-1)
      end do
#endif
      data = dtmp8*con(1)
      return
      end

