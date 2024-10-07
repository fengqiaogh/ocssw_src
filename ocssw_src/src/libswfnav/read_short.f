        subroutine read_short(data,loc,con,raw)

c $Header$
c $Log$
c
c  This subroutine converts the floating point value DATA from a signed 2-byte
c  integer using the conversion specification CON and the location specified 
c  by LOC in the array RAW
c
c  April 3, 1995 by Frederick S. Patt, GSC
c
c Modification History
c
c  Changed integer*1 to byte for Sun OS compatibility, B. A. Franz,
c  GAC, November 14, 1997.


        real*4 con(2)
        integer*4 loc(2)
        integer*2 itmp2
        byte temp2(2)
        byte raw(*)
        equivalence (itmp2,temp2)
        
        data = 0.0
#ifdef LINUX
        do j=1,2
          temp2(j) = raw(loc(1)+2-j) 
        end do
#else
        do j=1,2
          temp2(j) = raw(loc(1)+j-1) 
        end do
#endif
        data = itmp2*con(1) + con(2)
        return
        end
