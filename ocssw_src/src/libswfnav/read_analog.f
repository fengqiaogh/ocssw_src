        subroutine read_analog(data,loc,con,raw)
c  This subroutine converts the floating point value DATA from an integer
c  using the conversion specification CON and the location specified by LOC
c  in the array RAW
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
        data itmp2/0/
        
        data = 0.0
#ifdef LINUX
        do j=1,loc(2)
          temp2(1) = raw(loc(1)+j-1) 
          data = data*256.0 + itmp2
        end do
#else
        do j=1,loc(2)
          temp2(2) = raw(loc(1)+j-1) 
          data = data*256.0 + itmp2
        end do
#endif
        data = data*con(1) + con(2)
        return
        end
