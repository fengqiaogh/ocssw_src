        subroutine read_discrete(bdat,dis,raw)
c  This subroutine extracts the discrete value BDAT from a byte 
c  using the location specified by DIS in the array RAW
c
c Modification History
c
c  Changed integer*1 to byte for Sun OS compatibility, B. A. Franz,
c  GAC, November 14, 1997.

        integer*4 dis(3)
        integer*2 mask(8),itmp2,idat
        byte raw(*),bdat,temp2(2)
        equivalence(temp2,itmp2)
        data mask/128,64,32,16,8,4,2,1/
        
        idat = 0
        itmp2 = 0
#ifdef LINUX
        temp2(1) = raw(dis(1))
#else
        temp2(2) = raw(dis(1))
#endif
        do j=1,dis(3)
          idat = idat + iand(itmp2,mask(dis(2)+j-1))
        end do
        bdat = idat/mask(dis(2)+dis(3)-1)
        return
        end
