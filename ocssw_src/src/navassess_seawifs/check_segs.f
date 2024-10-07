      subroutine check_segs(ill,ilp,i1,i2,lflag1,lflag2,igood)

c  Subroutine to check for other land pixels contiguous to two segments
      integer*4 ill(*),ilp(2,*)
      byte lflag1(*),lflag2(*)
  
      igood = 0

c  Check at start of segments
      if (ilp(1,i2).gt.ilp(1,i1)) then
        do i=ilp(1,i1)-1,ilp(1,i2)-2
          if(lflag2(i).ge.1) igood=-1
        end do
      else if (ilp(1,i2).lt.ilp(1,i1)) then
        do i=ilp(1,i2)-1,ilp(1,i1)-2
          if(lflag1(i).ge.1) igood=-1
        end do
      end if

c  Check at end of segments
      if (ilp(2,i2).gt.ilp(2,i1)) then
        do i=ilp(2,i1)+2,ilp(2,i2)+1
          if(lflag1(i).ge.1) igood=-1
        end do
      else if (ilp(2,i2).lt.ilp(2,i1)) then
        do i=ilp(2,i2)+2,ilp(2,i1)+1
          if(lflag2(i).ge.1) igood=-1
        end do
      end if

      return
      end
