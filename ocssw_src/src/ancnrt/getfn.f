      subroutine getfn(filespec,lspec,name,lname)
c
c *** subroutine to return a filename given a full file specification
c
c     By Jim Firestone, GSC/SAIC, 5/11/94
c
      implicit none
      integer lspec, index, lnstrg, lname, slash, dot, i
      character filespec*(*), name*80

      do i = 1,80
        name(i:i) = ' '
      end do

      if (lspec .lt. 1) lspec = lnstrg(filespec)
      dot = index(filespec(1:lspec),'.')
c
c *** step backwards to preceding slash - everything between that slash and
c     '.' is part of the file name
c
      slash = 0
      do i = dot, 1, -1
        if (filespec(i:i) .eq. '/') then 
          slash = i
          goto 10
        end if
      end do

10    continue
      name = filespec(slash+1:dot-1)
      lname = lnstrg(name) 

      return
      end 
