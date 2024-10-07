        real*4 xlat(500),xlon(500),xlatc(8000),xlonc(8000)
        integer*2 ipix(500),inum(500)
        character*50 iname

        type *,'Enter island file name'
        read (*,'(a50)') iname  
        open(51,file=iname,readonly)
        do i=1,500
          read(51,*,end=500) inum(i),xlat(i),xlon(i),ilin,ipix(i)
          ni = i
        end do

  500   open(52,file='/u3/fred/bak/islands.all',readonly)
        do i=1,8000
          read(52,*,end=600) xlonc(i),xlatc(i),wlon,wlat
          nc = i
        end do

  600   continue

        do i=1,ni
          dlm = 999. 
          ic = 0
          nn = 0
          do j=1,nc
            dlon = abs(xlon(i)-xlonc(j))
            dlat = abs(xlat(i)-xlatc(j))
            if ((dlon.lt.0.2).and.(dlat.lt.0.2)) then
            nn = nn + 1
            dll = sqrt(dlon*dlon+dlat*dlat)
              if (dll.lt.dlm) then
                ic = j
                dlm = dll
              end if
            end if
          end do
          if (ic.ne.0) then
            write (*,1200) i,xlon(i),xlat(i),ic,xlonc(ic),xlatc(ic),
     *          dlm,ipix(i),inum(i),nn
 1200       format (2(i5,2f11.5),f9.5,3i5)
c         else
c           write (*,1200) i,xlon(i),xlat(i)
          end if
        end do
        stop
        end

