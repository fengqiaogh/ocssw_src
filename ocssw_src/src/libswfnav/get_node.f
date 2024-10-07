        subroutine get_node(nlines,timref,time,pos,vel,xnodel,tnode)
c
c
c  Purpose:  Determine longitude and time of descending node crossing
c  
c
c  Calling Arguments:
c
c  Name         Type    I/O     Description
c  --------     ----    ---     -----------
c  nlines       I*4      I      Number of scan line times
c  timref       R*8      I      Size 3 reference time at start line 
c                                of data:  year, day, sec
c  time         R*8      I      Array of time in seconds relative to 
c                                timref for every scan line
c  pos          R*4      I      Size 3 by nlines interpolated position
c  vel          R*4      I      Size 3 by nlines interpolated velocity
c  xnodel       R*4      O      Longitude of descending node
c  tnode        R*8      O      Time of node crossing in seconds of day
c
c
c  By: Frederick S. Patt, GSC, July 13, 1994
c
c  Modification History:
c
      implicit none
c
      integer*4 ind, nlines 
      real*8 timref(3), time(nlines), tnode
      real*4 pos(3,nlines), vel(3,nlines), xnodel
      real*4 x, y, rfac, pmag, vmag, on(3), oxymag, orbang

      real*8 pi,radeg,re,rem,f,omf2,omegae
      common /gconst/pi,radeg,re,rem,f,omf2,omegae


c  If orbit data span node crossing, find crossing
      if ((pos(3,1).gt.0.).and.(pos(3,nlines).lt.0.)) then
        ind = 1
        dowhile (pos(3,ind).gt.0.)
          ind = ind + 1
        end do

c   Interpolate adjacent vectors to crossing and compute node and time
        rfac = pos(3,ind-1)/(pos(3,ind-1)-pos(3,ind))
        x = pos(1,ind)*rfac + pos(1,ind-1)*(1.0-rfac)
        y = pos(2,ind)*rfac + pos(2,ind-1)*(1.0-rfac)
        xnodel = radeg*atan2(y,x)
        tnode = time(ind)*rfac + time(ind-1)*(1.0-rfac) + timref(3)

c  Else, need to estimate crossing from closest vector
c   If above equator
      else if (pos(3,nlines).gt.0) then
        pmag = sqrt(pos(1,nlines)**2 + pos(2,nlines)**2 
     *       + pos(3,nlines)**2)
        vmag = sqrt(vel(1,nlines)**2 + vel(2,nlines)**2 
     *       + vel(3,nlines)**2)
c   Compute orbit normal and node longitude
        call crossp(pos(1,nlines),vel(1,nlines),on)
        xnodel = radeg*atan2(-on(1),on(2))
c   Compute angle from position to crossing and time difference
        oxymag = sqrt(on(1)**2 + on(2)**2)
        orbang = acos((pos(1,nlines)*on(2)-pos(2,nlines)*on(1))/
     *                (oxymag*pmag))
        tnode = orbang*pmag/vmag + time(nlines) + timref(3)

c   Else if below equator
      else
        pmag = sqrt(pos(1,1)**2 + pos(2,1)**2 
     *       + pos(3,1)**2)
        vmag = sqrt(vel(1,1)**2 + vel(2,1)**2 
     *       + vel(3,1)**2)
c   Compute orbit normal and node longitude
        call crossp(pos(1,1),vel(1,1),on)
        xnodel = radeg*atan2(-on(1),on(2))
c   Compute angle from position to crossing and time difference
        oxymag = sqrt(on(1)**2 + on(2)**2)
        orbang = -acos((pos(1,1)*on(2)-pos(2,1)*on(1))/
     *                (oxymag*pmag))
        tnode = orbang*pmag/vmag + time(1) + timref(3)

      end if
      return
      end       
