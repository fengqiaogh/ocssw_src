      subroutine czcscal(iorbit,slopech,yintch,icnts,igain,rlt)
c
c  Calibrates CZCS radiances according to Evans and Gordon, 1994.
c
c     Name        Type    I/O   Description
c     ----------  ----    ---   --------------------------------
c     iorbit       int      I   orbit count
c     slopech     real      I   size 5 slope data for 5 vis bands
c     yintch      real      I   size 5  intercept for 5 vis bands
c     icnts      int*2      I   raw radiance counts for 5 vis bands
c     igain        int      O   derived gain
c     rlt         real      O   size 5 calibrated radiances
c
c     From Watson Gregg, modified be W. Robinson, 29 Apr 2004
c
c
c     rlt - O calibrated radiances for 5 vis bands
      save ifst
      parameter(npix=1968,nlt=5,nltc=4)
      integer*2 icnts(nlt)
      real slope(nltc,4),yint(nltc,4),rk(nltc,4)
      real slopech(nlt),yintch(nlt)
      real g(nltc)
      integer iorbdec(5)
      real dec(nltc,5)
      real S(nlt),rlt(nlt)
      data ifst /0/
c
c                 Band 1  Band 2  Band 3  Band 4
      data slope /0.04452,0.03103,0.02467,0.01136,   !Gain 1
     *            0.03598,0.02493,0.02015,0.00897,   !Gain 2
     *            0.02968,0.02032,0.01643,0.00741,   !Gain 3
     *            0.02113,0.01486,0.01181,0.00535/   !Gain 4
      data yint /0.03963,0.05361,0.08013,0.01136,
     *           0.03963,0.05361,0.08013,0.01136,
     *           0.03963,0.06461,0.09503,0.01136,
     *           0.03963,0.06361,0.09159,0.01136/
      data rk /1.018,0.982,0.974,1.008,
     *         1.021,0.983,0.963,1.020,
     *         1.011,0.988,0.947,1.016,
     *         1.020,0.972,0.950,1.010/
      data iorbdec /1,5000,7000,20000,39000/
      data dec /1.0, 1.0, 1.0, 1.0,
     *          0.92,0.98,0.99,1.0,
     *          0.84,0.97,0.99,1.0,
     *          0.73,0.86,0.91,1.0,
     *          0.59,0.78,0.84,0.91/
c
      if (ifst .eq. 0)then
c  Find gain
       call gain(slopech,slope,igain)
c       write(6,*)'Gain = ',igain
c  Find decay
c       write(6,*)'iorbit = ',iorbit
       iorbit = min(iorbit,39000)
       do i = 1,4
        if (iorbit .ge. iorbdec(i) .and. iorbit .lt. iorbdec(i+1))then
         idec = i
        endif
       enddo
c       write(6,*)'idec = ',idec
       i = idec
       x1 = float(iorbdec(i))
       x2 = float(iorbdec(i+1))
       rden = x2-x1
       do nl = 1,nltc
        rnum = dec(nl,i+1) - dec(nl,i)
        slp = rnum/rden
        gm1 = float(iorbit-iorbdec(i))*slp + dec(nl,i)
        g(nl) = 1.0/gm1
       enddo
c  Adjust calibration of Band 4 according to Evans and Antoine,
c  personal communication
       iorb4 = 6750
       if (iorbit .gt. iorb4)then
        x1 = float(iorb4)
        x2 = float(iorbdec(5))
        rden = x2-x1
        nl = nltc
        rnum = dec(nl,5) - dec(nl,4)
        slp = rnum/rden
        gm1 = float(iorbit-iorb4)*slp + dec(nl,4)
        g(nl) = 1.0/gm1
       endif
       do nl = 1,nltc
c        write(6,*)'i,nl,dec,dec,g = ',i,nl,dec(nl,i),dec(nl,i+1),g(nl)
       enddo
       S(nlt) = 0.0
       ifst = 1
       slope1 = slopech(1)
      endif
c
c  Slope
      if (slopech(1) .ne. slope1)call gain(slopech,slope,igain)
      ig = igain
      do nl = 1,nltc
       S(nl) = g(nl)*rk(nl,ig)*slope(nl,ig)
       rlt(nl) = S(nl)*float(icnts(nl)) + yint(nl,ig)
c       write(6,*)'nl,icnts,rlt = ',nl,icnts(nl),rlt(nl)
      enddo
      rlt(nlt) = slopech(nlt)*icnts(nlt) + yintch(nlt)
c
      return
      end
c
c ****************************************************************
      subroutine gain(slopech,slope,igain)
c
c  Computes gain from tabulated slopes given in the data file.
c
      parameter(nlt=5,nltc=4)
      real slope(nltc,4)
      real slopech(nlt)
c
      dif1 = 99999.0
      do ig = 1,4
       dif = abs(slopech(1) - slope(1,ig))
       if (dif .lt. dif1)then
        dif1 = dif
        igain = ig
       endif
      enddo
c
      return
      end
