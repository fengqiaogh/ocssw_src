      subroutine find_islands(prod_ID,npix,west,numi,xlon,xlat,
     *  wlon,wlat,nump,nill,nilp)

      save
      real*4 xlon(1000),xlat(1000),wlon(1000),wlat(1000)
      real*4 nav_orb_vec(3), nav_sun_ref(3), nav_sen_mat(3,3), 
     1     nav_scan_ell(6)
      real*4 plat(1285),plon(1285),sena(1285),senz(1285)
      real*4 sola(1285),solz(1285)
      real*4 b2(1285),b8(1285),b2p(1285),b82p(1285)
      integer*2 dat(8,1285)
      integer*4 ill(30000),ilp(2,30000),ilf(30000)
      integer*4 jll(1000),jlp(2,1000)
      integer*4 irec,prod_ID(2)
      integer*4 nump(1000),nill(1000),nilp(1000)
      byte lflag(1285,16000)
      logical land,end,first,west
      data b82_ice_min/20.0/
      data b2_cloud_min/18./
c      data b2_cloud_stray_min/25./
c      data b82_cloud_stray_max/-13./
      data b2_valid_min/6.0/
      data b82_valid_min/-20.0/
      data b82_valid_max/8.0/
      data b82_water_max/-3.0/
      data b82_fac/0.7/
      data a2/0.02/,a8/0.15/,tc/1400./
      data lp_maxg/15/,lp_maxl/40/
      data first/.true./

      logical(4) water
      real*8 pi,radeg,re,rem,f,omf2,omegae
      common /gconst/pi,radeg,re,rem,f,omf2,omegae

#include "navblk_s.fin"

      type(navblk_struct) :: navblk(16000)
  
      lp_max = lp_maxl
      ip = 1
      jp = 1
      ipst = 147
      ipen = 1135
      if (npix.eq.248) then
        lp_max = lp_maxg
        ip = 147
        jp = 4
        ipst = 1
        ipen = 248
      end if

      tr2 = a2*exp(-prod_id(2)/tc) + (1. - a2)
      tr8 = a8*exp(-prod_id(2)/tc) + (1. - a8)

      iline = 0
      nsegs = 0
      irec = 0
      ntpix = 0
      ncpix = 0
      nipix = 0
      nwpix = 0
      nlpix = 0
      nbpix = 0
      numi = 0

c  First read in data and find scan line land segments using Bands 1 and 6
      dowhile (.true.)
        lw = -1
        land = .false.
      call get_l1a_record( prod_ID, npix, irec, dat, nav_orb_vec, 
     1     nav_sun_ref, nav_sen_mat, nav_scan_ell, iret)

c      write(6, 500 ) iret, ((rec_arr(i,j),i=1,8 ),j=71,72),
c     1     (nav_orb_vec(i),i=1,3), (nav_sun_ref(i),i=1,3), 
c     1     (nav_sen_mat(i),i=1,9), (nav_scan_ell(i),i=1,6)
c  500 format( ' test, read: iret = ', i7,'  8 bands for 2 pixels:',/,
c     1     2(8(2x,i5),/),' nav_orb_vec: ',3e12.5,/,
c     1     ' nav_sun_ref: ',3e12.5,/,' nav_sen_mat: ',/,
c     1     3(3(2x,e12.5),/),' nav_scan_ell: ',/,
c     1     2(3(2x,e12.5),/) )

        if (iret.ne.0) go to 990

        iline = iline + 1

c  Store navigation data
        do i=1,3
          navblk(iline)%orb_vec(i) = nav_orb_vec(i)
          navblk(iline)%sun_ref(i) = nav_sun_ref(i)
          navblk(iline)%scan_ell(i) = nav_scan_ell(i)
          navblk(iline)%scan_ell(i+3) = nav_scan_ell(i+3)
          navblk(iline)%sen_mat(i,1) = nav_sen_mat(i,1)
          navblk(iline)%sen_mat(i,2) = nav_sen_mat(i,2)
          navblk(iline)%sen_mat(i,3) = nav_sen_mat(i,3)
        end do

c  Check if scene starts in Western hemisphere
        if (first) then
          first = .false.
          west = .false.
          ik = ip + jp*(npix-1)
          call geonav(navblk(iline)%orb_vec,navblk(iline)%sen_mat,
     *      navblk(iline)%scan_ell,navblk(iline)%sun_ref,ik,1,1,
     *      plat,plon,solz,sola,senz,sena)
          if (plon(1).lt.0.0) west = .true.
        end if

c  Get solar and sensor zenith angles
        call geonav(navblk(iline)%orb_vec,navblk(iline)%sen_mat,
     *    navblk(iline)%scan_ell,navblk(iline)%sun_ref,ip,jp,npix,
     *    plat,plon,solz,sola,senz,sena)

c       csz = cos(solz(1)/radeg)

c  Linearize counts
        call b28_lin(dat,npix,b2,b8)
        
c       write (71) (b2(i),i=1,npix),(b8(i),i=1,npix)

        do i=ipst,ipen

          csz = cos(solz(i)/radeg)
          cnz = cos(senz(i)/radeg)
          ssz = sin(solz(i)/radeg)
          snz = sin(senz(i)/radeg)
          ca = cos((sena(i)-sola(i))/radeg)

          fac2 = 2./(1.+csz/cnz)
c          fac8 = (csz*cnz + ssz*snz*ca)/csz
          b2(i)=b2(i)/tr2
          b8(i)=b8(i)/tr8
          b2p(i) = b2(i)/csz
c         b82p(i) = b8(i)*fac8-b82_fac*b2(i)*fac2
          b82p(i) = b8(i)*fac2-b82_fac*b2(i)*fac2
          ntpix = ntpix + 1

c  First check for invalid (corrupted) radiances
          if ((b2p(i).lt.b2_valid_min).or.
     *         (b82p(i).lt.b82_valid_min).or.
     *         (b82p(i).gt.b82_valid_max)) then
            land = .false.
            lw = -1
            lflag(i,iline) = 2
            nbpix = nbpix + 1

c  If ice, reset flags
c         if (b82p(i).gt.b82_ice_min) then
c           land = .false.
c           lw = -1
c           lflag(i,iline) = 2
c           nipix = nipix + 1

c  If clouds, reset flags
          else if (b2p(i).gt.b2_cloud_min) then
            land = .false.
            lw = -1
            lflag(i,iline) = 3
            ncpix = ncpix + 1

c  Check for stray light from clouds
c         else if ((b1p(i).gt.b1_cloud_stray_min).and.
c     *         (b61p(i).lt.b61_cloud_stray_max)) then
c           land = .false.
c           lw = -1
c           lflag(i,iline) = 3
c           ncpix = ncpix + 1

c  If water 
          else if (b82p(i).lt.b82_water_max) then
            water = .true.
            lflag(i,iline) = 0
            nwpix = nwpix + 1

c  If previous pixel is land 
            if (land) then

c  If land is preceded by water within maximum number of pixels 
              if ((lw.gt.0).and.((i-lw).lt.lp_max)) then
c               write(*,*)'Land at line, pixels',iline,lw,i

c  Store segment indices in table
                nsegs = nsegs + 1
                ill(nsegs) = iline
                ilp(1,nsegs) = lw + 1
                ilp(2,nsegs) = i - 1
                ilf(nsegs) = 0
              end if      

              land = .false.
            end if              
            lw = i

c  Else if land
          else 
            land = .true.
            lflag(i,iline) = 1
            nlpix = nlpix + 1

          end if
        end do
c       write (61) (lflag(i,iline),i=1,npix)
c       write (81) (b2p(i),i=1,npix),(b82p(i),i=1,npix)
      end do
  990 continue

      nipix = 100*nipix/ntpix
      ncpix = 100*ncpix/ntpix
      nwpix = 100*nwpix/ntpix
      nlpix = 100*nlpix/ntpix
      nbpix = 100*nbpix/ntpix
c      write(*,*)'Ice pixels  ',nipix,'%'
      write(*,*)'Cloud pixels',ncpix,'%'
      write(*,*)'Water pixels',nwpix,'%'
      write(*,*)'Land pixels ',nlpix,'%'
      write(*,*)'Corrupted pixels ',nbpix,'%'

      if (nsegs.lt.1) go to 998


c  Now try to identify islands in land segments by 
c   chaining segments if necessary

      il = 1

c  Check if segments are in first or last scan line
      dowhile(ill(il).lt.2)
        il = il + 1
        if (il.gt.nsegs) go to 998
      end do

      dowhile (ill(nsegs).ge.iline)
        nsegs = nsegs - 1
        if (nsegs.lt.il) go to 998
      end do

      dowhile (il.le.nsegs)     

        js = 0
c  Check next segment in table

        if (ilf(il).ne.0) then
          go to 997
        end if
        ilf(il) = 1

c   First check previous line for land or cloud pixels
        do j=ilp(1,il)-1,ilp(2,il)+1
          if (lflag(j,ill(il)-1).ge.1) then
            go to 997
          end if
        end do

c  Assume for the moment this is the first segment of a chain
        i1 = il
        jll(1) = ill(il)
        jlp(1,1) = ilp(1,il)
        jlp(2,1) = ilp(2,il)
        end = .false.

c  Look for contiguous land segments
        call find_all_segs(ill,ilp,ilf,i1,nsegs,jll(1),jlp(1,1),nj)

        if (nj.gt.1000) go to 997
c  If contiguous segments are found then
        if (nj.gt.1) then

c   Consolidate segment table
          nl = 1
          do i=2,nj
            jl = jll(i) - jll(1) + 1
            if (jl.gt.nl) then
              nl = jl
              jll(nl) = jll(i)
              jlp(1,nl) = jlp(1,i)
              jlp(2,nl) = jlp(2,i)
            else
              if (jlp(1,i).lt.jlp(1,jl)) jlp(1,jl) = jlp(1,i) 
              if (jlp(2,i).gt.jlp(2,jl)) jlp(2,jl) = jlp(2,i)
            end if
          end do              
               
c   Check previous line for land or cloud pixels
          do j=jlp(1,1)-1,jlp(2,1)+1
            if (lflag(j,jll(1)-1).ge.1) go to 997
          end do

c  Check for other contiguous land pixels
          do j=1,nl-1
            call check_segs(jll,jlp,j,j+1,lflag(1,jll(j)),
     *          lflag(1,jll(j)+1),igood)
            if (igood.ne.0) go to 997       
          enddo

        else
          nl = 1
        end if

c  Check after last line for land or cloud pixels
        do j=jlp(1,nl)-1,jlp(2,nl)+1
          if (lflag(j,jll(nl)+1).ge.1) go to 997
        end do

c       write(*,*)'Land segment chain'
        numi = numi + 1

c  Compute centroid of land pixels
        ni = 0
        xlati = 0.0
        xloni = 0.0
        xlomin = 999.
        xlomax = -999.
        xlamin = 999.
        xlamax = -999.
        do i=1,nl
c         write(*,*)jll(i),jlp(1,i),jlp(2,i)
          k = jll(i)
          ik = ip + jp*(jlp(1,i)-1)
          np = jlp(2,i) - jlp(1,i) + 1
          call geonav(navblk(k)%orb_vec,navblk(k)%sen_mat,
     *          navblk(k)%scan_ell,navblk(k)%sun_ref,ik,jp,np,
     *          plat,plon,solz,sola,senz,sena)
          do j=1,np
            if (lflag((j+jlp(1,i)-1),k).eq.1) then
              xlati = xlati + plat(j)
              if (west.and.(plon(j).lt.0.0)) plon(j) = plon(j) + 360.0
              xloni = xloni + plon(j)
              ni = ni + 1
              if (plon(j).lt.xlomin) xlomin = plon(j)
              if (plon(j).gt.xlomax) xlomax = plon(j)
              if (plat(j).lt.xlamin) xlamin = plat(j)
              if (plat(j).gt.xlamax) xlamax = plat(j)
            end if
          end do                
        end do
        nump(numi) = ni
        nill(numi) = jll(1)
        nilp(numi) = jlp(1,1)
        xlat(numi) = xlati/ni
        xlon(numi) = xloni/ni
        wlat(numi) = xlamax - xlamin
        wlon(numi) = xlomax - xlomin
        write(*,*)ni,xlat(numi),xlon(numi),wlat(numi),wlon(numi)
        write (51,1100) ni,xlat(numi),xlon(numi),jll(1),jlp(1,1),
     *    wlat(numi),wlon(numi)
 1100   format (i10,2f12.5,2i10,2f10.5)

  997   il = il + 1

        if (numi.ge.1000) go to 998

      end do

  998 return
      end
