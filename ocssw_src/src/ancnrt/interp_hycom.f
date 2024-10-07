      program interp_hycom

      use netcdf
      USE HDF5  
c
c     This program interpolate a global 1/12 degree hycom temperature
c     or salinity field at selected depth into a regular grid of
c     a quarter degree (1440x720).  
c
c     This program is based on isubaregion.f written by Alan J. Wallcraft
c     It takes two pre-generated NetCDF files:  hycom_weight.nc
c     that stores the interpolation weight table from the hycom grid
c     to the regular grid based on a bilinear interpolation
c     The weight table was generated using a modified program based 
c     on isuba_gmapi.f also written by Wallcraft.  
c
c     The second file is mask1.nc.  It provides the landmask for the
c     output grid.  This two files have to reside in the directory
c     where the program is executed.
c
c     The input file is a NetCDF file extracted from the HYCOM output
c     file via their OpenDAP server.  The variable field contains only
c     one depth only.  If you change the input file format, you will
c     have to modify the code to read it properly.
c
c     The output file is a NetCDF file containing 1D arrays for latitude
c     and longitude, two 3D arrays for temp and salt.  There are 
c     9 pre-defined depth at 0, 10, 20, 30, 50, 100, 250, 500, 1000 meter
c     below the sea surface.
c     
c     Original Author:
c     Alan J. Wallcraft,  NRL,  January 2005 and January 2009.c     
c     Modified by
c     Peggy Li, JPL, August 31, 2010
c     Modified by
c     Joel Gales, GSFC, September 5, 2010
c     Add HDF5 and binary output
c
c
      character*40         :: field, timetag, depth, varname
      character*128        :: dir_in, dir_out
      character*256        :: file_in, file_out
      integer              :: nargs,iargc
      integer              :: level
      integer              :: idm, jdm,idm_out,jdm_out,ddim
      integer              :: i,ii,ip,j,jj,jp
      integer              :: if_sm,il_sm,jf_sm,jl_sm
      integer              :: status,ncid,ncid2,vid
      integer              :: latid, lonid, depthid
      integer              :: varid, xoutid, youtid
      real, allocatable    :: lats(:), lons(:), dep(:)
      character*80 :: attstr
      integer              :: start_3d(4),count_3d(4)
!      logical              :: lmask
      integer, allocatable :: m_sm(:,:),  iv_sm(:,:)
      integer, allocatable :: m_out(:,:), m_osm(:,:)
      real,    allocatable :: a_in(:,:)
      real,    allocatable, target :: a_out(:,:)
c   
      integer, allocatable :: i_out(:,:),j_out(:,:)
      real,    allocatable :: x_out(:,:),y_out(:,:)
c
      real,    parameter   :: spval=2.0**100  ! spval
      real,    parameter   :: hspval=0.5*2.0**100  ! half spval
!      DOUBLE PRECISION                :: time
      logical              :: file_exists
      integer              :: count, count1

      INTEGER(HID_T) :: file_id   
      INTEGER(HID_T) :: dset_id
      INTEGER :: rank
      INTEGER(HSIZE_T), DIMENSION(2) :: dims_sal
      INTEGER(HSIZE_T), DIMENSION(3) :: offset=(/0,0,0/)
      INTEGER(HID_T) :: filespace
      INTEGER(HID_T) :: dataspace
      real(4), POINTER :: pt

      integer slash_pos
      CHARACTER datadir*255

      include "aquarius.fin"
c
c      define input/output dimension
c      parse input arguments
      idm=4500;jdm=3298
      idm_out=1440;jdm_out=720
      nargs=iargc()
      if (nargs.eq.2) then
!          call getarg(1,field)
!          call getarg(2,timetag)
!          call getarg(2,depth)
         call getarg(1,file_in)
         call getarg(2,dir_out)
      else
         print *,
!     &    "Command:get_hycom [temp|salt] yyyymmdd depth infile outdir"
     &        "Command:get_hycom [temp|salt] infile outdir"
!          print *,"   depth: 0,10,20,30,50,100,250,500,and 1000"
!         stop
         call exit(1)
      endif
      field = 'salt'
      depth = '0'

      slash_pos = INDEX(file_in,'/',BACK=.TRUE.)
      timetag = file_in(slash_pos+7:slash_pos+10)//file_in(slash_pos+12:slash_pos+14)

      if (field .eq. "temp") then 
         varname="temperature"
      elseif (field .eq. "salt") then
         varname="salinity"
      else
         print *, "field has to be either temp or salt"
!         stop
         call exit(1)
      endif
      level = -1
      if (depth .eq. '0') level=1
      if (depth .eq. '10') level=2
      if (depth .eq. '20') level=3
      if (depth .eq. '30') level=4
      if (depth .eq. '50') level=5
      if (depth .eq. '100') level=6
      if (depth .eq. '250') level=7
      if (depth .eq. '500') level=8
      if (depth .eq. '1000') level=9
      if (level < 0) then
         print *, "depth value is invalid"
!         stop
         call exit(1)
      endif

      allocate(  iv_sm(jdm,2) )
c
      allocate(   m_sm(idm,jdm),     m_osm(idm_out,jdm_out) )
      allocate(   m_out(idm_out,jdm_out) )
      allocate(   a_in(idm,jdm),     a_out(idm_out,jdm_out) )
c                                                   
      allocate( i_out(idm_out,jdm_out), j_out(idm_out,jdm_out) )
      allocate( x_out(idm_out,jdm_out), y_out(idm_out,jdm_out) )

      a_out(:,:)=spval
!      write(*,*) 'interpolating ',trim(file_in),' to ',
!     &  trim(file_out)

      call getenv('OCDATAROOT', datadir)
      
c     readin the land mask for the output grid
      status=nf90_open(trim(datadir)//'/aquarius/static/hycom_landmask.nc', 0, ncid)    
      if (status /= nf90_noerr) then
         write(*,*) 'Land mask "hycom_landmask.nc" cannot be opened.'
         call handle_err(status)
      endif

      status=nf90_inq_varid(ncid,"mask",vid)
      if (status /= nf90_noerr) then
         call handle_err(status)
      endif

      status=nf90_get_var(ncid,vid,m_out) 
      if (status /= nf90_noerr) then
         call handle_err(status)
      endif

      status=nf90_close(ncid)


c     stop if input file does not exit
c     and create the output file if not exist
      inquire(FILE=file_in, EXIST=file_exists)
      if (.not. file_exists) then
         print *, file_in, " does not exist"
!         stop
         call exit(1)
      endif
!      inquire(FILE=file_out, EXIST=file_exists)
!      if (.not. file_exists) then

      file_out=trim(dir_out)//'/N'//trim(timetag)//'00_SALINITY_HYCOM_24h.nc'

!      print *, "Create ",trim(file_out)
      ddim = 9

      status=nf90_create(file_out,0, ncid)
      if (status /= nf90_noerr) then
         write(*,*) 'NetCDF output file cannot be created.'
         call handle_err(status)
      endif

      status=nf90_def_dim(ncid,'lat', jdm_out, latid)
      if (status /= nf90_noerr) call handle_err(status)
      status=nf90_def_dim(ncid,'lon', idm_out, lonid)
      if (status /= nf90_noerr) call handle_err(status)

      status=nf90_def_var(ncid, 'salt', NF90_FLOAT, 
     &     (/lonid,latid/), varid)
      if (status /= nf90_noerr) call handle_err(status)
      attstr="Salinity"
      status=nf90_put_att(ncid,varid,'long_name',attstr)
      attstr="PS"
      status=nf90_put_att(ncid,varid,'units',attstr)
      if (status /= nf90_noerr) call handle_err(status)
      status=nf90_put_att(ncid,varid,'_FillValue',spval)
      if (status /= nf90_noerr) call handle_err(status)
      status=nf90_enddef(ncid)
      if (status /= nf90_noerr) call handle_err(status)

      allocate(lats(jdm_out), lons(idm_out), dep(ddim))
      
      do i=1,idm_out
         lons(i)=-180.125+i*0.25
      enddo
      do i=1,jdm_out
         lats(i)=-90.125+i*0.25
      enddo
      dep(1)=0; dep(2)=10; dep(3)=20; dep(4)=30
      dep(5)=50; dep(6)=100; dep(7)=250; dep(8)=500
      dep(9)=1000
      
      deallocate(lats,lons,dep)


c     Open input file and read in the field 
      status=nf90_open (file_in, nf90_nowrite, ncid)
      if (status /= nf90_noerr) call handle_err(status)
      status=nf90_inq_varid(ncid,varname,vid)
      status=nf90_get_var(ncid,vid,a_in)
      if (status /= nf90_noerr) call handle_err(status)
      status=nf90_close(ncid)
      if (status /= nf90_noerr) call handle_err(status)

      count = 0
      do j= 1,jdm
        do i= 1,idm
          if     (a_in(i,j).gt.hspval) then
            count = count+1
            m_sm(i,j) = 0
          else
            m_sm(i,j) = 1
          endif
        enddo
      enddo
      print *, "hspval", hspval
      print *, "Input data ", count, " land point"

c  read in the weight file and output mask

       status=nf90_open(trim(datadir)//'/aquarius/static/hycom_weight.nc',0, ncid)

       if (status /= nf90_noerr) then
          write(*,*) 'HYCOM weight file "hycom_weight.nc" cannot be opened.'
          call handle_err(status)
       endif

       status=nf90_inq_varid(ncid, 'x_out', xoutid)
       if (status /= nf90_noerr) call handle_err(status)
       status=nf90_inq_varid(ncid, 'y_out', youtid)
       if (status /= nf90_noerr) call handle_err(status)
       status=nf90_get_var(ncid,xoutid,x_out)     
       if (status /= nf90_noerr) call handle_err(status)
       status=nf90_get_var(ncid,youtid,y_out)     
       if (status /= nf90_noerr) call handle_err(status)
       status=nf90_close(ncid)
       if (status /= nf90_noerr) call handle_err(status)

c      m_out(:,:)=1
      do jj= 1,jdm_out
        do ii= 1,idm_out
          if     (x_out(ii,jj).lt.hspval) then
            i_out(ii,jj) = int(x_out(ii,jj))
            x_out(ii,jj) =     x_out(ii,jj) - i_out(ii,jj)
            j_out(ii,jj) = int(y_out(ii,jj))
            y_out(ii,jj) =     y_out(ii,jj) - j_out(ii,jj)
          else  !update output mask???
            m_out(ii,jj) = 0
          endif
        enddo !ii
      enddo !jj
c
c --- form the p-grid smoother mask,
c --- set to 2 if a land point is needed for interpolation.
c --- we are assuming that "2" is never needed outside the
c --- target subregion, which will be the case unless the
c --- subregion rectangle is poorly chosen.
c

      count = 0
      do jj= 1,jdm_out
        do ii= 1,idm_out
          count = 0
          if     (m_out(ii,jj).eq.1) then
            j  = j_out(ii,jj)
            i  = i_out(ii,jj)
            jp = min(j+1,jdm)
            if     (i.ne.idm) then
              ip = i+1            
            else  !lperiod        
              ip =   1    
            endif         
            if     (m_sm(i,j).eq.0) then
              count = count+1
              m_sm(i, j ) = 2
            endif
            if     (m_sm(i,jp).eq.0) then
              count = count+1
              m_sm(i, jp) = 2
            endif
            if     (m_sm(ip,j).eq.0) then
              count = count+1
              m_sm(ip,j) = 2
            endif
            if     (m_sm(ip,jp).eq.0) then
              count = count+1
              m_sm(ip,jp) = 2
            endif
            if (count > 0) count1=count1+1
          endif !sea-point
        enddo !ii
      enddo !jj
c
      print *, "Total ", count1, "ocean point hits land"

      do j= 1,jdm
        iv_sm(j,1) = idm
        do i= 1,idm
          if     (m_sm(i,j).eq.2) then
            iv_sm(j,1) = i
            exit
          endif
        enddo
        iv_sm(j,2) = 1
        do i= idm,iv_sm(j,1),-1
          if     (m_sm(i,j).eq.2) then
            iv_sm(j,2) = i
            exit
          endif
        enddo
      enddo
      jf_sm = jdm
      do j= 1,jdm
        if     (iv_sm(j,1).le.iv_sm(j,2)) then
          jf_sm = j
          exit
        endif
      enddo
      jl_sm = 1
      do j= jdm,jf_sm,-1
        if     (iv_sm(j,1).le.iv_sm(j,2)) then
          jl_sm = j
          exit
        endif
      enddo
      if_sm = minval(iv_sm(:,1))
      il_sm = maxval(iv_sm(:,2))
        
        call landfill(  a_in,m_sm,idm,    jdm,
     &                  iv_sm,if_sm,il_sm,jf_sm,jl_sm)
        call bilinear_p(a_in,     idm,    jdm,
     &                  a_out,    idm_out,jdm_out,
     &                  m_sm, m_out,i_out,j_out,x_out,y_out)


c   output the interpolated field
      status = nf90_open(file_out,nf90_write,ncid2)
      if (status /= nf90_noerr) call handle_err(status)
      start_3d(1)=1
      start_3d(2)=1
      start_3d(3)=level
      count_3d(1)=idm_out
      count_3d(2)=jdm_out
      count_3d(3)=1
      status=nf90_inq_varid(ncid2,field,vid)
      if (status /= nf90_noerr) call handle_err(status)
      status=nf90_put_var(ncid2, vid,a_out, 
     &     start=start_3d,count=count_3d)
      if (status /= nf90_noerr) call handle_err(status)
      status=nf90_close(ncid2)
      if (status /= nf90_noerr) call handle_err(status)


c     write out a land mask based on the output
      do j=1,jdm_out
        do i=1,idm_out
          if (a_out(i,j) .eq. spval) then
             m_out(i,j)=0
          else
             m_out(i,j)=1
          endif
        enddo
      enddo

C     Write out HD5 file
      file_out=trim(dir_out)//'/N'//trim(timetag)//'00_SALINITY_HYCOM_24h.h5'

      call create_hdf5(file_out, file_id)

      rank = 2
      dims_sal(1) = 360*4
      dims_sal(2) = 180*4

      call create_hdf5_real(file_id, 'salinity', rank, dims_sal, dset_id)
      call set_space_hdf5(dset_id, rank, dims_sal, filespace, dataspace)
      pt => a_out(1,1)
      call write_hdf5_real(dset_id, filespace, dataspace, offset, dims_sal, pt)
      call close_hdf5_ds(dset_id, dataspace, filespace)
      call close_hdf5_df(file_id)


C     Write binary file
      file_out=trim(dir_out)//'/N'//trim(timetag)//'00_SALINITY_HYCOM_24h.dat'
      open(unit=8,file=file_out, form='UNFORMATTED')
      write(8) a_out
      close(8)

      deallocate(iv_sm, m_sm, m_osm, m_out, a_in, a_out)
      deallocate(i_out, j_out, x_out, y_out)

      write(*,*) ' '
      write(*,*) 'NETCDF=N'//trim(timetag)//'00_SALINITY_HYCOM_24h.nc'
      write(*,*) 'HDF5  =N'//trim(timetag)//'00_SALINITY_HYCOM_24h.h5'
      write(*,*) 'BINARY=N'//trim(timetag)//'00_SALINITY_HYCOM_24h.dat'

      end program interp_hycom

      subroutine bilinear_p(a_in, idm_in, jdm_in,
     &                      a_out,idm_out,jdm_out,
     &                      m_in, m_out,i_out,j_out,x_out,y_out)
      implicit none
c
      integer idm_in, jdm_in,
     &        idm_out,jdm_out
      integer m_out(idm_out,jdm_out),
     &        m_in(idm_in,jdm_in),
     &        i_out(idm_out,jdm_out),
     &        j_out(idm_out,jdm_out)
      real    a_in( idm_in, jdm_in ),
     &        a_out(idm_out,jdm_out),
     &        x_out(idm_out,jdm_out),
     &        y_out(idm_out,jdm_out)
c
c --- interpolate from a_in to a_out.
c
      integer i,ii,ip,j,jj,jp
      integer ki,kj,iii,jjj
      real sa, ss
      real    sx,sy
      integer   count
      real      s(-1:1, -1:1)
      data      s / 1.0, 2.0, 1.0,
     &              2.0, 4.0, 2.0,
     &              1.0, 2.0, 1.0 /
      real,    parameter   :: spval=2.0**100  ! spval

c
      count = 0
      do jj= 1,jdm_out
        do ii= 1,idm_out
          if     (m_out(ii,jj).eq.1) then
            count = count+1
            sx = x_out(ii,jj)
            sy = y_out(ii,jj)
            i  = i_out(ii,jj)
            if     (i.ne.idm_in) then
              ip = i+1
            else
              ip =   1
            endif
            j  = j_out(ii,jj)
            jp = j+1
c
c   doing a 3x3 interpolation for jj<552
            if (jj<552) then
              sa = 0.0
              ss = 0.0
              do kj= -1,1
                jjj = j+kj
                if (jjj==0) jjj=1
                if (jjj>jdm_in) jjj=jdm_in
                do ki= -1,1
                  iii = i+ki
                  if (iii==0) iii=idm_in
                  if (iii>idm_in) iii=1
                     if (m_in(iii,jjj).eq.1) then
                       sa = sa + s(ki,kj)*a_in(iii,jjj)
                       ss = ss + s(ki,kj)
                     endif
                enddo
              enddo
              if     (ss.ne.0.0) then
c
c               at least one ocean point within stencil.
c
                a_out( ii,jj)     = sa/ss
              else
                a_out(ii,jj) = spval
              endif
            else
            a_out(ii,jj) = (1.0-sx)*(1.0-sy)*a_in(i, j ) +
     &                     (1.0-sx)*     sy *a_in(i, jp) +
     &                          sx *(1.0-sy)*a_in(ip,j ) +
     &                          sx *     sy *a_in(ip,jp)
*           if     (a_out(ii,jj).lt.0.0) then
*             write(6,'(a,6i5,2f7.3,5f9.2)')
*    &        'ii,jj,i,ip,j,jp,sx,sy,a_out,a_in',
*    &         ii,jj,i,ip,j,jp,
*    &         sx,sy,a_out(ii,jj),
*    &         a_in(i, j ), 
*    &         a_in(i, jp), 
*    &         a_in(ip,j ), 
*    &         a_in(ip,jp)
*           endif
            endif
          endif
        enddo
      enddo
      print *, "Output ocean point:", count
      return
      end

      subroutine landfill(a,mask,m,n, iv,if,il,jf,jl)
      implicit none
c
      integer m,n,mask(m,n), iv(n,2),if,il,jf,jl
      real    a(m,n)
c
c --- extrapolate a 1-grid cell into the land mask,
c ---   mask == 0 for land.  
c ---   mask == 1 for ocean.
c ---   mask == 2 for land to be extrapolated to ocean.
c
      integer, allocatable :: mm(:,:,:)
c
      integer i,ii,ip0,ip1,ipass,j,jj,ki,kj,nleft,nup
      real    sa,ss
c
      logical lfirst
c      real    s(-1:1,-1:1)
      real   s(-2:2,-2:2)
      save    lfirst,s
c
      data lfirst / .true. /
c      data      s / 1.0, 2.0, 1.0,
c     &              2.0, 4.0, 2.0,
c     &              1.0, 2.0, 1.0 /
      data   s /1.0, 1.5, 2.0, 1.5, 1.0,
     &          1.5, 2.0, 3.0, 2.0, 1.5,
     &          2.0, 3.0, 4.0, 3.0, 2.0, 
     &          1.5, 2.0, 3.0, 2.0, 1.5,
     &          1.0, 1.5, 2.0, 1.5, 1.0 /
c
c     adding a halo to mm simplifies ocean selection logic.
c
      allocate( mm(0:m+1,0:n+1,0:1) )
c
      mm( : , : ,0) = 0
      mm(1:m,1:n,0) = mask
      mm( : , : ,1) = mm(:,:,0)
c
c --- repeated passes of 9-point "smoother" to
c ---  convert all mask==2 points to mask==1.
c --- double-buffering mm allows in-place use of a.
c
      if     (lfirst) then
        write(6,'(/a,6i5/)')
     &    'landfill - m,n,if,il,jf,jl =',m,n,if,il,jf,jl
      endif
      do ipass= 1,n+m
        ip0   = mod(ipass+1,2)
        ip1   = mod(ipass,  2)
        nup   = 0
        nleft = 0
        do j= jf,jl
          do i= iv(j,1),iv(j,2)
            if     (mm(i,j,ip0).eq.2) then
              sa = 0.0
              ss = 0.0
              do kj= -2,2
                jj = j+kj
                do ki= -2,2
                  ii = i+ki
                  if ((ii>0) .and. (ii.le.m) .and. (jj>0) 
     &                   .and. (jj.le.n)) then
                     if (mm(ii,jj,ip0).eq.1) then
                       sa = sa + s(ki,kj)*a(ii,jj)
                       ss = ss + s(ki,kj)
                     endif
                  endif
                enddo
              enddo
              if     (ss.ne.0.0) then
c
c               at least one ocean point within stencil.
c
                a( i,j)     = sa/ss
                mm(i,j,ip1) = 1
                nup         = nup + 1
*               if     (mask(i,j).eq.1) then
*                 write(6,*) 'error - i,j,ip0,ip1,mask,mm = ',
*    &              i,j,ip0,ip1,mask(i,j),mm(i,j,ip0)
*                 stop
*               endif
*               if     (mod(nup,1000).eq.1) then
*                 write(6,'(a,2i5,f5.1,f10.3)') 
*    &              '   i,j,ss,a = ',i,j,ss,a(i,j)
*               endif
              else
                nleft = nleft + 1
              endif
            endif
          enddo
        enddo
        if     (lfirst) then
          write(6,'(a,i4,a,i6,a,i6,a)')
     &      'landfill: pass',ipass,
     &      ' filled in',nup,
     &      ' points, with',nleft,' still to fill'
          call flush(6)
        endif
        if     (nup.eq.0) then
          exit
        endif
        mm(if:il,jf:jl,ip0) = mm(if:il,jf:jl,ip1)
      enddo  ! ipass=1,...
      if     (lfirst) then
        write(6,*)
        lfirst = .false.
      endif
      if     (nleft.ne.0) then
        write(6,'(/a,i6,a/a/)')
     &    'error in landfill - ',
     &    nleft,' "mask==2" values are not fillable',
     &    'probably a mismatch between input and output land masks'
        call flush(6)
c        do j= jf,jl
c          do i= iv(j,1),iv(j,2)
c            if     (mm(i,j,ip1).eq.2) then
c              write(6,'(a,2i5)') 'mask==2 at (input) i,j = ',i,j
c            endif
c          enddo
c        enddo
c        write(6,*)
c        call flush(6)
c        stop
      endif
c
      mask=mm(1:m,1:n,ip1)
      deallocate( mm )
c
      return
      end subroutine landfill

      subroutine psmooth(a,mask,amn,amx,m,n)
      implicit none
c
      integer m,n,mask(m,n)
      real    a(m,n),amn(m,n),amx(m,n)
c
c --- smooth under mask and within amn,amx range.
c
      integer, allocatable :: mm(:,:)
      real,    allocatable :: aa(:,:)
c
      integer i,ii,j,jj,ki,kj
      real    rss,sa
c
      real    s(-1:1,-1:1)
      save    s
      data    s / 1.0, 2.0, 1.0,
     &            2.0, 4.0, 2.0,
     &            1.0, 2.0, 1.0 /
c
      rss = 1.0/sum(s(:,:))
c
c     local copy of a.
c
      allocate( aa(m,n) )
      aa = a
c
c     adding a halo to mm simplifies ocean selection logic.
c
      allocate( mm(0:m+1,0:n+1) )
      mm(  0, : ) = 0
      mm(m+1, : ) = 0
      mm( : ,  0) = 0
      mm( : ,n+1) = 0
      mm(1:m,1:n) = mask
c
      do j= 1,n
        do i= 1,m
          if     (mm(i,j).eq.1) then
            sa = 0.0
            do kj= -1,1
              jj = j+kj
              do ki= -1,1
                ii = i+ki
                if     (mm(ii,jj).eq.1) then
                  sa = sa + s(ki,kj)*aa(ii,jj)  ! must use local copy of a
                else
                  sa = sa + s(ki,kj)*aa(i ,j )
                endif
              enddo
            enddo
            a(i,j) = max( amn(i,j),
     &                    min( amx(i,j), sa*rss ) )
          endif
        enddo
      enddo
c
      deallocate( aa, mm )
c
      return
      end subroutine psmooth

c
      subroutine handle_err(status)
      integer   :: status
         print *, "NetCDF return error ", status
!         stop "Stopped"
         call exit(1)
      end subroutine handle_err


