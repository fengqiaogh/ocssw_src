module pfaast

	use FASCODE_routines
   use cld_tbl_names

	private
	
	public :: MODIS_fascode
	
	  integer, parameter :: Fnr = 17, Fnd = 10

      real ::	coefd(Fncd,Fnm,0:Fnd,Fnr),&
				coefo(Fnco,Fnm,0:Fnd,Fnr), &
                coefl(Fncl,Fnm,0:Fnd,Fnr),&
				coefs(Fncs,Fnm,0:Fnd,Fnr), &
                coefc(Fncc,Fnm,0:Fnd,Fnr)
	
contains


	subroutine MODIS_fascode(coeff_dir_path, year, jday, temp, wvmr, ozmr, theta, ang_2way, platform, &
							kban, jdet, taut, taut_2way, newang, newatm, new_2way, do_2way, iok, xxx, yyy ) 
	
! * MODIS band/detector 101-level fast transmittance routine
! .... version of 06.08.03

!    coeff_dir_path = directory that contains fast model coefficients
!        year = profile year
!        jday = profile day-of-year
!	 temp = temperature (Kelvin) profile
!	 wvmr = water-vapor mixing-ratio (g/kg) profile
!	 ozmr = ozone mixing-ratio (ppmv) profile
!	theta = local zenith angle
!	craft = TERRA, AQUA (either upper or lower case)
!	 kban = band number     (20...36)
!	 jdet = detector number (0...10) [ Product Order ]
!		 detector 0 is based on band-average response functions

!		 taut = total transmittance (see note below)
!		  iok = 0 if successful, 1 if I/O problem

! * NOTE: for kban = 26, return-arrays are filled with 1.0

! * PLOD/PFAAST regression model based on LBLRTM line-by-line transmittances.
! * Input temperatures, and water-vapor and ozone mixing ratios, must
! *	be defined at the pressure levels in array 'pstd'
! *    (see block data 'reference_atmosphere').
! * Units: temperature, deg-K; water vapor, g/kg; ozone, ppmv.
! * Logical units 31-35 are used for coefficient files.
! * Component taus are returned through common, product in 'taut'.
		character(*), intent(in) :: coeff_dir_path
		integer, intent(in) :: year, jday
      real, intent(inout) :: temp(*),wvmr(*),ozmr(*),taut(*), taut_2way(*)
		integer, intent(in) :: platform
		real, intent(in) :: theta, ang_2way
		integer, intent(in) :: kban, jdet, xxx, yyy
		integer, intent(inout) :: iok
		logical, intent(in) :: newang, newatm, new_2way, do_2way
		
      integer, parameter :: nd=10,nk=5,nl=101,nm=nl-1,koff=19,nr=17
	  integer, parameter :: nxc= 4,ncc=nxc+1,lencc=ncc*nm,lenccb=lencc*4
     integer, parameter :: nxd= 8,ncd=nxd+1,lencd=ncd*nm,lencdb=lencd*4
      integer, parameter :: nxo= 9,nco=nxo+1,lenco=nco*nm,lencob=lenco*4
      integer, parameter :: nxl= 2,ncl=nxl+1,lencl=ncl*nm,lenclb=lencl*4
      integer, parameter :: nxs=11,ncs=nxs+1,lencs=ncs*nm,lencsb=lencs*4
      integer, parameter :: ndt=nd+1,nrps=nr*ndt,nxw=nxl+nxs
      real, parameter :: slp=1.5/365.0
      real, parameter :: smag=3.0
      real, parameter :: pi=3.14159
	  real, parameter :: soff=0.41
      real, parameter :: coff=337.5

!      common/stdatm/pstd(nl),tstd(nl),wstd(nl),ostd(nl)
!     common/taudwo/taud(nl),tauw(nl),tauo(nl)
	  real :: taud(nl), tauw(nl), tauo(nl)
	  real :: taud_2way(nl), tauw_2way(nl), tauo_2way(nl)
	  
       real :: bufs(lencs)
!      real :: pavg(nm),tref(nm),wref(nm),oref(nm)
!      real ::  tavg(nm),wamt(nm),oamt(nm),secz(nm)
      real :: tauc(nl),tauc_2way(nl), tlas(nl),wlas(nl),olas(nl)

      real*4 x,rco2,ratio,tau_test
      character*28 xfile
      character*256 path
      character*256 cfile(nk),dfile
      character*6 craft
      character*6 cbt,cba
      character*6 cst,csa
      character*3 comp(nk)
      character*3 cbe,cle
      integer*4 lengcf(nk)
      integer*4 lengcx(nk)
      integer*4 iuc(nk)
      integer*2 jj
      logical big_endian

	  real ::  zlas, zlas2way
	  
     integer :: stat
     integer :: trans_id(5) = (/FAST_TRANS_COEFF_DRY, FAST_TRANS_COEFF_OZO, &
       FAST_TRANS_COEFF_WTS, FAST_TRANS_COEFF_WTL, FAST_TRANS_COEFF_WCO /)
	  integer :: ksat, iux, m, j, kk, ikrec, krec, krecx, k, lencx, l, i
	  integer :: lencf
	  integer*2 :: nsat
	  real :: dt, dw, fdo, datm, zsec, zsec_2way

     integer :: firsttime = 0 ! WDR To limit reads to first fill
     integer :: do_init ! WDR reduce computations in calpir
 
     ! WDR for test-out of new netcdf  file read
     integer*4 :: start(5) = (/ 1, 1, 1, 1, 1/), edge(5), nc_stat
     integer :: fid, sds_id, tbl_no
     real :: coef_dry( Fncd, Fnm, 0:Fnd, Fnr )
     real :: coef_oz( Fnco, Fnm, 0:Fnd, Fnr )
     real :: coef_wtr1( Fncs, Fnm, 0:Fnd, Fnr )
     real :: coef_wtr2( Fncl, Fnm, 0:Fnd, Fnr )
     real :: coef_cont( Fncc, Fnm, 0:Fnd, Fnr )
     character(100) :: tbl_nam
     integer :: ids(5), lendim

!      data cinit/'zzzzzz'/
	
     taud_2way = 0  ! WDR-UIV
     ratio = 1  ! WDR-UIV

		tlas = nl*0.
		wlas = nl*0.
		olas = nl*0.
		zlas = -999.
		zlas2way = -999.
		xfile = '/modisdet.com.101.xxx_end.v3'
		cbe = 'big'
		cle = 'lit'
		
		cbt = 'TERRA'
		cba = 'AQUA'
		cst = 'terra'
		csa = 'aqua'

		comp = (/'dry','ozo','wts','wtl','wco'/)
		lengcf = (/lencdb,lencob,lencsb,lenclb,lenccb/)
		lengcx = (/lencd,lenco,lencs,lencl,lencc/)

		iok = 0
		
		if (platform == 0) then 
!       Shifted bands 34-36, LBL 11.3
			craft = 'TERRA'
			path = coeff_dir_path
		else if (platform == 1) then 
			craft = 'AQUA'
			path = coeff_dir_path
		else
			write(*,'(''tran_modisd101- unknown spacecraft '',i2)') platform
			iok = 1
			return 
		endif
				
		if (craft /= cinit) then 
		
			if (craft == "TERRA") then	
				ksat = 1
			else if (craft == "AQUA") then 
				ksat = 2
			else
				write(*,'(''tran_modisd101- unknown spacecraft '',a6)') craft
				iok = 1
				return
			endif
			
! select big- or little-endian libraries.
!			if (big_endian()) then 
				xfile(19:21) = cbe
!			else
!				xfile(19:21) = cle
!			endif
			
! WDR if first time through, do the reads, otherwise rely on pre-read data
        if( firsttime == 0 ) then
          firsttime = 1
         !
         ! WDR replace with get_cld_tbl()
         ! we only have 1 set of tables and they are for modis aqua

         !  The following code is for the old binary yfast trans coeff files
         !  keep it behind a do_old_files variable so a return is easy
         do_old_files = 0
         coefd = 0
         coefo = 0
         coefs = 0
         coefl = 0
         coefc = 0
         if( do_old_files == 1 ) then
           iux = 30
           do m=1, nk
             iux = iux + 1
             call get_cld_tbl( 0, trans_id(m), dfile, stat )
             if( stat /= 0 ) then
               print *, __FILE__,__LINE__, &
                 'Unable to get the FAST_TRANS_COEFF file', m
               stop 123
             endif
             lencf=lengcf(m)
             open( iux, file=dfile, recl=lencf, access='direct', status='old', &
               convert='big_endian')
             iuc(m) = iux
             cfile(m)=dfile
           enddo
    
         ! first read each files fill-record for band 26/det 0
         ! and verify satellite number stored in word 1
         ! note: number of levels is in word 2, creation date 
         ! (yyyyddd) is in word 3
    	 	ikrec=nrps*(ksat-1)
    	   krecx=ikrec+7
    	   do k=1,nk
    	     lencx=lengcx(k)
    		  read(iuc(k),rec=krecx) (bufs(j),j=1,lencx)
    	     nsat=bufs(1)
    	     if(nsat /= ksat) then
    	       dfile=cfile(k)
    		    write(*,'(''In tran_modisd101 ... requested data for '', &
    & ''satellite '',i1/'' but read data for '', ''satellite '',i1,'' from file '',a80)') ksat,nsat,dfile
    	      iok = 1
    		   return
    		   endif
          enddo
    
        ! *     now read in the coefficients
        krec=ikrec
        do l=0,nd
          do k=1,nr
            krec=krec+1
            read(iuc(1),rec=krec) ((coefd(i,j,l,k),i=1,ncd),j=1,nm)
            read(iuc(2),rec=krec) ((coefo(i,j,l,k),i=1,nco),j=1,nm)
            read(iuc(3),rec=krec) ((coefs(i,j,l,k),i=1,ncs),j=1,nm)
            read(iuc(4),rec=krec) ((coefl(i,j,l,k),i=1,ncl),j=1,nm)
            read(iuc(5),rec=krec) ((coefc(i,j,l,k),i=1,ncc),j=1,nm)
          enddo
        enddo
      do k=1,nk
        close(iuc(k))
      enddo
    endif
    ! WDR I believe this is the end of first time read
    ! WDR now, read the new trans coeff arrays from the new file and compare
         call get_cld_tbl( 0, trans_id(1), dfile, stat )
         nc_stat = nf_open( dfile, NF_NOWRITE, fid )
         call cld_fchk( nc_stat, __FILE__, __LINE__ )
         start(5) = ksat
         edge(5) = 1
         edge(4) = Fnr
         edge(3) = Fnd + 1
         edge(2) = Fnm

         do tbl_no = 1,5
           select case (tbl_no)
             case (1)
               edge(1) = Fncd
               tbl_nam = "Dry_air_component"
             case (2)
               edge(1) = Fnco
               tbl_nam = "Ozone_component"
             case (3)
               edge(1) = Fncs
               tbl_nam = "Water_component_1"
             case (4)
               edge(1) = Fncl
               tbl_nam = "Water_component_2"
             case (5)
               edge(1) = Fncc
               tbl_nam = "Continuum_component"
           end select
           nc_stat = nf_inq_varid( fid, tbl_nam, sds_id )
           call cld_fchk( nc_stat, __FILE__, __LINE__ )
           ! find the dimids and their size
           nc_stat = nf_inq_vardimid(fid, sds_id, ids )
           !do l=1,5
           !  nc_stat = nf_inq_dimlen(fid, ids(l), lendim)
           !  print*,'dim len: ', lendim
           !enddo
           select case (tbl_no)
             case (1)
               nc_stat = nf_get_vara_real( fid, sds_id, start, edge, coef_dry )
               call cld_fchk( nc_stat, __FILE__, __LINE__ )
               !print*, 'in pfasst, dif min, max for: ', tbl_nam, &
               !  minval(coef_dry-coefd), maxval(coef_dry-coefd)
             case (2)
               nc_stat = nf_get_vara_real( fid, sds_id, start, edge, coef_oz )
               call cld_fchk( nc_stat, __FILE__, __LINE__ )
               !print*, 'in pfasst, dif min, max for: ', tbl_nam, &
               !  minval(coef_oz-coefo), maxval(coef_oz-coefo)
             case (3)
               nc_stat = nf_get_vara_real( fid, sds_id, start, edge, coef_wtr1 )
               call cld_fchk( nc_stat, __FILE__, __LINE__ )
              ! print*, 'in pfasst, dif min, max for: ', tbl_nam, &
              !   minval(coef_wtr1-coefs), maxval(coef_wtr1-coefs)
             case (4)
               nc_stat = nf_get_vara_real( fid, sds_id, start, edge, coef_wtr2 )
               call cld_fchk( nc_stat, __FILE__, __LINE__ )
               !print*, 'in pfasst, dif min, max for: ', tbl_nam, &
               !  minval(coef_wtr2-coefl), maxval(coef_wtr2-coefl)
             case (5)
               nc_stat = nf_get_vara_real( fid, sds_id, start, edge, coef_cont )
               call cld_fchk( nc_stat, __FILE__, __LINE__ )
               !print*, 'in pfasst, dif min, max for: ', tbl_nam, &
               !  minval(coef_cont-coefc), maxval(coef_cont-coefc)
           end select
         enddo
         !
         !  WDR - this will use the new nc file coefficient arrays
         coefd = coef_dry
         coefo = coef_oz
         coefs = coef_wtr1
         coefl = coef_wtr2
         coefc = coef_cont
         !
         !  close the file and be off
         nc_stat = nf_close( fid )
         call cld_fchk( nc_stat, __FILE__, __LINE__ )

        end if
		
			call conpir(pstd,tstd,wstd,ostd,nl,1,pavg,tref,wref,oref)

			cinit=craft
			iok=0
	  
!     End "first-time-through-only" block.	   
      endif
		

	  if(newatm) then
	   call conpir(pstd,temp,wvmr,ozmr,nl,1,pavg,tavg,wamt,oamt)
	  endif

! we have to enable doing two-way angles on an as-needed basis

	if(newang) then
		zsec = secant(theta)
		secz = zsec
	endif

	if (do_2way .and. new_2way) then 
		zsec_2way = secant(ang_2way)
		secz_2way = zsec_2way
	endif
! WDR try to avoid all init computations for same profile
   do_init = 1
	if(newang .or. newatm ) then
	   call calpir(tref,wref,oref,tavg,wamt,oamt,pavg,secz, &
     			 nm,nxd,nxw,nxo,nxc,xdry,xwet,xozo,xcon,do_init)
	endif

	if (do_2way .and. (new_2way .or. newatm)) then
		   call calpir(tref,wref,oref,tavg,wamt,oamt,pavg,secz_2way, &
     			 nm,nxd,nxw,nxo,nxc,xdry_2way,xwet_2way,xozo_2way,xcon_2way,do_init)
	endif


	if(kban == 26) then
	   do j=1,nl
	      taud(j)=1.0
	      tauo(j)=1.0
	      tauw(j)=1.0
	      taut(j)=1.0
	   enddo
	   return
	 endif

	  j=jdet
	  k=kban-koff
! *   dry
	  call taudoc(ncd,nxd,nm,coefd(:,:,j,k),xdry,taud)

!		do j=1, nl
!			print*, j, taud(j)
!		end do

!         Adjust dry tau for changes in CO2 concentration from model (R. Frey)

#ifdef CT_CODE
          x = (year - 1980) * 365.25 + jday
          rco2 = (slp*x - smag*sin(2*pi*(x/365.25 + soff))) + coff
          ratio=rco2/380.0
          do jj=1,nl
            tau_test = taud(jj)
            if(taud(jj) > 0.0 .and. taud(jj) < 1.0) then
              taud(jj)=taud(jj)**ratio
            endif
          enddo

! *   ozo
	  call taudoc(nco,nxo,nm,coefo(:,:,j,k),xozo,tauo)
#else
! we set the ozone to 0. if we're not doing CT, so transmittance due to ozone is 1.0
	   tauo = 1.0
#endif

! *   wet
	  call tauwtr(ncs,ncl,nxs,nxl,nxw,nm,coefs(:,:,j,k),coefl(:,:,j,k),xwet,tauw)
	  
	  call taudoc(ncc,nxc,nm,coefc(:,:,j,k),xcon,tauc)

	  do jj=1,nl
	   tauw(jj)=tauw(jj)*tauc(jj)
	  enddo
! * total
	  do jj=1,nl
	   taut(jj)=taud(jj)*tauo(jj)*tauw(jj)
	  enddo
	
	  if (do_2way) then 
	  
		j=jdet
	  
		call taudoc(ncd,nxd,nm,coefd(:,:,j,k),xdry_2way,taud_2way)
		do jj=1,nl
			tau_test = taud_2way(jj)
            if(taud_2way(jj) > 0.0 .and. taud_2way(jj) < 1.0) then
              taud_2way(jj)=taud_2way(jj)**ratio

            endif
		enddo
#ifdef CT_CODE 
		call taudoc(nco,nxo,nm,coefo(:,:,j,k),xozo_2way,tauo_2way)
#else
		tauo_2way = 1.0
#endif
		call tauwtr(ncs,ncl,nxs,nxl,nxw,nm,coefs(:,:,j,k),coefl(:,:,j,k),xwet_2way,tauw_2way)		
		call taudoc(ncc,nxc,nm,coefc(:,:,j,k),xcon_2way,tauc_2way)


		do jj=1,nl
			tauw_2way(jj)=tauw_2way(jj)*tauc_2way(jj)
		enddo
! * total
		do jj=1,nl
			taut_2way(jj)=taud_2way(jj)*tauo_2way(jj)*tauw_2way(jj)
		enddo
	  	  
	  
	  endif
	
	
	
	end subroutine MODIS_fascode	





end module pfaast
