 subroutine remove_edge_scenes_ext(&
                              CSR_results, &
                              xsize, ysize, &
                              status)

   use nonscience_parameters
   use modis_cloudstructure
   use mod06_run_settings
   use GeneralAuxType
   use core_arrays
   
   implicit none


!  for each scene, if the pixel is not surrounded by cloud we eliminate it
!  as a retreival candidate
!  In effect we remove all edge pixels

   integer, intent(in)  :: xsize, ysize
   integer(integer_onebyte), dimension(xsize,ysize), intent(inout) :: CSR_results
   integer, intent(out) :: status

   integer :: i, j, tempmask(3,3), lo_i, hi_i, lo_j, hi_j, sum_tempmask, main_tempmask, ii, jj
   integer, dimension(:,:), allocatable :: mask_solution


   status = 0

#if MAS_INST | EMAS_INST 
! this code is modified to be using the CSR results from the hi-res CSR product MAS_PR06CSR

	do i=1, xsize
		do j=1, ysize
			if (restoral_pos(i,j) == 2) then 
				CSR_results(i,j) = 1
			endif
		end do
	end do
#endif

  
 end subroutine remove_edge_scenes_ext



subroutine remove_edge_scenes(cloudmask, &
                              CSR_results, &
                              xsize, ysize, &
                              status)
!  WDR 
!   cloudmask   I   structure(?) with cloud state for each pixel in this chunk
!   CSR_results O   indicator (= 1) that this pixel has a cloud-free 
!                      pixel adjacent to it
!   xsize, ysize  I  Size of the chunk of data in pix, lin
!   status     O    result = 0 always

   use nonscience_parameters
   use modis_cloudstructure
   use mod06_run_settings
   use GeneralAuxType
   
   implicit none


!  for each scene, if the pixel is not surrounded by cloud we eliminate it
!  as a retreival candidate
!  In effect we remove all edge pixels

   integer, intent(in)  :: xsize, ysize
   type(processflag), dimension(xsize,ysize),  intent(in)         :: cloudmask
   integer(integer_onebyte), dimension(xsize,ysize), intent(inout) :: CSR_results
   integer, intent(out) :: status

   integer :: i, j, tempmask(3,3), lo_i, hi_i, lo_j, hi_j, sum_tempmask, main_tempmask, ii, jj
   integer, dimension(:,:), allocatable :: mask_solution


   status = 0

! this code is reworked to incorporate the tests for edge lines. Basically it's quite
! simple. If there is at least one point around the current point that is not cloudy, 
! we remove that point. 

   i=1
   j=1

   allocate(mask_solution(xsize, ysize))
   mask_solution = 0

   do i=1, xsize
      if (i==1) then 
         lo_i = i
         hi_i = i+1
      else if (i==xsize) then 
         lo_i = i-1
         hi_i = i
      else
         lo_i = i-1
         hi_i = i+1
      endif

      if (i==1 .or. i==xsize) then 
         sum_tempmask = 6
      else
         sum_tempmask = 9
      endif

      main_tempmask = sum_tempmask

      do j=1, ysize

         if (j==1) then 
            lo_j = j
            hi_j = j+1
         else if (j==ysize) then 
            lo_j = j-1
            hi_j = j
         else
            lo_j = j-1
            hi_j = j+1
         endif

         if (j==1 .or. j==ysize) then 
            if (sum_tempmask == 6) then 
               sum_tempmask = 4
            else
               sum_tempmask = 6
            endif
         endif


            tempmask = 0

         ! attempt edge removal only if a pixel is cloudy to begin with, otherwise, why would we 
         ! ever want to run the restoral in the first place. G.Wind 4.20.05
         if (cloudmask(i,j)%cloudobserved) then 

            do ii=lo_i, hi_i
               do jj=lo_j, hi_j
                  if (cloudmask(ii,jj)%cloudobserved .and. CSR_results(ii,jj) == 0) then 
                     tempmask( ii-lo_i+1, jj-lo_j+1) = 1
                  end if
               end do
             end do

             if (sum(tempmask) <  sum_tempmask) &
                  mask_solution(i,j) = 1
             
           endif

          sum_tempmask = main_tempmask

         enddo
      enddo
 
      do i=1, xsize
         do j=1, ysize
! don't overwrite the pixels with other CSR results
            if (mask_solution(i,j) == 1 .and. CSR_results(i,j) == 0) then 
               CSR_results(i,j) = 1
            endif

         end do
      end do


      deallocate(mask_solution)
  
end subroutine remove_edge_scenes
