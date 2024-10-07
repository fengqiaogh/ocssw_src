module GeneralAuxType

!f90
!
!description:
!           data type declarations and parameters for MODIS retrieval code
!
!
!input parameters: none
!
!
!output parameters: none
!revision history:
!team-unique header:
!             mark gray, lac
!             climate and radiation branch
!             nasa goddard space flight center
!             greenbelt, maryland, u.s.a.


!... symbolic names for data types:

integer, parameter :: integer_fourbyte= selected_int_kind(9)
integer, parameter :: integer_twobyte = selected_int_kind(4)
integer, parameter :: integer_onebyte = selected_int_kind(2)
integer, parameter :: single = kind(1.0)
integer, parameter :: double = kind(1.0d0)
integer, parameter :: singlecomplex = kind((1.0,1.0))
integer, parameter :: doublecomplex = kind((1.0d0,1.0d0))
integer, parameter :: logical = kind(.true.)

!... debug variables
integer :: iterationX, iterationY 
integer :: pixX, pixY
integer :: number_of_iterationsX, number_of_iterationsY
logical :: debugPRN

! metadata variables
integer :: total_number_of_pixels
integer :: IM_cloudy_count
integer :: IM_ice_cloud_count
integer :: IM_water_cloud_count
integer :: IM_successful_retrieval_count
integer :: IM_undet_count

integer, parameter :: OUTPUT_UNIT = 45

!... useful mathematical constants

real(single), parameter :: pi=3.141592653589793238462643383279502884197_single
real(single), parameter :: pio2=1.57079632679489661923132169163975144209858_single
real(single), parameter :: twopi=6.283185307179586476925286766559005768394_single

real(double), parameter :: pi_d=3.141592653589793238462643383279502884197_double
real(double), parameter :: pio2_d=1.57079632679489661923132169163975144209858_double
real(double), parameter :: twopi_d=6.283185307179586476925286766559005768394_double
real(double), parameter :: d2r_d = 0.017453292519943295_double

real(single),       parameter :: missingvalue = -1.0

contains

logical function real_s_equal(x,y)
   real :: x, y
   real_s_equal = (abs(x-y) <= epsilon(x)) 
end function real_s_equal

logical function realsingle_s_equal(x,y)
   real(single) :: x, y
   realsingle_s_equal = (abs(x-y) <= epsilon(x)) 
end function realsingle_s_equal

subroutine realsingle_s_where_equal(x,y)
   real(single) ,intent(inout) :: x(:)
   real(single) :: y

   where( abs(x - y) <= epsilon(y) )
      x = 1.
   elsewhere
      x = 0.
   endwhere

end subroutine realsingle_s_where_equal

end module GeneralAuxType
