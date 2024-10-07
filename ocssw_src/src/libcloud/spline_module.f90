module spline_module

  use GeneralAuxType
  implicit none

contains

SUBROUTINE spline(n, x,y,y2)
  IMPLICIT NONE
  REAL(single), DIMENSION(n), INTENT(IN) :: x,y
  REAL(single), DIMENSION(n), INTENT(OUT) :: y2
  INTEGER(integer_fourbyte), intent(in) :: n
  REAL(single)              :: yp1,ypn
  REAL(single), DIMENSION(size(x)) :: a,b,c,r
 
  yp1 = 9999.
  ypn = 9999.

  c(1:n-1)=x(2:n)-x(1:n-1)
  r(1:n-1)=6.0_single*((y(2:n)-y(1:n-1))/c(1:n-1))
  r(2:n-1)=r(2:n-1)-r(1:n-2)
  a(2:n-1)=c(1:n-2)
  b(2:n-1)=2.0_single*(c(2:n-1)+a(2:n-1))
  b(1)=1.0
  b(n)=1.0
  if (yp1 > 999.) then
    r(1)=0.0
    c(1)=0.0
  else
    r(1)=(3.0_single/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
    c(1)=0.5
  end if
  if (ypn > 999.) then
    r(n)=0.0
    a(n)=0.0
  else
    r(n)=(-3.0_single/(x(n)-x(n-1)))*((y(n)-y(n-1))/(x(n)-x(n-1))-ypn)
    a(n)=0.5
  end if
  call tridag(a(2:n),b(1:n),c(1:n-1),r(1:n),y2(1:n),n)


END SUBROUTINE spline


FUNCTION splint(n,x,xa,ya,y2a)
  use GeneralAuxType
  IMPLICIT NONE
  INTEGER(integer_fourbyte), intent(in) :: n
  REAL(single), DIMENSION(:), INTENT(IN) :: xa,ya,y2a
  REAL(single), INTENT(IN) :: x
  REAL(single) :: splint
  INTEGER(integer_fourbyte) :: khi,klo
  REAL(double) :: a,b,h, temp_1, temp_2, temp_3

  klo=max(min(locate(xa,x),n-1),1)
  khi=klo+1

  h=xa(khi)-xa(klo)

  a=(xa(khi)-x)/h
  b=(x-xa(klo))/h
 
  splint=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.0_single
END FUNCTION splint

FUNCTION locate(xx,x)
  use GeneralAuxType
  IMPLICIT NONE
  REAL(single), DIMENSION(:), INTENT(IN) :: xx
  REAL(single), INTENT(IN) :: x
  INTEGER(integer_fourbyte) :: locate
  INTEGER(integer_fourbyte) :: n,jl,jm,ju
  LOGICAL :: ascnd
  n=size(xx)
  ascnd = (xx(n) >= xx(1))
  jl=0
  ju=n+1
  do
    if (ju-jl <= 1) exit
      jm=(ju+jl)/2
      if (ascnd .eqv. (x >= xx(jm))) then
      jl=jm
    else
      ju=jm
    end if
  end do

  if (realsingle_s_equal(x,xx(1))) then
    locate=1
  else if (realsingle_s_equal(x,xx(n))) then
    locate=n-1
  else
    locate=jl
  end if
END FUNCTION locate

SUBROUTINE tridag(a,b,c,r,u,size)
  use GeneralAuxType, only: single, integer_fourbyte
  IMPLICIT NONE
  integer(integer_fourbyte), intent(in) :: size
  REAL(single),  INTENT(IN)  :: a(size-1),b(size),c(size-1),r(size)
  REAL(single),  INTENT(OUT) :: u(size)
  REAL(single)               :: gam(size)
  INTEGER(integer_fourbyte) :: n,j
  REAL(single) :: bet

  n=size
  bet=b(1)

  u(1)=r(1)/bet

  do j=2,n
    gam(j)=c(j-1)/bet
    bet=b(j)-a(j-1)*gam(j)
    u(j)=(r(j)-a(j-1)*u(j-1))/bet
  end do

  do j=n-1,1,-1
    u(j)=u(j)-gam(j+1)*u(j+1)
  end do
  END SUBROUTINE tridag

end module spline_module
