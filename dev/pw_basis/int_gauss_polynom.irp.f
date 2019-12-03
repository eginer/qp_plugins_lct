 subroutine num_int_gauss_pol(n,L,alpha,integral,nx)
  implicit none
!  BEGIN_DOC
!  computes numerically the following integral 
! 
!  \int(-L;+L) x^n exp(-alpha x^2)
!  END_DOC
  integer, intent(in) :: n
  double precision, intent(in) :: L, alpha
  double precision, intent(out):: integral
  integer, intent(in)         :: nx
  double precision :: pouet_tiny,x_min,x_max,domain,dx,x
  integer :: i
  pouet_tiny = dsqrt(40.d0/alpha)
  x_min = pouet_tiny
  x_min = - min(L,x_min)
  x_max = - x_min
  domain = x_max-x_min
  dx = domain/dble(nx)
  integral = 0.d0
  x = x_min
  do i = 1, nx
   x += dx
   integral += x**n * dexp(-alpha * x*x) 
  enddo
  integral *= dx
 
 end
 
 subroutine exact_int_gauss_pol(n,L,alpha,integral)
  implicit none
!  BEGIN_DOC
!  computes exactly the following integral 
! 
!  \int(-L;+L) x^n exp(-alpha x^2)
  END_DOC
  integer, intent(in) :: n
  double precision, intent(in) :: L, alpha
  double precision, intent(out):: integral
  include 'constants.include.F'
  double precision :: sqrt_alpha
  sqrt_alpha = dsqrt(alpha)
  integral = sqpi / sqrt_alpha * derf(sqrt_alpha * L)
  
 end
