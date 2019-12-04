 subroutine num_int_gauss_cos(n,alpha,integral,nx)
  implicit none
!  BEGIN_DOC
!  computes numerically the following integral 
! 
!  \int(-n*pi;+n*pi)  exp(-alpha x^2) * cos(x)
!  END_DOC
  integer, intent(in) :: n
  double precision, intent(in) :: alpha
  double precision, intent(out):: integral
  integer, intent(in)         :: nx
  double precision :: pouet_tiny,x_min,x_max,domain,dx,x
  double precision :: L
  integer :: i
  include 'constants.include.F'
  L = dble(n) * pi
  pouet_tiny = dsqrt(40.d0/alpha)
  x_min = pouet_tiny
  x_min = - min(L,x_min)
  x_max = - x_min
  domain = x_max-x_min
  dx = domain/dble(nx)
  integral = 0.d0
  x = x_min
  do i = 1, nx+1
   x += dx
   integral += dcos(x) * dexp(-alpha * x*x) 
  enddo
  integral *= dx
 
 end
 

 subroutine int_gauss_cos(m,n,alpha,integral,error)
  implicit none
!  BEGIN_DOC
!  computes the following integral 
! 
!  \int(-m*pi;+m*pi)  exp(-alpha x^2) * cos(x)
!
!  using the Taylor expansion of the gaussian up to order 2*n
!  END_DOC
  integer, intent(in) :: m,n
  double precision, intent(in) :: alpha
  double precision, intent(out):: integral
  double precision :: error,int_xn_cos_m

 integer :: i
 double precision :: pol(0:n),tmp
 call taylor_gauss_pol(alpha,n,pol)
 integral = 0.d0
 do i = 1, n
  tmp = int_xn_cos_m(m,2*i)
  integral += pol(i) * tmp 
  error = pol(i) * tmp 
 enddo

end
