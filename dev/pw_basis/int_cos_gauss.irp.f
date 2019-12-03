subroutine int_gauss_cos(n,alpha,nx,dx,integral)
 implicit none
 BEGIN_DOC
 ! int(0,+n*pi) cos(x) * exp(-alpha *x^2)
 END_DOC
 double precision, intent(in)  :: alpha,dx
 integer, intent(in) :: n
 double precision, intent(out) :: integral
 integer, intent(out) :: nx

 double precision :: domain,accuracy,xmin,x
 integer :: i
 accuracy = 1.d-18
 include 'constants.include.F'
 domain = dble(n) * pi
 xmin = dsqrt(dabs(dlog(accuracy))/alpha)
 domain = min(xmin,domain)
 nx = int(domain/dx)
 integral = 0.d0
 x = 0.d0
 do i = 1, nx
!  integral += dcos(x) * dexp(-alpha*x*x)
  integral +=  dexp(-alpha*x*x)
  x += dx
 enddo
 integral *= dx
end

