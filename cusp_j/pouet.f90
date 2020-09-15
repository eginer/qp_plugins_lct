program pouet
 implicit none
 double precision :: dx,int_gauss,gauss_norm,derf_deriv
 double precision :: ymin,ymax,dy,y,num_erf,dy2,mu,derf_deriv_mu
 integer :: i,ny
 mu = 3.d0
 dx = 1.d-4
 ymin = 0.d0
 ymax = 5.d0
 ny = 10000
 dy = (ymax - ymin)/dble(ny)
 int_gauss = num_erf(100.d0,dx)
 y = ymin
 print*,'int_gauss = ',int_gauss
 print*,'sqpi      = ',dsqrt(dacos(-1.d0))/2.d0
 dy2 = 1.d-4
 do i = 1, ny
!  write(33,'(100(F16.10,X))')y,derf(y), num_erf(y,dx)/int_gauss,derf_deriv(mu*y,dy2),gauss_norm(mu*y)
  write(33,'(100(F16.10,X))')y,derf_deriv_mu(y,mu,dy),mu * gauss_norm(mu*y)
  y = y + dy
 enddo




end

double precision function num_erf(y,dx)
 implicit none
 double precision, intent(in) :: y,dx
 double precision :: x,accu,xmin
 integer :: nx,i
 xmin = 0.d0
 nx = int((y - xmin)/dx)
 x = xmin
 accu = 0.d0
 do i = 1, nx
  accu = accu + dexp(-x*x)
  x = x + dx
 enddo
 accu = accu * dx
end

double precision function derf_deriv(y,dy)
 implicit none
 double precision, intent(in) :: y,dy
 derf_deriv = ( derf(y + dy) - derf(y) ) / dy
end

double precision function derf_deriv_mu(y,mu,dy)
 implicit none
 double precision, intent(in) :: y,mu,dy
 derf_deriv_mu = ( derf(mu*(y + dy)) - derf(mu*y) ) / dy
end

double precision function gauss_norm(y)
 implicit none
 double precision, intent(in) :: y
 double precision :: norm
 norm = dsqrt(dacos(-1.d0))/2.d0
 norm = 1.d0/norm
 gauss_norm = dexp(-y*y) * norm
end
