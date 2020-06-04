program ao_two_e_eff_ints
  implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC
  print *, 'Hello world'
  integer :: i,nx
  double precision :: x,dx,xmax,fit_1_erf_x,fit_1_erf_x_2
  double precision :: eff_pot_gauss,eff_pot_fit_gauss
  xmax = 10.d0
  nx = 1000
  dx = xmax/dble(nx)
  x = dx
  do i = 1, nx
   write(33,'(100(F16.10,X))')x,(1.d0 - derf(mu_erf*x)),fit_1_erf_x(x),(1.d0 - derf(mu_erf*x))**2.d0,fit_1_erf_x_2(x),eff_pot_gauss(x,mu_erf),eff_pot_fit_gauss(x)
   x += dx
  enddo

end
