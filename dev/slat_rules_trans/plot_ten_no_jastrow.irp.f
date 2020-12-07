program print_ten_no
 implicit none
 integer :: i,nx
 double precision :: x,dx,xmax
 double precision :: slater_fit_ten_no,full_jastrow_mu,c0,mu,f_mu
 xmax = 5.d0
 nx = 10000
 dx =xmax/dble(nx)
 c0 = slater_fit_ten_no(0.d0)
 print*,'slater_fit_ten_no(0) = ',c0
 mu = -1.d0/(2.d0 * dsqrt(dacos(-1.d0))*c0 )
 print*,'mu = ',mu
 x = 0.d0
 do i = 1, nx
  write(33,'(100(F16.10,x))')x,slater_fit_ten_no(x),f_mu(mu,x) , dexp(slater_fit_ten_no(x)), full_jastrow_mu(mu,x)
  x += dx
 enddo

end
