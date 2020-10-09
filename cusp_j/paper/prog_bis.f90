program prog_bis
 implicit none
 double precision :: x,dx,xmax,mu,w_ee_eff,derf_mu_x
 integer :: i,nx
 double precision :: mu0,mu1,mu2,mu3,mu4,mu5,mu6,mu_tmp_5,mu_tmp_10
 double precision :: mu7, mu8, mu9,mu10,mu11, mu12, mu13, mu14
 xmax = 15.d0
 nx = 1000
 dx = xmax/dble(nx)
 x = dx
 mu0 = 0.05d0
 mu1 = 0.10d0
 mu2 = 0.20d0
 mu3 = 0.3d0
 mu4 = 0.4d0
 mu5 = 0.5d0
 mu6 = 0.6d0
 mu7 = 0.7071d0
 mu8 = 0.8d0
 mu9 = 0.9d0
 mu10= 1.0d0
 mu11= 1.5d0
 mu12= 3.0d0
 mu13= 5.0d0
 mu_tmp_5 = 0.5284432686368106d0
 mu_tmp_10 = 1.27844d0
 do i = 1, nx
  write(38,'(100(F16.10,X))')x, w_ee_eff(x,mu0),w_ee_eff(x,mu1),w_ee_eff(x,mu2)& 
                              ,w_ee_eff(x,mu3),w_ee_eff(x,mu4),w_ee_eff(x,mu5),w_ee_eff(x,mu6)&
                              ,w_ee_eff(x,mu7),w_ee_eff(x,mu8),w_ee_eff(x,mu9),w_ee_eff(x,mu10)& 
                              ,w_ee_eff(x,mu11),w_ee_eff(x,mu12),w_ee_eff(x,mu13),w_ee_eff(x,mu14)  
  write(39,'(100(F16.10,X))')x, derf_mu_x(x,mu0),derf_mu_x(x,mu1),derf_mu_x(x,mu2)& 
                              ,derf_mu_x(x,mu3),derf_mu_x(x,mu4),derf_mu_x(x,mu5),derf_mu_x(x,mu6)&
                              ,derf_mu_x(x,mu7),derf_mu_x(x,mu8),derf_mu_x(x,mu9),derf_mu_x(x,mu10)& 
                              ,derf_mu_x(x,mu11),derf_mu_x(x,mu12),derf_mu_x(x,mu13),derf_mu_x(x,mu14) 
  write(40,'(100(F16.10,X))')x, derf_mu_x(x,mu_tmp_5),w_ee_eff(x,mu5),derf_mu_x(x,mu_tmp_10),w_ee_eff(x,mu10)
  x = x + dx
 enddo
end

double precision function w_ee_eff(x,mu)
 implicit none
 double precision, intent(in) :: x,mu 
 double precision :: derf_mu_x
 double precision :: inv_sq_pi 
 inv_sq_pi = 1.d0/dsqrt(dacos(-1.d0))
 w_ee_eff = derf_mu_x(x,mu)  + mu * inv_sq_pi * dexp(-(mu*x)**2.d0) - 0.25d0 * (1.d0 - derf(mu*x))**2.d0
end


double precision function derf_mu_x(x,mu)  
 implicit none
 double precision, intent(in) :: mu,x
 double precision :: inv_sq_pi 
 inv_sq_pi = 1.d0/dsqrt(dacos(-1.d0))
  if(dabs(x).gt.1.d-6)then
   derf_mu_x = derf(mu * x)/x
  else
   derf_mu_x =  inv_sq_pi * 2.d0 * mu * (1.d0 - mu*mu*x*x/3.d0)
  endif
end

double precision function mu_rs_from_mu_h(mu_h)
 implicit none
 double precision, intent(in) :: mu_h
 double precision :: sq_pi 
 sq_pi = dsqrt(dacos(-1.d0))
 mu_rs_from_mu_h = 1.5d0 * mu_h - sq_pi/8.d0
end

double precision function mu_h_from_mu_rs(mu_rs)
 implicit none
 double precision, intent(in) :: mu_rs
 double precision :: sq_pi 
 sq_pi = dsqrt(dacos(-1.d0))
 mu_h_from_mu_rs = 2.d0/3.d0 * mu_rs + sq_pi/12.d0
end
