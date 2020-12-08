program prog_bis
 implicit none
 double precision :: x,dx,xmax,mu,w_ee_eff,derf_mu_x
 integer :: i,nx
 double precision :: mu0,mu1,mu2,mu3,mu4,mu5,mu6,mu_tmp_5,mu_tmp_10
 double precision :: mu7, mu8, mu9,mu10,mu11, mu12, mu13, mu14
 double precision :: mu_ten_no,ten_no_square,ten_no_lapl,ten_no_pot
 xmax = 15.d0
 nx = 100000
 dx = xmax/dble(nx)
 x = dx
 mu_ten_no = 0.86975d0
 mu2 = 0.3d0
 mu3 = 0.5d0
 mu4= 1.0d0
 mu5 = 1.5d0
 mu6 = 3.0d0
 mu7 = 5.0d0
 mu8 = mu_ten_no

 do i = 1, nx
  call eff_pot_ten_no(x,ten_no_square,ten_no_lapl,ten_no_pot)
  write(38,'(100(F16.10,X))')x,w_ee_eff(x,mu2),w_ee_eff(x,mu3),w_ee_eff(x,mu4),w_ee_eff(x,mu5)& 
                              ,w_ee_eff(x,mu6),w_ee_eff(x,mu7),w_ee_eff(x,mu8)&
                              ,ten_no_pot+1.d0/x,1.d0/x
   write(39,'(100(F16.10,X))')x,w_ee_eff(x,0.3d0),w_ee_eff(x,0.5d0), w_ee_eff(x,1.d0)
   write(40,'(100(F16.10,X))')x,w_ee_eff(x,0.5d0), derf_mu_x(x,0.5d0)&
                               ,w_ee_eff(x,1.0d0), derf_mu_x(x,1.0d0)
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

subroutine eff_pot_ten_no(r12,ten_no_square,ten_no_lapl,ten_no_pot)
 implicit none
 double precision, intent(in) :: r12 
 double precision, intent(out):: ten_no_square,ten_no_lapl,ten_no_pot
 integer :: i,j
 double precision :: coef_fit_ten_no_slat_gauss(6)
 double precision :: expo_fit_ten_no_slat_gauss(6)
 double precision :: coef,alpha,beta,contrib
 coef_fit_ten_no_slat_gauss(1) =-0.078215d0
 coef_fit_ten_no_slat_gauss(2) =-0.132037d0
 coef_fit_ten_no_slat_gauss(3) =-0.068633d0
 coef_fit_ten_no_slat_gauss(4) =-0.029047d0
 coef_fit_ten_no_slat_gauss(5) =-0.012063d0
 coef_fit_ten_no_slat_gauss(6) =-0.004346d0

 expo_fit_ten_no_slat_gauss(1) = 0.621698d0
 expo_fit_ten_no_slat_gauss(2) = 3.371717d0
 expo_fit_ten_no_slat_gauss(3) = 14.27116d0
 expo_fit_ten_no_slat_gauss(4) = 82.76522d0
 expo_fit_ten_no_slat_gauss(5) = 605.5295d0
 expo_fit_ten_no_slat_gauss(6) = 6596.808d0

 ten_no_pot = 0.d0

 ten_no_square = 0.d0
 do i = 1, 6
  alpha = expo_fit_ten_no_slat_gauss(i)
  do j = 1, 6
   beta = expo_fit_ten_no_slat_gauss(j)
   coef = coef_fit_ten_no_slat_gauss(i) * coef_fit_ten_no_slat_gauss(j)
   ten_no_square = ten_no_square +4.d0 * alpha * beta * coef * dexp(-(alpha+beta)*r12*r12) *r12*r12
  enddo
 enddo
  
 ten_no_lapl = 0.d0
 do i = 1, 6
  alpha = expo_fit_ten_no_slat_gauss(i)
  coef = coef_fit_ten_no_slat_gauss(i) 
  contrib =  dexp(-alpha*r12*r12)
  ten_no_lapl = ten_no_lapl + coef*(- 6.d0 * alpha * contrib + 4.d0 * alpha**2.d0 * r12*r12 * contrib)
 enddo
 ten_no_pot = -ten_no_lapl - ten_no_square
end

