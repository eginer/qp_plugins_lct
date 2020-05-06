double precision function f_mu(mu,x)
 include 'constants.include.F'
 BEGIN_DOC
! correlation function corresponding to the erf(mu*x)/x interaction
 END_DOC
 implicit none
 double precision, intent(in) :: mu,x
 f_mu = 0.5d0 * x * (1.d0 - derf(mu*x)) - inv_sq_pi*0.5d0/mu * dexp(-x*mu*x*mu)
end

double precision function full_jastrow_mu(mu,x)
 implicit none
 double precision, intent(in) :: x,mu
 double precision :: f_mu
 full_jastrow_mu = dexp(f_mu(mu,x))
end


double precision function tilde_f_mu(mu,x)
 include 'constants.include.F'
 BEGIN_DOC
! rescaled correlation function corresponding to the erf(mu*x)/x interaction
 END_DOC
 implicit none
 double precision, intent(in) :: mu,x
 tilde_f_mu = sqpi * mu * x * (1.d0 - derf(mu*x)) - dexp(-mu*x*mu*x) 
end

double precision function exp_tilde_f_mu(mu,x)
 implicit none
 double precision, intent(in) :: x,mu
 double precision :: tilde_f_mu
 exp_tilde_f_mu = dexp(tilde_f_mu(x,mu))
end

double precision function exp_tilde_f_fit(mu,x)
 implicit none
 double precision, intent(in) :: x,mu
 double precision :: f_tilde_fit
 exp_tilde_f_fit = dexp(f_tilde_fit(x,mu))
end

double precision function exp_f_fit(mu,x)
 implicit none
 include 'constants.include.F'
 double precision, intent(in) :: x,mu
 double precision :: f_tilde_fit,a0
 a0 = 1.d0/(2.d0*dsqrt(pi)*mu)
 exp_f_fit = dexp(f_tilde_fit(x,mu)*a0)
end


subroutine give_param_fit_f_h(mu,gama,delta)
 include 'constants.include.F'
 BEGIN_DOC
! Returns the gama and delta parameters such that the function 
!
! f_tilde is fitted perfectly up to second order by the function h
 END_DOC
 implicit none
 double precision, intent(in) :: mu
 double precision, intent(out):: gama,delta
 gama = sqpi * mu
 delta = -mu*mu + 0.5d0 * gama * gama 
end

double precision function f_tilde_fit(mu,x)
 implicit none
 double precision, intent(in) :: mu,x
 double precision :: slater_fit_gam,gama,delta
 call give_param_fit_f_h(mu,gama,delta)
 f_tilde_fit = -slater_fit_gam(x,gama) * dexp(-delta *x*x)
! f_tilde_fit =  -dexp(-gama * x -delta *x*x)
end

double precision function ovlp_exp_f_phi_ij(mu,r1,A_center,B_center,power_A,power_B,alpha,beta,n_taylor)
 implicit none
 include 'constants.include.F'
 double precision, intent(in)    :: r1(3), mu        ! where the Jastrow factor is centered and its shape 
 integer, intent(in)             :: n_taylor         ! order of the Taylor expansion of the exponential
 double precision, intent(in)    :: A_center(3),B_center(3),alpha,beta ! gaussian/polynoms "A" and "B"
 integer, intent(in)             :: power_A(3),power_B(3)
 BEGIN_DOC
! integral of J(r1,r2)^2 phi_i(r2) phi_j(r2) erf(mu r_12)/r12
 END_DOC
 
 double precision :: gama, delta,shank_general
 double precision :: array_ints(0:n_taylor),integral_taylor,sum(0:n_taylor)
 double precision :: a0,ai
 integer :: i
 a0 = 1.d0/(dsqrt(pi)*mu)
 call give_param_fit_f_h(mu,gama,delta)
 call exp_dl_ovlp_stg_phi_ij(a0,r1,gama,delta,A_center,B_center,power_A,power_B,alpha,beta,n_taylor,array_ints,integral_taylor)
 sum(0) = array_ints(0) 
 do i = 1, n_taylor
  ai = a0**dble(i)
  sum(i) = sum(i-1) + array_ints(i) * (-1.d0)**dble(i) * fact_inv(i) * ai
 enddo
 ovlp_exp_f_phi_ij = shank_general(sum,n_taylor,n_taylor) 

end


double precision function erf_exp_f_phi_ij(muj,muc,r1,A_center,B_center,power_A,power_B,alpha,beta,n_taylor)
 implicit none
 BEGIN_DOC
! integral of J_muj(r1,r2)^2 phi_i(r2) phi_j(r2) erf(muc r_12)/r12
 END_DOC
 include 'constants.include.F'
 double precision, intent(in)    :: r1(3), muj,muc   ! where the Jastrow factor is centered and its shape 
 integer, intent(in)             :: n_taylor         ! order of the Taylor expansion of the exponential
 double precision, intent(in)    :: A_center(3),B_center(3),alpha,beta ! gaussian/polynoms "A" and "B"
 integer, intent(in)             :: power_A(3),power_B(3)
 
 double precision :: gama, delta,shank_general
 double precision :: array_ints(0:n_taylor),integral_taylor,sum(0:n_taylor)
 double precision :: a0,ai
 integer :: i

 a0 = 1.d0/(dsqrt(pi)*muj)
 call give_param_fit_f_h(muj,gama,delta)
 call exp_dl_erf_stg_phi_ij(a0, r1,gama,delta,A_center,B_center,power_A,power_B,alpha,beta,r1,muc,n_taylor,array_ints,integral_taylor)
 sum(0) = array_ints(0) 
 do i = 1, n_taylor
  ai = a0**dble(i)
  sum(i) = sum(i-1) + array_ints(i) * (-1.d0)**dble(i) * fact_inv(i) * ai
 enddo
 erf_exp_f_phi_ij = shank_general(sum,n_taylor,n_taylor) 

end

