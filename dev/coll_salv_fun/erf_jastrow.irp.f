double precision function f_mu(mu,x)
 include 'constants.include.F'
 BEGIN_DOC
! correlation function corresponding to the erf(mu*x)/x interaction
 END_DOC
 implicit none
 double precision, intent(in) :: mu,x
 f_mu = 0.5d0 * x * (1.d0 - derf(mu*x)) - inv_sq_pi*0.5d0/mu * dexp(-x*mu*x*mu)
end

double precision function full_jastrow_mu(x,mu)
 implicit none
 double precision, intent(in) :: x,mu
 double precision :: f_mu
 full_jastrow_mu = dexp(f_mu(x,mu))
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

double precision function jastrow_func_tilde(x,mu)
 implicit none
 double precision, intent(in) :: x,mu
 double precision :: tilde_f_mu
 jastrow_func_tilde = dexp(tilde_f_mu(x,mu))
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
 f_tilde_fit = slater_fit_gam(x,gama) * dexp(-delta *x*x)
end

double precision function ovlp_f_tilde_phi_ij(mu,r1,A_center,B_center,power_A,power_B,alpha,beta,n_taylor)
 implicit none
 double precision, intent(in)    :: r1(3), mu        ! where the Jastrow factor is centered and its shape 
 integer, intent(in)             :: n_taylor         ! order of the Taylor expansion of the exponential
 double precision, intent(in)    :: A_center(3),B_center(3),alpha,beta ! gaussian/polynoms "A" and "B"
 integer, intent(in)             :: power_A(3),power_B(3)
 
 double precision :: gama, delta
 double precision :: array_ints(0:n_taylor),integral_taylor

 call give_param_fit_f_h(mu,gama,delta)
 call exp_dl_ovlp_stg_phi_ij(r1,gama,delta,A_center,B_center,power_A,power_B,alpha,beta,n_taylor,array_ints,integral_taylor)
 ovlp_f_tilde_phi_ij = integral_taylor

end

