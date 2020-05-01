double precision function f_mu(mu,x)
 include 'constants.include.F'
 BEGIN_DOC
! correlation function corresponding to the erf(mu*x)/x interaction
 END_DOC
 implicit none
 double precision, intent(in) :: mu,x
 f_mu = 0.5d0 * x * (1.d0 - derf(mu*x)) - inv_sq_pi*0.5d0/mu * dexp(-x*mu*x*mu)
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

subroutine taylor_f_tilde(mu,x,g0,g1,g2,f_tilde)
 include 'constants.include.F'
 BEGIN_DOC
! Taylor expansion of the tilde_f_mu function
 END_DOC
 implicit none
 double precision, intent(in) :: mu,x
 double precision, intent(out):: g0,g1,g2,f_tilde
 g0 = -1.d0
 g1 = sqpi * mu
 g2 = -mu*mu + 0.5d0 * g1
 f_tilde = g0 + g1 * x + g2 * x*x 
end
