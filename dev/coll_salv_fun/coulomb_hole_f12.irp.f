double precision function gama_f12_mu_based(mu)
 implicit none
 double precision, intent(in) :: mu
 BEGIN_DOC
! you enter with a mu, and it returns the gama exponent of the F12 jastrow of the form 
!
! f(r12,gama) =exp(-1/gama * exp(-gama * r12) )
!
! such that f(0,gama) = 1./(1. + 2./(sqrt(pi) * mu)) 
!
! which is the extrapolation factor of the on-top pair density
 END_DOC
  include 'constants.include.F'
 gama_f12_mu_based = dlog(1.d0 + 2.d0/(sqpi * mu))
 gama_f12_mu_based = 1.0/gama_f12_mu_based
 
end
