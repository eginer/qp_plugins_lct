 BEGIN_PROVIDER [double precision, ecmd_large_mu, (N_states)]
 BEGIN_DOC
  ! Give the Ec_md energy with a good large mu behaviour in function of the on top pair density.
  ! Ec_md_on_top = (alpha/mu**3) * int n2(r,r) dr  where alpha = sqrt(2pi)*(-2+sqrt(2)) 
 END_DOC
 implicit none 
 integer :: istate
 double precision :: pi,mu
 mu = mu_erf_dft
 pi = 4.d0 * datan(1.d0)
 ecmd_large_mu = ((-2.d0+sqrt(2.d0))*sqrt(2.d0*pi)/(3.d0*(mu**3)))*integral_on_top/(1.d0 + 2.d0 / (dsqrt(pi)*mu ))
 END_PROVIDER


