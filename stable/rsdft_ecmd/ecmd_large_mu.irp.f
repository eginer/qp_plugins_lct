 BEGIN_PROVIDER [double precision, ecmd_large_mu, (N_states)]
 BEGIN_DOC
  ! Eq. (23) of JCP, 150, 084103 (2019) : 
  ! Ec_md energy at large mu behaviour in function of the on top pair density.
  !
  ! ecmd_large_mu = (alpha/mu**3) * int n2(r,r) dr  where alpha = sqrt(2pi)*(-2+sqrt(2)) 
  !
  ! n2(r,r) is supposed to be the exact on-top pair density 
  ! 
  ! Here we use the extrapolated on-top pair density based on the asymptotic expansion at large mu
  !
  ! of P. Gori-Giorgi and A. Savin, Phys. Rev. A73, 032506 (2006)
 END_DOC
 implicit none 
 integer :: istate
 double precision :: pi,mu,on_top_extrap,alpha
 mu = mu_erf_dft
 pi = 4.d0 * datan(1.d0)
 alpha = (-2.d0+sqrt(2.d0))*sqrt(2.d0*pi)/(3.d0)

 do istate = 1, N_states
 ! Eq. 29 of JCP extrapolation between the on-top pair density at mu and the exact one
 ! the factor "2" comes from a difference of normalization of the two-body density in the JCP paper and QP2
  on_top_extrap = 2.d0 * integral_on_top(istate)/(1.d0 + 2.d0 / (dsqrt(pi)*mu ))
  ! constant in the large mu behaviour
  ! alpha / mu^3 \int n2(r,r)
  ecmd_large_mu(istate) = alpha / (mu**3) * on_top_extrap
 enddo
 END_PROVIDER

