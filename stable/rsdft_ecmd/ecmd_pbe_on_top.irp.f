 BEGIN_PROVIDER [double precision, ecmd_pbe_on_top_at_mu, (N_states)]
 BEGIN_DOC
!
! Ecmd functional defined in Eq. (25) of JCP, 150, 084103 1-10 (2019) 
! 
! Such a functional is built by interpolating between two regimes : 
! 
!    +) the large mu behaviour in cst/(\mu^3) \int dr on-top(r) where on-top(r) is supposed to be the exact on-top of the system
!
!    +) mu= 0 with the usal ec_pbe(rho_a,rho_b,grad_rho_a,grad_rho_b) 
!
! Here the approximation to the exact on-top is done through the assymptotic expansion (in \mu) of the exact on-top pair density (see Eq. 29)
! 
! Such an asymptotic expansion was introduced in P. Gori-Giorgi and A. Savin, Phys. Rev. A73, 032506 (2006)
!
 END_DOC
 implicit none
 double precision :: weight
 double precision  :: eps_c_md_on_top_PBE,on_top_extrap,mu_correction_of_on_top
 integer :: ipoint,istate
 double precision :: eps_c_md_PBE,mu,rho_a,rho_b,grad_rho_a(3),grad_rho_b(3),on_top
 mu = mu_erf_dft
 ecmd_pbe_on_top_at_mu = 0.d0
  
 do istate = 1, N_states
  do ipoint = 1, n_points_final_grid
   weight = final_weight_at_r_vector(ipoint)
   
   rho_a = one_e_dm_and_grad_alpha_in_r(4,ipoint,istate)
   rho_b = one_e_dm_and_grad_beta_in_r(4,ipoint,istate)
   grad_rho_a(1:3) = one_e_dm_and_grad_alpha_in_r(1:3,ipoint,istate)
   grad_rho_b(1:3) = one_e_dm_and_grad_beta_in_r(1:3,ipoint,istate)
   
!  We take the extrapolated on-top pair density (Eq. 29)
   on_top = total_cas_on_top_density(ipoint,istate)
!  Multiplied by 2 because of difference of normalizations between the on_top of QP2 and that of JCP, 150, 084103 1-10 (2019)
   on_top_extrap = 2.d0 * mu_correction_of_on_top(mu,on_top)

   call ec_md_pbe_on_top_general(mu,rho_a,rho_b,grad_rho_a,grad_rho_b,on_top_extrap,eps_c_md_on_top_PBE)
 
   ecmd_pbe_on_top_at_mu(istate) += eps_c_md_on_top_PBE * weight
  enddo
 enddo
 END_PROVIDER


