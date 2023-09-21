
 BEGIN_PROVIDER [double precision, ecmd_pbe_ueg_prov, (N_states)]
 BEGIN_DOC
!
! Ecmd functional using the on-top pair density of the UEG. 
! 
! Obtained with the functional originally introduced in JCP, 150, 084103 1-10 (2019) which interpolates between 
! 
!    +) the large mu behaviour in cst/(\mu^3) \int dr on-top(r) where on-top(r) is supposed to be the exact on-top of the system
!
!    +) mu= 0 with the usal ec_pbe(rho_a,rho_b,grad_rho_a,grad_rho_b) 
!
! Here the approximation to the exact on-top is done through the UEG
!
 END_DOC
 implicit none
 double precision :: weight
 integer :: ipoint,istate
 double precision :: eps_c_md_PBE,rho_a,rho_b,grad_rho_a(3),grad_rho_b(3),on_top
 double precision :: g0_UEG_mu_inf
 ecmd_pbe_ueg_prov = 0.d0
  
 do istate = 1, N_states
  do ipoint = 1, n_points_final_grid
   weight = final_weight_at_r_vector(ipoint)
   
   rho_a = one_e_dm_and_grad_alpha_in_r(4,ipoint,istate)
   rho_b = one_e_dm_and_grad_beta_in_r(4,ipoint,istate)
   grad_rho_a(1:3) = one_e_dm_and_grad_alpha_in_r(1:3,ipoint,istate)
   grad_rho_b(1:3) = one_e_dm_and_grad_beta_in_r(1:3,ipoint,istate)
   
!  We take the on-top pair density of the UEG which is (1-zeta^2) rhoc^2 g0 = 4 rhoa * rhob * g0
   on_top = 4.d0 * rho_a * rho_b * g0_UEG_mu_inf(rho_a,rho_b)   

   double precision :: mu_local 
   mu_local = mu_of_r_dft(ipoint)
   call ec_md_pbe_on_top_general(mu_local,rho_a,rho_b,grad_rho_a,grad_rho_b,on_top,eps_c_md_PBE)

   ecmd_pbe_ueg_prov(istate) += eps_c_md_PBE * weight
  enddo
 enddo
 END_PROVIDER


