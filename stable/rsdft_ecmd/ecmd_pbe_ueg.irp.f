
 BEGIN_PROVIDER [double precision, ecmd_pbe_ueg_prov, (N_states)]
 BEGIN_DOC
  ! Give the Ec_md energy with a large mu behaviour in function of the UEG on top pair density coupled to the PBE correlation energy at mu=0
  ! Ec_md_PBE = Int epsilon_c_PBE_mu=0 / ( 1 + beta*mu**3 ) = Int eps_c_md_PBE  with beta chosen to recover the UEG large mu behaviour (JT)
 END_DOC
 implicit none
 double precision :: weight
 integer :: ipoint,istate
 double precision :: eps_c_md_PBE,mu,rho_a,rho_b,grad_rho_a(3),grad_rho_b(3),on_top
 double precision :: g0_UEG_mu_inf
 mu = mu_erf_dft
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

   call ec_md_pbe_on_top_general(mu,rho_a,rho_b,grad_rho_a,grad_rho_b,on_top,eps_c_md_PBE)

   ecmd_pbe_ueg_prov(istate) += eps_c_md_PBE * weight
  enddo
 enddo
 END_PROVIDER


