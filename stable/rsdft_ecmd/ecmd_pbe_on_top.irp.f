 BEGIN_PROVIDER [double precision, ecmd_pbe_on_top_at_mu, (N_states)]
 BEGIN_DOC
  ! Give the Ec_md energy with a good large mu behaviour in function of the on top pair density with mu correction based on the on-top of the UEG coupled to the PBE correlation energy at mu=0
  ! Ec_md_on_top_PBE = Int epsilon_c_PBE_mu=0 / ( 1 + beta*mu**3 ) = Int eps_c_md_on_top_PBE  with beta chosen to recover the good large mu behaviour of the Energy_c_md_on_top functional
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
   
!  We take the extrpolated on-top pair density * 2 because of normalization
   on_top = total_cas_on_top_density(ipoint,istate)
   on_top_extrap = 2.d0 * mu_correction_of_on_top(mu,on_top)

   call ec_md_pbe_on_top_general(mu,rho_a,rho_b,grad_rho_a,grad_rho_b,on_top_extrap,eps_c_md_on_top_PBE)
 
   ecmd_pbe_on_top_at_mu(istate) += eps_c_md_on_top_PBE * weight
  enddo
 enddo
 END_PROVIDER


