
 BEGIN_PROVIDER [double precision, Energy_c_md_on_top_PBE_mu_UEG_vector, (N_states)]
 BEGIN_DOC
  ! Give the Ec_md energy with a good large mu behaviour in function of the on top pair density with mu correction based on the on-top of the UEG coupled to the PBE correlation energy at mu=0
  ! Ec_md_on_top_PBE = Int epsilon_c_PBE_mu=0 / ( 1 + beta*mu**3 ) = Int eps_c_md_on_top_PBE  with beta chosen to recover the good large mu behaviour of the Energy_c_md_on_top functional
 END_DOC
 implicit none
 double precision ::  r(3)
 double precision :: weight,mu
 integer :: i,istate
 double precision,allocatable  :: eps_c_md_on_top_PBE(:),two_dm(:)
 allocate(eps_c_md_on_top_PBE(N_states),two_dm(N_states))
 mu = mu_erf_dft
 Energy_c_md_on_top_PBE_mu_UEG_vector = 0.d0
  
 do i = 1, n_points_final_grid
  r(1) = final_grid_points(1,i)
  r(2) = final_grid_points(2,i)
  r(3) = final_grid_points(3,i)
  weight = final_weight_at_r_vector(i)
  two_dm(:) = total_cas_on_top_density(i,:)
  call give_epsilon_c_md_on_top_PBE_mu_corrected_UEG_from_two_dm(mu,r,two_dm,eps_c_md_on_top_PBE)
  do istate = 1, N_states
   Energy_c_md_on_top_PBE_mu_UEG_vector(istate) += eps_c_md_on_top_PBE(istate) * weight
  enddo
 enddo
 END_PROVIDER

