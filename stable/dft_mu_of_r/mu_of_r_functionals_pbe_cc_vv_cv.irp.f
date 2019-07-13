BEGIN_PROVIDER [double precision, e_cmd_mu_of_r_pbe_vv, (N_states)]
 BEGIN_DOC
  ! e_cmd_mu_of_r_pbe_vv_vector           = PBE multi determinant functional with UEG on top for large mu using a mu(r) interaction (JT)
 END_DOC
 implicit none
 double precision ::  r(3)
 double precision :: weight,mu
 integer :: i,istate
 double precision              :: eps_c_md_PBE(N_states)
 double precision              :: rho_a(N_states),rho_b(N_states), grad_rho_a(3,N_states),grad_rho_b(3,N_states)

 e_cmd_mu_of_r_pbe_vv = 0.d0
  
 print*,'Providing e_cmd_mu_of_r_pbe_vv ...'
 call wall_time(wall0)
 do i = 1, n_points_final_grid
  r(1) = final_grid_points(1,i)
  r(2) = final_grid_points(2,i)
  r(3) = final_grid_points(3,i)
  weight = final_weight_at_r_vector(i)
  mu = mu_of_r_hf_coal_vv_vector(i)
  rho_a(:) = one_e_dm_no_core_and_grad_alpha_in_r(4,i,:)
  rho_b(:) = one_e_dm_no_core_and_grad_beta_in_r(4,i,:)
  grad_rho_a(1:3,:) = one_e_dm_no_core_and_grad_alpha_in_r(1:3,i,:)
  grad_rho_b(1:3,:) = one_e_dm_no_core_and_grad_beta_in_r(1:3,i,:)
  call give_epsilon_c_md_PBE_mu_grad_input(mu,rho_a,rho_b, grad_rho_a, grad_rho_b,eps_c_md_PBE)
  do istate = 1, N_states
   e_cmd_mu_of_r_pbe_vv(istate) += eps_c_md_PBE(istate) * weight
  enddo
 enddo
 double precision :: wall1, wall0
 call wall_time(wall1)
 print*,'Time for the e_cmd_mu_of_r_pbe_vv:',wall1-wall0

 END_PROVIDER
 
BEGIN_PROVIDER [double precision, e_cmd_mu_of_r_pbe_cc, (N_states)]
 BEGIN_DOC
  ! e_cmd_mu_of_r_pbe_cc_vector           = PBE multi determinant functional with UEG on top for large mu using a mu(r) interaction (JT)
 END_DOC
 implicit none
 double precision ::  r(3)
 double precision :: weight,mu
 integer :: i,istate
 double precision              :: eps_c_md_PBE(N_states)
 double precision              :: rho_a(N_states),rho_b(N_states), grad_rho_a(3,N_states),grad_rho_b(3,N_states)

 e_cmd_mu_of_r_pbe_cc = 0.d0
  
 print*,'Providing e_cmd_mu_of_r_pbe_cc ...'
 call wall_time(wall0)
 do i = 1, n_points_final_grid
  r(1) = final_grid_points(1,i)
  r(2) = final_grid_points(2,i)
  r(3) = final_grid_points(3,i)
  weight = final_weight_at_r_vector(i)
  mu = mu_of_r_hf_coal_cc_vector(i)
  rho_a(:) = one_e_dm_no_core_and_grad_alpha_in_r(4,i,:)
  rho_b(:) = one_e_dm_no_core_and_grad_beta_in_r(4,i,:)
  grad_rho_a(1:3,:) = one_e_dm_no_core_and_grad_alpha_in_r(1:3,i,:)
  grad_rho_b(1:3,:) = one_e_dm_no_core_and_grad_beta_in_r(1:3,i,:)
  call give_epsilon_c_md_PBE_mu_grad_input(mu,rho_a,rho_b, grad_rho_a, grad_rho_b,eps_c_md_PBE)
  do istate = 1, N_states
   e_cmd_mu_of_r_pbe_cc(istate) += eps_c_md_PBE(istate) * weight
  enddo
 enddo
 double precision :: wall1, wall0
 call wall_time(wall1)
 print*,'Time for the e_cmd_mu_of_r_pbe_cc:',wall1-wall0

 END_PROVIDER
