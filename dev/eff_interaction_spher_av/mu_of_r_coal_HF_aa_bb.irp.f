!!!!!!!!! alpha alpha part!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 BEGIN_PROVIDER [double precision, mu_of_r_hf_coal_vector_alpha, (n_points_final_grid) ]
&BEGIN_PROVIDER [double precision, mu_of_r_hf_coal_alpha, (n_points_integration_angular,n_points_radial_grid,nucl_num) ]
 implicit none
 BEGIN_DOC
 ! mu_of_r HF alpha alpha 
 END_DOC
 integer :: i_point,k,i,j
 double precision :: r(3)
 double precision :: mu
 double precision :: cpu0,cpu1,local_potential,two_bod,f_deriv2,f_deriv4,n2_deriv2,n2_deriv4,wee_paper
 print*,'providing the mu_of_r_hf_coal_vector_alpha ...'
 call wall_time(cpu0)
 r = 0.d0
 call give_eff_inter_alpha_alpha_hf_at_r1_r12(r,delta_r12,local_potential,two_bod)
!!$OMP PARALLEL DO &
!!$OMP DEFAULT (NONE)  &
!!$OMP PRIVATE (i_point,r,local_potential,two_bod,delta_r12,mu,f_deriv2,f_deriv4,n2_deriv2,n2_deriv4) & 
!!$OMP ShARED (n_points_final_grid,final_grid_points,mu_of_r_hf_coal_vector_alpha) 
 do i_point = 1, n_points_final_grid
  r(1) = final_grid_points(1,i_point)
  r(2) = final_grid_points(2,i_point)
  r(3) = final_grid_points(3,i_point)
  call give_f_paper_alpha_alpha_hf_at_r1_r12(r,delta_r12,local_potential,f_deriv2,f_deriv4)
  call give_n2_alpha_alpha_hf_at_r1_r12(r,delta_r12,two_bod,n2_deriv2,n2_deriv4)
  if(dabs(two_bod).lt.1.d-40)then
   mu_of_r_hf_coal_vector_alpha(i_point) = 1.d+20
  else if(local_potential/two_bod*delta_r12.ge.0.9999d0)then
   mu_of_r_hf_coal_vector_alpha(i_point) = 1.d+20
  else if(local_potential/two_bod.lt.0.d0)then
   mu_of_r_hf_coal_vector_alpha(i_point) = 1.d+20
  else 
   local_potential = local_potential /  two_bod
   call give_mu_r12(local_potential,delta_r12,mu)
   mu_of_r_hf_coal_vector_alpha(i_point) = mu 
  endif
 enddo
!!$OMP END PARALLEL DO
 do i_point = 1, n_points_final_grid
  k = index_final_points(1,i_point)
  i = index_final_points(2,i_point)
  j = index_final_points(3,i_point)
  mu_of_r_hf_coal_alpha(k,i,j) = mu_of_r_hf_coal_vector_alpha(i_point)
 enddo
 call wall_time(cpu1)
 print*,'Time to provide mu_of_r_hf_coal_vector_alpha = ',cpu1-cpu0
 END_PROVIDER 

 BEGIN_PROVIDER [double precision, mu_average_aa, (N_states)]
 implicit none
 BEGIN_DOC
 ! mu_average = \int dr mu(r) * n(r) /N_elec
 END_DOC
 integer :: i,istate
 double precision, allocatable :: rho_a(:), rho_b(:)
 double precision :: mu,weight
 mu_average_aa = 0.d0
 allocate(rho_a(N_states), rho_b(N_states))
 do i = 1, n_points_final_grid
  mu = mu_of_r_hf_coal_vector_alpha(i)
  if(mu.gt.1.d10)cycle
  weight=final_weight_at_r_vector(i)
  do istate = 1, N_states
   rho_a(istate) = one_e_dm_alpha_at_r(i,istate)
   rho_b(istate) = 0.d0
   mu_average_aa(istate) +=  weight *  mu * (rho_a(istate))
  enddo
 enddo
 mu_average_aa = mu_average_aa / dble(elec_alpha_num)
 END_PROVIDER

!!!!!!!!! beta beta part!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 BEGIN_PROVIDER [double precision, mu_of_r_hf_coal_vector_beta, (n_points_final_grid) ]
&BEGIN_PROVIDER [double precision, mu_of_r_hf_coal_beta, (n_points_integration_angular,n_points_radial_grid,nucl_num) ]
 implicit none
 BEGIN_DOC
 ! mu_of_r HF beta beta 
 END_DOC
 integer :: i_point,k,i,j
 double precision :: r(3)
 double precision :: mu
 double precision :: cpu0,cpu1,local_potential,two_bod,f_deriv2,f_deriv4,n2_deriv2,n2_deriv4,wee_paper
 print*,'providing the mu_of_r_hf_coal_vector_beta ...'
 call wall_time(cpu0)
 r = 0.d0
 call give_eff_inter_beta_beta_hf_at_r1_r12(r,delta_r12,local_potential,two_bod)
!!$OMP PARALLEL DO &
!!$OMP DEFAULT (NONE)  &
!!$OMP PRIVATE (i_point,r,local_potential,two_bod,delta_r12,mu,f_deriv2,f_deriv4,n2_deriv2,n2_deriv4) & 
!!$OMP ShARED (n_points_final_grid,final_grid_points,mu_of_r_hf_coal_vector_beta) 
 do i_point = 1, n_points_final_grid
  r(1) = final_grid_points(1,i_point)
  r(2) = final_grid_points(2,i_point)
  r(3) = final_grid_points(3,i_point)
  call give_f_paper_beta_beta_hf_at_r1_r12(r,delta_r12,local_potential,f_deriv2,f_deriv4)
  call give_n2_beta_beta_hf_at_r1_r12(r,delta_r12,two_bod,n2_deriv2,n2_deriv4)
  if(dabs(two_bod).lt.1.d-30)then
   mu_of_r_hf_coal_vector_beta(i_point) = 1.d+20
  else if(local_potential/two_bod *delta_r12.ge.0.9999d0)then
   mu_of_r_hf_coal_vector_beta(i_point) = 1.d+20
  else if(local_potential/two_bod .lt.0.d0)then
   mu_of_r_hf_coal_vector_beta(i_point) = 1.d+20
  else 
   local_potential = local_potential / two_bod
   call give_mu_r12(local_potential,delta_r12,mu)
   mu_of_r_hf_coal_vector_beta(i_point) = mu 
  endif
 enddo
!!$OMP END PARALLEL DO
 do i_point = 1, n_points_final_grid
  k = index_final_points(1,i_point)
  i = index_final_points(2,i_point)
  j = index_final_points(3,i_point)
  mu_of_r_hf_coal_beta(k,i,j) = mu_of_r_hf_coal_vector_beta(i_point)
 enddo
 call wall_time(cpu1)
 print*,'Time to provide mu_of_r_hf_coal_vector_beta = ',cpu1-cpu0
 END_PROVIDER 

 BEGIN_PROVIDER [double precision, mu_average_bb, (N_states)]
 implicit none
 BEGIN_DOC
 ! mu_average = \int dr mu(r) * n(r) /N_elec
 END_DOC
 integer :: i,istate
 double precision, allocatable :: rho_a(:), rho_b(:)
 double precision :: mu,weight
 mu_average_bb = 0.d0
 allocate(rho_a(N_states), rho_b(N_states))
 do i = 1, n_points_final_grid
  mu = mu_of_r_hf_coal_vector_beta(i)
  if(mu.gt.1.d10)cycle
  weight=final_weight_at_r_vector(i)
  do istate = 1, N_states
   rho_a(istate) = 0.d0 
   rho_b(istate) = one_e_dm_beta_at_r(i,istate)
   mu_average_bb(istate) +=  weight *  mu * (rho_b(istate))
  enddo
 enddo
 mu_average_bb = mu_average_bb / dble(elec_beta_num)
 END_PROVIDER

