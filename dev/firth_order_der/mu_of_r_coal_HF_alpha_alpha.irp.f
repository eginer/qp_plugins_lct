 BEGIN_PROVIDER [double precision, mu_of_r_hf_coal_vector_alpha, (n_points_final_grid) ]
&BEGIN_PROVIDER [double precision, mu_of_r_hf_coal_alpha, (n_points_integration_angular,n_points_radial_grid,nucl_num) ]
 implicit none
 BEGIN_DOC
 ! mu_of_r HF alpha alpha 
 END_DOC
 integer :: i_point,k,i,j
 double precision :: r(3)
 double precision :: r12,mu
 double precision :: cpu0,cpu1,local_potential,two_bod
 print*,'providing the mu_of_r_hf_coal_vector_alpha ...'
 call wall_time(cpu0)
 r = 0.d0
 r12 = 5.d-4
 call give_eff_inter_alpha_alpha_hf_at_r1_r12(r,r12,local_potential,two_bod)
!!$OMP PARALLEL DO &
!!$OMP DEFAULT (NONE)  &
!!$OMP PRIVATE (i_point,r,local_potential,two_bod,r12,mu) & 
!!$OMP ShARED (n_points_final_grid,final_grid_points,mu_of_r_hf_coal_vector_alpha) 
 do i_point = 1, n_points_final_grid
  r12 = 1.d-4
  r(1) = final_grid_points(1,i_point)
  r(2) = final_grid_points(2,i_point)
  r(3) = final_grid_points(3,i_point)
  call give_eff_inter_alpha_alpha_hf_at_r1_r12(r,r12,local_potential,two_bod)
  if(two_bod.le.1.d-12.or.local_potential.le.0.d0)then
    local_potential = 1.d-10
  else 
    local_potential = local_potential /  two_bod
  endif
  call give_mu_r12(local_potential,r12,mu)
  if(local_potential.lt.1.d-10)then
   print*,'r'
   print*, r 
   print*,'local_potential, dm, mu = '
   print*,local_potential,one_e_dm_alpha_at_r(i_point,1),mu
   pause
  endif
  mu_of_r_hf_coal_vector_alpha(i_point) = mu 
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
  weight=final_weight_at_r_vector(i)
  do istate = 1, N_states
   rho_a(istate) = one_e_dm_alpha_at_r(i,istate)
   rho_b(istate) = 0.d0
   mu_average_aa(istate) +=  weight *  mu * (rho_a(istate))
  enddo
 enddo
 mu_average_aa = mu_average_aa / dble(elec_alpha_num)
 END_PROVIDER

