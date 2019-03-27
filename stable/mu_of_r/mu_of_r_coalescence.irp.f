
 BEGIN_PROVIDER [double precision, mu_of_r_hf_coal_vector, (n_points_final_grid) ]
&BEGIN_PROVIDER [double precision, mu_of_r_hf_coal, (n_points_integration_angular,n_points_radial_grid,nucl_num) ]
 implicit none 
 BEGIN_DOC
 ! mu_of_r obtained from the FULL interaction and the FULL two body density of the HF wave function
 END_DOC
 integer :: i_point,k,i,j
 double precision :: r(3)
 double precision :: cpu0,cpu1,local_potential,two_bod
 print*,'providing the mu_of_r_hf_coal_vector ...'
 call wall_time(cpu0)
 r = 0.d0
 call f_HF_ab(r,r,local_potential,two_bod)
 !$OMP PARALLEL DO &
 !$OMP DEFAULT (NONE)  &
 !$OMP PRIVATE (i_point,r,local_potential,two_bod) & 
 !$OMP ShARED (n_points_final_grid,final_grid_points,mu_of_r_hf_coal_vector) 
 do i_point = 1, n_points_final_grid
  r(1) = final_grid_points(1,i_point)
  r(2) = final_grid_points(2,i_point)
  r(3) = final_grid_points(3,i_point)
  call f_HF_ab(r,r,local_potential,two_bod)
  if(two_bod.le.1.d-12.or.local_potential.le.0.d0.or.local_potential * two_bod.lt.0.d0)then
    local_potential = 1.d+10
  else 
    local_potential = local_potential /  two_bod
  endif
  mu_of_r_hf_coal_vector(i_point) =  local_potential * dsqrt(dacos(-1.d0)) * 0.5d0
 enddo
 !$OMP END PARALLEL DO
 do i_point = 1, n_points_final_grid
  k = index_final_points(1,i_point)
  i = index_final_points(2,i_point)
  j = index_final_points(3,i_point)
  mu_of_r_hf_coal(k,i,j) = mu_of_r_hf_coal_vector(i_point)
 enddo
 call wall_time(cpu1)
 print*,'Time to provide mu_of_r_hf_coal_vector = ',cpu1-cpu0
 END_PROVIDER 

 BEGIN_PROVIDER [double precision, mu_of_r_psi_coal_vector, (n_points_final_grid) ]
&BEGIN_PROVIDER [double precision, mu_of_r_psi_coal, (n_points_integration_angular,n_points_radial_grid,nucl_num) ]
 implicit none 
 BEGIN_DOC
 ! mu_of_r obtained from the FULL interaction and the FULL two body density of the wave function stored in psi_det/psi_coef
 END_DOC
 integer :: i_point,k,i,j
 double precision :: r(3)
 double precision :: cpu0,cpu1,local_potential,two_body_dm
 print*,'providing the mu_of_r_psi_coal_vector ...'
 call wall_time(cpu0)

 if(.True.)then
  provide on_top_of_r_vector 
  provide f_psi_ab
 endif
 !$OMP PARALLEL DO &
 !$OMP DEFAULT (NONE)  &
 !$OMP PRIVATE (i_point,r,local_potential,two_body_dm) & 
 !$OMP shARED (n_points_final_grid,final_grid_points,mu_of_r_psi_coal_vector,f_psi_ab,on_top_of_r_vector) 
 do i_point = 1, n_points_final_grid
  local_potential = f_psi_ab(i_point) / on_top_of_r_vector(i_point,1)
  if(on_top_of_r_vector(i_point,1).gt.1.d-12.and.f_psi_ab(i_point).gt.1.d-12)then
   local_potential = f_psi_ab(i_point)/on_top_of_r_vector(i_point,1)
  else 
   local_potential = 1.d-10
  endif
  mu_of_r_psi_coal_vector(i_point) =  local_potential * dsqrt(dacos(-1.d0)) * 0.5d0
 enddo
 !$OMP END PARALLEL DO


 do i_point = 1, n_points_final_grid
  k = index_final_points(1,i_point)
  i = index_final_points(2,i_point)
  j = index_final_points(3,i_point)
  mu_of_r_psi_coal(k,i,j) = mu_of_r_psi_coal_vector(i_point)
 enddo
 call wall_time(cpu1)
 print*,'Time to provide mu_of_r_psi_coal_vector = ',cpu1-cpu0
 END_PROVIDER 


 BEGIN_PROVIDER [double precision, mu_of_r_hf_coal_vv_vector, (n_points_final_grid) ]
&BEGIN_PROVIDER [double precision, mu_of_r_hf_coal_vv, (n_points_integration_angular,n_points_radial_grid,nucl_num) ]
&BEGIN_PROVIDER [double precision, mu_average_hf_coal_vv, ( N_states )  ]
 implicit none 
 BEGIN_DOC
 ! mu_of_r obtained from the valence-valence interaction and the valence-valence two body density of the HF wave function
 END_DOC
 integer :: i_point,k,i,j
 double precision :: r(3)
 double precision :: cpu0,cpu1,local_potential,two_bod
 print*,'providing the mu_of_r_hf_coal_vv_vector ...'
 call wall_time(cpu0)
 r = 0.d0
 call f_HF_ab(r,r,local_potential,two_bod)
 !$OMP PARALLEL DO &
 !$OMP DEFAULT (NONE)  &
 !$OMP PRIVATE (i_point,r,local_potential,two_bod) & 
 !$OMP ShARED (n_points_final_grid,final_grid_points,mu_of_r_hf_coal_vv_vector) 
 do i_point = 1, n_points_final_grid
  r(1) = final_grid_points(1,i_point)
  r(2) = final_grid_points(2,i_point)
  r(3) = final_grid_points(3,i_point)
  call f_HF_valence_ab(r,r,local_potential,two_bod)
  if(two_bod.le.1.d-12.or.local_potential.le.0.d0.or.local_potential * two_bod.lt.0.d0)then
    local_potential = 1.d+10
  else 
    local_potential = local_potential /  two_bod
  endif
  mu_of_r_hf_coal_vv_vector(i_point) =  local_potential * dsqrt(dacos(-1.d0)) * 0.5d0
 enddo
 !$OMP END PARALLEL DO
 do i_point = 1, n_points_final_grid
  k = index_final_points(1,i_point)
  i = index_final_points(2,i_point)
  j = index_final_points(3,i_point)
  mu_of_r_hf_coal_vv(k,i,j) = mu_of_r_hf_coal_vv_vector(i_point)
 enddo
 call wall_time(cpu1)
 print*,'Time to provide mu_of_r_hf_coal_vv_vector = ',cpu1-cpu0
 mu_average_hf_coal_vv = 0.d0
 double precision :: elec_num_val_tmp(N_states),weight
 integer :: istate 
 elec_num_val_tmp = 0.d0
 do i_point = 1, n_points_final_grid
  weight=final_weight_at_r_vector(i_point)
  do istate = 1, N_states 
   elec_num_val_tmp(istate) += ( one_e_dm_no_core_and_grad_alpha_in_r(4,i_point,istate) + one_e_dm_no_core_and_grad_beta_in_r(4,i_point,istate)) * weight 
   if( mu_of_r_hf_coal_vv_vector(i_point) .gt.1000d0)cycle
   mu_average_hf_coal_vv(istate) += mu_of_r_hf_coal_vv_vector(i_point) * & 
                                   ( one_e_dm_no_core_and_grad_alpha_in_r(4,i_point,istate) + one_e_dm_no_core_and_grad_beta_in_r(4,i_point,istate)) * weight 
  enddo
 enddo
 do istate = 1, N_states
  mu_average_hf_coal_vv(istate) = mu_average_hf_coal_vv(istate) / elec_num_val_tmp(istate)
 enddo
 END_PROVIDER 


 BEGIN_PROVIDER [double precision, mu_of_r_hf_coal_cc_vector, (n_points_final_grid) ]
&BEGIN_PROVIDER [double precision, mu_of_r_hf_coal_cc, (n_points_integration_angular,n_points_radial_grid,nucl_num) ]
&BEGIN_PROVIDER [double precision, mu_average_hf_coal_cc, ( N_states )  ]
 implicit none 
 BEGIN_DOC
 ! mu_of_r obtained from the core-core interaction and the core-core two body density of the HF wave function
 END_DOC
 integer :: i_point,k,i,j
 double precision :: r(3)
 double precision :: cpu0,cpu1,local_potential,two_bod
 print*,'providing the mu_of_r_hf_coal_cc_vector ...'
 call wall_time(cpu0)
 r = 0.d0
 call f_HF_ab(r,r,local_potential,two_bod)
 !$OMP PARALLEL DO &
 !$OMP DEFAULT (NONE)  &
 !$OMP PRIVATE (i_point,r,local_potential,two_bod) & 
 !$OMP ShARED (n_points_final_grid,final_grid_points,mu_of_r_hf_coal_cc_vector) 
 do i_point = 1, n_points_final_grid
  r(1) = final_grid_points(1,i_point)
  r(2) = final_grid_points(2,i_point)
  r(3) = final_grid_points(3,i_point)
  call f_HF_core_ab(r,r,local_potential,two_bod)
  if(two_bod.le.1.d-12.or.local_potential.le.0.d0.or.local_potential * two_bod.lt.0.d0)then
    local_potential = 1.d+10
  else 
    local_potential = local_potential /  two_bod
  endif
  mu_of_r_hf_coal_cc_vector(i_point) =  local_potential * dsqrt(dacos(-1.d0)) * 0.5d0
 enddo
 !$OMP END PARALLEL DO
 do i_point = 1, n_points_final_grid
  k = index_final_points(1,i_point)
  i = index_final_points(2,i_point)
  j = index_final_points(3,i_point)
  mu_of_r_hf_coal_cc(k,i,j) = mu_of_r_hf_coal_cc_vector(i_point)
 enddo
 call wall_time(cpu1)
 print*,'Time to provide mu_of_r_hf_coal_cc_vector = ',cpu1-cpu0

 mu_average_hf_coal_cc = 0.d0
 double precision :: elec_num_val_tmp(N_states),weight
 integer :: istate 
 elec_num_val_tmp = 0.d0
 do i_point = 1, n_points_final_grid
  weight=final_weight_at_r_vector(i_point)
  do istate = 1, N_states 
   elec_num_val_tmp(istate) += ( one_e_dm_no_core_and_grad_alpha_in_r(4,i_point,istate) + one_e_dm_no_core_and_grad_beta_in_r(4,i_point,istate)) * weight 
   if( mu_of_r_hf_coal_cc_vector(i_point) .gt.1000d0)cycle
   mu_average_hf_coal_cc(istate) += mu_of_r_hf_coal_cc_vector(i_point) * & 
                                   ( one_e_dm_no_core_and_grad_alpha_in_r(4,i_point,istate) + one_e_dm_no_core_and_grad_beta_in_r(4,i_point,istate)) * weight 
  enddo
 enddo
 do istate = 1, N_states
  mu_average_hf_coal_cc(istate) = mu_average_hf_coal_cc(istate) / elec_num_val_tmp(istate)
 enddo
 END_PROVIDER 


 BEGIN_PROVIDER [double precision, mu_of_r_hf_coal_cv_vector, (n_points_final_grid) ]
&BEGIN_PROVIDER [double precision, mu_of_r_hf_coal_cv, (n_points_integration_angular,n_points_radial_grid,nucl_num) ]
&BEGIN_PROVIDER [double precision, mu_average_hf_coal_cv, ( N_states )  ]
 implicit none 
 BEGIN_DOC
 ! mu_of_r obtained from the core-valence interaction and the core-valence two body density of the HF wave function
 END_DOC
 integer :: i_point,k,i,j
 double precision :: r(3)
 double precision :: cpu0,cpu1,local_potential,two_bod
 print*,'providing the mu_of_r_hf_coal_cv_vector ...'
 call wall_time(cpu0)
 r = 0.d0
 call f_HF_ab(r,r,local_potential,two_bod)
 !$OMP PARALLEL DO &
 !$OMP DEFAULT (NONE)  &
 !$OMP PRIVATE (i_point,r,local_potential,two_bod) & 
 !$OMP ShARED (n_points_final_grid,final_grid_points,mu_of_r_hf_coal_cv_vector) 
 do i_point = 1, n_points_final_grid
  r(1) = final_grid_points(1,i_point)
  r(2) = final_grid_points(2,i_point)
  r(3) = final_grid_points(3,i_point)
  call f_HF_core_valence_ab(r,r,local_potential,two_bod)
  if(two_bod.le.1.d-12.or.local_potential.le.0.d0.or.local_potential * two_bod.lt.0.d0)then
    local_potential = 1.d+10
  else 
    local_potential = local_potential /  two_bod
  endif
  mu_of_r_hf_coal_cv_vector(i_point) =  local_potential * dsqrt(dacos(-1.d0)) * 0.5d0
 enddo
 !$OMP END PARALLEL DO
 do i_point = 1, n_points_final_grid
  k = index_final_points(1,i_point)
  i = index_final_points(2,i_point)
  j = index_final_points(3,i_point)
  mu_of_r_hf_coal_cv(k,i,j) = mu_of_r_hf_coal_cv_vector(i_point)
 enddo
 call wall_time(cpu1)
 print*,'Time to provide mu_of_r_hf_coal_cv_vector = ',cpu1-cpu0

 mu_average_hf_coal_cv = 0.d0
 double precision :: elec_num_val_tmp(N_states),weight
 integer :: istate 
 elec_num_val_tmp = 0.d0
 do i_point = 1, n_points_final_grid
  weight=final_weight_at_r_vector(i_point)
  do istate = 1, N_states 
   elec_num_val_tmp(istate) += ( one_e_dm_no_core_and_grad_alpha_in_r(4,i_point,istate) + one_e_dm_no_core_and_grad_beta_in_r(4,i_point,istate)) * weight 
   if( mu_of_r_hf_coal_cv_vector(i_point) .gt.1000d0)cycle
   mu_average_hf_coal_cv(istate) += mu_of_r_hf_coal_cv_vector(i_point) * & 
                                   ( one_e_dm_no_core_and_grad_alpha_in_r(4,i_point,istate) + one_e_dm_no_core_and_grad_beta_in_r(4,i_point,istate)) * weight 
  enddo
 enddo
 do istate = 1, N_states
  mu_average_hf_coal_cv(istate) = mu_average_hf_coal_cv(istate) / elec_num_val_tmp(istate)
 enddo
 END_PROVIDER 

