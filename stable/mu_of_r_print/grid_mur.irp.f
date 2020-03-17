 BEGIN_PROVIDER [double precision, cas_full_mu_of_r_grid_mur_psi_coal_vector, (n_points_print_mur) ]
 implicit none 
 BEGIN_DOC
 ! cas_full_mu_of_r obtained from the FULL interaction and the FULL two body density of the wave function stored in psi_det/psi_coef
 END_DOC
 integer :: i_point,k,i,j
 double precision :: r(3)
 double precision :: cpu0,cpu1,local_potential,two_body_dm
 print*,'providing the cas_full_mu_of_r_grid_mur_psi_coal_vector ...'
 call wall_time(cpu0)
 if(.True.)then
  provide total_cas_on_top_density_grid_mur 
  provide core_inact_act_grid_mur_f_psi_ab
 endif
 !$OMP PARALLEL DO &
 !$OMP DEFAULT (NONE)  &
 !$OMP PRIVATE (i_point,r,local_potential,two_body_dm) & 
 !$OMP shARED (n_points_print_mur,cas_full_mu_of_r_grid_mur_psi_coal_vector,core_inact_act_grid_mur_f_psi_ab,total_cas_on_top_density_grid_mur) 
 do i_point = 1, n_points_print_mur
  local_potential = core_inact_act_grid_mur_f_psi_ab(i_point) / total_cas_on_top_density_grid_mur(i_point,1)
  if(total_cas_on_top_density_grid_mur(i_point,1).gt.1.d-12.and.core_inact_act_grid_mur_f_psi_ab(i_point).gt.1.d-12)then
   local_potential = core_inact_act_grid_mur_f_psi_ab(i_point)/total_cas_on_top_density_grid_mur(i_point,1)
  else 
   local_potential = 1.d+10
  endif
  cas_full_mu_of_r_grid_mur_psi_coal_vector(i_point) =  local_potential * dsqrt(dacos(-1.d0)) * 0.5d0
 enddo
 !$OMP END PARALLEL DO


 call wall_time(cpu1)
 print*,'Time to provide cas_full_mu_of_r_grid_mur_psi_coal_vector = ',cpu1-cpu0
 END_PROVIDER 

