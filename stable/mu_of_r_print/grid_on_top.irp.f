 double precision function total_cas_on_top_density_grid_mur_from_provider(ipoint,i_state)
 implicit none
 BEGIN_DOC
 ! on top pair density evaluated at a given point of the grid 
 END_DOC
 integer, intent(in) :: ipoint,i_state
 integer :: i,j,k,l
 total_cas_on_top_density_grid_mur_from_provider = 0.d0
 do l = 1, n_core_inact_act_orb
  do k = 1, n_core_inact_act_orb
    do j = 1, n_core_inact_act_orb
     do i = 1, n_core_inact_act_orb
     !                                                                                          1 2 1 2 
     total_cas_on_top_density_grid_mur_from_provider += core_inact_act_two_bod_alpha_beta_mo_physicist(i,j,k,l,i_state) * core_inact_act_mos_in_r_array_grid_mur(j,ipoint) * core_inact_act_mos_in_r_array_grid_mur(i,ipoint) * core_inact_act_mos_in_r_array_grid_mur(l,ipoint) * core_inact_act_mos_in_r_array_grid_mur(k,ipoint)
    enddo
   enddo
  enddo
 enddo
 end

 BEGIN_PROVIDER [double precision, total_cas_on_top_density_grid_mur,(n_points_print_mur,N_states) ]
&BEGIN_PROVIDER [double precision, wall_time_total_cas_on_top_density_grid_mur ]
 implicit none
 BEGIN_DOC
 ! on top pair density at each grid point computed using the full two-body density matrix 
 END_DOC
 integer :: i_point,i_state
 double precision :: wall_0,wall_1
 double precision :: total_cas_on_top_density_grid_mur_from_provider

 print*,'providing the total_cas_on_top_density_grid_mur'
 i_point = 1
 provide core_inact_act_two_bod_alpha_beta_mo_physicist
 i_state = 1
 total_cas_on_top_density_grid_mur(i_point,i_state) = total_cas_on_top_density_grid_mur_from_provider(i_point,i_state)
 call wall_time(wall_0)
 !$OMP PARALLEL DO &
 !$OMP DEFAULT (NONE)  &
 !$OMP PRIVATE (i_point,i_state) & 
 !$OMP SHARED(total_cas_on_top_density_grid_mur,n_points_print_mur,N_states)
 do i_point = 1, n_points_print_mur
  do i_state = 1, N_states
   total_cas_on_top_density_grid_mur(i_point,i_state) = total_cas_on_top_density_grid_mur_from_provider(i_point,i_state)
  enddo
 enddo
 !$OMP END PARALLEL DO
 call wall_time(wall_1)
 print*,'provided the total_cas_on_top_density_grid_mur'
 print*,'Time to provide :',wall_1 - wall_0
 wall_time_total_cas_on_top_density_grid_mur = wall_1 - wall_0

 END_PROVIDER 


 BEGIN_PROVIDER [double precision, extrapolated_core_inact_act_on_top_grid_mur,(n_points_print_mur,N_states) ]
 implicit none
 BEGIN_DOC
 ! on top pair density at each grid point computed using the full two-body density matrix 
 END_DOC
 integer :: i_point,i_state
 double precision :: two_dm,mu,on_top_two_dm_in_r_mu_corrected_from_two_dm
 do i_point = 1, n_points_print_mur
  do i_state = 1, N_states
   two_dm = total_cas_on_top_density_grid_mur(i_point,i_state) 
   mu = cas_full_mu_of_r_grid_mur_psi_coal_vector(i_point)
   extrapolated_core_inact_act_on_top_grid_mur(i_point,i_state) = on_top_two_dm_in_r_mu_corrected_from_two_dm(mu,i_state,two_dm)
  enddo
 enddo
 END_PROVIDER 

 BEGIN_PROVIDER [double precision, UEG_on_top_grid_mur,(n_points_print_mur,N_states) ]
 implicit none
 BEGIN_DOC
 ! on top pair density at each grid point computed using the full two-body density matrix 
 END_DOC

 integer :: i_point,i_state
 double precision :: rho_a(N_states),rho_b(N_states),g0_UEG_mu_inf
 do i_point = 1, n_points_print_mur
  call dm_dft_alpha_beta_at_r(grid_points_mur(1,i_point),rho_a,rho_b)
  do i_state = 1, N_states
   UEG_on_top_grid_mur(i_point,i_state) = 2.d0 * rho_a(i_state)*rho_b(i_state)*g0_UEG_mu_inf(rho_a(i_state),rho_b(i_state))
  enddo
 enddo

 END_PROVIDER 
