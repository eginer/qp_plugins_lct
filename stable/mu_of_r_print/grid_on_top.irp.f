 double precision function core_inact_act_on_top_of_r_grid_mur_from_provider(ipoint,istate)
 implicit none
 BEGIN_DOC
 ! on top pair density evaluated at a given point of the grid 
 END_DOC
 integer, intent(in) :: ipoint,istate
 integer :: i,j,k,l
 core_inact_act_on_top_of_r_grid_mur_from_provider = 0.d0
 do l = 1, n_core_inact_act_orb
  do k = 1, n_core_inact_act_orb
    do j = 1, n_core_inact_act_orb
     do i = 1, n_core_inact_act_orb
     !                                                                                          1 2 1 2 
     core_inact_act_on_top_of_r_grid_mur_from_provider += core_inact_act_two_bod_alpha_beta_mo_physicist(i,j,k,l,istate) * core_inact_act_mos_in_r_array_grid_mur(j,ipoint) * core_inact_act_mos_in_r_array_grid_mur(i,ipoint) * core_inact_act_mos_in_r_array_grid_mur(l,ipoint) * core_inact_act_mos_in_r_array_grid_mur(k,ipoint)
    enddo
   enddo
  enddo
 enddo
 end

 BEGIN_PROVIDER [double precision, core_inact_act_on_top_of_r_grid_mur,(n_points_final_grid,N_states) ]
&BEGIN_PROVIDER [double precision, wall_time_core_inact_act_on_top_of_r_grid_mur ]
 implicit none
 BEGIN_DOC
 ! on top pair density at each grid point computed using the full two-body density matrix 
 END_DOC
 integer :: i_point,i_state
 double precision :: wall_0,wall_1
 double precision :: core_inact_act_on_top_of_r_grid_mur_from_provider

 print*,'providing the core_inact_act_on_top_of_r_grid_mur'
 i_point = 1
 provide core_inact_act_two_bod_alpha_beta_mo_physicist
 i_state = 1
 core_inact_act_on_top_of_r_grid_mur(i_point,i_state) = core_inact_act_on_top_of_r_grid_mur_from_provider(i_point,i_state)
 call wall_time(wall_0)
 !$OMP PARALLEL DO &
 !$OMP DEFAULT (NONE)  &
 !$OMP PRIVATE (i_point,i_state) & 
 !$OMP SHARED(core_inact_act_on_top_of_r_grid_mur,n_points_final_grid,N_states)
 do i_point = 1, n_points_final_grid
  do i_state = 1, N_states
   core_inact_act_on_top_of_r_grid_mur(i_point,i_state) = core_inact_act_on_top_of_r_grid_mur_from_provider(i_point,i_state)
  enddo
 enddo
 !$OMP END PARALLEL DO
 call wall_time(wall_1)
 print*,'provided the core_inact_act_on_top_of_r_grid_mur'
 print*,'Time to provide :',wall_1 - wall_0
 wall_time_core_inact_act_on_top_of_r_grid_mur = wall_1 - wall_0

 END_PROVIDER 

