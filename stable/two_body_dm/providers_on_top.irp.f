
 BEGIN_PROVIDER [double precision, on_top_of_r_vector,(n_points_final_grid,N_states) ]
 implicit none
 BEGIN_DOC
 ! On top pair density as a function of r given as a single array of r points
 ! Can use the regular two-body tensor or approximated versions
 END_DOC
 integer :: i_point,i_state
 if(on_top_from_cas)then
  do i_point = 1, n_points_final_grid
   do i_state = 1, N_states
    on_top_of_r_vector(i_point,i_state) = core_inact_act_on_top_of_r(i_point,i_state)
   enddo
  enddo
  return
 endif
 if(ontop_approx)then
  do i_point = 1, n_points_final_grid
   do i_state = 1, N_states
    on_top_of_r_vector(i_point,i_state) = on_top_of_r_approx_svd(i_point,i_state)
   enddo
  enddo
 else
  do i_point = 1, n_points_final_grid
   do i_state = 1, N_states
    on_top_of_r_vector(i_point,i_state) = on_top_of_r_exact(i_point,i_state)
   enddo
  enddo
 endif
 END_PROVIDER


 BEGIN_PROVIDER [double precision, on_top_of_r,(n_points_integration_angular,n_points_radial_grid,nucl_num,N_states) ]
 implicit none
 BEGIN_DOC
 ! On top pair density as a function of r given by the angular and nuclear points attached per atoms 
 END_DOC
 integer :: i_point,i_state
 integer :: i,j,k
 do i_state = 1, N_states
  do i_point = 1, n_points_final_grid
   k = index_final_points(1,i_point)
   i = index_final_points(2,i_point)
   j = index_final_points(3,i_point)
   on_top_of_r(k,i,j,i_state) = on_top_of_r_vector(i_point,i_state)
  enddo  
 enddo
 END_PROVIDER
