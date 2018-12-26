
 BEGIN_PROVIDER [double precision, mu_of_r_vector, (n_points_final_grid) ]
&BEGIN_PROVIDER [double precision, mu_of_r, (n_points_integration_angular,n_points_radial_grid,nucl_num) ]
 implicit none 
 BEGIN_DOC
 ! mu_of_r and mu_average computation 
 END_DOC
 integer :: i_point,k,i,j
 double precision :: r(3)
 double precision :: cpu0,cpu1
 print*,'providing the mu_of_r ...'
 call wall_time(cpu0)
 r = 0.d0
 do i_point = 1, n_points_final_grid
  r(1) = final_grid_points(1,i_point)
  r(2) = final_grid_points(2,i_point)
  r(3) = final_grid_points(3,i_point)
  if(mu_of_r_potential.EQ."cusp_condition")then
   mu_of_r_vector(i_point) = mu_of_r_cusp_condition_vector(i_point,1)
  else if(mu_of_r_potential.EQ."hf_coallescence")then
   mu_of_r_vector(i_point) =  mu_of_r_hf_coalescence_vector(i_point)
  else if(mu_of_r_potential.EQ."psi_coallescence")then
   mu_of_r_vector(i_point) =  mu_of_r_psi_coalescence_vector(i_point)
  else if(mu_of_r_potential.EQ."hf_integral")then
   mu_of_r_vector(i_point) = mu_of_r_integral_hf_vector(i_point)
  else 
    print*,'you requested the following mu_of_r_potential'
    print*,mu_of_r_potential
    print*,'which does not correspond to any of the options for such keyword'
    stop
  endif
  k = index_final_points(1,i_point)
  i = index_final_points(2,i_point)
  j = index_final_points(3,i_point)
  mu_of_r(k,i,j) = mu_of_r_vector(i_point)
 enddo
 call wall_time(cpu1)
 print*,'Time to provide mu_of_r = ',cpu1-cpu0
 END_PROVIDER 
