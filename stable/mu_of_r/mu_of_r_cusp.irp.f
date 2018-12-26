

 BEGIN_PROVIDER [double precision, mu_of_r_cusp_condition, (n_points_integration_angular,n_points_radial_grid,nucl_num) ]
 implicit none 
 BEGIN_DOC
 ! mu_of_r and mu_average computation 
 END_DOC
 integer :: i,j,k,i_point
 double precision :: r(3),two_dm,two_dm_laplacian,total_dm
 do i_point = 1, n_points_final_grid
  k = index_final_points(1,i_point)
  i = index_final_points(2,i_point)
  j = index_final_points(3,i_point)
  mu_of_r_cusp_condition(k,i,j) = mu_of_r_cusp_condition_vector(i_point,1)
 enddo

 END_PROVIDER 
