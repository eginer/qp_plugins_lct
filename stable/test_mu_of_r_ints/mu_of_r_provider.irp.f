
BEGIN_PROVIDER [double precision, mu_of_r_for_ints_vector, (n_points_final_grid) ]
 implicit none
 BEGIN_DOC
 ! value of mu(r) in each point in space
 END_DOC
 integer :: i_point
 integer :: i_atom,k,l
 double precision :: r(3)
 do i_point = 1, n_points_final_grid
  l = index_final_points(1,i_point)
  k = index_final_points(2,i_point)
  i_atom = index_final_points(3,i_point)
  mu_of_r_for_ints_vector(i_point) = 1.d0
 enddo

END_PROVIDER 
