
BEGIN_PROVIDER [double precision, mu_of_r_prov, (n_points_integration_angular,n_points_radial_grid,nucl_num) ]
 implicit none
 BEGIN_DOC
 ! value of mu(r) in each point in space
 END_DOC
 integer :: i_atom,k,l
  do i_atom = 1,nucl_num
   do k = 1, n_points_radial_grid -1
    do l = 1 , n_points_integration_angular
     mu_of_r_prov(l,k,i_atom) =  mu_of_r(l,k,i_atom)
    enddo
   enddo
  enddo
END_PROVIDER 


BEGIN_PROVIDER [double precision, mu_of_r_prov_selected, (n_points_final_grid) ]
 implicit none
 BEGIN_DOC
 ! value of mu(r) in each point in space
 END_DOC
 integer :: i_point
 integer :: i_atom,k,l
 double precision :: r(3),distance
 do i_point = 1, n_points_final_grid
  l = index_final_points(1,i_point)
  k = index_final_points(2,i_point)
  i_atom = index_final_points(3,i_point)
  r(1:3) = final_grid_points(1:3,i_point)  
  distance = dsqrt(r(1)**2 + r(2)**2 + r(3)**2)
  mu_of_r_prov_selected(i_point) = mu_of_r_prov(l,k,i_atom) !* distance 
 enddo

END_PROVIDER 
