BEGIN_PROVIDER [double precision, mo_ij_at_r, (mo_num*mo_num,n_points_final_grid)]
 implicit none
 integer :: i,j,ij,ig
 do ig = 1, n_points_final_grid
  ij = 0
  do i = 1,mo_num
   do j=1,mo_num
    ij += 1
    mo_ij_at_r(ij,ig) = mos_in_r_array(i,ig) * mos_in_r_array(j,ig)
   enddo
  enddo
 enddo

END_PROVIDER 
