subroutine give_laplacian_2_all_aos(r,lapl_2_aos)
 implicit none
 double precision, intent(in) :: r(3)
 double precision, intent(out):: lapl_2_aos(ao_num)

 double precision :: aos_array(ao_num)
 double precision :: aos_grad_array(3,ao_num)
 double precision :: aos_lapl_array(3,ao_num)
 double precision :: aos_3rd_direct_array(3,ao_num)
 double precision :: aos_4th_direct_array(3,ao_num)
 double precision :: aos_2nd_cross_array(3,ao_num)
 double precision :: aos_3rd_cross_array(6,ao_num)
 double precision :: aos_4th_cross_array(3,ao_num)
 call give_all_aos_and_fourth_at_r(r,aos_array,aos_grad_array,aos_lapl_array,aos_3rd_direct_array,aos_4th_direct_array)
 call give_all_aos_and_fourth_order_cross_terms_at_r(r,aos_array,aos_grad_array,aos_2nd_cross_array,aos_3rd_cross_array,aos_4th_cross_array)
 integer :: i
 do i = 1, ao_num
  lapl_2_aos(i)  =       aos_4th_direct_array(1,i) + aos_4th_direct_array(2,i) + aos_4th_direct_array(3,i) 
  lapl_2_aos(i) += 2.d0*(aos_4th_cross_array(1,i) + aos_4th_cross_array(2,i) + aos_4th_cross_array(3,i)) 
 enddo

end

 BEGIN_PROVIDER[double precision, aos_lapl_2_in_r_array, (ao_num,n_points_final_grid)]
&BEGIN_PROVIDER[double precision, aos_lapl_2_in_r_array_transp, (n_points_final_grid,ao_num)]
 integer :: i,j
 double precision :: lapl_2_aos(ao_num), r(3)

 do i = 1, n_points_final_grid
  r(1) = final_grid_points(1,i)
  r(2) = final_grid_points(2,i)
  r(3) = final_grid_points(3,i)
  call give_laplacian_2_all_aos(r,lapl_2_aos)
  do j = 1, ao_num
   aos_lapl_2_in_r_array(j,i) = lapl_2_aos(j)
   aos_lapl_2_in_r_array_transp(i,j) = lapl_2_aos(j)
  enddo
 enddo

END_PROVIDER 

 BEGIN_PROVIDER[double precision, mos_lapl_2_in_r_array, (mo_num,n_points_final_grid)]
&BEGIN_PROVIDER[double precision, mos_lapl_2_in_r_array_transp, (n_points_final_grid,mo_num)]
 integer :: i,j,k
 double precision :: lapl_2_aos(ao_num), r(3)
 mos_lapl_2_in_r_array = 0.d0
 do i = 1, n_points_final_grid
  r(1) = final_grid_points(1,i)
  r(2) = final_grid_points(2,i)
  r(3) = final_grid_points(3,i)
  do k = 1, mo_num
   do j = 1, ao_num
    mos_lapl_2_in_r_array(k,i)  += aos_lapl_2_in_r_array(j,i) * mo_coef(j,k)
   enddo
  enddo
  do k = 1, mo_num
   mos_lapl_2_in_r_array_transp(i,k) = mos_lapl_2_in_r_array(k,i)
  enddo
 enddo

END_PROVIDER 

