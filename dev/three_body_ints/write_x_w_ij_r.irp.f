program write_x_w_ij_r
 implicit none
 my_grid_becke = .True.
 my_n_pt_r_grid = 30
 my_n_pt_a_grid = 50
 touch  my_grid_becke my_n_pt_r_grid my_n_pt_a_grid
 call routine_x_wij_r
 call routine_mos
end

subroutine routine_x_wij_r
 implicit none
 integer :: ipoint,m,i,j
 double precision :: wall0,wall1
 double precision,allocatable :: weight(:)
 character*(128) :: output
 integer :: i_unit_output
 integer :: getUnitAndOpen
 do ipoint = 1, n_points_final_grid
  final_weight_at_r_vector(ipoint)
 enddo
 output =trim(ezfio_filename)//'/x_w_ij_r'
 i_unit_output = getUnitAndOpen(output,'w')
 call wall_time(wall0)
 do i = 1, mo_num
  do j = 1, mo_num
   do m = 1, 3
    do ipoint = 1, n_points_final_grid
     write(i_unit_output,*)ipoint,m,j,i,x_W_ij_erf_rk(ipoint,m,j,i)
    enddo
   enddo
  enddo
 enddo
 call wall_time(wall1)
 print*,'Time to write x_W_ij_erf_rk = ',wall1 - wall0

end

subroutine routine_mos
 implicit none
 integer :: ipoint,i
 double precision :: wall0,wall1
 character*(128) :: output
 integer :: i_unit_output
 integer :: getUnitAndOpen
 output =trim(ezfio_filename)//'/mos_in_r'
 i_unit_output = getUnitAndOpen(output,'w')
 call wall_time(wall0)
 do i = 1, mo_num
  do ipoint = 1, n_points_final_grid
   write(i_unit_output,*)ipoint,i,mos_in_r_array_transp(ipoint,i)
  enddo
 enddo
 call wall_time(wall1)
 print*,'Time to write mos_in_r_array_transp = ',wall1 - wall0

end


