program write_x_w_ij_r
 implicit none
 my_grid_becke = .True.
 my_n_pt_r_grid = 30
 my_n_pt_a_grid = 50
 touch  my_grid_becke my_n_pt_r_grid my_n_pt_a_grid
 call routine_write_grid_pt_mos
 call routine_x_wij_r
 call routine_mos
end

subroutine routine_write_grid_pt_mos
 implicit none
 character*(128) :: output
 integer :: i_unit_output
 integer :: getUnitAndOpen
 output =trim(ezfio_filename)//'/n_grid_pts_mo_num'
 i_unit_output = getUnitAndOpen(output,'w')
 write(i_unit_output,*)n_points_final_grid
 write(i_unit_output,*)n_act_orb 
end

subroutine routine_x_wij_r
 implicit none
 integer :: ipoint,m,i,j,ii,jj
 double precision :: wall0,wall1
 double precision,allocatable :: weight(:)
 character*(128) :: output
 integer :: i_unit_output
 integer :: getUnitAndOpen
! do ipoint = 1, n_points_final_grid
!  final_weight_at_r_vector(ipoint)
! enddo
 output =trim(ezfio_filename)//'/x_w_ij_r'
 i_unit_output = getUnitAndOpen(output,'w')
 call wall_time(wall0)
 do ii = 1, n_act_orb
  i = list_act(ii)
  do jj = 1, n_act_orb 
   j = list_act(jj)
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
 integer :: ipoint,i,ii
 double precision :: wall0,wall1
 character*(128) :: output
 integer :: i_unit_output
 integer :: getUnitAndOpen
 output =trim(ezfio_filename)//'/mos_in_r'
 i_unit_output = getUnitAndOpen(output,'w')
 call wall_time(wall0)
 do ii = 1, n_act_orb
  i = list_act(ii)
  do ipoint = 1, n_points_final_grid
   write(i_unit_output,*)ipoint,i,mos_in_r_array_transp(ipoint,i) * sqrt_weight_at_r(ipoint)
  enddo
 enddo
 call wall_time(wall1)
 print*,'Time to write mos_in_r_array_transp = ',wall1 - wall0

end


