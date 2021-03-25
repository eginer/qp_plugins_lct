program plot_density_w
 implicit none
 read_wf = .True.
 touch read_wf
 call routine_print_density
end

subroutine routine_print_density
 implicit none
 integer :: ipoint,nx
 double precision :: x,dx,xmax,xmin,r(3)
 double precision :: dm_a,dm_b
 xmax = 5.d0
 nx = 10000
 dx = xmax / dble(nx)
 xmin = -0.5d0 * xmax
 r = 0.d0
 r(1) = xmin

 character*(128) :: output
 integer :: i_unit_output,getUnitAndOpen
 output=trim(ezfio_filename)//'.density'
 i_unit_output = getUnitAndOpen(output,'w')

 do ipoint = 1, nx
  call dm_dft_alpha_beta_at_r(r,dm_a,dm_b)
  write(i_unit_output,'(100(F16.10,X))')r(1),dm_a + dm_b
  r(1) += dx
 enddo

end


double precision function density_from_matrix(matrix,r)
 implicit none
 double precision, intent(in) :: matrix(mo_num, mo_num), r(3)
 double precision :: mos_array(mo_num)
 call give_all_mos_at_r(r,mos_array)
 density_from_matrix = 0.d0
 integer :: i,j
 do i = 1, mo_num
  do j = 1, mo_num
   density_from_matrix += mos_array(i) * mos_array(j) * matrix(j,i)
  enddo
 enddo
end
