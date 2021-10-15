program plot_density_tc
 implicit none
 read_wf = .True.
 touch read_wf
 call routine_check_tm
 call routine_print_density
 call print_z_dipole_tc_moment_only
end

subroutine routine_check_tm
 implicit none
 integer :: i
 print*,''
 print*,''
 print*,'one_e_tm_mo_norm   = ',one_e_tm_mo_norm
!print*,'left_right_overlap = ',left_right_overlap(1,1)
!print*,'tmp                = ',left_right_overlap_read
 print*,'one_e_tm_mo_norm   = ',one_e_tm_mo_norm 
 print*,'elec_num           = ',elec_num_tab(1) + elec_num_tab(2)
 print*,''
 print*,''
 print*,'right_overlap_read = ',right_overlap_read
 print*,'left_overlap_read  = ',left_overlap_read
 print*,''
 print*,'Printing the one_e_tm_mo'
 do i = 1, mo_num
  write(*,'(1000(F16.10,X))')one_e_tm_mo(i,:,1)
 enddo
 print*,''
 print*,''
 print*,''
end 

subroutine routine_print_density
 implicit none
 integer :: ipoint,i,j,m1,m2
 integer :: nx
 double precision :: x,dx,xmax,xmin,r(3)
 double precision :: density_from_matrix,dm_r,dm_l, tm_lr
 double precision, allocatable :: dm_mo(:,:,:)
 allocate(dm_mo(mo_num, mo_num, N_states)) ! state 1 == Right, state 2 == Left
 do i = 1, mo_num
  do j = 1, mo_num
   dm_mo(j,i,1) = one_e_dm_mo_alpha(j,i,1) + one_e_dm_mo_beta(j,i,1)
   dm_mo(j,i,2) = one_e_dm_mo_alpha(j,i,2) + one_e_dm_mo_beta(j,i,2)
  enddo
 enddo
 xmax = 5.d0
 nx = 10000
 dx = xmax / dble(nx)
 xmin = -0.5d0 * xmax
 r = 0.d0
 r(1) = xmin

 character*(128) :: output
 integer :: i_unit_output,getUnitAndOpen
 output=trim(ezfio_filename)//'.tc_density'
 i_unit_output = getUnitAndOpen(output,'w')

 do ipoint = 1, nx
  dm_r = density_from_matrix(dm_mo(1,1,1),r) 
  dm_l = density_from_matrix(dm_mo(1,1,2),r) 
  tm_lr= density_from_matrix(one_e_tm_mo(1,1,1),r)
  write(i_unit_output,'(100(F16.10,X))')r(1),dm_r,dm_l,tm_lr
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
