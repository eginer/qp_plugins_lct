program print_mat
 implicit none
 read_wf = .True.
 touch read_wf
 call routine

end

subroutine routine
 implicit none
 integer :: i,j,k,l
 j = 1
 l = 1
 call get_mo_two_e_integrals_i1j1(j,l,mo_num,integrals_matrix,mo_integrals_map)
 call compute_all_ijkl_for_jl_mu_of_r_int(j,l,eff_int_mat)
 if(.True.)then
  do i = 1, N_det
   write(*,'(100(F14.7,X))')H_matrix_all_dets(i,:)
  enddo
  print*,''
  print*,'Delta E = ',H_matrix_all_dets(3,3) - H_matrix_all_dets(1,1)
  double precision :: integrals_matrix(mo_num,mo_num),eff_int_mat(mo_num,mo_num)
  print*,'(11|11),eff(11|11) = ',integrals_matrix(1,1),eff_int_mat(1,1)
  print*,'(11|22),eff(11|22) = ',integrals_matrix(2,2),eff_int_mat(2,2)
 endif

end
