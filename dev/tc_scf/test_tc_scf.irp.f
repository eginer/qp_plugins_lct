program tc_scf

  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC

  implicit none

  print *, 'starting ...'

  my_grid_becke  = .True.
  my_n_pt_r_grid = 30
  my_n_pt_a_grid = 50
  touch my_grid_becke my_n_pt_r_grid my_n_pt_a_grid

  call routine_test

end

subroutine routine_test
 implicit none
 integer :: i,j 
 do i = 1, mo_num
  print*,eigval_f_mat_ortho_basis(i),eigval_fock_tc_mo(i),dabs(eigval_f_mat_ortho_basis(i)-eigval_fock_tc_mo(i))
 enddo
 do i = 1, mo_num
  do j = 1, mo_num
   if(i==j)cycle
   if(dabs(overlap_fock_tc_eigvec_mo(i,j)).gt.1.d-8)then
    print*,i,j,overlap_fock_tc_eigvec_mo(i,j)
   endif
  enddo
 enddo
 double precision :: accu
 accu = 0.d0
 do i = 1, mo_num
  do j = 1, mo_num
   accu += dabs(overlap_rleigv_f_mat_ortho_basis(j,i)) - dabs(overlap_fock_tc_eigvec_mo(j,i))
  enddo
 enddo
 print*,'accu = ',accu
 print*,'overlap_fock_tc_eigvec_mo'
 do i = 1, mo_num
  write(*,'(100(F8.4,X))')overlap_fock_tc_eigvec_mo(:,i)
 enddo
 print*,'overlap_rleigv_f_mat_ortho_basis'
 do i = 1, mo_num
  write(*,'(100(F8.4,X))')overlap_rleigv_f_mat_ortho_basis(:,i)
 enddo
end
