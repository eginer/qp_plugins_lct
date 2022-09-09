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

  bi_ortho = .True.
  touch bi_ortho
  call print_fock

end

subroutine print_fock
 implicit none
 integer :: i,j
 print*,'Diagonal of Fock matrix'
 do i = 1, mo_num
  print*,i,Fock_matrix_tc_mo_tot(i,i)
 enddo
 print*,'Fock matrix '
 do i = 1, mo_num
  write(*,'(100(F16.10,X))')Fock_matrix_tc_mo_tot(:,i)
 enddo
end

