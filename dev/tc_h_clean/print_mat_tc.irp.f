program print_mat_tc
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC

  implicit none

  my_grid_becke = .True.
  my_n_pt_r_grid = 30
  my_n_pt_a_grid = 50
  touch  my_grid_becke my_n_pt_r_grid my_n_pt_a_grid

  read_wf = .True.
  touch read_wf 

  call routine
end

subroutine routine
 implicit none
 integer :: i
 print*,''
 print*,'Total H_TC matrix '
 do i= 1, N_det
  write(*,'(100(F12.8,X))')htilde_matrix_elmt(i,:)
 enddo
 print*,''
 print*,'3-e   H_TC matrix '
 do i= 1, N_det
  write(*,'(100(F12.8,X))')htilde_matrix_elmt_hthree(i,:)
 enddo
end
