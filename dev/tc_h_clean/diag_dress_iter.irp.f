program diag_dress_iter
  implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC
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
 integer :: i,j
 print*,'eigval_right_tc = ',eigval_right_tc
 print*,'eigval_left _tc = ',eigval_left_tc
 print*,'Left, right and usual eigenvectors '
 do i = 1, N_det
  write(*,'(I5,X,(100(F9.5,X)))')i,leigvec_tc(i,1),reigvec_tc(i,1),psi_coef(i,1)
 enddo
end
