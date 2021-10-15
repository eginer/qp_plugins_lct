program write_htilde_mat
 implicit none
 read_wf = .True.
 touch read_wf
 my_grid_becke = .True. 
 my_n_pt_r_grid = 30
 my_n_pt_a_grid = 50
 touch  my_grid_becke my_n_pt_r_grid my_n_pt_a_grid 
 call routine
end

subroutine routine
 implicit none
 integer :: i,j
 do i = 1, N_det
  write(33,'(10000(F16.10,X))')htilde_matrix_elmt(i,:)
 enddo
end
