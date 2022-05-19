program bi_ortho_mos
  implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC
  my_grid_becke = .True.
  my_n_pt_r_grid = 30
  my_n_pt_a_grid = 50
  touch  my_grid_becke my_n_pt_r_grid my_n_pt_a_grid
  call test_overlap
end

subroutine test_overlap                                                                                                   
 implicit none
 integer :: i,j
 provide overlap_bi_ortho
end

