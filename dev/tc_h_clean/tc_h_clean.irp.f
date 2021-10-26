program tc_h_clean
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
  print*,''
  print*,''
  print*,''
  print*,''
  print*,''
  print*,'eigval_trans = ',eigval_tc(1)
end
