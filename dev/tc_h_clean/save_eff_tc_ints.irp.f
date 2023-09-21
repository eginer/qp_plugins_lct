program save_tc_ints
  implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC
  my_grid_becke = .True.
  my_n_pt_r_grid = 30
  my_n_pt_a_grid = 50
  touch  my_grid_becke my_n_pt_r_grid my_n_pt_a_grid
  zero_tc_eff_map = .True.
  touch zero_tc_eff_map 
  call save_eff_tc_two_e_ints_mo_into_ints_mo

end
