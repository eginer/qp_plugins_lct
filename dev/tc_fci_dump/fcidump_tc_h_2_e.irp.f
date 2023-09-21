program tc_fci_dump
  implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC
 my_grid_becke = .True. 
 my_n_pt_r_grid = 30
 my_n_pt_a_grid = 50
 touch  my_grid_becke my_n_pt_r_grid my_n_pt_a_grid 
 call fcidump_2_tc_cst_mu
end
