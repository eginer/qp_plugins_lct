program save_right_eigv
 implicit none
 my_grid_becke = .True.
 my_n_pt_r_grid = 30
 my_n_pt_a_grid = 50
 touch  my_grid_becke my_n_pt_r_grid my_n_pt_a_grid
 read_wf = .True.
 touch read_wf 
 print*,'Saving right (and maybe left) eigenvector in EZFIO as the first (second) states'
 call routine_save_left_right
end


