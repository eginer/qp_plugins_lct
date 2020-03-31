program write_mu_of_r
 implicit none
 read_wf = .True.
 touch read_wf
 call routine

end


subroutine routine
 implicit none
 provide mu_of_r_vector
 call ezfio_set_becke_numerical_grid_n_points_final_grid(n_points_final_grid)
 call ezfio_set_mu_of_r_mu_of_r_array(mu_of_r_vector)


end
