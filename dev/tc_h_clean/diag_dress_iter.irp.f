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
 double precision, allocatable :: energies_out(:),u_in(:,:)
 logical :: converged
 allocate(energies_out(N_states_diag),u_in(N_det,N_states_diag))
 u_in = 0.d0
 u_in(1:N_det,1:N_states) = reigvec_tc(1:N_det,1:N_states)
 call iterative_davidson_tc(psi_det,u_in,N_det,N_states,N_states_diag,1,energies_out,converged)
end
