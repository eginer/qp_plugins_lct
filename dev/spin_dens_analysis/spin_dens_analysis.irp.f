program spin_dens_analysis
  implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC
  print *, 'Hello world'
  read_wf = .True.
  touch read_wf
  call routine
end

subroutine routine
 implicit none
 provide integrated_delta_rho_all_points
 provide integrated_rho_tot_all_points
 call print_mulliken_sd
 call print_standard_deviation
end
