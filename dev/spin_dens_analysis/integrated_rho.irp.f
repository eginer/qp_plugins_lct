program integrated_rho
  implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC
  read_wf = .True.
  touch read_wf
  call routine
end

subroutine routine
 implicit none
 provide integrated_rho_tot_all_points
end
