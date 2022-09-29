program print_etcscf

  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC

  implicit none

  print *, 'starting ...'

  my_grid_becke  = .True.
  my_n_pt_r_grid = 30
  my_n_pt_a_grid = 50
!  my_n_pt_r_grid = 10 ! small grid for quick debug
!  my_n_pt_a_grid = 26 ! small grid for quick debug
  touch my_grid_becke my_n_pt_r_grid my_n_pt_a_grid

  call routine_etcscf()

end

! ---

subroutine routine_etcscf()

  implicit none

  print*,'***'
  print*,'TC HF total energy = ', TC_HF_energy
  print*,'TC HF 1 e   energy = ', TC_HF_one_electron_energy
  print*,'TC HF 2 e   energy = ', TC_HF_two_e_energy

end subroutine routine_etcscf

! ---

