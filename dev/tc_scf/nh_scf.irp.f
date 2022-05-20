program nh_scf

  BEGIN_DOC 
  ! non-hermitian scf
  END_DOC

  implicit none

  my_grid_becke = .True.
  my_n_pt_r_grid = 30
  my_n_pt_a_grid = 50
  touch  my_grid_becke my_n_pt_r_grid my_n_pt_a_grid

  call run

end

!__________________________________________________________________________________________________
!
subroutine run

  BEGIN_DOC
  ! run non-hermitian SCF calculation
  END_DOC

  implicit none

  mo_label = 'Orthonormalized'

  !thresh_SCF = 1d-2
  !TOUCH thresh_SCF
  !call Roothaan_Hall_SCF
  !call ezfio_set_hartree_fock_energy(SCF_energy)

  !thresh_SCF = 1d-6
  !TOUCH thresh_SCF

  call Roothaan_Hall_nhSCF
  !call ezfio_set_hartree_fock_energy(SCF_energy)

end subroutine run
!__________________________________________________________________________________________________


