program herm_Etc

  implicit none

  my_grid_becke = .True.
  my_n_pt_r_grid = 30
  my_n_pt_a_grid = 50
  touch  my_grid_becke my_n_pt_r_grid my_n_pt_a_grid

  read_wf = .True.
  touch read_wf

  PROVIDE N_int
  call print_energy

end


subroutine print_energy()

  implicit none
  integer          :: i, j
  double precision :: htilde_ij, hmono, heff, hderiv, hthree

  i = 1
  j = 1
  call htilde_mu_mat(psi_det(1,1,i), psi_det(1,1,j), N_int, hmono, heff, hderiv, hthree, htilde_ij)

  print *, '     hermitian energy: ', hmono + heff + hthree + nuclear_repulsion
  print *, ' non hermitian energy: ', hderiv
  print *, '         total energy: ', htilde_ij + nuclear_repulsion

  return
end subroutine print_energy
