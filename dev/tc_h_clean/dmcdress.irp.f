program dmcdress

  implicit none

  my_grid_becke = .True.
  my_n_pt_r_grid = 30
  my_n_pt_a_grid = 50
  touch  my_grid_becke my_n_pt_r_grid my_n_pt_a_grid

  read_wf = .True.
  touch read_wf

  PROVIDE N_int
  call delta_dmcdressing

end


subroutine delta_dmcdressing()

  implicit none
  double precision, allocatable :: delta(:) 

  print *, j1b_gauss
  print *, j1b_gauss_pen
  print *, mu_erf

  allocate( delta(N_det) )

  ! get < I | H_mu - H | psi > 
  call get_delta_tc_psi(psi_det, psi_coef, N_det, N_int, delta)

  ! order as QMCCHEM
  call dset_order(delta, psi_bilinear_matrix_order, N_det)

  call ezfio_set_dmc_dress_dmc_delta_h(delta)

  deallocate(delta)

  return
end subroutine delta_dmcdressing
