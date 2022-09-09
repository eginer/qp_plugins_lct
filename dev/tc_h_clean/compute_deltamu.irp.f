program compute_deltamu

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

! ---

subroutine delta_dmcdressing()

  implicit none
  integer                       :: k
  double precision, allocatable :: delta(:,:) 

  print *, j1b_gauss
  print *, j1b_gauss_pen
  print *, mu_erf

  allocate( delta(N_det,N_states) )
  delta = 0.d0

  do k = 1, N_states
  !do k = 1, 1

    ! get < I | H_mu - H | psi > 
    call get_delta_tc_psi(psi_det, psi_coef(:,k), N_det, N_int, delta(:,k))

    ! order as QMCCHEM
    call dset_order(delta(:,k), psi_bilinear_matrix_order, N_det)

  enddo

  call ezfio_set_dmc_dress_dmc_delta_h(delta)

  deallocate(delta)

  return
end subroutine delta_dmcdressing

! ---

