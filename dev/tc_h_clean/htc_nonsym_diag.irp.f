program htc_nonsym_diag

  BEGIN_DOC
    ! TODO : Put the documentation of the program here
  END_DOC

  implicit none

  my_grid_becke = .True.
  my_n_pt_r_grid = 30
  my_n_pt_a_grid = 50
  touch  my_grid_becke my_n_pt_r_grid my_n_pt_a_grid

  read_wf = .True.
  touch read_wf

  PROVIDE N_int
  PROVIDE n_states
  PROVIDE n_states_diag
  PROVIDE N_det
  PROVIDE psi_det
  PROVIDE psi_coef 

  !call check() 
  call nonsym_test() 

end 


subroutine nonsym_test()

  implicit none

  integer                       :: i, j
  logical                       :: converged
  double precision, allocatable :: u_in(:,:), H_jj(:), wr(:), wl(:)
  external                         htc_calc_tdav 
  external                         htcdag_calc_tdav

  allocate( H_jj(N_det) )
  H_jj = 0.d0
  i = 1
  call htilde_mu_mat_tot(psi_det(1,1,i), psi_det(1,1,i), N_int, H_jj(i))
  do i = 1, N_det
    call htilde_mu_mat_tot(psi_det(1,1,i), psi_det(1,1,i), N_int, H_jj(i))
  enddo

  allocate( u_in(N_det,n_states_diag) )
  u_in = 0.d0
  do i = 1, n_states 
    do j = 1, N_det
      u_in(j,i) = psi_coef(j,i)
    enddo
  enddo

  ! get right eigenvectors
  allocate( wr(n_states_diag) )
  print*, ' start right-eig calc'
  call davidson_general_ext_rout_nonsym_b1space(u_in, H_jj, wr, N_det, n_states, n_states_diag, converged, htc_calc_tdav)

  do i = 1, n_states
    print *, ' right eigen-energy = ', wr(i)
    do j = 1, N_det
      !reigvec_tc(j,i) = u_in(j,i)
      write(111, *) u_in(j,i)
    enddo
  enddo

  ! ---

  ! get left eigenvectors
  allocate( wl(n_states_diag) )
  print*, ' start left-eig calc'
  call davidson_general_ext_rout_nonsym_b1space(u_in, H_jj, wl, N_det, n_states, n_states_diag, converged, htcdag_calc_tdav)

  do i = 1, n_states
    print *, ' left eigen-energy = ', wl(i)
    do j = 1, N_det
      !leigvec_tc(j,i) = u_in(j,i)
      write(112, *) u_in(j,i)
    enddo
  enddo

  print *, ' | right - left | eigen-energy = ', dabs(wr(1)-wl(1))

  ! ---

  deallocate(wr, wl, u_in, H_jj) 


end subroutine nonsym_test

! ---

subroutine check()

  implicit none

  integer                       :: i, j, n_real_eigv
  double precision, allocatable :: H_tmp(:,:), reigvec(:,:), leigvec(:,:), eigval(:)

  allocate(H_tmp(N_det,N_det))
!  do i = 1, N_det
!    do j = 1, N_det
!      H_tmp(j,i) += H_matrix_all_dets(j,i)
!    enddo
!  enddo
!  do i = 1, N_det-1, 2
!    do j = i+1, N_det, 3
!      H_tmp(j,i) += 0.1d0 * dble(j)
!    enddo
!  enddo
  i = 1
  j = 1
  call htilde_mu_mat_tot(psi_det(1,1,i), psi_det(1,1,j), N_int, H_tmp(i,j))
  do i = 1, N_det
    do j = 1, N_det
      H_tmp(i,j) = 0.d0
      call htilde_mu_mat_tot(psi_det(1,1,i), psi_det(1,1,j), N_int, H_tmp(i,j))
    enddo
  enddo


  allocate( reigvec(N_det,N_det), leigvec(N_det,N_det), eigval(N_det) )
  call non_hrmt_real_diag(N_det, H_tmp, leigvec, reigvec, n_real_eigv, eigval)

  deallocate(H_tmp)

  do i = 1, n_states
    print *, ' exact eigen-energy = ', eigval(1)
    do j = 1, N_det
      write(211, *) reigvec_tc(j,1)
      write(212, *) leigvec_tc(j,1)
    enddo
  enddo
 
  deallocate( eigval, reigvec_tc, leigvec_tc )


end subroutine check

! ---


