program compute_Hmu

  implicit none

  my_grid_becke = .True.
  my_n_pt_r_grid = 30
  my_n_pt_a_grid = 50
  touch  my_grid_becke my_n_pt_r_grid my_n_pt_a_grid

  read_wf = .True.
  touch read_wf

  PROVIDE N_int

  !call Hmu_psi()

  !call Hmu_IJ()

  call routine_3()

end

! ---

subroutine Hmu_psi()

  implicit none
  integer                       :: i
  double precision, allocatable :: delta(:) 

  print *, j1b_gauss
  print *, j1b_gauss_pen
  print *, mu_erf

  allocate( delta(N_det) )

  ! get < I | H_mu | psi > 
  call get_Htc_psi(psi_det, psi_coef, N_det, N_int, delta)

  ! order as QMCCHEM
  call dset_order(delta, psi_bilinear_matrix_order, N_det)

  do i = 1, N_det
    print *, delta(i)
  enddo 

  deallocate(delta)

  return
end subroutine Hmu_psi

! ---

subroutine Hmu_IJ()

  implicit none
  integer          :: i, j
  double precision :: hmono, heff, hderiv, hthree, htilde_ij

  do i = 1, N_det
    do j = 1, N_det
      call htilde_mu_mat(psi_det(1,1,i), psi_det(1,1,j), N_int, hmono, heff, hderiv, hthree, htilde_ij)
      write(11, *) htilde_ij
    enddo
  enddo

  return
end subroutine Hmu_IJ

! ---

subroutine routine_3()

  use bitmasks ! you need to include the bitmasks_module.f90 features

  implicit none
  integer                        :: i, a, i_ok, s1
  double precision               :: hmono, heff, hderiv, hthree, htilde_ij
  double precision               :: err_ai, err_tot
  integer(bit_kind), allocatable :: det_i(:,:)

  allocate(det_i(N_int,2))

  err_tot = 0.d0

  s1 = 1

  det_i = ref_bitmask
  call debug_det(det_i, N_int)
  print*, ' HF det'
  call debug_det(det_i, N_int)

  do i = 1, elec_alpha_num ! occupied
    do a = elec_alpha_num+1, mo_num ! virtual 


      det_i = ref_bitmask
      call do_single_excitation(det_i, i, a, s1, i_ok)
      if(i_ok == -1) then
       print*, 'PB !!'
       print*, i, a
       stop
      endif
      !print*, ' excited det'
      !call debug_det(det_i, N_int)

      call htilde_mu_mat(det_i, ref_bitmask, N_int, hmono, heff, hderiv, hthree, htilde_ij)
      err_ai = dabs(htilde_ij)
      if(err_ai .gt. 1d-7) then
        print*, ' warning on', i, a
        print*, hmono, heff, hthree, htilde_ij
      endif
      err_tot += err_ai

      write(11, *) htilde_ij
    enddo
  enddo

  print *, ' err_tot = ', err_tot

  deallocate(det_i)

end subroutine routine_3

! ---


