program print_emu

  implicit none

  my_grid_becke = .True.
  my_n_pt_r_grid = 30
  my_n_pt_a_grid = 50
  touch  my_grid_becke my_n_pt_r_grid my_n_pt_a_grid

  read_wf = .True.
  touch read_wf

  PROVIDE N_det N_int
  call main()

end


subroutine main()

  use bitmasks

  implicit none
  integer                       :: i, j
  double precision              :: htilde_ij, hmono, heff, hderiv, hthree
  double precision, allocatable :: e_tmp(:), tmp, e_tot
  double precision, allocatable :: e_tmp_3b(:), tmp_3b, e_3b

  print *, j1b_gauss
  print *, j1b_gauss_pen
  print *, mu_erf

  i = 1
  j = 1
  call htilde_mu_mat(psi_det(1,1,i), psi_det(1,1,j), N_int, hmono, heff, hderiv, hthree, htilde_ij)
  print *, hthree
  print *, htilde_ij - hthree + nuclear_repulsion
  print *, htilde_ij + nuclear_repulsion

  allocate( e_tmp(N_det), e_tmp_3b(N_det) )
  e_tmp   (1:N_det) = 0.d0
  e_tmp_3b(1:N_det) = 0.d0

 !$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(dynamic,8)            &
 !$OMP SHARED(N_det, psi_det, psi_coef, N_int, e_tmp, e_tmp_3b) &
 !$OMP PRIVATE(i, j, tmp, tmp_3b, hmono, heff, hderiv, hthree, htilde_ij)
  do i = 1, N_det
    do j = 1, N_det
      call htilde_mu_mat(psi_det(1,1,i), psi_det(1,1,j), N_int, hmono, heff, hderiv, hthree, htilde_ij)
      tmp          = psi_coef(j,1) * htilde_ij
      tmp_3b       = psi_coef(j,1) * hthree
      e_tmp(i)    += tmp
      e_tmp_3b(i) += tmp_3b
    enddo
  enddo
 !$OMP END PARALLEL DO

  e_3b  = 0.d0
  e_tot = 0.d0
  do i = 1, N_det
    e_3b  += psi_coef(i,1) * e_tmp_3b(i)
    e_tot += psi_coef(i,1) * e_tmp(i)
  enddo

  print *, "e_3b   = ", e_3b
  print *, "e no3b = ", e_tot - e_3b + nuclear_repulsion
  print *, "e_tot  = ", e_tot + nuclear_repulsion

  deallocate( e_tmp, e_tmp_3b )

  return
end subroutine main()




