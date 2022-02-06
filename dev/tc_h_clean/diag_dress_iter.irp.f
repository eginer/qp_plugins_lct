program diag_dress_iter
  implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC
  my_grid_becke = .True.
  my_n_pt_r_grid = 30
  my_n_pt_a_grid = 50
  touch  my_grid_becke my_n_pt_r_grid my_n_pt_a_grid
  read_wf = .True.
  touch read_wf 
  call routine

end

subroutine routine_bis
 implicit none
 double precision, allocatable  :: delta(:) , psicoef(:),v(:),delta_av(:)
 allocate( delta(n_det) , psicoef(n_det), v(n_det) , delta_av(n_det))
 integer :: i
 do i = 1, N_det
  psicoef(i) = psi_coef(i,1)
 enddo
 call H_u_0_nstates_openmp(v,psicoef,1,N_det)
 call get_delta_tc_psi(psi_det,psicoef,n_det,delta)
 call get_delta_av_tc_psi(psi_det,psicoef,n_det,delta_av)
 double precision :: accu_h,accu_d,accu_av
 accu_h = 0.d0
 do i = 1, N_det
  accu_h += psicoef(i) * v(i)
 enddo
 print*,'accu_h = ',accu_h
 accu_d = 0.d0
 do i = 1, N_det
  accu_d += psicoef(i) * delta(i)
 enddo
 print*,'accu_d = ',accu_d
 accu_av = 0.d0
 do i = 1, N_det
  accu_av += psicoef(i) * delta_av(i)
 enddo
 print*,'accu_av= ',accu_av
 print*,'accutot= ',accu_h + accu_d
 print*,''

end

subroutine routine
 implicit none
 integer :: i,j
 print*,'eigval_right_tc = ',eigval_right_tc
 print*,'eigval_left_tc  = ',eigval_left_tc
 print*,'******************'
 print*,'< h_core >      = ',h_mono_comp_right_tc
 print*,'< h_eff_2e >    = ',h_eff_comp_right_tc
 print*,'< h_deriv_2_e > = ',h_deriv_comp_right_tc
 print*,'< h_eee >       = ',h_three_comp_right_tc
 print*,'< h_tot >       = ',h_tot_comp_right_tc
 print*,'******************'
 print*,'singles_hf_mat_ele',singles_hf_mat_elem
 print*,'e_tilde_00      = ',e_tilde_00
 print*,'E corr tc       = ',eigval_right_tc - e_tilde_00
 print*,'e_pt2_tc        = ',e_pt2_tc
 print*,'e_pt2_tc_single = ',e_pt2_tc_single
 print*,'e_pt2_tc_double = ',e_pt2_tc_double
 print*,'e_from_right_eigv ',e_from_right_eigv
 print*,'******************'
! print*,'******************'
! print*,'h_tilde_dagger_expect_right    = ',h_tilde_dagger_expect_right
! print*,'h_tilde_expect_right           = ',h_tilde_expect_right
! print*,'h_tilde_dagger_expect_psi_coef = ',h_tilde_dagger_expect_psi_coef
! print*,'h_tilde_expect_psi_coef        = ',h_tilde_dagger_expect_psi_coef
! print*,'eigval_sym_tc                  = ',eigval_sym_tc(1)
! print*,'exp_h_tilde_eigvec_sym_tc      = ',exp_h_tilde_eigvec_sym_tc
! print*,'******************'
 print*,'******************'
 if(comp_left_eigv.or.full_tc_h_solver)then
  print*,' norm_ground_left_right = ',norm_ground_left_right
  print*,' norm_ground_right       = ',norm_ground_right
  print*,' norm_ground_left       = ',norm_ground_left
 endif
 print*,'Left, right and usual eigenvectors '
 do i = 1, N_det
  double precision :: tmp
  tmp = (leigvec_tc(i,1)/leigvec_tc(1,1) + reigvec_tc(i,1)/reigvec_tc(1,1)) * 0.5d0
!  write(*,'(I5,X,(100(F12.7,X)))')i,leigvec_tc(i,1),reigvec_tc(i,1),psi_coef(i,1),eigvec_sym_tc(i,1)/eigvec_sym_tc(1,1),tmp
  write(*,'(I5,X,(100(F12.7,X)))')i,leigvec_tc(i,1),reigvec_tc(i,1),psi_coef(i,1)
 enddo
 call routine_save_right
end
