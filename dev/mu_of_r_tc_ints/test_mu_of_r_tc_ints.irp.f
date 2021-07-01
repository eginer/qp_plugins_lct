program test_mu_of_r_tc_ints
  implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC
 my_grid_becke = .True.
 my_n_pt_r_grid = 30
 my_n_pt_a_grid = 170
 touch  my_grid_becke my_n_pt_r_grid my_n_pt_a_grid
 extra_grid_type_sgn = 1
 touch extra_grid_type_sgn
 my_extra_grid_becke = .False.
 touch my_extra_grid_becke
 print*,'Warning : the Becke grid parameters are automatically set to '
 print*,'my_n_pt_a_grid = ',my_n_pt_a_grid
 print*,'my_n_pt_r_grid = ',my_n_pt_r_grid
 print*,'If you want to modify them, you have to modify the following file '
 print*,'qp2/plugins/qp_plugins_lct/dev/transcorr_h/transcorr_general.irp.f'
 print*,'and recompile doing ninja'
!  constant_mu = .True.
!  touch constant_mu
!  call test_gauss_ij_rk
! call test_erf_mu_squared_ij_rk
! call test_big_array
! call test_big_array_3
! call test_num_scal_pot
! call test_num_deriv_pot
! call test_big_array_mo
! call test_big_array_mo_scal
! call test_big_array_mo_deriv
! call test_big_array_ao
! call test_num_scal_pot_mo
! call test_num_deriv_pot_mo
 call test_grad_jastrow
! call test_grad_jastrow_sq
! call test_nabla_1_sq_ao
! call test_non_hermit_ao
! call test_non_hermit_mo
! call test_lapl_j_ao
! call test_grad_gamma_r1
end
