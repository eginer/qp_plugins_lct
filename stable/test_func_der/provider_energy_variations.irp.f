 BEGIN_PROVIDER  [double precision, ex_lda_at_n]
 &BEGIN_PROVIDER [double precision, ex_lda_at_n_plus_delta_n]
 &BEGIN_PROVIDER [double precision, int_vx_lda_at_n]
   implicit none
   BEGIN_DOC
   ! 
   END_DOC
   double precision :: dm_a(n_points_final_grid),dm_b(n_points_final_grid)
   integer :: i
   do i=1, n_points_final_grid
    dm_a(i) = one_e_dm_and_grad_alpha_in_r(4,i,1)
    dm_b(i) = one_e_dm_and_grad_beta_in_r(4,i,1)
   enddo 
   call energy_x_lda_test(dm_a, dm_b, ex_lda_at_n)
   call energy_x_lda_test(dm_a_plus_delta_rho_a_in_r, dm_b_plus_delta_rho_b_in_r, ex_lda_at_n_plus_delta_n)
!  call int_potential_x_lda_test (delta_rho_11, dm_a, dm_b, int_vx_lda_at_n)  
  call compute_func_der(delta_gamma_i_j_alpha, delta_gamma_i_j_beta, potential_x_alpha_mo_lda, potential_x_beta_mo_lda, int_vx_lda_at_n)
 END_PROVIDER


 BEGIN_PROVIDER  [ double precision, ex_pbe_at_n ]
 &BEGIN_PROVIDER [ double precision, ex_pbe_at_n_plus_delta_n ]
 &BEGIN_PROVIDER [ double precision, int_vx_pbe_at_n ]
 &BEGIN_PROVIDER [ double precision, ec_pbe_at_n ]
 &BEGIN_PROVIDER [ double precision, ec_pbe_at_n_plus_delta_n ]
 &BEGIN_PROVIDER [ double precision, int_vc_pbe_at_n ]
   implicit none
   BEGIN_DOC
   ! 
   END_DOC
   double precision :: dm_a(n_points_final_grid),dm_b(n_points_final_grid)
   double precision :: grad_dm_a(3,n_points_final_grid), grad_dm_b(3,n_points_final_grid)
   integer :: i
   do i=1, n_points_final_grid
    grad_dm_a(1:3,i) = one_e_dm_and_grad_alpha_in_r(1:3,i,1)
    grad_dm_b(1:3,i) = one_e_dm_and_grad_beta_in_r(1:3,i,1)
    dm_a(i) = one_e_dm_and_grad_alpha_in_r(4,i,1)
    dm_b(i) = one_e_dm_and_grad_beta_in_r(4,i,1)
   enddo
   call energy_xc_pbe_test (dm_a, dm_b, grad_dm_a, grad_dm_b, ex_pbe_at_n, ec_pbe_at_n)
   call energy_xc_pbe_test (dm_a_plus_delta_rho_a_in_r, dm_b_plus_delta_rho_b_in_r, & 
                            grad_plus_delta_grad_rho_a_in_r,grad_plus_delta_grad_rho_b_in_r,  &
                            ex_pbe_at_n_plus_delta_n, ec_pbe_at_n_plus_delta_n)
!  call int_potential_xc_pbe_test (delta_rho_11, delta_rho_11, int_vx_pbe_at_n, int_vc_pbe_at_n)
   call compute_func_der(delta_gamma_i_j_alpha, delta_gamma_i_j_beta, potential_c_alpha_mo_pbe, potential_x_beta_mo_pbe, int_vx_pbe_at_n)
   call compute_func_der(delta_gamma_i_j_alpha, delta_gamma_i_j_beta, potential_c_alpha_mo_pbe, potential_c_beta_mo_pbe, int_vc_pbe_at_n)
 END_PROVIDER

 BEGIN_PROVIDER  [ double precision , ex_pbeUEG_at_n ]
 &BEGIN_PROVIDER [ double precision , ex_pbeUEG_at_n_plus_delta_n ]
 &BEGIN_PROVIDER [ double precision , int_vx_pbeUEG_at_n ]
 &BEGIN_PROVIDER [ double precision , ec_pbeUEG_at_n ]
 &BEGIN_PROVIDER [ double precision , ec_pbeUEG_at_n_plus_delta_n ]
 &BEGIN_PROVIDER [ double precision , int_vc_pbeUEG_at_n ]
   implicit none
   BEGIN_DOC
   !
   END_DOC
!  double precision :: r(3) 
!  double precision :: dm_a(n_points_final_grid), dm_b(n_points_final_grid)
!  double precision :: grad_dm_a(3,n_points_final_grid), grad_dm_b(3,n_points_final_grid)
!  double precision :: aos_array(ao_num,n_points_final_grid), grad_aos_array(3,ao_num,n_points_final_grid)
!  double precision :: delta_dm_a(n_points_final_grid), delta_dm_b(n_points_final_grid), delta_grad_dm_a(3,n_points_final_grid), delta_grad_dm_b(3,n_points_final_grid)
!  double precision :: dm_a_plus_delta(n_points_final_grid), dm_b_plus_delta(n_points_final_grid), grad_dm_a_plus_delta(3,n_points_final_grid), grad_dm_b_plus_delta(3,n_points_final_grid)
!  double precision :: delta_rho_11, delta_grad_rho_11
!  integer :: i, m
!  delta_rho_11 = 1.d-5
!  delta_grad_rho_11 = 1.d-5
!  call delta_density_for_energy_test_general(delta_dm_a, delta_dm_b, delta_rho_11, delta_rho_11)
!  call delta_grad_density_for_energy_test_general(delta_grad_dm_a, delta_grad_dm_b, delta_grad_rho_11, delta_grad_rho_11)
!  do i=1, n_points_final_grid
!    r(1) = final_grid_points(1,i) 
!    r(2) = final_grid_points(2,i) 
!    r(3) = final_grid_points(3,i)
!    call density_and_grad_alpha_beta_and_all_aos_and_grad_aos_at_r(r,dm_a(i),dm_b(i), grad_dm_a(1,i), grad_dm_b(1,i), aos_array(1,i), grad_aos_array(1,1,i))
!    dm_a_plus_delta(i) = dm_a(i) + delta_dm_a(i)
!    dm_b_plus_delta(i) = dm_b(i) + delta_dm_b(i)
!    do m=1,3
!     grad_dm_a_plus_delta(m,i) = grad_dm_a(m,i) + delta_grad_dm_a(m,i)
!     grad_dm_b_plus_delta(m,i) = grad_dm_b(m,i) + delta_grad_dm_b(m,i)
!    enddo
!  enddo 
!  call energy_xc_pbeUEG_test (dm_a, dm_b, grad_dm_a, grad_dm_b, ex_pbeUEG_at_n, ec_pbeUEG_at_n)
!  call energy_xc_pbeUEG_test (dm_a_plus_delta, dm_b_plus_delta, grad_dm_a_plus_delta, grad_dm_b_plus_delta, ex_pbeUEG_at_n_plus_delta_n, ec_pbeUEG_at_n_plus_delta_n)
!  call int_potential_xc_pbeUEG_test (delta_rho_11, int_vx_pbeUEG_at_n, int_vc_pbeUEG_at_n)
 END_PROVIDER

 BEGIN_PROVIDER  [ double precision , ec_pben2_at_n_n2 ]
 &BEGIN_PROVIDER [ double precision , ec_pben2_at_n_plus_delta_n_n2_plus_delta_n2 ]
 &BEGIN_PROVIDER [ double precision , ec_pben2_at_n_plus_delta_n ]
 &BEGIN_PROVIDER [ double precision , ec_pben2_at_n2_plus_delta_n2 ]
 &BEGIN_PROVIDER [ double precision , int_vc_pben2_one_e_at_n_n2 ]
 &BEGIN_PROVIDER [ double precision , int_vc_pben2_two_e_at_n_n2 ]
   implicit none
   BEGIN_DOC
   !
   END_DOC
   double precision :: dm_a(n_points_final_grid), dm_b(n_points_final_grid), n2(n_points_final_grid)
   double precision :: grad_dm_a(3,n_points_final_grid), grad_dm_b(3,n_points_final_grid)
   double precision :: ex_pben2_at_n_n2  
   double precision :: ex_pben2_at_n_plus_delta_n_n2_plus_delta_n2  
   double precision :: ex_pben2_at_n_plus_delta_n  
   double precision :: ex_pben2_at_n2_plus_delta_n2  
  integer :: i, m

  do i=1, n_points_final_grid
    n2(i) = total_cas_on_top_density(i,1)
    grad_dm_a(1:3,i) = one_e_dm_and_grad_alpha_in_r(1:3,i,1)
    grad_dm_b(1:3,i) = one_e_dm_and_grad_beta_in_r(1:3,i,1)
    dm_a(i) = one_e_dm_and_grad_alpha_in_r(4,i,1)
    dm_b(i) = one_e_dm_and_grad_beta_in_r(4,i,1)
  enddo 
  call energy_xc_pben2_test (dm_a ,dm_b ,grad_dm_a ,grad_dm_b ,n2 ,ec_pben2_at_n_n2 ,ex_pben2_at_n_n2)
  call energy_xc_pben2_test (dm_a_plus_delta_rho_a_in_r,dm_b_plus_delta_rho_b_in_r, & 
                             grad_plus_delta_grad_rho_a_in_r,grad_plus_delta_grad_rho_b_in_r, & 
                             n2_plus_dn2,& 
                             ec_pben2_at_n_plus_delta_n_n2_plus_delta_n2,ex_pben2_at_n_plus_delta_n_n2_plus_delta_n2)
  call energy_xc_pben2_test (dm_a_plus_delta_rho_a_in_r,dm_b_plus_delta_rho_b_in_r, & 
                             grad_plus_delta_grad_rho_a_in_r,grad_plus_delta_grad_rho_b_in_r, & 
                             n2,          &
                             ec_pben2_at_n_plus_delta_n   ,ex_pben2_at_n_plus_delta_n)
  call energy_xc_pben2_test (dm_a, dm_b, grad_dm_a, grad_dm_b, & 
                             n2_plus_dn2,& 
                             ec_pben2_at_n2_plus_delta_n2  ,ex_pben2_at_n2_plus_delta_n2)
  call compute_func_der(delta_gamma_i_j_alpha, delta_gamma_i_j_beta, & 
                        pot_basis_alpha_mo_su_pbe_ot, pot_basis_beta_mo_su_pbe_ot, & 
                        int_vc_pben2_one_e_at_n_n2)
  call compute_func_der_two_e(delta_n2_ijkl, eff_two_e, int_vc_pben2_two_e_at_n_n2)
 END_PROVIDER


