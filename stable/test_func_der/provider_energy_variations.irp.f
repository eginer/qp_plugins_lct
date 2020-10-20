 BEGIN_PROVIDER  [double precision, ex_lda_at_n]
 &BEGIN_PROVIDER [double precision, ex_lda_at_n_plus_delta_n]
 &BEGIN_PROVIDER [double precision, int_vx_lda_at_n]
   implicit none
   BEGIN_DOC
   ! 
   END_DOC
   double precision :: delta_rho_11, r(3)
   double precision :: dm_a(n_points_final_grid),dm_b(n_points_final_grid)
   double precision :: dm_a_plus_delta(n_points_final_grid),dm_b_plus_delta(n_points_final_grid)
   double precision :: delta_dm_a(n_points_final_grid), delta_dm_b(n_points_final_grid) 
   integer :: i
   delta_rho_11 = 1.d-5
   call delta_density_for_energy_test(delta_dm_a, delta_dm_b, delta_rho_11, delta_rho_11)
   do i=1, n_points_final_grid
    r(1) = final_grid_points(1,i) 
    r(2) = final_grid_points(2,i) 
    r(3) = final_grid_points(3,i)
    call dm_dft_alpha_beta_at_r(r,dm_a(i),dm_b(i))
    dm_a_plus_delta(i) = dm_a(i) + delta_dm_a(i)
    dm_b_plus_delta(i) = dm_b(i) + delta_dm_b(i)
   enddo 
   call energy_x_lda_test(dm_a, dm_b, ex_lda_at_n)
   call energy_x_lda_test(dm_a_plus_delta, dm_b_plus_delta, ex_lda_at_n_plus_delta_n)
   call int_potential_x_lda_test (delta_rho_11, dm_a, dm_b, int_vx_lda_at_n)  
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
   double precision :: delta_rho_11, delta_grad_rho_11, r(3)
   double precision :: dm_a(n_points_final_grid),dm_b(n_points_final_grid)
   double precision :: grad_dm_a(3,n_points_final_grid), grad_dm_b(3,n_points_final_grid)
   double precision :: dm_a_plus_delta(n_points_final_grid),dm_b_plus_delta(n_points_final_grid)
   double precision :: grad_dm_a_plus_delta(3,n_points_final_grid), grad_dm_b_plus_delta(3,n_points_final_grid)
   double precision :: delta_dm_a(n_points_final_grid), delta_dm_b(n_points_final_grid) 
   double precision :: delta_grad_dm_a(3,n_points_final_grid), delta_grad_dm_b(3,n_points_final_grid) 
   double precision :: aos_array(ao_num,n_points_final_grid), grad_aos_array(3,ao_num,n_points_final_grid)
   integer :: i,m
   delta_rho_11 = 1.d-5
   delta_grad_rho_11 = 1.d-5
   call delta_density_for_energy_test(delta_dm_a, delta_dm_b, delta_rho_11, delta_rho_11)
   call delta_grad_density_for_energy_test(delta_grad_dm_a, delta_grad_dm_b, delta_grad_rho_11, delta_grad_rho_11)
   do i=1, n_points_final_grid
    r(1) = final_grid_points(1,i) 
    r(2) = final_grid_points(2,i) 
    r(3) = final_grid_points(3,i)
    call density_and_grad_alpha_beta_and_all_aos_and_grad_aos_at_r(r,dm_a(i),dm_b(i), grad_dm_a(1,i), grad_dm_b(1,i), aos_array(1,i), grad_aos_array(1,1,i))
    dm_a_plus_delta(i) = dm_a(i) + delta_dm_a(i)
    dm_b_plus_delta(i) = dm_b(i) + delta_dm_b(i)
    do m=1,3
     grad_dm_a_plus_delta(m,i) = grad_dm_a(m,i) + delta_grad_dm_a(m,i)
     grad_dm_b_plus_delta(m,i) = grad_dm_b(m,i) + delta_grad_dm_b(m,i)
    enddo
   enddo
   call energy_xc_pbe_test (dm_a, dm_b, grad_dm_a, grad_dm_b, ex_pbe_at_n, ec_pbe_at_n)
   call energy_xc_pbe_test (dm_a_plus_delta, dm_b_plus_delta, grad_dm_a_plus_delta, grad_dm_b_plus_delta, ex_pbe_at_n_plus_delta_n, ec_pbe_at_n_plus_delta_n)
   call int_potential_xc_pbe_test (delta_rho_11, int_vx_pbe_at_n, int_vc_pbe_at_n)
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
   double precision :: r(3) 
   double precision :: dm_a(n_points_final_grid), dm_b(n_points_final_grid)
   double precision :: grad_dm_a(3,n_points_final_grid), grad_dm_b(3,n_points_final_grid)
   double precision :: aos_array(ao_num,n_points_final_grid), grad_aos_array(3,ao_num,n_points_final_grid)
   double precision :: delta_dm_a(n_points_final_grid), delta_dm_b(n_points_final_grid), delta_grad_dm_a(3,n_points_final_grid), delta_grad_dm_b(3,n_points_final_grid)
   double precision :: dm_a_plus_delta(n_points_final_grid), dm_b_plus_delta(n_points_final_grid), grad_dm_a_plus_delta(3,n_points_final_grid), grad_dm_b_plus_delta(3,n_points_final_grid)
   double precision :: delta_rho_11, delta_grad_rho_11
   integer :: i, m
   delta_rho_11 = 1.d-5
   delta_grad_rho_11 = 1.d-5
   call delta_density_for_energy_test(delta_dm_a, delta_dm_b, delta_rho_11, delta_rho_11)
   call delta_grad_density_for_energy_test(delta_grad_dm_a, delta_grad_dm_b, delta_grad_rho_11, delta_grad_rho_11)
   do i=1, n_points_final_grid
     r(1) = final_grid_points(1,i) 
     r(2) = final_grid_points(2,i) 
     r(3) = final_grid_points(3,i)
     call density_and_grad_alpha_beta_and_all_aos_and_grad_aos_at_r(r,dm_a(i),dm_b(i), grad_dm_a(1,i), grad_dm_b(1,i), aos_array(1,i), grad_aos_array(1,1,i))
     dm_a_plus_delta(i) = dm_a(i) + delta_dm_a(i)
     dm_b_plus_delta(i) = dm_b(i) + delta_dm_b(i)
     do m=1,3
      grad_dm_a_plus_delta(m,i) = grad_dm_a(m,i) + delta_grad_dm_a(m,i)
      grad_dm_b_plus_delta(m,i) = grad_dm_b(m,i) + delta_grad_dm_b(m,i)
     enddo
   enddo 
   call energy_xc_pbeUEG_test (dm_a, dm_b, grad_dm_a, grad_dm_b, ex_pbeUEG_at_n, ec_pbeUEG_at_n)
   call energy_xc_pbeUEG_test (dm_a_plus_delta, dm_b_plus_delta, grad_dm_a_plus_delta, grad_dm_b_plus_delta, ex_pbeUEG_at_n_plus_delta_n, ec_pbeUEG_at_n_plus_delta_n)
   call int_potential_xc_pbeUEG_test (delta_rho_11, int_vx_pbeUEG_at_n, int_vc_pbeUEG_at_n)
 END_PROVIDER

 BEGIN_PROVIDER  [ double precision , ec_pben2_at_n_n2 ]
 &BEGIN_PROVIDER [ double precision , ec_pben2_at_n_plus_delta_n_n2_plus_delta_n2 ]
 &BEGIN_PROVIDER [ double precision , ec_pben2_at_n_plus_delta_n ]
 &BEGIN_PROVIDER [ double precision , ec_pben2_at_n2_plus_delta_n2 ]
 &BEGIN_PROVIDER [ double precision , ex_pben2_at_n_n2 ]
 &BEGIN_PROVIDER [ double precision , ex_pben2_at_n_plus_delta_n_n2_plus_delta_n2 ]
 &BEGIN_PROVIDER [ double precision , ex_pben2_at_n_plus_delta_n ]
 &BEGIN_PROVIDER [ double precision , ex_pben2_at_n2_plus_delta_n2 ]
 &BEGIN_PROVIDER [ double precision , int_vc_pben2_one_e_at_n_n2 ]
 &BEGIN_PROVIDER [ double precision , int_vc_pben2_two_e_at_n_n2 ]
   implicit none
   BEGIN_DOC
   !
   END_DOC
   double precision :: r(3) 
   double precision :: dm_a(n_points_final_grid), dm_b(n_points_final_grid), n2(n_points_final_grid)
   double precision :: grad_dm_a(3,n_points_final_grid), grad_dm_b(3,n_points_final_grid)
   double precision :: aos_array(ao_num,n_points_final_grid), grad_aos_array(3,ao_num,n_points_final_grid)
   double precision :: delta_dm_a(n_points_final_grid), delta_dm_b(n_points_final_grid), delta_grad_dm_a(3,n_points_final_grid), delta_grad_dm_b(3,n_points_final_grid), delta_n2(n_points_final_grid)
   double precision :: dm_a_plus_delta(n_points_final_grid), dm_b_plus_delta(n_points_final_grid), grad_dm_a_plus_delta(3,n_points_final_grid), grad_dm_b_plus_delta(3,n_points_final_grid), n2_plus_delta(n_points_final_grid)
   double precision :: delta_rho_11, delta_rho_11_alpha, delta_rho_11_beta, delta_grad_rho_11, delta_n2_11
  integer :: i, m
  delta_rho_11 = 1.d-5
  delta_rho_11_alpha = 1.d-5
  delta_rho_11_beta = 1.d-5
  delta_n2_11 = 1.d-5
  delta_grad_rho_11 = 1.d-5
  call delta_density_for_energy_test(delta_dm_a, delta_dm_b, delta_rho_11, delta_rho_11)
  call delta_grad_density_for_energy_test(delta_grad_dm_a, delta_grad_dm_b, delta_grad_rho_11, delta_grad_rho_11)
  call delta_n2_for_energy_test (delta_n2, delta_n2_11)
  do i=1, n_points_final_grid
    n2(i) = total_cas_on_top_density(i,1)
    r(1) = final_grid_points(1,i) 
    r(2) = final_grid_points(2,i) 
    r(3) = final_grid_points(3,i)
    call density_and_grad_alpha_beta_and_all_aos_and_grad_aos_at_r(r,dm_a(i),dm_b(i), grad_dm_a(1,i), grad_dm_b(1,i), aos_array(1,i), grad_aos_array(1,1,i))
    dm_a_plus_delta(i) = dm_a(i) + delta_dm_a(i)
    dm_b_plus_delta(i) = dm_b(i) + delta_dm_b(i)
    do m=1,3
     grad_dm_a_plus_delta(m,i) = grad_dm_a(m,i) + delta_grad_dm_a(m,i)
     grad_dm_b_plus_delta(m,i) = grad_dm_b(m,i) + delta_grad_dm_b(m,i)
    enddo
    n2_plus_delta(i) = n2(i) + delta_n2(i)
  enddo 
  call energy_xc_pben2_test (dm_a           ,dm_b           ,grad_dm_a           ,grad_dm_b           ,n2           ,ec_pben2_at_n_n2                ,ex_pben2_at_n_n2)
  call energy_xc_pben2_test (dm_a_plus_delta,dm_b_plus_delta,grad_dm_a_plus_delta,grad_dm_b_plus_delta,n2_plus_delta,ec_pben2_at_n_plus_delta_n_n2_plus_delta_n2,ex_pben2_at_n_plus_delta_n_n2_plus_delta_n2)
  call energy_xc_pben2_test (dm_a_plus_delta,dm_b_plus_delta,grad_dm_a_plus_delta,grad_dm_b_plus_delta,n2           ,ec_pben2_at_n_plus_delta_n   ,ex_pben2_at_n_plus_delta_n)
  call energy_xc_pben2_test (dm_a           ,dm_b           ,grad_dm_a           ,grad_dm_b           ,n2_plus_delta,ec_pben2_at_n2_plus_delta_n2  ,ex_pben2_at_n2_plus_delta_n2)
  call int_potential_c_pben2_two_e_test (delta_n2_11, int_vc_pben2_two_e_at_n_n2)
  call int_potential_c_pben2_one_e_test (delta_rho_11_alpha, delta_rho_11_beta, int_vc_pben2_one_e_at_n_n2)
 END_PROVIDER


