program test_energy_integral_lda
 implicit none
 BEGIN_DOC
 !
 END_DOC

 read_wf = .true.
 touch read_wf
 basis_cor_func = "su_pbe_ot"
 touch basis_cor_func

 ! total one-e integrals 
 io_mo_one_e_integrals = "None"
 touch io_mo_one_e_integrals  
 ! Vne integrals on the MO basis 
 io_mo_integrals_n_e = "None"
 touch io_mo_integrals_n_e
 ! kinetic integrals on the MO basis 
 io_mo_integrals_kinetic = "None"
 touch io_mo_integrals_kinetic 
 ! Vne integrals on the AO basis 
 io_ao_integrals_n_e = "None"
 touch io_ao_integrals_n_e 
 ! kinetic integrals on the AO basis 
 io_ao_integrals_kinetic = "None"
 touch io_ao_integrals_kinetic 

 double precision :: r(3) 
 double precision :: dm_a(n_points_final_grid), dm_b(n_points_final_grid), n2(n_points_final_grid)
 double precision :: grad_dm_a(3,n_points_final_grid), grad_dm_b(3,n_points_final_grid)
 double precision :: aos_array(ao_num,n_points_final_grid), grad_aos_array(3,ao_num,n_points_final_grid)
 double precision :: delta_dm_a(n_points_final_grid), delta_dm_b(n_points_final_grid), delta_grad_dm_a(3,n_points_final_grid), delta_grad_dm_b(3,n_points_final_grid), delta_n2(n_points_final_grid)
 double precision :: dm_a_plus_delta(n_points_final_grid), dm_b_plus_delta(n_points_final_grid), grad_dm_a_plus_delta(3,n_points_final_grid), grad_dm_b_plus_delta(3,n_points_final_grid), n2_plus_delta(n_points_final_grid)
 double precision :: ex_lda_test, ex_lda_test_plus_delta, int_vx_lda_test
 double precision :: ex_pbe_test, ex_pbe_test_plus_delta, int_vx_pbe_test
 double precision :: ec_pbe_test, ec_pbe_test_plus_delta, int_vc_pbe_test
 double precision :: ex_pbeUEG_test, ex_pbeUEG_test_plus_delta, int_vx_pbeUEG_test
 double precision :: ec_pbeUEG_test, ec_pbeUEG_test_plus_delta, int_vc_pbeUEG_test
 double precision :: ec_pben2_test, ec_pben2_test_plus_delta_n_n2 ,ec_pben2_test_plus_delta_n ,ec_pben2_test_plus_delta_n2, int_vc_pben2_one_e_test, int_vc_pben2_two_e_test
 double precision :: ex_pben2_test, ex_pben2_test_plus_delta_n_n2 ,ex_pben2_test_plus_delta_n ,ex_pben2_test_plus_delta_n2, int_vx_pben2_one_e_test, int_vx_pben2_two_e_test

 double precision :: pi
 double precision :: delta_rho_11, delta_grad_rho_11, delta_n2_11, integral_potential_lda, ex_lda_potential, weight
 integer :: i, m, istate
 pi = dacos(-1.d0)
 delta_rho_11 = 1.d-5
 delta_n2_11 = 1.d-5
 delta_grad_rho_11 = 1.d-5
 integral_potential_lda = 0.d0
 istate = 1
 call delta_density_for_energy_test(delta_dm_a, delta_dm_b, delta_rho_11)
 call delta_grad_density_for_energy_test(delta_grad_dm_a, delta_grad_dm_b, delta_grad_rho_11)
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
   
   ! exLDA
   call energy_x_lda_test(dm_a, dm_b, ex_lda_test)
   call energy_x_lda_test(dm_a_plus_delta, dm_b_plus_delta, ex_lda_test_plus_delta)
   call int_potential_x_lda_test (delta_rho_11, dm_a, dm_b, int_vx_lda_test) 

   !excPBE
   call energy_xc_pbe_test (dm_a, dm_b, grad_dm_a, grad_dm_b, ex_pbe_test, ec_pbe_test)
   call energy_xc_pbe_test (dm_a_plus_delta, dm_b_plus_delta, grad_dm_a_plus_delta, grad_dm_b_plus_delta, ex_pbe_test_plus_delta, ec_pbe_test_plus_delta)
   call int_potential_xc_pbe_test (delta_rho_11, int_vx_pbe_test, int_vc_pbe_test)
  
   !excPBEUEG
   call energy_xc_pbeUEG_test (dm_a, dm_b, grad_dm_a, grad_dm_b, ex_pbeUEG_test, ec_pbeUEG_test)
   call energy_xc_pbeUEG_test (dm_a_plus_delta, dm_b_plus_delta, grad_dm_a_plus_delta, grad_dm_b_plus_delta, ex_pbeUEG_test_plus_delta, ec_pbeUEG_test_plus_delta)
   call int_potential_xc_pbeUEG_test (delta_rho_11, int_vx_pbeUEG_test, int_vc_pbeUEG_test)

   !excPBEn2
   call energy_xc_pben2_test (dm_a           ,dm_b           ,grad_dm_a           ,grad_dm_b           ,n2           ,ec_pben2_test                ,ex_pben2_test)
   call energy_xc_pben2_test (dm_a_plus_delta,dm_b_plus_delta,grad_dm_a_plus_delta,grad_dm_b_plus_delta,n2_plus_delta,ec_pben2_test_plus_delta_n_n2,ex_pben2_test_plus_delta_n_n2)
   call energy_xc_pben2_test (dm_a_plus_delta,dm_b_plus_delta,grad_dm_a_plus_delta,grad_dm_b_plus_delta,n2           ,ec_pben2_test_plus_delta_n   ,ex_pben2_test_plus_delta_n)
   call energy_xc_pben2_test (dm_a           ,dm_b           ,grad_dm_a           ,grad_dm_b           ,n2_plus_delta,ec_pben2_test_plus_delta_n2  ,ex_pben2_test_plus_delta_n2)
   call int_potential_c_pben2_two_e_test (delta_n2_11, int_vc_pben2_two_e_test)
   call int_potential_c_pben2_one_e_test (delta_rho_11, int_vc_pben2_one_e_test)
   
   call write_energies_test_on_file(delta_rho_11, ex_lda_test, ex_lda_test_plus_delta,int_vx_lda_test, ex_pbe_test, ex_pbe_test_plus_delta,int_vx_pbe_test, ec_pbe_test, ec_pbe_test_plus_delta,int_vc_pbe_test,ex_pbeUEG_test, ex_pbeUEG_test_plus_delta,int_vx_pbeUEG_test, ec_pbeUEG_test, ec_pbeUEG_test_plus_delta,int_vc_pbeUEG_test,ec_pben2_test, ec_pben2_test_plus_delta_n_n2, ec_pben2_test_plus_delta_n, ec_pben2_test_plus_delta_n2, int_vc_pben2_one_e_test, int_vc_pben2_two_e_test)
 
end program

