program test_energy_integral_lda
 implicit none
 BEGIN_DOC
 !
 END_DOC

 read_wf = .true.
 touch read_wf
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
 double precision :: dm_a(n_points_final_grid), dm_b(n_points_final_grid)
 double precision :: grad_dm_a(3,n_points_final_grid), grad_dm_b(3,n_points_final_grid)
 double precision :: aos_array(ao_num,n_points_final_grid), grad_aos_array(3,ao_num,n_points_final_grid)
 double precision :: delta_dm_a(n_points_final_grid), delta_dm_b(n_points_final_grid), delta_grad_dm_a(3,n_points_final_grid), delta_grad_dm_b(3,n_points_final_grid)
 double precision :: dm_a_plus_delta(n_points_final_grid), dm_b_plus_delta(n_points_final_grid), grad_dm_a_plus_delta(3,n_points_final_grid), grad_dm_b_plus_delta(3,n_points_final_grid)
 double precision :: ex_lda_test, ex_lda_test_plus_delta, int_vx_lda_test
 double precision :: ex_pbe_test, ex_pbe_test_plus_delta, int_vx_pbe_test
 double precision :: ec_pbe_test, ec_pbe_test_plus_delta, int_vc_pbe_test
 double precision :: pi
 double precision :: delta_rho_11, delta_grad_rho_11, integral_potential_lda, ex_lda_potential, weight
 integer :: i, m, istate
 pi = dacos(-1.d0)
 delta_rho_11 = 1.d-5
 delta_grad_rho_11 = 1.d-5
 integral_potential_lda = 0.d0
 istate = 1
 call delta_density_for_energy_test(delta_dm_a, delta_dm_b, delta_rho_11)
 call delta_grad_density_for_energy_test(delta_grad_dm_a, delta_grad_dm_b, delta_grad_rho_11)
 
 do i=1, n_points_final_grid
   r(1) = final_grid_points(1,i) 
   r(2) = final_grid_points(2,i) 
   r(3) = final_grid_points(3,i)
   
  ! call dm_dft_alpha_beta_at_r(r, dm_a(i), dm_b(i))
   call density_and_grad_alpha_beta_and_all_aos_and_grad_aos_at_r(r,dm_a(i),dm_b(i), grad_dm_a(1,i), grad_dm_b(1,i), aos_array(1,i), grad_aos_array(1,1,i))
   dm_a_plus_delta(i) = dm_a(i) + delta_dm_a(i)
   dm_b_plus_delta(i) = dm_b(i) + delta_dm_b(i)
  
   do m=1,3
    grad_dm_a_plus_delta(m,i) = grad_dm_a(m,i) + delta_grad_dm_a(m,i)
    grad_dm_b_plus_delta(m,i) = grad_dm_b(m,i) + delta_grad_dm_b(m,i)
   enddo
   
 !  weight = final_weight_at_r_vector(i)
   
 !  ex_lda_potential = -((3.d0/pi)**(1.d0/3.d0))*(dm_a(i) + dm_b(i))**(1.d0/3.d0)
 !  integral_potential_lda += weight*(delta_dm_a(i) + delta_dm_b(i))*ex_lda_potential
 enddo 
   
   ! exLDA
   call energy_x_lda_test(dm_a, dm_b, ex_lda_test)
   call energy_x_lda_test(dm_a_plus_delta, dm_b_plus_delta, ex_lda_test_plus_delta)
   call int_potential_x_lda_test (delta_rho_11, dm_a, dm_b, int_vx_lda_test) 

   !excPBE
   call energy_xc_pbe_test (dm_a, dm_b, grad_dm_a, grad_dm_b, ex_pbe_test, ec_pbe_test)
   call energy_xc_pbe_test (dm_a_plus_delta, dm_b_plus_delta, grad_dm_a_plus_delta, grad_dm_b_plus_delta, ex_pbe_test_plus_delta, ec_pbe_test_plus_delta)
   call int_potential_xc_pbe_test (delta_rho_11, int_vx_pbe_test, int_vc_pbe_test)
   print*,'ex_lda_test            = ',ex_lda_test
   print*,'ex_lda_test_plus_delta = ',ex_lda_test_plus_delta
   print*,'delta e                = ',ex_lda_test_plus_delta - ex_lda_test
   print*,'new method             = ',int_vx_lda_test
!   print*,'delta_gamm_pot_x_lda   = ',delta_gamm_pot_x_lda
   print*,''
   print*,'ex_pbe_test            = ',ex_pbe_test
   print*,'ex_pbe_test_plus_delta = ',ex_pbe_test_plus_delta
   print*,'delta e                = ',ex_pbe_test_plus_delta - ex_pbe_test
   print*,'new method             = ',int_vx_pbe_test
!   print*,'delta_gamm_pot_x_pbe   = ',delta_gamm_pot_x_pbe
  ! call energy_xc_pbe_test (dm_a_plus_delta, dm_b_plus_delta, grad_dm_a, grad_dm_b, ex_pbe_test_plus_delta, ec_pbe_test_plus_delta)
   
   call write_energies_test_on_file(delta_rho_11, ex_lda_test, ex_lda_test_plus_delta,int_vx_lda_test, ex_pbe_test, ex_pbe_test_plus_delta,int_vx_pbe_test, ec_pbe_test, ec_pbe_test_plus_delta,int_vc_pbe_test)
 
end program

