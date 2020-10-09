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
 double precision :: delta_dm_a(n_points_final_grid), delta_dm_b(n_points_final_grid), dm_a_plus_delta(n_points_final_grid), dm_b_plus_delta(n_points_final_grid)
 double precision :: ex_lda_test, ex_lda_test_plus_delta, int_vx_lda_test
 double precision :: pi
 double precision :: delta_rho_11, integral_potential_lda, ex_lda_potential, weight
 integer :: i, istate
 pi = dacos(-1.d0)
 delta_rho_11 = 1.d-3
 integral_potential_lda = 0.d0
 istate = 1
 call delta_density_for_energy_test(delta_dm_a, delta_dm_b, delta_rho_11)
 
 do i=1, n_points_final_grid
   r(1) = final_grid_points(1,i) 
   r(2) = final_grid_points(2,i) 
   r(3) = final_grid_points(3,i)
   
   call dm_dft_alpha_beta_at_r(r, dm_a(i), dm_b(i))
   dm_a_plus_delta(i) = dm_a(i) + delta_dm_a(i)
   dm_b_plus_delta(i) = dm_b(i) + delta_dm_b(i)

   ex_lda_potential = -((3.d0/pi)**(1.d0/3.d0))*(dm_a(i) + dm_b(i))**(1.d0/3.d0)
   weight = final_weight_at_r_vector(i)
   integral_potential_lda += weight*(delta_dm_a(i) + delta_dm_b(i))*ex_lda_potential
 enddo 
   ! exLDA
   call energy_x_lda_test(dm_a, dm_b, ex_lda_test)
   call energy_x_lda_test(dm_a_plus_delta, dm_b_plus_delta, ex_lda_test_plus_delta)
   call int_potential_x_lda_test (delta_rho_11, int_vx_lda_test) 

   !exPBE
   call energy_x_pbe_test (dm_a, dm_b, grad_rho_a, grad_rho_b, ex_pbe_test)
   call energy_x_pbe_test (dm_a_plus_delta, dm_b_plus_delta, grad_rho_a, grad_rho_b, ex_pbe_test_plus_delta)
   call int_potential_x_pbe_test (delta_rho_11, int_vx_pbe_test)
 
   write(34,*) 'gamma_i_j = epsilon*delta_kronecker(i,1)*delta_kronecker(1,j)           =', delta_rho_11
   write(34,*) '----------------------------------------------LDA--------------------------------------------' 
   write(34,*) 'ex_lda_test[n]                                                          =', ex_lda_test
   write(34,*) 'ex_lda(provider)[n]                                                     =', energy_x_lda(1)
   write(34,*) 'ex_lda_test[n + delta_n]                                                =', ex_lda_test_plus_delta
   write(34,*) '------------------------------------------------'
   write(34,*) 'Theoretical value : int (dr delta_n * [delta (E) / delta(n(r))] )       =', integral_potential_lda
   write(34,*) 'We test the equality between the two following results :'
   write(34,*) '1st_method = ex_lda [n + delta_n] - ex_lda[n]                           =', ex_lda_test_plus_delta - ex_lda_test
   write(34,*) '2nd method = sum(i,j) delta_gamma_i_j * int (dr phi_i(r) v(r) phi_j(r)) =', int_vx_lda_test
   write(34,*) 'Relative error : (1st_method - 2nd_method)/1st_method                   =', (int_vx_lda_test - ex_lda_test_plus_delta + ex_lda_test)/(ex_lda_test_plus_delta - ex_lda_test)
   write(34,*) '------------------------------------------------'
   write(34,*) '------------------------------------------------'
   write(34,*) '----------------------------------------------PBE--------------------------------------------' 
end program

