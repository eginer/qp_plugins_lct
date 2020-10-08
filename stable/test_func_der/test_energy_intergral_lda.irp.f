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
 double precision :: ex_lda_test, ex_lda_test_plus_delta

 double precision :: mu_in
 double precision :: grad_rho_a_2(n_points_final_grid),grad_rho_b_2(n_points_final_grid),grad_rho_a_b(n_points_final_grid), rho2(n_points_final_grid)
 double precision :: ec_srmuPBE,decdrho_a(n_points_final_grid),decdrho_b(n_points_final_grid),decdgrad_rho_a_2(n_points_final_grid),decdgrad_rho_b_2(n_points_final_grid),decdgrad_rho_a_b(n_points_final_grid), decdgrad_rho_2(n_points_final_grid), decdrho(n_points_final_grid), decdrho2(n_points_final_grid), aos_array(ao_num), grad_aos_array(3,ao_num)!, decdrho2_a, decdrho2_b


 double precision :: delta_rho_11, integral_potential_lda, ex_lda_potential, weight
 double precision :: pi
 integer :: i, istate
 
 pi = dacos(-1.d0)
 delta_rho_11 = 1.d-3
 integral_potential_lda = 0.d0
 istate = 1
 call delta_density_for_energy_test(delta_dm_a, delta_dm_b, delta_rho_11)
 
 do i=1, n_points_final_grid
   mu_in = mu_of_r_prov(i,istate)
   r(1) = final_grid_points(1,i) 
   r(2) = final_grid_points(2,i) 
   r(3) = final_grid_points(3,i)
   
 !  call dm_dft_alpha_beta_at_r(r, dm_a(i), dm_b(i))
   call density_and_grad_alpha_beta_and_all_aos_and_grad_aos_at_r(r,dm_a(i),dm_b(i), grad_rho_a(i), grad_rho_b(i), aos_array, grad_aos_array)
   dm_a_plus_delta(i) = dm_a(i) + delta_dm_a(i)
   dm_b_plus_delta(i) = dm_b(i) + delta_dm_b(i)
!   grad_rho_a(1:3) = one_e_dm_and_grad_alpha_in_r(1:3,ipoint,istate)
!   grad_rho_b(1:3) = one_e_dm_and_grad_beta_in_r(1:3,ipoint,istate)
 !  grad_rho_a_2(i)    =
 !  grad_rho_b_2(i)    =
 !  grad_rho_a_b(i)    =
   rho2(i)            = on_top_cas_mu_r(i,istate)   

   ex_lda_potential = -((3.d0/pi)**(1.d0/3.d0))*(dm_a(i) + dm_b(i))**(1.d0/3.d0)
   weight = final_weight_at_r_vector(i)
   integral_potential_lda += weight*(delta_dm_a(i) + delta_dm_b(i))*ex_lda_potential
   write(35,*) 'delta_dm_a(', i, ') =', delta_dm_a(i)
   write(35,*) 'delta_dm_b(', i, ') =', delta_dm_b(i)
 enddo 
   ! exLDA
   call energy_x_lda_test(dm_a, dm_b, ex_lda_test)
   call energy_x_lda_test(dm_a_plus_delta, dm_b_plus_delta, ex_lda_test_plus_delta)
   
   !ecmdsrPBEn2
   call ecmdsrPBEn2_test(mu_in,dm_a,dm_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b,rho2,ec_srmuPBE,decdrho_a,decdrho_b, decdrho, decdgrad_rho_a_2,decdgrad_rho_b_2,decdgrad_rho_a_b, decdgrad_rho_2,decdrho2)

   write(34,*) 'ex_lda_test[n]                                =', ex_lda_test
   write(34,*) 'ex_lda(provider)[n]                           =', energy_x_lda(1)
   write(34,*) 'ex_lda_test[n + delta_n]                      =', ex_lda_test_plus_delta
   write(34,*) '------------------------------------------------' 
   write(34,*) 'We test the equality between the two following results (knowing the exact formulation for the LDA exchange potential) :'
   write(34,*) 'ex_lda [n + delta_n] - ex_lda[n]              =', ex_lda_test_plus_delta - ex_lda_test
   write(34,*) 'int (dr delta_n * [delta (E) / delta(n(r))] ) =', integral_potential_lda
   write(34,*) '  '
   write(34,*) '  '
   write(34,*) '  '
   write(34,*) 'Test on ECMD potentials :'

end program

