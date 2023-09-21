program test_d2ecmdLDAUEG  
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

 integer :: istate, ipoint, nx

 double precision :: mu
 double precision :: rho_a,rho_b,ec_srmuLDAn,decdrho_a,decdrho_b,d2ecdrho_a2,d2ecdrho_b2
 
 do istate = 1, N_states
  do ipoint = 1, n_points_final_grid

  rho_a =  one_e_dm_and_grad_alpha_in_r(4,ipoint,istate)
  rho_b =  one_e_dm_and_grad_beta_in_r(4,ipoint,istate)

  mu = mu_of_r_prov(ipoint,istate)
  call ecmdsrLDAn(mu,rho_a,rho_b,ec_srmuLDAn,decdrho_a,decdrho_b,d2ecdrho_a2,d2ecdrho_b2)

  enddo
 enddo
end program


