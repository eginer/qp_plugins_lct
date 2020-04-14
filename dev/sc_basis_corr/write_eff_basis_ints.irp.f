program write_effective_RSDFT_hamiltonian
 implicit none
 BEGIN_DOC
 ! This programs writes the effective RS-DFT Hamiltonian into the EZFIO folder. 
 ! The next programs that will run unto the EZFIO folder will, by default, 
 !
 ! have the one- and two-body integrals loaded from the EZFIO data. 
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

 ! regular 1/r12 integrals  on the MO basis
 io_mo_two_e_integrals = "None"
 touch io_mo_two_e_integrals
 ! regular 1/r12 integrals  on the AO basis
 io_ao_two_e_integrals = "None"
 touch io_ao_two_e_integrals
 ! integral of the effective potential 
 io_mo_int_mu_of_r = "None" 
 touch io_mo_int_mu_of_r
 call write_all_ints_basis
 call routines_compute_energy
end

subroutine routines_compute_energy
 implicit none
 BEGIN_DOC
! routine that computes the variational energy of the current wave function 
!
! and saves the current one-body density matrix into the aux_quantities folder the EZFIO folder
 END_DOC
 density_for_dft = "WFT" 
 touch density_for_dft 
 call print_variational_energy
 call ezfio_set_aux_quantities_data_one_e_dm_alpha_mo(one_e_dm_mo_alpha)
 call ezfio_set_aux_quantities_data_one_e_dm_beta_mo(one_e_dm_mo_beta)

end


