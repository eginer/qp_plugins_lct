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
 io_mo_one_e_integrals = "None"
 touch io_mo_one_e_integrals  
 io_mo_two_e_integrals = "None"
 touch io_mo_two_e_integrals
 io_ao_two_e_integrals = "None"
 touch io_ao_two_e_integrals
 io_mo_two_e_integrals_erf = "None" 
 touch io_mo_two_e_integrals_erf
 io_ao_two_e_integrals_erf = "None" 
 touch io_ao_two_e_integrals_erf
 call routines_write_int
 call routines_compute_energy
end

subroutine routines_write_int
 implicit none
 BEGIN_DOC
! routine that computes the effective RSDFT Hamiltoninan and writes it into the EZFIO folder
 END_DOC
 call write_all_integrals_for_mrdft
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
 call print_variational_energy_dft
 call ezfio_set_aux_quantities_data_one_body_alpha_dm_mo(one_body_dm_mo_alpha)
 call ezfio_set_aux_quantities_data_one_body_beta_dm_mo(one_body_dm_mo_beta)

end


