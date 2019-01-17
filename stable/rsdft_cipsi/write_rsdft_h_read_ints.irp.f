program write_rsdft_h_read_ints
 implicit none
 BEGIN_DOC
 ! This programs writes the effective RS-DFT Hamiltonian into the EZFIO folder. 
 ! !!! BUT !!! the integrals will be read from input, it must be used only very carefully 
 ! 
 END_DOC
 read_wf = .true.
 touch read_wf
 io_mo_one_e_integrals = "None"
 touch io_mo_one_e_integrals  

 io_mo_two_e_integrals = "Read"
 touch io_mo_two_e_integrals
 io_ao_two_e_integrals = "Read"
 touch io_ao_two_e_integrals
 io_mo_two_e_integrals_erf = "Read"
 touch io_mo_two_e_integrals_erf
 io_ao_two_e_integrals_erf = "Read"
 touch io_ao_two_e_integrals_erf


 io_mo_integrals_e_n = "None"
 touch io_mo_integrals_e_n
 io_mo_integrals_kinetic = "None"
 touch io_mo_integrals_kinetic 
 io_ao_integrals_e_n = "None"
 touch io_ao_integrals_e_n 
 io_ao_integrals_kinetic = "None"
 touch io_ao_integrals_kinetic 

 call routines_write_rsdft_h_read_ints
 call routines_compute_energy_h_read_ints
end

subroutine routines_write_rsdft_h_read_ints
 implicit none
 BEGIN_DOC
! routine that computes the effective RSDFT Hamiltoninan and writes it into the EZFIO folder
 END_DOC
 call write_all_integrals_for_mrdft_read_all
end 

subroutine routines_compute_energy_h_read_ints
 implicit none
 BEGIN_DOC
! routine that computes the variational energy of the current wave function 
!
! and saves the current one-body density matrix into the aux_quantities folder the EZFIO folder
 END_DOC
 density_for_dft = "WFT" 
 touch density_for_dft 
 call print_variational_energy_dft
 call ezfio_set_aux_quantities_data_one_e_dm_alpha_mo(one_e_dm_mo_alpha)
 call ezfio_set_aux_quantities_data_one_e_dm_beta_mo(one_e_dm_mo_beta)

end


