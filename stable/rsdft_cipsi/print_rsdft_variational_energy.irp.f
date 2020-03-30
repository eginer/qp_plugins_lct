program DFT_Utils_two_body_main
 implicit none
 read_wf = .true.
 touch read_wf
 io_mo_one_e_integrals = "None"
 touch io_mo_one_e_integrals
 io_mo_two_e_integrals = "None"
 touch io_mo_two_e_integrals
 io_ao_two_e_integrals = "None"
 touch io_ao_two_e_integrals

 io_mo_integrals_n_e = "None"
 touch io_mo_integrals_n_e
 io_mo_integrals_kinetic = "None"
 touch io_mo_integrals_kinetic 
 io_ao_integrals_n_e = "None"
 touch io_ao_integrals_n_e 
 io_ao_integrals_kinetic = "None"
 touch io_ao_integrals_kinetic 

 density_for_dft = "WFT" 
 touch density_for_dft 
 call print_variational_energy_dft

end
