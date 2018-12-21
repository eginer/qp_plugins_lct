program DFT_Utils_two_body_main
 implicit none
 read_wf = .true.
 touch read_wf
 disk_access_mo_one_integrals = "None"
 touch disk_access_mo_one_integrals
 disk_access_mo_integrals = "None"
 touch disk_access_mo_integrals
 disk_access_ao_integrals = "None"
 touch disk_access_ao_integrals
 density_for_dft = "WFT" 
 touch density_for_dft 
 call print_variational_energy_dft
 

end
