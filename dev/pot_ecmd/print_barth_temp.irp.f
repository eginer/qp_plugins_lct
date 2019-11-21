program DFT_Utils_ECMD
 implicit none
 read_wf = .true.
 touch read_wf
 no_core_density = "no_core_dm"
 touch no_core_density
 io_mo_one_e_integrals = "None"
 touch io_mo_one_e_integrals 
 io_mo_integrals_e_n = "None"
 touch io_mo_integrals_e_n 
 io_ao_integrals_e_n = "None"
 touch io_ao_integrals_e_n 
 call print_ecmd_var_energy_barth 

!call print_z_dipole_moment_only  
end


