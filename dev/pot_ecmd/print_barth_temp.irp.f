program DFT_Utils_ECMD
 implicit none
 read_wf = .true.
 touch read_wf
 no_core_density = .True.
 touch no_core_density
 io_mo_one_e_integrals = "None"
 touch io_mo_one_e_integrals 
 io_mo_integrals_n_e = "None"
 touch io_mo_integrals_n_e 
 io_ao_integrals_n_e = "None"
 touch io_ao_integrals_n_e 
 call print_ecmd_var_energy_barth 

!call print_z_dipole_moment_only  
end


