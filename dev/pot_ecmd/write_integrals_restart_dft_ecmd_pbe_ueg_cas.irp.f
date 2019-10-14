program write_integrals_for_dft_ecmd_pbe_ueg
 implicit none
 read_wf = .true.
 touch read_wf
 mu_of_r_potential = "psi_cas_truncated"
 touch mu_of_r_potential 
 on_top_from_cas = .True.
 touch on_top_from_cas
 no_core_density = "no_core_dm"
 touch no_core_density
 io_mo_one_e_integrals = "None"
 touch io_mo_one_e_integrals 
 io_mo_integrals_e_n = "None"
 touch io_mo_integrals_e_n 
 io_ao_integrals_e_n = "None"
 touch io_ao_integrals_e_n 
 call write_all_integrals_for_mrdft_ecmd_pbe_ueg

 !call print_contribution_dft_mu_of_r 
  call print_ecmd_var_energy_barth
end


