program write_integrals_for_dft_ecmd_lda
 implicit none
 read_wf = .true.
 touch read_wf
 io_mo_one_e_integrals = "None"
 touch io_mo_one_e_integrals 
 io_mo_integrals_e_n = "None"
 touch io_mo_integrals_e_n 
 io_ao_integrals_e_n = "None"
 touch io_ao_integrals_e_n 
 call write_all_integrals_for_mrdft_ecmd_lda

 !call print_contribution_dft_mu_of_r 
  call print_ecmd_var_energy_barth
end


