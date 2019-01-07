program write_integrals_for_dft
 implicit none
 read_wf = .true.
 touch read_wf
 disk_access_mo_one_integrals = "None"
 touch disk_access_only_mo_one_integrals
 disk_access_mo_integrals = "None"
 touch disk_access_mo_integrals
 disk_access_ao_integrals = "None"
 touch disk_access_ao_integrals


 print*,'**********************'
 print*,'**********************'
 print*,'LDA / HF coallescence'
 md_correlation_functional ="basis_set_LDA"
 touch md_correlation_functional
 mu_of_r_potential = "hf_coallescence"
 touch mu_of_r_potential 
 call print_contribution_dft_mu_of_r
 print*,'**********************'
 print*,'**********************'
 print*,'LDA and PBE / HF coallescence'
 md_correlation_functional ="basis_set_on_top_PBE"
 touch md_correlation_functional
 mu_of_r_potential = "hf_coallescence" 
 touch mu_of_r_potential 
 call print_contribution_dft_mu_of_r

 print*,'**********************'
 print*,'**********************'
 print*,'LDA and PBE / PSI coallescence'
 md_correlation_functional ="basis_set_on_top_PBE"
 touch md_correlation_functional
 mu_of_r_potential = "psi_coallescence"
 touch mu_of_r_potential 
 call print_contribution_dft_mu_of_r

end

