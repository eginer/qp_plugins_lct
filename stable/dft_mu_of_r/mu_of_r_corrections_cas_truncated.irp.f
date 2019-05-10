program all_mu_of_r_corrections
 implicit none
 read_wf = .true.
 touch read_wf
 mu_of_r_potential = "psi_cas_truncated"
 touch mu_of_r_potential 
 no_core_density = "no_core_dm"
 touch no_core_density
 on_top_from_cas = .True.
 touch on_top_from_cas
 call routine_print
end

subroutine routine_print
 implicit none
 print*,'LDA, PBE and PBE-on-top / mu(r) PSI coallescence with frozen core interaction'
 mu_of_r_functional ="basis_set_on_top_PBE"
 touch mu_of_r_functional
 call print_contribution_dft_mu_of_r

end
