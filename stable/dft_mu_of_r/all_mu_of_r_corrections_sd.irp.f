program all_mu_of_r_corrections
 implicit none
 read_wf = .true.
 touch read_wf

!print*,'**********************'
!print*,'**********************'
!print*,'LDA and PBE / mu(r) HF coallescence with full density and interaction '
!print*,'**********************'
!mu_of_r_functional ="basis_set_PBE"
!touch mu_of_r_functional
!mu_of_r_potential = "hf_coallescence"
!touch mu_of_r_potential 
!no_core_density = "pouet"
!touch no_core_density
!call print_contribution_dft_mu_of_r
 print*,'**********************'
 print*,'**********************'
 print*,'**********************'
 print*,'LDA and PBE / mu(r) HF coallescence with frozen core density and interaction '
 print*,'**********************'
 mu_of_r_functional ="basis_set_PBE"
 touch mu_of_r_functional
 mu_of_r_potential = "hf_valence_coallescence"
 touch mu_of_r_potential 
 no_core_density = "no_core_dm"
 touch no_core_density
 call print_contribution_dft_mu_of_r
 print*,'**********************'

!print*,'LDA and PBE / mu(r) HF/AO coallescence with frozen core density and interaction '
!print*,'**********************'
!mu_of_r_functional ="basis_set_PBE"
!touch mu_of_r_functional
!mu_of_r_potential = "hf_valence_coallescence_ao"
!touch mu_of_r_potential 
!no_core_density = "no_core_dm"
!touch no_core_density
!call print_contribution_dft_mu_of_r
!print*,'**********************'


end

