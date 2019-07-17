 subroutine print_ecmd_var_energy
 implicit none
 BEGIN_DOC
! routines that prints the variational energy with the ECMD correction
 END_DOC
 print*,'/////////////////////////'
  print*,  '****************************************'
  print*,'///////////////////'
  print*,  ' Regular range separated DFT energy '
  write(*, '(A22,X,F32.10)') 'mu_erf_dft              = ',mu_erf_dft          
  write(*, '(A22,X,F16.10)') 'TOTAL ENERGY        = ',psi_energy + Energy_c_md_on_top_PBE_mu_vector  + nuclear_repulsion 
  print*, ''
  print*, 'Component of the energy ....'
  print*, ''
 
  write(*, '(A22,X,F16.10)') '<Psi| H | Psi>      = ',psi_energy
  write(*, '(A22,X,F16.10)') 'Energy_c_md_on_top_PBE_mu_vector  = ',Energy_c_md_on_top_PBE_mu_vector
  print*,  '****************************************'
  write(*, '(A22,X,F16.10)') 'psi_energy_two_e    = ',psi_energy_two_e
  write(*, '(A22,X,F16.10)') 'psi_energy_h_core   = ',psi_energy_h_core
  print*,  '****************************************'
  print*,'Test for the coherence between density and wave function used'
  print*,'psi_energy_h_core - psi_dft_energy_h_core = ',psi_energy_h_core - psi_dft_energy_h_core
 end


!subroutine print_ecmd_var_energy_mu_of_r
!implicit none
!BEGIN_DOC
! routines that prints the variational energy with the ECMD correction
!END_DOC
!print*,'/////////////////////////'
! print*,  '****************************************'
! print*,'///////////////////'
! print*,  ' Regular range separated DFT energy '
! write(*, '(A22,X,F32.10)') 'mu_erf_dft              = ',mu_erf_dft          
! write(*, '(A22,X,F16.10)') 'TOTAL ENERGY mu of r        = ',psi_energy + Energy_c_md_on_top_PBE_mu_of_r_rsdft + nuclear_repulsion 
! print*, ''
! print*, 'Component of the energy ....'
! print*, ''
!
! write(*, '(A22,X,F16.10)') '<Psi| H | Psi>      = ',psi_energy
! write(*, '(A22,X,F16.10)') 'Energy_c_md_on_top_PBE_mu_of_r_rsdft  = ', Energy_c_md_on_top_PBE_mu_of_r_rsdft
! print*,  '****************************************'
! write(*, '(A22,X,F16.10)') 'psi_energy_two_e    = ',psi_energy_two_e
! write(*, '(A22,X,F16.10)') 'psi_energy_h_core   = ',psi_energy_h_core
! print*,  '****************************************'
! print*,'Test for the coherence between density and wave function used'
! print*,'psi_energy_h_core - psi_dft_energy_h_core = ',psi_energy_h_core - psi_dft_energy_h_core
!end



!subroutine print_ecmd_var_energy_barth
!implicit none
!BEGIN_DOC
! routines that prints the variational energy with the ECMD correction
!END_DOC
!print*,'/////////////////////////'
! print*,  '****************************************'
! print*,'///////////////////'
! print*,  ' Regular range separated DFT energy '
! print*, 'mu_erf_dft                 = ',mu_erf_dft          
! print*, 'TOTAL ENERGY PBE           = ',psi_energy + Energy_c_md_on_top_PBE_mu_vector  + nuclear_repulsion 
! print*, 'TOTAL ENERGY mu_of_r PBE   = ',psi_energy + Energy_c_md_on_top_PBE_mu_of_r_rsdft  + nuclear_repulsion 
! print*, '********* FOR LDA Calculation*****'
! print*, 'TOTAL ENERGY LDA           = ',psi_energy + Energy_c_md_LDA + nuclear_repulsion
! print*, 'TOTAL ENERGY mu_of_r LDA   = ',psi_energy + Energy_c_md_mu_of_r_LDA_rsdft  + nuclear_repulsion 
! print*, ''
! print*, 'Component of the energy ....'
! print*, ''
!
! print*, '<Psi| H | Psi>      = ',psi_energy
!!print*, 'Energy_c_md_on_top_PBE_mu_vector  = ',Energy_c_md_on_top_PBE_mu_vector
!!print*, 'Energy_c_md_on_top_PBE_mu_of_r    = ',Energy_c_md_on_top_PBE_mu_of_r_rsdft
! print*,  '****************************************'
! write(*, '(A22,X,F16.10)') 'psi_energy_two_e    = ',psi_energy_two_e
! write(*, '(A22,X,F16.10)') 'psi_energy_h_core   = ',psi_energy_h_core
! print*,  '****************************************'
! print*,'Test for the coherence between density and wave function used'
! print*,'psi_energy_h_core - psi_dft_energy_h_core = ',psi_energy_h_core - psi_dft_energy_h_core
!end
