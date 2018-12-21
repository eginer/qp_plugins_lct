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
  write(*, '(A22,X,F16.10)') 'TOTAL ENERGY        = ',psi_energy + Energy_c_md_on_top  + nuclear_repulsion 
  print*, ''
  print*, 'Component of the energy ....'
  print*, ''
 
  write(*, '(A22,X,F16.10)') '<Psi| H | Psi>      = ',psi_energy
  write(*, '(A22,X,F16.10)') 'Energy_c_md_on_top  = ',Energy_c_md_on_top
  print*,  '****************************************'
  write(*, '(A22,X,F16.10)') 'psi_energy_bielec   = ',psi_energy_bielec
  write(*, '(A22,X,F16.10)') 'psi_energy_h_core   = ',psi_energy_h_core
  print*,  '****************************************'
  print*,'Test for the coherence between density and wave function used'
  print*,'psi_energy_h_core - psi_dft_energy_h_core = ',psi_energy_h_core - psi_dft_energy_h_core
 end

