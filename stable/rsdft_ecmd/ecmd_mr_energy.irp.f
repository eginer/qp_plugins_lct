 subroutine print_ecmd_var_energy
 implicit none
 BEGIN_DOC
! routines that prints the variational energy with the ECMD correction
 END_DOC
 print*,'/////////////////////////'
  print*,  '****************************************'
  print*,'///////////////////'
  print*,  ' Regular range separated DFT energy '
  write(*, '(A22,X,F32.10)') 'mu_erf_dft               = ',mu_erf_dft          
  write(*, '(A22,X,F16.10)') 'TOTAL ENERGY             = ',psi_energy + ecmd_pbe_on_top_at_mu  + nuclear_repulsion 
  print*,  ''
  print*,  'Component of the energy ....'
  print*,  ''
 
  write(*, '(A22,X,F16.10)') '<Psi| H | Psi>           = ',psi_energy + nuclear_repulsion
  write(*, '(A22,X,F16.10)') 'ecmd_pbe_on_top_at_mu    = ',ecmd_pbe_on_top_at_mu
  print*,  '****************************************'
  write(*, '(A22,X,F16.10)') 'psi_energy_two_e         = ',psi_energy_two_e
  write(*, '(A22,X,F16.10)') 'psi_energy_h_core        = ',psi_energy_h_core
  print*,  '****************************************'
  print*,  'Test for the coherence between density and wave function used'
  print*,  'psi_energy_h_core - psi_dft_energy_h_core  = ',psi_energy_h_core - psi_dft_energy_h_core
 end

