 subroutine print_ecmd_var_energy
 implicit none
 BEGIN_DOC
! routines that prints the variational energy with the ECMD correction
 END_DOC
 provide psi_energy ecmd_pbe_on_top_at_mu ecmd_pbe_ueg_prov ecmd_pbe_ueg_prov 
 print*,'/////////////////////////'
  print*,  '****************************************'
  print*,'///////////////////'
  print*,  ' Regular range separated DFT energy '
  write(*, '(A22,X,F32.10)') 'mu_erf_dft               = ',mu_erf_dft          
  write(*, '(A22,X,F16.10)') 'TOTAL ENERGY             = ',psi_energy + ecmd_pbe_on_top_at_mu  + nuclear_repulsion 
  write(*, '(A22,X,F16.10)') 'Ecmd PBE-on-top          = ',ecmd_pbe_on_top_at_mu
  print*,  ''
  print*,  'Component of the energy ....'
  print*,  ''
 
  write(*, '(A22,X,F16.10)') '<Psi| H | Psi>           = ',psi_energy + nuclear_repulsion
  print*,'Correlation part '
  print*,'On-top functional '
  write(*, '(A22,X,F16.10)') 'ecmd_pbe_on_top_at_mu    = ',ecmd_pbe_on_top_at_mu
  print*,'On-top-UEG functional '
  write(*, '(A22,X,F16.10)') 'Ecmd PBE-UEG             = ',ecmd_pbe_ueg_prov
  write(*, '(A22,X,F16.10)') 'energy_c_md_sr_pbe       = ',energy_c_md_sr_pbe
  print*,  '****************************************'
  print*,'Exchange    part '
  write(*, '(A22,X,F16.10)') ' Exact EXmd energy       = ',psi_energy_wee_sr - short_range_Hartree
  write(*, '(A22,X,F16.10)') 'energy_x_md_sr_pbe       = ',energy_x_md_sr_pbe
  print*,  '****************************************'
  write(*, '(A22,X,F16.10)') ' psi_energy_erf          = ',psi_energy_erf 
  print*,''
  write(*, '(A22,X,F16.10)') 'psi_energy_two_e         = ',psi_energy_two_e
  write(*, '(A22,X,F16.10)') 'psi_energy_h_core        = ',psi_energy_h_core
  print*,  '****************************************'
  print*,  'Test for the coherence between density and wave function used'
  print*,  'psi_energy_h_core - psi_dft_energy_h_core  = ',psi_energy_h_core - psi_dft_energy_h_core
 end


 subroutine print_many_energy_component
 implicit none
 BEGIN_DOC
! routines that prints the variational energy with the ECMD correction
 END_DOC
 double precision :: energy(N_states)
 provide psi_energy ecmd_pbe_on_top_at_mu ecmd_pbe_ueg_prov ecmd_pbe_ueg_prov energy_x_md_sr_pbe
 print*,'/////////////////////////'
  print*,  '****************************************'
  write(*, '(A22,X,F32.10)') 'mu_erf_dft               = ',mu_erf_dft          
  print*,'///////////////////'
  !       <Psi | T + V_ne        +  W_ee^lr      | Psi >  + E_x^sr  +  E_c^sr         + E_H^sr 
  energy = psi_dft_energy_h_core + psi_energy_erf + energy_x_sr_pbe + energy_c_sr_pbe +  short_range_Hartree
  write(*, '(A22,X,F16.10)') 'Usual sr_pbe             = ',energy + nuclear_repulsion 
  energy = psi_energy + ecmd_pbe_on_top_at_mu  
  write(*, '(A22,X,F16.10)') 'PBE-OT energy            = ',energy + nuclear_repulsion 
  write(*, '(A22,X,F16.10)') 'Test                     = ',ecmd_pbe_on_top_at_mu + nuclear_repulsion 
  !       <Psi | T + V_ne        +  W_ee^lr      | Psi >  + E_xmd^sr   +  E_cmd^sr         + E_H^sr 
  energy = psi_dft_energy_h_core + psi_energy_erf + energy_x_md_sr_pbe + energy_c_md_sr_pbe +  short_range_Hartree
  write(*, '(A22,X,F16.10)') 'PBE-UEG MD               = ',energy + nuclear_repulsion
  print*,  ''
  print*,  'Component of the energy ....'
  print*,  ''
 
  write(*, '(A22,X,F16.10)') '<Psi| H | Psi>           = ',psi_energy + nuclear_repulsion
  write(*, '(A22,X,F16.10)') 'Short range Hartree      = ',short_range_Hartree
  print*,'Correlation part '
  print*,'On-top functional '
  write(*, '(A22,X,F16.10)') 'ecmd_pbe_on_top_at_mu    = ',ecmd_pbe_on_top_at_mu
  print*,'On-top-UEG functional '
  write(*, '(A22,X,F16.10)') 'Ecmd PBE-UEG             = ',ecmd_pbe_ueg_prov
  write(*, '(A22,X,F16.10)') 'energy_c_md_sr_pbe       = ',energy_c_md_sr_pbe
  print*,  '****************************************'
  print*,'Exchange    part '
  write(*, '(A22,X,F16.10)') ' Exact EXmd energy       = ',psi_energy_wee_sr - short_range_Hartree
  write(*, '(A22,X,F16.10)') 'energy_x_md_sr_pbe       = ',energy_x_md_sr_pbe
  print*,  '****************************************'
  write(*, '(A22,X,F16.10)') ' psi_energy_erf          = ',psi_energy_erf 
  print*,''
  write(*, '(A22,X,F16.10)') 'psi_energy_two_e         = ',psi_energy_two_e
  write(*, '(A22,X,F16.10)') 'psi_energy_h_core        = ',psi_energy_h_core
  print*,  '****************************************'
  print*,  'Test for the coherence between density and wave function used'
  print*,  'psi_energy_h_core - psi_dft_energy_h_core  = ',psi_energy_h_core - psi_dft_energy_h_core
 end

