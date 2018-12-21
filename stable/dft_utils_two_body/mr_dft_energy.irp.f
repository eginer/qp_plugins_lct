BEGIN_PROVIDER [double precision, electronic_energy_mr_dft, (N_states)]
 implicit none
 BEGIN_DOC
 ! Energy for the multi determinantal DFT calculation
 END_DOC
 
  print*,'You are using a variational method which uses the wave function stored in the EZFIO folder'
  electronic_energy_mr_dft = total_range_separated_electronic_energy


END_PROVIDER 

 subroutine print_variational_energy_dft
 implicit none
 print*,'/////////////////////////'
  print*,  '****************************************'
  print*,'///////////////////'
  print*,  ' Regular range separated DFT energy '
  write(*, '(A22,X,F32.10)') 'mu_erf              = ',mu_erf          
  write(*, '(A22,X,F16.10)') 'TOTAL ENERGY        = ',electronic_energy_mr_dft+nuclear_repulsion
  print*, ''
  print*, 'Component of the energy ....'
  print*, ''
  write(*, '(A22,X,F16.10)') 'nuclear_repulsion   = ',nuclear_repulsion
  write(*, '(A22,X,F16.10)') 'psi_energy_erf      = ',psi_energy_erf      
  write(*, '(A22,X,F16.10)') 'psi_dft_energy_h_cor= ',psi_dft_energy_h_core    
  write(*, '(A22,X,F16.10)') 'short_range_Hartree = ',short_range_Hartree
  write(*, '(A22,X,F16.10)') 'two_elec_energy     = ',two_elec_energy_dft
  write(*, '(A22,X,F16.10)') 'energy_x            = ',energy_x         
  write(*, '(A22,X,F16.10)') 'energy_c            = ',energy_c          
  write(*, '(A22,X,F16.10)') 'E_xc                = ',energy_x + energy_c
  write(*, '(A22,X,F16.10)') 'E_Hxc               = ',energy_x + energy_c + short_range_Hartree
  print*, ''
  print*,  '****************************************'
  print*, ''
  write(*, '(A22,X,F16.10)') 'Approx eigenvalue   = ',electronic_energy_mr_dft+nuclear_repulsion + Trace_v_Hxc - (short_range_Hartree + energy_x + energy_c)
  write(*, '(A22,X,F16.10)') 'Trace_v_xc          = ',Trace_v_xc
  write(*, '(A22,X,F16.10)') 'Trace_v_Hxc         = ',Trace_v_Hxc
 
  write(*, '(A22,X,F16.10)') '<Psi| H | Psi>      = ',psi_energy
  write(*, '(A22,X,F16.10)') 'psi_energy_bielec   = ',psi_energy_bielec
  write(*, '(A22,X,F16.10)') 'psi_energy_h_core   = ',psi_energy_h_core
 end

