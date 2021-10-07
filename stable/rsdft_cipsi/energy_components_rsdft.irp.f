
BEGIN_PROVIDER [ double precision, var_rsdft_energy, (N_states) ]
  implicit none
  BEGIN_DOC
! Total_range_separated_electronic_energy = <Psi| h_{core} |Psi> + (1/2) <Psi| v_{H}^{sr} |Psi> + <i|W_{ee}^{lr}|i> + E_{x} + E_{c}
  END_DOC
   var_rsdft_energy = psi_dft_energy_h_core + short_range_Hartree + psi_energy_erf + energy_x + energy_c
END_PROVIDER


BEGIN_PROVIDER [ double precision, two_elec_energy_dft, (N_states) ]
  implicit none
  BEGIN_DOC
! two_elec_energy_dft = (1/2) <Psi| v_{H}^{sr} |Psi> + <Psi|W_{ee}^{lr}|Psi>
  END_DOC
   two_elec_energy_dft = short_range_Hartree + psi_energy_erf
END_PROVIDER


BEGIN_PROVIDER [double precision, electronic_energy_mr_dft, (N_states)]
 implicit none
 BEGIN_DOC
 ! Energy for the multi determinantal DFT calculation
 END_DOC

 print*,'You are using a variational method which uses the wave function stored in the EZFIO folder'
 electronic_energy_mr_dft = var_rsdft_energy
END_PROVIDER


 subroutine print_variational_energy_dft
 implicit none
 BEGIN_DOC
! Routines that prints the variational energy, and many more quantities
 END_DOC
 integer :: i
 print*,'/////////////////////////'
 print*,'mu_erf_dft              = ',mu_erf_dft
 print*,'mu_of_r_dft_average     = ',mu_of_r_dft_average
 print*,'nuclear_repulsion   = ',nuclear_repulsion
 do i=1,N_states
  print*,  '****************************************'
  print*,'///////////////////'
  print*,  ' Regular range separated DFT energy '
  print*,'State ',i
  print*,'TOTAL ENERGY        = ',electronic_energy_mr_dft(i)+nuclear_repulsion
  print*, ''
  print*, 'Component of the energy ....'
  print*, ''
  print*,'psi_energy_erf      = ',psi_energy_erf(i)
  print*,'psi_dft_energy_h_cor= ',psi_dft_energy_h_core(i)
  print*,'short_range_Hartree = ',short_range_Hartree(i)
  print*,'two_elec_energy     = ',two_elec_energy_dft(i)
  print*,'energy_x            = ',energy_x(i)
  print*,'energy_c            = ',energy_c(i)
  print*,'E_xc                = ',energy_x(i) + energy_c(i)
  print*,'E_Hxc               = ',energy_x(i) + energy_c(i) + short_range_Hartree(i)
  print*, ''
  print*,  '****************************************'
  print*, ''
  write(*, '(A22,X,F16.10)') 'Trace_v_Hxc         = ',Trace_v_Hxc(i)
 enddo
 end

