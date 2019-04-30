 subroutine print_variational_energy_dft_boltzmann
 implicit none
 BEGIN_DOC
! Routines that prints the variational energy, and many more quantities
 END_DOC
 integer :: i
 double precision :: states_ave_psi_ener,total_boltz_energy
 double precision, allocatable :: psi_boltz_dft_energy(:)
 allocate(psi_boltz_dft_energy(N_states))

 do i=1,N_states
  psi_boltz_dft_energy(i) =  psi_energy_erf(i)  + psi_dft_energy_h_core(i) 
 !print*,'psi_boltz_dft_energy(i)',psi_boltz_dft_energy(i)
 enddo

 print*,  '**************************************        **'
 write(*, '(A22,X,3(F16.10,X))') 'Weight normal           = ',state_average_weight(1),state_average_weight(2),state_average_weight(3)
 print*,  '**************************************        **'

 state_average_weight = state_average_weight_boltzmann
 soft_touch state_average_weight 

 density_for_dft = "state_average_dens"
 soft_touch density_for_dft 
 provide energy_x energy_c short_range_Hartree
 
 states_ave_psi_ener= 0.d0

 do i=1,N_states
  states_ave_psi_ener += psi_boltz_dft_energy(i)*state_average_weight(i) 
 enddo 

 total_boltz_energy = states_ave_psi_ener + energy_x(1) + energy_c(1) + short_range_Hartree(1) 

 print*,  '**************************************        **'
 write(*, '(A22,X,3(F16.10,X))') 'Weight boltz            = ',state_average_weight(1),state_average_weight(2),state_average_weight(3)
 print*,  '**************************************        **'

 print*,'/////////////////////////'
  print*,  '****************************************'
  print*,'///////////////////'
  print*,  ' Regular range separated DFT energy Bolztmann average'
  write(*, '(A22,X,F32.10)') 'mu_erf_dft                           = ',mu_erf_dft          
  write(*, '(A22,X,F16.10)') 'TOTAL ENERGY   Bolztmann average     = ',total_boltz_energy+nuclear_repulsion
  print*, ''
  print*, 'Component of the energy ....'
  print*, ''
  write(*, '(A22,X,F16.10)') 'nuclear_repulsion           = ',nuclear_repulsion
  write(*, '(A22,X,F16.10)') 'states_ave_psi_ener         = ',states_ave_psi_ener
  write(*, '(A22,X,F16.10)') 'short_range_Hartree[n_ave]  = ',short_range_Hartree(1)
  write(*, '(A22,X,F16.10)') 'energy_x[n_ave]             = ',energy_x(1)         
  write(*, '(A22,X,F16.10)') 'energy_c[n_ave]             = ',energy_c(1)          
  write(*, '(A22,X,F16.10)') 'E_xc[n_ave]                 = ',energy_x(1) + energy_c(1)
  write(*, '(A22,X,F16.10)') 'E_Hxc[n_ave]                = ',energy_x(1) + energy_c(1) + short_range_Hartree(1)
  print*, ''                                               
  print*,  '**************************************        **'
  print*, ''                                               
  write(*, '(A22,X,F16.10)') 'Trace_v_xc[n_ave]           = ',Trace_v_xc(1)
  write(*, '(A22,X,F16.10)') 'Trace_v_Hxc[n_ave]          = ',Trace_v_Hxc(1)
  print*,  '****************************************'
 end

