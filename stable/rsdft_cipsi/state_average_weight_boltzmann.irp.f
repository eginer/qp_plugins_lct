BEGIN_PROVIDER [ double precision, state_average_weight_boltzmann, (N_states) ]
   implicit none
   BEGIN_DOC
   ! Weights in the state-average blabla 
   END_DOC
   integer :: i
   double precision :: beta_state_av,delta_ei,state_average_weight_tot,state_e_lowest
   double precision, allocatable :: state_energy(:)
   allocate(state_energy(N_states))
   beta_state_av = 10.d0

   do i=1,N_states
    if(weight_boltz_type .EQ. 'cipsi')then
     state_energy(i) = psi_energy(i)+nuclear_repulsion 
    !state_energy(i) = psi_energy(i) !+ psi_energy_h_core(i) 
     print*,'cipsi_energy(',i,') =',state_energy(i)
    else if (weight_boltz_type .EQ. 'rsdft')then
     state_energy(i) = nuclear_repulsion + psi_energy_erf(i) + psi_dft_energy_h_core(i) + energy_x(i) + energy_c(i) + short_range_Hartree(i)
     print*,'rs_dft_energy(',i,') =',state_energy(i)
    endif
   enddo

   state_e_lowest = minval(state_energy)

   do i=1,N_states
    delta_ei = state_energy(i)-state_e_lowest
    !rint*,'rs_dft_energy(',i,') =',delta_ei
    state_average_weight_boltzmann(i) = exp(-beta_state_av*delta_ei)
    state_average_weight_tot += state_average_weight_boltzmann(i)
   enddo

   state_average_weight_boltzmann = state_average_weight_boltzmann/state_average_weight_tot

END_PROVIDER

