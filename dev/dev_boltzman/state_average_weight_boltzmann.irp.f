
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
    !state_energy(i) = psi_energy(i)+nuclear_repulsion 
     state_energy(i) = psi_energy(i)
    !state_energy(i) = psi_energy(i) !+ psi_energy_h_core(i) 
     print*,'cipsi_energy(',i,') =',state_energy(i)
    else if (weight_boltz_type .EQ. 'rsdft')then
    !state_energy(i) = nuclear_repulsion + psi_energy_erf(i) + psi_dft_energy_h_core(i) + energy_x(i) + energy_c(i) + short_range_Hartree(i)
     state_energy(i) = psi_energy_erf(i) + psi_energy_h_core(i) + potential_hxc_exp_value_ave(i)
     print*,'rs_dft_energy(',i,') =',state_energy(i)
    endif
   enddo

   state_e_lowest = minval(state_energy)
   state_average_weight_tot = 0.d0
   do i=1,N_states
    delta_ei = state_energy(i)-state_e_lowest
    !rint*,'rs_dft_energy(',i,') =',delta_ei
    state_average_weight_boltzmann(i) = exp(-beta_state_av*delta_ei)
    state_average_weight_tot += state_average_weight_boltzmann(i)
   enddo

   state_average_weight_boltzmann = state_average_weight_boltzmann/state_average_weight_tot

END_PROVIDER



 BEGIN_PROVIDER [double precision, effective_one_e_pot_hxc_state_average, (mo_num,mo_num)]
 implicit none
 integer :: i,j,istate
 effective_one_e_pot_hxc_state_average = 0.d0
 BEGIN_DOC
 ! Effective_one_e_pot_hxc(i,j) = $\rangle i_{MO}| v_{H}^{sr} |j_{MO}\rangle  +  \rangle i_{MO}|v_{xc} |j_{MO}\rangle$
 ! density for DFT must me set at states average blabal
 END_DOC
 do j = 1, mo_num
  do i = 1, mo_num
   effective_one_e_pot_hxc_state_average(i,j) = short_range_Hartree_operator(i,j,1)  &
                                  + 0.5d0 * (potential_x_alpha_mo(i,j,1) + potential_c_alpha_mo(i,j,1)                          &
                                  +          potential_x_beta_mo(i,j,1)  + potential_c_beta_mo(i,j,1)   )
  enddo
 enddo
 END_PROVIDER


 BEGIN_PROVIDER [ double precision, potential_hxc_exp_value_ave, (N_states) ]
 implicit none
 integer :: i,j,k
 BEGIN_DOC
  !Doc
 END_DOC
 potential_hxc_exp_value_ave = 0.d0
 do i = 1, N_states
  do j = 1, mo_num
   do k = 1, mo_num
    potential_hxc_exp_value_ave(i) += effective_one_e_pot_hxc_state_average(k,j) * (one_e_dm_mo_alpha(k,j,i) + one_e_dm_mo_beta(k,j,i))
   enddo
  enddo
 enddo
 END_PROVIDER
