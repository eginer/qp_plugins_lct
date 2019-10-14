!BEGIN_PROVIDER [double precision, short_range_Hartree_operator, (mo_tot_num,mo_tot_num)]
!BEGIN_PROVIDER [double precision, short_range_Hartree]
!implicit none
!BEGIN_DOC
! short_range_Hartree_operator(i,j) = \int dr i(r)j(r) \int r' \rho(r') W_{ee}^{sr}
! short_range_Hartree = 0.5 * \sum_{i,j} \rho_{ij} short_range_Hartree_operator(i,j) 
!                     = 0.5 * \int dr \int r' \rho(r) \rho(r') W_{ee}^{sr}
!END_DOC
!integer :: i,j,k,l,m,n
!double precision :: get_mo_bielec_integral,get_mo_bielec_integral_erf
!double precision :: integral, integral_erf, contrib
!short_range_Hartree_operator = 0.d0
!short_range_Hartree = 0.d0
!do i = 1, mo_tot_num
! do j = 1, mo_tot_num
!  if(dabs(one_body_dm_average_mo_for_dft(i,j)).le.1.d-10)cycle
!  do k = 1, mo_tot_num
!   do l = 1, mo_tot_num
!    integral = get_mo_bielec_integral(i,k,j,l,mo_integrals_map) ! <ik|jl> = (ij|kl)
!    integral_erf = get_mo_bielec_integral_erf(i,k,j,l,mo_integrals_erf_map)
!    contrib = one_body_dm_average_mo_for_dft(i,j) * (integral  - integral_erf)
!    short_range_Hartree_operator(l,k) += contrib 
!    short_range_Hartree += contrib * one_body_dm_average_mo_for_dft(k,l) 
!   enddo
!  enddo
! enddo
!enddo
!short_range_Hartree = short_range_Hartree * 0.5d0
!print*, 'short_range_Hartree',short_range_Hartree
!ND_PROVIDER


!BEGIN_PROVIDER [double precision, effective_one_e_potential, (mo_tot_num, mo_tot_num,N_states)]
!BEGIN_PROVIDER [double precision, effective_one_e_potential_without_kin, (mo_tot_num, mo_tot_num,N_states)]
!implicit none
!integer :: i,j,i_state
!effective_one_e_potential = 0.d0
!BEGIN_DOC 
! effective_one_e_potential(i,j) = <i| v_{H}^{sr} |j> + <i| h_{core} |j> + <i|v_{xc} |j> 
! Taking the expectation value does not provide any energy
! but effective_one_e_potential(i,j) is the potential coupling DFT and WFT part to be used in any WFT calculation
!END_DOC
!do i_state = 1, N_states
! do i = 1, mo_tot_num
!  do j = 1, mo_tot_num
!   effective_one_e_potential(i,j,i_state) = short_range_Hartree_operator(i,j) + mo_nucl_elec_integral(i,j) + mo_kinetic_integral(i,j) & 
!                                  + 0.5d0 * (potential_x_alpha_mo(i,j,i_state) + potential_c_alpha_mo(i,j,i_state)                               &
!                                  +          potential_x_beta_mo(i,j,i_state) + potential_c_beta_mo(i,j,i_state)   )
!   effective_one_e_potential_without_kin(i,j,i_state) = short_range_Hartree_operator(i,j) + mo_nucl_elec_integral(i,j)  & 
!                                  + 0.5d0 * (potential_x_alpha_mo(i,j,i_state) + potential_c_alpha_mo(i,j,i_state)                               &
!                                  +          potential_x_beta_mo(i,j,i_state)  + potential_c_beta_mo(i,j,i_state)   )
!  enddo
! enddo
!enddo
!ND_PROVIDER 

 BEGIN_PROVIDER [double precision, effective_one_e_potential_ecmd_lda, (mo_num,mo_num,N_states)]
&BEGIN_PROVIDER [double precision, effective_one_e_potential_without_kin_ecmd_lda, (mo_num, mo_num,N_states)]
 implicit none
 integer :: i,j,i_state
 effective_one_e_potential_ecmd_lda = 0.d0
 BEGIN_DOC
! effective_one_e_potential(i,j) = <i| h_{core} |j> + <i| v_{ecmd,LDA} |j>
! Taking the expectation value does not provide any energy but effective_one_e_potential(i,j) is the potential coupling DFT and WFT part to be used in any WFT calculation
 END_DOC
 do i_state = 1, N_states
  do i = 1, mo_num
   do j = 1, mo_num
    effective_one_e_potential_ecmd_lda(i,j,i_state) = mo_integrals_n_e(i,j) + mo_kinetic_integrals(i,j) &
                                   + 0.5d0 * (potential_deltarho_ecmd_alpha_mo(i,j,i_state) + potential_deltarho_ecmd_beta_mo(i,j,i_state)  &
                                   + potential_e_c_lda_ecmd_alpha_mo(i,j,i_state) + potential_e_c_lda_ecmd_beta_mo(i,j,i_state)   )
    effective_one_e_potential_without_kin_ecmd_lda(i,j,i_state) =  mo_integrals_n_e(i,j)  &
                                   + 0.5d0 * (potential_deltarho_ecmd_alpha_mo(i,j,i_state) + potential_deltarho_ecmd_beta_mo(i,j,i_state)  &
                                   + potential_e_c_lda_ecmd_alpha_mo(i,j,i_state) + potential_e_c_lda_ecmd_beta_mo(i,j,i_state)   ) 
   enddo
  enddo
 enddo
END_PROVIDER



 BEGIN_PROVIDER [double precision, effective_one_e_potential_ecmd_pbe_ueg, (mo_num,mo_num,N_states)]
&BEGIN_PROVIDER [double precision, effective_one_e_potential_without_kin_ecmd_pbe_ueg, (mo_num, mo_num,N_states)]
 implicit none
 integer :: i,j,i_state
 effective_one_e_potential_ecmd_pbe_ueg = 0.d0
 BEGIN_DOC
! effective_one_e_potential(i,j) = <i| h_{core} |j> + <i| v_{ecmd,LDA} |j>
! Taking the expectation value does not provide any energy but effective_one_e_potential(i,j) is the potential coupling DFT and WFT part to be used in any WFT calculation
 END_DOC
 do i_state = 1, N_states
  do i = 1, mo_num
   do j = 1, mo_num
    effective_one_e_potential_ecmd_pbe_ueg(i,j,i_state) = mo_integrals_n_e(i,j) + mo_kinetic_integrals(i,j) &
                                   + 0.5d0 * (potential_c_alpha_mo_sr_pbe_ueg(i,j,i_state) + potential_c_beta_mo_sr_pbe_ueg(i,j,i_state))
    effective_one_e_potential_without_kin_ecmd_pbe_ueg(i,j,i_state) =  mo_integrals_n_e(i,j)  &
                                   + 0.5d0 * (potential_c_alpha_mo_sr_pbe_ueg(i,j,i_state) + potential_c_beta_mo_sr_pbe_ueg(i,j,i_state))
   enddo
  enddo
 enddo
END_PROVIDER



!EGIN_PROVIDER [double precision, Fock_matrix_expectation_value]
!implicit none
! call get_average(effective_one_e_potential,one_body_dm_average_mo_for_dft,Fock_matrix_expectation_value)

!ND_PROVIDER 

!BEGIN_PROVIDER [double precision, Trace_v_xc, (N_states)]
!BEGIN_PROVIDER [double precision, Trace_v_Hxc, (N_states)]
!implicit none
!integer :: i,j,istate
!double precision :: tmp(mo_tot_num,mo_tot_num)
!BEGIN_DOC 
! Trace_v_xc  = \sum_{i,j} rho_{ij} v^{xc}_{ij} 
! Trace_v_Hxc = \sum_{i,j} rho_{ij} v^{Hxc}_{ij} 
!END_DOC
! WARNING: I think there is a bug it potential_alpha should be contracted with density_matrix_alpha and 
! potential_beta should be contracted with density_matrix_beta for opne-shell systems ! JT
!print *,'WARNING: Trace_v_xc is wrong for open-shell systems: MUST BE CORRECTED'
!do istate = 1, N_states
! tmp = 0.d0
! do i = 1, mo_tot_num
!  do j = 1, mo_tot_num
!    tmp(i,j) =   + 0.5d0 * (potential_x_alpha_mo(i,j,istate) + potential_c_alpha_mo(i,j,istate)&
!                 +          potential_x_beta_mo(i,j,istate)  + potential_c_beta_mo(i,j,istate)   )
!  enddo
! enddo
! call get_average(tmp,one_body_dm_mo_for_dft(1,1,istate),Trace_v_xc(istate))

! tmp = 0.d0
! do i = 1, mo_tot_num
!  do j = 1, mo_tot_num
!    tmp(i,j) =    short_range_Hartree_operator(j,i) & 
!                + 0.5d0 * (potential_x_alpha_mo(j,i,istate) + potential_c_alpha_mo(j,i,istate)&
!                     +     potential_x_beta_mo(j,i,istate)  + potential_c_beta_mo(j,i,istate)   )
!  enddo
! enddo
! call get_average(tmp,one_body_dm_mo_for_dft(1,1,istate),Trace_v_Hxc(istate))
!enddo

!ND_PROVIDER 

!EGIN_PROVIDER [double precision, one_e_energy_potential, (mo_tot_num, mo_tot_num)]
!implicit none
!integer :: i,j,i_state
!BEGIN_DOC 
! one_e_energy_potential(i,j) = <i|h_{core}|j> + \int dr i(r)j(r) \int r' \rho(r') W_{ee}^{sr}
! If one take the expectation value over Psi, one gets the total one body energy
!END_DOC
! do i = 1, mo_tot_num
!  do j = 1, mo_tot_num
!   one_e_energy_potential(i,j) = mo_nucl_elec_integral(i,j) + mo_kinetic_integral(i,j) + short_range_Hartree_operator(i,j) * 0.5d0
!  enddo
! enddo

!ND_PROVIDER 

!!!!!!!!!!!!!!test pour manu RS DFT !!!

!BEGIN_PROVIDER [double precision, exp_value_V_SR_mu, (N_states)]
!implicit none
!integer :: i,j,i_state
!exp_value_V_SR_mu = 0.d0
!do i_state = 1, N_states
! do j = 1, mo_tot_num
!  do i = 1, mo_tot_num
!   exp_value_V_SR_mu(i_state) += (one_body_dm_mo_alpha(i,j,i_state)+one_body_dm_mo_beta(i,j,i_state))*(short_range_Hartree_operator(i,j)+ 0.5d0 * (potential_x_alpha_mo(i,j,i_state)+potential_c_alpha_mo(i,j,i_state)+potential_x_beta_mo(i,j,i_state)+potential_c_beta_mo(i,j,i_state)))


!  enddo
! enddo
!enddo
!ND_PROVIDER
