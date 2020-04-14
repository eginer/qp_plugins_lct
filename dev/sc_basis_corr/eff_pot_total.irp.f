 BEGIN_PROVIDER [double precision, mo_tot_eff_pot_basis, (mo_num, mo_num,N_states)]
&BEGIN_PROVIDER [double precision, mo_tot_eff_pot_basis_no_kin, (mo_num, mo_num,N_states)]
 implicit none
 integer :: i,j,istate
 mo_tot_eff_pot_basis = 0.d0
 BEGIN_DOC
! mo_tot_eff_pot_basis(i,j) = $\rangle i_{MO}| v_{H}^{sr} |j_{MO}\rangle  + \rangle i_{MO}| h_{core} |j_{MO}\rangle  + \rangle i_{MO}|v_{xc} |j_{MO}\rangle$
!
! on the |MO| basis
! 
! Taking the expectation value does not provide any energy, but
!
! mo_tot_eff_pot_basis(i,j) is the potential coupling DFT and WFT parts 
!
! and it is used in any RS-DFT based calculations  
 END_DOC
 do istate = 1, N_states
  do j = 1, mo_num
   do i = 1, mo_num

    mo_tot_eff_pot_basis(i,j,istate) =  mo_integrals_n_e(i,j) + mo_kinetic_integrals(i,j)   & 
    + 0.5d0 * ( pot_basis_alpha_mo(i,j,istate) + pot_basis_alpha_mo(i,j,istate) )

    mo_tot_eff_pot_basis_no_kin(i,j,istate) = mo_integrals_n_e(i,j) &
    + 0.5d0 * ( pot_basis_alpha_mo(i,j,istate) + pot_basis_alpha_mo(i,j,istate) )
   enddo
  enddo
 enddo
END_PROVIDER


 BEGIN_PROVIDER [double precision, ao_tot_eff_pot_basis, (ao_num, ao_num,N_states)]
&BEGIN_PROVIDER [double precision, ao_tot_eff_pot_basis_no_kin, (ao_num, ao_num,N_states)]
 implicit none
 BEGIN_DOC
! ao_tot_eff_pot_basis(i,j) = $\rangle i_{AO}| v_{H}^{sr} |j_{AO}\rangle  + \rangle i_{AO}| h_{core} |j_{AO}\rangle  + \rangle i_{AO}|v_{xc} |j_{AO}\rangle$
!
 END_DOC

 integer :: istate

 do istate = 1, N_states
  call mo_to_ao_no_overlap(mo_tot_eff_pot_basis(1,1,istate),size(mo_tot_eff_pot_basis,1),ao_tot_eff_pot_basis(1,1,istate),size(ao_tot_eff_pot_basis,1))

  call mo_to_ao_no_overlap(mo_tot_eff_pot_basis_no_kin(1,1,istate),size(mo_tot_eff_pot_basis_no_kin,1),ao_tot_eff_pot_basis_no_kin(1,1,istate),size(ao_tot_eff_pot_basis_no_kin,1))
 enddo

END_PROVIDER 
