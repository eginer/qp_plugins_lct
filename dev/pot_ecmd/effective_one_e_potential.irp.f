 BEGIN_PROVIDER [double precision, effective_one_e_potential_ecmd_lda, (mo_num,mo_num,N_states)]
&BEGIN_PROVIDER [double precision, effective_one_e_potential_without_kin_ecmd_lda, (mo_num, mo_num,N_states)]
 implicit none
 integer :: i,j,istate
 effective_one_e_potential_ecmd_lda = 0.d0
 BEGIN_DOC
! effective_one_e_potential(i,j) = <i| h_{core} |j> + <i| v_{ecmd,LDA} |j>
! Taking the expectation value does not provide any energy but effective_one_e_potential(i,j) is the potential coupling DFT and WFT part to be used in any WFT calculation
 END_DOC
 do istate = 1, N_states
  do i = 1, mo_num
   do j = 1, mo_num
    effective_one_e_potential_ecmd_lda(i,j,istate) = mo_integrals_n_e(i,j) + mo_kinetic_integrals(i,j) &
                                   + 0.5d0 * (pot_basis_alpha_mo_lda(j,i,istate) + pot_basis_beta_mo_lda(j,i,istate) )
    effective_one_e_potential_without_kin_ecmd_lda(i,j,istate) =  mo_integrals_n_e(i,j)  &
                                   + 0.5d0 * (pot_basis_alpha_mo_lda(j,i,istate) + pot_basis_beta_mo_lda(j,i,istate) )
   enddo
  enddo
 enddo
END_PROVIDER



 BEGIN_PROVIDER [double precision, effective_one_e_potential_ecmd_pbe_ueg, (mo_num,mo_num,N_states)]
&BEGIN_PROVIDER [double precision, effective_one_e_potential_without_kin_ecmd_pbe_ueg, (mo_num, mo_num,N_states)]
 implicit none
 integer :: i,j,istate
 effective_one_e_potential_ecmd_pbe_ueg = 0.d0
 BEGIN_DOC
! effective_one_e_potential(i,j) = <i| h_{core} |j> + <i| v_{ecmd,LDA} |j>
! Taking the expectation value does not provide any energy but effective_one_e_potential(i,j) is the potential coupling DFT and WFT part to be used in any WFT calculation
 END_DOC
 do istate = 1, N_states
  do i = 1, mo_num
   do j = 1, mo_num
    effective_one_e_potential_ecmd_pbe_ueg(i,j,istate) = mo_integrals_n_e(i,j) + mo_kinetic_integrals(i,j) &
                                   + 0.5d0 * (pot_basis_alpha_mo_pbe_ueg(i,j,istate) + pot_basis_beta_mo_pbe_ueg(i,j,istate))
    effective_one_e_potential_without_kin_ecmd_pbe_ueg(i,j,istate) =  mo_integrals_n_e(i,j)  &
                                   + 0.5d0 * (pot_basis_alpha_mo_pbe_ueg(i,j,istate) + pot_basis_beta_mo_pbe_ueg(i,j,istate))
   enddo
  enddo
 enddo
END_PROVIDER
