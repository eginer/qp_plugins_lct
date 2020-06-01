
 BEGIN_PROVIDER [double precision, full_occ_2_rdm_ab_chemist_mo, (n_core_inact_act_orb,n_core_inact_act_orb,n_core_inact_act_orb,n_core_inact_act_orb,N_states)]
 implicit none
 BEGIN_DOC
! full_occ_2_rdm_ab_chemist_mo(i,j,k,l,istate) = (1 1|2 2) 2-Rdm == Chemist notation
 END_DOC
 integer :: i,j,k,l,istate
 do istate = 1, N_states
  do i = 1, n_core_inact_act_orb ! 1 
   do j = 1, n_core_inact_act_orb ! 2 
    do k = 1, n_core_inact_act_orb ! 1
     do l = 1, n_core_inact_act_orb ! 2 
      full_occ_2_rdm_ab_chemist_mo(i,k,j,l,istate) = full_occ_2_rdm_ab_mo(i,j,k,l,istate)
     enddo
    enddo
   enddo
  enddo
 enddo
 END_PROVIDER 
