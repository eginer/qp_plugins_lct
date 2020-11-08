
subroutine routine_active_only_test(two_rdm)
 implicit none
 double precision, intent(in) :: two_rdm(n_act_orb, n_act_orb, n_act_orb, n_act_orb)
 integer :: i,j,k,l,iorb,jorb,korb,lorb,istate
 BEGIN_DOC
! This routine computes the two electron repulsion within the active space using various providers 
! 
 END_DOC

 double precision :: vijkl,get_two_e_integral
 double precision :: wee_ab(N_states),rdmab

 wee_ab  = 0.d0

 iorb = 1
 jorb = 1
 korb = 1
 lorb = 1
 vijkl = get_two_e_integral(lorb,korb,jorb,iorb,mo_integrals_map)                                 
! call diagonalize_ci
 print*,'**************************'
 print*,'**************************'
 do istate = 1, N_states
   !! PURE ACTIVE PART 
   !! 
   do i = 1, n_act_orb
    iorb = list_act(i)
    do j = 1, n_act_orb
     jorb = list_act(j)
     do k = 1, n_act_orb
      korb = list_act(k)
      do l = 1, n_act_orb
       lorb = list_act(l)

       vijkl = get_two_e_integral(lorb,korb,jorb,iorb,mo_integrals_map)                                 


       rdmab   =  two_rdm(l,k,j,i)
!      if(dabs(rdmab).gt.1.d-10)then
!        print*,i,j,k,l
!        print*,iorb,jorb,korb,lorb
!        print*,rdmab
!        print*,vijkl
!      endif


       wee_ab(istate)     += vijkl * rdmab

      enddo
     enddo
    enddo
   enddo
   print*,''
   print*,''
   print*,'Active space only energy with usual 2RDM'
   print*,'wee_ab(istate)          = ',wee_ab(istate)
   print*,''
  enddo
end

subroutine routine_full_mos_test
 implicit none
 integer :: i,j,k,l,iorb,jorb,korb,lorb,istate
 BEGIN_DOC
! This routine computes the two electron repulsion using various providers 
! 
 END_DOC

 double precision :: vijkl,rdmaa,get_two_e_integral,rdmab,rdmbb,rdmtot
 double precision :: wee_ab(N_states),wee_tot(N_states)
 double precision :: wee_ab_st_av, rdm_ab_st_av
 double precision :: wee_tot_st_av, rdm_tot_st_av
 double precision :: ab_norm(N_states),tot_norm(N_states)

 ab_norm  = 0.d0
 tot_norm = 0.d0

 wee_ab = 0.d0
 wee_tot = 0.d0

 iorb = 1
 jorb = 1
 korb = 1
 lorb = 1
 vijkl = get_two_e_integral(lorb,korb,jorb,iorb,mo_integrals_map)                                 
 provide full_occ_2_rdm_ab_mo  full_occ_2_rdm_aa_mo  full_occ_2_rdm_bb_mo  full_occ_2_rdm_spin_trace_mo 
 print*,'**************************'
 print*,'**************************'
 do istate = 1, N_states
   do i = 1, n_core_inact_act_orb
    iorb = list_core_inact_act(i)
    do j = 1, n_core_inact_act_orb
     jorb = list_core_inact_act(j)
     do k = 1, n_core_inact_act_orb
      korb = list_core_inact_act(k)
      do l = 1, n_core_inact_act_orb
       lorb = list_core_inact_act(l)
       vijkl = get_two_e_integral(lorb,korb,jorb,iorb,mo_integrals_map)                                 

       rdmab  =  full_occ_2_rdm_ab_mo(l,k,j,i,istate)
!      rdmtot =  full_occ_2_rdm_spin_trace_mo(l,k,j,i,istate)

       wee_ab(istate) += vijkl * rdmab
!      wee_tot(istate)+= vijkl * rdmtot
      enddo
     enddo
     ab_norm(istate) += full_occ_2_rdm_ab_mo(j,i,j,i,istate)
!    tot_norm(istate)+= full_occ_2_rdm_spin_trace_mo(j,i,j,i,istate)
    enddo
   enddo
   print*,''
   print*,''
   print*,'Full energy for state ',istate
   print*,'wee_ab(istate)          = ',wee_ab(istate)
   print*,''
!  print*,'wee_tot(istate)         = ',wee_tot(istate)
   print*,''
   print*,'Normalization of two-rdms '
   print*,''
   print*,'ab_norm(istate)         = ',ab_norm(istate)
   print*,'N_alpha * N_beta        = ',elec_num_tab(1) * elec_num_tab(2)
   print*,''
!  print*,'tot_norm(istate)        = ',tot_norm(istate)
!  print*,'N(N-1)/2                = ',elec_num*(elec_num - 1)/2
  enddo

end
