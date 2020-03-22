 BEGIN_PROVIDER [ double precision , integral_on_top , (N_states) ]
 implicit none
 BEGIN_DOC
 ! Integral in real space of the on-top pair density
 ! 
 ! !!!!! WARNING !!!!! if no_core_density = .True. then all contributions from the core orbitals will be removed
 END_DOC
 integer :: i,j,k,l,istate,iorb,jorb,korb,lorb
 double precision :: get_mo_bielec_integral_ijkl_r3
 double precision, allocatable :: integrals_ij(:,:)
 allocate(integrals_ij(mo_num,mo_num))
 double precision :: wall0,wall1
 integral_on_top = 0.d0
 provide full_occ_2_rdm_ab_mo mo_integrals_ijkl_r3_map
 call wall_time(wall0)
  do istate = 1, N_states
   do i = 1, n_core_inact_act_orb
    iorb = list_core_inact_act(i)
    do j = 1, n_core_inact_act_orb
     jorb = list_core_inact_act(j)
     call get_mo_bielec_integrals_ijkl_r3_ij(iorb,jorb,mo_num,integrals_ij,mo_integrals_ijkl_r3_map)
     do l = 1, n_core_inact_act_orb
      lorb = list_core_inact_act(l)
      do k = 1, n_core_inact_act_orb
       korb = list_core_inact_act(k)
       integral_on_top(istate) += full_occ_2_rdm_ab_mo(k,l,i,j,istate) * integrals_ij(korb,lorb)
       enddo
      enddo
     enddo
    enddo
   enddo
 deallocate(integrals_ij)
 call wall_time(wall1)
 print*,'Time to provide integral_on_top   = ',wall1-wall0
 END_PROVIDER

