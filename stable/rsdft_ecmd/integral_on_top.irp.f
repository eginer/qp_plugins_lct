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


subroutine routine_test_cas_based_on_top_density
 implicit none
 BEGIN_DOC
! This routine computes the integral in real-space of the on-top pair density in two ways:
!
! +) using the PROVIDER "total_cas_on_top_density" 
!
! +) using the ROUTINE  "give_on_top_in_r_one_state" 
!
! +) using the PROVIDER "integral_on_top" which is the analytical integral 
 END_DOC
 integer :: i_point,istate
 double precision :: r(3),prov_dm,on_top_in_r
 double precision :: accu(N_states),accu_2(N_states)
 accu = 0.d0
 accu_2 = 0.d0
 do i_point = 1, n_points_final_grid
  r(1) = final_grid_points(1,i_point)
  r(2) = final_grid_points(2,i_point)
  r(3) = final_grid_points(3,i_point)
  do istate = 1, N_states
   ! provider to get the on-top on the becke grid 
   prov_dm = total_cas_on_top_density(i_point,istate)
   ! subroutine to get the on-top in any points in space
   call give_on_top_in_r_one_state(r,istate,on_top_in_r)
   accu(istate) += dabs(prov_dm - on_top_in_r) * final_weight_at_r_vector(i_point)
   accu_2(istate) += prov_dm * final_weight_at_r_vector(i_point)
  enddo
 enddo
 print*,'difference between provider and routine = ',accu
 print*,'integral of the on-top                  = ',accu_2
 print*,'integral_on_top                         = ',integral_on_top(:)
end

