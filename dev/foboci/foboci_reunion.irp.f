program osoci_program
implicit none
   do_it_perturbative = .True.
   touch do_it_perturbative
   call FOBOCI_lmct_mlct_old_thr
   call provide_all_the_rest
end
subroutine provide_all_the_rest
implicit none
integer :: i
   call update_one_e_dm_mo
   call set_lmct_mlct_to_psi_det
   call diagonalize_CI
   call save_wavefunction


end

subroutine set_lmct_mlct_to_psi_det
 implicit none
 integer :: i,j,m,n,i_hole,i_particle
 integer :: hole_particle(1000,2), n_couples
 integer(bit_kind) :: key_tmp(N_int,2)
 integer :: N_det_total,i_ok

 call collect_lmct_mlct(hole_particle,n_couples)
 call set_psi_det_to_generators_restart
 N_det_total = N_det_generators_restart 
 do i = 1, n_couples
  i_hole     = hole_particle(i,1)
  i_particle = hole_particle(i,2)
  do m = 1, N_det_generators_restart
   do n = 1, N_int
    key_tmp(n,1) = psi_det_generators_restart(n,1,m)
    key_tmp(n,2) = psi_det_generators_restart(n,2,m)
   enddo
   ! You excite the beta electron from i_hole to i_particle
   call do_single_excitation(key_tmp,i_hole,i_particle,2,i_ok)
   if(i_ok==1)then
    N_det_total +=1
    do n = 1, N_int
     psi_det(n,1,N_det_total) = key_tmp(n,1)
     psi_det(n,2,N_det_total) = key_tmp(n,2)
    enddo
   endif

   do n = 1, N_int
    key_tmp(n,1) = psi_det_generators_restart(n,1,m)
    key_tmp(n,2) = psi_det_generators_restart(n,2,m)
   enddo

   ! You excite the alpha electron from i_hole to i_particle
   call do_single_excitation(key_tmp,i_hole,i_particle,1,i_ok)
   if(i_ok==1)then
    N_det_total +=1
    do n = 1, N_int
     psi_det(n,1,N_det_total) = key_tmp(n,1)
     psi_det(n,2,N_det_total) = key_tmp(n,2)
    enddo
   endif
  enddo
 enddo
 
 N_det = N_det_total
 integer :: k
 do k = 1, N_states
  do i = 1, N_det
   psi_coef(i,k) = 1.d0/dsqrt(dble(N_det_total))
  enddo
 enddo
 SOFT_TOUCH N_det psi_det psi_coef
 logical :: found_duplicates
 call remove_duplicates_in_psi_det(found_duplicates)
end

subroutine collect_lmct_mlct(hole_particle,n_couples)
 implicit none
 integer, intent(out) :: hole_particle(1000,2), n_couples
 BEGIN_DOC
 ! Collect all the couple holes/particles of the important LMCT
 ! hole_particle(i,1) = ith hole
 ! hole_particle(i,2) = ith particle
 ! n_couples is the number of important excitations
 END_DOC
 double precision, allocatable :: tmp(:,:)
 print*,'COLLECTING THE PERTINENT LMCT (1h)'
 print*,'AND THE PERTINENT MLCT (1p)'
 allocate(tmp(size(one_e_dm_mo_alpha_osoci,1),size(one_e_dm_mo_alpha_osoci,2)))
 tmp =  one_e_dm_mo_alpha_osoci + one_e_dm_mo_beta_osoci
 integer :: i,j,iorb,jorb
 n_couples = 0
 do i = 1,n_act_orb
  iorb = list_act(i)
  do j = 1, n_inact_orb
   jorb = list_inact(j)
   if(dabs(tmp(iorb,jorb)).gt.threshold_lmct)then
    n_couples +=1
    hole_particle(n_couples,1) = jorb
    hole_particle(n_couples,2) = iorb
    print*,'DM'
    print*,hole_particle(n_couples,1),hole_particle(n_couples,2),tmp(iorb,jorb)
   endif
  enddo
  do j = 1, n_virt_orb
   jorb = list_virt(j)
   if(dabs(tmp(iorb,jorb)).gt.threshold_mlct)then
    n_couples +=1
    hole_particle(n_couples,1) = iorb
    hole_particle(n_couples,2) = jorb
    print*,'DM'
    print*,hole_particle(n_couples,1),hole_particle(n_couples,2),tmp(iorb,jorb)
   endif
  enddo
 enddo
 deallocate(tmp)
 print*,'number of meaning full couples of holes/particles '
 print*,'n_couples = ',n_couples


end


