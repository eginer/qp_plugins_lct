
BEGIN_PROVIDER [double precision, f_psi_cas_ab, (n_points_final_grid,N_states)]
 implicit none
!
! Function f_{\Psi^B}(r,r) of Eq. (22) of J. Chem. Phys. 149, 194301 (2018) on each point of the grid and for all states
! 
! Assumes that the wave function in psi_det is developped within an active space defined
 integer :: ipoint,k,l,istate 
 double precision :: wall0,wall1
 provide full_occ_2_rdm_ab_mo 
 print*,'Providing  f_psi_cas_ab ..... '
 call wall_time(wall0)
 provide core_inact_act_V_kl_contracted full_occ_2_rdm_cntrctd
 !$OMP PARALLEL        &
 !$OMP DEFAULT (NONE)  &
 !$OMP PRIVATE (ipoint,k,l,istate) & 
 !$OMP SHARED  (n_core_inact_act_orb, n_points_final_grid, full_occ_2_rdm_cntrctd, core_inact_act_V_kl_contracted, f_psi_cas_ab,N_states)
 !$OMP DO              
 do istate = 1, N_states
  do ipoint = 1, n_points_final_grid
   f_psi_cas_ab(ipoint,istate) = 0.d0
   do l = 1, n_core_inact_act_orb ! 2 
    do k = 1, n_core_inact_act_orb ! 1
     f_psi_cas_ab(ipoint,istate) += core_inact_act_V_kl_contracted(k,l,ipoint) * full_occ_2_rdm_cntrctd(k,l,ipoint,istate)
    enddo
   enddo
  enddo
 enddo
 !$OMP END DO
 !$OMP END PARALLEL
 call wall_time(wall1)
 print*,'Time to provide f_psi_cas_ab = ',wall1 - wall0
END_PROVIDER 


 BEGIN_PROVIDER [double precision, f_psi_hf_ab, (n_points_final_grid)]
&BEGIN_PROVIDER [double precision, on_top_psi_hf, (n_points_final_grid)]
 implicit none
 BEGIN_DOC
!
! Function f_{\Psi^B}(r,r) of Eq. (22) of J. Chem. Phys. 149, 194301 (2018) on each point of the grid for a HF wave function
!
 END_DOC 
 integer :: ipoint
 double precision :: wall0,wall1,r(3),f_HF_val_ab,two_bod_dens
 f_psi_hf_ab = 0.d0
 r = 0.d0
 ! To initialize parallelization
 call f_HF_valence_ab(r,r,f_HF_val_ab,two_bod_dens)

 call wall_time(wall0)
 !$OMP PARALLEL        &
 !$OMP DEFAULT (NONE)  &
 !$OMP PRIVATE (ipoint,r,f_HF_val_ab,two_bod_dens) & 
 !$OMP SHARED  (n_points_final_grid,f_psi_hf_ab,on_top_psi_hf,final_grid_points)
 !$OMP DO              
  do ipoint = 1, n_points_final_grid
   r(1)   = final_grid_points(1,ipoint)
   r(2)   = final_grid_points(2,ipoint)
   r(3)   = final_grid_points(3,ipoint)
   call f_HF_valence_ab(r,r,f_HF_val_ab,two_bod_dens)
   f_psi_hf_ab(ipoint)   = f_HF_val_ab
   on_top_psi_hf(ipoint) = two_bod_dens
  enddo
 !$OMP END DO
 !$OMP END PARALLEL
 call wall_time(wall1)
 print*,'Time to provide f_psi_hf_ab = ',wall1 - wall0

END_PROVIDER 
