 BEGIN_PROVIDER [double precision, f_psi_ii_ab, (n_points_final_grid)]
&BEGIN_PROVIDER [double precision, on_top_psi_ii, (n_points_final_grid)]
 implicit none
 BEGIN_DOC
!
! Function f_{\Psi^B}(r,r) of Eq. (22) of J. Chem. Phys. 149, 194301 (2018) on each point of the grid for a HF wave function
!
 END_DOC 
 integer :: ipoint
 double precision :: wall0,wall1,r(3),f_ii_val_ab,two_bod_dens
 f_psi_ii_ab = 0.d0
 r = 0.d0
 ! To initialize parallelization
 call f_ii_valence_ab(r,r,f_ii_val_ab,two_bod_dens)

 call wall_time(wall0)
 !$OMP PARALLEL        &
 !$OMP DEFAULT (NONE)  &
 !$OMP PRIVATE (ipoint,r,f_ii_val_ab,two_bod_dens) & 
 !$OMP SHARED  (n_points_final_grid,f_psi_ii_ab,on_top_psi_ii,final_grid_points)
 !$OMP DO              
  do ipoint = 1, n_points_final_grid
   r(1)   = final_grid_points(1,ipoint)
   r(2)   = final_grid_points(2,ipoint)
   r(3)   = final_grid_points(3,ipoint)
   call f_ii_valence_ab(r,r,f_ii_val_ab,two_bod_dens)
   f_psi_ii_ab(ipoint)   = f_ii_val_ab
   on_top_psi_ii(ipoint) = two_bod_dens
  enddo
 !$OMP END DO
 !$OMP END PARALLEL
 call wall_time(wall1)
 print*,'Time to provide f_psi_ii_ab = ',wall1 - wall0

END_PROVIDER 
