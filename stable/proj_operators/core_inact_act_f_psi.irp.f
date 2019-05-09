
BEGIN_PROVIDER [double precision, core_inact_act_rho2_kl_contracted_transposed, (n_points_final_grid,n_core_inact_act_orb,n_core_inact_act_orb)]
 implicit none
 BEGIN_DOC
! core_inact_act_rho2_kl_contracted_transposed(k,l,ipoint) = \sum_{ij} core_inact_act_rho2_{ij}^{kl} phi_i(r_ipoint) phi_j(r_ipoint)
 END_DOC
 integer :: i,j,k,l
 integer :: ipoint
 double precision, allocatable :: mos_array_r(:),r(:)
 provide core_inact_act_two_bod_alpha_beta_mo_physicist
 double precision :: wall0,wall1
 print*,'Providing  core_inact_act_rho2_kl_contracted_transposed ..... '
 call wall_time(wall0)
 core_inact_act_rho2_kl_contracted_transposed = 0.d0
 !$OMP PARALLEL        &
 !$OMP DEFAULT (NONE)  &
 !$OMP PRIVATE (ipoint,k,l,i,j) & 
 !$OMP SHARED  (n_core_inact_act_orb, n_points_final_grid, core_inact_act_rho2_kl_contracted_transposed, final_grid_points,core_inact_act_two_bod_alpha_beta_mo_physicist,mos_in_r_array )
 !$OMP DO              
 do l = 1, n_core_inact_act_orb ! 2 
  do k = 1, n_core_inact_act_orb ! 1 
   do ipoint = 1, n_points_final_grid
    do j = 1, n_core_inact_act_orb
     do i = 1, n_core_inact_act_orb
                                          !1 2                                     1 2 1 2 
      core_inact_act_rho2_kl_contracted_transposed(ipoint,k,l) += core_inact_act_two_bod_alpha_beta_mo_physicist(i,j,k,l,1) * mos_in_r_array(j,ipoint) * mos_in_r_array(i,ipoint)
     enddo
    enddo
   enddo
  enddo
 enddo
 !$OMP END DO
 !$OMP END PARALLEL
 call wall_time(wall1)
 print*,'Time to provide core_inact_act_rho2_kl_contracted_transposed = ',wall1 - wall0

END_PROVIDER 

BEGIN_PROVIDER [double precision, core_inact_act_rho2_kl_contracted, (n_core_inact_act_orb,n_core_inact_act_orb,n_points_final_grid)]
 implicit none
 BEGIN_DOC
! core_inact_act_rho2_kl_contracted(k,l,ipoint) = \sum_{ij} V_{ij}^{kl} phi_i(r_ipoint) phi_j(r_ipoint)
 END_DOC
 integer :: ipoint,k,l
 do k = 1, n_core_inact_act_orb
  do l = 1, n_core_inact_act_orb
   do ipoint = 1, n_points_final_grid
    core_inact_act_rho2_kl_contracted(k,l,ipoint) = core_inact_act_rho2_kl_contracted_transposed(ipoint,k,l)
   enddo
  enddo
 enddo
 free core_inact_act_rho2_kl_contracted_transposed 
END_PROVIDER 




BEGIN_PROVIDER [double precision, core_inact_act_V_kl_contracted_transposed, (n_points_final_grid,n_core_inact_act_orb,n_core_inact_act_orb)]
 implicit none
 BEGIN_DOC
! core_inact_act_V_kl_contracted_transposed(ipoint,k,l) = \sum_{ij} V_{ij}^{kl} phi_i(r_ipoint) phi_j(r_ipoint)
 END_DOC
 integer :: i,j,k,l,kk,ll
 integer :: ipoint
 double precision, allocatable :: integrals_array(:,:), mos_array_r(:),r(:)
 ! just not to mess with parallelization
 allocate(integrals_array(mo_num,mo_num))
  k = 1
  l = 1
  call get_mo_two_e_integrals_ij(k,l,mo_num,integrals_array,mo_integrals_map) 
 deallocate(integrals_array)
 double precision :: wall0,wall1
 call wall_time(wall0)

 core_inact_act_V_kl_contracted_transposed = 0.d0
 print*,'Providing  core_inact_act_V_kl_contracted_transposed ..... '
 !$OMP PARALLEL        &
 !$OMP DEFAULT (NONE)  &
 !$OMP PRIVATE (ipoint,kk,ll,k,l,i,j,integrals_array) & 
 !$OMP SHARED (mo_num, n_points_final_grid, core_inact_act_V_kl_contracted_transposed, mo_integrals_map,final_grid_points,mos_in_r_array, n_core_inact_act_orb, list_core_inact_act,n_orb_max_basis)
 allocate(integrals_array(mo_num,mo_num))
 !$OMP DO              
  do l = 1, n_core_inact_act_orb! 2 
   ll = list_core_inact_act(l)
   do k = 1, n_core_inact_act_orb ! 1 
    kk = list_core_inact_act(k)
    call get_mo_two_e_integrals_ij(kk,ll,mo_num,integrals_array,mo_integrals_map) 
    do ipoint = 1, n_points_final_grid
     do j = 1, n_orb_max_basis ! condition on mo_num in order to ensure the correct CBS limit 
      do i = 1, n_orb_max_basis ! 
                                        !1 2                     1 2 
       core_inact_act_V_kl_contracted_transposed(ipoint,k,l) += integrals_array(i,j) * mos_in_r_array(j,ipoint) * mos_in_r_array(i,ipoint)
      enddo
     enddo
    enddo
   enddo
  enddo
 !$OMP END DO
 deallocate(integrals_array)
 !$OMP END PARALLEL

 call wall_time(wall1)
 print*,'Time to provide core_inact_act_V_kl_contracted_transposed = ',wall1 - wall0
END_PROVIDER 

BEGIN_PROVIDER [double precision, core_inact_act_V_kl_contracted, (n_core_inact_act_orb,n_core_inact_act_orb,n_points_final_grid)]
 implicit none
 BEGIN_DOC
! core_inact_act_V_kl_contracted(k,l,ipoint) = \sum_{ij} V_{ij}^{kl} phi_i(r_ipoint) phi_j(r_ipoint)
 END_DOC
 integer :: ipoint,k,l
 do k = 1, n_core_inact_act_orb
  do l = 1, n_core_inact_act_orb
   do ipoint = 1, n_points_final_grid
    core_inact_act_V_kl_contracted(k,l,ipoint) = core_inact_act_V_kl_contracted_transposed(ipoint,k,l)
   enddo
  enddo
 enddo
 free core_inact_act_V_kl_contracted_transposed 

END_PROVIDER 


BEGIN_PROVIDER [double precision, core_inact_act_f_psi_ab, (n_points_final_grid)]
 implicit none
 integer :: ipoint
 integer :: k,l 
 double precision :: wall0,wall1
 provide core_inact_act_two_bod_alpha_beta_mo_physicist 
 print*,'Providing  core_inact_act_f_psi_ab ..... '
 call wall_time(wall0)
 provide core_inact_act_V_kl_contracted core_inact_act_rho2_kl_contracted
 !$OMP PARALLEL        &
 !$OMP DEFAULT (NONE)  &
 !$OMP PRIVATE (ipoint,k,l) & 
 !$OMP SHARED  (n_core_inact_act_orb, n_points_final_grid, core_inact_act_rho2_kl_contracted, core_inact_act_V_kl_contracted, core_inact_act_f_psi_ab)
 !$OMP DO              
 do ipoint = 1, n_points_final_grid
  core_inact_act_f_psi_ab(ipoint) = 0.d0
  do l = 1, n_core_inact_act_orb ! 2 
   do k = 1, n_core_inact_act_orb ! 1
    core_inact_act_f_psi_ab(ipoint) += core_inact_act_V_kl_contracted(k,l,ipoint) * core_inact_act_rho2_kl_contracted(k,l,ipoint)
   enddo
  enddo
 enddo
 !$OMP END DO
 !$OMP END PARALLEL
 call wall_time(wall1)
 print*,'Time to provide core_inact_act_f_psi_ab = ',wall1 - wall0
END_PROVIDER 

