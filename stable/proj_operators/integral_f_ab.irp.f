BEGIN_PROVIDER [double precision, V_kl_contracted_r6, (n_points_final_grid,n_points_final_grid,mo_num,mo_num)]
 implicit none
 BEGIN_DOC
! V_kl_contracted(ipoint,k,l) = \sum_{ij} V_{ij}^{kl} phi_i(r_ipoint) phi_j(r_jpoint)
 END_DOC
 integer :: i,j,k,l
 integer :: ipoint,jpoint
 double precision, allocatable :: integrals_array(:,:), mos_array_r(:),r(:)
 ! just not to mess with parallelization
 allocate(integrals_array(mo_num,mo_num))
  k = 1
  l = 1
  call get_mo_two_e_integrals_ij(k,l,mo_num,integrals_array,mo_integrals_map) 
 provide mos_in_r_array
 deallocate(integrals_array)
 double precision :: wall0,wall1
 call wall_time(wall0)

 print*,'Providing  V_kl_contracted ..... '
 V_kl_contracted_r6(:,:,:,:) = 0.d0
 !$OMP PARALLEL        &
 !$OMP DEFAULT (NONE)  &
 !$OMP PRIVATE (ipoint,jpoint,k,l,i,j,integrals_array) & 
 !$OMP SHARED (mo_num, n_points_final_grid, V_kl_contracted_r6, mo_integrals_map,final_grid_points,mos_in_r_array)
 allocate(integrals_array(mo_num,mo_num))
 !$OMP DO              
  do l = 1, mo_num ! 2 
   do k = 1, mo_num ! 1 
    call get_mo_two_e_integrals_ij(k,l,mo_num,integrals_array,mo_integrals_map) 
    do jpoint = 1, n_points_final_grid
     do ipoint = 1, n_points_final_grid
      do j = 1, mo_num
       do i = 1, mo_num
                           !1       2   !1 2                     1 2                     2                          1
        V_kl_contracted_r6(ipoint,jpoint,k,l) += integrals_array(i,j) * mos_in_r_array(j,jpoint) * mos_in_r_array(i,ipoint)
       enddo
      enddo
     enddo
    enddo
   enddo
  enddo
 !$OMP END DO
 deallocate(integrals_array)
 !$OMP END PARALLEL

 call wall_time(wall1)
 print*,'Time to provide V_kl_contracted_r6 = ',wall1 - wall0
END_PROVIDER 


BEGIN_PROVIDER [double precision, rho2_kl_contracted_r6, (n_points_final_grid,n_points_final_grid,mo_num,mo_num)]
 implicit none
 BEGIN_DOC
! rho2_kl_contracted_r6(k,l,ipoint) = \sum_{ij} rho2_{ij}^{kl} phi_i(r_ipoint) phi_j(r_ipoint)
 END_DOC
 integer :: i,j,k,l
 integer :: ipoint,jpoint
 double precision, allocatable :: mos_array_r(:),r(:)
 provide two_bod_alpha_beta_mo_physicist
 provide mos_in_r_array
 double precision :: wall0,wall1
 print*,'Providing  rho2_kl_contracted_r6 ..... '
 call wall_time(wall0)
 rho2_kl_contracted_r6(:,:,:,:) = 0.d0
 !$OMP PARALLEL        &
 !$OMP DEFAULT (NONE)  &
 !$OMP PRIVATE (ipoint,jpoint,k,l,i,j) & 
 !$OMP SHARED  (mo_num, n_points_final_grid, rho2_kl_contracted_r6, final_grid_points,two_bod_alpha_beta_mo_physicist,mos_in_r_array )
 !$OMP DO              
 do l = 1, mo_num ! 2 
  do k = 1, mo_num ! 1 
   do jpoint = 1, n_points_final_grid
    do ipoint = 1, n_points_final_grid
     do j = 1, mo_num
      do i = 1, mo_num
                            ! 1       2    1 2                                     1 2 1 2                      2                         1
       rho2_kl_contracted_r6(ipoint,jpoint,k,l) += two_bod_alpha_beta_mo_physicist(i,j,k,l,1) * mos_in_r_array(j,jpoint) * mos_in_r_array(i,ipoint)
      enddo
     enddo
    enddo
   enddo
  enddo
 enddo
 !$OMP END DO
 !$OMP END PARALLEL
 call wall_time(wall1)
 print*,'Time to provide rho2_kl_contracted_r6 = ',wall1 - wall0

END_PROVIDER 


BEGIN_PROVIDER [double precision, integral_r1r2_f_psi_ab]
 implicit none
 integer :: ipoint,jpoint
 integer :: k,l 
 double precision :: wall0,wall1,weight1,weight2
 provide two_bod_alpha_beta_mo_physicist 
 print*,'Providing  f_psi_ab ..... '
 call wall_time(wall0)
 provide V_kl_contracted_r6 rho2_kl_contracted_r6
 integral_r1r2_f_psi_ab = 0.d0
 !$OMP PARALLEL        &
 !$OMP DEFAULT (NONE)  &
 !$OMP PRIVATE (ipoint,jpoint,k,l,weight1,weight2) & 
 !$OMP SHARED  (mo_num, n_points_final_grid, rho2_kl_contracted_r6, V_kl_contracted_r6, final_weight_at_r_vector,integral_r1r2_f_psi_ab)
 !$OMP DO              
 do l = 1, mo_num ! 2 
  do k = 1, mo_num ! 1
   do ipoint = 1, n_points_final_grid
    weight1=final_weight_at_r_vector(ipoint)
    do jpoint = 1, n_points_final_grid
     weight2=final_weight_at_r_vector(jpoint)
     integral_r1r2_f_psi_ab += weight1 * weight2 * V_kl_contracted_r6(jpoint,ipoint,k,l) * rho2_kl_contracted_r6(jpoint,ipoint,k,l)
    enddo
   enddo
  enddo
 enddo
 !$OMP END DO
 !$OMP END PARALLEL
 call wall_time(wall1)
 print*,'Time to provide integral_r1r2_f_psi_ab = ',wall1 - wall0
END_PROVIDER 

BEGIN_PROVIDER [double precision, integral_r1r2_f_HF_aa]
 implicit none
 integer :: ipoint,jpoint
 integer :: k,l 
 double precision :: wall0,wall1,weight1,weight2,r1(3),r2(3),f_HF_aa
 provide two_bod_alpha_beta_mo_physicist 
 print*,'Providing  integral_r1r2_f_HF_aa ..... '
 call wall_time(wall0)
 integral_r1r2_f_HF_aa = 0.d0
 !$OMP PARALLEL        &
 !$OMP DEFAULT (NONE)  &
 !$OMP PRIVATE (ipoint,jpoint,k,l,weight1,weight2,r1,r2) & 
 !$OMP SHARED  (mo_num, n_points_final_grid, final_weight_at_r_vector,integral_r1r2_f_HF_aa,final_grid_points)
 !$OMP DO              
 do ipoint = 1, n_points_final_grid
  weight1=final_weight_at_r_vector(ipoint)
  r1(1) = final_grid_points(1,ipoint)
  r1(2) = final_grid_points(2,ipoint)
  r1(3) = final_grid_points(3,ipoint)
  do jpoint = 1, n_points_final_grid
   r2(1) = final_grid_points(1,jpoint)
   r2(2) = final_grid_points(2,jpoint)
   r2(3) = final_grid_points(3,jpoint)
   weight2=final_weight_at_r_vector(jpoint)
   integral_r1r2_f_HF_aa += weight1 * weight2 * f_HF_aa(r1,r2)
  enddo
 enddo
 !$OMP END DO
 !$OMP END PARALLEL
 call wall_time(wall1)
 print*,'Time to provide integral_r1r2_f_HF_aa = ',wall1 - wall0
END_PROVIDER 

BEGIN_PROVIDER [double precision, integral_r1r2_f_HF_ab]
 implicit none
 integer :: ipoint,jpoint
 integer :: k,l 
 double precision :: wall0,wall1,weight1,weight2,r1(3),r2(3),f_HF_aa,integral_psi,two_bod
 provide two_bod_alpha_beta_mo_physicist 
 print*,'Providing  integral_r1r2_f_HF_ab ..... '
 call wall_time(wall0)
 integral_r1r2_f_HF_ab = 0.d0
 !$OMP PARALLEL        &
 !$OMP DEFAULT (NONE)  &
 !$OMP PRIVATE (ipoint,jpoint,k,l,weight1,weight2,r1,r2,integral_psi,two_bod) & 
 !$OMP SHARED  (mo_num, n_points_final_grid, final_weight_at_r_vector,integral_r1r2_f_HF_ab,final_grid_points)
 !$OMP DO              
 do ipoint = 1, n_points_final_grid
  weight1=final_weight_at_r_vector(ipoint)
  r1(1) = final_grid_points(1,ipoint)
  r1(2) = final_grid_points(2,ipoint)
  r1(3) = final_grid_points(3,ipoint)
  do jpoint = 1, n_points_final_grid
   r2(1) = final_grid_points(1,jpoint)
   r2(2) = final_grid_points(2,jpoint)
   r2(3) = final_grid_points(3,jpoint)
   weight2=final_weight_at_r_vector(jpoint)
   call f_HF_ab(r1,r2,integral_psi,two_bod)
   integral_r1r2_f_HF_ab += weight1 * weight2 * integral_psi
  enddo
 enddo
 !$OMP END DO
 !$OMP END PARALLEL
 call wall_time(wall1)
 print*,'Time to provide integral_r1r2_f_HF_ab = ',wall1 - wall0
END_PROVIDER 

