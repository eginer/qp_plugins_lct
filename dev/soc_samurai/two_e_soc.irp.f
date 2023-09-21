BEGIN_PROVIDER [ double precision, v_c_ij_grid, ( ao_num, ao_num,n_points_final_grid)]
 implicit none
 BEGIN_DOC
! int dr phi_i(r) phi_j(r) / |r - R|
 END_DOC
 integer :: i,j,ipoint
 double precision :: r(3),NAI_pol_mult_erf_ao
 double precision :: int_coulomb
 provide final_grid_points 
 double precision :: wall0, wall1
 call wall_time(wall0)
 !$OMP PARALLEL                  &
 !$OMP DEFAULT (NONE)            &
 !$OMP PRIVATE (i,j,ipoint,r,int_coulomb) & 
 !$OMP SHARED (ao_num,n_points_final_grid,v_c_ij_grid,final_grid_points)
 !$OMP DO SCHEDULE (dynamic)
 do ipoint = 1, n_points_final_grid
   do i = 1, ao_num
    do j = i, ao_num
     r(1) = final_grid_points(1,ipoint)
     r(2) = final_grid_points(2,ipoint)
     r(3) = final_grid_points(3,ipoint)
     int_coulomb = NAI_pol_mult_erf_ao(i,j,1.d+9,r)
     v_c_ij_grid(j,i,ipoint)= ( int_coulomb )
   enddo
  enddo
 enddo
 !$OMP END DO
 !$OMP END PARALLEL

 do ipoint = 1, n_points_final_grid
  do i = 1, ao_num
   do j = 1, i-1
     v_c_ij_grid(j,i,ipoint)= v_c_ij_grid(i,j,ipoint)
    enddo
   enddo
  enddo

 call wall_time(wall1)
 print*,'wall time for v_c_ij_grid  ',wall1 - wall0
END_PROVIDER

 BEGIN_PROVIDER[double precision, v_c_ij_grid_transp, (n_points_final_grid,ao_num, ao_num)]
 implicit none
 BEGIN_DOC
! Transposed v_c_ij_grid 
 END_DOC
 integer :: i,j,ipoint
 do i = 1, ao_num
  do j = 1, ao_num
   do ipoint = 1,n_points_final_grid 
    v_c_ij_grid_transp(ipoint,j,i) = v_c_ij_grid(j,i,ipoint)
   enddo
  enddo
 enddo
 END_PROVIDER

BEGIN_PROVIDER [ double precision, rdm_one_e_soc, ( ao_num, ao_num,n_states)]
 implicit none
 BEGIN_DOC
! One electron (reduced??) density matrix for SOC
!!!!!!!!! WARNING!!!!!!!!!!!!!!
! It does not depend on istate 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 END_DOC
 integer :: i,j,istate
 do istate = 1, n_states
  do i = 1, ao_num
   do j = 1, ao_num 
    rdm_one_e_soc(j,i,istate) =   one_e_dm_ao_alpha(j,i) + one_e_dm_ao_beta(j,i) 
   enddo
  enddo
 enddo

END_PROVIDER


BEGIN_PROVIDER [ double precision, D_mu_nu_find_a_better_name, (3,ao_num, ao_num,n_states)]
 implicit none
 BEGIN_DOC
 !D_mu_nu_find_a_better_name(ao1,ao2,nstate) = 
 END_DOC
 integer :: i,j,ipoint,istate,k,l,coord
 double precision :: grid_sum(3)
 D_mu_nu_find_a_better_name = 0.d0
 do istate = 1, n_states
  do i= 1, ao_num
   do j = 1, ao_num
    do k = 1, ao_num
     do l = 1, ao_num

      grid_sum = 0d0
      do ipoint = 1,n_points_final_grid 
       do coord = 1,3 
         grid_sum(coord) += vec_prod_grad_ao_at_r(coord,ipoint,j,i) * v_c_ij_grid_transp(ipoint,l,k) * final_weight_at_r_vector(ipoint)
       enddo
      enddo

      do coord = 1,3 
       D_mu_nu_find_a_better_name(coord,j,i,istate) += rdm_one_e_soc(l,k,istate) * grid_sum(coord) 
      enddo

     enddo
    enddo
   enddo
  enddo
 enddo
END_PROVIDER

BEGIN_PROVIDER [ double precision, E1_mu_nu_find_a_better_name, (3,ao_num, ao_num,n_states)]
 implicit none
 BEGIN_DOC
 !E1_mu_nu_find_a_better_name(ao1,ao2,nstate) = 
 END_DOC
 integer :: i,j,ipoint,istate,k,l,coord
 double precision :: grid_sum(3)
 E1_mu_nu_find_a_better_name = 0.d0
 do istate = 1, n_states
  do i= 1, ao_num
   do j = 1, ao_num
    do k = 1, ao_num
     do l = 1, ao_num

      grid_sum = 0d0
      do ipoint = 1, n_points_final_grid 
       do coord = 1,3 
        grid_sum(coord) += vec_prod_grad_ao_at_r(coord,ipoint,k,i) * v_c_ij_grid_transp(ipoint,l,j) * final_weight_at_r_vector_extra(ipoint)
       enddo
      enddo

      do coord = 1,3 
       E1_mu_nu_find_a_better_name(coord,j,i,istate) += rdm_one_e_soc(l,k,istate) * grid_sum(coord) 
      enddo

     enddo
    enddo
   enddo
  enddo
 enddo
 E1_mu_nu_find_a_better_name = -3.d0/2.d0 * E1_mu_nu_find_a_better_name 
END_PROVIDER

BEGIN_PROVIDER [ double precision, E2_mu_nu_find_a_better_name, (3,ao_num, ao_num,n_states)]
 implicit none
 BEGIN_DOC
 !E1_mu_nu_find_a_better_name(ao1,ao2,nstate) = 
 END_DOC
 integer :: i,j,ipoint,istate,k,l,coord
 double precision :: grid_sum(3)
 E2_mu_nu_find_a_better_name = 0.d0
 do istate = 1, n_states
  do i= 1, ao_num
   do j = 1, ao_num
    do k = 1, ao_num
     do l = 1, ao_num

      grid_sum = 0d0
      do ipoint = 1, n_points_final_grid 
       do coord = 1,3 
        grid_sum(coord) += vec_prod_grad_ao_at_r(coord,ipoint,j,l) * v_c_ij_grid_transp(ipoint,i,k) * final_weight_at_r_vector(ipoint)
       enddo
      enddo

      do coord = 1,3 
       E2_mu_nu_find_a_better_name(coord,j,i,istate) += rdm_one_e_soc(l,k,istate) * grid_sum(coord) 
      enddo

     enddo
    enddo
   enddo
  enddo
 enddo
 E2_mu_nu_find_a_better_name = -3.d0/2.d0 * E2_mu_nu_find_a_better_name 
END_PROVIDER
