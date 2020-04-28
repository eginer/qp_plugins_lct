 BEGIN_PROVIDER [logical, needs_eff_two_e_ints_lda]
 implicit none
 BEGIN_DOC
 ! needs_eff_two_e_ints_lda = False because the self consistent procedure for this functional DOES NOT require to write effective two electron integrals
 END_DOC
 needs_eff_two_e_ints_lda = .False.

 END_PROVIDER

BEGIN_PROVIDER[double precision, e_c_md_basis_lda, (N_states) ]
 implicit none
 BEGIN_DOC
 ! Bla bla bla
 END_DOC 
 e_c_md_basis_lda = ecmd_lda_mu_of_r
END_PROVIDER

 BEGIN_PROVIDER [double precision, pot_basis_alpha_ao_lda,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, pot_basis_beta_ao_lda,(ao_num,ao_num,N_states)]
   implicit none
 BEGIN_DOC
 ! blablabla bis 
 !
 ! 
 END_DOC 
   integer :: i,j,istate
   do istate = 1, n_states 
    do i = 1, ao_num
     do j = 1, ao_num
      pot_basis_alpha_ao_lda(j,i,istate) = pot_basis_alpha_ao_ecmd_delta_lda(j,i,istate) + pot_basis_alpha_ao_ec_lda(j,i,istate) 
      pot_basis_beta_ao_lda(j,i,istate)  = pot_basis_beta_ao_ecmd_delta_lda(j,i,istate)  + pot_basis_beta_ao_ec_lda(j,i,istate) 
     enddo
    enddo
   enddo

END_PROVIDER 

 BEGIN_PROVIDER [double precision, d_dn2_e_cmd_lda, (n_points_final_grid,N_states)]
 implicit none
 d_dn2_e_cmd_lda = 0.d0
 END_PROVIDER 

 BEGIN_PROVIDER [double precision, pot_basis_alpha_mo_lda,(mo_num,mo_num,N_states)]
&BEGIN_PROVIDER [double precision, pot_basis_beta_mo_lda,(mo_num,mo_num,N_states)]
   implicit none
 call ao_to_mo(pot_basis_alpha_ao_lda,ao_num,pot_basis_alpha_mo_lda,mo_num)
 call ao_to_mo(pot_basis_beta_ao_lda,ao_num,pot_basis_beta_mo_lda,mo_num)
END_PROVIDER 



 BEGIN_PROVIDER [double precision, pot_basis_alpha_ao_ecmd_delta_lda,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, pot_basis_beta_ao_ecmd_delta_lda,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, pot_basis_alpha_ao_ec_lda,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, pot_basis_beta_ao_ec_lda,(ao_num,ao_num,N_states)]
 implicit none
 integer :: istate
 double precision :: wall_1,wall_2
 call wall_time(wall_1)
 do istate = 1, N_states 
  call dgemm('N','N',ao_num,ao_num,n_points_final_grid,1.d0,aos_in_r_array,ao_num,aos_deltarho_w_alpha(1,1,istate),n_points_final_grid,0.d0,pot_basis_alpha_ao_ecmd_delta_lda(1,1,istate),ao_num)
  call dgemm('N','N',ao_num,ao_num,n_points_final_grid,1.d0,aos_in_r_array,ao_num,aos_deltarho_w_beta(1,1,istate),n_points_final_grid,0.d0,pot_basis_beta_ao_ecmd_delta_lda(1,1,istate),ao_num)
  call dgemm('N','N',ao_num,ao_num,n_points_final_grid,1.d0,aos_in_r_array,ao_num,aos_e_c_w_alpha(1,1,istate),n_points_final_grid,0.d0,pot_basis_alpha_ao_ec_lda(1,1,istate),ao_num)
  call dgemm('N','N',ao_num,ao_num,n_points_final_grid,1.d0,aos_in_r_array,ao_num,aos_e_c_w_beta(1,1,istate),n_points_final_grid,0.d0,pot_basis_beta_ao_ec_lda(1,1,istate),ao_num)
 enddo
 call wall_time(wall_2)
 print*,'time to provide pot_basis_alpha_ao_ecmd_delta_lda_2 = ',wall_2 - wall_1
 END_PROVIDER




 BEGIN_PROVIDER[double precision, aos_deltarho_w_alpha, (n_points_final_grid,ao_num,N_states)]
&BEGIN_PROVIDER[double precision, aos_deltarho_w_beta, (n_points_final_grid,ao_num,N_states)]
&BEGIN_PROVIDER[double precision, aos_e_c_w_alpha, (n_points_final_grid,ao_num,N_states)]
&BEGIN_PROVIDER[double precision, aos_e_c_w_beta, (n_points_final_grid,ao_num,N_states)]
 implicit none
 integer :: istate,i,j
 double precision :: mu,weight
 double precision :: d_total_deltarho_rhoa,d_total_deltarho_rhob , e_c,vc_a,vc_b
 double precision, allocatable :: aos_array(:), r(:),rhoa(:),rhob(:)
 allocate(aos_array(ao_num),r(3),rhoa(N_states),rhob(N_states))
 double precision :: threshold
 threshold = 1d-15

 do istate = 1, N_states
  do i = 1, n_points_final_grid
   r(1) = final_grid_points(1,i)
   r(2) = final_grid_points(2,i)
   r(3) = final_grid_points(3,i)
   mu =mu_of_r_prov(i,istate)
   weight=final_weight_at_r_vector(i)

   rhoa(istate) = one_e_dm_and_grad_alpha_in_r(4,i,istate)
   rhob(istate) = one_e_dm_and_grad_beta_in_r(4,i,istate)
   call dm_dft_alpha_beta_and_all_aos_at_r(r,rhoa(istate),rhob(istate),aos_array)
   call ec_lda_sr(mu,rhoa(istate),rhob(istate),e_c,vc_a,vc_b)
   if(dabs(rhoa(istate)+rhob(istate)).lt.threshold) cycle 
   do j = 1, ao_num
    aos_deltarho_w_alpha(i,j,istate) = d_total_deltarho_rhoa(rhoa(istate),rhob(istate),mu)*aos_array(j)*weight
    aos_deltarho_w_beta(i,j,istate)  = d_total_deltarho_rhob(rhoa(istate),rhob(istate),mu)*aos_array(j)*weight
    aos_e_c_w_alpha(i,j,istate)      = vc_a*aos_array(j)*weight
    aos_e_c_w_beta(i,j,istate)       = vc_b*aos_array(j)*weight
   enddo
  enddo
 enddo
 END_PROVIDER

