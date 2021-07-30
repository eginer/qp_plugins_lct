 BEGIN_PROVIDER [logical, needs_eff_two_e_ints_pbe_ueg]
 implicit none
 BEGIN_DOC
 ! needs_eff_two_e_ints_pbe_ueg = False because the self consistent procedure for this functional DOES NOT require to write effective two electron integrals
 END_DOC
 needs_eff_two_e_ints_pbe_ueg = .False.

 END_PROVIDER


BEGIN_PROVIDER[double precision, e_c_md_basis_pbe_ueg, (N_states) ]
 implicit none
 BEGIN_DOC
 ! Bla bla bla
 END_DOC 
 e_c_md_basis_pbe_ueg = ecmd_pbe_ueg_mu_of_r
END_PROVIDER


 BEGIN_PROVIDER [double precision, pot_basis_alpha_ao_sp_pbe_ueg,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, pot_basis_beta_ao_sp_pbe_ueg,(ao_num,ao_num,N_states)]
   implicit none
 BEGIN_DOC
 ! Correlation potential for alpha / beta electrons  with the PBE UEG functional 
 END_DOC 
   integer :: i,j,istate
   do istate = 1, n_states 
    do i = 1, ao_num
     do j = 1, ao_num
      pot_basis_alpha_ao_sp_pbe_ueg(j,i,istate) = pot_scal_alpha_ao_sp_pbe_ueg(j,i,istate) + pot_grad_alpha_ao_sp_pbe_ueg(j,i,istate) 
      pot_basis_beta_ao_sp_pbe_ueg(j,i,istate)  = pot_scal_beta_ao_sp_pbe_ueg(j,i,istate)  + pot_grad_beta_ao_sp_pbe_ueg(j,i,istate) 
     enddo
    enddo
   enddo

END_PROVIDER



 BEGIN_PROVIDER [double precision, pot_basis_alpha_mo_sp_pbe_ueg,(mo_num,mo_num,N_states)]
&BEGIN_PROVIDER [double precision, pot_basis_beta_mo_sp_pbe_ueg, (mo_num,mo_num,N_states)]
 implicit none
 BEGIN_DOC
! Providers for the alpha/beta correlation potentials of the PBE UEG on the MO basis
 END_DOC
 integer :: istate
 do istate = 1, N_states
    call ao_to_mo(                                                   &
        pot_basis_alpha_ao_sp_pbe_ueg(1,1,istate),                                 &
        size(pot_basis_alpha_ao_sp_pbe_ueg,1),                                &
        pot_basis_alpha_mo_sp_pbe_ueg(1,1,istate),                                 &
        size(pot_basis_alpha_mo_sp_pbe_ueg,1)                                 &
        )

    call ao_to_mo(                                                   &
        pot_basis_beta_ao_sp_pbe_ueg(1,1,istate),                                  &
        size(pot_basis_beta_ao_sp_pbe_ueg,1),                                 &
        pot_basis_beta_mo_sp_pbe_ueg(1,1,istate),                                  &
        size(pot_basis_beta_mo_sp_pbe_ueg,1)                                  &
        )
 enddo

END_PROVIDER


 BEGIN_PROVIDER[double precision, aos_vc_alpha_pbe_ueg_w  , (ao_num,n_points_final_grid,N_states)]
&BEGIN_PROVIDER[double precision, aos_vc_beta_pbe_ueg_w   , (ao_num,n_points_final_grid,N_states)]
&BEGIN_PROVIDER[double precision, aos_dvc_alpha_pbe_ueg_w  , (ao_num,n_points_final_grid,N_states)]
&BEGIN_PROVIDER[double precision, aos_dvc_beta_pbe_ueg_w   ,  (ao_num,n_points_final_grid,N_states)]
 implicit none
 BEGIN_DOC
  ! Intermediate quantities for the calculation of the vc potentials for the PBE UEG functional
 END_DOC
 integer :: istate,i,j,m
 double precision :: r(3)
 double precision :: mu,weight
 double precision :: d_ec_pbeueg_d_ec_pbe(N_states),d_ec_pbeueg_d_grad_n_alpha(3,N_states),d_ec_pbeueg_d_grad_n_beta(3,N_states)
 double precision :: d_ec_pbeueg_rhoa(N_states),d_ec_pbeueg_rhob(N_states)

 aos_dvc_alpha_pbe_ueg_w= 0.d0
 aos_dvc_beta_pbe_ueg_w = 0.d0
 do istate = 1, N_states
  do i = 1, n_points_final_grid
   r(1) = final_grid_points(1,i)
   r(2) = final_grid_points(2,i)
   r(3) = final_grid_points(3,i)
   mu = mu_of_r_prov(i,istate)
   weight = final_weight_at_r_vector(i)

   call give_d_Ec_pbeueg_rho(r,mu,d_ec_pbeueg_rhoa,d_ec_pbeueg_rhob)
   call give_d_Ec_pbeueg_d_grad_n(r,mu,d_ec_pbeueg_d_ec_pbe,d_ec_pbeueg_d_grad_n_alpha,d_ec_pbeueg_d_grad_n_beta)

   do j = 1, ao_num
    aos_vc_alpha_pbe_ueg_w(j,i,istate) = d_ec_pbeueg_rhoa(istate) * weight * aos_in_r_array(j,i)
    aos_vc_beta_pbe_ueg_w (j,i,istate) = d_ec_pbeueg_rhob(istate) * weight * aos_in_r_array(j,i)
   enddo
   do j = 1, ao_num
    do m = 1,3
     aos_dvc_alpha_pbe_ueg_w(j,i,istate) += weight * d_ec_pbeueg_d_grad_n_alpha(m,istate) * aos_grad_in_r_array_transp(m,j,i)
     aos_dvc_beta_pbe_ueg_w (j,i,istate) += weight * d_ec_pbeueg_d_grad_n_beta(m,istate)  * aos_grad_in_r_array_transp(m,j,i)
    enddo

   enddo
  enddo
 enddo

 END_PROVIDER


 BEGIN_PROVIDER [double precision, pot_scal_alpha_ao_sp_pbe_ueg, (ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, pot_scal_beta_ao_sp_pbe_ueg, (ao_num,ao_num,N_states)]
 implicit none
 integer                        :: istate
   BEGIN_DOC
   ! Intermediate quantities for the calculation of the vc potentials related to the scalar part for the PBE UEG functional
   END_DOC
   pot_scal_alpha_ao_sp_pbe_ueg = 0.d0
   pot_scal_beta_ao_sp_pbe_ueg = 0.d0
   double precision               :: wall_1,wall_2
   call wall_time(wall_1)
   do istate = 1, N_states
     ! correlation alpha
     call dgemm('N','T',ao_num,ao_num,n_points_final_grid,1.d0,                                                                                       &
                 aos_vc_alpha_pbe_ueg_w(1,1,istate),size(aos_vc_alpha_pbe_ueg_w,1),                                                                   &
                 aos_in_r_array,size(aos_in_r_array,1),1.d0,                                                                                          &
                 pot_scal_alpha_ao_sp_pbe_ueg(1,1,istate),size(pot_scal_alpha_ao_sp_pbe_ueg,1))
     ! correlation beta
     call dgemm('N','T',ao_num,ao_num,n_points_final_grid,1.d0,                                                                                         &
                 aos_vc_beta_pbe_ueg_w(1,1,istate),size(aos_vc_beta_pbe_ueg_w,1),                                                                       &
                 aos_in_r_array,size(aos_in_r_array,1),1.d0,                                                                                            &
                 pot_scal_beta_ao_sp_pbe_ueg(1,1,istate),size(pot_scal_beta_ao_sp_pbe_ueg,1))
 
   enddo
 call wall_time(wall_2)

END_PROVIDER



 
 BEGIN_PROVIDER [double precision, pot_grad_alpha_ao_sp_pbe_ueg,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, pot_grad_beta_ao_sp_pbe_ueg,(ao_num,ao_num,N_states)]
   implicit none
   BEGIN_DOC
   ! Intermediate quantity for the calculation of the vc potentials for the PBA UEG functional related to the gradienst of the density and orbitals 
   END_DOC
   integer                        :: istate
   double precision               :: wall_1,wall_2
   call wall_time(wall_1)
   pot_grad_alpha_ao_sp_pbe_ueg = 0.d0
   pot_grad_beta_ao_sp_pbe_ueg = 0.d0
   do istate = 1, N_states
       ! correlation alpha
       call dgemm('N','N',ao_num,ao_num,n_points_final_grid,1.d0,                                                                                             &
                  aos_dvc_alpha_pbe_ueg_w(1,1,istate),size(aos_dvc_alpha_pbe_ueg_w,1),                                                                      &
                  aos_in_r_array_transp,size(aos_in_r_array_transp,1),1.d0,                                                                                &
                  pot_grad_alpha_ao_sp_pbe_ueg(1,1,istate),size(pot_grad_alpha_ao_sp_pbe_ueg,1))
       ! correlation beta
       call dgemm('N','N',ao_num,ao_num,n_points_final_grid,1.d0,                                                                                             &
                  aos_dvc_beta_pbe_ueg_w(1,1,istate),size(aos_dvc_beta_pbe_ueg_w,1),                                                                      &
                  aos_in_r_array_transp,size(aos_in_r_array_transp,1),1.d0,                                                                                &
                  pot_grad_beta_ao_sp_pbe_ueg(1,1,istate),size(pot_grad_beta_ao_sp_pbe_ueg,1))
   enddo
   
 call wall_time(wall_2)

END_PROVIDER
