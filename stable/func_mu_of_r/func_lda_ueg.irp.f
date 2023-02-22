
!-----------------------------------------------------------------Integrales------------------------------------------------------------------
 BEGIN_PROVIDER [logical, needs_eff_two_e_ints_lda_ueg]
 implicit none
 BEGIN_DOC
 ! needs_eff_two_e_ints_lda_ueg = False because the self consistent procedure for this functional DOES NOT require to write effective two electron integrals
 END_DOC
 needs_eff_two_e_ints_lda_ueg = .False.

 END_PROVIDER
 
 BEGIN_PROVIDER [double precision, d_dn2_e_cmd_lda_ueg, (n_points_final_grid,N_states)]
 implicit none
 d_dn2_e_cmd_lda_ueg = 0.d0
 END_PROVIDER 


 BEGIN_PROVIDER[double precision, e_c_md_basis_lda_ueg, (N_states) ]
 implicit none
 BEGIN_DOC
 ! Bla bla bla
 ! 
 END_DOC 
 integer :: istate,m,ipoint
 double precision :: weight,mu
 double precision :: ec_srmuLDAn
 double precision :: rho_a,rho_b
 double precision :: decdrho_a, decdrho_b, d2ecdrho_a2,d2ecdrho_b2

 do istate = 1, N_states
  do ipoint = 1, n_points_final_grid

   weight = final_weight_at_r_vector(ipoint)
   rho_a =  one_e_dm_and_grad_alpha_in_r(4,ipoint,istate)
   rho_b =  one_e_dm_and_grad_beta_in_r(4,ipoint,istate)
   
   mu = mu_of_r_prov(ipoint,istate)
!   call ecmdsrPBE(mu,rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b,ec_srmuPBE,decdrho_a,decdrho_b, decdgrad_rho_a_2,decdgrad_rho_b_2,decdgrad_rho_a_b)
    call ecmdsrLDAn(mu,rho_a,rho_b,ec_srmuLDAn,decdrho_a,decdrho_b,d2ecdrho_a2,d2ecdrho_b2)

   e_c_md_basis_lda_ueg(istate) += ec_srmuLDAn * weight
  enddo
 enddo

END_PROVIDER


  BEGIN_PROVIDER [double precision, pot_basis_alpha_mo_basis_lda_ueg,(mo_num,mo_num,N_states)]
 &BEGIN_PROVIDER [double precision, pot_basis_beta_mo_basis_lda_ueg, (mo_num,mo_num,N_states)]
  implicit none
  integer :: istate

  do istate = 1, N_states
     call ao_to_mo(                                                   &
         pot_basis_alpha_ao_lda_ueg(1,1,istate),                                 &
         size(pot_basis_alpha_ao_lda_ueg,1),                                &
         pot_basis_alpha_mo_basis_lda_ueg(1,1,istate),                                 &
         size(pot_basis_alpha_mo_basis_lda_ueg,1)                                 &
         )

     call ao_to_mo(                                                   &
         pot_basis_beta_ao_lda_ueg(1,1,istate),                                  &
         size(pot_basis_beta_ao_lda_ueg,1),                                 &
         pot_basis_beta_mo_basis_lda_ueg(1,1,istate),                                  &
         size(pot_basis_beta_mo_basis_lda_ueg,1)                                  &
         )
  enddo
                                                                                                                                             
 END_PROVIDER


 BEGIN_PROVIDER [double precision, pot_basis_alpha_ao_lda_ueg,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, pot_basis_beta_ao_lda_ueg,(ao_num,ao_num,N_states)]
   implicit none
   integer :: i,j,istate
   do istate = 1, n_states 
    do i = 1, ao_num
     do j = 1, ao_num

      pot_basis_alpha_ao_lda_ueg(j,i,istate) = pot_scal_c_alpha_ao_lda_ueg(j,i,istate) 
      pot_basis_beta_ao_lda_ueg(j,i,istate) = pot_scal_c_beta_ao_lda_ueg(j,i,istate) 
     enddo
    enddo
   enddo

END_PROVIDER 

 BEGIN_PROVIDER[double precision, aos_vc_alpha_basis_lda_ueg_w    , (ao_num,n_points_final_grid,N_states)]
&BEGIN_PROVIDER[double precision, aos_vc_beta_basis_lda_ueg_w    , (ao_num,n_points_final_grid,N_states)]
 implicit none
 BEGIN_DOC
! intermediates to compute the sr_pbe potentials 
! 
! aos_vxc_alpha_pbe_w(j,i) = ao_i(r_j) * (v^x_alpha(r_j) + v^c_alpha(r_j)) * W(r_j)
 END_DOC
 integer :: istate,i,j,m
 double precision :: weight,mu
 double precision :: ec_srmuLDAn
 double precision :: rho_a,rho_b
 double precision :: decdrho_a, decdrho_b

 do istate = 1, N_states
  do i = 1, n_points_final_grid
   weight = final_weight_at_r_vector(i)
   rho_a =  one_e_dm_and_grad_alpha_in_r(4,i,istate)
   rho_b =  one_e_dm_and_grad_beta_in_r(4,i,istate)

   mu = mu_of_r_prov(i,istate)
!   call ecmdsrPBE(mu,rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b,ec_srmuPBE,decdrho_a,decdrho_b, decdgrad_rho_a_2,decdgrad_rho_b_2,decdgrad_rho_a_b)
   call ecmdsrLDAn(mu,rho_a,rho_b,ec_srmuLDAn,decdrho_a,decdrho_b)

   decdrho_a *= weight
   decdrho_b *= weight

   do j = 1, ao_num
    aos_vc_alpha_basis_lda_ueg_w(j,i,istate) = decdrho_a * aos_in_r_array(j,i)
    aos_vc_beta_basis_lda_ueg_w (j,i,istate) = decdrho_b * aos_in_r_array(j,i)
   enddo

  enddo
 enddo

 END_PROVIDER

 BEGIN_PROVIDER [double precision, pot_scal_c_alpha_ao_lda_ueg, (ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, pot_scal_c_beta_ao_lda_ueg, (ao_num,ao_num,N_states)]
 implicit none
! intermediates to compute the sr_pbe potentials 
! 
 integer                        :: istate
   BEGIN_DOC
   ! intermediate quantity for the calculation of the vxc potentials for the GGA functionals  related to the scalar part of the potential 
   END_DOC
   pot_scal_c_alpha_ao_lda_ueg = 0.d0
   pot_scal_c_beta_ao_lda_ueg = 0.d0
   double precision               :: wall_1,wall_2
   call wall_time(wall_1)
   do istate = 1, N_states
     ! correlation alpha
     call dgemm('N','T',ao_num,ao_num,n_points_final_grid,1.d0,                   &
                 aos_vc_alpha_basis_lda_ueg_w(1,1,istate),size(aos_vc_alpha_basis_lda_ueg_w,1), &
                 aos_in_r_array,size(aos_in_r_array,1),1.d0,                      &
                 pot_scal_c_alpha_ao_lda_ueg(1,1,istate),size(pot_scal_c_alpha_ao_lda_ueg,1))
     ! correlation beta
     call dgemm('N','T',ao_num,ao_num,n_points_final_grid,1.d0,                   &
                 aos_vc_beta_basis_lda_ueg_w(1,1,istate),size(aos_vc_beta_basis_lda_ueg_w,1),   &
                 aos_in_r_array,size(aos_in_r_array,1),1.d0,                      &
                 pot_scal_c_beta_ao_lda_ueg(1,1,istate),size(pot_scal_c_beta_ao_lda_ueg,1))
 
   enddo
 call wall_time(wall_2)

END_PROVIDER 

! BEGIN_PROVIDER [double precision, pot_grad_c_alpha_ao_lda_ueg,(ao_num,ao_num,N_states)]
!&BEGIN_PROVIDER [double precision, pot_grad_c_beta_ao_lda_ueg,(ao_num,ao_num,N_states)]
!  implicit none
!  BEGIN_DOC
!  ! intermediate quantity for the calculation of the vxc potentials for the GGA functionals  related to the gradienst of the density and orbitals 
!  END_DOC
!  integer                        :: istate
!  double precision               :: wall_1,wall_2
!  call wall_time(wall_1)
!  pot_grad_c_alpha_ao_lda_ueg = 0.d0
!  pot_grad_c_beta_ao_lda_ueg = 0.d0
!  do istate = 1, N_states
!      ! correlation alpha
!      call dgemm('N','N',ao_num,ao_num,n_points_final_grid,1.d0,                     &
!                 aos_d_vc_alpha_basis_lda_ueg_w(1,1,istate),size(aos_d_vc_alpha_basis_lda_ueg_w,1),  &
!                 aos_in_r_array_transp,size(aos_in_r_array_transp,1),1.d0,           &
!                 pot_grad_c_alpha_ao_lda_ueg(1,1,istate),size(pot_grad_c_alpha_ao_lda_ueg,1))
!      ! correlation beta
!      call dgemm('N','N',ao_num,ao_num,n_points_final_grid,1.d0,                     &
!                 aos_d_vc_beta_basis_lda_ueg_w(1,1,istate),size(aos_d_vc_beta_basis_lda_ueg_w,1),    &
!                 aos_in_r_array_transp,size(aos_in_r_array_transp,1),1.d0,           &
!                 pot_grad_c_beta_ao_lda_ueg(1,1,istate),size(pot_grad_c_beta_ao_lda_ueg,1))
!  enddo
!  
!call wall_time(wall_2)

!END_PROVIDER

