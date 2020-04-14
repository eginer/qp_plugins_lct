!-----------------------------------------------------------------Integrales------------------------------------------------------------------

 BEGIN_PROVIDER [logical, needs_eff_two_e_ints_su_pbe_ot_no_extrap]
 implicit none
 BEGIN_DOC
 ! needs_eff_two_e_ints_su_pbe_ot_no_extrap = True because the self consistent procedure for this functional requires to write effective two electron integrals
 END_DOC
! print*,'needs_eff_two_e_ints_su_pbe_ot_no_extrap = ',needs_eff_two_e_ints_su_pbe_ot_no_extrap
 needs_eff_two_e_ints_su_pbe_ot_no_extrap = .True.

 END_PROVIDER

 BEGIN_PROVIDER[double precision, e_c_md_basis_su_pbe_ot_no_extrap, (N_states) ]
 implicit none
 BEGIN_DOC
 ! Bla bla bla
 ! 
 END_DOC 
 integer :: istate,ipoint,j,m
 double precision :: two_dm_in_r_exact
 double precision :: weight, r(3)
 double precision :: ec_srmuPBE, mu
 double precision :: rho2, rho_a,rho_b,grad_rho_a(3),grad_rho_b(3),grad_rho_a_2,grad_rho_b_2,grad_rho_a_b
 double precision :: decdrho_a, decdrho_b, decdrho, decdrho2
 double precision :: decdgrad_rho_a_2,decdgrad_rho_b_2,decdgrad_rho_a_b, decdgrad_rho_2

 e_c_md_basis_su_pbe_ot_no_extrap = 0.d0
 do istate = 1, N_states
  do ipoint = 1, n_points_final_grid
   r(1) = final_grid_points(1,ipoint)
   r(2) = final_grid_points(2,ipoint)
   r(3) = final_grid_points(3,ipoint)

   weight = final_weight_at_r_vector(ipoint)
   rho_a =  one_e_dm_and_grad_alpha_in_r(4,ipoint,istate)
   rho_b =  one_e_dm_and_grad_beta_in_r(4,ipoint,istate)
   
   call give_on_top_in_r_one_state(r,istate,rho2)
   grad_rho_a(1:3) =  one_e_dm_and_grad_alpha_in_r(1:3,ipoint,istate)
   grad_rho_b(1:3) =  one_e_dm_and_grad_beta_in_r(1:3,ipoint,istate)
   grad_rho_a_2 = 0.d0
   grad_rho_b_2 = 0.d0
   grad_rho_a_b = 0.d0
   do m = 1, 3
    grad_rho_a_2 += grad_rho_a(m) * grad_rho_a(m)
    grad_rho_b_2 += grad_rho_b(m) * grad_rho_b(m)
    grad_rho_a_b += grad_rho_a(m) * grad_rho_b(m)
   enddo

   mu = mu_of_r_prov(ipoint,istate)
   rho2 = rho2*2.d0 ! normalization
   call ecmdsrPBEn2(mu,rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b,rho2,ec_srmuPBE,decdrho_a,decdrho_b, decdrho, decdgrad_rho_a_2,decdgrad_rho_b_2,decdgrad_rho_a_b, decdgrad_rho_2,decdrho2)

   decdrho2 = 2.d0*decdrho2 ! normalization
   e_c_md_basis_su_pbe_ot_no_extrap(istate) += ec_srmuPBE * weight
  enddo
 enddo

END_PROVIDER


 BEGIN_PROVIDER [double precision, pot_basis_alpha_mo_su_pbe_ot_no_extrap,(mo_num,mo_num,N_states)]
&BEGIN_PROVIDER [double precision, pot_basis_beta_mo_su_pbe_ot_no_extrap,(mo_num,mo_num,N_states)]
   implicit none
 call ao_to_mo(pot_basis_alpha_ao_su_pbe_ot_no_extrap,ao_num,pot_basis_alpha_mo_su_pbe_ot_no_extrap,mo_num)
 call ao_to_mo(pot_basis_beta_ao_su_pbe_ot_no_extrap,ao_num,pot_basis_beta_mo_su_pbe_ot_no_extrap,mo_num)
END_PROVIDER 

!1
 BEGIN_PROVIDER [double precision, pot_basis_alpha_ao_su_pbe_ot_no_extrap,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, pot_basis_beta_ao_su_pbe_ot_no_extrap,(ao_num,ao_num,N_states)]
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
      pot_basis_alpha_ao_su_pbe_ot_no_extrap(j,i,istate) = pot_scal_alpha_ao_su_pbe_ot_no_extrap(j,i,istate) + pot_grad_alpha_ao_su_pbe_ot_no_extrap(j,i,istate) + pot_grad_alpha_ao_su_pbe_ot_no_extrap(i,j,istate)
      pot_basis_beta_ao_su_pbe_ot_no_extrap(j,i,istate) = pot_scal_beta_ao_su_pbe_ot_no_extrap(j,i,istate) + pot_grad_beta_ao_su_pbe_ot_no_extrap(j,i,istate) + pot_grad_beta_ao_su_pbe_ot_no_extrap(i,j,istate)
     enddo
    enddo
   enddo

END_PROVIDER 

!3
 BEGIN_PROVIDER[double precision, aos_vc_alpha_pbe_n2_w  , (ao_num,n_points_final_grid,N_states)]
&BEGIN_PROVIDER[double precision, aos_vc_beta_pbe_n2_w   , (ao_num,n_points_final_grid,N_states)]
&BEGIN_PROVIDER[double precision, aos_d_vc_alpha_pbe_n2_w  , (ao_num,n_points_final_grid,N_states)]
&BEGIN_PROVIDER[double precision, aos_d_vc_beta_pbe_n2_w   ,  (ao_num,n_points_final_grid,N_states)]
&BEGIN_PROVIDER[double precision, d_dn2_e_cmd_su_pbe_ot_no_extrap, (n_points_final_grid,N_states) ]
 implicit none
 BEGIN_DOC
! intermediates to compute the pbe potentials 
 END_DOC
 integer :: istate,ipoint,j,m
 double precision :: weight, r(3)
 double precision :: ec_srmuPBE,mu
 double precision :: rho_a,rho_b,grad_rho_a(3),grad_rho_b(3),grad_rho_a_2,grad_rho_b_2,grad_rho_a_b,rho2
 double precision :: contrib_grad_ca(3),contrib_grad_cb(3)
 double precision :: decdrho_a, decdrho_b, decdrho, decdrho2
 double precision :: decdgrad_rho_a_2,decdgrad_rho_b_2,decdgrad_rho_a_b, decdgrad_rho_2

 aos_d_vc_alpha_pbe_n2_w = 0.d0
 aos_d_vc_beta_pbe_n2_w = 0.d0
 do istate = 1, N_states
  do ipoint = 1, n_points_final_grid
   r(1) = final_grid_points(1,ipoint)
   r(2) = final_grid_points(2,ipoint)
   r(3) = final_grid_points(3,ipoint)
   weight = final_weight_at_r_vector(ipoint)
   rho_a =  one_e_dm_and_grad_alpha_in_r(4,ipoint,istate)
   rho_b =  one_e_dm_and_grad_beta_in_r(4,ipoint,istate)
   call give_on_top_in_r_one_state(r,istate,rho2)
   grad_rho_a(1:3) =  one_e_dm_and_grad_alpha_in_r(1:3,ipoint,istate)
   grad_rho_b(1:3) =  one_e_dm_and_grad_beta_in_r(1:3,ipoint,istate)
   grad_rho_a_2 = 0.d0
   grad_rho_b_2 = 0.d0
   grad_rho_a_b = 0.d0
   do m = 1, 3
    grad_rho_a_2 += grad_rho_a(m) * grad_rho_a(m)
    grad_rho_b_2 += grad_rho_b(m) * grad_rho_b(m)
    grad_rho_a_b += grad_rho_a(m) * grad_rho_b(m)
   enddo
   
   rho2 = 2.d0*rho2
   ! mu_erf_dft -> mu_b
   mu = mu_of_r_prov(ipoint,istate)
   call ecmdsrPBEn2(mu,rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b,rho2,ec_srmuPBE,decdrho_a,decdrho_b, decdrho, decdgrad_rho_a_2,decdgrad_rho_b_2,decdgrad_rho_a_b, decdgrad_rho_2,decdrho2)
   d_dn2_e_cmd_su_pbe_ot_no_extrap(ipoint,istate) = 2.d0 * decdrho2
   
   decdrho_a *= weight
   decdrho_b *= weight

   do m= 1,3
    contrib_grad_ca(m) = weight * (2.d0 * decdgrad_rho_a_2 *  grad_rho_a(m) + decdgrad_rho_a_b  * grad_rho_b(m) )
    contrib_grad_cb(m) = weight * (2.d0 * decdgrad_rho_b_2 *  grad_rho_b(m) + decdgrad_rho_a_b  * grad_rho_a(m) )
   enddo

   do j = 1, ao_num
    aos_vc_alpha_pbe_n2_w(j,ipoint,istate) = decdrho_a * aos_in_r_array(j,ipoint)
    aos_vc_beta_pbe_n2_w (j,ipoint,istate) = decdrho_b * aos_in_r_array(j,ipoint)
   enddo
   do j = 1, ao_num
    do m = 1,3
     aos_d_vc_alpha_pbe_n2_w(j,ipoint,istate) += contrib_grad_ca(m) * aos_grad_in_r_array_transp(m,j,ipoint)
     aos_d_vc_beta_pbe_n2_w (j,ipoint,istate) += contrib_grad_cb(m) * aos_grad_in_r_array_transp(m,j,ipoint)
    enddo
   enddo
  enddo
 enddo

 END_PROVIDER

!4
 BEGIN_PROVIDER [double precision, pot_scal_alpha_ao_su_pbe_ot_no_extrap, (ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, pot_scal_beta_ao_su_pbe_ot_no_extrap, (ao_num,ao_num,N_states)]
 implicit none
! intermediates to compute the pbe potentials 
! 
 integer                        :: istate
   BEGIN_DOC
   ! intermediate quantity for the calculation of the vxc potentials for the GGA functionals  related to the scalar part of the potential 
   END_DOC
   pot_scal_alpha_ao_su_pbe_ot_no_extrap = 0.d0
   pot_scal_beta_ao_su_pbe_ot_no_extrap = 0.d0
   double precision               :: wall_1,wall_2
   call wall_time(wall_1)
   do istate = 1, N_states
     ! correlation alpha
     call dgemm('N','T',ao_num,ao_num,n_points_final_grid,1.d0,                   &
                 aos_vc_alpha_pbe_n2_w(1,1,istate),size(aos_vc_alpha_pbe_n2_w,1), &
                 aos_in_r_array,size(aos_in_r_array,1),1.d0,                      &
                 pot_scal_alpha_ao_su_pbe_ot_no_extrap(1,1,istate),size(pot_scal_alpha_ao_su_pbe_ot_no_extrap,1))
     ! correlation beta
     call dgemm('N','T',ao_num,ao_num,n_points_final_grid,1.d0,                   &
                 aos_vc_beta_pbe_n2_w(1,1,istate),size(aos_vc_beta_pbe_n2_w,1),   &
                 aos_in_r_array,size(aos_in_r_array,1),1.d0,                      &
                 pot_scal_beta_ao_su_pbe_ot_no_extrap(1,1,istate),size(pot_scal_beta_ao_su_pbe_ot_no_extrap,1))
 
   enddo
 call wall_time(wall_2)

END_PROVIDER 

!5
 BEGIN_PROVIDER [double precision, pot_grad_alpha_ao_su_pbe_ot_no_extrap,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, pot_grad_beta_ao_su_pbe_ot_no_extrap,(ao_num,ao_num,N_states)]
   implicit none
   BEGIN_DOC
   ! intermediate quantity for the calculation of the vxc potentials for the GGA functionals  related to the gradienst of the density and orbitals 
   END_DOC
   integer                        :: istate
   double precision               :: wall_1,wall_2
   call wall_time(wall_1)
   pot_grad_alpha_ao_su_pbe_ot_no_extrap = 0.d0
   pot_grad_beta_ao_su_pbe_ot_no_extrap = 0.d0
   do istate = 1, N_states
       ! correlation alpha
       call dgemm('N','N',ao_num,ao_num,n_points_final_grid,1.d0,                     &
                  aos_d_vc_alpha_pbe_n2_w(1,1,istate),size(aos_d_vc_alpha_pbe_n2_w,1),  &
                  aos_in_r_array_transp,size(aos_in_r_array_transp,1),1.d0,           &
                  pot_grad_alpha_ao_su_pbe_ot_no_extrap(1,1,istate),size(pot_grad_alpha_ao_su_pbe_ot_no_extrap,1))
       ! correlation beta
       call dgemm('N','N',ao_num,ao_num,n_points_final_grid,1.d0,                     &
                  aos_d_vc_beta_pbe_n2_w(1,1,istate),size(aos_d_vc_beta_pbe_n2_w,1),    &
                  aos_in_r_array_transp,size(aos_in_r_array_transp,1),1.d0,           &
                  pot_grad_beta_ao_su_pbe_ot_no_extrap(1,1,istate),size(pot_grad_beta_ao_su_pbe_ot_no_extrap,1))
   enddo
   
 call wall_time(wall_2)

END_PROVIDER
