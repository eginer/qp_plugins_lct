!-----------------------------------------------------------------Integrales------------------------------------------------------------------
 BEGIN_PROVIDER[double precision, energy_x_md_sr_pbe, (N_states) ]
&BEGIN_PROVIDER[double precision, energy_c_md_sr_pbe, (N_states) ]
 implicit none
 BEGIN_DOC
! exchange/correlation energy with the short range pbe functional
 END_DOC
 integer :: istate,i,j,m
 double precision :: weight
 double precision :: ex_srmuPBE, ec_srmuPBE
 double precision :: rho_a,rho_b,grad_rho_a(3),grad_rho_b(3),grad_rho_a_2,grad_rho_b_2,grad_rho_a_b
 double precision :: decdrho_a, decdrho_b, dexdrho_a, dexdrho_b
 double precision :: dexdgrad_rho_a_2,dexdgrad_rho_b_2,dexdgrad_rho_a_b,decdgrad_rho_a_2,decdgrad_rho_b_2,decdgrad_rho_a_b



 energy_x_md_sr_pbe = 0.d0
 energy_c_md_sr_pbe = 0.d0
 do istate = 1, N_states
  do i = 1, n_points_final_grid
   weight = final_weight_at_r_vector(i)
   rho_a =  one_e_dm_and_grad_alpha_in_r(4,i,istate)
   rho_b =  one_e_dm_and_grad_beta_in_r(4,i,istate)
   grad_rho_a(1:3) =  one_e_dm_and_grad_alpha_in_r(1:3,i,istate)
   grad_rho_b(1:3) =  one_e_dm_and_grad_beta_in_r(1:3,i,istate)
   grad_rho_a_2 = 0.d0
   grad_rho_b_2 = 0.d0
   grad_rho_a_b = 0.d0
   do m = 1, 3
    grad_rho_a_2 += grad_rho_a(m) * grad_rho_a(m)
    grad_rho_b_2 += grad_rho_b(m) * grad_rho_b(m)
    grad_rho_a_b += grad_rho_a(m) * grad_rho_b(m)
   enddo
   double precision :: mu_local
   mu_local = mu_of_r_dft(i)
   call exc_dexc_md_sr_PBE(mu_local,rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b, &
       ex_srmuPBE,dexdrho_a,dexdrho_b,dexdgrad_rho_a_2,dexdgrad_rho_b_2,dexdgrad_rho_a_b, &
       ec_srmuPBE,decdrho_a,decdrho_b,decdgrad_rho_a_2,decdgrad_rho_b_2,decdgrad_rho_a_b)

   energy_x_md_sr_pbe(istate) += ex_srmuPBE * weight
   energy_c_md_sr_pbe(istate) += ec_srmuPBE * weight
  enddo
 enddo

END_PROVIDER


!1
 BEGIN_PROVIDER [double precision, potential_x_alpha_ao_md_sr_pbe,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, potential_x_beta_ao_md_sr_pbe,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, potential_c_alpha_ao_md_sr_pbe,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, potential_c_beta_ao_md_sr_pbe,(ao_num,ao_num,N_states)]
   implicit none
 BEGIN_DOC
 ! exchange / correlation potential for alpha / beta electrons  with the short-range version Perdew-Burke-Ernzerhof GGA functional 
 !
 ! defined in Chem. Phys.329, 276 (2006)
 END_DOC 
   integer :: i,j,istate
   do istate = 1, n_states 
    do i = 1, ao_num
     do j = 1, ao_num
      potential_x_alpha_ao_md_sr_pbe(j,i,istate) = pot_scal_x_alpha_ao_md_sr_pbe(j,i,istate) + pot_grad_x_alpha_ao_md_sr_pbe(j,i,istate) + pot_grad_x_alpha_ao_md_sr_pbe(i,j,istate)
      potential_x_beta_ao_md_sr_pbe(j,i,istate) = pot_scal_x_beta_ao_md_sr_pbe(j,i,istate) + pot_grad_x_beta_ao_md_sr_pbe(j,i,istate) + pot_grad_x_beta_ao_md_sr_pbe(i,j,istate)

      potential_c_alpha_ao_md_sr_pbe(j,i,istate) = pot_scal_c_alpha_ao_md_sr_pbe(j,i,istate) + pot_grad_c_alpha_ao_md_sr_pbe(j,i,istate) + pot_grad_c_alpha_ao_md_sr_pbe(i,j,istate)
      potential_c_beta_ao_md_sr_pbe(j,i,istate) = pot_scal_c_beta_ao_md_sr_pbe(j,i,istate) + pot_grad_c_beta_ao_md_sr_pbe(j,i,istate) + pot_grad_c_beta_ao_md_sr_pbe(i,j,istate)
     enddo
    enddo
   enddo

END_PROVIDER 

!2
 BEGIN_PROVIDER [double precision, potential_xc_alpha_ao_md_sr_pbe,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, potential_xc_beta_ao_md_sr_pbe,(ao_num,ao_num,N_states)]
   implicit none
 BEGIN_DOC
 ! exchange / correlation potential for alpha / beta electrons  with the Perdew-Burke-Ernzerhof GGA functional 
 END_DOC 
   integer :: i,j,istate
   do istate = 1, n_states 
    do i = 1, ao_num
     do j = 1, ao_num
      potential_xc_alpha_ao_md_sr_pbe(j,i,istate) = pot_scal_xc_alpha_ao_md_sr_pbe(j,i,istate) + pot_grad_xc_alpha_ao_md_sr_pbe(j,i,istate) + pot_grad_xc_alpha_ao_md_sr_pbe(i,j,istate)
      potential_xc_beta_ao_md_sr_pbe(j,i,istate)  = pot_scal_xc_beta_ao_md_sr_pbe(j,i,istate)  + pot_grad_xc_beta_ao_md_sr_pbe(j,i,istate)  + pot_grad_xc_beta_ao_md_sr_pbe(i,j,istate)
     enddo
    enddo
   enddo

END_PROVIDER 


!3
 BEGIN_PROVIDER[double precision, aos_vc_alpha_md_sr_pbe_w  , (ao_num,n_points_final_grid,N_states)]
&BEGIN_PROVIDER[double precision, aos_vc_beta_md_sr_pbe_w   , (ao_num,n_points_final_grid,N_states)]
&BEGIN_PROVIDER[double precision, aos_vx_alpha_md_sr_pbe_w  , (ao_num,n_points_final_grid,N_states)]
&BEGIN_PROVIDER[double precision, aos_vx_beta_md_sr_pbe_w   , (ao_num,n_points_final_grid,N_states)]
&BEGIN_PROVIDER[double precision, aos_d_vc_alpha_md_sr_pbe_w  , (ao_num,n_points_final_grid,N_states)]
&BEGIN_PROVIDER[double precision, aos_d_vc_beta_md_sr_pbe_w   ,  (ao_num,n_points_final_grid,N_states)]
&BEGIN_PROVIDER[double precision, aos_d_vx_alpha_md_sr_pbe_w  , (ao_num,n_points_final_grid,N_states)]
&BEGIN_PROVIDER[double precision, aos_d_vx_beta_md_sr_pbe_w   ,  (ao_num,n_points_final_grid,N_states)]
 implicit none
 BEGIN_DOC
! intermediates to compute the sr_pbe potentials 
! 
! aos_vxc_alpha_pbe_w(j,i) = ao_i(r_j) * (v^x_alpha(r_j) + v^c_alpha(r_j)) * W(r_j)
 END_DOC
 integer :: istate,i,j,m
 double precision :: weight
 double precision :: ex_srmuPBE, ec_srmuPBE
 double precision :: rho_a,rho_b,grad_rho_a(3),grad_rho_b(3),grad_rho_a_2,grad_rho_b_2,grad_rho_a_b
 double precision :: contrib_grad_xa(3),contrib_grad_xb(3),contrib_grad_ca(3),contrib_grad_cb(3)
 double precision :: decdrho_a, decdrho_b, dexdrho_a, dexdrho_b
 double precision :: dexdgrad_rho_a_2,dexdgrad_rho_b_2,dexdgrad_rho_a_b,decdgrad_rho_a_2,decdgrad_rho_b_2,decdgrad_rho_a_b

 aos_d_vc_alpha_md_sr_pbe_w= 0.d0
 aos_d_vc_beta_md_sr_pbe_w = 0.d0
 aos_d_vx_alpha_md_sr_pbe_w= 0.d0
 aos_d_vx_beta_md_sr_pbe_w = 0.d0
 do istate = 1, N_states
  do i = 1, n_points_final_grid
   weight = final_weight_at_r_vector(i)

   rho_a =  one_e_dm_and_grad_alpha_in_r(4,i,istate)
   rho_b =  one_e_dm_and_grad_beta_in_r(4,i,istate)
   grad_rho_a(1:3) =  one_e_dm_and_grad_alpha_in_r(1:3,i,istate)
   grad_rho_b(1:3) =  one_e_dm_and_grad_beta_in_r(1:3,i,istate)
   grad_rho_a_2 = 0.d0
   grad_rho_b_2 = 0.d0
   grad_rho_a_b = 0.d0
   do m = 1, 3
    grad_rho_a_2 += grad_rho_a(m) * grad_rho_a(m)
    grad_rho_b_2 += grad_rho_b(m) * grad_rho_b(m)
    grad_rho_a_b += grad_rho_a(m) * grad_rho_b(m)
   enddo

   double precision :: mu_local
   mu_local = mu_of_r_dft(i)
   call exc_dexc_md_sr_PBE(mu_local,rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b, &
       ex_srmuPBE,dexdrho_a,dexdrho_b,dexdgrad_rho_a_2,dexdgrad_rho_b_2,dexdgrad_rho_a_b, &
       ec_srmuPBE,decdrho_a,decdrho_b,decdgrad_rho_a_2,decdgrad_rho_b_2,decdgrad_rho_a_b)

   dexdrho_a *= weight
   decdrho_a *= weight
   dexdrho_b *= weight
   decdrho_b *= weight

   do m= 1,3
    contrib_grad_ca(m) = weight * (2.d0 * decdgrad_rho_a_2 *  grad_rho_a(m) + decdgrad_rho_a_b  * grad_rho_b(m) )
    contrib_grad_xa(m) = weight * (2.d0 * dexdgrad_rho_a_2 *  grad_rho_a(m) + dexdgrad_rho_a_b  * grad_rho_b(m) )
    contrib_grad_cb(m) = weight * (2.d0 * decdgrad_rho_b_2 *  grad_rho_b(m) + decdgrad_rho_a_b  * grad_rho_a(m) )
    contrib_grad_xb(m) = weight * (2.d0 * dexdgrad_rho_b_2 *  grad_rho_b(m) + dexdgrad_rho_a_b  * grad_rho_a(m) )    
   enddo

   do j = 1, ao_num
    aos_vc_alpha_md_sr_pbe_w(j,i,istate) = decdrho_a * aos_in_r_array(j,i)
    aos_vc_beta_md_sr_pbe_w (j,i,istate) = decdrho_b * aos_in_r_array(j,i)
    aos_vx_alpha_md_sr_pbe_w(j,i,istate) = dexdrho_a * aos_in_r_array(j,i)
    aos_vx_beta_md_sr_pbe_w (j,i,istate) = dexdrho_b * aos_in_r_array(j,i)
   enddo
   do j = 1, ao_num
    do m = 1,3
     aos_d_vc_alpha_md_sr_pbe_w(j,i,istate) += contrib_grad_ca(m) * aos_grad_in_r_array_transp(m,j,i)
     aos_d_vc_beta_md_sr_pbe_w (j,i,istate) += contrib_grad_cb(m) * aos_grad_in_r_array_transp(m,j,i)
     aos_d_vx_alpha_md_sr_pbe_w(j,i,istate) += contrib_grad_xa(m) * aos_grad_in_r_array_transp(m,j,i)
     aos_d_vx_beta_md_sr_pbe_w (j,i,istate) += contrib_grad_xb(m) * aos_grad_in_r_array_transp(m,j,i)
    enddo
   enddo
  enddo
 enddo

 END_PROVIDER

!4
 BEGIN_PROVIDER [double precision, pot_scal_x_alpha_ao_md_sr_pbe, (ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, pot_scal_c_alpha_ao_md_sr_pbe, (ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, pot_scal_x_beta_ao_md_sr_pbe, (ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, pot_scal_c_beta_ao_md_sr_pbe, (ao_num,ao_num,N_states)]
 implicit none
! intermediates to compute the sr_pbe potentials 
! 
 integer                        :: istate
   BEGIN_DOC
   ! intermediate quantity for the calculation of the vxc potentials for the GGA functionals  related to the scalar part of the potential 
   END_DOC
   pot_scal_c_alpha_ao_md_sr_pbe = 0.d0
   pot_scal_x_alpha_ao_md_sr_pbe = 0.d0
   pot_scal_c_beta_ao_md_sr_pbe = 0.d0
   pot_scal_x_beta_ao_md_sr_pbe = 0.d0
   double precision               :: wall_1,wall_2
   call wall_time(wall_1)
   do istate = 1, N_states
     ! correlation alpha
     call dgemm('N','T',ao_num,ao_num,n_points_final_grid,1.d0,                   &
                 aos_vc_alpha_md_sr_pbe_w(1,1,istate),size(aos_vc_alpha_md_sr_pbe_w,1), &
                 aos_in_r_array,size(aos_in_r_array,1),1.d0,                      &
                 pot_scal_c_alpha_ao_md_sr_pbe(1,1,istate),size(pot_scal_c_alpha_ao_md_sr_pbe,1))
     ! correlation beta
     call dgemm('N','T',ao_num,ao_num,n_points_final_grid,1.d0,                   &
                 aos_vc_beta_md_sr_pbe_w(1,1,istate),size(aos_vc_beta_md_sr_pbe_w,1),   &
                 aos_in_r_array,size(aos_in_r_array,1),1.d0,                      &
                 pot_scal_c_beta_ao_md_sr_pbe(1,1,istate),size(pot_scal_c_beta_ao_md_sr_pbe,1))
     ! exchange alpha
     call dgemm('N','T',ao_num,ao_num,n_points_final_grid,1.d0,                   &
                 aos_vx_alpha_md_sr_pbe_w(1,1,istate),size(aos_vx_alpha_md_sr_pbe_w,1), &
                 aos_in_r_array,size(aos_in_r_array,1),1.d0,                      &
                 pot_scal_x_alpha_ao_md_sr_pbe(1,1,istate),size(pot_scal_x_alpha_ao_md_sr_pbe,1))
     ! exchange beta
     call dgemm('N','T',ao_num,ao_num,n_points_final_grid,1.d0,                   &
                 aos_vx_beta_md_sr_pbe_w(1,1,istate),size(aos_vx_beta_md_sr_pbe_w,1),   &
                 aos_in_r_array,size(aos_in_r_array,1),1.d0,                      &
                 pot_scal_x_beta_ao_md_sr_pbe(1,1,istate), size(pot_scal_x_beta_ao_md_sr_pbe,1))
 
   enddo
 call wall_time(wall_2)

END_PROVIDER 

!5
 BEGIN_PROVIDER [double precision, pot_grad_x_alpha_ao_md_sr_pbe,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, pot_grad_x_beta_ao_md_sr_pbe,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, pot_grad_c_alpha_ao_md_sr_pbe,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, pot_grad_c_beta_ao_md_sr_pbe,(ao_num,ao_num,N_states)]
   implicit none
   BEGIN_DOC
   ! intermediate quantity for the calculation of the vxc potentials for the GGA functionals  related to the gradienst of the density and orbitals 
   END_DOC
   integer                        :: istate
   double precision               :: wall_1,wall_2
   call wall_time(wall_1)
   pot_grad_c_alpha_ao_md_sr_pbe = 0.d0
   pot_grad_x_alpha_ao_md_sr_pbe = 0.d0
   pot_grad_c_beta_ao_md_sr_pbe = 0.d0
   pot_grad_x_beta_ao_md_sr_pbe = 0.d0
   do istate = 1, N_states
       ! correlation alpha
       call dgemm('N','N',ao_num,ao_num,n_points_final_grid,1.d0,                     &
                  aos_d_vc_alpha_md_sr_pbe_w(1,1,istate),size(aos_d_vc_alpha_md_sr_pbe_w,1),  &
                  aos_in_r_array_transp,size(aos_in_r_array_transp,1),1.d0,           &
                  pot_grad_c_alpha_ao_md_sr_pbe(1,1,istate),size(pot_grad_c_alpha_ao_md_sr_pbe,1))
       ! correlation beta
       call dgemm('N','N',ao_num,ao_num,n_points_final_grid,1.d0,                     &
                  aos_d_vc_beta_md_sr_pbe_w(1,1,istate),size(aos_d_vc_beta_md_sr_pbe_w,1),    &
                  aos_in_r_array_transp,size(aos_in_r_array_transp,1),1.d0,           &
                  pot_grad_c_beta_ao_md_sr_pbe(1,1,istate),size(pot_grad_c_beta_ao_md_sr_pbe,1))
       ! exchange alpha
       call dgemm('N','N',ao_num,ao_num,n_points_final_grid,1.d0,                     &
                  aos_d_vx_alpha_md_sr_pbe_w(1,1,istate),size(aos_d_vx_alpha_md_sr_pbe_w,1),  &
                  aos_in_r_array_transp,size(aos_in_r_array_transp,1),1.d0,           &
                  pot_grad_x_alpha_ao_md_sr_pbe(1,1,istate),size(pot_grad_x_alpha_ao_md_sr_pbe,1))
       ! exchange beta
       call dgemm('N','N',ao_num,ao_num,n_points_final_grid,1.d0,                     &
                  aos_d_vx_beta_md_sr_pbe_w(1,1,istate),size(aos_d_vx_beta_md_sr_pbe_w,1),    &
                  aos_in_r_array_transp,size(aos_in_r_array_transp,1),1.d0,           &
                  pot_grad_x_beta_ao_md_sr_pbe(1,1,istate),size(pot_grad_x_beta_ao_md_sr_pbe,1))
   enddo
   
 call wall_time(wall_2)

END_PROVIDER

!6
 BEGIN_PROVIDER[double precision, aos_vxc_alpha_md_sr_pbe_w  , (ao_num,n_points_final_grid,N_states)]  ! sr_pbe
&BEGIN_PROVIDER[double precision, aos_vxc_beta_md_sr_pbe_w   , (ao_num,n_points_final_grid,N_states)]
&BEGIN_PROVIDER[double precision, aos_d_vxc_alpha_md_sr_pbe_w  , (ao_num,n_points_final_grid,N_states)]
&BEGIN_PROVIDER[double precision, aos_d_vxc_beta_md_sr_pbe_w   ,  (ao_num,n_points_final_grid,N_states)]
 implicit none
 BEGIN_DOC
! aos_vxc_alpha_pbe_w(j,i) = ao_i(r_j) * (v^x_alpha(r_j) + v^c_alpha(r_j)) * W(r_j)
 END_DOC
 integer :: istate,i,j,m
 double precision :: weight
 double precision :: ex_srmuPBE, ec_srmuPBE
 double precision :: rho_a,rho_b,grad_rho_a(3),grad_rho_b(3),grad_rho_a_2,grad_rho_b_2,grad_rho_a_b
 double precision :: contrib_grad_xa(3),contrib_grad_xb(3),contrib_grad_ca(3),contrib_grad_cb(3)
 double precision :: decdrho_a, decdrho_b, dexdrho_a, dexdrho_b
 double precision :: dexdgrad_rho_a_2,dexdgrad_rho_b_2,dexdgrad_rho_a_b,decdgrad_rho_a_2,decdgrad_rho_b_2,decdgrad_rho_a_b
 
 aos_d_vxc_alpha_md_sr_pbe_w = 0.d0
 aos_d_vxc_beta_md_sr_pbe_w = 0.d0

 do istate = 1, N_states
  do i = 1, n_points_final_grid
   weight = final_weight_at_r_vector(i)
   rho_a =  one_e_dm_and_grad_alpha_in_r(4,i,istate)
   rho_b =  one_e_dm_and_grad_beta_in_r(4,i,istate)
   grad_rho_a(1:3) =  one_e_dm_and_grad_alpha_in_r(1:3,i,istate)
   grad_rho_b(1:3) =  one_e_dm_and_grad_beta_in_r(1:3,i,istate)
   grad_rho_a_2 = 0.d0
   grad_rho_b_2 = 0.d0
   grad_rho_a_b = 0.d0
   do m = 1, 3
    grad_rho_a_2 += grad_rho_a(m) * grad_rho_a(m)
    grad_rho_b_2 += grad_rho_b(m) * grad_rho_b(m)
    grad_rho_a_b += grad_rho_a(m) * grad_rho_b(m)
   enddo

   double precision :: mu_local
   mu_local = mu_of_r_dft(i)
   call exc_dexc_md_sr_PBE(mu_local,rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b, &
       ex_srmuPBE,dexdrho_a,dexdrho_b,dexdgrad_rho_a_2,dexdgrad_rho_b_2,dexdgrad_rho_a_b, &
       ec_srmuPBE,decdrho_a,decdrho_b,decdgrad_rho_a_2,decdgrad_rho_b_2,decdgrad_rho_a_b)
 
   dexdrho_a *= weight
   decdrho_a *= weight
   dexdrho_b *= weight
   decdrho_b *= weight
   do m= 1,3
    contrib_grad_ca(m) = weight * (2.d0 * decdgrad_rho_a_2 *  grad_rho_a(m) + decdgrad_rho_a_b  * grad_rho_b(m) )
    contrib_grad_xa(m) = weight * (2.d0 * dexdgrad_rho_a_2 *  grad_rho_a(m) + dexdgrad_rho_a_b  * grad_rho_b(m) )
    contrib_grad_cb(m) = weight * (2.d0 * decdgrad_rho_b_2 *  grad_rho_b(m) + decdgrad_rho_a_b  * grad_rho_a(m) )
    contrib_grad_xb(m) = weight * (2.d0 * dexdgrad_rho_b_2 *  grad_rho_b(m) + dexdgrad_rho_a_b  * grad_rho_a(m) )
   enddo
   do j = 1, ao_num
    aos_vxc_alpha_md_sr_pbe_w(j,i,istate) = ( decdrho_a + dexdrho_a ) * aos_in_r_array(j,i)
    aos_vxc_beta_md_sr_pbe_w (j,i,istate) = ( decdrho_b + dexdrho_b ) * aos_in_r_array(j,i)
   enddo
   do j = 1, ao_num
    do m = 1,3
     aos_d_vxc_alpha_md_sr_pbe_w(j,i,istate) += ( contrib_grad_ca(m) + contrib_grad_xa(m) ) * aos_grad_in_r_array_transp(m,j,i)
     aos_d_vxc_beta_md_sr_pbe_w (j,i,istate) += ( contrib_grad_cb(m) + contrib_grad_xb(m) ) * aos_grad_in_r_array_transp(m,j,i)
    enddo
   enddo
  enddo
 enddo

 END_PROVIDER

!7
 BEGIN_PROVIDER [double precision, pot_scal_xc_alpha_ao_md_sr_pbe, (ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, pot_scal_xc_beta_ao_md_sr_pbe, (ao_num,ao_num,N_states)]
 implicit none
 integer                        :: istate
   BEGIN_DOC
   ! intermediate quantity for the calculation of the vxc potentials for the GGA functionals  related to the scalar part of the potential 
   END_DOC
   pot_scal_xc_alpha_ao_md_sr_pbe = 0.d0
   pot_scal_xc_beta_ao_md_sr_pbe = 0.d0
   double precision               :: wall_1,wall_2
   call wall_time(wall_1)
   do istate = 1, N_states
     ! exchange - correlation alpha
     call dgemm('N','T',ao_num,ao_num,n_points_final_grid,1.d0,                     &
                 aos_vxc_alpha_md_sr_pbe_w(1,1,istate),size(aos_vxc_alpha_md_sr_pbe_w,1), &
                 aos_in_r_array,size(aos_in_r_array,1),1.d0,                        &
                 pot_scal_xc_alpha_ao_md_sr_pbe(1,1,istate),size(pot_scal_xc_alpha_ao_md_sr_pbe,1))
     ! exchange - correlation beta
     call dgemm('N','T',ao_num,ao_num,n_points_final_grid,1.d0,                     &
                 aos_vxc_beta_md_sr_pbe_w(1,1,istate),size(aos_vxc_beta_md_sr_pbe_w,1),   &
                 aos_in_r_array,size(aos_in_r_array,1),1.d0,                        &
                 pot_scal_xc_beta_ao_md_sr_pbe(1,1,istate),size(pot_scal_xc_beta_ao_md_sr_pbe,1))
   enddo
 call wall_time(wall_2)

END_PROVIDER 

!8
 BEGIN_PROVIDER [double precision, pot_grad_xc_alpha_ao_md_sr_pbe,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, pot_grad_xc_beta_ao_md_sr_pbe,(ao_num,ao_num,N_states)]
   implicit none
   BEGIN_DOC
   ! intermediate quantity for the calculation of the vxc potentials for the GGA functionals  related to the gradienst of the density and orbitals 
   END_DOC
   integer                        :: istate
   double precision               :: wall_1,wall_2
   call wall_time(wall_1)
   pot_grad_xc_alpha_ao_md_sr_pbe = 0.d0
   pot_grad_xc_beta_ao_md_sr_pbe = 0.d0
   do istate = 1, N_states
       ! exchange - correlation alpha
       call dgemm('N','N',ao_num,ao_num,n_points_final_grid,1.d0,                      &
                  aos_d_vxc_alpha_md_sr_pbe_w(1,1,istate),size(aos_d_vxc_alpha_md_sr_pbe_w,1), &
                  aos_in_r_array_transp,size(aos_in_r_array_transp,1),1.d0,            &
                  pot_grad_xc_alpha_ao_md_sr_pbe(1,1,istate),size(pot_grad_xc_alpha_ao_md_sr_pbe,1))
       ! exchange - correlation beta
       call dgemm('N','N',ao_num,ao_num,n_points_final_grid,1.d0,                      &
                  aos_d_vxc_beta_md_sr_pbe_w(1,1,istate),size(aos_d_vxc_beta_md_sr_pbe_w,1),   &
                  aos_in_r_array_transp,size(aos_in_r_array_transp,1),1.d0,            &
                  pot_grad_xc_beta_ao_md_sr_pbe(1,1,istate),size(pot_grad_xc_beta_ao_md_sr_pbe,1))
   enddo
   
 call wall_time(wall_2)

END_PROVIDER

