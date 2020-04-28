 BEGIN_PROVIDER [logical, needs_eff_two_e_ints_su_pbe_ot]
 implicit none
 BEGIN_DOC
 ! needs_eff_two_e_ints_su_pbe_ot = True because the self consistent procedure for this functional requires to write effective two electron integrals
 END_DOC
 needs_eff_two_e_ints_su_pbe_ot = .True.

 END_PROVIDER

 BEGIN_PROVIDER[double precision, e_c_md_basis_su_pbe_ot, (N_states) ]
 implicit none
 BEGIN_DOC
!
! Ecmd functional evaluated with mu(r) and depending on 
!    +) the EXTRAPOLATED on-top pair density 
! 
!    +) the total density, density gradients 
!
!    +) !!!!! NO SPIN POLAIRIZATION !!!!! 
!
!   See the provider of ecmd_pbe_on_top_su_mu_of_r for more details
 END_DOC 

 e_c_md_basis_su_pbe_ot = ecmd_pbe_on_top_su_mu_of_r
END_PROVIDER


 BEGIN_PROVIDER [double precision, pot_basis_alpha_mo_su_pbe_ot,(mo_num,mo_num,N_states)]
&BEGIN_PROVIDER [double precision, pot_basis_beta_mo_su_pbe_ot,(mo_num,mo_num,N_states)]
   implicit none
 BEGIN_DOC
 ! Effective one electron potential (on the MO basis) coming from the functional derivative 
 !
 ! with respect to the total density of the SU-PBE-OT functional evaluated with mu(r). 
 !
 ! As the SU-PBE-OT is S_z independent, the alpha and beta potentials are equal
 ! 
 END_DOC 
 call ao_to_mo(pot_basis_alpha_ao_su_pbe_ot,ao_num,pot_basis_alpha_mo_su_pbe_ot,mo_num)
 call ao_to_mo(pot_basis_beta_ao_su_pbe_ot,ao_num,pot_basis_beta_mo_su_pbe_ot,mo_num)
END_PROVIDER 

!1
 BEGIN_PROVIDER [double precision, pot_basis_alpha_ao_su_pbe_ot,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, pot_basis_beta_ao_su_pbe_ot,(ao_num,ao_num,N_states)]
   implicit none
 BEGIN_DOC
 ! Effective one electron potential (on the AO basis) coming from the functional derivative 
 !
 ! with respect to the total density of the SU-PBE-OT functional evaluated with mu(r). 
 !
 ! As the SU-PBE-FULL-OT is S_z independent, the alpha and beta potentials are equal
 ! 
 END_DOC 
   integer :: i,j,istate
   do istate = 1, n_states 
    do i = 1, ao_num
     do j = 1, ao_num
      pot_basis_alpha_ao_su_pbe_ot(j,i,istate) = pot_scal_alpha_ao_su_pbe_ot(j,i,istate) + pot_grad_alpha_ao_su_pbe_ot(j,i,istate) + pot_grad_alpha_ao_su_pbe_ot(i,j,istate)
      pot_basis_beta_ao_su_pbe_ot(j,i,istate) = pot_scal_beta_ao_su_pbe_ot(j,i,istate) + pot_grad_beta_ao_su_pbe_ot(j,i,istate) + pot_grad_beta_ao_su_pbe_ot(i,j,istate)
     enddo
    enddo
   enddo

END_PROVIDER 

!3
 BEGIN_PROVIDER[double precision, aos_vc_alpha_su_pbe_ot_w  , (ao_num,n_points_final_grid,N_states)]
&BEGIN_PROVIDER[double precision, aos_vc_beta_su_pbe_ot_w   , (ao_num,n_points_final_grid,N_states)]
&BEGIN_PROVIDER[double precision, aos_d_vc_alpha_su_pbe_ot_w  , (ao_num,n_points_final_grid,N_states)]
&BEGIN_PROVIDER[double precision, aos_d_vc_beta_su_pbe_ot_w   ,  (ao_num,n_points_final_grid,N_states)]
&BEGIN_PROVIDER[double precision, d_dn2_e_cmd_su_pbe_ot, (n_points_final_grid,N_states) ]
 implicit none
 BEGIN_DOC
! Intermediates to compute the su-pbe-ot potentials 
!
! An important quantity is the "d_dn2_e_cmd_su_pbe_ot" which is the functional derivative with respect to the on-top 
!
! evaluated on a given grid point. 
!
! NOTICE THAT SUCH QUANTITY IS RETURNED WITH NORMALIZATION N(N-1) WHICH IS SUITED FOR DIRECT MODIFICATION OF V_ijkl
 END_DOC
 integer :: istate,ipoint,j,m
 double precision :: weight, r(3)
 double precision :: ec_srmuPBE,mu
 double precision :: rho_a,rho_b,grad_rho_a(3),grad_rho_b(3),grad_rho_a_2,grad_rho_b_2,grad_rho_a_b,rho2
 double precision :: contrib_grad_ca(3),contrib_grad_cb(3)
 double precision :: decdrho_a, decdrho_b, decdrho, decdrho2
 double precision :: decdgrad_rho_a_2,decdgrad_rho_b_2,decdgrad_rho_a_b, decdgrad_rho_2
 double precision :: mu_correction_of_on_top,rho_tot

 aos_d_vc_alpha_su_pbe_ot_w = 0.d0
 aos_d_vc_beta_su_pbe_ot_w = 0.d0
 do istate = 1, N_states
  do ipoint = 1, n_points_final_grid
   r(1) = final_grid_points(1,ipoint)
   r(2) = final_grid_points(2,ipoint)
   r(3) = final_grid_points(3,ipoint)
   weight = final_weight_at_r_vector(ipoint)

   rho_tot = one_e_dm_and_grad_alpha_in_r(4,ipoint,istate) + one_e_dm_and_grad_beta_in_r(4,ipoint,istate)
   ! This ensures that the spin polarization is zero
   rho_a =  rho_tot * 0.5d0 
   rho_b =  rho_tot * 0.5d0 

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
   
   if(mu_of_r_potential == "cas_ful")then
    ! You take the on-top of the CAS wave function which is computed with mu(r) 
    rho2 = on_top_cas_mu_r(ipoint,istate)
   else
    ! You take the on-top of the CAS wave function computed separately
    rho2 = total_cas_on_top_density(ipoint,istate)
   endif

   ! The factor 2 is because the on-top is normalized to N(N-1)/2 
   ! whereas the routine ecmdsrPBEn2 assumes a on-top normalized to N(N-1)
   rho2 = 2.d0*rho2
   ! mu = mu(r)
   mu = mu_of_r_prov(ipoint,istate)
   ! extrapolation toward the exact on-top based on mu(r)
   rho2 = mu_correction_of_on_top(mu,rho2) 

   call ecmdsrPBEn2(mu,rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b,rho2,ec_srmuPBE,decdrho_a,decdrho_b, decdrho, decdgrad_rho_a_2,decdgrad_rho_b_2,decdgrad_rho_a_b, decdgrad_rho_2,decdrho2)

   ! The factor 1 (and not 2) is because the output is with the correct normalization factor (i.e. N(N-1)) 
   ! for the computation of the effective operator of type 1/2  \sum_{ijkl} <ij|kl> a^k a^l a_j a_i
   d_dn2_e_cmd_su_pbe_ot(ipoint,istate) = 1.d0 * decdrho2
   
   decdrho_a *= weight
   decdrho_b *= weight

   do m= 1,3
    contrib_grad_ca(m) = weight * (2.d0 * decdgrad_rho_a_2 *  grad_rho_a(m) + decdgrad_rho_a_b  * grad_rho_b(m) )
    contrib_grad_cb(m) = weight * (2.d0 * decdgrad_rho_b_2 *  grad_rho_b(m) + decdgrad_rho_a_b  * grad_rho_a(m) )
   enddo

   do j = 1, ao_num
    aos_vc_alpha_su_pbe_ot_w(j,ipoint,istate) = decdrho_a * aos_in_r_array(j,ipoint)
    aos_vc_beta_su_pbe_ot_w (j,ipoint,istate) = decdrho_b * aos_in_r_array(j,ipoint)
   enddo
   do j = 1, ao_num
    do m = 1,3
     aos_d_vc_alpha_su_pbe_ot_w(j,ipoint,istate) += contrib_grad_ca(m) * aos_grad_in_r_array_transp(m,j,ipoint)
     aos_d_vc_beta_su_pbe_ot_w (j,ipoint,istate) += contrib_grad_cb(m) * aos_grad_in_r_array_transp(m,j,ipoint)
    enddo
   enddo
  enddo
 enddo

 END_PROVIDER

!4
 BEGIN_PROVIDER [double precision, pot_scal_alpha_ao_su_pbe_ot, (ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, pot_scal_beta_ao_su_pbe_ot, (ao_num,ao_num,N_states)]
 implicit none
 integer                        :: istate
   BEGIN_DOC
   ! Scalar part of the functional derivatives with respect to the density of the su-pbe-ot
   !
   END_DOC
   pot_scal_alpha_ao_su_pbe_ot = 0.d0
   pot_scal_beta_ao_su_pbe_ot = 0.d0
   double precision               :: wall_1,wall_2
   call wall_time(wall_1)
   do istate = 1, N_states
     ! correlation alpha
     call dgemm('N','T',ao_num,ao_num,n_points_final_grid,1.d0,                   &
                 aos_vc_alpha_su_pbe_ot_w(1,1,istate),size(aos_vc_alpha_su_pbe_ot_w,1), &
                 aos_in_r_array,size(aos_in_r_array,1),1.d0,                      &
                 pot_scal_alpha_ao_su_pbe_ot(1,1,istate),size(pot_scal_alpha_ao_su_pbe_ot,1))
     ! correlation beta
     call dgemm('N','T',ao_num,ao_num,n_points_final_grid,1.d0,                   &
                 aos_vc_beta_su_pbe_ot_w(1,1,istate),size(aos_vc_beta_su_pbe_ot_w,1),   &
                 aos_in_r_array,size(aos_in_r_array,1),1.d0,                      &
                 pot_scal_beta_ao_su_pbe_ot(1,1,istate),size(pot_scal_beta_ao_su_pbe_ot,1))
 
   enddo
 call wall_time(wall_2)

END_PROVIDER 

!5
 BEGIN_PROVIDER [double precision, pot_grad_alpha_ao_su_pbe_ot,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, pot_grad_beta_ao_su_pbe_ot,(ao_num,ao_num,N_states)]
   implicit none
   BEGIN_DOC
   ! Gradient part of the functional derivatives with respect to the density of the su-pbe-ot
   !
   END_DOC
   END_DOC
   integer                        :: istate
   double precision               :: wall_1,wall_2
   call wall_time(wall_1)
   pot_grad_alpha_ao_su_pbe_ot = 0.d0
   pot_grad_beta_ao_su_pbe_ot = 0.d0
   do istate = 1, N_states
       ! correlation alpha
       call dgemm('N','N',ao_num,ao_num,n_points_final_grid,1.d0,                     &
                  aos_d_vc_alpha_su_pbe_ot_w(1,1,istate),size(aos_d_vc_alpha_su_pbe_ot_w,1),  &
                  aos_in_r_array_transp,size(aos_in_r_array_transp,1),1.d0,           &
                  pot_grad_alpha_ao_su_pbe_ot(1,1,istate),size(pot_grad_alpha_ao_su_pbe_ot,1))
       ! correlation beta
       call dgemm('N','N',ao_num,ao_num,n_points_final_grid,1.d0,                     &
                  aos_d_vc_beta_su_pbe_ot_w(1,1,istate),size(aos_d_vc_beta_su_pbe_ot_w,1),    &
                  aos_in_r_array_transp,size(aos_in_r_array_transp,1),1.d0,           &
                  pot_grad_beta_ao_su_pbe_ot(1,1,istate),size(pot_grad_beta_ao_su_pbe_ot,1))
   enddo
   
 call wall_time(wall_2)

END_PROVIDER
