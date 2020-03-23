 subroutine give_nabla2_n2_alpha_alpha_hf_at_r1(r1,n2_deriv2)
 implicit none
 double precision, intent(in)  :: r1(3)
 double precision, intent(out) :: n2_deriv2
 integer :: i,j

 double precision :: mos_array_r1(mo_num)
 double precision :: nabla_2_at_r_mo(mo_num,mo_num)
 call give_all_mos_at_r(r1,mos_array_r1)
 call give_nabla_2_at_r_mo(r1,nabla_2_at_r_mo)


 n2_deriv2 = 0.d0

 do j = 1, elec_alpha_num
  do i = 1, elec_alpha_num
   n2_deriv2 += 0.3333333333333d0 * ( mos_array_r1(i)**2 * nabla_2_at_r_mo(j,j) - mos_array_r1(i)*mos_array_r1(j) * nabla_2_at_r_mo(i,j) )
  enddo
 enddo

 n2_deriv2 = 0.5d0*n2_deriv2
end

BEGIN_PROVIDER [double precision, nabla2_n2_hf_alpha_alpha, (n_points_final_grid)]
 implicit none
 integer :: i,j,i_point

 nabla2_n2_hf_alpha_alpha = 0.d0

 do i_point = 1, n_points_final_grid
  do j = 1, elec_alpha_num
   do i = 1, elec_alpha_num
    nabla2_n2_hf_alpha_alpha(i_point) += 0.3333333333333d0 * ( mos_in_r_array(i,i_point)**2 * mos_nabla_2_in_r_array_2(j,j,i_point) - mos_in_r_array(i,i_point)*mos_in_r_array(j,i_point) * mos_nabla_2_in_r_array_2(i,j,i_point) )
   enddo
  enddo
 enddo

 nabla2_n2_hf_alpha_alpha = 0.5d0 * nabla2_n2_hf_alpha_alpha

END_PROVIDER 

BEGIN_PROVIDER [double precision, Ecmd_pbe_n2_hf_aa, (N_states)]
 implicit none
 integer :: i_state,i_point
 double precision :: weight
 double precision :: constant,beta_tmp,extrapol
 constant = (-3.d0 + 2.d0 * dsqrt(2.d0)) * dsqrt(dacos(-1.d0)) / (10.d0*dsqrt(2.d0) ) 
 if(mu_erf_dft.lt.1.d-10)then
  extrapol = 1.d-20
 else
  extrapol = 1.d0/(1.d0 + 2.d0 / (3.d0 * dsqrt(dacos(-1.d0))*mu_erf_dft )) ! extrapolation factor for the two body density 
 endif
 Ecmd_pbe_n2_hf_aa = 0.d0
 do i_state = 1, N_states
  do i_point = 1, n_points_final_grid
   weight = final_weight_at_r_vector(i_point)
   if(dabs(e_c_pbe_on_grid(1,i_point,i_state)).lt.1.d-13.or.dabs(nabla2_n2_hf_alpha_alpha(i_point)).lt.1.d-13)cycle

   beta_tmp = e_c_pbe_on_grid(1,i_point,i_state)/(constant * nabla2_n2_hf_alpha_alpha(i_point) * extrapol ) 
   Ecmd_pbe_n2_hf_aa += e_c_pbe_on_grid(1,i_point,i_state) / (1.d0 + beta_tmp * mu_erf_dft **5) * weight 
  enddo
 enddo

END_PROVIDER 


BEGIN_PROVIDER [double precision, e_c_pbe_on_grid, (3,n_points_final_grid,n_states)]
 implicit none
 BEGIN_DOC
 ! e_c_pbe_on_grid(1,i,istate) = PBE correlation energy at grid point i for the alpha density of the state i_state
 !
 ! e_c_pbe_on_grid(2,i,istate) = PBE correlation energy at grid point i for the beta density  of the state i_state
 !
 ! e_c_pbe_on_grid(3,i,istate) = PBE correlation energy at grid point i for the TOTAL density of the state i_state
 END_DOC
 integer :: i_point,i_state,j
 double precision :: rhoc,rhoo,sigmacc,sigmaco,sigmaoo
 double precision :: rho_a,rho_b,grad_rho_ab,grad_rho_a2,grad_rho_b2
  
 do i_state = 1, N_states
  do i_point = 1, n_points_final_grid
   ! Total density 
   rho_a = one_e_dm_and_grad_alpha_in_r(4,i_point,i_state)
   rho_b = one_e_dm_and_grad_beta_in_r(4,i_point,i_state)
   call rho_ab_to_rho_oc(rho_a,rho_b,rhoo,rhoc)
   grad_rho_ab = 0.d0
   do j = 1, 3
    grad_rho_ab += one_e_dm_and_grad_alpha_in_r(j,i_point,i_state ) * one_e_dm_and_grad_beta_in_r(j,i_point,i_state )
   enddo
   grad_rho_a2 = one_e_grad_2_dm_alpha_at_r(i_point,i_state)
   grad_rho_b2 = one_e_grad_2_dm_beta_at_r(i_point,i_state)
   call grad_rho_ab_to_grad_rho_oc(grad_rho_a2,grad_rho_b2,grad_rho_ab,sigmaoo,sigmacc,sigmaco)
   call ec_pbe_only(0.d0,rhoc,rhoo,sigmacc,sigmaco,sigmaoo,e_c_pbe_on_grid(3,i_point,i_state))

   ! alpha density 
   rho_b = 0.d0
   call rho_ab_to_rho_oc(rho_a,rho_b,rhoo,rhoc)
   grad_rho_ab = 0.d0
   grad_rho_b2 = 0.d0
   call grad_rho_ab_to_grad_rho_oc(grad_rho_a2,grad_rho_b2,grad_rho_ab,sigmaoo,sigmacc,sigmaco)
   call ec_pbe_only(0.d0,rhoc,rhoo,sigmacc,sigmaco,sigmaoo,e_c_pbe_on_grid(1,i_point,i_state))
  
  ! beta density 
   rho_a = one_e_dm_and_grad_beta_in_r(4,i_point,i_state)
   rho_b = 0.d0
   call rho_ab_to_rho_oc(rho_a,rho_b,rhoo,rhoc)
   grad_rho_ab = 0.d0
   grad_rho_b2 = 0.d0
   grad_rho_a2 = one_e_grad_2_dm_beta_at_r(i_point,i_state)
   call grad_rho_ab_to_grad_rho_oc(grad_rho_a2,grad_rho_b2,grad_rho_ab,sigmaoo,sigmacc,sigmaco)
   call ec_pbe_only(0.d0,rhoc,rhoo,sigmacc,sigmaco,sigmaoo,e_c_pbe_on_grid(2,i_point,i_state))
   
  enddo
 enddo

END_PROVIDER 

 BEGIN_PROVIDER [double precision, Energy_c_md_on_top_spin_pol]
 BEGIN_DOC
! new ECMD functional for spin polarized electrons
 END_DOC
 implicit none
 integer :: i_point
 double precision :: weight
 double precision :: constant
 constant = (-3.d0 + 2.d0 * dsqrt(2.d0)) * dsqrt(dacos(-1.d0)) / (10.d0*dsqrt(2.d0) * mu_erf_dft**5)
 Energy_c_md_on_top_spin_pol = 0.d0
 do i_point = 1, n_points_final_grid
  weight = final_weight_at_r_vector(i_point)
  Energy_c_md_on_top_spin_pol += nabla2_n2_hf_alpha_alpha(i_point) * weight 
 enddo
 ! extrapolation of the on top
 Energy_c_md_on_top_spin_pol *= constant / (1.d0 + 2.d0 / (3.d0 * dsqrt(dacos(-1.d0))*mu_erf_dft ))
 END_PROVIDER 

