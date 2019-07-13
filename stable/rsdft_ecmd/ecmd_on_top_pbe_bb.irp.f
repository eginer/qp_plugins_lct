 subroutine give_nabla2_n2_beta_beta_hf_at_r1(r1,n2_deriv2)
 implicit none
 double precision, intent(in)  :: r1(3)
 double precision, intent(out) :: n2_deriv2
 integer :: i,j

 double precision :: mos_array_r1(mo_num)
 double precision :: nabla_2_at_r_mo(mo_num,mo_num)
 call give_all_mos_at_r(r1,mos_array_r1)
 call give_nabla_2_at_r_mo(r1,nabla_2_at_r_mo)


 n2_deriv2 = 0.d0

 do j = 1, elec_beta_num
  do i = 1, elec_beta_num
   n2_deriv2 += 0.3333333333333d0 * ( mos_array_r1(i)**2 * nabla_2_at_r_mo(j,j) - mos_array_r1(i)*mos_array_r1(j) * nabla_2_at_r_mo(i,j) )
  enddo
 enddo

 n2_deriv2 = 0.5d0*n2_deriv2
end

BEGIN_PROVIDER [double precision, nabla2_n2_hf_beta_beta, (n_points_final_grid)]
 implicit none
 integer :: i,j,i_point

 nabla2_n2_hf_beta_beta = 0.d0

 do i_point = 1, n_points_final_grid
  do j = 1, elec_beta_num
   do i = 1, elec_beta_num
    nabla2_n2_hf_beta_beta(i_point) += 0.3333333333333d0 * ( mos_in_r_array(i,i_point)**2 * mos_nabla_2_in_r_array_2(j,j,i_point) - mos_in_r_array(i,i_point)*mos_in_r_array(j,i_point) * mos_nabla_2_in_r_array_2(i,j,i_point) )
   enddo
  enddo
 enddo

 nabla2_n2_hf_beta_beta = 0.5d0 * nabla2_n2_hf_beta_beta

END_PROVIDER 

BEGIN_PROVIDER [double precision, Ecmd_pbe_n2_hf_bb, (N_states)]
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

 Ecmd_pbe_n2_hf_bb = 0.d0
 do i_state = 1, N_states
  do i_point = 1, n_points_final_grid
   weight = final_weight_at_r_vector(i_point)
   if(dabs(e_c_pbe_on_grid(2,i_point,i_state)).lt.1.d-13.or.dabs(nabla2_n2_hf_beta_beta(i_point)).lt.1.d-13)cycle
   beta_tmp = e_c_pbe_on_grid(2,i_point,i_state)/(constant * nabla2_n2_hf_beta_beta(i_point) * extrapol ) 
   Ecmd_pbe_n2_hf_bb += e_c_pbe_on_grid(2,i_point,i_state) / (1.d0 + beta_tmp * mu_erf_dft **5) * weight 
  enddo
 enddo

END_PROVIDER 

