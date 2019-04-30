 
 subroutine give_h_sr_component_at_r1_r12(r1,h_0,h_deriv2,h_deriv4)
 implicit none
 double precision, intent(in)  :: r1(3)
 double precision, intent(out) :: h_0(n_states)
 double precision, intent(out) :: h_deriv2(n_states),h_deriv4(n_states)
 integer :: i,j,istate

 double precision :: dm_a(n_states),dm_b(n_states)
 double precision :: aos_array_r1(ao_num)
 double precision :: nabla_2_at_r(ao_num,ao_num)
 double precision :: nabla_2_cop(3,ao_num,ao_num)
 double precision :: nabla_4_at_r(ao_num,ao_num)

 call give_all_aos_at_r(r1,aos_array_r1)
 call give_nabla_2_at_r(r1,nabla_2_cop,nabla_2_at_r)
 call give_nabla_4_at_r(r1,nabla_4_at_r)
 call dm_dft_alpha_beta_at_r(r1,dm_a,dm_b)

 h_0      = 0.d0 
 h_deriv2 = 0.d0
 h_deriv4 = 0.d0
 do istate = 1, n_states
  do j = 1, ao_num 
   do i = 1, ao_num 
   !h_0(istate)      += 2.d0*(dm_a(istate)+dm_b(istate))*((one_e_dm_alpha_ao_for_dft(i,j,istate)+one_e_dm_beta_ao_for_dft(i,j,istate))* aos_array_r1(i)* aos_array_r1(j)) 
   !h_deriv2(istate) += 2.d0*(dm_a(istate)+dm_b(istate))*(0.3333333333333d0 *((one_e_dm_alpha_ao_for_dft(i,j,istate)+one_e_dm_beta_ao_for_dft(i,j,istate)) * nabla_2_at_r(i,j)))
   !h_deriv4(istate) += 2.d0*(dm_a(istate)+dm_b(istate))*( 0.2d0             *  ((one_e_dm_alpha_ao_for_dft(i,j,istate)+one_e_dm_beta_ao_for_dft(i,j,istate)) * nabla_4_at_r(i,j))) 
    h_0(istate)      += (dm_a(istate)+dm_b(istate))*((one_e_dm_alpha_ao_for_dft(i,j,istate)+one_e_dm_beta_ao_for_dft(i,j,istate))* aos_array_r1(i)* aos_array_r1(j)) 
    h_deriv2(istate) += (dm_a(istate)+dm_b(istate))*(0.3333333333333d0 *((one_e_dm_alpha_ao_for_dft(i,j,istate)+one_e_dm_beta_ao_for_dft(i,j,istate)) * nabla_2_at_r(i,j)))
    h_deriv4(istate) += (dm_a(istate)+dm_b(istate))*( 0.2d0             *  ((one_e_dm_alpha_ao_for_dft(i,j,istate)+one_e_dm_beta_ao_for_dft(i,j,istate)) * nabla_4_at_r(i,j))) 
   enddo
  enddo
 enddo  
 end

 BEGIN_PROVIDER[double precision, e_h_sr_exp0,(n_states)]
&BEGIN_PROVIDER[double precision, e_h_sr_exp1,(n_states)]
&BEGIN_PROVIDER[double precision, e_h_sr_exp2,(n_states)]
&BEGIN_PROVIDER[double precision, e_h_sr_tot,(n_states)]
 implicit none
 include 'utils/constants.include.F'
 integer :: i,istate
 double precision, allocatable :: r(:),h_0(:),h_deriv2(:),h_deriv4(:),int_h_0(:),int_h_1(:),int_h_2(:)  
 allocate(r(3),h_0(n_states),h_deriv2(n_states),h_deriv4(n_states),int_h_0(n_states),int_h_1(n_states),int_h_2(n_states))

 int_h_0 = 0.d0
 int_h_1 = 0.d0
 int_h_2 = 0.d0

 do i = 1, n_points_final_grid
  r(1) = final_grid_points(1,i)
  r(2) = final_grid_points(2,i)
  r(3) = final_grid_points(3,i)
  call give_h_sr_component_at_r1_r12(r,h_0,h_deriv2,h_deriv4) 
  int_h_0 += h_0 * final_weight_at_r_vector(i)
  int_h_1 += h_deriv2 * final_weight_at_r_vector(i)
  int_h_2 += h_deriv4 * final_weight_at_r_vector(i)

 enddo

 e_h_sr_exp0 = 0.5d0 * pi * int_h_0 
 e_h_sr_exp1 = 3.d0 * pi * int_h_1 /16.d0
 e_h_sr_exp2 = 15.d0 * pi * int_h_2 /576.d0 

 e_h_sr_tot = e_h_sr_exp0/mu_erf **2.d0 + e_h_sr_exp1/ mu_erf **4.d0 + e_h_sr_exp2/ mu_erf **6.d0

 END_PROVIDER


 BEGIN_PROVIDER[double precision, pade_1_h,(n_states)]
&BEGIN_PROVIDER[double precision, pade_2_h,(n_states)]
&BEGIN_PROVIDER[double precision, pade_3_h,(n_states)]
 implicit none
 double precision :: xpad
 double precision :: c2(n_states),c4(n_states),c6(n_states)

 xpad = 1/mu_erf**2.d0
 c2 = e_h_sr_exp0
 c4 = e_h_sr_exp1
 c6 = e_h_sr_exp2

 pade_1_h = (c2 * xpad)/(1.d0 - c4*xpad/c2 + (c4**2.d0-c2*c6)*xpad**2.d0/c2**2.d0 )
 pade_2_h = (c2*xpad +(c4**2.d0-c2*c6)*xpad**2.d0/c4 )/(1.d0 - c6*xpad/c4 )  
 pade_3_h = (c2*xpad +(c4*((c4**2.d0)-2.d0*c2*c6)*(xpad**2.d0))/((c4**2.d0)-c2*c6) )/(1.d0 - c4*c6*xpad/((c4**2.d0)-c2*c6)+ (c6**2.d0)*(xpad**2.d0)/((c4**2.d0)-c2*c6))  
 END_PROVIDER


 BEGIN_PROVIDER[double precision, e_h_local,(n_states)]
&BEGIN_PROVIDER[double precision, e_h_pade_local_1,(n_states)]
&BEGIN_PROVIDER[double precision, e_h_pade_local_2,(n_states)]
&BEGIN_PROVIDER[double precision, e_h_pade_local_3,(n_states)]
 implicit none
 include 'utils/constants.include.F'
 integer :: i,istate
 double precision :: xpad  
 double precision, allocatable :: r(:),h_0(:),h_deriv2(:),h_deriv4(:),c2(:),c4(:),c6(:) 
 allocate(r(3),h_0(n_states),h_deriv2(n_states),h_deriv4(n_states),c2(n_states),c4(n_states),c6(n_states))

 e_h_local = 0.d0
 e_h_pade_local_1 = 0.d0
 e_h_pade_local_2 = 0.d0
 e_h_pade_local_3 = 0.d0
 xpad = 1/mu_erf**2.d0

 do i = 1, n_points_final_grid
  r(1) = final_grid_points(1,i)
  r(2) = final_grid_points(2,i)
  r(3) = final_grid_points(3,i)
  call give_h_sr_component_at_r1_r12(r,h_0,h_deriv2,h_deriv4)
  c2 = 0.5d0 * pi * h_0
  c4 = 3.d0 * pi * h_deriv2 /16.d0
  c6 = 15.d0 * pi * h_deriv4 / 576.d0
  e_h_local += (c2*xpad + c4*xpad**2 + c6*xpad**3) * final_weight_at_r_vector(i)
  e_h_pade_local_1 += ((c2 * xpad)/(1.d0 - c4*xpad/c2 + (c4**2.d0-c2*c6)*xpad**2.d0/c2**2.d0 )) * final_weight_at_r_vector(i)
  e_h_pade_local_2 += ((c2*xpad +(c4**2.d0-c2*c6)*xpad**2.d0/c4 )/(1.d0 -c6*xpad/c4)) * final_weight_at_r_vector(i)
  e_h_pade_local_3 += ((c2*xpad +(c4*((c4**2.d0)-2.d0*c2*c6)*(xpad**2.d0))/((c4**2.d0)-c2*c6) )/(1.d0 - c4*c6*xpad/((c4**2.d0)-c2*c6)+ (c6**2.d0)*(xpad**2.d0)/((c4**2.d0)-c2*c6)))* final_weight_at_r_vector(i) 
 enddo


 END_PROVIDER
