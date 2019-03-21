 
 subroutine give_n2_alpha_beta_component_at_r1_r12(r1,n2_0,n2_deriv2,n2_deriv4)
 implicit none
 double precision, intent(in)  :: r1(3)
 double precision, intent(out) :: n2_0(n_states)
 double precision, intent(out) :: n2_deriv2(n_states),n2_deriv4(n_states)
 integer :: i,j,k,l,istate

 double precision :: mos_array_r1(mo_num)
 double precision :: nabla_2_at_r_mo(mo_num,mo_num)
 double precision :: nabla_4_at_r_mo(mo_num,mo_num)
 call give_all_mos_at_r(r1,mos_array_r1)
 call give_nabla_2_at_r_mo(r1,nabla_2_at_r_mo)
 call give_nabla_4_at_r_mo(r1,nabla_4_at_r_mo)
 
 n2_0      = 0.d0 
 n2_deriv2 = 0.d0
 n2_deriv4 = 0.d0
 do istate = 1, n_states
  do l = 1, mo_num       
   do k = 1, mo_num        
    do j = 1, mo_num 
     do i = 1, mo_num 
      n2_0(istate)      += 2.d0*( two_bod_alpha_beta_mo_physicist(i,j,k,l,istate)* mos_array_r1(i)* mos_array_r1(k) * mos_array_r1(j) * mos_array_r1(l))
      n2_deriv2(istate) += 2.d0*(0.3333333333333d0 * two_bod_alpha_beta_mo_physicist(i,j,k,l,istate) * mos_array_r1(i)* mos_array_r1(k) * nabla_2_at_r_mo(j,l))
      n2_deriv4(istate) += 2.d0*( 0.2d0             * two_bod_alpha_beta_mo_physicist(i,j,k,l,istate) * mos_array_r1(i)* mos_array_r1(k) * nabla_4_at_r_mo(j,l)) 
     enddo
    enddo
   enddo
  enddo
 enddo  
 end

 BEGIN_PROVIDER[double precision, e_hx_sr_exp0,(n_states)]
&BEGIN_PROVIDER[double precision, e_hx_sr_exp1,(n_states)]
&BEGIN_PROVIDER[double precision, e_hx_sr_exp2,(n_states)]
&BEGIN_PROVIDER[double precision, e_hx_sr_tot,(n_states)]
 implicit none
 include 'utils/constants.include.F'
 integer :: i,istate
 double precision, allocatable :: r(:),n2_0(:),n2_deriv2(:),n2_deriv4(:),int_n_2_0(:),int_n_2_1(:),int_n_2_2(:)  
 allocate(r(3),n2_0(n_states),n2_deriv2(n_states),n2_deriv4(n_states),int_n_2_0(n_states),int_n_2_1(n_states),int_n_2_2(n_states))

 int_n_2_0 = 0.d0
 int_n_2_1 = 0.d0
 int_n_2_2 = 0.d0

 do i = 1, n_points_final_grid
  r(1) = final_grid_points(1,i)
  r(2) = final_grid_points(2,i)
  r(3) = final_grid_points(3,i)
  call give_n2_alpha_beta_component_at_r1_r12(r,n2_0,n2_deriv2,n2_deriv4) 
  int_n_2_0 += n2_0 * final_weight_at_r_vector(i)
  int_n_2_1 += n2_deriv2 * final_weight_at_r_vector(i)
  int_n_2_2 += n2_deriv4 * final_weight_at_r_vector(i)
 !do istate = 1, n_states
 ! int_n_2_0(istate) += n2_0(istate) * final_weight_at_r_vector(i)
 ! int_n_2_1(istate) += n2_deriv2(istate) * final_weight_at_r_vector(i)
 ! int_n_2_2(istate) += n2_deriv4(istate) * final_weight_at_r_vector(i)
 !enddo

 enddo

!e_hx_sr_exp0 = 0.5d0 * pi * int_n_2_0 /(mu_erf **2.d0) 
!e_hx_sr_exp1 = 3.d0 * pi * int_n_2_1 /(16.d0 * mu_erf **4.d0)
!e_hx_sr_exp2 = 15.d0 * pi * int_n_2_2 /(576.d0 * mu_erf **6.d0)


 e_hx_sr_exp0 = 0.5d0 * pi * int_n_2_0 
 e_hx_sr_exp1 = 3.d0 * pi * int_n_2_1 /16.d0
 e_hx_sr_exp2 = 15.d0 * pi * int_n_2_2 /576.d0 

 e_hx_sr_tot = e_hx_sr_exp0/mu_erf **2.d0 + e_hx_sr_exp1/ mu_erf **4.d0 + e_hx_sr_exp2/ mu_erf **6.d0

 END_PROVIDER


 BEGIN_PROVIDER[double precision, pade_1,(n_states)]
&BEGIN_PROVIDER[double precision, pade_2,(n_states)]
&BEGIN_PROVIDER[double precision, pade_3,(n_states)]
 implicit none
 double precision :: xpad
 double precision :: c2(n_states),c4(n_states),c6(n_states)

 xpad = 1/mu_erf**2.d0
 c2 = e_hx_sr_exp0
 c4 = e_hx_sr_exp1
 c6 = e_hx_sr_exp2

 pade_1 = (c2 * xpad)/(1.d0 - c4*xpad/c2 + (c4**2.d0-c2*c6)*xpad**2.d0/c2**2.d0 )
 pade_2 = (c2*xpad +(c4**2.d0-c2*c6)*xpad**2.d0/c4 )/(1.d0 - c6*xpad/c4 )  
 pade_3 = (c2*xpad +(c4*((c4**2.d0)-2.d0*c2*c6)*(xpad**2.d0))/((c4**2.d0)-c2*c6) )/(1.d0 - c4*c6*xpad/((c4**2.d0)-c2*c6)+ (c6**2.d0)*(xpad**2.d0)/((c4**2.d0)-c2*c6))  
 END_PROVIDER


 BEGIN_PROVIDER[double precision, e_hx_local,(n_states)]
&BEGIN_PROVIDER[double precision, e_hx_pade_local_1,(n_states)]
&BEGIN_PROVIDER[double precision, e_hx_pade_local_2,(n_states)]
&BEGIN_PROVIDER[double precision, e_hx_pade_local_3,(n_states)]
 implicit none
 include 'utils/constants.include.F'
 integer :: i,istate
 double precision :: xpad  
 double precision, allocatable :: r(:),n2_0(:),n2_deriv2(:),n2_deriv4(:),c2(:),c4(:),c6(:) 
 allocate(r(3),n2_0(n_states),n2_deriv2(n_states),n2_deriv4(n_states),c2(n_states),c4(n_states),c6(n_states))

 e_hx_local = 0.d0
 e_hx_pade_local_1 = 0.d0
 e_hx_pade_local_2 = 0.d0
 e_hx_pade_local_3 = 0.d0
 xpad = 1/mu_erf**2.d0

 do i = 1, n_points_final_grid
  r(1) = final_grid_points(1,i)
  r(2) = final_grid_points(2,i)
  r(3) = final_grid_points(3,i)
  call give_n2_alpha_beta_component_at_r1_r12(r,n2_0,n2_deriv2,n2_deriv4) 
  c2 = 0.5d0 * pi * n2_0
  c4 = 3.d0 * pi * n2_deriv2 /16.d0
  c6 = 15.d0 * pi * n2_deriv4 / 576.d0
  e_hx_local += (c2*xpad + c4*xpad**2 + c6*xpad**3) * final_weight_at_r_vector(i)
  e_hx_pade_local_1 += ((c2 * xpad)/(1.d0 - c4*xpad/c2 + (c4**2.d0-c2*c6)*xpad**2.d0/c2**2.d0 )) * final_weight_at_r_vector(i)
  e_hx_pade_local_2 += ((c2*xpad +(c4**2.d0-c2*c6)*xpad**2.d0/c4 )/(1.d0 -c6*xpad/c4)) * final_weight_at_r_vector(i)
  e_hx_pade_local_3 += ((c2*xpad +(c4*((c4**2.d0)-2.d0*c2*c6)*(xpad**2.d0))/((c4**2.d0)-c2*c6) )/(1.d0 - c4*c6*xpad/((c4**2.d0)-c2*c6)+ (c6**2.d0)*(xpad**2.d0)/((c4**2.d0)-c2*c6)))* final_weight_at_r_vector(i) 
 enddo


 END_PROVIDER
