 BEGIN_PROVIDER[double precision, Vijkl_eff_int_hf_alpha_beta, (mo_num,mo_num,elec_alpha_num,elec_beta_num)]
 implicit none
 BEGIN_DOC
 ! blablablabla
 END_DOC
 integer :: i,j,k,l
 double precision :: get_two_e_integral

 do l = 1, elec_beta_num
  do k = 1, elec_alpha_num
   do j = 1,mo_num
    do i = 1,mo_num
     Vijkl_eff_int_hf_alpha_beta(i,j,k,l) = get_two_e_integral(i,j,k,l,mo_integrals_map)
    enddo
   enddo
  enddo
 enddo 

 END_PROVIDER

 subroutine give_f_alpha_beta_hf_at_r1_r2(r1,r2,f_hf_alpha_beta)
 implicit none
 double precision, intent(in) :: r1(3),r2(3)
 double precision, intent(out) :: f_hf_alpha_beta 
 integer :: i,j,k,l

 double precision :: mos_array_r1(mo_num),mos_array_r2(mo_num)

 call give_all_mos_at_r(r1,mos_array_r1)
 call give_all_mos_at_r(r2,mos_array_r2)
 f_hf_alpha_beta = 0.d0

 do l = 1, elec_beta_num 
  do k = 1, elec_alpha_num
   do j = 1,mo_num
    do i = 1,mo_num
     f_hf_alpha_beta += Vijkl_eff_int_hf_alpha_beta(i,j,k,l)*mos_array_r1(k)*mos_array_r2(l)*mos_array_r1(i)*mos_array_r2(j)
    enddo
   enddo
  enddo
 enddo

 end



 subroutine give_f_alpha_beta_hf_at_r(r,f_hf_alpha_beta)
 implicit none
 double precision, intent(in) :: r(3)
 double precision, intent(out) :: f_hf_alpha_beta 
 double precision :: get_two_e_integral 
 integer :: i,j,k,l

 double precision :: mos_array_r(mo_num)

 call give_all_mos_at_r(r,mos_array_r)
 f_hf_alpha_beta  = 0.d0

 do l = 1, elec_beta_num 
  do k = 1, elec_alpha_num
   do i = 1,mo_num
    f_hf_alpha_beta += get_two_e_integral(i,l,k,l,mo_integrals_map)*mos_array_r(k)*mos_array_r(i)
   enddo
  enddo
 enddo


 end



 subroutine give_n2_alpha_beta_hf_at_r1_r2(r1,r2,n2_hf_alpha_beta)
 implicit none
 double precision, intent(in) :: r1(3),r2(3)
 double precision, intent(out) :: n2_hf_alpha_beta  
 integer :: i,j

 double precision :: mos_array_r1(mo_num),mos_array_r2(mo_num)

 call give_all_mos_at_r(r1,mos_array_r1)
 call give_all_mos_at_r(r2,mos_array_r2)
 n2_hf_alpha_beta = 0.d0

 do j = 1, elec_beta_num 
  do i = 1, elec_alpha_num
   n2_hf_alpha_beta += mos_array_r1(i)**2 *mos_array_r2(j)**2
  enddo
 enddo

 end


 subroutine give_n2_alpha_beta_hf_at_r(r,n2_hf_alpha_beta)
 implicit none
 double precision, intent(in) :: r(3)
 double precision, intent(out) :: n2_hf_alpha_beta 
 integer :: i,j

 double precision :: mos_array_r(mo_num)

 call give_all_mos_at_r(r,mos_array_r)
 n2_hf_alpha_beta = 0.d0

 do j = 1, elec_beta_num 
  do i = 1, elec_alpha_num
   n2_hf_alpha_beta += mos_array_r(i)**2 
  enddo
 enddo

 end


 
 subroutine give_n2_alpha_beta_hf_at_r1_r12(r1,r12,n2_hf_alpha_beta,n2_0,n2_deriv2,n2_deriv4)
 implicit none
 double precision, intent(in)  :: r1(3)
 double precision, intent(in)  :: r12
 double precision, intent(out) :: n2_0
 double precision, intent(out) :: n2_hf_alpha_beta 
 double precision, intent(out) :: n2_deriv2,n2_deriv4
 integer :: i,j

 double precision :: mos_array_r1(mo_num)
 double precision :: nabla_2_at_r_mo(mo_num,mo_num)
 double precision :: nabla_4_at_r_mo(mo_num,mo_num)
 call give_all_mos_at_r(r1,mos_array_r1)
 call give_nabla_2_at_r_mo(r1,nabla_2_at_r_mo)
 call give_nabla_4_at_r_mo(r1,nabla_4_at_r_mo)
 
 n2_hf_alpha_beta = 0.d0

 n2_0      = 0.d0 
 n2_deriv2 = 0.d0
 n2_deriv4 = 0.d0

 do j = 1, elec_beta_num
  do i = 1, elec_alpha_num
   n2_0      += mos_array_r1(i)**2 * mos_array_r1(j)**2 
   n2_deriv2 += 0.3333333333333d0 * mos_array_r1(i)**2 * nabla_2_at_r_mo(j,j)
   n2_deriv4 += 0.2d0             * mos_array_r1(i)**2 * nabla_4_at_r_mo(j,j) 
  enddo
 enddo
   
 n2_hf_alpha_beta = n2_0 + (0.5d0)*n2_deriv2*(r12**2.d0)+(fact_inv(4))*n2_deriv4*(r12**4.d0)
 end


 subroutine give_f_paper_alpha_beta_hf_at_r1_r12(r1,r12,f_paper_hf_alpha_beta,n2_0,n2_deriv2,n2_deriv4)
 implicit none
 double precision, intent(in)  :: r1(3)
 double precision, intent(in)  :: r12
 double precision, intent(out) :: f_paper_hf_alpha_beta
 double precision, intent(out) :: n2_0,n2_deriv2,n2_deriv4
 integer :: i,j,k,l

 double precision :: mos_array_r1(mo_num)
 double precision :: nabla_2_at_r_mo(mo_num,mo_num)
 double precision :: nabla_4_at_r_mo(mo_num,mo_num)
 call give_all_mos_at_r(r1,mos_array_r1)
 call give_nabla_2_at_r_mo(r1,nabla_2_at_r_mo)
 call give_nabla_4_at_r_mo(r1,nabla_4_at_r_mo)
 
 f_paper_hf_alpha_beta = 0.d0

 n2_0      = 0.d0  
 n2_deriv2 = 0.d0
 n2_deriv4 = 0.d0

 do l = 1, elec_beta_num 
  do k = 1, elec_alpha_num
   do j = 1,mo_num
    do i = 1,mo_num
     n2_0      += Vijkl_eff_int_hf_alpha_beta(i,j,k,l) *  mos_array_r1(i)* mos_array_r1(k) * mos_array_r1(j)* mos_array_r1(l) 
     n2_deriv2 += 0.3333333333333d0 * Vijkl_eff_int_hf_alpha_beta(i,j,k,l) *  mos_array_r1(i)* mos_array_r1(k) * nabla_2_at_r_mo(j,l) 
     n2_deriv4 += 0.2d0             * Vijkl_eff_int_hf_alpha_beta(i,j,k,l) *  mos_array_r1(i)* mos_array_r1(k) * nabla_4_at_r_mo(j,l) 
    enddo
   enddo
  enddo
 enddo
 
 f_paper_hf_alpha_beta = n2_0 + (0.5d0)*n2_deriv2*(r12**2.d0)+(fact_inv(4))*n2_deriv4*(r12**4.d0)
 end

 subroutine give_eff_inter_alpha_beta_hf_at_r1_r12(r1,r12,f_paper_hf,n2_hf)
 implicit none
 double precision, intent(in)  :: r1(3)
 double precision, intent(in)  :: r12
 double precision, intent(out) :: f_paper_hf
 double precision, intent(out) :: n2_hf

 double precision :: f_0,f_deriv2,f_deriv4,n2_0,n2_deriv2,n2_deriv4
 
 call give_f_paper_alpha_beta_hf_at_r1_r12(r1,r12,f_paper_hf,f_0,f_deriv2,f_deriv4)  
 call give_n2_alpha_beta_hf_at_r1_r12(r1,r12,n2_hf,n2_0,n2_deriv2,n2_deriv4) 

 end
