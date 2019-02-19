 BEGIN_PROVIDER[double precision, Vijkl_eff_int_hf_alpha_alpha, (mo_num,mo_num,elec_alpha_num,elec_alpha_num)]
 implicit none
 BEGIN_DOC
 ! blablablabla
 END_DOC
 integer :: i,j,k,l
 double precision :: get_two_e_integral

 do l = 1, elec_alpha_num
  do k = 1, elec_alpha_num
   do j = 1,mo_num
    do i = 1,mo_num
     Vijkl_eff_int_hf_alpha_alpha(i,j,k,l) = get_two_e_integral(i,j,k,l,mo_integrals_map)
    enddo
   enddo
  enddo
 enddo 

 END_PROVIDER

 subroutine give_f_alpha_alpha_hf_at_r1_r2(r1,r2,f_hf_alpha_alpha)
 implicit none
 double precision, intent(in) :: r1(3),r2(3)
 double precision, intent(out) :: f_hf_alpha_alpha 
 integer :: i,j,k,l

 double precision :: mos_array_r1(mo_num),mos_array_r2(mo_num)

 call give_all_mos_at_r(r1,mos_array_r1)
 call give_all_mos_at_r(r2,mos_array_r2)
 f_hf_alpha_alpha = 0.d0

 do l = 1, elec_alpha_num 
  do k = 1, elec_alpha_num
   do j = 1,mo_num
    do i = 1,mo_num
     f_hf_alpha_alpha += Vijkl_eff_int_hf_alpha_alpha(i,j,k,l)*(mos_array_r1(k)*mos_array_r2(l)-mos_array_r2(k)*mos_array_r1(l))*mos_array_r1(i)*mos_array_r2(j)
    enddo
   enddo
  enddo
 enddo

 end



 subroutine give_f_alpha_alpha_hf_at_r(r,f_hf_alpha_alpha)
 implicit none
 double precision, intent(in) :: r(3)
 double precision, intent(out) :: f_hf_alpha_alpha 
 double precision :: get_two_e_integral 
 integer :: i,j,k,l

 double precision :: mos_array_r(mo_num)

 call give_all_mos_at_r(r,mos_array_r)
 f_hf_alpha_alpha = 0.d0

 do l = 1, elec_alpha_num 
  do k = 1, elec_alpha_num
   do i = 1,mo_num
    f_hf_alpha_alpha += get_two_e_integral(i,l,k,l,mo_integrals_map)*mos_array_r(k)*mos_array_r(i)-get_two_e_integral(i,k,k,l,mo_integrals_map)*mos_array_r(l)*mos_array_r(i)
   enddo
  enddo
 enddo


 end



 subroutine give_n2_alpha_alpha_hf_at_r1_r2(r1,r2,n2_hf_alpha_alpha)
 implicit none
 double precision, intent(in) :: r1(3),r2(3)
 double precision, intent(out) :: n2_hf_alpha_alpha 
 integer :: i,j

 double precision :: mos_array_r1(mo_num),mos_array_r2(mo_num)

 call give_all_mos_at_r(r1,mos_array_r1)
 call give_all_mos_at_r(r2,mos_array_r2)
 n2_hf_alpha_alpha = 0.d0

 do j = 1, elec_alpha_num 
  do i = 1, elec_alpha_num
   n2_hf_alpha_alpha += mos_array_r1(i)**2 *mos_array_r2(j)**2-mos_array_r1(i)*mos_array_r2(i)*mos_array_r1(j)*mos_array_r2(j)
  enddo
 enddo

 end


 subroutine give_n2_alpha_alpha_hf_at_r(r,n2_hf_alpha_alpha)
 implicit none
 double precision, intent(in) :: r(3)
 double precision, intent(out) :: n2_hf_alpha_alpha 
 integer :: i,j

 double precision :: mos_array_r(mo_num)

 call give_all_mos_at_r(r,mos_array_r)
 n2_hf_alpha_alpha = 0.d0

 do j = 1, elec_alpha_num 
  do i = 1, elec_alpha_num
   n2_hf_alpha_alpha += mos_array_r(i)**2 
  enddo
  n2_hf_alpha_alpha -= mos_array_r(j)**2
 enddo

 end


 
 subroutine give_n2_alpha_alpha_hf_at_r1_r12(r1,r12,n2_hf_alpha_alpha,n2_deriv2,n2_deriv4)
 implicit none
 double precision, intent(in)  :: r1(3)
 double precision, intent(in)  :: r12
 double precision, intent(out) :: n2_hf_alpha_alpha
 double precision, intent(out) :: n2_deriv2,n2_deriv4
 integer :: i,j

 double precision :: mos_array_r1(mo_num)
 double precision :: nabla_2_at_r_mo(mo_num,mo_num)
 double precision :: nabla_4_at_r_mo(mo_num,mo_num)
 call give_all_mos_at_r(r1,mos_array_r1)
 call give_nabla_2_at_r_mo(r1,nabla_2_at_r_mo)
 call give_nabla_4_at_r_mo(r1,nabla_4_at_r_mo)
 
 n2_hf_alpha_alpha = 0.d0
 
 n2_deriv2 = 0.d0
 n2_deriv4 = 0.d0

 do j = 1, elec_alpha_num
  do i = 1, elec_alpha_num
   n2_deriv2 += 0.3333333333333d0 * ( mos_array_r1(i)**2 * nabla_2_at_r_mo(j,j) - mos_array_r1(i)*mos_array_r1(j) * nabla_2_at_r_mo(i,j) )
   n2_deriv4 += 0.2d0             * ( mos_array_r1(i)**2 * nabla_4_at_r_mo(j,j) - mos_array_r1(i)*mos_array_r1(j) * nabla_4_at_r_mo(i,j) )
  enddo
 enddo
   
 n2_hf_alpha_alpha = (0.5d0)*n2_deriv2*(r12**2.d0)+(fact_inv(4))*n2_deriv4*(r12**4.d0)
 end


 subroutine give_f_paper_alpha_alpha_hf_at_r1_r12(r1,r12,f_paper_hf_alpha_alpha,n2_deriv2,n2_deriv4)
 implicit none
 double precision, intent(in)  :: r1(3)
 double precision, intent(in)  :: r12
 double precision, intent(out) :: f_paper_hf_alpha_alpha
 double precision, intent(out) :: n2_deriv2,n2_deriv4
 integer :: i,j,k,l

 double precision :: mos_array_r1(mo_num)
 double precision :: nabla_2_at_r_mo(mo_num,mo_num)
 double precision :: nabla_4_at_r_mo(mo_num,mo_num)
 call give_all_mos_at_r(r1,mos_array_r1)
 call give_nabla_2_at_r_mo(r1,nabla_2_at_r_mo)
 call give_nabla_4_at_r_mo(r1,nabla_4_at_r_mo)
 
 f_paper_hf_alpha_alpha = 0.d0
 
 n2_deriv2 = 0.d0
 n2_deriv4 = 0.d0

 do l = 1, elec_alpha_num 
  do k = 1, elec_alpha_num
   do j = 1,mo_num
    do i = 1,mo_num
     n2_deriv2 += 0.3333333333333d0 * Vijkl_eff_int_hf_alpha_alpha(i,j,k,l) * ( mos_array_r1(i)* mos_array_r1(k) * nabla_2_at_r_mo(j,l) - mos_array_r1(i)*mos_array_r1(l) * nabla_2_at_r_mo(j,k) )
     n2_deriv4 += 0.2d0             * Vijkl_eff_int_hf_alpha_alpha(i,j,k,l) * ( mos_array_r1(i)* mos_array_r1(k) * nabla_4_at_r_mo(j,l) - mos_array_r1(i)*mos_array_r1(l) * nabla_4_at_r_mo(j,k) ) 
    enddo
   enddo
  enddo
 enddo
 
 f_paper_hf_alpha_alpha = (0.5d0)*n2_deriv2*(r12**2.d0)+(fact_inv(4))*n2_deriv4*(r12**4.d0)
 end

 subroutine give_eff_inter_alpha_alpha_hf_at_r1_r12(r1,r12,f_paper_hf,n2_hf)
 implicit none
 double precision, intent(in)  :: r1(3)
 double precision, intent(in)  :: r12
 double precision, intent(out) :: f_paper_hf
 double precision, intent(out) :: n2_hf

 double precision :: f_deriv2,f_deriv4,n2_deriv2,n2_deriv4
 
 call give_f_paper_alpha_alpha_hf_at_r1_r12(r1,r12,f_paper_hf,f_deriv2,f_deriv4)  
 call give_n2_alpha_alpha_hf_at_r1_r12(r1,r12,n2_hf,n2_deriv2,n2_deriv4) 

 end


double precision  function spherical_av_wee_paper(r1,r12)
 implicit none
 include 'utils/constants.include.F'
 double precision, intent(in) :: r1(3),r12
 integer :: k
 double precision :: f_hf_alpha_alpha,r2(3),n2_hf_alpha_alpha
 spherical_av_wee_paper = 0.d0
 do k = 1,n_points_integration_angular
  r2(1)= r1(1)+r12*angular_quadrature_points(k,1)
  r2(2)= r1(2)+r12*angular_quadrature_points(k,2)
  r2(3)= r1(3)+r12*angular_quadrature_points(k,3)
  call give_f_alpha_alpha_hf_at_r1_r2(r1,r2,f_hf_alpha_alpha) 
  call give_n2_alpha_alpha_hf_at_r1_r2(r1,r2,n2_hf_alpha_alpha) 
  spherical_av_wee_paper += f_hf_alpha_alpha/n2_hf_alpha_alpha * weights_angular_points(k)
 enddo
 spherical_av_wee_paper = spherical_av_wee_paper * 0.25d0 / pi 
end


double precision  function spherical_av_n2_ab(r1,r12)
 implicit none
 include 'utils/constants.include.F'
 double precision, intent(in) :: r1(3),r12
 integer :: k
 double precision :: n2_hf_alpha_beta,r2(3)
 spherical_av_n2_ab = 0.d0
 do k = 1,n_points_integration_angular
  r2(1)= r1(1)+r12*angular_quadrature_points(k,1)
  r2(2)= r1(2)+r12*angular_quadrature_points(k,2)
  r2(3)= r1(3)+r12*angular_quadrature_points(k,3)
  call give_n2_alpha_beta_hf_at_r1_r2(r1,r2,n2_hf_alpha_beta) 
  spherical_av_n2_ab += n2_hf_alpha_beta*weights_angular_points(k)
 enddo
 spherical_av_n2_ab = spherical_av_n2_ab * 0.25d0 / pi 
end

double precision  function spherical_av_f_paper_ab(r1,r12)
 implicit none
 include 'utils/constants.include.F'
 double precision, intent(in) :: r1(3),r12
 integer :: k
 double precision :: f_hf_alpha_beta,r2(3)
 spherical_av_f_paper_ab = 0.d0
 do k = 1,n_points_integration_angular
  r2(1)= r1(1)+r12*angular_quadrature_points(k,1)
  r2(2)= r1(2)+r12*angular_quadrature_points(k,2)
  r2(3)= r1(3)+r12*angular_quadrature_points(k,3)
  call give_f_alpha_beta_hf_at_r1_r2(r1,r2,f_hf_alpha_beta) 
  spherical_av_f_paper_ab += f_hf_alpha_beta * weights_angular_points(k)
 enddo
 spherical_av_f_paper_ab = spherical_av_f_paper_ab * 0.25d0 / pi 
end

double precision  function spherical_av_n2(r1,r12)
 implicit none
 include 'utils/constants.include.F'
 double precision, intent(in) :: r1(3),r12
 integer :: k
 double precision :: n2_hf_alpha_alpha,r2(3)
 spherical_av_n2 = 0.d0
 do k = 1,n_points_integration_angular
  r2(1)= r1(1)+r12*angular_quadrature_points(k,1)
  r2(2)= r1(2)+r12*angular_quadrature_points(k,2)
  r2(3)= r1(3)+r12*angular_quadrature_points(k,3)
  call give_n2_alpha_alpha_hf_at_r1_r2(r1,r2,n2_hf_alpha_alpha) 
  spherical_av_n2 += n2_hf_alpha_alpha*weights_angular_points(k)
 enddo
 spherical_av_n2 = spherical_av_n2 * 0.25d0 / pi 
end

double precision  function spherical_av_f(r1,r12)
 implicit none
 include 'utils/constants.include.F'
 double precision, intent(in) :: r1(3),r12
 integer :: k
 double precision :: f_hf_alpha_alpha,r2(3)
 spherical_av_f= 0.d0
 do k = 1,n_points_integration_angular
  r2(1)= r1(1)+r12*angular_quadrature_points(k,1)
  r2(2)= r1(2)+r12*angular_quadrature_points(k,2)
  r2(3)= r1(3)+r12*angular_quadrature_points(k,3)
  call give_f_alpha_alpha_hf_at_r1_r2(r1,r2,f_hf_alpha_alpha) 
  spherical_av_f+= f_hf_alpha_alpha * weights_angular_points(k)
 enddo
 spherical_av_f= spherical_av_f * 0.25d0 / pi 
end

