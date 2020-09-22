 BEGIN_PROVIDER [double precision, ao_erf_mu_r_inv_r, (ao_num, ao_num)]
&BEGIN_PROVIDER [double precision, ao_gauss_mu_pot_one_e, (ao_num, ao_num)]
&BEGIN_PROVIDER [double precision, ao_deriv_mu_pot_one_e, (ao_num, ao_num)]
&BEGIN_PROVIDER [double precision, ao_grad_one_e_j_squared , (ao_num, ao_num)]
 implicit none
 include 'utils/constants.include.F'
 BEGIN_DOC
! ao_erf_mu_r_inv_r(i,j) = \sum_A (-Z_A) \int dr phi_i(r) erf(mu * r_ia) / r_ia phi_j(r) 
! 
! ao_gauss_mu_pot_one_e(i,j) = \sum_A (-Z_A) \int dr phi_i(r) mu/sqrt(pi)exp(-(mu*r_ia)^2) phi_j(r) 
!
! ao_deriv_mu_pot_one_e(i,j) = \sum_A Z_A \int phi_i(r) (1 - erf(mu * r_ia))/r_ia vec(r_ia) . grad phi_j(r) 
!
! ao_grad_one_e_j_squared(i,j) = \sum_A Z_A 
 END_DOC
 double precision :: mu_in , r1(3) , aos_array_r1(ao_num) , aos_grad_array_r1(3,ao_num) , weight1
 double precision :: derf_mu_x, r_ia(3,nucl_num), dist_r_ia(nucl_num),erf_mu_ria(nucl_num),dist_vec(3,nucl_num)
 double precision :: poly_inv_r(3, nucl_num), poly_one_m_erf_ria_inv_ria(3,nucl_num)
 integer :: ipoint,i,j,k,l
 mu_in = mu_one_e_j
 ao_erf_mu_r_inv_r = 0.d0
 ao_gauss_mu_pot_one_e = 0.d0
 ao_deriv_mu_pot_one_e = 0.d0
 ao_grad_one_e_j_squared = 0.d0
 do ipoint = 1, n_points_final_grid
  r1(1) = final_grid_points(1,ipoint)
  r1(2) = final_grid_points(2,ipoint)
  r1(3) = final_grid_points(3,ipoint)
  call give_all_aos_and_grad_at_r(r1,aos_array_r1,aos_grad_array_r1)
  weight1 = final_weight_at_r_vector(ipoint)
  do i = 1, nucl_num
   r_ia(1,i) = r1(1) - nucl_coord_transp(1,i) 
   r_ia(2,i) = r1(2) - nucl_coord_transp(2,i) 
   r_ia(3,i) = r1(3) - nucl_coord_transp(3,i) 

   dist_r_ia(i)  = r_ia(1,i)*r_ia(1,i) + r_ia(2,i)*r_ia(2,i) + r_ia(3,i)*r_ia(3,i) 
   dist_r_ia(i)  = dsqrt(dist_r_ia(i))

   dist_vec(1,i) = dsqrt(r_ia(2,i)*r_ia(2,i) + r_ia(3,i)*r_ia(3,i))
   dist_vec(2,i) = dsqrt(r_ia(1,i)*r_ia(1,i) + r_ia(3,i)*r_ia(3,i))
   dist_vec(3,i) = dsqrt(r_ia(1,i)*r_ia(1,i) + r_ia(2,i)*r_ia(2,i))
   ! poly_inv_r(1) = x / sqrt(x^2+y^2+z^2), poly_inv_r(2) = y / sqrt(x^2+y^2+z^2), poly_inv_r(3) = z / sqrt(x^2+y^2+z^2)
   call inv_r_times_poly(r_ia(1,i), dist_r_ia(i), dist_vec(1,i), poly_inv_r(1,i))
   ! erf_mu_ria(i) = erf(mu * r_ia)/r_ia 
   erf_mu_ria(i) = derf_mu_x(mu_in,dist_r_ia(i))
   ! poly_one_m_erf_ria_inv_ria(1,i) = (x_i - x_A) * (1 - erf(mu * r_ia))/r_ia 
   poly_one_m_erf_ria_inv_ria(1,i) = (poly_inv_r(1,i) - erf_mu_ria(i) * r_ia(1,i))
   poly_one_m_erf_ria_inv_ria(2,i) = (poly_inv_r(2,i) - erf_mu_ria(i) * r_ia(2,i))
   poly_one_m_erf_ria_inv_ria(3,i) = (poly_inv_r(3,i) - erf_mu_ria(i) * r_ia(3,i))
   do l = 1, ao_num
    do k = 1, ao_num
     !  phi_k (r) * (-Z_A) * erf(mu * r_ia)/r_ia * phi_l(r) 
     ao_erf_mu_r_inv_r(k,l)     += - nucl_charge(i) * aos_array_r1(k) * aos_array_r1(l) * weight1 * erf_mu_ria(i) 
     ! effective gaussian potential integral 
     ! phi_k (r) * (-Z_a * mu/sqrt(pi) ) * exp(-(mu*r_ia)^2) * phi_l(r) 
     ao_gauss_mu_pot_one_e(k,l) += - nucl_charge(i) * mu_in * inv_sq_pi * aos_array_r1(k) * aos_array_r1(l) * weight1 * dexp(-mu_in*mu_in*dist_r_ia(i) * dist_r_ia(i))
     ! phi_k(r) * (Z_A) * (1 - erf(mu r_ia))/r_ia vec(r_ia) . grad phi_l 
     ao_deriv_mu_pot_one_e(k,l) += nucl_charge(i) *  weight1 * & 
    (poly_one_m_erf_ria_inv_ria(1,i) * aos_array_r1(k) * aos_grad_array_r1(1,l) & 
    +poly_one_m_erf_ria_inv_ria(2,i) * aos_array_r1(k) * aos_grad_array_r1(2,l) & 
    +poly_one_m_erf_ria_inv_ria(3,i) * aos_array_r1(k) * aos_grad_array_r1(3,l) ) 
    enddo
   enddo
  enddo
  do l = 1, ao_num
   do k = 1, ao_num
    do i = 1, nucl_num
     do j = 1, nucl_num
      ao_grad_one_e_j_squared(k,l) += -0.5d0 * nucl_charge(i) * nucl_charge(j) * aos_array_r1(k) * aos_array_r1(l) * weight1 & 
                                  *(poly_one_m_erf_ria_inv_ria(1,i) * poly_one_m_erf_ria_inv_ria(1,j) & 
                                   +poly_one_m_erf_ria_inv_ria(2,i) * poly_one_m_erf_ria_inv_ria(2,j) & 
                                   +poly_one_m_erf_ria_inv_ria(3,i) * poly_one_m_erf_ria_inv_ria(3,j) )

     enddo
    enddo
   enddo
  enddo

 enddo
END_PROVIDER 

 BEGIN_PROVIDER [double precision, mo_erf_mu_r_inv_r, (mo_num, mo_num)]
&BEGIN_PROVIDER [double precision, mo_gauss_mu_pot_one_e, (mo_num, mo_num)]
&BEGIN_PROVIDER [double precision, mo_deriv_mu_pot_one_e, (mo_num, mo_num)]
&BEGIN_PROVIDER [double precision, mo_grad_one_e_j_squared , (mo_num, mo_num)]
 implicit none
    call ao_to_mo(ao_erf_mu_r_inv_r,size(ao_erf_mu_r_inv_r,1),mo_erf_mu_r_inv_r,size(mo_erf_mu_r_inv_r,1))
    call ao_to_mo(ao_gauss_mu_pot_one_e,size(ao_gauss_mu_pot_one_e,1),mo_gauss_mu_pot_one_e,size(mo_gauss_mu_pot_one_e,1))
    call ao_to_mo(ao_deriv_mu_pot_one_e,size(ao_deriv_mu_pot_one_e,1),mo_deriv_mu_pot_one_e,size(mo_deriv_mu_pot_one_e,1))
    call ao_to_mo(ao_grad_one_e_j_squared,size(ao_grad_one_e_j_squared,1),mo_grad_one_e_j_squared,size(mo_grad_one_e_j_squared,1))

END_PROVIDER 

 BEGIN_PROVIDER [double precision, mo_eff_pot_one_e_j, (mo_num, mo_num)]
 implicit none
 integer :: i,j
 do i = 1, mo_num
  do j = 1, mo_num
!   mo_eff_pot_one_e_j(j,i) = mo_erf_mu_r_inv_r(j,i) + mo_gauss_mu_pot_one_e(j,i) + mo_deriv_mu_pot_one_e(j,i) + mo_grad_one_e_j_squared(j,i)
   mo_eff_pot_one_e_j(j,i) = mo_erf_integrals_n_e(j,i) + mo_gauss_mu_pot_one_e(j,i) + mo_deriv_mu_pot_one_e(j,i) + mo_grad_one_e_j_squared(j,i)
  enddo
 enddo

 END_PROVIDER 

BEGIN_PROVIDER [double precision, h_tilde_one_j, (mo_num, mo_num)]
 implicit none
 integer :: i,j
 do i = 1, mo_num
  do j = 1, mo_num
   h_tilde_one_j(j,i) = mo_kinetic_integrals(j,i) + mo_eff_pot_one_e_j(j,i)
  enddo
 enddo
END_PROVIDER 


subroutine inv_r_times_poly(r, dist_r, dist_vec, poly)
 implicit none
 BEGIN_DOC
! returns 
!
! poly(1) = x / sqrt(x^2+y^2+z^2), poly(2) = y / sqrt(x^2+y^2+z^2), poly(3) = z / sqrt(x^2+y^2+z^2)
!
! with the arguments  
!
! r(1)  = x, r(2) = y, r(3) = z, dist_r = sqrt(x^2+y^2+z^2)
!
! dist_vec(1) = sqrt(y^2+z^2), dist_vec(2) = sqrt(x^2+z^2), dist_vec(3) = sqrt(x^2+y^2)
 END_DOC
 double precision, intent(in) :: r(3), dist_r, dist_vec(3)
 double precision, intent(out):: poly(3)
 double precision :: inv_dist
 integer :: i
 if (dist_r.gt. 1.d-8)then
  inv_dist = 1.d0/dist_r
  do i = 1, 3
   poly(i) = r(i) * inv_dist 
  enddo
 else
  do i = 1, 3
   if(dabs(r(i)).lt.dist_vec(i))then
    inv_dist = 1.d0/dist_r
    poly(i) = r(i) * inv_dist 
   else !if(dabs(r(i)))then
    poly(i) = 1.d0 
!    poly(i) = 0.d0 
   endif
  enddo
 endif
end
