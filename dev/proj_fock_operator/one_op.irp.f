double precision function kin_proj_op_in_r_on_phi_k(k,r)
 implicit none
 double precision, intent(in) :: r(3)
 integer, intent(in) :: k
 double precision :: mos_array(mo_num)
 integer :: j
 kin_proj_op_in_r_on_phi_k = 0.d0
 call give_all_mos_at_r(r,mos_array) 
 do j = 1, mo_num
  kin_proj_op_in_r_on_phi_k += mo_kinetic_integrals(j,k) * mos_array(j)
 enddo

end

double precision function v_ne_proj_op_in_r_on_phi_k(k,r)
 implicit none
 double precision, intent(in) :: r(3)
 integer, intent(in) :: k
 double precision :: mos_array(mo_num)
 integer :: j
 v_ne_proj_op_in_r_on_phi_k = 0.d0
 call give_all_mos_at_r(r,mos_array) 
 do j = 1, mo_num
  v_ne_proj_op_in_r_on_phi_k += mo_integrals_n_e(j,k) * mos_array(j)
 enddo
end

double precision function kin_op_in_r_on_phi_k(k,r)
 implicit none
 double precision, intent(in) :: r(3)
 integer, intent(in) :: k
 double precision :: aos_array(ao_num)
 double precision :: aos_grad_array(3,ao_num)
 double precision :: aos_lapl_array(3,ao_num)
 integer :: j
 kin_op_in_r_on_phi_k = 0.d0 
 call give_all_aos_and_grad_and_lapl_at_r(r,aos_array,aos_grad_array,aos_lapl_array)
 do j = 1, ao_num 
  kin_op_in_r_on_phi_k  +=  mo_coef(j,k) * ( aos_lapl_array(1,j) + aos_lapl_array(2,j) + aos_lapl_array(3,j) )
 enddo
 kin_op_in_r_on_phi_k *= -0.5d0 
end


double precision function v_ne_op_in_r_on_phi_k(k,r)
 implicit none
 double precision, intent(in) :: r(3)
 integer, intent(in) :: k
 integer :: j 
 double precision :: r_j
 double precision :: mos_array(mo_num)
 v_ne_op_in_r_on_phi_k = 0.d0
 call give_all_mos_at_r(r,mos_array) 
 do j = 1, nucl_num
  r_j = (nucl_coord_transp(1,j) - r(1))**2 + (nucl_coord_transp(2,j) - r(2))**2 + (nucl_coord_transp(3,j) - r(3))**2
  r_j = dsqrt(r_j)
  v_ne_op_in_r_on_phi_k -= nucl_charge(j) /r_j 
 enddo
  v_ne_op_in_r_on_phi_k *= mos_array(k)

end
