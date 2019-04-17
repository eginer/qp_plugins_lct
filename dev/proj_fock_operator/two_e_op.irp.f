BEGIN_PROVIDER [ double precision, mo_F_two_e_operator, (mo_num,mo_num)]
 implicit none
 integer :: i,j,k
 double precision :: mo_two_e_integral
 mo_F_two_e_operator = 0.d0
 do i = 1, mo_num
  do j = 1, mo_num
   do k = 1, elec_alpha_num
    mo_F_two_e_operator(j,i) += 2.d0 * mo_two_e_integral(i,k,j,k) - mo_two_e_integral(i,k,k,j)
   enddo
  enddo
 enddo
 mo_F_two_e_operator = mo_F_two_e_operator * 0.5d0
END_PROVIDER 

double precision function coulomb_exch_proj_op_in_r_on_phi_k(k,r)
 implicit none
 double precision, intent(in) :: r(3)
 integer, intent(in) :: k
 double precision :: mos_array(mo_num)
 integer :: j
 coulomb_exch_proj_op_in_r_on_phi_k = 0.d0
 call give_all_mos_at_r(r,mos_array) 
 do j = 1, mo_num
  coulomb_exch_proj_op_in_r_on_phi_k += mo_F_two_e_operator(j,k) * mos_array(j)
 enddo

end
