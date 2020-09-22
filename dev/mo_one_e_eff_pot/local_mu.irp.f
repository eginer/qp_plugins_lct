double precision function local_v_ne(r)
 implicit none
 double precision, intent(in) :: r(3)
 integer :: i,j,k
 double precision :: rho, mos_array(mo_num)
 call give_all_mos_at_r(r,mos_array)
 rho = mos_array(1) * mos_array(1)
 local_v_ne = 0.d0
 do i = 1, mo_num
  local_v_ne += mo_integrals_n_e(i,1) * mos_array(i) * mos_array(1) 
 enddo
 local_v_ne = local_v_ne / rho
end
