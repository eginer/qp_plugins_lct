
program test
 implicit none
 call test_gauss_ints_mos

end


subroutine test_gauss_ints_mos
 implicit none 
 integer :: imo,jmo,kmo,lmo
 double precision :: integral_new, integral_1, integral_2,accu_abs, accu_relat, integral_ref
 double precision :: get_mo_tc_sym_two_e_pot, mo_tc_sym_two_e_pot, mo_two_e_integral_erf 
 double precision :: delta
 accu_abs = 0.d0
 accu_relat = 0.d0
 provide 
 do imo = 1, mo_num ! r1
  do jmo = 1, mo_num ! r2
   do kmo = 1, mo_num ! r1
    do lmo = 1, mo_num ! r2
     integral_new = get_mo_tc_sym_two_e_pot(imo,jmo,kmo,lmo,mo_tc_sym_two_e_pot_map)
     integral_1   = mo_two_e_integral_erf(imo,jmo,jmo,lmo)
     integral_2   = mo_tc_sym_two_e_pot(imo,kmo,jmo,lmo)
     integral_ref = integral_2 + integral_1
     delta = dabs(integral_ref - integral_new)
     accu_abs += delta
     if(delta.gt.1.d-7)then
      print*,imo,jmo,kmo,lmo
      print*,integral_ref,integral_new,delta
     endif
    enddo
   enddo
  enddo
 enddo
 print*,'accu_abs = ',accu_abs/dble(mo_num)**4

end
