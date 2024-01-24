
program test
 implicit none
 call test_gauss_ints_aos

end


subroutine test_gauss_ints_aos
 implicit none 
 integer :: iao,jao,kao,lao
 double precision :: integral_new, integral_1, integral_2,accu_abs, accu_relat, integral_ref
 double precision :: get_ao_tc_sym_two_e_pot, ao_tc_sym_two_e_pot, ao_two_e_integral_erf 
 double precision :: delta
 double precision :: j1b_gauss_coul_modifdebug, integral_3

 accu_abs = 0.d0
 accu_relat = 0.d0
! provide 
 do iao = 1, ao_num ! r1
  do jao = 1, ao_num ! r2
   do kao = 1, ao_num ! r1
    do lao = 1, ao_num ! r2

     integral_new = get_ao_tc_sym_two_e_pot(iao,jao,kao,lao,ao_tc_sym_two_e_pot_map)

     integral_1   = ao_two_e_integral_erf(iao,kao,jao,lao)
     integral_2   = ao_tc_sym_two_e_pot(iao,kao,jao,lao)
     integral_3   = j1b_gauss_coul_modifdebug(iao,kao,jao,lao)
     integral_ref = integral_2 + integral_1 + integral_3

     delta = dabs(integral_ref - integral_new)
     accu_abs += delta
     if(delta.gt.1.d-7)then
      print*,iao,jao,kao,lao
      print*,integral_ref,integral_new,delta
     endif
    enddo
   enddo
  enddo
 enddo
 print*,'accu_abs = ',accu_abs/dble(ao_num)**4

end
