
subroutine test_gauss_ints_aos
 implicit none 
 double precision :: weight,r(3),weightj,r2(3),int_ao,int_gauss_num_2
 double precision :: ao_two_e_integral_eff_pot
 double precision :: int_r2,r_12,int_gauss_num,alpha_r12,coef
 integer :: ipoint,i,j,n_pt_in,jpoint
 double precision :: aos_array_r1(ao_num),aos_array_r2(ao_num)
 double precision :: accu_relat,accu_abs,err_relat,err_abs,exact
 double precision :: overlap_gauss_r12_ao
 include 'utils/constants.include.F'
 integer :: iao,jao,kao,lao
 accu_abs = 0.d0
 accu_relat = 0.d0
 do iao = 1, ao_num ! r1
  do jao = 1, ao_num ! r2
   do kao = 1, ao_num ! r1
    do lao = 1, ao_num ! r2
     int_ao = ao_two_e_integral_eff_pot(iao,kao,jao,lao)
     int_gauss_num = 0.d0
     int_gauss_num_2 = 0.d0
     do ipoint = 1, n_points_final_grid
      r(1) = final_grid_points(1,ipoint)
      r(2) = final_grid_points(2,ipoint)
      r(3) = final_grid_points(3,ipoint)
      call give_all_aos_at_r(r,aos_array_r1)
      weight = final_weight_at_r_vector(ipoint)
      int_r2 = 0.d0
      do i = 1,n_gauss_eff_pot
       alpha_r12 = expo_gauss_eff_pot(i)
       coef      = coef_gauss_eff_pot(i)
       exact = overlap_gauss_r12_ao(r,alpha_r12,jao,lao)
       int_r2   +=  coef * exact
      enddo
      int_gauss_num += weight * int_r2 * aos_array_r1(iao) * aos_array_r1(kao)
!      int_r2 = 0.d0
!      do jpoint = 1, n_points_final_grid
!       r2(1) = final_grid_points(1,jpoint)
!       r2(2) = final_grid_points(2,jpoint)
!       r2(3) = final_grid_points(3,jpoint)
!       call give_all_aos_at_r(r2,aos_array_r2)
!       weightj = final_weight_at_r_vector(jpoint)
!       r_12 = (r(1) - r2(1))**2 + (r(2) - r2(2))**2 + (r(3) - r2(3))**2 
!       do i = 1,n_gauss_eff_pot
!        alpha_r12 = expo_gauss_eff_pot(i)
!        if(alpha_r12 * r_12.gt.20.d0)cycle
!        coef      = coef_gauss_eff_pot(i)
!        int_r2   += aos_array_r2(jao) * aos_array_r2(lao) * dexp(-alpha_r12*r_12) * coef * weightj
!       enddo
!      enddo
!      int_gauss_num_2 += weight * int_r2 * aos_array_r1(iao) * aos_array_r1(kao)
     enddo
     err_abs = dabs(int_gauss_num - int_ao)
     if(dabs(int_gauss_num).gt.1.d-10)then
      err_relat = err_abs/dabs(int_gauss_num)
     else
      err_relat = 0.d0
     endif
     if(err_relat.gt.1.d-6)then
      print*,'AHAHAHAH'
      print*,'<ij|kl> = ',iao,jao,kao,lao
      print*,'int_ao        = ',int_ao
      print*,'int_gauss_num = ',int_gauss_num
      print*,'abs error     = ',err_abs
      print*,'err_relat     = ',err_relat
!      print*,'int_gauss_num2= ',int_gauss_num_2
      stop
     endif
     accu_abs += err_abs
     accu_relat += err_abs
    enddo
   enddo
  enddo
 enddo
 print*,'accu_abs   = ',accu_abs/dble(ao_num**4)
 print*,'accu_relat = ',accu_relat/dble(ao_num**4)

end


subroutine test_gauss_ints_map
 implicit none
 include 'utils/constants.include.F'
 integer :: iao,jao,kao,lao
 double precision :: int_ao,ao_two_e_integral_eff_pot,map_ao
 double precision :: get_ao_two_e_integral_eff_pot

 print*,'Beginning the AO map test'
 do iao = 1, ao_num ! r1
  do jao = 1, ao_num ! r2
   do kao = 1, ao_num ! r1
    do lao = 1, ao_num ! r2
     int_ao = ao_two_e_integral_eff_pot(iao,kao,jao,lao)
     map_ao = get_ao_two_e_integral_eff_pot(iao,jao,kao,lao,ao_integrals_eff_pot_map)
     if(dabs(int_ao - map_ao).gt.1.d-10)then
      print*,'ahahah'
      print*,iao,jao,kao,lao
      print*,int_ao,map_ao
     endif
    enddo
   enddo
  enddo
 enddo
 print*,'passed the AO map test'

end
