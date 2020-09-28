
subroutine test_gauss_ints_mos
 implicit none 
 double precision :: weight,r(3),weightj,r2(3),int_mo
 double precision :: mo_two_e_integral_eff_pot
 double precision :: int_r2,r12,int_gauss_num,alpha_r12,coef
 integer :: ipoint,i,j,n_pt_in,jpoint
 double precision :: mos_array_r1(mo_num),mos_array_r2(mo_num)
 double precision :: accu_relat,accu_abs,err_relat,err_abs
 double precision :: mos_ints_exact(mo_num,mo_num),mo_ints_tmp(mo_num,mo_num),mo_ints_num(mo_num,mo_num)
 include 'utils/constants.include.F'
 integer :: imo,jmo,kmo,lmo
 do imo = 1, mo_num
  do kmo = 1, mo_num
   mo_ints_num = 0.d0
   do ipoint = 1, n_points_final_grid
    r(1) = final_grid_points(1,ipoint)
    r(2) = final_grid_points(2,ipoint)
    r(3) = final_grid_points(3,ipoint)
    call give_all_mos_at_r(r,mos_array_r1)
    weight = final_weight_at_r_vector(ipoint)
    do i = 1,n_gauss_eff_pot                                                                                            
     alpha_r12 = expo_gauss_eff_pot(i)
     call overlap_gauss_r12_all_mo(r, alpha_r12, mo_ints_tmp)
     coef      = coef_gauss_eff_pot(i)
     do jmo = 1, mo_num
      do lmo = 1, mo_num
       mo_ints_num(lmo,jmo) +=  coef * mo_ints_tmp(lmo,jmo) * weight * mos_array_r1(imo) * mos_array_r1(kmo)
      enddo
     enddo
    enddo
   enddo
   do jmo = 1, mo_num
    do lmo = 1, mo_num
     int_mo = mo_two_e_integral_eff_pot(imo,jmo,kmo,lmo)
     int_gauss_num = mo_ints_num(lmo,jmo)
     err_abs = dabs(int_gauss_num - int_mo)
     if(dabs(int_gauss_num).gt.1.d-10)then
      err_relat = err_abs/dabs(int_gauss_num)
     else
      err_relat = 0.d0
     endif
     if(err_relat.gt.1.d-1)then
      print*,'AHAHHAHAH'
      print*,'imo,jmo,kmo,lmo'
      print*, imo,jmo,kmo,lmo 
      print*,int_gauss_num,int_mo
      print*,'AHAHHAHAH'
      print*,'int_gauss_num = ',int_gauss_num
      print*,'abs error     = ',err_abs
      print*,'err_relat     = ',err_relat
      stop
     endif
     accu_abs += err_abs
     accu_relat += err_relat
    enddo
   enddo
  enddo
 enddo
 print*,'accu_abs   = ',accu_abs/dble(mo_num**4)
 print*,'accu_relat = ',accu_relat/dble(mo_num**4)
end

subroutine test_gauss_ints_mos_exchange
 implicit none 
 double precision :: weight,r(3),weightj,r2(3),int_mo
 double precision :: mo_two_e_integral_eff_pot
 double precision :: int_r2,r12,int_gauss_num,alpha_r12,coef
 integer :: ipoint,i,j,n_pt_in,jpoint
 double precision :: mos_array_r1(mo_num),mos_array_r2(mo_num)
 double precision :: accu_relat,accu_abs,err_relat,err_abs
 double precision :: mos_ints_exact(mo_num,mo_num),mo_ints_tmp(mo_num,mo_num),mo_ints_num(mo_num,mo_num)
 include 'utils/constants.include.F'
 integer :: imo,jmo,kmo,lmo
 do imo = 1, mo_num
  do kmo = 1, mo_num
   print*,'imo,kmo',imo,kmo
   int_gauss_num = 0.d0
   do ipoint = 1, n_points_final_grid
    r(1) = final_grid_points(1,ipoint)
    r(2) = final_grid_points(2,ipoint)
    r(3) = final_grid_points(3,ipoint)
    call give_all_mos_at_r(r,mos_array_r1)
    weight = final_weight_at_r_vector(ipoint)
    int_r2 = 0.d0
    do jpoint = 1, n_points_final_grid
     r2(1) = final_grid_points(1,jpoint)
     r2(2) = final_grid_points(2,jpoint)
     r2(3) = final_grid_points(3,jpoint)
     call give_all_mos_at_r(r2,mos_array_r2)
     r12 = (r(1) - r2(1))**2.d0 + (r(2) - r2(2))**2.d0 + (r(3) - r2(3))**2.d0
     weightj = final_weight_at_r_vector(jpoint)
     do i = 1,n_gauss_eff_pot                                                                                            
      alpha_r12 = expo_gauss_eff_pot(i)
      coef      = coef_gauss_eff_pot(i)
      int_r2   += coef * weightj * mos_array_r2(imo) * mos_array_r2(kmo) * dexp(-alpha_r12 * r12)
     enddo
    enddo
    int_gauss_num += weight * int_r2 * mos_array_r1(imo) * mos_array_r1(kmo) 
   enddo
   int_mo = mo_two_e_integral_eff_pot(imo,imo,kmo,kmo)
   err_abs = dabs(int_gauss_num - int_mo)
   if(dabs(int_gauss_num).gt.1.d-10)then
    err_relat = err_abs/dabs(int_gauss_num)
   else
    err_relat = 0.d0
   endif
   print*,'int_gauss_num = ',int_gauss_num
   print*,'int_mo        = ',int_mo
   if(err_relat.gt.1.d-1)then
    print*,'AHAHHAHAH'
    print*,'imo,jmo,kmo,lmo'
    print*, imo,jmo,kmo,lmo 
    print*,int_gauss_num,int_mo
    print*,'AHAHHAHAH'
    print*,'int_gauss_num = ',int_gauss_num
    print*,'abs error     = ',err_abs
    print*,'err_relat     = ',err_relat
    stop
   endif
   accu_abs += err_abs
   accu_relat += err_relat
  enddo
 enddo
 print*,'accu_abs   = ',accu_abs/dble(mo_num**4)
 print*,'accu_relat = ',accu_relat/dble(mo_num**4)

end

subroutine test_fit_pouet
 implicit none
 double precision :: dx,xmax,x,eff_pot_gauss,eff_pot_fit_gauss
 integer :: nx,i
 nx = 1000
 xmax = 10.d0
 dx = xmax/dble(nx)

 x=0.d0
 do i = 1, nx
  write(33,*)x,eff_pot_gauss(x,mu_erf),eff_pot_fit_gauss(x)
  x += dx
 enddo

end

subroutine test_coulomb_exchange
 implicit none
 integer :: imo,jmo
 double precision :: int_mo,mo_two_e_integral_eff_pot
 do imo = 1, mo_num
  do jmo = 1, mo_num
     int_mo = mo_two_e_integral_eff_pot(imo,jmo,imo,jmo)
     print*,'*******'
     print*,'imo,jmo = ',imo,jmo
     print*,'coulomb = ',int_mo
     int_mo = mo_two_e_integral_eff_pot(imo,jmo,jmo,imo)
     print*,'exchang = ',int_mo
  enddo
 enddo


end

subroutine overlap_gauss_r12_all_mo(D_center, delta, mo_ints)
 implicit none
 double precision, intent(in) :: D_center(3), delta
 double precision, intent(out):: mo_ints(mo_num,mo_num)
 double precision :: ao_ints(ao_num,ao_num)
 call overlap_gauss_r12_all_ao(D_center,delta,ao_ints)
 call ao_to_mo(ao_ints,ao_num,mo_ints,mo_num)
end
