program pouet
 implicit none
! provide x_v_ij_erf_rk v_ij_erf_rk
! call test_ints_semi
! call test_semi_num_bis
! call test_new_ints

 call test_new_old_ints
! provide ao_two_e_eff_dr12_pot_array
end

subroutine test_ints_semi
 implicit none
 integer :: ipoint,i,j,m
 double precision :: r1(3), aos_grad_array_r1(3, ao_num), aos_array_r1(ao_num)
 double precision :: C_center(3), weight1,mu_in,r12,derf_mu_x,ints(3),integral,NAI_pol_mult_erf_ao
 double precision :: ao_mat(ao_num,ao_num),ao_xmat(3,ao_num,ao_num),accu1, accu2(3)
 mu_in = mu_erf 
 C_center = 0.d0
 C_center(1) = 0.25d0
 C_center(3) = 1.12d0
 C_center(2) = -1.d0
 ao_mat = 0.d0
 ao_xmat = 0.d0
 do ipoint = 1, n_points_final_grid
  r1(1) = final_grid_points(1,ipoint)
  r1(2) = final_grid_points(2,ipoint)
  r1(3) = final_grid_points(3,ipoint)
  call give_all_aos_and_grad_at_r(r1,aos_array_r1,aos_grad_array_r1)
  weight1 = final_weight_at_r_vector(ipoint)
  r12 = (r1(1) - C_center(1))**2.d0 + (r1(2) - C_center(2))**2.d0 + (r1(3) - C_center(3))**2.d0 
  r12 = dsqrt(r12)
  do i = 1, ao_num
   do j = 1, ao_num
    ao_mat(j,i)  += aos_array_r1(i) * aos_array_r1(j) * weight1 * derf_mu_x(mu_in,r12)
    do m = 1, 3
     ao_xmat(m,j,i) += r1(m) * aos_array_r1(i) * aos_array_r1(j) * weight1 * derf_mu_x(mu_in,r12)
    enddo
   enddo
  enddo
 enddo

 accu1 = 0.d0
 accu2 = 0.d0
 accu1relat = 0.d0
 accu2relat = 0.d0
 double precision :: accu1relat, accu2relat(3)
 do i = 1, ao_num
  do j = 1, ao_num
   integral =  NAI_pol_mult_erf_ao(i,j,mu_in,C_center)
   if(dabs(ao_mat(j,i)).gt.1.d-10)then
    accu1relat = dabs(integral - ao_mat(j,i))/dabs(ao_mat(j,i))
   endif
   if(dabs(integral - ao_mat(j,i)).gt.1.d-5)then
    print*,'i,j',i,j
    print*,integral,ao_mat(j,i),dabs(integral - ao_mat(j,i))
   endif
   accu1 += dabs(integral - ao_mat(j,i))
   call NAI_pol_x_mult_erf_ao(i,j,mu_in,C_center,ints)
   do m = 1, 3
    accu2(m) += dabs(ao_xmat(m,j,i) - ints(m))
    if(dabs(ao_xmat(m,j,i)).gt.1.d-10)then
     accu2relat(m) += dabs(ao_xmat(m,j,i) - ints(m))/dabs(ao_xmat(m,j,i))
    endif
   enddo
  enddo
 enddo
 print*,'accu1      = ',accu1/dble(ao_num * ao_num)
 print*,'accu2      = '
 print*,accu2(:)/dble(ao_num * ao_num)
 print*,''
 print*,'accu1relat = ',accu1relat/dble(ao_num * ao_num)
 print*,'accu2relat = '
 print*, accu2relat /dble(ao_num * ao_num)

end


subroutine test_semi_num_bis
 implicit none 
 double precision :: weight1,r1(3),weight2,r2(3),int_ao
 double precision :: ao_two_e_integral_schwartz_accel_gauss
 double precision :: int_r2(3),int_gauss_num,alpha_r12,coef,r12,int_r2_bis(3)
 integer :: ipoint,i,j,n_pt_in,jpoint,m
 double precision :: aos_array_r1(ao_num),aos_array_r2(ao_num),aos_grad_array_r1(3,ao_num),aos_grad_array_r2(3,ao_num)
 double precision :: accu_relat,accu_abs,err_relat,err_abs
 double precision :: accu_tmp(3),mu_in,derf_mu_x,accu_tmp_bis(3)
 double precision :: d_dr12_large(3),d_dr12(3)
 double precision :: d_dr12_large_2(3),d_dr12_2(3),total_int
 include 'utils/constants.include.F'
 integer :: iao,jao,kao,lao
 double precision :: test, ao_two_e_integral_schwartz_accel_erf,num_int
 mu_in = mu_erf
 accu_abs = 0.d0
 accu_relat = 0.d0
 num_int = 0.d0
 do jao = 1, ao_num ! r1
  do lao = 1, ao_num ! r2
   do iao = 1, ao_num ! r1
    do kao = 1, ao_num ! r2
    num_int += 1.d0
     print*,'<ij|kl> = ',jao,lao,iao,kao
     call ao_two_e_d_dr12_int(iao,jao,kao,lao,mu_in,d_dr12)
     int_ao = ao_two_e_eff_dr12_pot_array_new(jao,iao,lao,kao)
     print*,'int_ao = ',int_ao
     int_gauss_num = 0.d0
     accu_tmp = 0.d0
     accu_tmp_bis = 0.d0
     do ipoint = 1, n_points_final_grid
      r1(1) = final_grid_points(1,ipoint)
      r1(2) = final_grid_points(2,ipoint)
      r1(3) = final_grid_points(3,ipoint)
      call give_all_aos_and_grad_at_r(r1,aos_array_r1,aos_grad_array_r1)
      weight1 = final_weight_at_r_vector(ipoint)
      int_r2 = 0.d0
      int_r2_bis = 0.d0
      do jpoint = 1, n_points_final_grid
       r2(1) = final_grid_points(1,jpoint)
       r2(2) = final_grid_points(2,jpoint)
       r2(3) = final_grid_points(3,jpoint)
       weight2 = final_weight_at_r_vector(jpoint)
       call give_all_aos_and_grad_at_r(r2,aos_array_r2,aos_grad_array_r2)
       r12 = (r1(1) - r2(1))**2.d0 + (r1(2) - r2(2))**2.d0 + (r1(3) - r2(3))**2.d0 
       r12 = dsqrt(r12)
       do m = 1, 3
        int_r2_bis(m) += weight2 * 0.5d0 * derf_mu_x(mu_in,r12) & 
                             * (r1(m) - r2(m)) *  aos_array_r1(jao) * aos_array_r2(lao)  & 
                             * (aos_grad_array_r1(m,iao) * aos_array_r2(kao) - aos_grad_array_r2(m,kao) * aos_array_r1(iao))

!       !!!!!!!!! x1 dx1 i(r1)
!        int_r2(m) += weight2 * 0.5d0 * derf_mu_x(mu_in,r12)      & 
!                             * r1(m) * aos_grad_array_r1(m,iao)  & 
!                             *  aos_array_r2(kao)                & 
!                             *  aos_array_r1(jao) * aos_array_r2(lao)  
       !!!!!!!!! x2 dx2 k(r2)
!        int_r2(m) += weight2 * 0.5d0 * derf_mu_x(mu_in,r12)      & 
!                             * r2(m) * aos_grad_array_r2(m,kao)  & 
!                             *  aos_array_r1(iao)                & 
!                             *  aos_array_r1(jao) * aos_array_r2(lao)  
       !!!!!!!!! x1 i(r1) dx2 k(r1)
!        int_r2(m) -= weight2 * 0.5d0 * derf_mu_x(mu_in,r12)      & 
!                             * r1(m) * aos_array_r1(iao)         & 
!                             * aos_grad_array_r2(m,kao)          & 
!                             * aos_array_r1(jao) * aos_array_r2(lao)  
       !!!!!!!!! x2 k(r2) dx1 i(r1)
!        int_r2(m) -= weight2 * 0.5d0 * derf_mu_x(mu_in,r12)      & 
!                             * r2(m) * aos_array_r2(kao)  & 
!                             * aos_grad_array_r1(m,iao)                & 
!                             * aos_array_r1(jao) * aos_array_r2(lao)  
       enddo
      enddo
      do m = 1, 3
       accu_tmp(m) += weight1 * int_r2_bis(m)
!       accu_tmp_bis(m) += weight1 * int_r2(m)
      enddo
     enddo
     int_gauss_num = 0.d0
     do m = 1, 3
      int_gauss_num += accu_tmp(m)
     enddo
      err_abs = dabs(int_gauss_num - int_ao)
      if(dabs(int_gauss_num).gt.1.d-10)then
       err_relat = err_abs/dabs(int_gauss_num)
      else
       err_relat = 0.d0
      endif
      print*,'int_gauss_num = ',int_gauss_num
      print*,'int_ao        = ',int_ao
      print*,'abs error     = ',err_abs
      print*,'err_relat     = ',err_relat
      if(err_relat .gt. 1.d-1)then
       print*,'AHAHAHAAH'
       stop
      endif
      accu_abs += err_abs
      accu_relat += err_relat
    enddo
   enddo
  enddo
 enddo
 print*,''
 print*,''
 print*,''
 print*,'Summary'
 print*,''
 print*,''
 print*,''
 print*,'accu_abs   = ',accu_abs/num_int
 print*,'accu_relat = ',accu_relat/num_int

end

subroutine test_new_ints
 implicit none
 double precision :: accu,err
 integer :: i,j,k,l
 provide ao_two_e_eff_dr12_pot_array_new_3
 accu = 0.d0
 do j = 1, ao_num
  do l = 1, ao_num
   do i = 1, ao_num
    do k = 1, ao_num
     err = dabs(ao_two_e_eff_dr12_pot_array_new(k,i,l,j) - ao_two_e_eff_dr12_pot_array_new_3(k,i,l,j))
     if(err.gt.1.d-10)then
      print*,'k,i,l,j',k,i,l,j
      print*,err,ao_two_e_eff_dr12_pot_array_new(k,i,l,j),ao_two_e_eff_dr12_pot_array_new_3(k,i,l,j)
     endif
     accu += err
    enddo
   enddo
  enddo
 enddo
 print*,'accu = ',accu

end

subroutine test_new_old_ints
 implicit none
 integer :: i,j,k,l
 double precision :: accu,err
 accu = 0.d0
 do j = 1, mo_num
  do i = 1, mo_num
   do l = 1, mo_num
    do k = 1, mo_num
     err = dabs(mo_two_e_eff_dr12_pot_array_physicist(k,l,i,j) -  &
                mo_two_e_eff_dr12_pot_array(k,l,i,j) ) 
     if(err.gt.1.d-10)then
      print*,'k,l,i,j',k,l,i,j
      print*,err,mo_two_e_eff_dr12_pot_array_physicist(k,l,i,j),mo_two_e_eff_dr12_pot_array(k,l,i,j)
     endif
     accu += err
    enddo
   enddo
  enddo
 enddo
 print*,'accu = ',accu

end
