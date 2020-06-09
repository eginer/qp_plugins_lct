
subroutine test_gauss_ints_aos
 implicit none 
 double precision :: weight,r(3),weightj,r2(3),int_ao
 double precision :: ao_two_e_integral_schwartz_accel_gauss
 double precision :: int_r2,r_12,int_gauss_num,alpha_r12,coef
 integer :: ipoint,i,j,n_pt_in,jpoint
 double precision :: aos_array_r1(ao_num),aos_array_r2(ao_num)
 double precision :: accu_relat,accu_abs,err_relat,err_abs
 include 'utils/constants.include.F'
 integer :: iao,jao,kao,lao
 do iao = 1, ao_num ! r1
  do jao = 1, ao_num ! r2
   do kao = 1, ao_num ! r1
    do lao = 1, ao_num ! r2
 

     print*,'<ij|kl> = ',iao,jao,kao,lao
     int_ao = ao_two_e_integral_schwartz_accel_gauss(iao,kao,jao,lao)
     print*,'int_ao        = ',int_ao
     int_gauss_num = 0.d0
     do ipoint = 1, n_points_final_grid
      r(1) = final_grid_points(1,ipoint)
      r(2) = final_grid_points(2,ipoint)
      r(3) = final_grid_points(3,ipoint)
      call give_all_aos_at_r(r,aos_array_r1)
      weight = final_weight_at_r_vector(ipoint)
      int_r2 = 0.d0
      do jpoint = 1, n_points_final_grid
       r2(1) = final_grid_points(1,jpoint)
       r2(2) = final_grid_points(2,jpoint)
       r2(3) = final_grid_points(3,jpoint)
       call give_all_aos_at_r(r2,aos_array_r2)
       weightj = final_weight_at_r_vector(jpoint)
       r_12 = (r(1) - r2(1))**2 + (r(2) - r2(2))**2 + (r(3) - r2(3))**2 
       do i = 1,n_gauss_eff_pot
        alpha_r12 = expo_gauss_eff_pot(i)
        if(alpha_r12 * r_12.gt.20.d0)cycle
        coef      = coef_gauss_eff_pot(i)
        int_r2   += aos_array_r2(jao) * aos_array_r2(lao) * dexp(-alpha_r12*r_12) * coef * weightj
       enddo
      enddo
      int_gauss_num += weight * int_r2 * aos_array_r1(iao) * aos_array_r1(kao)
     enddo
     err_abs = dabs(int_gauss_num - int_ao)
     if(int_gauss_num.gt.1.d-10)then
      err_relat = err_abs/dabs(int_gauss_num)
     else
      err_relat = 0.d0
     endif
     print*,'int_gauss_num = ',int_gauss_num
     print*,'abs error     = ',err_abs
     print*,'err_relat     = ',err_relat
     accu_abs += err_abs
     accu_relat += accu_relat
    enddo
   enddo
  enddo
 enddo
 print*,'accu_abs   = ',accu_abs/dble(ao_num**4)
 print*,'accu_relat = ',accu_relat/dble(ao_num**4)

end


subroutine test_extra_basis
 implicit none
 integer :: ipoint,i_ao,m
 double precision :: r(3),accu(3,ao_num),weight,aos_array(ao_num)
 double precision :: xyz_phi(3),aos_grad_array(3,ao_num),grad_xyz_phi(3),xyz_grad_phi(3)
 accu = 0.d0
 do ipoint = 1, n_points_final_grid
  r(1) = final_grid_points(1,ipoint)
  r(2) = final_grid_points(2,ipoint)
  r(3) = final_grid_points(3,ipoint)
  call give_all_aos_and_grad_at_r(r,aos_array,aos_grad_array)
  weight = final_weight_at_r_vector(ipoint)
  do i_ao = 1, ao_num
   call xyz_phi_ao(r,i_ao,xyz_phi)
   call xyz_grad_phi_ao(r,i_ao,xyz_grad_phi)
   do m = 1, 3
!    if(dabs(aos_array(i_ao) * r(m) - xyz_phi(m)).gt.1.d-10)then
!    if(dabs(aos_grad_array(m,i_ao) - grad_xyz_phi(m)).gt.1.d-10)then
    if(dabs(r(m) * aos_grad_array(m,i_ao) - xyz_grad_phi(m)).gt.1.d-10)then
!     print*,i_ao,m,r(m)
     print*,i_ao,m,ao_power(i_ao,m)
     print*,r
     print*,aos_grad_array(m,i_ao)*r(m),xyz_grad_phi(m)
     pause
    endif
    accu(m,i_ao) += dabs(aos_grad_array(m,i_ao) * r(m) - xyz_grad_phi(m)) * weight
   enddo
  enddo
 enddo
 print*,''
 print*,'errors '
 print*,''
 do i_ao = 1, ao_num
  write(*,'(100(F16.10,X))')accu(:,i_ao)
 enddo



end
