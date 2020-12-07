subroutine full_num_deriv_ao
 implicit none
 integer :: ipoint,jpoint,i,j,k,l,m,pp
 double precision :: r1(3), r2(3), weight1, weight2,r12 ,accu_tmp
 double precision :: alpha, coef,contrib,accu
 double precision, allocatable :: array_tmp(:,:,:,:)
 allocate(array_tmp(ao_num,ao_num,ao_num,ao_num))
 array_tmp = 0.d0
 do ipoint = 1, n_points_final_grid
  r1(1) = final_grid_points(1,ipoint)
  r1(2) = final_grid_points(2,ipoint)
  r1(3) = final_grid_points(3,ipoint)
  weight1 = final_weight_at_r_vector(ipoint)
  do jpoint = 1, n_points_final_grid
   r2(1) = final_grid_points(1,jpoint)
   r2(2) = final_grid_points(2,jpoint)
   r2(3) = final_grid_points(3,jpoint)
   weight2 = final_weight_at_r_vector(jpoint)
   r12 = (r1(1) - r2(1))**2.d0 + (r1(2) - r2(2))**2.d0 + (r1(3) - r2(3))**2.d0
!   print*,ipoint,jpoint
   do m = 1, 3
    do j = 1, ao_num
     do l = 1, ao_num
      do i = 1, ao_num
       do k = 1, ao_num
        contrib =  weight1 * weight2 * aos_in_r_array_transp(ipoint,k) * aos_in_r_array_transp(jpoint,l)  & 
                           !   x1 (dx1 phi_i(1)) phi_j(2)
                           * ( r1(m) * aos_grad_in_r_array(i,ipoint,m) * aos_in_r_array_transp(jpoint,j)              &
                           !   x2 (dx2 phi_j(2)) phi_i(2)
                              +r2(m) * aos_grad_in_r_array(j,jpoint,m) * aos_in_r_array_transp(ipoint,i)              & 
                           !   x1 phi_i(1) (dx2 phi_j(2))
                              -r1(m) * aos_in_r_array_transp(ipoint,i) * aos_grad_in_r_array(j,jpoint,m)              & 
                           !   x2 phi_j(2) (dx1 phi_i(1))
                              -r2(m) * aos_in_r_array_transp(jpoint,j) * aos_grad_in_r_array(i,ipoint,m)              ) 
        do pp = 1, n_max_fit_ten_no_slat
         alpha = expo_fit_ten_no_slat_gauss(pp)
         coef = coef_fit_ten_no_slat_gauss(pp)
         array_tmp(k,i,l,j) += contrib *  (-2.d0 * alpha)   * dexp(-alpha * r12) * coef 
        enddo
       enddo
      enddo
     enddo
    enddo
   enddo
  enddo
 enddo

 accu = 0.d0
 do j = 1, ao_num
  do l = 1, ao_num
   do i = 1, ao_num
    do k = 1, ao_num
     accu_tmp = dabs(ao_ten_no_dr12_pot(k,i,l,j) - array_tmp(k,i,l,j))
     accu += accu_tmp
     print*,k,i,l,j
     print*,'ao_ten_no_dr12_pot, array_tmp, accu_tmp'
     print*, ao_ten_no_dr12_pot(k,i,l,j), array_tmp(k,i,l,j), accu_tmp
    enddo
   enddo
  enddo
 enddo
 print*,''
 print*,''
 print*,'accu/ao_num**4',accu/dble(ao_num)**4.d0

end

subroutine full_num_deriv_mo
 implicit none
 integer :: ipoint,jpoint,i,j,k,l,m,pp
 double precision :: r1(3), r2(3), weight1, weight2,r12 ,accu_tmp
 double precision :: alpha, coef,contrib,accu
 double precision, allocatable :: array_tmp(:,:,:,:)
 allocate(array_tmp(ao_num,ao_num,ao_num,ao_num))
 array_tmp = 0.d0
 do ipoint = 1, n_points_final_grid
  r1(1) = final_grid_points(1,ipoint)
  r1(2) = final_grid_points(2,ipoint)
  r1(3) = final_grid_points(3,ipoint)
  weight1 = final_weight_at_r_vector(ipoint)
  do jpoint = 1, n_points_final_grid
   r2(1) = final_grid_points(1,jpoint)
   r2(2) = final_grid_points(2,jpoint)
   r2(3) = final_grid_points(3,jpoint)
   weight2 = final_weight_at_r_vector(jpoint)
   r12 = (r1(1) - r2(1))**2.d0 + (r1(2) - r2(2))**2.d0 + (r1(3) - r2(3))**2.d0
   do m = 1, 3
    do j = 1, mo_num
     do l = 1, mo_num
      do i = 1, mo_num
       do k = 1, mo_num
        contrib =  weight1 * weight2 * mos_in_r_array_transp(ipoint,k) * mos_in_r_array_transp(jpoint,l)  & 
                           !   x1 (dx1 phi_i(1)) phi_j(2)
                           * ( r1(m) * mos_grad_in_r_array(i,ipoint,m) * mos_in_r_array_transp(jpoint,j)              &
                           !   x2 (dx2 phi_j(2)) phi_i(2)
                              +r2(m) * mos_grad_in_r_array(j,jpoint,m) * mos_in_r_array_transp(ipoint,i)              & 
                           !   x1 phi_i(1) (dx2 phi_j(2))
                              -r1(m) * mos_in_r_array_transp(ipoint,i) * mos_grad_in_r_array(j,jpoint,m)              & 
                           !   x2 phi_j(2) (dx1 phi_i(1))
                              -r2(m) * mos_in_r_array_transp(jpoint,j) * mos_grad_in_r_array(i,ipoint,m)              ) 
                                                
        do pp = 1, n_max_fit_ten_no_slat
         alpha = expo_fit_ten_no_slat_gauss(pp)
         coef = coef_fit_ten_no_slat_gauss(pp)
         array_tmp(k,i,l,j) += contrib *  (-2.d0 * alpha)   * dexp(-alpha * r12) * coef 
        enddo
       enddo
      enddo
     enddo
    enddo
   enddo
  enddo
 enddo

 accu = 0.d0
 do j = 1, mo_num
  do l = 1, mo_num
   do i = 1, mo_num
    do k = 1, mo_num
     accu_tmp = dabs(mo_ten_no_dr12_pot_chemist(k,i,l,j) - array_tmp(k,i,l,j))
     accu += accu_tmp
     print*,k,i,l,j
     print*,'mo_ten_no_dr12_pot_chemist, array_tmp, accu_tmp'
     print*, mo_ten_no_dr12_pot_chemist(k,i,l,j), array_tmp(k,i,l,j), accu_tmp
    enddo
   enddo
  enddo
 enddo
 print*,''
 print*,''
 print*,'accu/mo_num**4',accu/dble(mo_num)**4.d0

end
