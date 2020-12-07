subroutine full_num_lapl_ao
 implicit none
 integer :: ipoint,jpoint,i,j,k,l,m,pp
 double precision :: r1(3), r2(3), weight1, weight2,r12 ,accu_tmp
 double precision :: alpha, coef,contrib,accu
 double precision, allocatable :: array_tmp(:,:,:,:)
! double precision, allocatable :: array_tmp_2(:,:,:,:)
 allocate(array_tmp(ao_num,ao_num,ao_num,ao_num))
! allocate(array_tmp_2(ao_num,ao_num,ao_num,ao_num))
 array_tmp = 0.d0
 array_tmp_2 = 0.d0
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
   do j = 1, ao_num
    do l = 1, ao_num
     do i = 1, ao_num
      do k = 1, ao_num
       do pp = 1, n_max_fit_ten_no_slat
        alpha = expo_fit_ten_no_slat_gauss(pp)
        coef = coef_fit_ten_no_slat_gauss(pp) 
        contrib =  weight1 * weight2 * aos_in_r_array_transp(ipoint,k) * aos_in_r_array_transp(ipoint,i)        &
                                     * aos_in_r_array_transp(jpoint,l) * aos_in_r_array_transp(jpoint,j)        & 
                                     * (-6.d0 * alpha + 4.d0 * alpha*alpha * r12)
        array_tmp(k,i,l,j) += contrib* dexp(-alpha * r12) * coef 
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
     accu_tmp = dabs(ao_ten_no_lapl_pot_old(k,i,l,j) - array_tmp(k,i,l,j))
     accu += accu_tmp
     print*,k,i,l,j
     print*,'ao_ten_no_lapl_pot, array_tmp, accu_tmp'
     print*, ao_ten_no_lapl_pot(k,i,l,j), array_tmp(k,i,l,j), accu_tmp
    enddo
   enddo
  enddo
 enddo
 print*,''
 print*,''
 print*,'accu/ao_num**4',accu/dble(ao_num)**4.d0
!
end

