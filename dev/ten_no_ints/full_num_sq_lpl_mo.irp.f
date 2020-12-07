
subroutine full_num_square_lapl_mo
 implicit none
 integer :: ipoint,jpoint,i,j,k,l,m,pp,qq
 double precision :: r1(3), r2(3), weight1, weight2,r12 ,accu_tmp
 double precision :: alpha, coef,contrib,accu,beta,contrib2
 double precision, allocatable :: array_tmp(:,:,:,:)
 allocate(array_tmp(mo_num,mo_num,mo_num,mo_num))
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
   do j = 1, mo_num
    do l = 1, mo_num
     do i = 1, mo_num
      do k = 1, mo_num
       contrib =  weight1 * weight2 * r12 * mos_in_r_array_transp(ipoint,k) * mos_in_r_array_transp(ipoint,i)        &
                                          * mos_in_r_array_transp(jpoint,l) * mos_in_r_array_transp(jpoint,j)       
        ! square term 
        do qq =1, n_max_fit_ten_no_slat
         beta = expo_fit_ten_no_slat_gauss(qq)
         do pp = 1, n_max_fit_ten_no_slat
          alpha = expo_fit_ten_no_slat_gauss(pp)
          coef = coef_fit_ten_no_slat_gauss(pp) * coef_fit_ten_no_slat_gauss(qq)
          array_tmp(k,i,l,j) += contrib *  (4.d0 * alpha*beta)   * dexp(-(alpha+beta) * r12) * coef 
          enddo
         enddo

        ! laplacian term
        do pp = 1, n_max_fit_ten_no_slat
         alpha = expo_fit_ten_no_slat_gauss(pp)
         coef = coef_fit_ten_no_slat_gauss(pp) 
         contrib =  weight1 * weight2 * mos_in_r_array_transp(ipoint,k) * mos_in_r_array_transp(ipoint,i)        &
                                      * mos_in_r_array_transp(jpoint,l) * mos_in_r_array_transp(jpoint,j)        & 
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
 do j = 1, mo_num
  do l = 1, mo_num
   do i = 1, mo_num
    do k = 1, mo_num
     accu_tmp = dabs(mo_ten_no_eff_sq_lpl_pot_chemist(k,i,l,j) - array_tmp(k,i,l,j))
     accu += accu_tmp
     print*,k,i,l,j
     print*,'mo_ten_no_eff_sq_lpl_pot_chemist, array_tmp, accu_tmp'
     print*, mo_ten_no_eff_sq_lpl_pot_chemist(k,i,l,j), array_tmp(k,i,l,j), accu_tmp
    enddo
   enddo
  enddo
 enddo
 print*,''
 print*,''
 print*,'accu/mo_num**4',accu/dble(mo_num)**4.d0
!
end

