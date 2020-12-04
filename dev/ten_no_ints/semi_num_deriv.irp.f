
BEGIN_PROVIDER [double precision, ao_two_e_eff_dr12_pot_array_new_bis, (ao_num, ao_num, ao_num, ao_num)]
 implicit none
BEGIN_DOC
!                             1 1 2 2      1 2                                1 2 
!
! ao_two_e_eff_dr12_pot_array(k,i,l,j) = < k l | [erf( mu r12) - 1] d/d_r12 | i j > on the AO basis
END_DOC
 integer :: i,j,k,l,ipoint,m,pp
 double precision :: weight1,thr,r(3)
 thr = 1.d-8
 double precision, allocatable :: b_mat(:,:,:,:),ac_mat(:,:,:,:)
 double precision :: alpha,coef
 provide v_ij_gauss_rk x_v_ij_gauss_rk
 call wall_time(wall0)
 allocate(b_mat(ao_num,ao_num,n_points_final_grid,3),ac_mat(ao_num, ao_num, ao_num, ao_num))
 do m = 1, 3
  do ipoint = 1, n_points_final_grid
   r(1) = final_grid_points(1,ipoint)
   r(2) = final_grid_points(2,ipoint)
   r(3) = final_grid_points(3,ipoint)
   weight1 = final_weight_at_r_vector(ipoint)
   do i = 1, ao_num
    do k = 1, ao_num
     ! x * phi_k * d/dx phi_i
     b_mat(k,i,ipoint,m) = aos_in_r_array_transp(ipoint,k) * r(m) * weight1 * aos_grad_in_r_array(i,ipoint,m) 
    enddo
   enddo
  enddo
 enddo

 ac_mat = 0.d0
 do pp = 1, n_max_fit_ten_no_slat
  alpha = expo_fit_ten_no_slat_gauss(pp)
  coef = coef_fit_ten_no_slat_gauss(pp)
  do m = 1, 3
   do ipoint = 1, n_points_final_grid
    do j = 1, ao_num ! 2
     do l = 1, ao_num ! 2
      do i = 1, ao_num ! 1
       do k = 1, ao_num ! 1
        !      1 1 2 2                          [k*i](1)     * [l * x * grad_x j](2)
        ac_mat(k,i,l,j) += -coef * 2.d0 * alpha * v_ij_gauss_rk(k,i,ipoint,pp) * b_mat(l,j,ipoint,m) 
       enddo
      enddo
     enddo
    enddo
   enddo
  enddo
 enddo
!
! do m = 1, 3
!  do ipoint = 1, n_points_final_grid
!   r(1) = final_grid_points(1,ipoint)
!   r(2) = final_grid_points(2,ipoint)
!   r(3) = final_grid_points(3,ipoint)
!   weight1 = final_weight_at_r_vector(ipoint)
!   do i = 1, ao_num
!    do k = 1, ao_num
!     ! phi_k * d/dx phi_i
!     b_mat(k,i,ipoint,m) = aos_in_r_array_transp(ipoint,k) * weight1 * aos_grad_in_r_array(i,ipoint,m) 
!    enddo
!   enddo
!  enddo
! enddo

! do pp = 1, n_max_fit_ten_no_slat
!  alpha = expo_fit_ten_no_slat_gauss(pp)
!  coef = coef_fit_ten_no_slat_gauss(pp)
!  do m = 1, 3
!   do ipoint = 1, n_points_final_grid
!    do j = 1, ao_num ! 2
!     do l = 1, ao_num ! 2
!      do i = 1, ao_num ! 1
!       do k = 1, ao_num ! 1
!        !      1 1 2 2                                     [k*x*i](1)       * [l *  grad_x j](2)
!        ac_mat(k,i,l,j) += + coef * 2.d0 * alpha * x_v_ij_gauss_rk(k,i,ipoint,m,pp) * b_mat(l,j,ipoint,m) 
!       enddo
!      enddo
!     enddo
!    enddo
!   enddo
!  enddo
! enddo

 do j = 1, ao_num
  do l = 1, ao_num
   do i = 1, ao_num
    do k = 1, ao_num
      ao_two_e_eff_dr12_pot_array_new_bis(k,i,l,j) = ac_mat(k,i,l,j) !+ ac_mat(l,j,k,i)    
    enddo
   enddo
  enddo
 enddo
 double precision :: wall0, wall1
 call wall_time(wall1)
 print*,'time to provide ao_two_e_eff_dr12_pot_array_new_bis = ',wall1 - wall0
END_PROVIDER 

