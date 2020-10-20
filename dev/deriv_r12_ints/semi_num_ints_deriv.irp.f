BEGIN_PROVIDER [double precision, ao_two_e_eff_dr12_pot_array_new, (ao_num, ao_num, ao_num, ao_num)]
 implicit none
                                       !   1 2                                1 2 
! ao_two_e_eff_dr12_pot_array(k,i,l,j) = < k l | [erf( mu r12) - 1] d/d_r12 | i j > on the AO basis
 integer :: i,j,k,l,ipoint,m
 double precision :: weight1,thr,r(3)
 thr = 1.d-8
 ao_two_e_eff_dr12_pot_array_new = 0.d0
! double precision, allocatable :: a_mat(:,:,:,:)
! allocate(a_mat(n_points_final_grid,ao_num,ao_num,3))
! do m = 1, 3
!  do i = 1, ao_num
!   do k = 1, ao_num
!    do ipoint = 1, n_points_final_grid
!     r(1) = final_grid_points(1,ipoint)
!     r(2) = final_grid_points(2,ipoint)
!     r(3) = final_grid_points(3,ipoint)
!     weight1 = final_weight_at_r_vector(ipoint)
!     !   x*k grad i 
!     a_mat(ipoint,k,i,m) = aos_grad_in_r_array_transp_bis(ipoint,i,m) * aos_in_r_array_transp(ipoint,k) * r(m) * weight1
!    enddo
!   enddo
!  enddo
! enddo
 double precision, allocatable :: b_mat(:,:,:,:)
 allocate(b_mat(ao_num,ao_num,n_points_final_grid,3))
 do m = 1, 3
  do ipoint = 1, n_points_final_grid
   r(1) = final_grid_points(1,ipoint)
   r(2) = final_grid_points(2,ipoint)
   r(3) = final_grid_points(3,ipoint)
   weight1 = final_weight_at_r_vector(ipoint)
   do i = 1, ao_num
    do k = 1, ao_num
     b_mat(k,i,ipoint,m) = 0.5d0 * aos_in_r_array_transp(ipoint,k) * r(m) * weight1 * aos_grad_in_r_array(i,ipoint,m) 
    enddo
   enddo
  enddo
 enddo
  do m = 1, 3
   do ipoint = 1, n_points_final_grid
    do j = 1, ao_num ! 2
!     if(dabs(aos_grad_in_r_array(j,ipoint,m)).lt.thr)cycle
     do l = 1, ao_num ! 2
!      if(dabs(b_mat(l,j,ipoint,m)).lt.thr)cycle
      do i = 1, ao_num ! 1
       do k = 1, ao_num ! 1
        !                               1 1 2 2                [k*i](1)     * [l * x * grad_x j](2)
        ao_two_e_eff_dr12_pot_array_new(k,i,l,j) += v_ij_erf_rk(k,i,ipoint) * b_mat(l,j,ipoint,m) 
       enddo
      enddo
     enddo
    enddo
   enddo
  enddo

  do m = 1, 3
   do ipoint = 1, n_points_final_grid
    do j = 1, ao_num ! 2
     do l = 1, ao_num ! 2
      if(dabs(v_ij_erf_rk(l,j,ipoint)).lt.thr) cycle
      do i = 1, ao_num ! 1
       if(dabs(aos_grad_in_r_array(i,ipoint,m)*v_ij_erf_rk(l,j,ipoint)).lt.thr)cycle
       do k = 1, ao_num ! 1
        !                               1 1 2 2                [l*j](2)     * [k * x * grad_x i](1)
        ao_two_e_eff_dr12_pot_array_new(k,i,l,j) += v_ij_erf_rk(l,j,ipoint) * b_mat(k,i,ipoint,m) 
       enddo
      enddo
     enddo
    enddo
   enddo
  enddo

 do m = 1, 3
  do ipoint = 1, n_points_final_grid
   r(1) = final_grid_points(1,ipoint)
   r(2) = final_grid_points(2,ipoint)
   r(3) = final_grid_points(3,ipoint)
   weight1 = final_weight_at_r_vector(ipoint)
   do i = 1, ao_num
    do k = 1, ao_num
     b_mat(k,i,ipoint,m) = 0.5d0 * aos_in_r_array_transp(ipoint,k) * weight1 * aos_grad_in_r_array(i,ipoint,m) 
    enddo
   enddo
  enddo
 enddo

  do m = 1, 3
   do ipoint = 1, n_points_final_grid
    do j = 1, ao_num ! 2
     if(dabs(aos_grad_in_r_array(j,ipoint,m)).lt.thr)cycle
     do l = 1, ao_num ! 2
      if(dabs(b_mat(l,j,ipoint,m)).lt.thr)cycle
      do i = 1, ao_num ! 1
       do k = 1, ao_num ! 1
        !                               1 1 2 2                        [k*x*i](1)       * [l *  grad_x j](2)
        ao_two_e_eff_dr12_pot_array_new(k,i,l,j) -= x_v_ij_erf_rk_trans(k,i,ipoint,m) * b_mat(l,j,ipoint,m) 
       enddo
      enddo
     enddo
    enddo
   enddo
  enddo

  do m = 1, 3
   do ipoint = 1, n_points_final_grid
    do j = 1, ao_num ! 2
     do l = 1, ao_num ! 2
      if(dabs(v_ij_erf_rk(l,j,ipoint)).lt.thr) cycle
      do i = 1, ao_num ! 1
       if(dabs(aos_grad_in_r_array(i,ipoint,m)*v_ij_erf_rk(l,j,ipoint)).lt.thr)cycle
       do k = 1, ao_num ! 1
        !                               1 1 2 2                [l*x*j](2)     * [k * grad_x i](1)
        ao_two_e_eff_dr12_pot_array_new(k,i,l,j) -= x_v_ij_erf_rk_trans(l,j,ipoint,m) * b_mat(k,i,ipoint,m) 
       enddo
      enddo
     enddo
    enddo
   enddo
  enddo

END_PROVIDER 
