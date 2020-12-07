
BEGIN_PROVIDER [double precision, ao_ten_no_dr12_pot_old, (ao_num, ao_num, ao_num, ao_num)]
 implicit none
BEGIN_DOC
! Integral of the non hermite terms for the TenNo potential 
!
!                        1 1 2 2      1 2                                             1 2 
!
! ao_ten_no_dr12_pot_old(k,i,l,j) = < k l \grad_1 u_Tenno(r12) . (\grad_1 - \grad2) | i j > on the AO basis
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
     b_mat(k,i,ipoint,m) = aos_in_r_array(k,ipoint) * r(m) * weight1 * aos_grad_in_r_array(i,ipoint,m) 
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

 do m = 1, 3
  do ipoint = 1, n_points_final_grid
   r(1) = final_grid_points(1,ipoint)
   r(2) = final_grid_points(2,ipoint)
   r(3) = final_grid_points(3,ipoint)
   weight1 = final_weight_at_r_vector(ipoint)
   do i = 1, ao_num
    do k = 1, ao_num
     ! phi_k * d/dx phi_i
     b_mat(k,i,ipoint,m) = aos_in_r_array(k,ipoint) * weight1 * aos_grad_in_r_array(i,ipoint,m) 
    enddo
   enddo
  enddo
 enddo
!
 do pp = 1, n_max_fit_ten_no_slat
  alpha = expo_fit_ten_no_slat_gauss(pp)
  coef = coef_fit_ten_no_slat_gauss(pp)
  do m = 1, 3
   do ipoint = 1, n_points_final_grid
    do j = 1, ao_num ! 2
     do l = 1, ao_num ! 2
      do i = 1, ao_num ! 1
       do k = 1, ao_num ! 1
        !      1 1 2 2                                     [k*x*i](1)       * [l *  grad_x j](2)
        ac_mat(k,i,l,j) += + coef * 2.d0 * alpha * x_v_ij_gauss_rk(k,i,ipoint,m,pp) * b_mat(l,j,ipoint,m) 
       enddo
      enddo
     enddo
    enddo
   enddo
  enddo
 enddo

 do j = 1, ao_num
  do l = 1, ao_num
   do i = 1, ao_num
    do k = 1, ao_num
      ao_ten_no_dr12_pot_old(k,i,l,j) = ac_mat(k,i,l,j) + ac_mat(l,j,k,i)    
    enddo
   enddo
  enddo
 enddo
 double precision :: wall0, wall1
 call wall_time(wall1)
 print*,'time to provide ao_ten_no_dr12_pot_old = ',wall1 - wall0
END_PROVIDER 


BEGIN_PROVIDER [double precision, ao_ten_no_dr12_pot, (ao_num, ao_num, ao_num, ao_num)]
 implicit none
BEGIN_DOC
! Integral of the non hermite terms for the TenNo potential 
!
!                        1 1 2 2      1 2                                             1 2 
!
! ao_ten_no_dr12_pot_old(k,i,l,j) = < k l \grad_1 u_Tenno(r12) . (\grad_1 - \grad2) | i j > on the AO basis
END_DOC
 integer :: i,j,k,l,ipoint,m,pp
 double precision :: weight1,thr,r(3),alpha,coef,coeftmp
 thr = 1.d-8
 double precision, allocatable :: b_mat(:,:,:,:),ac_mat(:,:,:,:)
 provide v_ij_gauss_rk x_v_ij_gauss_rk
  call wall_time(wall0)
 allocate(b_mat(n_points_final_grid,ao_num,ao_num,3),ac_mat(ao_num, ao_num, ao_num, ao_num))
 !$OMP PARALLEL                  &
 !$OMP DEFAULT (NONE)            &
 !$OMP PRIVATE (i,k,m,ipoint,r,weight1) & 
 !$OMP SHARED (aos_in_r_array_transp,aos_grad_in_r_array_transp_bis,b_mat)& 
 !$OMP SHARED (ao_num,n_points_final_grid,final_grid_points,final_weight_at_r_vector)
 !$OMP DO SCHEDULE (static)
 do m = 1, 3
  do i = 1, ao_num
   do k = 1, ao_num
    do ipoint = 1, n_points_final_grid
     r(1) = final_grid_points(1,ipoint)
     r(2) = final_grid_points(2,ipoint)
     r(3) = final_grid_points(3,ipoint)
     weight1 = final_weight_at_r_vector(ipoint)
     b_mat(ipoint,k,i,m) = aos_in_r_array_transp(ipoint,k) * r(m) * weight1 * aos_grad_in_r_array_transp_bis(ipoint,i,m) 
    enddo
   enddo
  enddo
 enddo
 !$OMP END DO
 !$OMP END PARALLEL

 ! -2 * alpha * x * phi_k * d/dx phi_i \int phi_l phi_j exp(-alpha*r12*r12)
  ac_mat = 0.d0
  do pp = 1, n_max_fit_ten_no_slat
   alpha = expo_fit_ten_no_slat_gauss(pp)
   coef = coef_fit_ten_no_slat_gauss(pp)
   coeftmp = -coef * 2.d0 * alpha
   do m = 1, 3
    !           A   B^T  dim(A,1)       dim(B,2)       dim(A,2)        alpha * A                LDA 
    call dgemm("N","N",ao_num*ao_num,ao_num*ao_num,n_points_final_grid,coeftmp,v_ij_gauss_rk(1,1,1,pp),ao_num*ao_num & 
                      ,b_mat(1,1,1,m),n_points_final_grid,1.d0,ac_mat,ao_num*ao_num)
   enddo
  enddo

 !$OMP PARALLEL                  &
 !$OMP DEFAULT (NONE)            &
 !$OMP PRIVATE (i,k,m,ipoint,weight1) & 
 !$OMP SHARED (aos_in_r_array_transp,aos_grad_in_r_array_transp_bis,b_mat,ao_num,n_points_final_grid,final_weight_at_r_vector)
 !$OMP DO SCHEDULE (static)
 do m = 1, 3
  do i = 1, ao_num
   do k = 1, ao_num
    do ipoint = 1, n_points_final_grid
     weight1 = final_weight_at_r_vector(ipoint)
     b_mat(ipoint,k,i,m) =  aos_in_r_array_transp(ipoint,k) * weight1 * aos_grad_in_r_array_transp_bis(ipoint,i,m) 
    enddo
   enddo
  enddo
 enddo
 !$OMP END DO
 !$OMP END PARALLEL
 do pp = 1, n_max_fit_ten_no_slat
  do m = 1, 3
   alpha = expo_fit_ten_no_slat_gauss(pp)
   coef = coef_fit_ten_no_slat_gauss(pp)
   coeftmp = coef * 2.d0 * alpha
   !           A   B^T  dim(A,1)       dim(B,2)       dim(A,2)        alpha * A                LDA 
   call dgemm("N","N",ao_num*ao_num,ao_num*ao_num,n_points_final_grid,coeftmp,x_v_ij_gauss_rk(1,1,1,m,pp),ao_num*ao_num & 
                     ,b_mat(1,1,1,m),n_points_final_grid,1.d0,ac_mat,ao_num*ao_num)
  enddo
 enddo

 !$OMP PARALLEL                  &
 !$OMP DEFAULT (NONE)            &
 !$OMP PRIVATE (i,k,j,l) & 
 !$OMP SHARED (ac_mat,ao_ten_no_dr12_pot,ao_num)
 !$OMP DO SCHEDULE (static)
 do j = 1, ao_num
  do l = 1, ao_num
   do i = 1, ao_num
    do k = 1, ao_num
      ao_ten_no_dr12_pot(k,i,l,j) = ac_mat(k,i,l,j)  + ac_mat(l,j,k,i)    
    enddo
   enddo
  enddo
 enddo
 !$OMP END DO
 !$OMP END PARALLEL
  double precision :: wall1, wall0
  call wall_time(wall1)
  print*,'wall time dgemm ',wall1 - wall0
END_PROVIDER 
