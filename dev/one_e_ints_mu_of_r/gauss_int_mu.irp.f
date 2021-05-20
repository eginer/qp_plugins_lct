BEGIN_PROVIDER [ double precision, gauss_ij_rk,  (n_points_final_grid ,ao_num, ao_num)]
 implicit none
 BEGIN_DOC
! mu_of_r_gauss(R,j,l) = int dr phi_l(r) phi_j(r) exp(-(\mu(R) |r - R|)^2)
 END_DOC
 integer :: l,j,ipoint
 double precision :: mu,r(3),overlap_gauss_r12_ao
 double precision :: int_mu, delta,thr
 provide mu_erf final_grid_points 
 double precision :: wall0, wall1
 thr = 1.d-12
 call wall_time(wall0)
  provide mu_of_r_for_ints
 !$OMP PARALLEL                  &
 !$OMP DEFAULT (NONE)            &
 !$OMP PRIVATE (l,j,ipoint,mu,r,int_mu,delta) & 
 !$OMP SHARED (ao_num,n_points_final_grid,mu_of_r_for_ints,gauss_ij_rk,final_grid_points,aos_in_r_array_transp,thr)
 !$OMP DO SCHEDULE (dynamic)
 do l = 1, ao_num
  do j = l, ao_num
   do ipoint = 1, n_points_final_grid
     gauss_ij_rk(ipoint,j,l) = 0.d0
     if(dabs( aos_in_r_array_transp(ipoint,j) * aos_in_r_array_transp(ipoint,l)).lt.thr)cycle
     mu = mu_of_r_for_ints(ipoint,1)
     delta = mu * mu
     r(1) = final_grid_points(1,ipoint)
     r(2) = final_grid_points(2,ipoint)
     r(3) = final_grid_points(3,ipoint)
     int_mu = overlap_gauss_r12_ao(r,delta,l,j)
     gauss_ij_rk(ipoint,j,l)= int_mu 
   enddo
  enddo
 enddo
 !$OMP END DO
 !$OMP END PARALLEL


 do l = 1, ao_num
  do j = 1, l-1
   do ipoint = 1, n_points_final_grid
     gauss_ij_rk(ipoint,j,l)= gauss_ij_rk(ipoint,l,j)
   enddo
  enddo
 enddo

 call wall_time(wall1)
 print*,'wall time for gauss_ij_rk  ',wall1 - wall0

END_PROVIDER 

BEGIN_PROVIDER [ double precision, gauss_2_ij_rk,  (n_points_final_grid ,ao_num, ao_num)]
 implicit none
 BEGIN_DOC
! mu_of_r_gauss(R,j,i) = int dr phi_i(r) phi_j(r) exp(-2(\mu(R) |r - R|)^2)
 END_DOC
 integer :: l,j,ipoint
 double precision :: mu,r(3),overlap_gauss_r12_ao
 double precision :: int_mu, delta,thr
 thr = 1.d-12
 provide mu_erf final_grid_points 
 double precision :: wall0, wall1
 call wall_time(wall0)
  provide mu_of_r_for_ints
 !$OMP PARALLEL                  &
 !$OMP DEFAULT (NONE)            &
 !$OMP PRIVATE (l,j,ipoint,mu,r,int_mu,delta) & 
 !$OMP SHARED (ao_num,n_points_final_grid,mu_of_r_for_ints,gauss_2_ij_rk,final_grid_points,thr,aos_in_r_array_transp)
 !$OMP DO SCHEDULE (dynamic)
 do l = 1, ao_num
  do j = l, ao_num
   do ipoint = 1, n_points_final_grid
     gauss_2_ij_rk(ipoint,j,l) = 0.d0
     if(dabs( aos_in_r_array_transp(ipoint,j) * aos_in_r_array_transp(ipoint,l)).lt.thr)cycle
     mu = mu_of_r_for_ints(ipoint,1)
     delta = 2.d0 * mu * mu
     r(1) = final_grid_points(1,ipoint)
     r(2) = final_grid_points(2,ipoint)
     r(3) = final_grid_points(3,ipoint)
     int_mu = overlap_gauss_r12_ao(r,delta,l,j)
     gauss_2_ij_rk(ipoint,j,l)= int_mu 
   enddo
  enddo
 enddo
 !$OMP END DO
 !$OMP END PARALLEL


 do l = 1, ao_num
  do j = 1, l-1
   do ipoint = 1, n_points_final_grid
     gauss_2_ij_rk(ipoint,j,l)= gauss_2_ij_rk(ipoint,l,j)
    enddo
   enddo
  enddo

 call wall_time(wall1)
 print*,'wall time for gauss_2_ij_rk  ',wall1 - wall0

END_PROVIDER 

BEGIN_PROVIDER [ double precision, gauss_ij_xyz_rk_tmp,  (3,n_points_final_grid, ao_num, ao_num)]
 implicit none
 BEGIN_DOC
! gauss_ij_xyz_rk_tmp(m,j,i,R) = int dr x/y/z phi_i(r) phi_j(r) exp(-(\mu(R) |r - R|)^2)
!
! with m == 1 ==> x, m == 2 ==> y, m == 3 ==> z
!
! only temporary --> will be reshaped in gauss_ij_xyz_rk with x/y/z index as the last index 
 END_DOC
 integer :: l,j,ipoint,m
 double precision :: mu,r(3),overlap_gauss_r12_ao
 double precision :: int_mu(3), delta
 provide mu_erf final_grid_points 
 double precision :: wall0, wall1,thr
 thr = 1.d-12
 call wall_time(wall0)
  provide mu_of_r_for_ints
 !$OMP PARALLEL                  &
 !$OMP DEFAULT (NONE)            &
 !$OMP PRIVATE (l,j,ipoint,mu,r,int_mu,delta,m) & 
 !$OMP SHARED (ao_num,n_points_final_grid,mu_of_r_for_ints,gauss_ij_xyz_rk_tmp,final_grid_points,thr,aos_in_r_array_transp)
 !$OMP DO SCHEDULE (dynamic)
 do l = 1, ao_num
  do j = l, ao_num
   do ipoint = 1, n_points_final_grid
     gauss_ij_xyz_rk_tmp(:,ipoint,j,l) = 0.d0
     if(dabs( aos_in_r_array_transp(ipoint,j) * aos_in_r_array_transp(ipoint,l)).lt.thr)cycle
     mu = mu_of_r_for_ints(ipoint,1)
     delta = mu * mu
     r(1) = final_grid_points(1,ipoint)
     r(2) = final_grid_points(2,ipoint)
     r(3) = final_grid_points(3,ipoint)
     call overlap_gauss_xyz_r12_ao(r,delta,l,j,int_mu)
     do m = 1, 3
      gauss_ij_xyz_rk_tmp(m,ipoint,j,l)= int_mu(m) 
     enddo
   enddo
  enddo
 enddo
 !$OMP END DO
 !$OMP END PARALLEL


  do l = 1, ao_num
   do j = 1, l-1
    do ipoint = 1, n_points_final_grid
     do m = 1, 3
      gauss_ij_xyz_rk_tmp(m,ipoint,j,l)= gauss_ij_xyz_rk_tmp(m,ipoint,l,j)
     enddo
    enddo
   enddo
  enddo

 call wall_time(wall1)
 print*,'wall time for gauss_ij_xyz_rk_tmp  ',wall1 - wall0

END_PROVIDER 

BEGIN_PROVIDER [ double precision, gauss_ij_xyz_rk,  (n_points_final_grid,ao_num, ao_num,3)]
 implicit none
 BEGIN_DOC
! gauss_ij_xyz_rk(j,i,R,m) = int dr x/y/z phi_i(r) phi_j(r) exp(-(\mu(R) |r - R|)^2)
!
! with m == 1 ==> x, m == 2 ==> y, m == 3 ==> z
 END_DOC

 integer :: ipoint,l,j,m
 provide gauss_ij_xyz_rk_tmp

 double precision :: wall0, wall1
 call wall_time(wall0)
 do l = 1, ao_num
  do j = 1, ao_num
   do ipoint = 1, n_points_final_grid
    do m = 1, 3
     gauss_ij_xyz_rk(ipoint,j,l,m)= gauss_ij_xyz_rk_tmp(m,ipoint,j,l)
    enddo
   enddo
  enddo
 enddo

 call wall_time(wall1)
 print*,'wall time for gauss_ij_xyz_rk = ',wall1 - wall0
 FREE gauss_ij_xyz_rk_tmp


END_PROVIDER

subroutine test_gauss_ij_rk
 implicit none
 integer :: ipoint,i,j,m,jpoint
 double precision :: r1(3)
 double precision :: C_center(3), weight1,mu,r12,integral,delta, num_int
 double precision, allocatable :: ao_mat(:,:,:), aos_array_r1(:),ao_mat_xyz(:,:,:,:)
 allocate(ao_mat(ao_num,ao_num,n_points_final_grid), aos_array_r1(ao_num))
 allocate(ao_mat_xyz(3,ao_num,ao_num,n_points_final_grid))
 do jpoint = 1, n_points_final_grid
  mu = mu_of_r_for_ints(jpoint,1)
  delta = mu * mu
  C_center(:) = final_grid_points(:,jpoint)
  ao_mat(:,:,jpoint) = 0.d0
  do ipoint = 1, n_points_final_grid
   r1(1) = final_grid_points(1,ipoint)
   r1(2) = final_grid_points(2,ipoint)
   r1(3) = final_grid_points(3,ipoint)
   call give_all_aos_at_r(r1,aos_array_r1)
   weight1 = final_weight_at_r_vector(ipoint)
   r12 = (r1(1) - C_center(1))**2.d0 + (r1(2) - C_center(2))**2.d0 + (r1(3) - C_center(3))**2.d0 
   do i = 1, ao_num
    do j = 1, ao_num
     ao_mat(j,i,jpoint)  += aos_array_r1(i) * aos_array_r1(j) * weight1 * dexp(-delta*r12)
     do m = 1, 3
     ao_mat_xyz(m,j,i,jpoint)  += aos_array_r1(i) * aos_array_r1(j) * weight1 * dexp(-delta*r12) * r1(m)
     enddo
    enddo
   enddo
  enddo
 enddo

 double precision :: accu1relat,accu1
 do jpoint = 1, n_points_final_grid
  accu1 = 0.d0
  accu1relat = 0.d0
  print*,'jpoint = ',jpoint
  do i = 1, ao_num
   do j = 1, ao_num
    integral =  gauss_ij_rk(j,i,jpoint)
    num_int  =  ao_mat(j,i,jpoint)
    if(dabs(num_int).gt.1.d-10)then
     accu1relat = dabs(integral - num_int )/dabs(num_int)
    endif
    if(dabs(integral - num_int).gt.1.d-5)then
     print*,'i,j,jpoint',i,j,jpoint
     print*,'prov, num, difference'
     print*,integral,num_int,dabs(integral - num_int)
    endif
    accu1 += dabs(integral - num_int)
   enddo
  enddo
  print*,'accu1      = ',accu1/dble(ao_num * ao_num)
  print*,''
  print*,'accu1relat = ',accu1relat/dble(ao_num * ao_num)
 enddo

 print*,'m = ',m
 do jpoint = 1, n_points_final_grid
  accu1 = 0.d0
  accu1relat = 0.d0
  print*,'jpoint = ',jpoint
  do i = 1, ao_num
   do j = 1, ao_num
    do m = 1, 3
     integral =  gauss_ij_xyz_rk_tmp(m,j,i,jpoint)
     num_int  = ao_mat_xyz(m,j,i,jpoint) 
     if(dabs(num_int).gt.1.d-10)then
      accu1relat = dabs(integral - num_int)/dabs(num_int)
     endif
     if(dabs(integral - num_int).gt.1.d-5)then
      print*,'i,j,jpoint',i,j,jpoint
      print*,'prov, num, difference'
      print*,integral,num_int,dabs(integral - num_int)
     endif
     accu1 += dabs(integral - num_int)
    enddo
   enddo
  enddo
  print*,'accu1      = ',accu1/dble(ao_num * ao_num)
  print*,''
  print*,'accu1relat = ',accu1relat/dble(ao_num * ao_num)
 enddo

end

subroutine test_gauss_2_ij
 implicit none
 integer :: ipoint,i,j,m,jpoint
 double precision :: r2(3)
 double precision :: r1(3), weight1,mu,r12,integral,delta, num_int,ao_prod_r1,weight2
 double precision, allocatable :: ao_mat(:,:,:), aos_array_r2(:), ao_mat_2(:,:)
 allocate(ao_mat(ao_num,ao_num,n_points_final_grid), aos_array_r2(ao_num),ao_mat_2(ao_num, ao_num))
 ao_mat_2 = 0.d0
 do ipoint = 1, n_points_final_grid
  mu = mu_of_r_for_ints(ipoint,1)
  delta = 2.d0 * mu * mu
  r1(:) = final_grid_points(:,ipoint)
  ao_mat(:,:,ipoint) = 0.d0
  weight1 = final_weight_at_r_vector(ipoint)
  do jpoint = 1, n_points_final_grid
   r2(:) = final_grid_points(:,jpoint)
   weight2 = final_weight_at_r_vector(jpoint)
   r12 = (r1(1) - r2(1))**2.d0 + (r1(2) - r2(2))**2.d0 + (r1(3) - r2(3))**2.d0 
   do i = 1, ao_num
    do j = 1, ao_num
     ao_prod_r1 = aos_in_r_array(i,ipoint) * aos_in_r_array(j,ipoint) * weight1
     ao_mat_2(j,i) += ao_prod_r1 * aos_in_r_array(i,jpoint) * aos_in_r_array(j,jpoint) * weight2 * dexp(-delta*r12) 
     ao_mat(j,i,ipoint)  += aos_in_r_array(i,jpoint) * aos_in_r_array(j,jpoint) * weight2 * dexp(-delta*r12) 
    enddo
   enddo
  enddo
 enddo
 do i = 1, ao_num
  do j = 1, ao_num
   write(33,*)i,j,ao_mat_2(j,i)
  enddo
 enddo

 double precision :: accu1relat,accu1
 do ipoint = 1, n_points_final_grid
  accu1 = 0.d0
  accu1relat = 0.d0
  print*,'ipoint = ',ipoint
  do i = 1, ao_num
   do j = 1, ao_num
    integral =  gauss_2_ij_rk(j,i,ipoint)
    num_int  =  ao_mat(j,i,ipoint)
    if(dabs(num_int).gt.1.d-10)then
     accu1relat = dabs(integral - num_int )/dabs(num_int)
    endif
    if(dabs(integral - num_int).gt.1.d-5)then
     print*,'i,j,ipoint',i,j,ipoint
     print*,'prov, num, difference'
     print*,integral,num_int,dabs(integral - num_int)
    endif
    accu1 += dabs(integral - num_int)
   enddo
  enddo
  print*,'accu1      = ',accu1/dble(ao_num * ao_num)
  print*,''
  print*,'accu1relat = ',accu1relat/dble(ao_num * ao_num)
 enddo


end
