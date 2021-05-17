BEGIN_PROVIDER [ double precision, gauss_ij_rk,  ( ao_num, ao_num,n_points_final_grid)]
 implicit none
 BEGIN_DOC
! mu_of_r_gauss(j,i,R) = int dr phi_i(r) phi_j(r) exp(-(\mu(R) |r - R|)^2)
 END_DOC
 integer :: i,j,ipoint
 double precision :: mu,r(3),overlap_gauss_r12_ao
 double precision :: int_mu, delta
 provide mu_erf final_grid_points constant_mu 
 double precision :: wall0, wall1
 call wall_time(wall0)
 if(constant_mu)then
 !$OMP PARALLEL                  &
 !$OMP DEFAULT (NONE)            &
 !$OMP PRIVATE (i,j,ipoint,mu,r,int_mu,delta) & 
 !$OMP SHARED (ao_num,n_points_final_grid,mu_erf,gauss_ij_rk,final_grid_points)
 !$OMP DO SCHEDULE (dynamic)
 do ipoint = 1, n_points_final_grid
   do i = 1, ao_num
    do j = i, ao_num
     mu = mu_erf 
     r(1) = final_grid_points(1,ipoint)
     r(2) = final_grid_points(2,ipoint)
     r(3) = final_grid_points(3,ipoint)
     delta = mu * mu
     int_mu = overlap_gauss_r12_ao(r,delta,i,j)
     gauss_ij_rk(j,i,ipoint)= int_mu 
   enddo
  enddo
 enddo
 !$OMP END DO
 !$OMP END PARALLEL

 else 
  provide mu_of_r_for_ints
 !$OMP PARALLEL                  &
 !$OMP DEFAULT (NONE)            &
 !$OMP PRIVATE (i,j,ipoint,mu,r,int_mu,delta) & 
 !$OMP SHARED (ao_num,n_points_final_grid,mu_of_r_for_ints,gauss_ij_rk,final_grid_points)
 !$OMP DO SCHEDULE (dynamic)
 do ipoint = 1, n_points_final_grid
   do i = 1, ao_num
    do j = i, ao_num
     mu = mu_of_r_for_ints(ipoint,1)
     delta = mu * mu
     r(1) = final_grid_points(1,ipoint)
     r(2) = final_grid_points(2,ipoint)
     r(3) = final_grid_points(3,ipoint)
     int_mu = overlap_gauss_r12_ao(r,delta,i,j)
     gauss_ij_rk(j,i,ipoint)= int_mu 
   enddo
  enddo
 enddo
 !$OMP END DO
 !$OMP END PARALLEL

 endif

 do ipoint = 1, n_points_final_grid
  do i = 1, ao_num
   do j = 1, i-1
     gauss_ij_rk(j,i,ipoint)= gauss_ij_rk(i,j,ipoint)
    enddo
   enddo
  enddo

 call wall_time(wall1)
 print*,'wall time for gauss_ij_rk  ',wall1 - wall0

END_PROVIDER 

BEGIN_PROVIDER [ double precision, gauss_2_ij_rk,  ( ao_num, ao_num,n_points_final_grid)]
 implicit none
 BEGIN_DOC
! mu_of_r_gauss(j,i,R) = int dr phi_i(r) phi_j(r) exp(-2(\mu(R) |r - R|)^2)
 END_DOC
 integer :: i,j,ipoint
 double precision :: mu,r(3),overlap_gauss_r12_ao
 double precision :: int_mu, delta
 provide mu_erf final_grid_points constant_mu 
 double precision :: wall0, wall1
 call wall_time(wall0)
 if(constant_mu)then
 !$OMP PARALLEL                  &
 !$OMP DEFAULT (NONE)            &
 !$OMP PRIVATE (i,j,ipoint,mu,r,int_mu,delta) & 
 !$OMP SHARED (ao_num,n_points_final_grid,mu_erf,gauss_2_ij_rk,final_grid_points)
 !$OMP DO SCHEDULE (dynamic)
 do ipoint = 1, n_points_final_grid
   do i = 1, ao_num
    do j = i, ao_num
     mu = mu_erf 
     r(1) = final_grid_points(1,ipoint)
     r(2) = final_grid_points(2,ipoint)
     r(3) = final_grid_points(3,ipoint)
     delta = 2.d0 * mu * mu
     int_mu = overlap_gauss_r12_ao(r,delta,i,j)
     gauss_2_ij_rk(j,i,ipoint)= int_mu 
   enddo
  enddo
 enddo
 !$OMP END DO
 !$OMP END PARALLEL

 else 
  provide mu_of_r_for_ints
 !$OMP PARALLEL                  &
 !$OMP DEFAULT (NONE)            &
 !$OMP PRIVATE (i,j,ipoint,mu,r,int_mu,delta) & 
 !$OMP SHARED (ao_num,n_points_final_grid,mu_of_r_for_ints,gauss_2_ij_rk,final_grid_points)
 !$OMP DO SCHEDULE (dynamic)
 do ipoint = 1, n_points_final_grid
   do i = 1, ao_num
    do j = i, ao_num
     mu = mu_of_r_for_ints(ipoint,1)
     delta = 2.d0 * mu * mu
     r(1) = final_grid_points(1,ipoint)
     r(2) = final_grid_points(2,ipoint)
     r(3) = final_grid_points(3,ipoint)
     int_mu = overlap_gauss_r12_ao(r,delta,i,j)
     gauss_2_ij_rk(j,i,ipoint)= int_mu 
   enddo
  enddo
 enddo
 !$OMP END DO
 !$OMP END PARALLEL

 endif

 do ipoint = 1, n_points_final_grid
  do i = 1, ao_num
   do j = 1, i-1
     gauss_2_ij_rk(j,i,ipoint)= gauss_2_ij_rk(i,j,ipoint)
    enddo
   enddo
  enddo

 call wall_time(wall1)
 print*,'wall time for gauss_2_ij_rk  ',wall1 - wall0

END_PROVIDER 

BEGIN_PROVIDER [ double precision, gauss_ij_xyz_rk_tmp_bis,  ( ao_num, ao_num,n_points_final_grid,3)]
 implicit none
 BEGIN_DOC
! mu_of_r_gauss(j,i,R,m) = int dr x/y/z phi_i(r) phi_j(r) exp(-(\mu(R) |r - R|)^2)
!
! with m == 1 ==> x, m == 2 ==> y, m == 3 ==> z
 END_DOC
 integer :: i,j,ipoint,m
 double precision :: mu,r(3),overlap_gauss_r12_ao
 double precision :: int_mu, delta
 double precision :: overlap_gauss_xyz_r12_ao_specific
 provide mu_erf final_grid_points constant_mu 
 double precision :: wall0, wall1
 call wall_time(wall0)
 if(constant_mu)then
 !$OMP PARALLEL                  &
 !$OMP DEFAULT (NONE)            &
 !$OMP PRIVATE (i,j,ipoint,mu,r,int_mu,delta,m) & 
 !$OMP SHARED (ao_num,n_points_final_grid,mu_erf,gauss_ij_xyz_rk_tmp_bis,final_grid_points)
 !$OMP DO SCHEDULE (dynamic)
 do m = 1, 3
  do ipoint = 1, n_points_final_grid
    do i = 1, ao_num
     do j = i, ao_num
      mu = mu_erf 
      r(1) = final_grid_points(1,ipoint)
      r(2) = final_grid_points(2,ipoint)
      r(3) = final_grid_points(3,ipoint)
      delta = mu * mu
      int_mu = overlap_gauss_xyz_r12_ao_specific(r,delta,i,j,m)
      gauss_ij_xyz_rk_tmp_bis(j,i,ipoint,m)= int_mu 
    enddo
   enddo
  enddo
 enddo
 !$OMP END DO
 !$OMP END PARALLEL

 else 
  provide mu_of_r_for_ints
 !$OMP PARALLEL                  &
 !$OMP DEFAULT (NONE)            &
 !$OMP PRIVATE (i,j,ipoint,mu,r,int_mu,delta) & 
 !$OMP SHARED (ao_num,n_points_final_grid,mu_of_r_for_ints,gauss_ij_xyz_rk_tmp_bis,final_grid_points)
 !$OMP DO SCHEDULE (dynamic)
 do m = 1, 3
  do ipoint = 1, n_points_final_grid
    do i = 1, ao_num
     do j = i, ao_num
      mu = mu_of_r_for_ints(ipoint,1)
      delta = mu * mu
      r(1) = final_grid_points(1,ipoint)
      r(2) = final_grid_points(2,ipoint)
      r(3) = final_grid_points(3,ipoint)
      int_mu = overlap_gauss_xyz_r12_ao_specific(r,delta,i,j,m)
      gauss_ij_xyz_rk_tmp_bis(j,i,ipoint,m)= int_mu 
    enddo
   enddo
  enddo
 enddo
 !$OMP END DO
 !$OMP END PARALLEL

 endif

 do m = 1, 3
  do ipoint = 1, n_points_final_grid
   do i = 1, ao_num
    do j = 1, i-1
      gauss_ij_xyz_rk_tmp_bis(j,i,ipoint,m)= gauss_ij_xyz_rk_tmp_bis(i,j,ipoint,m)
    enddo
   enddo
  enddo
 enddo

 call wall_time(wall1)
 print*,'wall time for gauss_ij_xyz_rk_tmp_bis  ',wall1 - wall0

END_PROVIDER 

BEGIN_PROVIDER [ double precision, gauss_ij_xyz_rk_tmp,  (3, ao_num, ao_num,n_points_final_grid)]
 implicit none
 BEGIN_DOC
! gauss_ij_xyz_rk_tmp(m,j,i,R) = int dr x/y/z phi_i(r) phi_j(r) exp(-(\mu(R) |r - R|)^2)
!
! with m == 1 ==> x, m == 2 ==> y, m == 3 ==> z
!
! only temporary --> will be reshaped in gauss_ij_xyz_rk with x/y/z index as the last index 
 END_DOC
 integer :: i,j,ipoint,m
 double precision :: mu,r(3),overlap_gauss_r12_ao
 double precision :: int_mu(3), delta
 provide mu_erf final_grid_points constant_mu 
 double precision :: wall0, wall1
 call wall_time(wall0)
 if(constant_mu)then
 !$OMP PARALLEL                  &
 !$OMP DEFAULT (NONE)            &
 !$OMP PRIVATE (i,j,ipoint,mu,r,int_mu,delta,m) & 
 !$OMP SHARED (ao_num,n_points_final_grid,mu_erf,gauss_ij_xyz_rk_tmp,final_grid_points)
 !$OMP DO SCHEDULE (dynamic)
 do ipoint = 1, n_points_final_grid
   do i = 1, ao_num
    do j = i, ao_num
     mu = mu_erf 
     r(1) = final_grid_points(1,ipoint)
     r(2) = final_grid_points(2,ipoint)
     r(3) = final_grid_points(3,ipoint)
     delta = mu * mu
     call overlap_gauss_xyz_r12_ao(r,delta,i,j,int_mu)
     do m = 1, 3
      gauss_ij_xyz_rk_tmp(m,j,i,ipoint)= int_mu(m) 
     enddo
   enddo
  enddo
 enddo
 !$OMP END DO
 !$OMP END PARALLEL

 else 
  provide mu_of_r_for_ints
 !$OMP PARALLEL                  &
 !$OMP DEFAULT (NONE)            &
 !$OMP PRIVATE (i,j,ipoint,mu,r,int_mu,delta,m) & 
 !$OMP SHARED (ao_num,n_points_final_grid,mu_of_r_for_ints,gauss_ij_xyz_rk_tmp,final_grid_points)
 !$OMP DO SCHEDULE (dynamic)
 do ipoint = 1, n_points_final_grid
   do i = 1, ao_num
    do j = i, ao_num
     mu = mu_of_r_for_ints(ipoint,1)
     delta = mu * mu
     r(1) = final_grid_points(1,ipoint)
     r(2) = final_grid_points(2,ipoint)
     r(3) = final_grid_points(3,ipoint)
     call overlap_gauss_xyz_r12_ao(r,delta,i,j,int_mu)
     do m = 1, 3
      gauss_ij_xyz_rk_tmp(m,j,i,ipoint)= int_mu(m) 
     enddo
   enddo
  enddo
 enddo
 !$OMP END DO
 !$OMP END PARALLEL

 endif

  do ipoint = 1, n_points_final_grid
   do i = 1, ao_num
    do j = 1, i-1
     do m = 1, 3
      gauss_ij_xyz_rk_tmp(m,j,i,ipoint)= gauss_ij_xyz_rk_tmp(m,i,j,ipoint)
     enddo
    enddo
   enddo
  enddo

 call wall_time(wall1)
 print*,'wall time for gauss_ij_xyz_rk_tmp  ',wall1 - wall0

END_PROVIDER 

BEGIN_PROVIDER [ double precision, gauss_ij_xyz_rk,  (ao_num, ao_num,n_points_final_grid,3)]
 implicit none
 BEGIN_DOC
! gauss_ij_xyz_rk(j,i,R,m) = int dr x/y/z phi_i(r) phi_j(r) exp(-(\mu(R) |r - R|)^2)
!
! with m == 1 ==> x, m == 2 ==> y, m == 3 ==> z
 END_DOC

 integer :: ipoint,i,j,m
 provide gauss_ij_xyz_rk_tmp

 double precision :: wall0, wall1
 call wall_time(wall0)
 do ipoint = 1, n_points_final_grid
  do i = 1, ao_num
   do j = 1, ao_num
    do m = 1, 3
     gauss_ij_xyz_rk(j,i,ipoint,m)= gauss_ij_xyz_rk_tmp(m,j,i,ipoint)
    enddo
   enddo
  enddo
 enddo

 call wall_time(wall1)
 print*,'wall time for gauss_ij_xyz_rk = ',wall1 - wall0


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
