
BEGIN_PROVIDER [double precision, erf_mu_r12_inv_r12_rk,( ao_num, ao_num,n_points_final_grid)]
 implicit none
 BEGIN_DOC
! erf_mu_r12_inv_r12_rk(i,j,r) = \int dr' phi_i(r') phi_j(r') erf(mu(r) |r - r'|)/|r-r'|
 END_DOC
 double precision :: mu,r(3),int_mu,delta,wall0,wall1,phi_j_erf_mu_r_phi
 integer :: i,j,ipoint
  provide mu_of_r_for_ints
 call wall_time(wall0)
  !$OMP PARALLEL                  &
  !$OMP DEFAULT (NONE)            &
  !$OMP PRIVATE (i,j,ipoint,mu,r,int_mu,delta) & 
  !$OMP SHARED (ao_num,n_points_final_grid,mu_of_r_for_ints,erf_mu_r12_inv_r12_rk,final_grid_points)
  !$OMP DO SCHEDULE (dynamic)
  do ipoint = 1, n_points_final_grid
   do i = 1, ao_num
    do j = i, ao_num
     mu = mu_of_r_for_ints(ipoint,1)
     r(1) = final_grid_points(1,ipoint)
     r(2) = final_grid_points(2,ipoint)
     r(3) = final_grid_points(3,ipoint)
     int_mu = phi_j_erf_mu_r_phi(i,j,mu, r)
     erf_mu_r12_inv_r12_rk(j,i,ipoint) = int_mu 
    enddo
   enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

 do ipoint = 1, n_points_final_grid
  do i = 1, ao_num
   do j = 1, i-1
    erf_mu_r12_inv_r12_rk(j,i,ipoint)= erf_mu_r12_inv_r12_rk(i,j,ipoint)
   enddo
  enddo
 enddo

 call wall_time(wall1)
 print*,'wall time for erf_mu_r12_inv_r12_rk  ',wall1 - wall0

END_PROVIDER 


BEGIN_PROVIDER [double precision, coulomb_rk,( ao_num, ao_num,n_points_final_grid)]
 implicit none
 BEGIN_DOC
! coulomb_rk(i,j,r) = \int dr' phi_i(r') phi_j(r') 1/|r-r'|
 END_DOC
 double precision :: r(3),int_mu,delta,wall0,wall1,mu,phi_j_erf_mu_r_phi
 integer :: i,j,ipoint
 call wall_time(wall0)
 !$OMP PARALLEL                  &
 !$OMP DEFAULT (NONE)            &
 !$OMP PRIVATE (i,j,ipoint,r,int_mu,delta,mu) & 
 !$OMP SHARED (ao_num,n_points_final_grid,coulomb_rk,final_grid_points)
 !$OMP DO SCHEDULE (dynamic)
 do ipoint = 1, n_points_final_grid
  do i = 1, ao_num
   do j = i, ao_num
    r(1) = final_grid_points(1,ipoint)
    r(2) = final_grid_points(2,ipoint)
    r(3) = final_grid_points(3,ipoint)
    mu = 1.d+10
    int_mu = phi_j_erf_mu_r_phi(i,j,mu, r)
    coulomb_rk(j,i,ipoint) = int_mu
   enddo
  enddo
 enddo
 !$OMP END DO
 !$OMP END PARALLEL


 do ipoint = 1, n_points_final_grid
  do i = 1, ao_num
   do j = 1, i-1
    coulomb_rk(j,i,ipoint)= coulomb_rk(i,j,ipoint)
   enddo
  enddo
 enddo

 call wall_time(wall1)
 print*,'wall time for coulomb_rk  ',wall1 - wall0

END_PROVIDER 


BEGIN_PROVIDER [ double precision, gauss_erfc_mu_r12_inv_r12_rk_tmp,  (4,n_points_final_grid, ao_num, ao_num)]
 implicit none
 BEGIN_DOC
! gauss_erfc_mu_r12_inv_r12_rk_tmp(m,j,i,R) = int dr x/y/z phi_i(r) phi_j(r) exp(-(\mu(R) |r - R|)^2) * (1 - erf(mu(R) |r-R|))/|r-R|
!
! with m == 1 ==> x, m == 2 ==> y, m == 3 ==> z, m == 4 x^0/y^0/z^0
!
! only temporary --> will be reshaped in gauss_ij_xyz_rk with x/y/z index as the last index 
 END_DOC
 integer :: l,j,ipoint,m
 double precision :: mu,r(3),overlap_gauss_r12_ao
 double precision :: int_mu(3), delta, thr
 provide mu_erf final_grid_points 
 double precision :: wall0, wall1
 thr = 0.d0
 call wall_time(wall0)
  provide mu_of_r_for_ints
  print*,'Providing gauss_erfc_mu_r12_inv_r12_rk_tmp ...'
 !$OMP PARALLEL                  &
 !$OMP DEFAULT (NONE)            &
 !$OMP PRIVATE (l,j,ipoint,mu,r,int_mu,delta,m) & 
 !$OMP SHARED (ao_num,n_points_final_grid,mu_of_r_for_ints,gauss_erfc_mu_r12_inv_r12_rk_tmp,final_grid_points,thr,aos_in_r_array_transp)
 !$OMP DO SCHEDULE (dynamic)
 do l = 1, ao_num
  do j = l, ao_num
   do ipoint = 1, n_points_final_grid
    gauss_erfc_mu_r12_inv_r12_rk_tmp(:,ipoint,j,l) = 0.d0
    if(dabs( aos_in_r_array_transp(ipoint,j) * aos_in_r_array_transp(ipoint,l)).lt.thr)cycle
    mu = mu_of_r_for_ints(ipoint,1)
    delta = mu * mu
    r(1) = final_grid_points(1,ipoint)
    r(2) = final_grid_points(2,ipoint)
    r(3) = final_grid_points(3,ipoint)
    call erfc_mu_gauss_xyz_ij_ao(l,j,mu, r, delta,int_mu)
    do m = 1, 4
     gauss_erfc_mu_r12_inv_r12_rk_tmp(m,ipoint,j,l)= int_mu(m) 
    enddo
   enddo
  enddo
 enddo
 !$OMP END DO
 !$OMP END PARALLEL


  do l = 1, ao_num
   do j = 1, l-1
    do ipoint = 1, n_points_final_grid
     do m = 1, 4
      gauss_erfc_mu_r12_inv_r12_rk_tmp(m,ipoint,j,l)= gauss_erfc_mu_r12_inv_r12_rk_tmp(m,ipoint,l,j)
     enddo
    enddo
   enddo
  enddo
 call wall_time(wall1)
 print*,'wall time for gauss_erfc_mu_r12_inv_r12_rk_tmp ',wall1 - wall0


END_PROVIDER 

BEGIN_PROVIDER [ double precision, gauss_erfc_mu_r12_inv_r12_rk,  (n_points_final_grid,ao_num, ao_num,4)]
 implicit none
 BEGIN_DOC
! gauss_erfc_mu_r12_inv_r12_rk(j,i,R,m) = int dr x/y/z phi_i(r) phi_j(r) exp(-(\mu(R) |r - R|)^2) * (1 - erf(mu(R) |r-R|))/|r-R|
!
! with m == 1 ==> x, m == 2 ==> y, m == 3 ==> z
 END_DOC

 integer :: ipoint,i,j,m
 provide gauss_erfc_mu_r12_inv_r12_rk_tmp

 double precision :: wall0, wall1
 call wall_time(wall0)
 do i = 1, ao_num
  do j = 1, ao_num
   do ipoint = 1, n_points_final_grid
    do m = 1, 4
     gauss_erfc_mu_r12_inv_r12_rk(ipoint,j,i,m)= gauss_erfc_mu_r12_inv_r12_rk_tmp(m,ipoint,j,i)
    enddo
   enddo
  enddo
 enddo

 call wall_time(wall1)
 print*,'wall time for gauss_erfc_mu_r12_inv_r12_rk = ',wall1 - wall0
 FREE gauss_erfc_mu_r12_inv_r12_rk_tmp


END_PROVIDER


subroutine test_erfc_gauss
 implicit none
 integer :: ipoint,i,j,m,jpoint
 double precision :: r1(3),r12_sq
 double precision :: C_center(3), weight1,mu,r12,integral,delta, num_int,contrib
 double precision, allocatable :: ao_mat(:,:,:), aos_array_r1(:),ao_mat_xyz(:,:,:,:)
 allocate(ao_mat(ao_num,ao_num,n_points_final_grid), aos_array_r1(ao_num))
 allocate(ao_mat_xyz(4,ao_num,ao_num,n_points_final_grid))
 ao_mat_xyz = 0.d0
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
   r12_sq = (r1(1) - C_center(1))**2.d0 + (r1(2) - C_center(2))**2.d0 + (r1(3) - C_center(3))**2.d0 
   r12 = dsqrt(r12)
   if(r12.lt.1.d-10)cycle
   do i = 1, ao_num
    do j = i, i
     contrib = aos_array_r1(i) * aos_array_r1(j) * weight1 * dexp(-delta*r12_sq) * (1.d0 - derf(mu * r12))/r12
     do m = 1, 3
      ao_mat_xyz(m,j,i,jpoint)  += contrib * r1(m)
     enddo
     ao_mat_xyz(4,j,i,jpoint) += contrib
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
   do j = i,i
    do m = 1, 4
     integral =  gauss_erfc_mu_r12_inv_r12_rk_tmp(m,j,i,jpoint)
     num_int  = ao_mat_xyz(m,j,i,jpoint) 
     if(dabs(num_int).gt.1.d-10)then
      accu1relat = dabs(integral - num_int)/dabs(num_int)
     endif
     if(dabs(integral - num_int).gt.1.d-5)then
      print*,'m = ',m
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
