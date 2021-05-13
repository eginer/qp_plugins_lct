BEGIN_PROVIDER [ double precision, gauss_ij_rk,  ( ao_num, ao_num,n_points_final_grid)]
 implicit none
 BEGIN_DOC
! mu_of_r_gauss(j,i,R) = int dr phi_i(r) phi_j(r) exp(-\mu(R) |r - R|^2)
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
  provide mu_of_r_prov
 !$OMP PARALLEL                  &
 !$OMP DEFAULT (NONE)            &
 !$OMP PRIVATE (i,j,ipoint,mu,r,int_mu,delta) & 
 !$OMP SHARED (ao_num,n_points_final_grid,mu_of_r_prov,gauss_ij_rk,final_grid_points)
 !$OMP DO SCHEDULE (dynamic)
 do ipoint = 1, n_points_final_grid
   do i = 1, ao_num
    do j = i, ao_num
     mu = mu_of_r_prov(ipoint,1)
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

subroutine test_gauss_ij_rk
 implicit none
 integer :: ipoint,i,j,m,jpoint
 double precision :: r1(3)
 double precision :: C_center(3), weight1,mu,r12,integral,delta
 double precision, allocatable :: ao_mat(:,:,:), aos_array_r1(:)
 allocate(ao_mat(ao_num,ao_num,n_points_final_grid), aos_array_r1(ao_num))
 do jpoint = 1, n_points_final_grid
  mu = mu_of_r_prov(jpoint,1)
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
    if(dabs(ao_mat(j,i,jpoint)).gt.1.d-10)then
     accu1relat = dabs(integral - ao_mat(j,i,jpoint))/dabs(ao_mat(j,i,jpoint))
    endif
    if(dabs(integral - ao_mat(j,i,jpoint)).gt.1.d-5)then
     print*,'i,j,jpoint',i,j,jpoint
     print*,'prov, num, difference'
     print*,integral,ao_mat(j,i,jpoint),dabs(integral - ao_mat(j,i,jpoint))
    endif
    accu1 += dabs(integral - ao_mat(j,i,jpoint))
   enddo
  enddo
  print*,'accu1      = ',accu1/dble(ao_num * ao_num)
  print*,''
  print*,'accu1relat = ',accu1relat/dble(ao_num * ao_num)
 enddo

end
