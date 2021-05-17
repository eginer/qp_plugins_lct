subroutine fit_erf_mu_squared(mu, expo_fit, coef_fit)
 implicit none
 BEGIN_DOC
! fit of (erf(mu * |r - r'|) - 1)^2 with a linear combination of Gaussians
 END_DOC
 double precision, intent(in) :: mu
 double precision, intent(out):: expo_fit(n_max_fit_slat)
 double precision, intent(out):: coef_fit(n_max_fit_slat)
 integer :: i
 double precision :: expos(n_max_fit_slat),alpha,beta
 alpha = 2.d0 * expos_slat_gauss_1_erf_x(1) * mu
 call expo_fit_slater_gam(alpha,expos)
 beta = 2.d0 * expos_slat_gauss_1_erf_x(2) * mu**2.d0
 do i = 1, n_max_fit_slat
  expo_fit(i) = expos(i) + beta
  coef_fit(i) = coef_fit_slat_gauss(i)
 enddo
end

BEGIN_PROVIDER [double precision, erf_mu_squared_ij_rk,( ao_num, ao_num,n_points_final_grid)]
 implicit none
 BEGIN_DOC
! erf_mu_squared_ij_rk(i,j,r) = \int dr' phi_i(r') phi_j(r') (erf(mu(r) * |r - r'|) - 1)^2
!
! it is fitted with fit_erf_mu_squared
 END_DOC
 double precision :: expo_fit(n_max_fit_slat)
 double precision :: coef_fit(n_max_fit_slat),overlap_gauss_r12_ao
 double precision :: mu,r(3),int_mu,delta,wall0,wall1
 integer :: i,j,ipoint,ifit
  provide mu_of_r_for_ints
  call fit_erf_mu_squared(1.d0, expo_fit, coef_fit)
 call wall_time(wall0)
  !$OMP PARALLEL                  &
  !$OMP DEFAULT (NONE)            &
  !$OMP PRIVATE (i,j,ipoint,mu,r,int_mu,delta,ifit,expo_fit,coef_fit) & 
  !$OMP SHARED (ao_num,n_points_final_grid,mu_of_r_for_ints,erf_mu_squared_ij_rk,final_grid_points,n_max_fit_slat)
  !$OMP DO SCHEDULE (dynamic)
  do ipoint = 1, n_points_final_grid
   do i = 1, ao_num
    do j = i, ao_num
     mu = mu_of_r_for_ints(ipoint,1)
     call fit_erf_mu_squared(mu, expo_fit, coef_fit)
     r(1) = final_grid_points(1,ipoint)
     r(2) = final_grid_points(2,ipoint)
     r(3) = final_grid_points(3,ipoint)
     erf_mu_squared_ij_rk(j,i,ipoint) = 0.d0
     do ifit = 1, n_max_fit_slat
      delta = expo_fit(ifit)
      int_mu = overlap_gauss_r12_ao(r,delta,i,j)
      erf_mu_squared_ij_rk(j,i,ipoint) += coef_fit(ifit) * int_mu 
     enddo
    enddo
   enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL


 do ipoint = 1, n_points_final_grid
  do i = 1, ao_num
   do j = 1, i-1
    erf_mu_squared_ij_rk(j,i,ipoint)= erf_mu_squared_ij_rk(i,j,ipoint)
   enddo
  enddo
 enddo

 call wall_time(wall1)
 print*,'wall time for erf_mu_squared_ij_rk  ',wall1 - wall0

END_PROVIDER 

subroutine test_erf_mu_squared_ij_rk
 implicit none
 integer :: ipoint,i,j,m,jpoint
 double precision :: r1(3)
 double precision :: C_center(3), weight1,mu,r12,integral,delta, num_int
 double precision, allocatable :: ao_mat(:,:,:), aos_array_r1(:)
 allocate(ao_mat(ao_num,ao_num,n_points_final_grid), aos_array_r1(ao_num))
 do jpoint = 1, n_points_final_grid
  mu = mu_of_r_for_ints(jpoint,1)
  C_center(:) = final_grid_points(:,jpoint)
  ao_mat(:,:,jpoint) = 0.d0
  do ipoint = 1, n_points_final_grid
   r1(1) = final_grid_points(1,ipoint)
   r1(2) = final_grid_points(2,ipoint)
   r1(3) = final_grid_points(3,ipoint)
   call give_all_aos_at_r(r1,aos_array_r1)
   weight1 = final_weight_at_r_vector(ipoint)
   r12 = (r1(1) - C_center(1))**2.d0 + (r1(2) - C_center(2))**2.d0 + (r1(3) - C_center(3))**2.d0 
   r12 = dsqrt(r12)
   do i = 1, ao_num
    do j = 1, ao_num
     ao_mat(j,i,jpoint)  += aos_array_r1(i) * aos_array_r1(j) * weight1 * (1.d0 - derf(mu * r12))**2.d0
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
    integral =  erf_mu_squared_ij_rk(j,i,jpoint)
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

end
