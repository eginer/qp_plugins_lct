BEGIN_PROVIDER [ double precision, ao_prod_on_grid, (n_points_final_grid,ao_num,ao_num)]
 implicit none
 integer :: ipoint, i,j,k,l,m
 double precision :: weight
 do i = 1, ao_num
  do j = 1, ao_num
   do ipoint = 1, n_points_final_grid
    weight = final_weight_at_r_vector(ipoint)
    ao_prod_on_grid(ipoint,j,i) = aos_in_r_array_transp(ipoint,j) * aos_in_r_array_transp(ipoint,i) * weight
   enddo
  enddo
 enddo

END_PROVIDER 


BEGIN_PROVIDER [ double precision, ao_prod_on_grid_transp, (ao_num,ao_num,n_points_final_grid)]
 implicit none
 integer :: ipoint, i,j,k,l,m
 double precision :: weight
 do ipoint = 1, n_points_final_grid
  weight = final_weight_at_r_vector(ipoint)
  do i = 1, ao_num
   do j = 1, ao_num
    ao_prod_on_grid_transp(j,i,ipoint) = aos_in_r_array(j,ipoint) * aos_in_r_array(i,ipoint) * weight
   enddo
  enddo
 enddo

END_PROVIDER 


BEGIN_PROVIDER [ double precision, big_array_naive, (ao_num, ao_num, ao_num, ao_num)]
 implicit none
 include 'constants.include.F'
 integer :: ipoint, i,j,k,l,m
 double precision :: mu
 big_array_naive = 0.d0

!
! ! erf(mu(r) r12)/r12
! do ipoint = 1, n_points_final_grid
!  do l = 1, ao_num
!   do j = 1, ao_num
!    do k = 1, ao_num
!     do i = 1, ao_num
!      big_array_naive(i,k,j,l) += ao_prod_on_grid_transp(i,k,ipoint) * erf_mu_r12_inv_r12_rk(j,l,ipoint)
!     enddo
!    enddo
!   enddo
!  enddo
! enddo
!
! !  mu(r) / sqpi * exp(-(mu(r) r12)^2)
! do ipoint = 1, n_points_final_grid
!  mu = mu_of_r_for_ints(ipoint,1)
!  do l = 1, ao_num
!   do j = 1, ao_num
!    do k = 1, ao_num
!     do i = 1, ao_num
!      big_array_naive(i,k,j,l) += ao_prod_on_grid_transp(i,k,ipoint) * gauss_ij_rk(j,l,ipoint) * mu * inv_sq_pi
!     enddo
!    enddo
!   enddo
!  enddo
! enddo
!
!  ! - 1/4 (1 - erf(mu(r1)r12))^2
!  do ipoint = 1, n_points_final_grid
!   mu = mu_of_r_for_ints(ipoint,1)
!   do l = 1, ao_num
!    do j = 1, ao_num
!     do k = 1, ao_num
!      do i = 1, ao_num
!       big_array_naive(i,k,j,l) += ao_prod_on_grid_transp(i,k,ipoint) * erf_mu_squared_ij_rk(j,l,ipoint) * (-0.25d0 )
!      enddo
!     enddo
!    enddo
!   enddo
!  enddo
!
 if(.not.constant_mu)then
 !  -1/2 * (grad mu(r1))^2/(4 * pi * mu(r1)^4)
  double precision :: min_inv_8_pi
  min_inv_8_pi = -1.d0/(8.d0 * pi)
  do ipoint = 1, n_points_final_grid
   mu = mu_of_r_for_ints(ipoint,1)
   do l = 1, ao_num
    do j = 1, ao_num
     do k = 1, ao_num
      do i = 1, ao_num
       big_array_naive(i,k,j,l) += ao_prod_on_grid_transp(i,k,ipoint) * gauss_2_ij_rk(j,l,ipoint) & 
                                 * min_inv_8_pi * grad_sq_mu_of_r_for_ints(ipoint,1) * inv_4_mu_of_r_for_ints(ipoint,1)
!       big_array_naive(i,k,j,l) += ao_prod_on_grid_transp(i,k,ipoint) * gauss_2_ij_rk(j,l,ipoint)  
!                                 * min_inv_8_pi * grad_sq_mu_of_r_for_ints(ipoint,1) * inv_4_mu_of_r_for_ints(ipoint,1)
      enddo
     enddo
    enddo
   enddo
  enddo
 endif



END_PROVIDER 



