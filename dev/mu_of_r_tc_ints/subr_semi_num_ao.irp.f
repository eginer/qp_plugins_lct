BEGIN_PROVIDER [ double precision, scalar_mu_r_pot_chemist_ao, (ao_num, ao_num, ao_num, ao_num)]
 implicit none
 BEGIN_DOC
! scalar_mu_r_pot_chemist_ao(i,k,j,l) = \int dr1 \int dr2 \phi_i(r1) \phi_k(r1) W_ee(r12,\mu(r1)) \phi_j(r2) \phi_l(r2)
!
! NOTICE THAT : because of mu(r1) and the non symmetric form of the Jastrow factor, the integrals ARE NOT SYMMETRIC in r1, r2 
!
! scalar_mu_r_pot_chemist_ao(i,k,j,l) NOT EQUAL TO scalar_mu_r_pot_chemist_ao(j,l,i,k) for instance 
 END_DOC

 print*,''
 print*,'providing scalar_mu_r_pot_chemist_ao ...'
 print*,''
 print*,''
 call wall_time(wall0)
 double precision :: wall0,wall1


 scalar_mu_r_pot_chemist_ao = 0.d0

 call all_erf_mu_r1_lr_int_big_mat(scalar_mu_r_pot_chemist_ao)
 call all_gauss_r12_big_mat(scalar_mu_r_pot_chemist_ao)
 call all_erf_r12_sq_big_mat(scalar_mu_r_pot_chemist_ao)
 call all_nabla_mu_r1_sq_big_mat(scalar_mu_r_pot_chemist_ao)
 call all_nabla_mu_r1_cdot_r12_big_mat(scalar_mu_r_pot_chemist_ao)
 call all_nabla_mu_r1_cdot_r12_erfc_r12_big_mat(scalar_mu_r_pot_chemist_ao)
 call wall_time(wall1)
 print*,''
 print*,''
 print*,'wall time for calar_mu_r_pot_chemist_ao ', wall1 - wall0   
 print*,''
 print*,''
END_PROVIDER 




subroutine all_gauss_r12_big_mat(big_mat)
 implicit none
 BEGIN_DOC
! enters with the array big_mat and add the following integrals 
!
! big_mat(i,k,j,l) += \int dr1 \phi_i(r1) \phi_k(r1) \mu(r1)/sqpi \int dr2 exp[-(mu(r1)r12)^2] \phi_j(r2) \phi_l(r2)
 END_DOC
 include 'constants.include.F'
 integer :: i,k,ipoint,m
 double precision :: weight,mu
 double precision, allocatable :: a_mat(:,:,:)
 double precision, intent(inout) :: big_mat(ao_num,ao_num,ao_num,ao_num)

 print*,'computing all_gauss_r12_big_mat ...'
 call wall_time(wall0)
 double precision :: wall0,wall1

 allocate(a_mat(n_points_final_grid,ao_num,ao_num))

 ! mu(r) / sqpi * exp[-(mu(r1)r12)^2]
  do i = 1, ao_num
   do k = 1, ao_num
    do ipoint = 1, n_points_final_grid
     mu = mu_of_r_for_ints(ipoint,1)
     weight = final_weight_at_r_vector(ipoint)
     a_mat(ipoint,k,i) = aos_in_r_array_transp(ipoint,k) * aos_in_r_array_transp(ipoint,i) * weight * mu * inv_sq_pi
    enddo
   enddo
  enddo
 call dgemm("N","N",ao_num*ao_num,ao_num*ao_num,n_points_final_grid,1.d0,gauss_ij_rk(1,1,1),ao_num*ao_num & 
                   ,a_mat(1,1,1),n_points_final_grid,1.d0,big_mat,ao_num*ao_num)

  call wall_time(wall1)
  print*,'wall time for all_gauss_r12_big_mat ',wall1 - wall0   
end

subroutine all_erf_r12_sq_big_mat(big_mat)
 implicit none
 BEGIN_DOC
! enters with the array big_mat and add the following integrals 
!
! big_mat(i,k,j,l) += -1/4 \int dr1 \phi_i(r1) \phi_k(r1) \int dr2 [ 1 - erf(\mu(r1)r12)]^2 \phi_j(r2) \phi_l(r2)
 END_DOC
 include 'constants.include.F'
 integer :: i,k,ipoint,m
 double precision :: weight,mu
 double precision, allocatable :: a_mat(:,:,:)
 double precision, intent(inout) :: big_mat(ao_num,ao_num,ao_num,ao_num)

 print*,'computing all_erf_r12_sq_big_mat ...'
 call wall_time(wall0)
 double precision :: wall0,wall1

 allocate(a_mat(n_points_final_grid,ao_num,ao_num))

  ! - 1/4 (1 - erf(mu(r1)r12))^2
  do i = 1, ao_num
   do k = 1, ao_num
    do ipoint = 1, n_points_final_grid
     mu = mu_of_r_for_ints(ipoint,1)
     weight = final_weight_at_r_vector(ipoint)
     a_mat(ipoint,k,i) = aos_in_r_array_transp(ipoint,k) * aos_in_r_array_transp(ipoint,i) * weight * (-0.25d0)
    enddo
   enddo
  enddo
 call dgemm("N","N",ao_num*ao_num,ao_num*ao_num,n_points_final_grid,1.d0,erf_mu_squared_ij_rk(1,1,1),ao_num*ao_num & 
                   ,a_mat(1,1,1),n_points_final_grid,1.d0,big_mat,ao_num*ao_num)

  call wall_time(wall1)
  print*,'wall time for all_erf_r12_sq_big_mat ',wall1 - wall0   
end

subroutine all_nabla_mu_r1_sq_big_mat(big_mat)
 implicit none
 BEGIN_DOC
! enters with the array big_mat and add the following integrals 
!
! big_mat(i,k,j,l) += -1/2 \int dr1 \phi_i(r1) \phi_k(r1) (grad mu(r1))^2/(4 * pi * mu(r1)^4) \int dr2 exp[-2 (\mu(r1) r12)^2] \phi_j(r2) \phi_l(r2)
 END_DOC
 include 'constants.include.F'
 integer :: i,k,ipoint,m
 double precision :: weight,mu
 double precision, allocatable :: a_mat(:,:,:)
 double precision, intent(inout) :: big_mat(ao_num,ao_num,ao_num,ao_num)

 print*,'computing all_nabla_mu_r1_sq_big_mat ...'
 call wall_time(wall0)
 double precision :: wall0,wall1

 allocate(a_mat(n_points_final_grid,ao_num,ao_num))

 !  -1/2 * (grad mu(r1))^2/(4 * pi * mu(r1)^4)
  double precision :: min_inv_8_pi
  min_inv_8_pi = -1.d0/(8.d0 * pi)
  do i = 1, ao_num
   do k = 1, ao_num
    do ipoint = 1, n_points_final_grid
     mu = mu_of_r_for_ints(ipoint,1)
     weight = final_weight_at_r_vector(ipoint)
     a_mat(ipoint,k,i) = aos_in_r_array_transp(ipoint,k) * aos_in_r_array_transp(ipoint,i) * weight & 
                       * grad_sq_mu_of_r_for_ints(ipoint,1) * inv_4_mu_of_r_for_ints(ipoint,1) * min_inv_8_pi
    enddo
   enddo
  enddo
  call dgemm("N","N",ao_num*ao_num,ao_num*ao_num,n_points_final_grid,1.d0,gauss_2_ij_rk(1,1,1),ao_num*ao_num & 
                    ,a_mat(1,1,1),n_points_final_grid,1.d0,big_mat,ao_num*ao_num)

  call wall_time(wall1)
  print*,'wall time for all_nabla_mu_r1_sq_big_mat ',wall1 - wall0   
end

subroutine all_nabla_mu_r1_cdot_r12_big_mat(big_mat)
 implicit none
 BEGIN_DOC
! enters with the array big_mat and add the following integrals 
!
! big_mat(i,k,j,l) += 1/(2 * sqpi) \int dr1 \int dr2 \phi_i(r1) \phi_k(r1) \nabla_1 . cdot  (r1 - r2) exp[-2 (\mu(r1) r12)^2] \phi_j(r2) \phi_l(r2)
 END_DOC
 include 'constants.include.F'
 integer :: i,k,ipoint,m
 double precision :: weight,mu
 double precision, allocatable :: a_mat(:,:,:)
 double precision, intent(inout) :: big_mat(ao_num,ao_num,ao_num,ao_num)
 print*,'computing all_nabla_mu_r1_cdot_r12_big_mat ...'
 call wall_time(wall0)
 double precision :: wall0,wall1

 allocate(a_mat(n_points_final_grid,ao_num,ao_num))

! ! + 1/(2 * sqrt(pi)) e^{-(\mu(r1)r12)^2} (r_1 - r_2) . \nabla_1 \mu(r1)
  double precision :: inv_2_sqpi,r1(3)
  inv_2_sqpi = 0.5d0 * inv_sq_pi
  do m = 1, 3 
   ! first part : 1/(2 * sqrt(pi)) e^{-(\mu(r1)r12)^2} * x1 \deriv{x1} \mu(r1) 
   do i = 1, ao_num
    do k = 1, ao_num
     do ipoint = 1, n_points_final_grid
      r1(:) = final_grid_points(:,ipoint)
      weight = final_weight_at_r_vector(ipoint)
      a_mat(ipoint,k,i) = aos_in_r_array_transp(ipoint,k) * aos_in_r_array_transp(ipoint,i) * weight & 
                        * inv_2_sqpi * r1(m) * grad_mu_of_r_for_ints(m,ipoint,1) 
     enddo
    enddo
   enddo

   call dgemm("N","N",ao_num*ao_num,ao_num*ao_num,n_points_final_grid,1.d0,gauss_ij_rk(1,1,1),ao_num*ao_num & 
                     ,a_mat(1,1,1),n_points_final_grid,1.d0,big_mat,ao_num*ao_num)

   ! second part : - 1/(2 * sqrt(pi)) e^{-(\mu(r1)r12)^2} * x2 \deriv{x1} \mu(r1) 
   do i = 1, ao_num
    do k = 1, ao_num
     do ipoint = 1, n_points_final_grid
      r1(:) = final_grid_points(:,ipoint)
      weight = final_weight_at_r_vector(ipoint)
      a_mat(ipoint,k,i) = aos_in_r_array_transp(ipoint,k) * aos_in_r_array_transp(ipoint,i) * weight & 
                        * (-1.d0) * inv_2_sqpi * grad_mu_of_r_for_ints(m,ipoint,1) 
     enddo
    enddo
   enddo

   call dgemm("N","N",ao_num*ao_num,ao_num*ao_num,n_points_final_grid,1.d0,gauss_ij_xyz_rk(1,1,1,m),ao_num*ao_num & 
                     ,a_mat(1,1,1),n_points_final_grid,1.d0,big_mat,ao_num*ao_num)
  enddo

  call wall_time(wall1)
  print*,'wall time for all_nabla_mu_r1_cdot_r12_big_mat ',wall1 - wall0   
end

subroutine all_nabla_mu_r1_cdot_r12_erfc_r12_big_mat(big_mat)
 implicit none
 BEGIN_DOC
! enters with the array big_mat and add the following integrals 
!
! big_mat(i,k,j,l) += - 1/(4 * sqpi) \int dr1 \int dr2 \phi_i(r1) \phi_k(r1) 1 / mu(r1)^2 * \nabla_1 . cdot  (r1 - r2) exp[-2 (\mu(r1) r12)^2] (1 - erf(mu*r12))/r12 * \phi_j(r2) \phi_l(r2)
 END_DOC
 include 'constants.include.F'
 integer :: i,k,ipoint,m
 double precision :: weight,mu
 double precision, allocatable :: a_mat(:,:,:)
 double precision, intent(inout) :: big_mat(ao_num,ao_num,ao_num,ao_num)
 print*,'computing all_nabla_mu_r1_cdot_r12_erfc_r12_big_mat ...'
 call wall_time(wall0)
 double precision :: wall0,wall1
 allocate(a_mat(n_points_final_grid,ao_num,ao_num))

! ! - 1/(4 * sqrt(pi)) (1 - erf(mu * r12))/r12 * e^{-(\mu(r1)r12)^2} (r_1 - r_2) . \nabla_1 \mu(r1)
  double precision :: inv_2_sqpi,r1(3)
  inv_2_sqpi = -0.25d0 * inv_sq_pi
  do m = 1, 3 
   ! first part : 1/(4 * sqrt(pi)) * (1 - erf(mu * r12))/r12 e^{-(\mu(r1)r12)^2} * x1 \deriv{x1} \mu(r1) 
   do i = 1, ao_num
    do k = 1, ao_num
     do ipoint = 1, n_points_final_grid
      r1(:) = final_grid_points(:,ipoint)
      weight = final_weight_at_r_vector(ipoint)
      a_mat(ipoint,k,i) = aos_in_r_array_transp(ipoint,k) * aos_in_r_array_transp(ipoint,i) * weight & 
                        * inv_2_sqpi * inv_2_mu_of_r_for_ints(ipoint,1) *  r1(m) * grad_mu_of_r_for_ints(m,ipoint,1) 
     enddo
    enddo
   enddo

   call dgemm("N","N",ao_num*ao_num,ao_num*ao_num,n_points_final_grid,1.d0,gauss_erfc_mu_r12_inv_r12_rk(1,1,1,4),ao_num*ao_num & 
                     ,a_mat(1,1,1),n_points_final_grid,1.d0,big_mat,ao_num*ao_num)

   ! second part : - 1/(4 * sqrt(pi)) e^{-(\mu(r1)r12)^2} * (1 - erf(mu * r12))/r12 * x2 \deriv{x1} \mu(r1) 
   do i = 1, ao_num
    do k = 1, ao_num
     do ipoint = 1, n_points_final_grid
      r1(:) = final_grid_points(:,ipoint)
      weight = final_weight_at_r_vector(ipoint)
      a_mat(ipoint,k,i) = aos_in_r_array_transp(ipoint,k) * aos_in_r_array_transp(ipoint,i) * weight & 
                        * (-1.d0) * inv_2_sqpi * inv_2_mu_of_r_for_ints(ipoint,1) * grad_mu_of_r_for_ints(m,ipoint,1) 
     enddo
    enddo
   enddo

   call dgemm("N","N",ao_num*ao_num,ao_num*ao_num,n_points_final_grid,1.d0,gauss_erfc_mu_r12_inv_r12_rk(1,1,1,m),ao_num*ao_num & 
                     ,a_mat(1,1,1),n_points_final_grid,1.d0,big_mat,ao_num*ao_num)
  enddo
  call wall_time(wall1)
  print*,'wall time for all_nabla_mu_r1_cdot_r12_erfc_r12_big_mat ',wall1 - wall0   

end
