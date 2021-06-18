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

!!! Laplacian of Jastrow: 
!! !call lr_int_mu_r1(scalar_mu_r_pot_chemist_ao) ! other approch for erf(mu(r1) r12)/r12 
 call lr_int_mu_r1_bis(scalar_mu_r_pot_chemist_ao)! Long range interaction erf(mu(r1) r12)/r12
 call gauss_int_mu_r1(scalar_mu_r_pot_chemist_ao) ! Gaussian e^{(-mu(r1) r12)^2}
 call gauss_grad_scal_r12_int_mu_r1(scalar_mu_r_pot_chemist_ao) ! e^{(-mu(r1) r12)^2} r12 . \grad mu(r1)
 call lapl_gamm_r1(scalar_mu_r_pot_chemist_ao) ! \grad . gamma(r1) 

!!! (Nabla of Jastrow)^2
 call erf_sq_int_mu_r1(scalar_mu_r_pot_chemist_ao) ! erf(mu(r1) r12)^2
 call gauss_sq_int_mu_r1(scalar_mu_r_pot_chemist_ao) ! e^{-2(mu(r1)r12)^2} \grad mu(r1)^2
 call erfc_gauss_int_mu_r1(scalar_mu_r_pot_chemist_ao) ! erfc(mu(r1)r12)/r12 * e^{(-mu(r1) r12)^2}
!!! 
 call wall_time(wall1)
 print*,''
 print*,''
 print*,'wall time for calar_mu_r_pot_chemist_ao ', wall1 - wall0   
 print*,''
 print*,''
END_PROVIDER 




subroutine gauss_int_mu_r1(big_mat)
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

 print*,'computing gauss_int_mu_r1 ...'
 call wall_time(wall0)
 double precision :: wall0,wall1

 allocate(a_mat(ao_num,ao_num,n_points_final_grid))

 ! mu(r) / sqpi * exp[-(mu(r1)r12)^2]
  do ipoint = 1, n_points_final_grid
   mu = mu_of_r_for_ints(ipoint,1)
   weight = final_weight_at_r_vector(ipoint)
   do k = 1, ao_num
    do i = 1, ao_num
     a_mat(i,k,ipoint) = aos_in_r_array(k,ipoint) * aos_in_r_array(i,ipoint) * weight * mu * inv_sq_pi
    enddo
   enddo
  enddo
 call dgemm("N","N",ao_num*ao_num,ao_num*ao_num,n_points_final_grid,1.d0,a_mat(1,1,1),ao_num*ao_num & 
                   ,gauss_ij_rk(1,1,1),n_points_final_grid,1.d0,big_mat,ao_num*ao_num)

  call wall_time(wall1)
  print*,'wall time for gauss_int_mu_r1 ',wall1 - wall0   
end

subroutine erf_sq_int_mu_r1(big_mat)
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

 print*,'computing erf_sq_int_mu_r1 ...'
 call wall_time(wall0)
 double precision :: wall0,wall1

 allocate(a_mat(ao_num,ao_num,n_points_final_grid))

  ! - 1/4 (1 - erf(mu(r1)r12))^2
 do ipoint = 1, n_points_final_grid
  weight = final_weight_at_r_vector(ipoint)
  mu = mu_of_r_for_ints(ipoint,1)
  do k = 1, ao_num
   do i = 1, ao_num
     a_mat(i,k,ipoint) = aos_in_r_array(k,ipoint) * aos_in_r_array(i,ipoint) * weight * (-0.25d0)
    enddo
   enddo
  enddo
 call dgemm("N","N",ao_num*ao_num,ao_num*ao_num,n_points_final_grid,1.d0,a_mat(1,1,1),ao_num*ao_num & 
                   ,erf_mu_squared_ij_rk(1,1,1),n_points_final_grid,1.d0,big_mat,ao_num*ao_num)

  call wall_time(wall1)
  print*,'wall time for erf_sq_int_mu_r1 ',wall1 - wall0   
end

subroutine gauss_sq_int_mu_r1(big_mat)
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

 print*,'computing gauss_sq_int_mu_r1 ...'
 call wall_time(wall0)
 double precision :: wall0,wall1

 allocate(a_mat(ao_num,ao_num,n_points_final_grid))

 !  -1/2 * (grad mu(r1))^2/(4 * pi * mu(r1)^4)
  double precision :: min_inv_8_pi
  min_inv_8_pi = -1.d0/(8.d0 * pi)
  do ipoint = 1, n_points_final_grid
   mu = mu_of_r_for_ints(ipoint,1)
   weight = final_weight_at_r_vector(ipoint)
   do k = 1, ao_num
    do i = 1, ao_num
     a_mat(i,k,ipoint) = aos_in_r_array(k,ipoint) * aos_in_r_array(i,ipoint) * weight & 
                       * grad_sq_mu_of_r_for_ints(ipoint,1) * inv_4_mu_of_r_for_ints(ipoint,1) * min_inv_8_pi
    enddo
   enddo
  enddo
  call dgemm("N","N",ao_num*ao_num,ao_num*ao_num,n_points_final_grid,1.d0,a_mat(1,1,1),ao_num*ao_num & 
                    ,gauss_2_ij_rk(1,1,1),n_points_final_grid,1.d0,big_mat,ao_num*ao_num)

  call wall_time(wall1)
  print*,'wall time for gauss_sq_int_mu_r1 ',wall1 - wall0   
end

subroutine gauss_grad_scal_r12_int_mu_r1(big_mat)
 implicit none
 BEGIN_DOC
! enters with the array big_mat and add the following integrals 
!
! big_mat(i,k,j,l) += 1/(2 * sqpi) \int dr1 \int dr2 \phi_i(r1) \phi_k(r1) \nabla_1 \mu(r1). cdot  (r1 - r2) exp[- (\mu(r1) r12)^2] \phi_j(r2) \phi_l(r2)
 END_DOC
 include 'constants.include.F'
 integer :: i,k,ipoint,m
 double precision :: weight,mu
 double precision, allocatable :: a_mat(:,:,:)
 double precision, intent(inout) :: big_mat(ao_num,ao_num,ao_num,ao_num)
 print*,'computing gauss_grad_scal_r12_int_mu_r1 ...'
 call wall_time(wall0)
 double precision :: wall0,wall1

 allocate(a_mat(ao_num,ao_num,n_points_final_grid))

! ! + 1/(2 * sqrt(pi)) e^{-(\mu(r1)r12)^2} (r_1 - r_2) . \nabla_1 \mu(r1)
  double precision :: inv_2_sqpi,r1(3)
  inv_2_sqpi = 0.5d0 * inv_sq_pi
  do m = 1, 3 
   ! first part : 1/(2 * sqrt(pi)) e^{-(\mu(r1)r12)^2} * x1 \deriv{x1} \mu(r1) 
   do ipoint = 1, n_points_final_grid
    r1(:) = final_grid_points(:,ipoint)
    weight = final_weight_at_r_vector(ipoint)
    do k = 1, ao_num
     do i = 1, ao_num
      a_mat(i,k,ipoint) = aos_in_r_array(k,ipoint) * aos_in_r_array(i,ipoint) * weight & 
                        * inv_2_sqpi * r1(m) * grad_mu_of_r_for_ints(m,ipoint,1) 
     enddo
    enddo
   enddo

   call dgemm("N","N",ao_num*ao_num,ao_num*ao_num,n_points_final_grid,1.d0,a_mat(1,1,1),ao_num*ao_num & 
                     ,gauss_ij_rk(1,1,1),n_points_final_grid,1.d0,big_mat,ao_num*ao_num)

   ! second part : - 1/(2 * sqrt(pi)) e^{-(\mu(r1)r12)^2} * x2 \deriv{x1} \mu(r1) 
   do ipoint = 1, n_points_final_grid
    r1(:) = final_grid_points(:,ipoint)
    weight = final_weight_at_r_vector(ipoint)
    do k = 1, ao_num
     do i = 1, ao_num
      a_mat(i,k,ipoint) = aos_in_r_array(k,ipoint) * aos_in_r_array(i,ipoint) * weight & 
                        * (-1.d0) * inv_2_sqpi * grad_mu_of_r_for_ints(m,ipoint,1) 
     enddo
    enddo
   enddo

   call dgemm("N","N",ao_num*ao_num,ao_num*ao_num,n_points_final_grid,1.d0,a_mat(1,1,1),ao_num*ao_num & 
                     ,gauss_ij_xyz_rk(1,1,1,m),n_points_final_grid,1.d0,big_mat,ao_num*ao_num)
  enddo

  call wall_time(wall1)
  print*,'wall time for gauss_grad_scal_r12_int_mu_r1 ',wall1 - wall0   
end

subroutine erfc_gauss_int_mu_r1(big_mat)
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
 print*,'computing erfc_gauss_int_mu_r1 ...'
 call wall_time(wall0)
 double precision :: wall0,wall1
 allocate(a_mat(ao_num,ao_num,n_points_final_grid))

! ! - 1/(4 * sqrt(pi)) (1 - erf(mu * r12))/r12 * e^{-(\mu(r1)r12)^2} (r_1 - r_2) . \nabla_1 \mu(r1)
  double precision :: inv_2_sqpi,r1(3)
  inv_2_sqpi = -0.25d0 * inv_sq_pi
  do m = 1, 3 
   ! first part : 1/(4 * sqrt(pi)) * (1 - erf(mu * r12))/r12 e^{-(\mu(r1)r12)^2} * x1 \deriv{x1} \mu(r1) 
   do ipoint = 1, n_points_final_grid
    r1(:) = final_grid_points(:,ipoint)
    weight = final_weight_at_r_vector(ipoint)
    do i = 1, ao_num
     do k = 1, ao_num
      a_mat(k,i,ipoint) = aos_in_r_array(k,ipoint) * aos_in_r_array(i,ipoint) * weight & 
                        * inv_2_sqpi * inv_2_mu_of_r_for_ints(ipoint,1) *  r1(m) * grad_mu_of_r_for_ints(m,ipoint,1) 
     enddo
    enddo
   enddo

   call dgemm("N","N",ao_num*ao_num,ao_num*ao_num,n_points_final_grid,1.d0,a_mat(1,1,1),ao_num*ao_num & 
                     ,gauss_erfc_mu_r12_inv_r12_rk(1,1,1,4),n_points_final_grid,1.d0,big_mat,ao_num*ao_num)

   ! second part : - 1/(4 * sqrt(pi)) e^{-(\mu(r1)r12)^2} * (1 - erf(mu * r12))/r12 * x2 \deriv{x1} \mu(r1) 
   do ipoint = 1, n_points_final_grid
    r1(:) = final_grid_points(:,ipoint)
    weight = final_weight_at_r_vector(ipoint)
    do i = 1, ao_num
     do k = 1, ao_num
      a_mat(k,i,ipoint) = aos_in_r_array(k,ipoint) * aos_in_r_array(i,ipoint) * weight & 
                        * (-1.d0) * inv_2_sqpi * inv_2_mu_of_r_for_ints(ipoint,1) * grad_mu_of_r_for_ints(m,ipoint,1) 
     enddo
    enddo
   enddo

   call dgemm("N","N",ao_num*ao_num,ao_num*ao_num,n_points_final_grid,1.d0,a_mat(1,1,1),ao_num*ao_num & 
                     ,gauss_erfc_mu_r12_inv_r12_rk(1,1,1,m),n_points_final_grid,1.d0,big_mat,ao_num*ao_num)
  enddo
  call wall_time(wall1)
  print*,'wall time for erfc_gauss_int_mu_r1 ',wall1 - wall0   

end


subroutine test_num_scal_pot
 implicit none
 include 'constants.include.F'
 integer :: ipoint,i,j,k,l,m,jpoint
 double precision :: r1(3),r2(3),weight1,weight2,weight_tot
 double precision :: nabla_sq_term,cst_nabla,ao_prod_r1,ao_prod_r2,nabla_r12_1,erf_mu_sq,nabla_r12_2
 double precision :: thr, sq_thr,cst_nabla_r12_1,cst_gauss_r12,gauss_r12_mu_r1,erf_mur1,cst_nabla_r12_2
 double precision, allocatable :: accu(:,:,:,:)
 allocate(accu(ao_num, ao_num, ao_num, ao_num))
 thr = 1.d-15
 sq_thr = dsqrt(thr)
 cst_nabla = -1.d0/(8.d0 * pi)
 cst_nabla_r12_1 = 0.5d0 * inv_sq_pi
 cst_nabla_r12_2 = -0.25d0 * inv_sq_pi
 cst_gauss_r12 = inv_sq_pi
 accu = 0.d0

 do jpoint = 1, n_points_final_grid ! r2
  weight2 = final_weight_at_r_vector(jpoint)
  do l = 1, ao_num
!  do l = 1, 1
!   do j = 1, ao_num
   do j = l, l
    ao_prod_r2 = aos_in_r_array(j,jpoint) * aos_in_r_array(l,jpoint) * weight2
    if(dabs(ao_prod_r2).lt.sq_thr)cycle
    do ipoint = 1, n_points_final_grid ! r1 
     weight1 = final_weight_at_r_vector(ipoint)
     weight_tot = weight1 * weight2
     if(dabs(weight_tot).lt.thr)cycle
!     do i = 1, ao_num
!      do k = 1, ao_num
     do k = l, l
      do i = j, j
       ao_prod_r1 = aos_in_r_array(i,ipoint) * aos_in_r_array(k,ipoint) * weight1
       if(dabs(ao_prod_r1).lt.sq_thr)cycle
        accu(i,k,j,l) += erf_mur1(ipoint,jpoint) * ao_prod_r1 * ao_prod_r2
        accu(i,k,j,l) += gauss_r12_mu_r1(ipoint,jpoint,cst_gauss_r12) * ao_prod_r1 * ao_prod_r2
        accu(i,k,j,l) += erf_mu_sq(ipoint,jpoint) * ao_prod_r1 * ao_prod_r2
        accu(i,k,j,l) += nabla_sq_term(ipoint,jpoint) * ao_prod_r1 * ao_prod_r2
        accu(i,k,j,l) += nabla_r12_1(ipoint,jpoint,cst_nabla_r12_1) * ao_prod_r1 * ao_prod_r2
        accu(i,k,j,l) += nabla_r12_2(ipoint,jpoint,cst_nabla_r12_2) * ao_prod_r1 * ao_prod_r2
      enddo
     enddo
    enddo

   enddo
  enddo

 enddo

 double precision :: num_int,contrib,accu_naive
 accu_naive = 0.d0
! do l = 1, 1
 do l = 1, ao_num
!  do j = 1, ao_num
  do j = l, l
   do k = l, l
    do i = 1, ao_num
     num_int = accu(i,k,j,l)
     contrib = dabs(num_int - scalar_mu_r_pot_chemist_ao(i,k,j,l) )
     accu_naive += contrib
     if(contrib .gt. 1.d-10)then
      print*,''
      print*,'i,k,j,l',i,k,j,l
      print*,contrib,num_int, scalar_mu_r_pot_chemist_ao(i,k,j,l)
     endif
     
    enddo
   enddo
  enddo
 enddo

 print*,'accu naive  = ',accu_naive/dble(ao_num**2)
! print*,'accu dgemm  = ',accu2/dble(ao_num**4)
end


double precision function nabla_sq_term(ipoint,jpoint)
 implicit none
 BEGIN_DOC
! -1/(8 pi * mu(r1)^4) (\nabla_1 \mu(r1))^2 exp(-2 (\mu(r1)r12)^2)
 END_DOC
 integer, intent(in) :: ipoint ! r1 
 integer, intent(in) :: jpoint ! r1 
 double precision:: cst
 include 'constants.include.F'
 double precision :: r12,r1(3),r2(3),mu
 cst = -0.125d0 * inv_pi
 r1(:) = final_grid_points(:,ipoint)

 r2(:) = final_grid_points(:,jpoint)

 r12  = (r1(1) - r2(1))*(r1(1) - r2(1))
 r12 += (r1(2) - r2(2))*(r1(2) - r2(2))
 r12 += (r1(3) - r2(3))*(r1(3) - r2(3))
 mu = mu_of_r_for_ints(ipoint,1)
 nabla_sq_term = cst * inv_4_mu_of_r_for_ints(ipoint,1) * dexp(-2.d0 * mu * mu * r12) * grad_sq_mu_of_r_for_ints(ipoint,1)

end

double precision function erf_mu_sq(ipoint,jpoint)
 implicit none
 BEGIN_DOC
! -1/4 (erf(mu(r1) r12) - 1)^2
 END_DOC
 integer, intent(in) :: ipoint ! r1 
 integer, intent(in) :: jpoint ! r1 
 double precision :: cst
 include 'constants.include.F'
 double precision :: r12,r1(3),r2(3),mu,contrib
 cst = -0.25d0
 r1(:) = final_grid_points(:,ipoint)

 r2(:) = final_grid_points(:,jpoint)

 r12  = (r1(1) - r2(1))*(r1(1) - r2(1))
 r12 += (r1(2) - r2(2))*(r1(2) - r2(2))
 r12 += (r1(3) - r2(3))*(r1(3) - r2(3))
 r12 = dsqrt(r12)
 mu = mu_of_r_for_ints(ipoint,1)
 contrib = ( 1.d0 - derf(mu * r12)) 
 erf_mu_sq = cst * contrib * contrib 

end


double precision function nabla_r12_1(ipoint,jpoint,cst)
 implicit none 
 integer, intent(in) :: ipoint
 integer, intent(in) :: jpoint
 double precision, intent(in) :: cst
 include 'constants.include.F'
 double precision :: r12,r12_vec(3),r1(3),r2(3),mu
 r1(:) = final_grid_points(:,ipoint)

 r2(:) = final_grid_points(:,jpoint)

 r12_vec(1) = (r1(1) - r2(1))
 r12_vec(2) = (r1(2) - r2(2))
 r12_vec(3) = (r1(3) - r2(3))
 r12 = r12_vec(1)*r12_vec(1) + r12_vec(2)*r12_vec(2) + r12_vec(3)*r12_vec(3)
 mu = mu_of_r_for_ints(ipoint,1)
 nabla_r12_1 = cst * dexp(- mu * mu * r12) * ( grad_mu_of_r_for_ints(1,ipoint,1) * r12_vec(1) & 
                                             + grad_mu_of_r_for_ints(2,ipoint,1) * r12_vec(2) &
                                             + grad_mu_of_r_for_ints(3,ipoint,1) * r12_vec(3) )
end

double precision function nabla_r12_2(ipoint,jpoint,cst)
 implicit none 
 BEGIN_DOC
! - 1/(4 sqrt(pi) (\mu(r1))^2) \nabla_1 \mu(r1) . r_12  e^{(-\mu(r1)r12)^2} (1 - erf(\mu(r1)r12))/r12
 END_DOC
 integer, intent(in) :: ipoint
 integer, intent(in) :: jpoint
 double precision :: cst
 include 'constants.include.F'
 double precision :: r12,r12_sq,r12_vec(3),r1(3),r2(3),mu
 cst = -0.25 * inv_sq_pi

 nabla_r12_2 = 0.d0
 r1(:) = final_grid_points(:,ipoint)

 r2(:) = final_grid_points(:,jpoint)

 r12_vec(1) = (r1(1) - r2(1))
 r12_vec(2) = (r1(2) - r2(2))
 r12_vec(3) = (r1(3) - r2(3))
 r12_sq = r12_vec(1)*r12_vec(1) + r12_vec(2)*r12_vec(2) + r12_vec(3)*r12_vec(3)
 r12 = dsqrt(r12_sq)
 if(r12.lt.1.d-10)return
 mu = mu_of_r_for_ints(ipoint,1)
 nabla_r12_2 = cst * inv_2_mu_of_r_for_ints(ipoint,1) * dexp(- mu * mu * r12_sq) & 
             * ( grad_mu_of_r_for_ints(1,ipoint,1) * r12_vec(1) & 
               + grad_mu_of_r_for_ints(2,ipoint,1) * r12_vec(2) &
               + grad_mu_of_r_for_ints(3,ipoint,1) * r12_vec(3) ) & 
             * (1.d0 - derf(mu * r12))/r12
end

double precision function gauss_r12_mu_r1(ipoint,jpoint,cst)
 implicit none 
 integer, intent(in) :: ipoint
 integer, intent(in) :: jpoint
 double precision, intent(in) :: cst
 include 'constants.include.F'
 double precision :: r12,r1(3),r2(3),mu
 r1(:) = final_grid_points(:,ipoint)

 r2(:) = final_grid_points(:,jpoint)

 r12  = (r1(1) - r2(1))*(r1(1) - r2(1))
 r12 += (r1(2) - r2(2))*(r1(2) - r2(2))
 r12 += (r1(3) - r2(3))*(r1(3) - r2(3))
 mu = mu_of_r_for_ints(ipoint,1)

 gauss_r12_mu_r1 = cst * dexp(- mu * mu * r12) * mu 
end

double precision function erf_mur1(ipoint,jpoint)
 implicit none 
 integer, intent(in) :: ipoint
 integer, intent(in) :: jpoint
 include 'constants.include.F'
 double precision :: r12,r1(3),r2(3),mu
 double precision :: derf_mu_x
 r1(:) = final_grid_points(:,ipoint)

 r2(:) = final_grid_points(:,jpoint)

 r12  = (r1(1) - r2(1))*(r1(1) - r2(1))
 r12 += (r1(2) - r2(2))*(r1(2) - r2(2))
 r12 += (r1(3) - r2(3))*(r1(3) - r2(3))
 r12 = dsqrt(r12)
 mu = mu_of_r_for_ints(ipoint,1)
 erf_mur1 = derf_mu_x(mu,r12)

end

subroutine test_big_array_ao
 implicit none
 integer :: i,j,k,l
 double precision :: accu1, accu2,contrib,exact,ao_two_e_integral_erf,ao_two_e_integral_eff_pot
 double precision :: num,accu_n,thr,contrib_relat
 thr = 1.d-10
 accu1 = 0.d0
 accu2 = 0.d0
 accu_n = 0.d0
  do l = 1, ao_num
   do j = 1, ao_num
    do k = 1, ao_num
     do i = 1, ao_num
      exact = ao_two_e_integral_erf(i,k,j,l) 
      exact += ao_two_e_integral_eff_pot(i,k,j,l)

      num = scalar_mu_r_pot_chemist_ao(i,k,j,l)
      if(dabs(exact).lt.thr)cycle
       contrib = dabs(num - exact)
       accu_n += 1.d0
       contrib_relat = contrib/dabs(exact)
       if(contrib .gt. thr.or.contrib_relat .gt. thr )then
        print*,'naive'
        print*,'i,k,j,l',i,k,j,l
        print*,exact, num
        print*,contrib,contrib_relat
       endif

      accu2 += contrib
      accu1 += contrib_relat
     enddo
    enddo
   enddo
  enddo
  print*,'accu absolute ',accu2/accu_n
  print*,'accu relat    ',accu1/accu_n

end
