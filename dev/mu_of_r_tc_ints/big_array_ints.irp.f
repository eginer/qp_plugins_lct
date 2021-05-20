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


subroutine test_big_array
 implicit none
 integer :: i,j,k,l
 double precision :: accu1, accu2,contrib,exact,ao_two_e_integral_erf,ao_two_e_integral_eff_pot
 accu1 = 0.d0
 accu2 = 0.d0
  do l = 1, ao_num
   do j = 1, ao_num
    do k = 1, ao_num
     do i = 1, ao_num
      exact = ao_two_e_integral_erf(i,k,j,l) 
      exact += ao_two_e_integral_eff_pot(i,k,j,l)

      contrib = dabs(exact - big_array_naive(i,k,j,l))
      accu1 += contrib
      if(contrib .gt. 1.d-10)then
       print*,'naive'
       print*,'i,k,j,l',i,k,j,l
       print*,contrib,exact, big_array_naive(i,k,j,l)
      endif

      contrib = dabs(scalar_mu_r_pot_chemist_ao(i,k,j,l) - exact)
      accu2 += contrib
      if(contrib .gt. 1.d-10)then
       print*,'dgemm'
       print*,'i,k,j,l',i,k,j,l
       print*,contrib,exact, scalar_mu_r_pot_chemist_ao(i,k,j,l)
      endif
     enddo
    enddo
   enddo
  enddo
  print*,'accu naive  = ',accu1/dble(ao_num**4)
  print*,'accu dgemm  = ',accu2/dble(ao_num**4)

end

subroutine test_int_r6
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
!        accu(i,k,j,l) += erf_mur1(ipoint,jpoint) * ao_prod_r1 * ao_prod_r2
!        accu(i,k,j,l) += gauss_r12_mu_r1(ipoint,jpoint,cst_gauss_r12) * ao_prod_r1 * ao_prod_r2
!        accu(i,k,j,l) += erf_mu_sq(ipoint,jpoint) * ao_prod_r1 * ao_prod_r2
!        accu(i,k,j,l) += nabla_sq_term(ipoint,jpoint,cst_nabla) * ao_prod_r1 * ao_prod_r2
!        accu(i,k,j,l) += nabla_r12_1(ipoint,jpoint,cst_nabla_r12_1) * ao_prod_r1 * ao_prod_r2
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
    do i = j, j
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

 print*,'accu naive  = ',accu_naive/dble(ao_num**1)
! print*,'accu dgemm  = ',accu2/dble(ao_num**4)


end

double precision function nabla_sq_term(ipoint,jpoint,cst)
 implicit none
 integer, intent(in) :: ipoint ! r1 
 integer, intent(in) :: jpoint ! r1 
 double precision, intent(in) :: cst
 include 'constants.include.F'
 double precision :: r12,r1(3),r2(3),mu
 r1(:) = final_grid_points(:,ipoint)

 r2(:) = final_grid_points(:,jpoint)

 r12  = (r1(1) - r2(1))*(r1(1) - r2(1))
 r12 += (r1(2) - r2(2))*(r1(2) - r2(2))
 r12 += (r1(3) - r2(3))*(r1(3) - r2(3))
 mu = mu_of_r_for_ints(ipoint,1)
 nabla_sq_term = cst * inv_4_mu_of_r_for_ints(ipoint,1) * dexp(-2.d0 * mu * mu * r12) * grad_sq_mu_of_r_for_ints(ipoint,1)

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
 integer, intent(in) :: ipoint
 integer, intent(in) :: jpoint
 double precision, intent(in) :: cst
 include 'constants.include.F'
 double precision :: r12,r12_sq,r12_vec(3),r1(3),r2(3),mu

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

double precision function erf_mu_sq(ipoint,jpoint)
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
 erf_mu_sq = -0.25d0 * (1.d0 - derf(mu*r12))**2.d0

end
