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


 ! erf(mu(r) r12)/r12
 do ipoint = 1, n_points_final_grid
  do l = 1, ao_num
   do j = 1, ao_num
    do k = 1, ao_num
     do i = 1, ao_num
      big_array_naive(i,k,j,l) += ao_prod_on_grid_transp(i,k,ipoint) * erf_mu_r12_inv_r12_rk(j,l,ipoint)
     enddo
    enddo
   enddo
  enddo
 enddo

 !  mu(r) / sqpi * exp(-(mu(r) r12)^2)
 do ipoint = 1, n_points_final_grid
  if(constant_mu)then
   mu = mu_erf
  else
   mu = mu_of_r_for_ints(ipoint,1)
  endif
  do l = 1, ao_num
   do j = 1, ao_num
    do k = 1, ao_num
     do i = 1, ao_num
      big_array_naive(i,k,j,l) += ao_prod_on_grid_transp(i,k,ipoint) * gauss_ij_rk(j,l,ipoint) * mu * inv_sq_pi
     enddo
    enddo
   enddo
  enddo
 enddo

  ! - 1/4 (1 - erf(mu(r1)r12))^2
  do ipoint = 1, n_points_final_grid
   if(constant_mu)then
    mu = mu_erf
   else
    mu = mu_of_r_for_ints(ipoint,1)
   endif
   do l = 1, ao_num
    do j = 1, ao_num
     do k = 1, ao_num
      do i = 1, ao_num
       big_array_naive(i,k,j,l) += ao_prod_on_grid_transp(i,k,ipoint) * erf_mu_squared_ij_rk(j,l,ipoint) * (-0.25d0 )
      enddo
     enddo
    enddo
   enddo
  enddo


END_PROVIDER 

BEGIN_PROVIDER [ double precision, big_array_dgemm, (ao_num, ao_num, ao_num, ao_num)]
 implicit none
 include 'constants.include.F'
 integer :: i,k,ipoint
 double precision :: weight,mu
 double precision, allocatable :: a_mat(:,:,:)

 big_array_dgemm= 0.d0

 allocate(a_mat(n_points_final_grid,ao_num,ao_num))
 do i = 1, ao_num
  do k = 1, ao_num
   do ipoint = 1, n_points_final_grid
    weight = final_weight_at_r_vector(ipoint)
    a_mat(ipoint,k,i) = aos_in_r_array_transp(ipoint,k) * aos_in_r_array_transp(ipoint,i) * weight
   enddo
  enddo
 enddo

 ! erf(mu(r) * r12)/r12
 call dgemm("N","N",ao_num*ao_num,ao_num*ao_num,n_points_final_grid,1.d0,erf_mu_r12_inv_r12_rk(1,1,1),ao_num*ao_num & 
                   ,a_mat(1,1,1),n_points_final_grid,1.d0,big_array_dgemm,ao_num*ao_num)
!
! ! mu(r) / sqpi * exp[-(mu(r1)r12)^2]
 if(constant_mu)then
  do i = 1, ao_num
   do k = 1, ao_num
    do ipoint = 1, n_points_final_grid
     mu = mu_erf
     weight = final_weight_at_r_vector(ipoint)
     a_mat(ipoint,k,i) = aos_in_r_array_transp(ipoint,k) * aos_in_r_array_transp(ipoint,i) * weight * mu * inv_sq_pi
    enddo
   enddo
  enddo
 else
  do i = 1, ao_num
   do k = 1, ao_num
    do ipoint = 1, n_points_final_grid
     mu = mu_of_r_for_ints(ipoint,1)
     weight = final_weight_at_r_vector(ipoint)
     a_mat(ipoint,k,i) = aos_in_r_array_transp(ipoint,k) * aos_in_r_array_transp(ipoint,i) * weight * mu * inv_sq_pi
    enddo
   enddo
  enddo
 endif
 call dgemm("N","N",ao_num*ao_num,ao_num*ao_num,n_points_final_grid,1.d0,gauss_ij_rk(1,1,1),ao_num*ao_num & 
                   ,a_mat(1,1,1),n_points_final_grid,1.d0,big_array_dgemm,ao_num*ao_num)

  ! - 1/4 (1 - erf(mu(r1)r12))^2
 if(constant_mu)then
  do i = 1, ao_num
   do k = 1, ao_num
    do ipoint = 1, n_points_final_grid
     mu = mu_erf
     weight = final_weight_at_r_vector(ipoint)
     a_mat(ipoint,k,i) = aos_in_r_array_transp(ipoint,k) * aos_in_r_array_transp(ipoint,i) * weight * (-0.25d0)
    enddo
   enddo
  enddo
 else
  do i = 1, ao_num
   do k = 1, ao_num
    do ipoint = 1, n_points_final_grid
     mu = mu_of_r_for_ints(ipoint,1)
     weight = final_weight_at_r_vector(ipoint)
     a_mat(ipoint,k,i) = aos_in_r_array_transp(ipoint,k) * aos_in_r_array_transp(ipoint,i) * weight * (-0.25d0)
    enddo
   enddo
  enddo
 endif
 call dgemm("N","N",ao_num*ao_num,ao_num*ao_num,n_points_final_grid,1.d0,erf_mu_squared_ij_rk(1,1,1),ao_num*ao_num & 
                   ,a_mat(1,1,1),n_points_final_grid,1.d0,big_array_dgemm,ao_num*ao_num)

! call dgemm("N","N",ao_num*ao_num,ao_num*ao_num,n_points_final_grid,1.d0,erf_mu_r12_inv_r12_rk(1,1,1),ao_num*ao_num & 
!                     ,ao_prod_on_grid(1,1,1,m),n_points_final_grid,1.d0,ac_mat,ao_num*ao_num)
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

      contrib = dabs(big_array_dgemm(i,k,j,l) - exact)
      accu2 += contrib
      if(contrib .gt. 1.d-10)then
       print*,'dgemm'
       print*,'i,k,j,l',i,k,j,l
       print*,contrib,exact, big_array_dgemm(i,k,j,l)
      endif
     enddo
    enddo
   enddo
  enddo
  print*,'accu naive  = ',accu1/dble(ao_num**4)
  print*,'accu dgemm  = ',accu2/dble(ao_num**4)

end
