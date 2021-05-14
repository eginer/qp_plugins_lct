program test_stufff
 implicit none
 call test_erf_mu_gauss_xyz
end

subroutine test_erf_mu_gauss_xyz
 implicit none
 integer :: ipoint,i,j,m,jpoint
 double precision :: r1(3),derf_mu_x
 double precision :: C_center(3), weight1,mu,r12,integrals(4),integral_tmp(3,ao_num,ao_num),delta, num_int
 double precision :: integral,integral_scal(ao_num,ao_num)
 double precision, allocatable :: aos_array_r1(:),ao_mat_xyz(:,:,:),ao_mat(:,:)
 allocate(ao_mat_xyz(3,ao_num,ao_num),aos_array_r1(ao_num), ao_mat(ao_num, ao_num))
 do jpoint = 1, n_points_final_grid
  mu = mu_of_r_prov(jpoint,1)
  delta = mu * mu
  C_center(:) = final_grid_points(:,jpoint)
  do i = 1, ao_num
   do j = 1, ao_num
    call erf_mu_gauss_xyz_ij_ao(i,j,mu, C_center, delta,integrals)
    integral_scal(j,i) = integrals(4) 
    do m = 1, 3
     integral_tmp(m,j,i) = integrals(4) * C_center(m) - integrals(m) 
!     integral_tmp(m,j,i) = integrals(m) 
    enddo
   enddo
  enddo
  ao_mat_xyz = 0.d0
  ao_mat = 0.d0
  do ipoint = 1, n_points_final_grid
   r1(1) = final_grid_points(1,ipoint)
   r1(2) = final_grid_points(2,ipoint)
   r1(3) = final_grid_points(3,ipoint)
   call give_all_aos_at_r(r1,aos_array_r1)
   weight1 = final_weight_at_r_vector(ipoint)
   r12 = (r1(1) - C_center(1))**2.d0 + (r1(2) - C_center(2))**2.d0 + (r1(3) - C_center(3))**2.d0 
   r12 = dsqrt(r12)
   if(r12.lt.1.d-6)cycle
   do i = 1, ao_num
    do j = 1, ao_num
     ao_mat(j,i)  += aos_array_r1(i) * aos_array_r1(j) * weight1 * dexp(-mu*mu*r12*r12) * derf_mu_x(mu,r12)
     do m = 1, 3
      ao_mat_xyz(m,j,i)  += aos_array_r1(i) * aos_array_r1(j) * weight1 * dexp(-mu*mu*r12*r12) * (1.d0 / r12 - derf_mu_x(mu,r12)) * (C_center(m) - r1(m) )
!      ao_mat_xyz(m,j,i)  += aos_array_r1(i) * aos_array_r1(j) * weight1 * dexp(-mu*mu*r12*r12) * derf_mu_x(mu,r12) * r1(m) 
     enddo
    enddo
   enddo
  enddo


 double precision :: accu1relat,accu1
  accu1 = 0.d0
  accu1relat = 0.d0
  print*,'jpoint = ',jpoint
  do i = 1, ao_num
   do j = 1, ao_num

!    integral = integral_scal(j,i)
!    num_int  = ao_mat(j,i) 
!    if(dabs(num_int).gt.1.d-10)then
!     accu1relat = dabs(integral - num_int)/dabs(num_int)
!    endif
!    if(dabs(integral - num_int).gt.1.d-5)then
!     print*,'i,j,jpoint',i,j,jpoint
!     print*,'prov, num, difference'
!     print*,integral,num_int,dabs(integral - num_int)
!    endif
!    accu1 += dabs(integral - num_int)

    do m = 1, 3
     integral = integral_tmp(m,j,i)
     num_int  = ao_mat_xyz(m,j,i) 
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
    print*,'accu1      = ',accu1/dble(ao_num * ao_num)
    print*,''
    print*,'accu1relat = ',accu1relat/dble(ao_num * ao_num)
!    pause
   enddo
  enddo

 enddo


end
