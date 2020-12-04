
subroutine test_overlap_x_gauss_r12_ao
 implicit none
 integer :: i,j,m
 double precision :: D_center(3), delta
 double precision :: weight, r(3),accu,overlap_gauss_r12_ao,exact,aos_array(ao_num)
 double precision :: r12,int_gauss
 integer :: ipoint

 D_center = 0.d0
 D_center(1) =  0.d0
 D_center(3) =  0.0
 delta = 0.2d0
! do i = 1, ao_num
!  do j = 1, ao_num
 do i = 106,106
  do j = 106, 106
   do m = 3, 3
    call gauss_int_x_ao(i,j,delta,D_center,m,exact)
!    exact = overlap_gauss_r12_ao(D_center,delta,i,j)
    accu = 0.d0
    do ipoint = 1, n_points_final_grid
     r(1) = final_grid_points(1,ipoint)
     r(2) = final_grid_points(2,ipoint)
     r(3) = final_grid_points(3,ipoint)
     weight = final_weight_at_r_vector(ipoint)
     r12 = (r(1) - D_center(1))**2.d0 + (r(2) - D_center(2))**2.d0 + (r(3) - D_center(3))**2.d0
     call give_all_aos_at_r(r,aos_array)
     accu += weight * aos_array(i) * aos_array(j) * dexp(-delta * r12) * r(m)
!     accu += weight * aos_array(i) * aos_array(j) * dexp(-delta * r12)
    enddo
    print*,'i,j,m',i,j,m
    print*,accu,exact,dabs(accu - exact)
   enddo
  enddo
 enddo
end

subroutine test_overlap_x_gauss_ten_no
 implicit none
 integer :: i,j,m,pp,jpoint
 double precision :: D_center(3), delta
 double precision :: weight, r(3),accu(3),overlap_gauss_r12_ao,exact,aos_array(ao_num)
 double precision :: r12
 integer :: ipoint

 do pp = 1, n_max_fit_ten_no_slat
  do jpoint = 1, n_points_final_grid
   D_center(1) = final_grid_points(1,jpoint)
   D_center(2) = final_grid_points(2,jpoint)
   D_center(3) = final_grid_points(3,jpoint)
   delta = expo_fit_ten_no_slat_gauss(pp)
   do i = 1, ao_num
    do j = 1, ao_num
     accu = 0.d0
     do m = 1, 3
      exact = x_v_ij_gauss_rk(j,i,jpoint,m,pp)
      do ipoint = 1, n_points_final_grid
       r(1) = final_grid_points(1,ipoint)
       r(2) = final_grid_points(2,ipoint)
       r(3) = final_grid_points(3,ipoint)
       weight = final_weight_at_r_vector(ipoint)
       r12 = (r(1) - D_center(1))**2.d0 + (r(2) - D_center(2))**2.d0 + (r(3) - D_center(3))**2.d0
       call give_all_aos_at_r(r,aos_array)
       accu(m) += weight * aos_array(i) * aos_array(j) * dexp(-delta * r12) * r(m)
      enddo
      print*,'delta = ',delta
      print*,'D_center'
      print*, D_center 
      print*,'i,j,m',i,j,m
      print*,accu(m),exact,dabs(accu(m) - exact)
     enddo
    enddo
   enddo
  enddo
 enddo


end


subroutine test_deriv_ints
 implicit none
 integer :: i,j,k,l,m,pp,jpoint
 double precision :: alpha,integral,integral_x,coef
 double precision :: weight, r(3), accu,overlap_gauss_r12_ao,exact,aos_array(ao_num)
 double precision :: r12
 integer :: ipoint
 double precision, allocatable :: array_tmp(:,:,:,:)
 allocate(array_tmp(ao_num, ao_num, ao_num, ao_num))
 array_tmp = 0.d0
 do ipoint = 1, n_points_final_grid
  r(1) = final_grid_points(1,ipoint)
  r(2) = final_grid_points(2,ipoint)
  r(3) = final_grid_points(3,ipoint)
  weight = final_weight_at_r_vector(ipoint)
  do m = 1, 3
   do i = 1, ao_num
    do k = 1, ao_num
     do pp = 1, n_max_fit_ten_no_slat
      alpha = expo_fit_ten_no_slat_gauss(pp)
      coef = coef_fit_ten_no_slat_gauss(pp)
!      call gauss_int_x_ao(i,k,alpha,r,m,integral_x)
      integral = overlap_gauss_r12_ao(r,alpha,i,k)
      do j = 1, ao_num
       do l = 1, ao_num 
        ! -2 alpha * x_1 * phi_k * d/dx phi_i * w_jl(r)
        array_tmp(k,i,l,j) += -2.d0 * weight * coef * alpha * integral * aos_in_r_array_transp(ipoint,l) * r(m) * aos_grad_in_r_array(j,ipoint,m) * integral
!        ! +2 alpha phi_k * d/dx phi_i * w_jl^x(r)
!        array_tmp(l,j,k,i) +=  2.d0 * weight * coef * alpha * integral * aos_in_r_array_transp(ipoint,l)        * aos_grad_in_r_array(j,ipoint,m) * integral_x
       enddo
      enddo
     enddo
    enddo 
   enddo
  enddo
 enddo

! do m = 1, 3
!  do i = 1, ao_num
!    do k = 1, ao_num
!      do ipoint = 1, n_points_final_grid
!       r(1) = final_grid_points(1,ipoint)
!       r(2) = final_grid_points(2,ipoint)
!       r(3) = final_grid_points(3,ipoint)
!       weight = final_weight_at_r_vector(ipoint)
!       do pp = 1, n_max_fit_ten_no_slat
!        alpha = expo_fit_ten_no_slat_gauss(pp)
!        coef = coef_fit_ten_no_slat_gauss(pp)
!        call gauss_int_x_ao(i,k,alpha,r,m,integral_x)
!        integral = overlap_gauss_r12_ao(r,alpha,i,k)
!        do j = 1, ao_num
!         do l = 1, ao_num 
!          ! -2 alpha * x_1 * phi_l * d/dx phi_j * w_ki(r)
!          array_tmp(l,j,k,i) += -2.d0 * weight * coef * alpha * integral * aos_in_r_array_transp(ipoint,l) * r(m) * aos_grad_in_r_array(j,ipoint,m) * integral
!          ! +2 alpha phi_l * d/dx phi_j * w_ik^x(r)
!          array_tmp(l,j,k,i) +=  2.d0 * weight * coef * alpha * integral * aos_in_r_array_transp(ipoint,l)        * aos_grad_in_r_array(j,ipoint,m) * integral_x
!         enddo
!        enddo
!       enddo
!      enddo 
!    enddo
!   enddo
! enddo

 do j = 1, ao_num
  do l = 1, ao_num
   do i = 1, ao_num
    do k = 1, ao_num
     accu = dabs(array_tmp(k,i,l,j) - ao_two_e_eff_dr12_pot_array_new_bis(k,i,l,j))
     print*,k,i,l,j
     print*,'array_tmp, ao_two_e_eff_dr12_pot_array_new_bis, accu'
     print*, array_tmp(k,i,l,j), ao_two_e_eff_dr12_pot_array_new_bis(k,i,l,j), accu 
    enddo
   enddo
  enddo
 enddo
end
