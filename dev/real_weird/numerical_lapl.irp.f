 BEGIN_PROVIDER [ double precision, mo_kin_plus, (3,mo_num,n_points_final_grid)]
 implicit none
 integer :: i,j,k
 double precision :: mos_array(mo_num),r(3),r_tmp(3)
 do i = 1, n_points_final_grid
  r(1) = final_grid_points(1,i)
  r(2) = final_grid_points(2,i)
  r(3) = final_grid_points(3,i)
  do k = 1, 3
   r_tmp = r
   r_tmp(k) += dx_deriv
   call give_all_mos_at_r(r_tmp,mos_array)
   do j = 1, mo_num
    mo_kin_plus(k,j,i) = mos_array(j)
   enddo
  enddo 
 enddo
 END_PROVIDER 

 BEGIN_PROVIDER [ double precision, mo_kin_minus, (3,mo_num,n_points_final_grid)]
 implicit none
 integer :: i,j,k
 double precision :: mos_array(mo_num),r(3),r_tmp(3)
 do i = 1, n_points_final_grid
  r(1) = final_grid_points(1,i)
  r(2) = final_grid_points(2,i)
  r(3) = final_grid_points(3,i)
  do k = 1, 3
   r_tmp = r
   r_tmp(k) -= dx_deriv
   call give_all_mos_at_r(r_tmp,mos_array)
   do j = 1, mo_num
    mo_kin_minus(k,j,i) = mos_array(j)
   enddo
  enddo 
 enddo
 END_PROVIDER 

 BEGIN_PROVIDER [ double precision, mo_num_lapl, (3,mo_num, n_points_final_grid)]
 implicit none
 integer :: i,j,k
 do i = 1, n_points_final_grid
  do j = 1, mo_num
   do k = 1, 3
    mo_num_lapl(k,j,i) = (mo_kin_minus(k,j,i) + mo_kin_plus(k,j,i) - 2.d0 * mos_in_r_array(j,i)) * inv_dx_deriv_2
   enddo
  enddo
 enddo
 END_PROVIDER 

 BEGIN_PROVIDER [ double precision, kin_num_ints, (mo_num, mo_num)]
 implicit none
 integer :: i,j,k,l
 kin_num_ints = 0.d0
 do i = 1, n_points_final_grid
  do j = 1, mo_num
   do l = 1, mo_num
    do k = 1, 3
     kin_num_ints(l,j) += mo_num_lapl(k,j,i) * mos_in_r_array(l,i) * final_weight_at_r_vector(i)
    enddo
   enddo
  enddo
 enddo
 kin_num_ints *= -0.5d0
 END_PROVIDER 

 BEGIN_PROVIDER [ double precision, kin_num_ints_bis, (mo_num, mo_num)]
 implicit none
 integer :: i,j,k,l
 kin_num_ints_bis = 0.d0
 do i = 1, n_points_final_grid
  do j = 1, mo_num
   do l = 1, mo_num
    do k = 1, 3
     kin_num_ints_bis(l,j) += mos_lapl_in_r_array_tranp(k,j,i) * mos_in_r_array(l,i) * final_weight_at_r_vector(i)
    enddo
   enddo
  enddo
 enddo
 kin_num_ints_bis *= -0.5d0
 END_PROVIDER 

subroutine test_num_lapl
 implicit none
 integer :: i,j
 double precision :: accu,accu2
 accu = 0.d0
 accu2 = 0.d0
 do i = 1, mo_num
  do j = 1, mo_num
!   accu2 += dabs(mo_kinetic_integrals(j,i) - kin_num_ints_bis(j,i))
!   print*,mo_kinetic_integrals(j,i),kin_num_ints_bis(j,i),dabs(mo_kinetic_integrals(j,i) - kin_num_ints_bis(j,i))
   if(dabs(kin_num_ints_bis(j,i) - kin_num_ints(j,i)).gt.1.d-3)then
    print*,kin_num_ints_bis(j,i),kin_num_ints(j,i),dabs(kin_num_ints_bis(j,i) - kin_num_ints(j,i))
   endif
   accu += dabs(kin_num_ints_bis(j,i) - kin_num_ints(j,i))
  enddo
 enddo
 print*,'accu = ',accu/dble(mo_num*mo_num)
! print*,'accu2= ',accu2/dble(mo_num * mo_num)

end
