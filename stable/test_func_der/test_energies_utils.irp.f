subroutine delta_gamma_i_j_for_energy_test (delta_rho_11,delta_gamma_i_j)
 implicit none
 BEGIN_DOC
!
 END_DOC

 double precision, intent(in)  :: delta_rho_11
 double precision, intent(out) :: delta_gamma_i_j(mo_num,mo_num)
 
 integer :: istate,i,j
 double precision :: r(3), mo_i(mo_num)


!  do i = 1, n_points_final_grid
!   r(1) = final_grid_points(1,i)   
!   r(2) = final_grid_points(2,i) 
!   r(3) = final_grid_points(3,i)
!   call give_all_mos_at_r(r, mo_i)
   do i=1, mo_num
    do j=1, mo_num
     delta_gamma_i_j(i,j) = 0.d0
    enddo
   enddo
   delta_gamma_i_j(1,1) = delta_rho_11
end

subroutine delta_gamma_i_j_k_l_for_energy_test (delta_rho_11,delta_gamma_i_j_k_l)
 implicit none
 BEGIN_DOC
!
 END_DOC

 double precision, intent(in)  :: delta_rho_11
 double precision, intent(out) :: delta_gamma_i_j_k_l(mo_num,mo_num,mo_num,mo_num)
 
 integer :: istate,j,k,l,m
 double precision :: r(3), mo_i(mo_num)


!  do i = 1, n_points_final_grid
!   r(1) = final_grid_points(1,i)   
!   r(2) = final_grid_points(2,i) 
!   r(3) = final_grid_points(3,i)
!   call give_all_mos_at_r(r, mo_i)
   do j=1, mo_num
    do k=1, mo_num
     do l=1, mo_num
      do m=1, mo_num
       delta_gamma_i_j_k_l(j,k,l,m) = 0.d0
      enddo
     enddo
    enddo
   enddo
   delta_gamma_i_j_k_l(1,1,1,1) = delta_rho_11
end

subroutine delta_density_for_energy_test (delta_rho_a, delta_rho_b, delta_rho_11)
 implicit none
 BEGIN_DOC
!
 END_DOC

 double precision, intent(in)  :: delta_rho_11
 double precision, intent(out) :: delta_rho_a(n_points_final_grid), delta_rho_b(n_points_final_grid)
 
 integer :: istate,i,j
 double precision :: r(3), mo_i(mo_num)


  do i = 1, n_points_final_grid
   r(1) = final_grid_points(1,i)   
   r(2) = final_grid_points(2,i) 
   r(3) = final_grid_points(3,i)
   call give_all_mos_at_r(r, mo_i)

   delta_rho_a(i) = delta_rho_11  * mo_i(1) * mo_i(1)   !delta_rho_ij = delta_rho_11 = kro(i1)kro(j1)*epsilon
   delta_rho_b(i) = delta_rho_11  * mo_i(1) * mo_i(1)  

  enddo
end

subroutine delta_grad_density_for_energy_test (delta_grad_rho_a, delta_grad_rho_b, delta_grad_rho_11)
 implicit none

 double precision, intent(in)  :: delta_grad_rho_11
 double precision, intent(out) :: delta_grad_rho_a(3,n_points_final_grid), delta_grad_rho_b(3,n_points_final_grid)
 integer :: i,m
 double precision :: r(3), mo_i(mo_num), mo_grad_i(mo_num,3)


  do i = 1, n_points_final_grid
   r(1) = final_grid_points(1,i)   
   r(2) = final_grid_points(2,i)   
   r(3) = final_grid_points(3,i)   
 
   call give_all_mos_and_grad_at_r(r,mo_i,mo_grad_i)

   do m=1,3
     delta_grad_rho_a(m,i) = 2.d0 * delta_grad_rho_11 * mo_i(1) * mo_grad_i(1,m) 
     delta_grad_rho_b(m,i) = 2.d0 * delta_grad_rho_11 * mo_i(1) * mo_grad_i(1,m)
   enddo

  enddo

end

