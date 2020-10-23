subroutine delta_gamma_i_j_for_energy_test (delta_rho_11,delta_gamma_i_j)
 implicit none
 BEGIN_DOC
! Considering the small variation of density : delta_n(r) = sum_i,j delta_gamma_i,j * mo_i(r) * mo_j(r),
! we choose delta_gamma_i,j = epsilon * kronecker(i,1) * kronecker(1,j)
! so that delta_n(r) = delta_gamma_11 * mo_1(r) * mo_1(r)
!
! To distinguish the different epsilon we'll need, we note epsilon -> delta_rho_11
 END_DOC
 double precision, intent(in)  :: delta_rho_11
 double precision, intent(out) :: delta_gamma_i_j(mo_num,mo_num)
 integer :: i,j
   do i=1, mo_num
    do j=1, mo_num
     delta_gamma_i_j(i,j) = 0.d0
    enddo
   enddo
   delta_gamma_i_j(1,1) = delta_rho_11
end

subroutine delta_gamma_i_j_for_energy_test_general (delta_rho_11,delta_gamma_i_j)
 implicit none
 BEGIN_DOC
! Considering the small variation of density : delta_n(r) = sum_i,j delta_gamma_i,j * mo_i(r) * mo_j(r),
! we choose delta_gamma_i,j = epsilon/(mo_num*mo_num)
!
! To distinguish the different epsilon we'll need, we note epsilon -> delta_rho_11
 END_DOC
 double precision, intent(in)  :: delta_rho_11
 double precision, intent(out) :: delta_gamma_i_j(mo_num,mo_num)
 integer :: i,j
   do i=1, mo_num
    do j=1, mo_num
     delta_gamma_i_j(i,j) = delta_rho_11 / (mo_num)
    enddo
   enddo
end

subroutine delta_gamma_i_j_k_l_for_energy_test (delta_rho2_1111,delta_gamma_i_j_k_l)
 implicit none
 BEGIN_DOC
! Considering the small variation of on-top density : delta_n2(r) = sum_i,j,k,l delta_gamma_i,j,k,l * mo_i(r) * mo_j(r) * mo_k(r) * mo_l(r),
! we choose delta_gamma_i,j,k,l = epsilon if 1=i=j=k=l
! so that delta_n2(r) = delta_gamma_1111 * mo_1(r) * mo_1(r) * mo_1(r) * mo_1(r)
!
! To distinguish the different epsilon we'll need, we note epsilon -> delta_rho2_1111
 END_DOC
 double precision, intent(in)  :: delta_rho2_1111
 double precision, intent(out) :: delta_gamma_i_j_k_l(mo_num,mo_num,mo_num,mo_num)
 integer :: j,k,l,m
   do j=1, mo_num
    do k=1, mo_num
     do l=1, mo_num
      do m=1, mo_num
       delta_gamma_i_j_k_l(j,k,l,m) = 0.d0
      enddo
     enddo
    enddo
   enddo
   delta_gamma_i_j_k_l(1,1,1,1) = delta_rho2_1111
end

subroutine delta_gamma_i_j_k_l_for_energy_test_general (delta_rho2_1111,delta_gamma_i_j_k_l)
 implicit none
 BEGIN_DOC
! Considering the small variation of on-top density : delta_n2(r) = sum_i,j,k,l delta_gamma_i,j,k,l * mo_i(r) * mo_j(r) * mo_k(r) * mo_l(r),
! we choose delta_gamma_i,j,k,l = epsilon / (mo_num**4)
!
! To distinguish the different epsilon we'll need, we note epsilon -> delta_rho2_1111
 END_DOC
 double precision, intent(in)  :: delta_rho2_1111
 double precision, intent(out) :: delta_gamma_i_j_k_l(mo_num,mo_num,mo_num,mo_num)
 integer :: j,k,l,m
   do j=1, mo_num
    do k=1, mo_num
     do l=1, mo_num
      do m=1, mo_num
       delta_gamma_i_j_k_l(j,k,l,m) = delta_rho2_1111 / (mo_num**2)
      enddo
     enddo
    enddo
   enddo
end

subroutine delta_gamma_i_j_for_energy_test_alpha_beta (delta_rho_11_alpha,delta_rho_11_beta,delta_gamma_i_j_alpha, delta_gamma_i_j_beta)
 implicit none
 BEGIN_DOC
! Same as the subroutine delta_gamma_i_j_for_energy_test : we distinguish the contribution of alpha and beta spins.
 END_DOC
 double precision, intent(in)  :: delta_rho_11_alpha, delta_rho_11_beta
 double precision, intent(out) :: delta_gamma_i_j_alpha(mo_num,mo_num) ,delta_gamma_i_j_beta(mo_num,mo_num)
 integer :: i,j

   do i=1, mo_num
    do j=1, mo_num
     delta_gamma_i_j_alpha(i,j) = 0.d0
     delta_gamma_i_j_beta(i,j) = 0.d0
    enddo
   enddo
   delta_gamma_i_j_alpha(1,1) = delta_rho_11_alpha
   delta_gamma_i_j_beta(1,1) = delta_rho_11_beta
end

subroutine delta_gamma_i_j_for_energy_test_alpha_beta_general (delta_rho_11_alpha,delta_rho_11_beta,delta_gamma_i_j_alpha, delta_gamma_i_j_beta)
 implicit none
 BEGIN_DOC
! Same as the subroutine delta_gamma_i_j_for_energy_test_general : we distinguish the contribution of alpha and beta spins.
 END_DOC
 double precision, intent(in)  :: delta_rho_11_alpha, delta_rho_11_beta
 double precision, intent(out) :: delta_gamma_i_j_alpha(mo_num,mo_num) ,delta_gamma_i_j_beta(mo_num,mo_num)
 integer :: i,j

   do i=1, mo_num
    do j=1, mo_num
     delta_gamma_i_j_alpha(i,j) = delta_rho_11_alpha / (mo_num)
     delta_gamma_i_j_beta(i,j) = delta_rho_11_beta / (mo_num)
    enddo
   enddo
end


subroutine delta_density_for_energy_test (delta_rho_a, delta_rho_b, delta_rho_11_alpha, delta_rho_11_beta)
 implicit none
 BEGIN_DOC
! delta_n(r) = sum_i,j delta_gamma_i,j * mo_i(r) * mo_j(r)
! where delta_gamma_i,j = epsilon * kronecker(i,1) * kronecker(1,j) (from delta_gamma_i_j_for_energy_test)
! then  delta_n(r) = epsilon * mo_i(r) * mo_i(r)
 END_DOC
 double precision, intent(in)  :: delta_rho_11_alpha, delta_rho_11_beta
 double precision, intent(out) :: delta_rho_a(n_points_final_grid), delta_rho_b(n_points_final_grid)
 integer :: i,j
 double precision :: r(3), mo_i(mo_num)
  do i = 1, n_points_final_grid
   r(1) = final_grid_points(1,i)   
   r(2) = final_grid_points(2,i) 
   r(3) = final_grid_points(3,i)
   call give_all_mos_at_r(r, mo_i)
   delta_rho_a(i) = delta_rho_11_alpha  * mo_i(1) * mo_i(1)
   delta_rho_b(i) = delta_rho_11_beta  * mo_i(1) * mo_i(1)  
  enddo
end

subroutine delta_density_for_energy_test_general (delta_rho_a, delta_rho_b, delta_rho_11_alpha, delta_rho_11_beta)
 implicit none
 BEGIN_DOC
! delta_n(r) = sum_i,j delta_gamma_i,j * mo_i(r) * mo_j(r)
! where delta_gamma_i,j = epsilon' * mo_i(r) * mo_j(r)
! with epsilon' = epsilon / (mo_num**2)
 END_DOC
 double precision, intent(in)  :: delta_rho_11_alpha, delta_rho_11_beta
 double precision, intent(out) :: delta_rho_a(n_points_final_grid), delta_rho_b(n_points_final_grid)
 integer :: i,j,k
 double precision :: r(3), mo_i(mo_num)
  do i = 1, n_points_final_grid
   r(1) = final_grid_points(1,i)   
   r(2) = final_grid_points(2,i) 
   r(3) = final_grid_points(3,i)
   call give_all_mos_at_r(r, mo_i)
   delta_rho_a(i)=0.d0
   delta_rho_b(i)=0.d0
   do j=1, mo_num
    do k=1, mo_num
     delta_rho_a(i) += (delta_rho_11_alpha/(mo_num))  * mo_i(j) * mo_i(k)
     delta_rho_b(i) += (delta_rho_11_beta /(mo_num))  * mo_i(j) * mo_i(k)  
    enddo
   enddo
  enddo
end

subroutine delta_grad_density_for_energy_test (delta_grad_rho_a, delta_grad_rho_b, delta_grad_rho_11_alpha, delta_grad_rho_11_beta)
 implicit none
 BEGIN_DOC
 ! Small variation of the gradient of the density
 ! delta_grad_n(r) = sum_i,j delta_gamma_i,j * (grad_mo_i(r) * mo_j(r) + mo_i(r) * grad_mo_j(r))
 ! where delta_gamma_i,j = epsilon * kronecker(i,1) * kronecker(1,j) (from delta_gamma_i_j_for_energy_test)
 ! then  delta_grad_n(r) = epsilon * (grad_mo_i(r) * mo_i(r) + mo_i(r) * grad_mo_i(r))
 END_DOC
 double precision, intent(in)  :: delta_grad_rho_11_alpha, delta_grad_rho_11_beta
 double precision, intent(out) :: delta_grad_rho_a(3,n_points_final_grid), delta_grad_rho_b(3,n_points_final_grid)
 integer :: i,m
 double precision :: r(3), mo_i(mo_num), mo_grad_i(3,mo_num)
  do i = 1, n_points_final_grid
   r(1) = final_grid_points(1,i)   
   r(2) = final_grid_points(2,i)   
   r(3) = final_grid_points(3,i)   
   call give_all_mos_and_grad_at_r(r,mo_i,mo_grad_i)
   do m=1,3
     delta_grad_rho_a(m,i) = 2.d0 * delta_grad_rho_11_alpha * mo_i(1) * mo_grad_i(m,1) 
     delta_grad_rho_b(m,i) = 2.d0 * delta_grad_rho_11_beta * mo_i(1) * mo_grad_i(m,1)
   enddo
  enddo

end

subroutine delta_grad_density_for_energy_test_general (delta_grad_rho_a, delta_grad_rho_b, delta_grad_rho_11_alpha, delta_grad_rho_11_beta)
 implicit none
 BEGIN_DOC
 ! Small variation of the gradient of the density
 ! delta_grad_n(r) = sum_i,j delta_gamma_i,j * (grad_mo_i(r) * mo_j(r) + mo_i(r) * grad_mo_j(r))
 ! where delta_gamma_i,j = epsilon' (from delta_gamma_i_j_for_energy_test_general)
 ! with epsilon' = epsilon / (mo_num**2)
 END_DOC
 double precision, intent(in)  :: delta_grad_rho_11_alpha, delta_grad_rho_11_beta
 double precision, intent(out) :: delta_grad_rho_a(3,n_points_final_grid), delta_grad_rho_b(3,n_points_final_grid)
 integer :: i,m,j,k
 double precision :: r(3), mo_i(mo_num), mo_grad_i(3,mo_num)
  do i = 1, n_points_final_grid
   r(1) = final_grid_points(1,i)   
   r(2) = final_grid_points(2,i)   
   r(3) = final_grid_points(3,i)   
   call give_all_mos_and_grad_at_r(r,mo_i,mo_grad_i)
   do m=1,3
    delta_grad_rho_a(m,i) = 0.d0
    delta_grad_rho_b(m,i) = 0.d0
    do j=1, mo_num
     do k=1, mo_num
      delta_grad_rho_a(m,i) += delta_grad_rho_11_alpha/(mo_num) * (mo_i(j) * mo_grad_i(m,k) + mo_i(k) * mo_grad_i(m,j))
      delta_grad_rho_b(m,i) += delta_grad_rho_11_beta/(mo_num) *  (mo_i(j) * mo_grad_i(m,k) + mo_i(k) * mo_grad_i(m,j))
     enddo
    enddo
   enddo
  enddo

end

subroutine delta_n2_for_energy_test (delta_n2, delta_n2_11) 
 implicit none
 BEGIN_DOC
! Small variation of the on-top pair density
! delta_n2 (r) = sum_i,j,k,l delta_gamma_i,j,k,l * mo_i(r) * mo_j(r) * mo_k(r) * mo_l(r) 
! where delta_gamma_i,j,k,l is non-zero if i=j=k=l=1
! then delta_n2(r) = epsilon * mo_i(r) * mo_i(r) * mo_i(r) * mo_i(r)
 END_DOC
 double precision, intent(in)  :: delta_n2_11
 double precision, intent(out) :: delta_n2(n_points_final_grid) 
 integer :: i,j
 double precision :: r(3), mo_i(mo_num)
  do i = 1, n_points_final_grid
   r(1) = final_grid_points(1,i)   
   r(2) = final_grid_points(2,i) 
   r(3) = final_grid_points(3,i)
   call give_all_mos_at_r(r, mo_i)
   delta_n2(i) = delta_n2_11  * mo_i(1) * mo_i(1) * mo_i(1) * mo_i(1)
  enddo
end

subroutine delta_n2_for_energy_test_general (delta_n2, delta_n2_11) 
 implicit none
 BEGIN_DOC
! Small variation of the on-top pair density
! delta_n2 (r) = sum_i,j,k,l delta_gamma_i,j,k,l * mo_i(r) * mo_j(r) * mo_k(r) * mo_l(r) 
! where delta_gamma_i,j,k,l = epsilon / (mo_num**4)
 END_DOC
 double precision, intent(in)  :: delta_n2_11
 double precision, intent(out) :: delta_n2(n_points_final_grid) 
 integer :: i,j,k,l,m
 double precision :: r(3), mo_i(mo_num)
  do i = 1, n_points_final_grid
   r(1) = final_grid_points(1,i)   
   r(2) = final_grid_points(2,i) 
   r(3) = final_grid_points(3,i)
   call give_all_mos_at_r(r, mo_i)
   delta_n2(i) = 0.d0
   do j = 1, mo_num
    do k = 1, mo_num
     do l = 1, mo_num
      do m = 1, mo_num
       delta_n2(i) += (delta_n2_11/(mo_num**2)) * mo_i(j) * mo_i(k) * mo_i(l) * mo_i(m)
      enddo
     enddo
    enddo
   enddo
  enddo
end

