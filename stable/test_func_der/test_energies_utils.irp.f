 BEGIN_PROVIDER [ double precision, delta_gamma_i_j_alpha, (mo_num, mo_num)]
&BEGIN_PROVIDER [ double precision, delta_gamma_i_j_beta, (mo_num, mo_num)]
 implicit none
 integer :: i,j
   do i=1, mo_num
    do j=1, mo_num
     delta_gamma_i_j_alpha(i,j) = delta_rho_11_alpha / dble(mo_num)
     delta_gamma_i_j_beta(i,j) = delta_rho_11_beta / dble(mo_num)
    enddo
   enddo
END_PROVIDER 

 BEGIN_PROVIDER [ double precision, delta_rho_11_alpha ]
&BEGIN_PROVIDER [ double precision, delta_rho_11_beta ]
 implicit none
 delta_rho_11_alpha = 1.d-5
 delta_rho_11_beta  = 1.d-5
END_PROVIDER 


subroutine delta_density_for_energy_test_general (delta_rho_a, delta_rho_b)
 implicit none
 BEGIN_DOC
! delta_n(r) = sum_i,j delta_gamma_i,j * mo_i(r) * mo_j(r)
! where delta_gamma_i,j = epsilon' * mo_i(r) * mo_j(r)
! with epsilon' = epsilon / (mo_num**2)
 END_DOC
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
     delta_rho_a(i) += delta_gamma_i_j_alpha(k,j)  * mo_i(j) * mo_i(k)
     delta_rho_b(i) += delta_gamma_i_j_beta(k,j)   * mo_i(j) * mo_i(k)  
    enddo
   enddo
  enddo
end

subroutine delta_grad_density_for_energy_test_general (delta_grad_rho_a, delta_grad_rho_b)
 implicit none
 BEGIN_DOC
 ! Small variation of the gradient of the density
 ! delta_grad_n(r) = sum_i,j delta_gamma_i,j * (grad_mo_i(r) * mo_j(r) + mo_i(r) * grad_mo_j(r))
 ! where delta_gamma_i,j = epsilon' (from delta_gamma_i_j_for_energy_test_general)
 ! with epsilon' = epsilon / (mo_num**2)
 END_DOC
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
      delta_grad_rho_a(m,i) += delta_gamma_i_j_alpha(k,j)   * (mo_i(j) * mo_grad_i(m,k) + mo_i(k) * mo_grad_i(m,j))
      delta_grad_rho_b(m,i) += delta_gamma_i_j_beta(k,j)    * (mo_i(j) * mo_grad_i(m,k) + mo_i(k) * mo_grad_i(m,j))
     enddo
    enddo
   enddo
  enddo

end

subroutine delta_n2_for_energy_test_general (delta_n2) 
 implicit none
 include 'constants.include.F'
 BEGIN_DOC
! Small variation of the on-top pair density
! delta_n2 (r) = sum_i,j,k,l delta_gamma_i,j,k,l * mo_i(r) * mo_j(r) * mo_k(r) * mo_l(r) 
! where delta_gamma_i,j,k,l = epsilon / (mo_num**4)
 END_DOC
 double precision, intent(out) :: delta_n2(n_points_final_grid) 
 integer :: i,j,k,l,m
 double precision :: r(3), mo_i(mo_num)
  do i = 1, n_points_final_grid
   r(1) = final_grid_points(1,i)   
   r(2) = final_grid_points(2,i) 
   r(3) = final_grid_points(3,i)
   call give_all_mos_at_r(r, mo_i)
!  mos_in_r_array(m,i)
   delta_n2(i) = 0.d0
   do j = 1, mo_num
    do k = 1, mo_num
     do l = 1, mo_num
      do m = 1, mo_num
       delta_n2(i) += delta_n2_ijkl(m,l,k,j) * mo_i(j) * mo_i(k) * mo_i(l) * mo_i(m)! * 2.d0 / ((1.d0 + 2.d0/(sqpi*mu))**2)
      enddo
     enddo
    enddo
   enddo
  enddo
end

BEGIN_PROVIDER [double precision, delta_n2_ijkl, (mo_num, mo_num, mo_num, mo_num)]
 implicit none
 integer :: i,j,k,l
   do j = 1, mo_num
    do k = 1, mo_num
     do l = 1, mo_num
      do i = 1, mo_num
       delta_n2_ijkl(i,l,k,j) = (delta_n2_11/(mo_num**2)) 
      enddo
     enddo
    enddo
   enddo
END_PROVIDER 

BEGIN_PROVIDER [double precision, delta_n2_11]
 implicit none
 delta_n2_11 = 1.d-5
END_PROVIDER 


 BEGIN_PROVIDER [ double precision, delta_rho_a_in_r, (n_points_final_grid)]
&BEGIN_PROVIDER [ double precision, delta_rho_b_in_r, (n_points_final_grid)]
 implicit none
 call delta_density_for_energy_test_general(delta_rho_a_in_r, delta_rho_b_in_r)
END_PROVIDER 

 BEGIN_PROVIDER [ double precision, dm_a_plus_delta_rho_a_in_r, (n_points_final_grid)]
&BEGIN_PROVIDER [ double precision, dm_b_plus_delta_rho_b_in_r, (n_points_final_grid)]
 implicit none
 integer :: i
 do i = 1, n_points_final_grid
  dm_a_plus_delta_rho_a_in_r(i) = one_e_dm_and_grad_alpha_in_r(4,i,1) + delta_rho_a_in_r(i)
  dm_b_plus_delta_rho_b_in_r(i) = one_e_dm_and_grad_beta_in_r(4,i,1)  + delta_rho_b_in_r(i)
 enddo
END_PROVIDER 


subroutine compute_func_der(delta_dm_a, delta_dm_b, pot_a, pot_b, delta_e)
 implicit none 
 double precision, intent(in) :: delta_dm_a(mo_num, mo_num), delta_dm_b(mo_num, mo_num)
 double precision, intent(in) :: pot_a(mo_num, mo_num), pot_b(mo_num, mo_num)
 double precision, intent(out):: delta_e
 integer :: i,j
 delta_e = 0.d0
 do i = 1, mo_num
  do j = 1, mo_num
   delta_e += delta_dm_a(j,i) * pot_a(j,i) + delta_dm_b(j,i) * pot_b(j,i)
  enddo 
 enddo
end

subroutine compute_func_der_two_e(delta_gamma_ijkl, delta_ijkl, delta_e)
 implicit none
 double precision, intent(in) :: delta_gamma_ijkl(mo_num,mo_num,mo_num,mo_num)
 double precision, intent(in) :: delta_ijkl(mo_num,mo_num,mo_num,mo_num)
 double precision, intent(out):: delta_e
 delta_e = 0.d0
 integer :: i,j,k,l
 do i = 1, mo_num
  do j = 1, mo_num
   do k = 1, mo_num
    do l = 1, mo_num
     delta_e += delta_gamma_ijkl(l,k,j,i) * delta_ijkl(l,k,j,i)
    enddo
   enddo
  enddo
 enddo
end

 BEGIN_PROVIDER [ double precision, delta_grad_rho_a_in_r, (3,n_points_final_grid)]
&BEGIN_PROVIDER [ double precision, delta_grad_rho_b_in_r, (3,n_points_final_grid)]
 implicit none
 call delta_grad_density_for_energy_test_general (delta_grad_rho_a_in_r, delta_grad_rho_b_in_r)
END_PROVIDER 

 BEGIN_PROVIDER [ double precision, grad_plus_delta_grad_rho_a_in_r, (3,n_points_final_grid)]
&BEGIN_PROVIDER [ double precision, grad_plus_delta_grad_rho_b_in_r, (3,n_points_final_grid)]
 implicit none
 integer :: i
 do i = 1, n_points_final_grid
  grad_plus_delta_grad_rho_a_in_r(1:3,i) = delta_grad_rho_a_in_r(1:3,i) + one_e_dm_and_grad_alpha_in_r(1:3,i,1)
  grad_plus_delta_grad_rho_b_in_r(1:3,i) = delta_grad_rho_b_in_r(1:3,i) + one_e_dm_and_grad_beta_in_r(1:3,i,1)
 enddo
END_PROVIDER 

 BEGIN_PROVIDER [ double precision, n2_plus_dn2, (n_points_final_grid)
 implicit none
 integer :: i
 call delta_n2_for_energy_test_general (n2_plus_dn2) 
 do i = 1, n_points_final_grid
  n2_plus_dn2(i) += total_cas_on_top_density(i,1)
 enddo

END_PROVIDER 
