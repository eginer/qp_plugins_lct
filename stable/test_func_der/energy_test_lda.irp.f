subroutine energy_x_lda_test (rho_a, rho_b,ex_lda_test)
 implicit none
 BEGIN_DOC
! exchange energy with the lda functional
 END_DOC
 integer :: istate,i,j
 double precision, intent(in) :: rho_a(n_points_final_grid), rho_b(n_points_final_grid)
 double precision, intent(out) :: ex_lda_test
 double precision :: weight
 double precision :: e_x,vx_a,vx_b
 
 ex_lda_test = 0.d0
  do i = 1, n_points_final_grid
   weight = final_weight_at_r_vector(i)
   call ex_lda(rho_a(i),rho_b(i),e_x,vx_a,vx_b)
   ex_lda_test += weight * e_x
  enddo
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

   delta_rho_a(i) = delta_rho_11 * 0.5d0 * mo_i(1) * mo_i(1)
   delta_rho_b(i) = delta_rho_11 * 0.5d0 * mo_i(1) * mo_i(1)

  enddo
end

BEGIN_PROVIDER [double precision, potential_x_alpha_mo_lda, (mo_num,mo_num)]
 implicit none
 call ao_to_mo(potential_x_alpha_ao_lda,ao_num,potential_x_alpha_mo_lda,mo_num)
END_PROVIDER 

BEGIN_PROVIDER [double precision, potential_x_beta_mo_lda, (mo_num,mo_num)]
 implicit none
 call ao_to_mo(potential_x_beta_ao_lda,ao_num,potential_x_beta_mo_lda,mo_num)
END_PROVIDER 
