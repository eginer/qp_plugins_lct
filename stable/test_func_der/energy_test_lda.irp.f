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


 BEGIN_PROVIDER [double precision, potential_x_alpha_mo_lda, (mo_num,mo_num)]
 implicit none
 BEGIN_DOC
 ! Representation of vx_alpha in the molecular orbitals basis 
 END_DOC
 call ao_to_mo(potential_x_alpha_ao_lda,ao_num,potential_x_alpha_mo_lda,mo_num)
 END_PROVIDER 

 BEGIN_PROVIDER [double precision, potential_x_beta_mo_lda, (mo_num,mo_num)]
 implicit none
 BEGIN_DOC
 ! Representation of vx_beta in the molecular orbitals basis 
 END_DOC
 call ao_to_mo(potential_x_beta_ao_lda,ao_num,potential_x_beta_mo_lda,mo_num)
 END_PROVIDER 
