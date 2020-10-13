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

subroutine int_potential_x_lda_test (delta_rho_11,rho_a,rho_b,int_vx_lda_test)
 implicit none
 BEGIN_DOC
! delta_Ex[n] = int_dr delta_n(r) dEx/dn(r) 
! delta_n(r) <- subroutine delta_density_for_energy_test
! d_e_d_n    <- subroutine de_dn_at_r
 END_DOC
 integer :: i
 double precision, intent(in) :: delta_rho_11
 double precision, intent(out):: int_vx_lda_test
 double precision, intent(in) :: rho_a(n_points_final_grid), rho_b(n_points_final_grid)
 double precision :: vx_i_j
 double precision :: delta_gamma_i_j(mo_num, mo_num)
 integer :: k,l
 
 int_vx_lda_test = 0.d0
 call delta_gamma_i_j_for_energy_test (delta_rho_11,delta_gamma_i_j)
 do k=1, mo_num
  do l=1, mo_num
   call potential_x_lda(k,l,rho_a,rho_b,vx_i_j)
   int_vx_lda_test += delta_gamma_i_j(k,l)*vx_i_j
  enddo
 enddo
end

 subroutine potential_x_lda(k,l,rho_a,rho_b,vx_i_j)
 implicit none
 BEGIN_DOC
 ! dEx/dn(r) = sum(k,l) {(vx_lda*mo_k*mo_l)}
 ! vx_lda = vx_lda_alpha + vx_lda_beta
 END_DOC
 double precision, intent(out) :: vx_i_j
 double precision, intent(in) :: rho_a(n_points_final_grid), rho_b(n_points_final_grid)
 integer, intent(in) :: k,l
 integer :: i
 double precision  :: r(3), weight
 double precision :: mo_k(mo_num), mo_l(mo_num)
 double precision :: e_x, vx_a, vx_b
 vx_i_j = 0.d0
  
 do i = 1, n_points_final_grid
   r(1) = final_grid_points(1,i)   
   r(2) = final_grid_points(2,i) 
   r(3) = final_grid_points(3,i)
   weight = final_weight_at_r_vector(i)

   call give_all_mos_at_r(r, mo_k)
   call ex_lda(rho_a(i),rho_b(i),e_x,vx_a,vx_b)
   
   vx_i_j += weight*(vx_a + vx_b)* mo_k(k) * mo_k(l)
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
