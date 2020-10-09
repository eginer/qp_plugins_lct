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

subroutine int_potential_x_lda_test (delta_dm_a, delta_dm_b, int_vx_lda_test)
 implicit none
 BEGIN_DOC
! delta_Ex[n] = int_dr delta_n(r) dEx/dn(r) 
! delta_n(r) <- subroutine delta_density_for_energy_test
! d_e_d_n    <- subroutine de_dn_at_r
 END_DOC
 integer :: i
 double precision, intent(in) :: delta_dm_a(n_points_final_grid), delta_dm_b(n_points_final_grid)
 double precision, intent(out):: int_vx_lda_test
 double precision :: r(3), d_e_d_n, weight
 
 int_vx_lda_test = 0.d0
  do i = 1, n_points_final_grid
   r(1) = final_grid_points(1,i)   
   r(2) = final_grid_points(2,i) 
   r(3) = final_grid_points(3,i)
   weight = final_weight_at_r_vector(i)

   call de_dn_at_r(r, d_e_d_n)
   int_vx_lda_test += (delta_dm_a(i) + delta_dm_b(i))*weight*d_e_d_n*0.5d0
  enddo
end

subroutine de_dn_at_r(r, d_e_d_n)
 implicit none
 BEGIN_DOC
! dEx/dn(r) = sum(k,l) {(vx_lda*mo_k*mo_l)}
! vx_lda = vx_lda_alpha + vx_lda_beta
 END_DOC
 double precision, intent(out) :: d_e_d_n
 double precision, intent(in)  :: r(3)
 integer :: k,l
 double precision :: mo_k(mo_num), mo_l(mo_num)

 call give_all_mos_at_r(r, mo_k)
 call give_all_mos_at_r(r, mo_l)
 d_e_d_n = 0.d0

 do k = 1, mo_num
  do l = 1, mo_num
   d_e_d_n += (potential_x_alpha_mo_lda(l,k) + potential_x_beta_mo_lda(l,k)) * mo_k(k) * mo_l(l)! * 0.5d0
  enddo
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

   delta_rho_a(i) = delta_rho_11 * mo_i(1) * mo_i(1)  *0.5d0 !delta_rho_ij = delta_rho_11 = kro(i1)kro(j1)*epsilon
   delta_rho_b(i) = delta_rho_11 * mo_i(1) * mo_i(1)  *0.5d0

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
