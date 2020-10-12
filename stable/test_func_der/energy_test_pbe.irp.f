subroutine energy_xc_pbe_test (rho_a, rho_b, grad_rho_a, grad_rho_b, ex_pbe_test, ec_pbe_test)
 implicit none
 BEGIN_DOC
! exchange energy with the lda functional
 END_DOC
 integer :: i,j,m
 double precision, intent(in) :: rho_a(n_points_final_grid), rho_b(n_points_final_grid), grad_rho_a(3,n_points_final_grid), grad_rho_b(3,n_points_final_grid)
 double precision, intent(out) :: ex_pbe_test, ec_pbe_test
 double precision :: mu,weight
 double precision :: ex, ec
 double precision :: grad_rho_a_2,grad_rho_b_2,grad_rho_a_b
 double precision :: vc_rho_a, vc_rho_b, vx_rho_a, vx_rho_b
 double precision :: vx_grad_rho_a_2, vx_grad_rho_b_2, vx_grad_rho_a_b, vc_grad_rho_a_2, vc_grad_rho_b_2, vc_grad_rho_a_b

 ex_pbe_test = 0.d0
 ec_pbe_test = 0.d0
 mu = 0.d0
  do i = 1, n_points_final_grid
   weight = final_weight_at_r_vector(i)
   grad_rho_a_2 = 0.d0
   grad_rho_b_2 = 0.d0
   grad_rho_a_b = 0.d0
   do m = 1, 3
    grad_rho_a_2 += grad_rho_a(m,i) * grad_rho_a(m,i)
    grad_rho_b_2 += grad_rho_b(m,i) * grad_rho_b(m,i)
    grad_rho_a_b += grad_rho_a(m,i) * grad_rho_b(m,i)
   enddo

                             ! inputs
   call GGA_sr_type_functionals(mu,rho_a(i),rho_b(i),grad_rho_a_2,grad_rho_b_2,grad_rho_a_b,                 &  ! outputs exchange
                             ex,vx_rho_a,vx_rho_b,vx_grad_rho_a_2,vx_grad_rho_b_2,vx_grad_rho_a_b, &  ! outputs correlation
                             ec,vc_rho_a,vc_rho_b,vc_grad_rho_a_2,vc_grad_rho_b_2,vc_grad_rho_a_b  )
   ex_pbe_test += ex * weight
   ec_pbe_test += ec * weight
 
  enddo
end

subroutine int_potential_xc_pbe_test (delta_rho_11, int_vx_pbe_test, int_vc_pbe_test)
 implicit none
 BEGIN_DOC
! delta_Ex[n] = int_dr delta_n(r) dEx/dn(r) 
! delta_n(r) <- subroutine delta_density_for_energy_test
! d_e_d_n    <- subroutine de_dn_at_r
 END_DOC
 integer :: i
 double precision, intent(in) :: delta_rho_11
 double precision, intent(out):: int_vx_pbe_test, int_vc_pbe_test
 double precision :: vx_i_j, vc_i_j
 double precision :: delta_gamma_i_j(mo_num, mo_num)
 integer :: k,l
 
 int_vx_pbe_test = 0.d0
 int_vc_pbe_test = 0.d0
 call delta_gamma_i_j_for_energy_test (delta_rho_11,delta_gamma_i_j)
 do k=1, mo_num
  do l=1, mo_num
   call potential_xc_pbe(k,l,vx_i_j, vc_i_j)
   int_vx_pbe_test += delta_gamma_i_j(k,l)*vx_i_j
   int_vc_pbe_test += delta_gamma_i_j(k,l)*vc_i_j
  enddo
 enddo
end

subroutine potential_xc_pbe(k,l,vx_i_j, vc_i_j)
 implicit none
 BEGIN_DOC
! dEx/dn(r) = sum(k,l) {(vx_lda*mo_k*mo_l)}
! vx_lda = vx_lda_alpha + vx_lda_beta
 END_DOC
 double precision, intent(out) :: vx_i_j, vc_i_j
 integer, intent(in) :: k,l
 integer :: i
 double precision  :: r(3), weight
 double precision :: mo_k(mo_num), mo_l(mo_num)
 vx_i_j = 0.d0
 vc_i_j = 0.d0
  
 do i = 1, n_points_final_grid
   r(1) = final_grid_points(1,i)   
   r(2) = final_grid_points(2,i) 
   r(3) = final_grid_points(3,i)
   weight = final_weight_at_r_vector(i)

   call give_all_mos_at_r(r, mo_k)
   call give_all_mos_at_r(r, mo_l)

   vx_i_j += weight*(potential_x_alpha_mo_pbe(k,l) + potential_x_beta_mo_pbe(k,l)) * mo_k(k) * mo_l(l)! * 0.5d0
   vc_i_j += weight*(potential_c_alpha_mo_pbe(k,l) + potential_c_beta_mo_pbe(k,l)) * mo_k(k) * mo_l(l)! * 0.5d0
 enddo
end


BEGIN_PROVIDER [double precision, potential_x_alpha_mo_pbe, (mo_num,mo_num)]
 implicit none
 BEGIN_DOC
 ! Representation of vx_alpha in the molecular orbitals basis 
 END_DOC
 call ao_to_mo(potential_x_alpha_ao_pbe,ao_num,potential_x_alpha_mo_pbe,mo_num)
END_PROVIDER 

BEGIN_PROVIDER [double precision, potential_x_beta_mo_pbe, (mo_num,mo_num)]
 implicit none
 BEGIN_DOC
 ! Representation of vx_beta in the molecular orbitals basis 
 END_DOC
 call ao_to_mo(potential_x_beta_ao_pbe,ao_num,potential_x_beta_mo_pbe,mo_num)
END_PROVIDER

 
BEGIN_PROVIDER [double precision, potential_c_alpha_mo_pbe, (mo_num,mo_num)]
 implicit none
 BEGIN_DOC
 ! Representation of vx_alpha in the molecular orbitals basis 
 END_DOC
 call ao_to_mo(potential_c_alpha_ao_pbe,ao_num,potential_c_alpha_mo_pbe,mo_num)
END_PROVIDER 

BEGIN_PROVIDER [double precision, potential_c_beta_mo_pbe, (mo_num,mo_num)]
 implicit none
 BEGIN_DOC
 ! Representation of vx_beta in the molecular orbitals basis 
 END_DOC
 call ao_to_mo(potential_c_beta_ao_pbe,ao_num,potential_c_beta_mo_pbe,mo_num)
END_PROVIDER

 
