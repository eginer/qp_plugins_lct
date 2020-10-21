subroutine energy_xc_pbe_test (rho_a, rho_b, grad_rho_a, grad_rho_b, ex_pbe_test, ec_pbe_test, vx_pbe_test, vc_pbe_test)
 implicit none
 BEGIN_DOC
! exchange energy with the lda functional
 END_DOC
 integer :: i,j,m
 double precision, intent(in) :: rho_a(n_points_final_grid), rho_b(n_points_final_grid), grad_rho_a(3,n_points_final_grid), grad_rho_b(3,n_points_final_grid)
 double precision, intent(out) :: ex_pbe_test, ec_pbe_test, vx_pbe_test, vc_pbe_test
 double precision :: mu,weight, r(3)
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

subroutine int_potential_xc_pbe_test (delta_rho_11_alpha, delta_rho_11_beta,int_vx_pbe_test, int_vc_pbe_test)
 implicit none
 BEGIN_DOC
! delta_Ex[n] = int_dr delta_n(r) dEx/dn(r) 
! delta_n(r) <- subroutine delta_density_for_energy_test
! d_e_d_n    <- subroutine de_dn_at_r
 END_DOC
 integer :: i
 double precision, intent(in) :: delta_rho_11_alpha, delta_rho_11_beta
 double precision, intent(out):: int_vx_pbe_test, int_vc_pbe_test
 double precision :: vx_i_j, vc_i_j
 double precision :: delta_gamma_i_j_alpha(mo_num, mo_num)
 double precision :: delta_gamma_i_j_beta(mo_num, mo_num)
 integer :: k,l
! double precision, intent(in) :: rho_a(n_points_final_grid), rho_b(n_points_final_grid), grad_rho_a(3,n_points_final_grid), grad_rho_b(3,n_points_final_grid)
 double precision :: ex_pbe_test, ec_pbe_test, vx_pbe_test, vc_pbe_test

 int_vx_pbe_test = 0.d0
 int_vc_pbe_test = 0.d0
 call delta_gamma_i_j_for_energy_test_alpha_beta_general (delta_rho_11_alpha,delta_rho_11_beta,delta_gamma_i_j_alpha,delta_gamma_i_j_beta)

 do k=1, mo_num
  do l=1, mo_num
!   call energy_xc_pbe_test (rho_a, rho_b, grad_rho_a, grad_rho_b, k, l, ex_pbe_test, ec_pbe_test, vx_pbe_test, vc_pbe_test)
   int_vx_pbe_test += delta_gamma_i_j_alpha(k,l)*potential_x_alpha_mo_pbe(k,l) + delta_gamma_i_j_beta(k,l)*potential_x_beta_mo_pbe(k,l)
   int_vc_pbe_test += delta_gamma_i_j_alpha(k,l)*potential_c_alpha_mo_pbe(k,l) + delta_gamma_i_j_beta(k,l)*potential_c_beta_mo_pbe(k,l)
  enddo
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

 
