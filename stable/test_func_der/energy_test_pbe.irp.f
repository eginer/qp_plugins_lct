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

 
