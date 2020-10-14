subroutine energy_xc_pben2_test (rho_a,rho_b,grad_rho_a,grad_rho_b,rho2,ecmd_pben2_test,exmd_pben2_test)

 implicit none
 BEGIN_DOC
! func_ecmd_utils/rout_pbe_ueg.irp.f
 END_DOC
 double precision, intent(in) :: rho_a(n_points_final_grid),rho_b(n_points_final_grid),grad_rho_a(3, n_points_final_grid),grad_rho_b(3, n_points_final_grid), rho2(n_points_final_grid)
 double precision, intent(out):: ecmd_pben2_test, exmd_pben2_test
 double precision :: grad_rho_a_2,grad_rho_b_2,grad_rho_a_b
 double precision :: dexdrho_a,dexdrho_b,dexdgrad_rho_a_2,dexdgrad_rho_b_2,dexdgrad_rho_a_b,ex, decdgrad_rho_2,decdrho2, decdrho
 double precision :: decdrho_a,decdrho_b,decdgrad_rho_a_2,decdgrad_rho_b_2,decdgrad_rho_a_b,ec, dexdgrad_rho_2,dexdrho2, dexdrho
 integer :: i, m 
 double precision :: weight, mu
 
 ecmd_pben2_test = 0.d0
 exmd_pben2_test = 0.d0
 do i = 1, n_points_final_grid
   mu = mu_erf_dft
   weight = final_weight_at_r_vector(i)
   grad_rho_a_2 = 0.d0
   grad_rho_b_2 = 0.d0
   grad_rho_a_b = 0.d0
   do m = 1, 3
    grad_rho_a_2 += grad_rho_a(m,i) * grad_rho_a(m,i)
    grad_rho_b_2 += grad_rho_b(m,i) * grad_rho_b(m,i)
    grad_rho_a_b += grad_rho_a(m,i) * grad_rho_b(m,i)
   enddo

!!!!!Dans le main
!  We take the extrapolated on-top pair density (Eq. 29)
!   on_top = total_cas_on_top_density(1,1) !! C'EST PAS LE BON MU ICI
!  Multiplied by 2 because of difference of normalizations between the on_top of QP2 and that of JCP, 150, 084103 1-10 (2019)
!   on_top_extrap = 2.d0 * mu_correction_of_on_top(mu,on_top)

   call ecmdsrPBEn2(mu,rho_a(i),rho_b(i),grad_rho_a_2,grad_rho_b_2,grad_rho_a_b,rho2(i),ec,decdrho_a,decdrho_b, decdrho, decdgrad_rho_a_2,decdgrad_rho_b_2,decdgrad_rho_a_b, decdgrad_rho_2,decdrho2)
   call exmdsrPBEn2(mu,rho_a(i),rho_b(i),grad_rho_a_2,grad_rho_b_2,grad_rho_a_b,rho2(i),ex,dexdrho_a,dexdrho_b, dexdrho, dexdgrad_rho_a_2,dexdgrad_rho_b_2,dexdgrad_rho_a_b, dexdgrad_rho_2,dexdrho2)

   exmd_pben2_test += weight*ex
   ecmd_pben2_test += weight*ec
 enddo   

end

 subroutine int_potential_c_pben2_two_e_test (delta_n2_11,int_vc_pben2_test)
 implicit none
 BEGIN_DOC
 ! delta_Ex[n] = int_dr delta_n(r) dEx/dn(r) 
 ! delta_n(r) <- subroutine delta_density_for_energy_test
 ! d_e_d_n    <- subroutine de_dn_at_r
 END_DOC
 double precision, intent(in) :: delta_n2_11
 double precision, intent(out):: int_vc_pben2_test
 double precision :: delta_gamma_i_j_k_l(mo_num, mo_num, mo_num, mo_num)
 integer :: i,j,k,l

 int_vc_pben2_test = 0.d0
 call delta_gamma_i_j_k_l_for_energy_test (delta_n2_11,delta_gamma_i_j_k_l)

  do i=1, mo_num
   do j=1, mo_num
    do k=1, mo_num
     do l=1, mo_num
      int_vc_pben2_test += delta_gamma_i_j_k_l(i,j,k,l)*eff_two_e
     enddo
    enddo
   enddo
  enddo
 end

 subroutine int_potential_c_pben2_one_e_test (delta_rho_11,int_vc_pben2_test)
 implicit none
 BEGIN_DOC
 ! 
 END_DOC
 integer :: i
 double precision, intent(in) :: delta_rho_11
 double precision, intent(out):: int_vc_pben2_test
 double precision :: delta_gamma_i_j(mo_num, mo_num)
 integer :: k,l

 int_vc_pben2_test = 0.d0
 call delta_gamma_i_j_for_energy_test (delta_rho_11,delta_gamma_i_j)

 do k=1, mo_num
  do l=1, mo_num
   int_vc_pben2_test += delta_gamma_i_j(k,l)*(pot_basis_alpha_mo_su_pbe_ot(k,l,1) + pot_basis_beta_mo_su_pbe_ot(k,l,1))
  enddo
 enddo
 end




!BEGIN_PROVIDER [double precision, potential_x_alpha_mo_pben2, (mo_num,mo_num)]
!implicit none
!BEGIN_DOC
!! Representation of vx_alpha in the molecular orbitals basis 
!END_DOC
!call ao_to_mo(potential_x_alpha_ao_pben2,ao_num,potential_x_alpha_mo_pben2,mo_num)
!END_PROVIDER 

!BEGIN_PROVIDER [double precision, potential_x_beta_mo_pben2, (mo_num,mo_num)]
!implicit none
!BEGIN_DOC
!! Representation of vx_beta in the molecular orbitals basis 
!END_DOC
!call ao_to_mo(potential_x_beta_ao_pben2,ao_num,potential_x_beta_mo_pben2,mo_num)
!END_PROVIDER

!
!BEGIN_PROVIDER [double precision, potential_c_alpha_mo_pben2, (mo_num,mo_num)]
!implicit none
!BEGIN_DOC
!! Representation of vx_alpha in the molecular orbitals basis 
!END_DOC
!call ao_to_mo(potential_c_alpha_ao_pben2,ao_num,potential_c_alpha_mo_pben2,mo_num)
!END_PROVIDER 

!BEGIN_PROVIDER [double precision, potential_c_beta_mo_pben2, (mo_num,mo_num)]
!implicit none
!BEGIN_DOC
!! Representation of vx_beta in the molecular orbitals basis 
!END_DOC
!call ao_to_mo(potential_c_beta_ao_pben2,ao_num,potential_c_beta_mo_pben2,mo_num)
!END_PROVIDER

 
