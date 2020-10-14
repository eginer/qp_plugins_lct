subroutine energy_xc_pbeUEG_test (rho_a,rho_b,grad_rho_a,grad_rho_b,ex_pbeUEG_test, ec_pbeUEG_test)

 implicit none
 BEGIN_DOC
! func_ecmd_utils/rout_pbe_ueg.irp.f
 END_DOC
 double precision, intent(in) :: rho_a(n_points_final_grid),rho_b(n_points_final_grid),grad_rho_a(3, n_points_final_grid),grad_rho_b(3, n_points_final_grid)
 double precision, intent(out):: ex_pbeUEG_test, ec_pbeUEG_test
 double precision :: grad_rho_a_2,grad_rho_b_2,grad_rho_a_b
 double precision :: dexdrho_a,dexdrho_b,dexdgrad_rho_a_2,dexdgrad_rho_b_2,dexdgrad_rho_a_b,ex
 double precision :: decdrho_a,decdrho_b,decdgrad_rho_a_2,decdgrad_rho_b_2,decdgrad_rho_a_b,ec
 integer :: i, m
 double precision :: weight

 ex_pbeUEG_test = 0.d0
 ec_pbeUEG_test = 0.d0

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

   call exc_dexc_md_sr_PBE(mu_erf_dft,rho_a(i),rho_b(i),grad_rho_a_2,grad_rho_b_2,grad_rho_a_b, &
       ex,dexdrho_a,dexdrho_b,dexdgrad_rho_a_2,dexdgrad_rho_b_2,dexdgrad_rho_a_b, &
       ec,decdrho_a,decdrho_b,decdgrad_rho_a_2,decdgrad_rho_b_2,decdgrad_rho_a_b)

   ec_pbeUEG_test += ec*weight
   ex_pbeUEG_test += ex*weight
 enddo
end

subroutine int_potential_xc_pbeUEG_test (delta_rho_11,int_vx_pbeUEG_test, int_vc_pbeUEG_test)
 implicit none
 BEGIN_DOC
! delta_Ex[n] = int_dr delta_n(r) dEx/dn(r) 
! delta_n(r) <- subroutine delta_density_for_energy_test
! d_e_d_n    <- subroutine de_dn_at_r
 END_DOC
 integer :: i
 double precision, intent(in) :: delta_rho_11
 double precision, intent(out):: int_vx_pbeUEG_test, int_vc_pbeUEG_test
 double precision :: vx_i_j, vc_i_j
 double precision :: delta_gamma_i_j(mo_num, mo_num)
 integer :: k,l
! double precision, intent(in) :: rho_a(n_points_final_grid), rho_b(n_points_final_grid), grad_rho_a(3,n_points_final_grid), grad_rho_b(3,n_points_final_grid)
 double precision :: ex_pbeUEG_test, ec_pbeUEG_test, vx_pbeUEG_test, vc_pbeUEG_test

 int_vx_pbeUEG_test = 0.d0
 int_vc_pbeUEG_test = 0.d0
 call delta_gamma_i_j_for_energy_test (delta_rho_11,delta_gamma_i_j)

 do k=1, mo_num
  do l=1, mo_num
!   call energy_xc_pbeUEG_test (rho_a, rho_b, grad_rho_a, grad_rho_b, k, l, ex_pbeUEG_test, ec_pbeUEG_test, vx_pbeUEG_test, vc_pbeUEG_test)
   int_vx_pbeUEG_test += delta_gamma_i_j(k,l)*(potential_x_alpha_mo_md_sr_pbe(k,l) + potential_x_beta_mo_md_sr_pbe(k,l))
   int_vc_pbeUEG_test += delta_gamma_i_j(k,l)*(potential_c_alpha_mo_md_sr_pbe(k,l) + potential_c_beta_mo_md_sr_pbe(k,l))
  enddo
 enddo
end


BEGIN_PROVIDER [double precision, potential_x_alpha_mo_md_sr_pbe, (mo_num,mo_num)]
 implicit none
 BEGIN_DOC
 ! Representation of vx_alpha in the molecular orbitals basis 
 END_DOC
 call ao_to_mo(potential_x_alpha_ao_md_sr_pbe,ao_num,potential_x_alpha_mo_md_sr_pbe,mo_num)
END_PROVIDER 

BEGIN_PROVIDER [double precision, potential_x_beta_mo_md_sr_pbe, (mo_num,mo_num)]
 implicit none
 BEGIN_DOC
 ! Representation of vx_beta in the molecular orbitals basis 
 END_DOC
 call ao_to_mo(potential_x_beta_ao_md_sr_pbe,ao_num,potential_x_beta_mo_md_sr_pbe,mo_num)
END_PROVIDER

 
BEGIN_PROVIDER [double precision, potential_c_alpha_mo_md_sr_pbe, (mo_num,mo_num)]
 implicit none
 BEGIN_DOC
 ! Representation of vx_alpha in the molecular orbitals basis 
 END_DOC
 call ao_to_mo(potential_c_alpha_ao_md_sr_pbe,ao_num,potential_c_alpha_mo_md_sr_pbe,mo_num)
END_PROVIDER 

BEGIN_PROVIDER [double precision, potential_c_beta_mo_md_sr_pbe, (mo_num,mo_num)]
 implicit none
 BEGIN_DOC
 ! Representation of vx_beta in the molecular orbitals basis 
 END_DOC
 call ao_to_mo(potential_c_beta_ao_md_sr_pbe,ao_num,potential_c_beta_mo_md_sr_pbe,mo_num)
END_PROVIDER

 
