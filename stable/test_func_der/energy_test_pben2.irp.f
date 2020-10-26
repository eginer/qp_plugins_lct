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
 double precision :: on_top_extrap, mu_correction_of_on_top 
 ecmd_pben2_test = 0.d0
 exmd_pben2_test = 0.d0
 do i = 1, n_points_final_grid
  mu = mu_of_r_prov(i,1)
!   mu = mu_of_r_hf(i)
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
   on_top_extrap = 2.d0 * mu_correction_of_on_top(mu,rho2(i))
!   on_top_extrap = rho2(i)
!   on_top_extrap_alpha = rho2_alpha(i)
!   on_top_extrap_beta = rho2_beta(i)

   call ecmdsrPBEn2(mu,rho_a(i),rho_b(i),grad_rho_a_2,grad_rho_b_2,grad_rho_a_b,on_top_extrap,ec,decdrho_a,decdrho_b, decdrho, decdgrad_rho_a_2,decdgrad_rho_b_2,decdgrad_rho_a_b, decdgrad_rho_2,decdrho2)
   call exmdsrPBEn2(mu,rho_a(i),rho_b(i),grad_rho_a_2,grad_rho_b_2,grad_rho_a_b,on_top_extrap,ex,dexdrho_a,dexdrho_b, dexdrho, dexdgrad_rho_a_2,dexdgrad_rho_b_2,dexdgrad_rho_a_b, dexdgrad_rho_2,dexdrho2)

   exmd_pben2_test += weight*ex
   ecmd_pben2_test += weight*ec
 enddo   

end

