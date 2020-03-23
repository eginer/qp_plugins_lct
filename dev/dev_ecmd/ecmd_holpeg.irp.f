
subroutine give_epsilon_holpeg_provider(mu,i_point,eps_c_md_holpeg)
  implicit none
  double precision, intent(in)  :: mu 
  double precision, intent(out) :: eps_c_md_holpeg(N_states)
  integer, intent(in) :: i_point
  double precision :: two_dm, pi, e_holpeg,beta,on_top_two_dm_in_r_mu_corrected_from_two_dm
  double precision :: grad_rho_a(3),grad_rho_b(3)
  double precision :: grad_rho_a_2,grad_rho_b_2,grad_rho_a_b
  double precision :: rhoc,rhoo,ec_holpeg_88
  double precision :: delta,two_dm_corr,rho_a,rho_b
  double precision :: grad_rho_2,denom,g0_UEG_mu_inf
  double precision :: sigmacc,sigmaco,sigmaoo
  integer :: m, istate

  pi = 4.d0 * datan(1.d0)

  eps_c_md_holpeg = 0.d0
  do istate = 1, N_states
   ! total and spin density 
   double precision :: on_top,lapl_aa,lapl_bb,on_top_aa_lapl,on_top_bb_lapl,ec_ab,ec_aa,ec_bb
   rho_a = one_e_dm_and_grad_alpha_in_r(4,i_point,istate) 
   rho_b = one_e_dm_and_grad_beta_in_r(4,i_point,istate) 
   rhoc = rho_a + rho_b
   on_top = rho_a * rho_b
   lapl_aa   = on_top_aa_lapl(i_point) 
   lapl_bb   = on_top_bb_lapl(i_point) 
   call give_ec_holpeg(rhoc,on_top,lapl_aa,lapl_bb,ec_ab,ec_aa,ec_bb)
   e_holpeg = ec_ab+ec_aa+ec_bb
   if(mu == 0.d0) then
    eps_c_md_holpeg(istate)=e_holpeg
   else
!   note: the on-top pair density is (1-zeta^2) rhoc^2 g0 = 4 rhoa * rhob * g0
    denom = (-2.d0+dsqrt(2d0))*sqrt(2.d0*pi) * 4.d0*rho_a*rho_b*g0_UEG_mu_inf(rho_a,rho_b)
    if (dabs(denom) > 1.d-12) then
     beta = (3.d0*e_holpeg)/denom
     ! Ecmd functional with the UEG ontop pair density when mu -> infty 
     ! and the usual holpeg correlation energy when mu = 0
     eps_c_md_holpeg(istate)=e_holpeg/(1.d0+beta*mu**3)
    else
     eps_c_md_holpeg(istate)=0.d0
    endif
   endif
  enddo
end
