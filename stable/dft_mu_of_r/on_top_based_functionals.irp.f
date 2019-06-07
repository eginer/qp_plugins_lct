
subroutine give_epsilon_pbe_ontop_effective_spin_dens_provider(mu,i_point,eps_c_md_on_top_PBE)
  implicit none
  double precision, intent(in)  :: mu 
  double precision, intent(out) :: eps_c_md_on_top_PBE(N_states)
  integer, intent(in) :: i_point
  double precision :: two_dm, pi, beta,on_top_two_dm_in_r_mu_corrected_from_two_dm
  double precision :: grad_rho_a(3),grad_rho_b(3)
  double precision :: grad_rho_a_2,grad_rho_b_2,grad_rho_a_b
  double precision :: rhoc,rhoo,sigmacc,sigmaco,sigmaoo,e_PBE
  integer :: m, istate

  pi = 4.d0 * datan(1.d0)

  eps_c_md_on_top_PBE = 0.d0
  do istate = 1, N_states
   ! total density 
   rhoc = one_e_dm_and_grad_alpha_in_r(4,i_point,istate) + one_e_dm_and_grad_beta_in_r(4,i_point,istate) 
   ! effective spin density 
   rhoo = effective_spin_dm(i_point,istate)
   grad_rho_a_2 = 0.D0
   grad_rho_b_2 = 0.D0
   grad_rho_a_b = 0.D0
   ! gradients of the effective spin density 
   do m = 1, 3
    grad_rho_a_2 += one_e_dm_and_grad_alpha_in_r(m,i_point,istate)**2.d0
    grad_rho_b_2 += one_e_dm_and_grad_beta_in_r(m,i_point,istate) **2.d0
    grad_rho_a_b += one_e_dm_and_grad_alpha_in_r(m,i_point,istate) * one_e_dm_and_grad_beta_in_r(m,i_point,istate)
   enddo
   sigmacc = grad_rho_a_2 + grad_rho_b_2 + 2.d0 * grad_rho_a_b
   sigmaco = 0.d0
   sigmaoo = 0.d0
   ! Note EG: for the pbe correlation energy, we don't need the gradients of the effective spin density 
   !          as it reduces to the usual square of the gradient of the total density 
   call ec_pbe_only(0.d0,rhoc,rhoo,sigmacc,sigmaco,sigmaoo,e_PBE) !  PBE correlation energy with effective spin density 
   double precision :: delta,two_dm_corr
   two_dm = core_inact_act_on_top_of_r(i_point,istate) ! on top of the wave function 
   two_dm_corr = on_top_two_dm_in_r_mu_corrected_from_two_dm(mu,istate,two_dm) ! extrapolated "exact" on top
   if(dabs(( (-2.d0+sqrt(2d0))*sqrt(2.d0*pi)*2.d0*two_dm_corr )).lt.1.d-12)cycle
   beta = (3.d0*e_PBE)/( (-2.d0+sqrt(2d0))*sqrt(2.d0*pi)*2.d0*two_dm_corr ) 
   ! Ecmd functional with the extrapolated exact on top when mu -> infty 
   ! and the usual PBE functional when mu = 0
   eps_c_md_on_top_PBE(istate)=e_PBE/(1.d0+beta*mu**3)
  enddo
end


subroutine give_epsilon_pbe_effective_spin_dens_provider(mu,i_point,eps_c_md_PBE)
  implicit none
  double precision, intent(in)  :: mu 
  double precision, intent(out) :: eps_c_md_PBE(N_states)
  integer, intent(in) :: i_point
  double precision :: two_dm, pi, beta,on_top_two_dm_in_r_mu_corrected_from_two_dm
  double precision :: grad_rho_a(3),grad_rho_b(3)
  double precision :: grad_rho_a_2,grad_rho_b_2,grad_rho_a_b
  double precision :: rhoc,rhoo,sigmacc,sigmaco,sigmaoo
  double precision :: delta,two_dm_corr,rho_a,rho_b,denom,e_PBE
  double precision :: g0_UEG_mu_inf
  integer :: m, istate

  pi = 4.d0 * datan(1.d0)

  eps_c_md_PBE = 0.d0
  do istate = 1, N_states
   ! total density 
   rhoc = one_e_dm_and_grad_alpha_in_r(4,i_point,istate) + one_e_dm_and_grad_beta_in_r(4,i_point,istate) 
   ! effective spin density 
   rhoo = effective_spin_dm(i_point,istate)
   rho_a  = effective_alpha_dm(i_point,istate)
   rho_b  = effective_beta_dm(i_point,istate)
   grad_rho_a_2 = 0.D0
   grad_rho_b_2 = 0.D0
   grad_rho_a_b = 0.D0
   ! gradients of the effective spin density 
   do m = 1, 3
    grad_rho_a_2 += one_e_dm_and_grad_alpha_in_r(m,i_point,istate)**2.d0
    grad_rho_b_2 += one_e_dm_and_grad_beta_in_r(m,i_point,istate) **2.d0
    grad_rho_a_b += one_e_dm_and_grad_alpha_in_r(m,i_point,istate) * one_e_dm_and_grad_beta_in_r(m,i_point,istate)
   enddo
   sigmacc = grad_rho_a_2 + grad_rho_b_2 + 2.d0 * grad_rho_a_b
   sigmaco = 0.d0
   sigmaoo = 0.d0
   ! Note EG: for the pbe correlation energy, we don't need the gradients of the effective spin density 
   !          as it reduces to the usual square of the gradient of the total density 
   call ec_pbe_only(0.d0,rhoc,rhoo,sigmacc,sigmaco,sigmaoo,e_PBE) !  PBE correlation energy with effective spin density 
   if(mu == 0.d0) then
    eps_c_md_PBE(istate)=e_PBE
   else
!   note: the on-top pair density is (1-zeta^2) rhoc^2 g0 = 4 rhoa * rhob * g0
    denom = (-2.d0+dsqrt(2d0))*sqrt(2.d0*pi) * 4.d0*rho_a*rho_b*g0_UEG_mu_inf(rho_a,rho_b)
    if (dabs(denom) > 1.d-12) then
     beta = (3.d0*e_PBE)/denom
     ! Ecmd functional with the UEG ontop pair density when mu -> infty 
     ! and the usual PBE correlation energy when mu = 0
     eps_c_md_PBE(istate)=e_PBE/(1.d0+beta*mu**3)
    else
     eps_c_md_PBE(istate)=0.d0
    endif
   endif
  enddo

end

subroutine give_epsilon_lyp_ontop_effective_spin_dens_provider(mu,i_point,eps_c_md_ontop_LYP)
  implicit none
  double precision, intent(in)  :: mu 
  double precision, intent(out) :: eps_c_md_ontop_LYP(N_states)
  integer, intent(in) :: i_point
  double precision :: two_dm, pi, e_lyp,beta,on_top_two_dm_in_r_mu_corrected_from_two_dm
  double precision :: grad_rho_a(3),grad_rho_b(3)
  double precision :: grad_rho_a_2,grad_rho_b_2,grad_rho_a_b
  double precision :: rhoc,rhoo,ec_lyp_88
  double precision :: delta,two_dm_corr,rho_a,rho_b
  double precision :: grad_rho_2
  integer :: m, istate

  pi = 4.d0 * datan(1.d0)

  eps_c_md_ontop_LYP = 0.d0
  do istate = 1, N_states
   ! total density 
   rhoc = one_e_dm_and_grad_alpha_in_r(4,i_point,istate) + one_e_dm_and_grad_beta_in_r(4,i_point,istate) 
   ! effective spin density 
   rhoo = effective_spin_dm(i_point,istate)
   ! gradients of the effective spin density 
   grad_rho_a_2 = 0.D0
   grad_rho_b_2 = 0.D0
   grad_rho_a_b = 0.D0
   do m = 1, 3
    grad_rho_a_2 += grad_effective_alpha_dm(m,i_point,istate)**2.d0
    grad_rho_b_2 += grad_effective_beta_dm(m,i_point,istate) **2.d0
    grad_rho_a_b += grad_effective_alpha_dm(m,i_point,istate) * grad_effective_beta_dm(m,i_point,istate)
   enddo
   grad_rho_2 = grad_rho_a_2 + grad_rho_b_2 + 2.D0 * grad_rho_a_b
   rho_a = effective_alpha_dm(i_point,istate)
   rho_b = effective_beta_dm(i_point,istate)


   e_LYP = ec_lyp_88(rhoc,rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_2)
   two_dm = core_inact_act_on_top_of_r(i_point,istate) ! on top of the wave function 
   two_dm_corr = on_top_two_dm_in_r_mu_corrected_from_two_dm(mu,istate,two_dm) ! extrapolated "exact" on top
   if(dabs(( (-2.d0+sqrt(2d0))*sqrt(2.d0*pi)*2.d0*two_dm_corr )).lt.1.d-12)cycle
   beta = (3.d0*e_LYP)/( (-2.d0+sqrt(2d0))*sqrt(2.d0*pi)*2.d0*two_dm_corr )
   ! Ecmd functional with the extrapolated exact on top when mu -> infty 
   ! and the usual PBE functional when mu = 0
   eps_c_md_ontop_LYP(istate)=e_LYP/(1.d0+beta*mu**3)
  enddo

end

subroutine give_epsilon_lyp_effective_spin_dens_provider(mu,i_point,eps_c_md_LYP)
  implicit none
  double precision, intent(in)  :: mu 
  double precision, intent(out) :: eps_c_md_LYP(N_states)
  integer, intent(in) :: i_point
  double precision :: two_dm, pi, e_lyp,beta,on_top_two_dm_in_r_mu_corrected_from_two_dm
  double precision :: grad_rho_a(3),grad_rho_b(3)
  double precision :: grad_rho_a_2,grad_rho_b_2,grad_rho_a_b
  double precision :: rhoc,rhoo,ec_lyp_88
  double precision :: delta,two_dm_corr,rho_a,rho_b
  double precision :: grad_rho_2,denom,g0_UEG_mu_inf
  integer :: m, istate

  pi = 4.d0 * datan(1.d0)

  eps_c_md_LYP = 0.d0
  do istate = 1, N_states
   ! total density 
   rhoc = one_e_dm_and_grad_alpha_in_r(4,i_point,istate) + one_e_dm_and_grad_beta_in_r(4,i_point,istate) 
   ! effective spin density 
   rhoo = effective_spin_dm(i_point,istate)
!  rhoo = one_e_dm_and_grad_alpha_in_r(4,i_point,istate) - one_e_dm_and_grad_beta_in_r(4,i_point,istate) 
   ! gradients of the effective spin density 
   grad_rho_a_2 = 0.D0
   grad_rho_b_2 = 0.D0
   grad_rho_a_b = 0.D0
   do m = 1, 3
    grad_rho_a_2 += grad_effective_alpha_dm(m,i_point,istate)**2.d0
    grad_rho_b_2 += grad_effective_beta_dm(m,i_point,istate) **2.d0
    grad_rho_a_b += grad_effective_alpha_dm(m,i_point,istate) * grad_effective_beta_dm(m,i_point,istate)
   enddo
   grad_rho_2 = grad_rho_a_2 + grad_rho_b_2 + 2.D0 * grad_rho_a_b
   rho_a = effective_alpha_dm(i_point,istate)
   rho_b = effective_beta_dm(i_point,istate)


   e_LYP = ec_lyp_88(rhoc,rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_2)
   if(mu == 0.d0) then
    eps_c_md_LYP(istate)=e_LYP
   else
!   note: the on-top pair density is (1-zeta^2) rhoc^2 g0 = 4 rhoa * rhob * g0
    denom = (-2.d0+dsqrt(2d0))*sqrt(2.d0*pi) * 4.d0*rho_a*rho_b*g0_UEG_mu_inf(rho_a,rho_b)
    if (dabs(denom) > 1.d-12) then
     beta = (3.d0*e_LYP)/denom
     ! Ecmd functional with the UEG ontop pair density when mu -> infty 
     ! and the usual LYP correlation energy when mu = 0
     eps_c_md_LYP(istate)=e_LYP/(1.d0+beta*mu**3)
    else
     eps_c_md_LYP(istate)=0.d0
    endif
   endif
  enddo

end

