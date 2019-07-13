
subroutine give_epsilon_scan_provider(mu,i_point,eps_c_md_SCAN)
  implicit none
  double precision, intent(in)  :: mu 
  double precision, intent(out) :: eps_c_md_SCAN(N_states)
  integer, intent(in) :: i_point
  double precision :: two_dm, pi, e_scan,beta,on_top_two_dm_in_r_mu_corrected_from_two_dm
  double precision :: grad_rho_a(3),grad_rho_b(3)
  double precision :: grad_rho_a_2,grad_rho_b_2,grad_rho_a_b,grad_rho_2
  double precision :: ec_scan,tau
  double precision :: delta,two_dm_corr,rho_a,rho_b
  double precision :: denom,g0_UEG_mu_inf
  integer :: m, istate

  pi = 4.d0 * datan(1.d0)

  eps_c_md_SCAN = 0.d0
  do istate = 1, N_states
   ! total and spin density 
   ! gradients of the effective spin density 
   grad_rho_a_2 = 0.D0
   grad_rho_b_2 = 0.D0
   grad_rho_a_b = 0.D0
   do m = 1, 3
    grad_rho_a_2 += one_e_dm_and_grad_alpha_in_r(m,i_point,istate)**2.d0
    grad_rho_b_2 += one_e_dm_and_grad_beta_in_r(m,i_point,istate) **2.d0
    grad_rho_a_b += one_e_dm_and_grad_alpha_in_r(m,i_point,istate) * one_e_dm_and_grad_beta_in_r(m,i_point,istate)
   enddo
   grad_rho_2 = grad_rho_a_2 + grad_rho_b_2 + 2.d0 * grad_rho_a_b
   rho_a = one_e_dm_and_grad_alpha_in_r(4,i_point,istate)
   rho_b = one_e_dm_and_grad_beta_in_r(4,i_point,istate)
   tau = kinetic_density_generalized(i_point)
   e_SCAN = (rho_a + rho_b) * ec_scan(rho_a,rho_b,tau,grad_rho_2)
   if(isnan(e_SCAN))then
    print*,'e_SCAN is NAN !!'
    write(*,'(100(F10.5,X))')rho_a,rho_b,tau,grad_rho_2
   endif
   if(mu == 0.d0) then
    eps_c_md_SCAN(istate)=e_SCAN
   else
!   note: the on-top pair density is (1-zeta^2) rhoc^2 g0 = 4 rhoa * rhob * g0
    denom = (-2.d0+dsqrt(2d0))*sqrt(2.d0*pi) * 4.d0*rho_a*rho_b*g0_UEG_mu_inf(rho_a,rho_b)
    if (dabs(denom) > 1.d-12) then
     beta = (3.d0*e_SCAN)/denom
     ! Ecmd functional with the UEG ontop pair density when mu -> infty 
     ! and the usual SCAN correlation energy when mu = 0
     eps_c_md_SCAN(istate)=e_SCAN/(1.d0+beta*mu**3)
    else
     eps_c_md_SCAN(istate)=0.d0
    endif
   endif
  enddo
end

subroutine give_epsilon_scan_ontop_provider(mu,i_point,eps_c_md_ontop_SCAN)
  implicit none
  double precision, intent(in)  :: mu 
  double precision, intent(out) :: eps_c_md_ontop_SCAN(N_states)
  integer, intent(in) :: i_point
  double precision :: two_dm, pi, e_scan,beta,on_top_two_dm_in_r_mu_corrected_from_two_dm
  double precision :: grad_rho_a(3),grad_rho_b(3)
  double precision :: grad_rho_a_2,grad_rho_b_2,grad_rho_a_b,grad_rho_2
  double precision :: ec_scan,tau
  double precision :: delta,two_dm_corr,rho_a,rho_b
  integer :: m, istate

  pi = 4.d0 * datan(1.d0)

  eps_c_md_ontop_SCAN = 0.d0
  do istate = 1, N_states
   ! gradients of the density 
   grad_rho_a_2 = 0.D0
   grad_rho_b_2 = 0.D0
   grad_rho_a_b = 0.D0
   do m = 1, 3
    grad_rho_a_2 += one_e_dm_and_grad_alpha_in_r(m,i_point,istate)**2.d0
    grad_rho_b_2 += one_e_dm_and_grad_beta_in_r(m,i_point,istate) **2.d0
    grad_rho_a_b += one_e_dm_and_grad_alpha_in_r(m,i_point,istate) * one_e_dm_and_grad_beta_in_r(m,i_point,istate)
   enddo
   grad_rho_2 = grad_rho_a_2 + grad_rho_b_2 + 2.d0 * grad_rho_a_b
   rho_a = one_e_dm_and_grad_alpha_in_r(4,i_point,istate)
   rho_b = one_e_dm_and_grad_beta_in_r(4,i_point,istate)
   tau = kinetic_density_generalized(i_point)
   e_SCAN = (rho_a + rho_b) * ec_scan(rho_a,rho_b,tau,grad_rho_2)
   if(isnan(e_SCAN))then
    print*,'e_SCAN is NAN !!'
    write(*,'(100(F10.5,X))')rho_a,rho_b,tau,grad_rho_2
   endif

   two_dm = core_inact_act_on_top_of_r(i_point,istate) ! on top of the wave function 
   two_dm_corr = on_top_two_dm_in_r_mu_corrected_from_two_dm(mu,istate,two_dm) ! extrapolated "exact" on top
   if(dabs(( (-2.d0+sqrt(2d0))*sqrt(2.d0*pi)*2.d0*two_dm_corr )).lt.1.d-12)cycle
   beta = dabs((3.d0*e_SCAN)/( (-2.d0+sqrt(2d0))*sqrt(2.d0*pi)*2.d0*two_dm_corr ))
   ! Ecmd functional with the extrapolated exact on top when mu -> infty 
   ! and the usual PBE functional when mu = 0
   eps_c_md_ontop_SCAN(istate)=e_SCAN/(1.d0+beta*mu**3)
  enddo

end

