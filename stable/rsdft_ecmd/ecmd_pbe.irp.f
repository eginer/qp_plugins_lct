 BEGIN_PROVIDER [double precision, Energy_c_md_PBE_mu_vector, (N_states)]
 BEGIN_DOC
  ! Give the Ec_md energy with a large mu behaviour in function of the UEG on top pair density coupled to the PBE correlation energy at mu=0
  ! Ec_md_PBE = Int epsilon_c_PBE_mu=0 / ( 1 + beta*mu**3 ) = Int eps_c_md_PBE  with beta chosen to recover the UEG large mu behaviour (JT)
 END_DOC
 implicit none
 double precision ::  r(3)
 double precision :: weight,mu
 integer :: i,istate
 double precision,allocatable  :: eps_c_md_PBE(:)
 allocate(eps_c_md_PBE(N_states))
 mu = mu_erf_dft
 Energy_c_md_PBE_mu_vector = 0.d0
  
 do i = 1, n_points_final_grid
  r(1) = final_grid_points(1,i)
  r(2) = final_grid_points(2,i)
  r(3) = final_grid_points(3,i)
  weight = final_weight_at_r_vector(i)
  call give_epsilon_c_md_PBE_mu(mu,r,eps_c_md_PBE)
  do istate = 1, N_states
   Energy_c_md_PBE_mu_vector(istate) += eps_c_md_PBE(istate) * weight
  enddo
 enddo
 END_PROVIDER



subroutine give_epsilon_c_md_PBE_mu(mu,r,eps_c_md_PBE) ! JT
  implicit none
  double precision, intent(in)  :: mu , r(3)
  double precision, intent(out) :: eps_c_md_PBE(N_states)
  double precision :: pi, e_PBE, beta
  double precision :: aos_array(ao_num), grad_aos_array(3,ao_num)
  double precision :: rho_a(N_states),rho_b(N_states)
  double precision :: grad_rho_a(3,N_states),grad_rho_b(3,N_states)
  double precision :: grad_rho_a_2(N_states),grad_rho_b_2(N_states),grad_rho_a_b(N_states)
  double precision :: rhoc,rhoo,sigmacc,sigmaco,sigmaoo,vrhoc,vrhoo,vsigmacc,vsigmaco,vsigmaoo
  double precision :: g0_UEG_mu_inf, denom
  integer :: m, istate

  pi = 4.d0 * datan(1.d0)

  eps_c_md_PBE = 0.d0
  call density_and_grad_alpha_beta_and_all_aos_and_grad_aos_at_r(r,rho_a,rho_b, grad_rho_a, grad_rho_b, aos_array, grad_aos_array)
  grad_rho_a_2 = 0.d0
  grad_rho_b_2 = 0.d0
  grad_rho_a_b = 0.d0
  do istate = 1, N_states
   do m = 1, 3
    grad_rho_a_2(istate) += grad_rho_a(m,istate)*grad_rho_a(m,istate)
    grad_rho_b_2(istate) += grad_rho_b(m,istate)*grad_rho_b(m,istate)
    grad_rho_a_b(istate) += grad_rho_a(m,istate)*grad_rho_b(m,istate)
   enddo
  enddo
  do istate = 1, N_states
   ! convertion from (alpha,beta) formalism to (closed, open) formalism
   call rho_ab_to_rho_oc(rho_a(istate),rho_b(istate),rhoo,rhoc)
   call grad_rho_ab_to_grad_rho_oc(grad_rho_a_2(istate),grad_rho_b_2(istate),grad_rho_a_b(istate),sigmaoo,sigmacc,sigmaco)
   call ec_pbe_only(0.d0,rhoc,rhoo,sigmacc,sigmaco,sigmaoo,e_PBE)

   if(mu == 0.d0) then
    eps_c_md_PBE(istate)=e_PBE
   else
!   note: the on-top pair density is (1-zeta^2) rhoc^2 g0 = 4 rhoa * rhob * g0
    denom = (-2.d0+sqrt(2d0))*sqrt(2.d0*pi) * 4.d0*rho_a(istate)*rho_b(istate)*g0_UEG_mu_inf(rho_a(istate),rho_b(istate)) 
    if (dabs(denom) > 1.d-12) then
     beta = (3.d0*e_PBE)/denom
     eps_c_md_PBE(istate)=e_PBE/(1.d0+beta*mu**3)
    else
     eps_c_md_PBE(istate)=0.d0
    endif
   endif
  enddo
 end
 

 
subroutine give_epsilon_c_md_PBE_mu_grad_input(mu,rho_a,rho_b, grad_rho_a, grad_rho_b,eps_c_md_PBE) ! EG
  implicit none
  double precision, intent(in)  :: mu , rho_a(N_states),rho_b(N_states), grad_rho_a(3,N_states),grad_rho_b(3,N_states)
  double precision, intent(out) :: eps_c_md_PBE(N_states)
  double precision :: pi, e_PBE, beta
  double precision :: aos_array(ao_num), grad_aos_array(3,ao_num)
  double precision :: grad_rho_a_2(N_states),grad_rho_b_2(N_states),grad_rho_a_b(N_states)
  double precision :: rhoc,rhoo,sigmacc,sigmaco,sigmaoo,vrhoc,vrhoo,vsigmacc,vsigmaco,vsigmaoo
  double precision :: g0_UEG_mu_inf, denom
  integer :: m, istate

  pi = 4.d0 * datan(1.d0)

  eps_c_md_PBE = 0.d0
  grad_rho_a_2 = 0.d0
  grad_rho_b_2 = 0.d0
  grad_rho_a_b = 0.d0
  do istate = 1, N_states
   do m = 1, 3
    grad_rho_a_2(istate) += grad_rho_a(m,istate)*grad_rho_a(m,istate)
    grad_rho_b_2(istate) += grad_rho_b(m,istate)*grad_rho_b(m,istate)
    grad_rho_a_b(istate) += grad_rho_a(m,istate)*grad_rho_b(m,istate)
   enddo
  enddo
  do istate = 1, N_states
   ! convertion from (alpha,beta) formalism to (closed, open) formalism
   call rho_ab_to_rho_oc(rho_a(istate),rho_b(istate),rhoo,rhoc)
   call grad_rho_ab_to_grad_rho_oc(grad_rho_a_2(istate),grad_rho_b_2(istate),grad_rho_a_b(istate),sigmaoo,sigmacc,sigmaco)
   call ec_pbe_only(0.d0,rhoc,rhoo,sigmacc,sigmaco,sigmaoo,e_PBE)

   if(mu == 0.d0) then
    eps_c_md_PBE(istate)=e_PBE
   else
!   note: the on-top pair density is (1-zeta^2) rhoc^2 g0 = 4 rhoa * rhob * g0
    denom = (-2.d0+sqrt(2d0))*sqrt(2.d0*pi) * 4.d0*rho_a(istate)*rho_b(istate)*g0_UEG_mu_inf(rho_a(istate),rho_b(istate)) 
    if (dabs(denom) > 1.d-12) then
     beta = (3.d0*e_PBE)/denom
     eps_c_md_PBE(istate)=e_PBE/(1.d0+beta*mu**3)
    else
     eps_c_md_PBE(istate)=0.d0
    endif
   endif
  enddo
 end
 

 


subroutine give_epsilon_pbe_provider(mu,i_point,eps_c_md_PBE)
  implicit none
  double precision, intent(in)  :: mu 
  double precision, intent(out) :: eps_c_md_PBE(N_states)
  integer, intent(in) :: i_point
  double precision :: two_dm, pi, e_pbe,beta,on_top_two_dm_in_r_mu_corrected_from_two_dm
  double precision :: grad_rho_a(3),grad_rho_b(3)
  double precision :: grad_rho_a_2,grad_rho_b_2,grad_rho_a_b
  double precision :: rhoc,rhoo,ec_pbe_88
  double precision :: delta,two_dm_corr,rho_a,rho_b
  double precision :: grad_rho_2,denom,g0_UEG_mu_inf
  double precision :: sigmacc,sigmaco,sigmaoo
  integer :: m, istate

  pi = 4.d0 * datan(1.d0)

  eps_c_md_PBE = 0.d0
  do istate = 1, N_states
   ! total and spin density 
   rhoc = one_e_dm_and_grad_alpha_in_r(4,i_point,istate) + one_e_dm_and_grad_beta_in_r(4,i_point,istate) 
   rhoo = one_e_dm_and_grad_alpha_in_r(4,i_point,istate) - one_e_dm_and_grad_beta_in_r(4,i_point,istate) 
   ! gradients of the effective spin density 
   grad_rho_a_2 = 0.D0
   grad_rho_b_2 = 0.D0
   grad_rho_a_b = 0.D0
   do m = 1, 3
    grad_rho_a_2 += one_e_dm_and_grad_alpha_in_r(m,i_point,istate)**2.d0
    grad_rho_b_2 += one_e_dm_and_grad_beta_in_r(m,i_point,istate) **2.d0
    grad_rho_a_b += one_e_dm_and_grad_alpha_in_r(m,i_point,istate) * one_e_dm_and_grad_beta_in_r(m,i_point,istate)
   enddo
   sigmacc = grad_rho_a_2 + grad_rho_b_2 + 2.d0 * grad_rho_a_b
   sigmaco = 0.d0
   sigmaoo = 0.d0
   rho_a = one_e_dm_and_grad_alpha_in_r(4,i_point,istate)
   rho_b = one_e_dm_and_grad_beta_in_r(4,i_point,istate)

   call ec_pbe_only(0.d0,rhoc,rhoo,sigmacc,sigmaco,sigmaoo,e_PBE)
   if(e_PBE.gt.0.d0)then
    print*,'PBE gt 0 with regular dens'
   endif
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
