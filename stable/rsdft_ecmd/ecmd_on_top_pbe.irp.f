 BEGIN_PROVIDER [double precision, Energy_c_md_on_top_PBE_mu_UEG_vector, (N_states)]
 BEGIN_DOC
  ! Give the Ec_md energy with a good large mu behaviour in function of the on top pair density with mu correction based on the on-top of the UEG coupled to the PBE correlation energy at mu=0
  ! Ec_md_on_top_PBE = Int epsilon_c_PBE_mu=0 / ( 1 + beta*mu**3 ) = Int eps_c_md_on_top_PBE  with beta chosen to recover the good large mu behaviour of the Energy_c_md_on_top functional
 END_DOC
 implicit none
 double precision ::  r(3)
 double precision :: weight,mu
 integer :: i,istate
 double precision,allocatable  :: eps_c_md_on_top_PBE(:),two_dm(:)
 allocate(eps_c_md_on_top_PBE(N_states),two_dm(N_states))
 mu = mu_erf_dft
 Energy_c_md_on_top_PBE_mu_UEG_vector = 0.d0
  
 do i = 1, n_points_final_grid
  r(1) = final_grid_points(1,i)
  r(2) = final_grid_points(2,i)
  r(3) = final_grid_points(3,i)
  weight = final_weight_at_r_vector(i)
  two_dm(:) = on_top_of_r_vector(i,:)
  call give_epsilon_c_md_on_top_PBE_mu_corrected_UEG_from_two_dm(mu,r,two_dm,eps_c_md_on_top_PBE)
  do istate = 1, N_states
   Energy_c_md_on_top_PBE_mu_UEG_vector(istate) += eps_c_md_on_top_PBE(istate) * weight
  enddo
 enddo
 END_PROVIDER

 BEGIN_PROVIDER [double precision, Energy_c_md_on_top_PBE_mu_vector, (N_states)]
 BEGIN_DOC
  ! Give the Ec_md energy with a good large mu behaviour in function of the on top pair density with mu correction based on the on-top of the UEG coupled to the PBE correlation energy at mu=0
  ! Ec_md_on_top_PBE = Int epsilon_c_PBE_mu=0 / ( 1 + beta*mu**3 ) = Int eps_c_md_on_top_PBE  with beta chosen to recover the good large mu behaviour of the Energy_c_md_on_top functional
 END_DOC
 implicit none
 double precision ::  r(3)
 double precision :: weight,mu
 integer :: i,istate
 double precision,allocatable  :: eps_c_md_on_top_PBE(:),two_dm(:)
 allocate(eps_c_md_on_top_PBE(N_states),two_dm(N_states))
 mu = mu_erf_dft
 Energy_c_md_on_top_PBE_mu_vector = 0.d0
  
 do i = 1, n_points_final_grid
  r(1) = final_grid_points(1,i)
  r(2) = final_grid_points(2,i)
  r(3) = final_grid_points(3,i)
  weight = final_weight_at_r_vector(i)
  two_dm(:) = on_top_of_r_vector(i,:)
  call give_epsilon_c_md_on_top_PBE_mu_corrected_from_two_dm(mu,r,two_dm,eps_c_md_on_top_PBE)
  do istate = 1, N_states
   Energy_c_md_on_top_PBE_mu_vector(istate) += eps_c_md_on_top_PBE(istate) * weight
  enddo
 enddo
 END_PROVIDER



subroutine give_epsilon_c_md_on_top_PBE_mu_corrected_from_two_dm(mu,r,two_dm,eps_c_md_on_top_PBE)
  implicit none
  double precision, intent(in)  :: mu , r(3), two_dm(N_states)
  double precision, intent(out) :: eps_c_md_on_top_PBE(N_states)
  double precision :: two_dm_in_r, pi, e_pbe(N_states),beta(N_states),on_top_two_dm_in_r_mu_corrected_from_two_dm
  double precision :: aos_array(ao_num), grad_aos_array(3,ao_num)
  double precision :: rho_a(N_states),rho_b(N_states)
  double precision :: grad_rho_a(3,N_states),grad_rho_b(3,N_states)
  double precision :: grad_rho_a_2(N_states),grad_rho_b_2(N_states),grad_rho_a_b(N_states)
  double precision :: rhoc,rhoo,sigmacc,sigmaco,sigmaoo,vrhoc,vrhoo,vsigmacc,vsigmaco,vsigmaoo
  integer :: m, istate

  pi = 4.d0 * datan(1.d0)

  eps_c_md_on_top_PBE = 0.d0
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
   call ec_pbe_only(0.d0,rhoc,rhoo,sigmacc,sigmaco,sigmaoo,e_PBE(istate))
   beta(istate) = (3.d0*e_PBE(istate))/( (-2.d0+sqrt(2d0))*sqrt(2.d0*pi)*2.d0*on_top_two_dm_in_r_mu_corrected_from_two_dm(mu,r,istate,two_dm) )
   eps_c_md_on_top_PBE(istate)=e_PBE(istate)/(1.d0+beta(istate)*mu**3)
  enddo
 end


subroutine give_epsilon_c_md_n_and_on_top_PBE_mu_corrected_from_two_dm(mu,r,two_dm,eps_c_md_on_top_PBE)
  implicit none
  double precision, intent(in)  :: mu , r(3), two_dm(N_states)
  double precision, intent(out) :: eps_c_md_on_top_PBE(N_states)
  double precision :: two_dm_in_r, pi, e_pbe(N_states),beta(N_states),on_top_two_dm_in_r_mu_corrected_from_two_dm
  double precision :: aos_array(ao_num), grad_aos_array(3,ao_num)
  double precision :: rho_a(N_states),rho_b(N_states)
  double precision :: grad_rho_a(3,N_states),grad_rho_b(3,N_states)
  double precision :: grad_rho_a_2(N_states),grad_rho_b_2(N_states),grad_rho_a_b(N_states)
  double precision :: rhoc,rhoo,sigmacc,sigmaco,sigmaoo,vrhoc,vrhoo,vsigmacc,vsigmaco,vsigmaoo
  integer :: m, istate

  pi = 4.d0 * datan(1.d0)

  eps_c_md_on_top_PBE = 0.d0
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
   sigmaco = 0.d0
   sigmaoo = 0.d0
   double precision :: delta,two_dm_corr,rhoo_2
   rhoo_2 = rhoo
   two_dm_corr = on_top_two_dm_in_r_mu_corrected_from_two_dm(mu,r,istate,two_dm)
   if(rhoc*rhoc - 4.d0 * two_dm_corr .ge.0.d0)then
    rhoo =  dsqrt(rhoc*rhoc - 4.d0 * two_dm_corr) ! effective spin polarization from the on-top pair density and total density
   else 
    if(dabs(rhoc*rhoc).gt.1.d-10.and.dabs(two_dm_corr).gt.1.d-10)then
     print*,'Imaginary effective spin polarization !'
     print*,'r = '
     print*,r
     print*,rhoc*rhoc , 4.d0 * two_dm_corr
    endif
 
   endif
   call ec_pbe_only(0.d0,rhoc,rhoo,sigmacc,sigmaco,sigmaoo,e_PBE(istate))
   beta(istate) = (3.d0*e_PBE(istate))/( (-2.d0+sqrt(2d0))*sqrt(2.d0*pi)*2.d0*two_dm_corr )
   eps_c_md_on_top_PBE(istate)=e_PBE(istate)/(1.d0+beta(istate)*mu**3)
  enddo
 end
 
 
 double precision function on_top_two_dm_in_r_mu_corrected_from_two_dm(mu,r,istate,two_dm)
 implicit none
 double precision, intent(in) :: mu,r(3),two_dm(N_states)
 integer, intent(in) :: istate
 double precision :: pi
 pi = 4.d0 * datan(1.d0)
 on_top_two_dm_in_r_mu_corrected_from_two_dm = two_dm(istate)  / ( 1.d0 + 2.d0/(dsqrt(pi)*mu) )
 on_top_two_dm_in_r_mu_corrected_from_two_dm = max(on_top_two_dm_in_r_mu_corrected_from_two_dm ,1.d-15)
 end

 double precision function on_top_two_dm_in_r_mu_corrected_UEG_from_two_dm(mu,r,istate,two_dm)
 implicit none
 double precision, intent(in) :: mu,r(3),two_dm(N_states)
 integer, intent(in) :: istate
 
 double precision :: correction_to_on_top_from_UEG
 on_top_two_dm_in_r_mu_corrected_UEG_from_two_dm = two_dm(istate) * correction_to_on_top_from_UEG(mu,r,istate)
 on_top_two_dm_in_r_mu_corrected_UEG_from_two_dm = max(on_top_two_dm_in_r_mu_corrected_UEG_from_two_dm,1.d-15)
 end



subroutine give_epsilon_c_md_on_top_PBE_mu_corrected_UEG_from_two_dm(mu,r,two_dm,eps_c_md_on_top_PBE)
  implicit none
  double precision, intent(in)  :: mu , r(3), two_dm(N_states)
  double precision, intent(out) :: eps_c_md_on_top_PBE(N_states)
  double precision :: two_dm_in_r, pi, e_pbe(N_states),beta(N_states),on_top_two_dm_in_r_mu_corrected_UEG_from_two_dm
  double precision :: aos_array(ao_num), grad_aos_array(3,ao_num)
  double precision :: rho_a(N_states),rho_b(N_states)
  double precision :: grad_rho_a(3,N_states),grad_rho_b(3,N_states)
  double precision :: grad_rho_a_2(N_states),grad_rho_b_2(N_states),grad_rho_a_b(N_states)
  double precision :: rhoc,rhoo,sigmacc,sigmaco,sigmaoo,vrhoc,vrhoo,vsigmacc,vsigmaco,vsigmaoo
  integer :: m, istate
  pi = 4d0 * datan(1d0)

  eps_c_md_on_top_PBE = 0d0
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
   call ec_pbe_only(0.d0,rhoc,rhoo,sigmacc,sigmaco,sigmaoo,e_PBE(istate))
   beta(istate) = (3d0*e_PBE(istate))/( (-2d0+sqrt(2d0))*sqrt(2d0*pi)*2d0*on_top_two_dm_in_r_mu_corrected_UEG_from_two_dm(mu,r,istate,two_dm) )
   eps_c_md_on_top_PBE(istate)=e_PBE(istate)/(1d0+beta(istate)*mu**3d0)
  enddo
 end
 
