subroutine give_epsilon_c_md_n_and_on_top_PBE_mu_corrected_from_two_dm(mu,r,two_dm,eps_c_md_on_top_PBE)
  implicit none
  BEGIN_DOC
! enter with "r(3)", and "two_dm(N_states)" which is the on-top pair density at "r" for each states
!
! you get out with the energy density defined in J. Chem. Phys.150, 084103 (2019); doi: 10.1063/1.508263
!
! by Eq. (26), which includes the correction of the on-top pair density of Eq. (29). 
!
! !!!!! BUT !!!!! with the effective spin polarization defined in Phys. Rev. A 51, 4531 (1995).
!
! !!!!! Such an effective spin polarization depends only on the on-top and the total density.
!
! !!!!! Therefore the effective spin polarization is S_z independent
  END_DOC
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
   ! correction of the on-top pair density according to Eq. (29) 
   two_dm_corr = on_top_two_dm_in_r_mu_corrected_from_two_dm(mu,istate,two_dm)

   ! effective spin polarization
   if(rhoc*rhoc - 4.d0 * two_dm_corr .ge.0.d0)then
    rhoo =  dsqrt(rhoc*rhoc - 4.d0 * two_dm_corr) ! effective spin polarization from the on-top pair density and total density
   else 
    if(dabs(rhoc*rhoc).gt.1.d-10.and.dabs(two_dm_corr).gt.1.d-10)then
     print*,'Imaginary effective spin polarization !'
     print*,'r = '
     print*,r
     print*,rhoc*rhoc , 4.d0 * two_dm_corr
    endif
    rhoo = 0.d0 
   endif
   ! PBE correlation energy using the effective spin polarization as input 
   call ec_pbe_only(0.d0,rhoc,rhoo,sigmacc,sigmaco,sigmaoo,e_PBE(istate))
   ! quantity of Eq. (27)
   beta(istate) = (3.d0*e_PBE(istate))/( (-2.d0+sqrt(2d0))*sqrt(2.d0*pi)*2.d0*two_dm_corr )
   ! quantity of Eq. (26)
   eps_c_md_on_top_PBE(istate)=e_PBE(istate)/(1.d0+beta(istate)*mu**3)
  enddo
 end
 
subroutine give_epsilon_c_md_n_and_on_top_LYP_mu_corrected_from_two_dm(mu,r,two_dm,eps_c_md_on_top_LYP)
  implicit none
  BEGIN_DOC
! enter with "r(3)", and "two_dm(N_states)" which is the on-top pair density at "r" for each states
!
! you get out with the energy density defined in J. Chem. Phys.150, 084103 (2019); doi: 10.1063/1.508263
!
! by Eq. (26), which includes the correction of the on-top pair density of Eq. (29). 
!
! !!!!! BUT !!!!! with the effective spin polarization defined in Phys. Rev. A 51, 4531 (1995).
!
! !!!!! Such an effective spin polarization depends only on the on-top and the total density.
!
! !!!!! Therefore the effective spin polarization is S_z independent
  END_DOC
  double precision, intent(in)  :: mu , r(3), two_dm(N_states)
  double precision, intent(out) :: eps_c_md_on_top_LYP(N_states)
  double precision :: two_dm_in_r, pi, e_lyp(N_states),beta(N_states),on_top_two_dm_in_r_mu_corrected_from_two_dm
  double precision :: rho_a(N_states),rho_b(N_states),grad_rho_a_2(N_states),grad_rho_b_2(N_states),grad_rho_2(N_states),rho(N_states)
  double precision :: ec_lyp_88
  integer :: m, istate

  pi = 4.d0 * datan(1.d0)

  eps_c_md_on_top_LYP = 0.d0
  call give_all_stuffs_in_r_for_lyp_88(r,rho,rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_2)
  do istate = 1, N_states
   ! convertion from (alpha,beta) formalism to (closed, open) formalism
  !call rho_ab_to_rho_oc(rho_a(istate),rho_b(istate),rhoo,rhoc)
  !call grad_rho_ab_to_grad_rho_oc(grad_rho_a_2(istate),grad_rho_b_2(istate),grad_rho_a_b(istate),sigmaoo,sigmacc,sigmaco)
  !sigmaco = 0.d0
  !sigmaoo = 0.d0
   double precision :: delta,two_dm_corr,rhoo_2
  !rhoo_2 = rhoo
   two_dm_corr = on_top_two_dm_in_r_mu_corrected_from_two_dm(mu,istate,two_dm)
  !if(rhoc*rhoc - 4.d0 * two_dm_corr .ge.0.d0)then
  ! rhoo =  dsqrt(rhoc*rhoc - 4.d0 * two_dm_corr) ! effective spin polarization from the on-top pair density and total density
  !else 
  ! if(dabs(rhoc*rhoc).gt.1.d-10.and.dabs(two_dm_corr).gt.1.d-10)then
  !  print*,'Imaginary effective spin polarization !'
  !  print*,'r = '
  !  print*,r
  !  print*,rhoc*rhoc , 4.d0 * two_dm_corr
  ! endif
  !endif
   e_lyp(istate) = ec_lyp_88(rho(istate),rho_a(istate),rho_b(istate),grad_rho_a_2(istate),grad_rho_b_2(istate),grad_rho_2(istate))
   beta(istate) = (3.d0*e_LYP(istate))/( (-2.d0+sqrt(2d0))*sqrt(2.d0*pi)*2.d0*two_dm_corr )
   eps_c_md_on_top_LYP(istate)=e_LYP(istate)/(1.d0+beta(istate)*mu**3)
  enddo
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
 
