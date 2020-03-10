subroutine give_d_Ec_supbeot_d_rho_and_d_gradrho(rho,on_top,grad_rho2,mu,d_ec_pbeueg_d_rho,d_ec_pbeueg_d_grad_rho2)
 BEGIN_DOC
  !Derivative of the epsilon_{c,md}^{sr,SU-PBE-ot} function (eq 14a of ref J.Chem.Lett 2019, 2931-2937) with respect to total density
 END_DOC
 implicit none
 double precision, intent(in)  :: rho,on_top,grad_rho2,mu
 double precision, intent(out) :: d_ec_pbeueg_d_rho,d_ec_pbeueg_d_grad_rho2
 double precision :: rhoo,sigmacc,sigmaco,sigmaoo,e_PBE,vrhoc,vrhoo,vsigmacc,vsigmaco,vsigmaoo
 double precision :: denom,beta,thr,d_beta_d_rho,d_beta_d_grad_rho2,one_beta_mu3,mu3,c_beta,pi

 thr = 1.d-12
 pi=dacos(-1.d0) 
 c_beta=3.d0/(2.d0*sqrt(pi)*(1.d0-sqrt(2.d0)))
 mu3 = mu * mu * mu

  ! the functional is SPIN-UNPOLAIRIZED !!!!
  rhoo = 0.d0
  sigmaoo = 0.d0
  sigmaco = 0.d0

  call ec_pbe_sr(0.d0,rho,rhoo,grad_rho2,sigmaco,sigmaoo,e_PBE,vrhoc,vrhoo,vsigmacc,vsigmaco,vsigmaoo)


  if(dabs(on_top).gt.thr)then
   beta = c_beta * e_PBE / on_top
   d_beta_d_rho = c_beta / on_top * vrhoc 
   d_beta_d_grad_rho2 = c_beta / on_top * vsigmacc
  else
   beta = c_beta * e_PBE / thr 
   d_beta_d_rho = c_beta / thr * vrhoc 
   d_beta_d_grad_rho2 = c_beta / thr * vsigmacc
  endif

  one_beta_mu3 = 1.d0 + beta * mu3
  if(dabs(one_beta_mu3*one_beta_mu3).gt.thr)then
   denom = 1.d0/(one_beta_mu3*one_beta_mu3)
  else
   denom = 1.d0/thr
  endif
  d_ec_pbeueg_d_rho  = vrhoc * one_beta_mu3 - mu3 * e_PBE * d_beta_d_rho 
  d_ec_pbeueg_d_rho *= denom 

  d_ec_pbeueg_d_grad_rho2  = vsigmacc * one_beta_mu3 - mu3 * e_PBE * d_beta_d_grad_rho2 
  d_ec_pbeueg_d_grad_rho2 *= denom 


 end


! BEGIN_PROVIDER[double precision, aos_sr_vc_alpha_pbe_ueg_w  , (ao_num,n_points_final_grid,N_states)]
!&BEGIN_PROVIDER[double precision, aos_sr_vc_beta_pbe_ueg_w   , (ao_num,n_points_final_grid,N_states)]
!&BEGIN_PROVIDER[double precision, aos_dsr_vc_alpha_pbe_ueg_w  , (ao_num,n_points_final_grid,N_states)]
!&BEGIN_PROVIDER[double precision, aos_dsr_vc_beta_pbe_ueg_w   ,  (ao_num,n_points_final_grid,N_states)]
! implicit none
! BEGIN_DOC
!  !! Intermediate quantities for the calculation of the vc potentials for the PBE UEG functional
! END_DOC
! integer :: istate,i,j,m
! double precision :: r(3)
! double precision :: mu,weight
! double precision :: d_ec_pbeueg_d_ec_pbe(N_states),d_ec_pbeueg_d_grad_n_alpha(3,N_states),d_ec_pbeueg_d_grad_n_beta(3,N_states)
! double precision :: d_ec_pbeueg_rhoa(N_states),d_ec_pbeueg_rhob(N_states)
!
! aos_dsr_vc_alpha_pbe_ueg_w= 0.d0
! aos_dsr_vc_beta_pbe_ueg_w = 0.d0
! do istate = 1, N_states
!  do i = 1, n_points_final_grid
!   r(1) = final_grid_points(1,i)
!   r(2) = final_grid_points(2,i)
!   r(3) = final_grid_points(3,i)
!   mu = mu_of_r_vector(i)
!   weight = final_weight_at_r_vector(i)
!
!   call give_d_Ec_pbeueg_rho(r,mu,d_ec_pbeueg_rhoa,d_ec_pbeueg_rhob)
!   call give_d_Ec_pbeueg_d_grad_n(r,mu,d_ec_pbeueg_d_ec_pbe,d_ec_pbeueg_d_grad_n_alpha,d_ec_pbeueg_d_grad_n_beta)
!
!   do j = 1, ao_num
!    aos_sr_vc_alpha_pbe_ueg_w(j,i,istate) = d_ec_pbeueg_rhoa(istate) * weight * aos_in_r_array(j,i)
!    aos_sr_vc_beta_pbe_ueg_w (j,i,istate) = d_ec_pbeueg_rhob(istate) * weight * aos_in_r_array(j,i)
!   enddo
!   do j = 1, ao_num
!    do m = 1,3
!     aos_dsr_vc_alpha_pbe_ueg_w(j,i,istate) += weight * d_ec_pbeueg_d_grad_n_alpha(m,istate) * aos_grad_in_r_array_transp_xyz(m,j,i)
!     aos_dsr_vc_beta_pbe_ueg_w (j,i,istate) += weight * d_ec_pbeueg_d_grad_n_beta(m,istate)  * aos_grad_in_r_array_transp_xyz(m,j,i)
!    enddo
!
!   enddo
!  enddo
! enddo
!
! END_PROVIDER
