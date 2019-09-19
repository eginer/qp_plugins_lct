program write_integrals_for_dft_ecmd_lda
 implicit none
 read_wf = .true.
 touch read_wf
 no_core_density = "no_core_dm"
 mu_of_r_potential = "hf_valence_coallescence"
 touch mu_of_r_potential
 touch no_core_density
 call truc_a_faire_lda
!call truc_a_faire_sr_pbe
 call truc_a_faire_pbeueg    
end


 subroutine truc_a_faire_lda 
 implicit none
 integer :: i,istate
 double precision :: r(3)
 double precision :: rhoa(n_states),rhob(n_states),mos_array(mo_num)
 double precision :: d_total_deltarho_rhoa,d_total_deltarho_rhob 
 double precision :: accu,accu_2,accu_3,mu 
 print*,'\\\\\\\\\\\\\\\\\'
 print*,' '


!!VErrrrifffff
!logical :: dospin
 double precision :: e_c,sr_vc_a,sr_vc_b
!dospin = .true. ! JT dospin have to be set to true for open shell
!poto =0.d0

 double precision :: threshold
 double precision :: homo
 threshold = 1d-15

 do istate = 1, n_states

  accu = 0.d0 
  accu_2 = 0.d0 
  accu_3 = 0.d0 
  do i = 1, n_points_final_grid
   r(1) = final_grid_points(1,i)
   r(2) = final_grid_points(2,i)
   r(3) = final_grid_points(3,i)
   mu = mu_of_r_vector(i)

   rhoa(istate) = one_e_dm_alpha_at_r(i,istate)
   rhob(istate) = one_e_dm_beta_at_r(i,istate)
   call ec_lda_sr(mu,rhoa(istate),rhob(istate),e_c,sr_vc_a,sr_vc_b)
   call give_all_mos_at_r(r,mos_array)
   homo = mos_in_r_array_transp(i,elec_alpha_num)
   if(dabs(rhoa(istate)+rhob(istate)).lt.threshold) cycle 
   accu += (d_total_deltarho_rhoa(rhoa(istate),rhob(istate),mu)+d_total_deltarho_rhob(rhoa(istate),rhob(istate),mu)+sr_vc_a+sr_vc_b)*homo**2.d0*final_weight_at_r_vector(i)
   accu_2 += (d_total_deltarho_rhoa(rhoa(istate),rhob(istate),mu)+sr_vc_a)*homo**2.d0*final_weight_at_r_vector(i)
   accu_3 += sr_vc_a*homo**2.d0*final_weight_at_r_vector(i)

  !print*,d_total_deltarho_rhoa(rhoa(istate),rhob(istate),mu),homo**2.d0 
  !!VErrrrifffff
  !call ESRC_MD_LDAERF (mu,rhoa(istate),rhob(istate),dospin,ec(istate))
  !if(isnan(ec(istate)))then
  ! print*,'ec is nan'
  ! stop
  !endif
  !poto(istate) += ec(istate)*final_weight_at_r_vector(i)

 
  enddo
  print*,'< HOMO | dE_c,md/dn_alpha+dE_c/dn_beta | HOMO >          =',accu
  print*,'< HOMO | dE_c,md/dn_alpha | HOMO >                       =',accu_2
  print*,'< HOMO | dE_c/dn_alpha | HOMO >                          =',accu_3
 !print*,'Vrai Ecmd LDA                =',Energy_c_md_LDA_mu_of_r(istate)
 !print*,'Ecmd LDA bart                =',poto(istate)
 enddo
end


 subroutine truc_a_faire_pbeueg
 implicit none
 integer :: istate
 print*,'\\\\\\\\\\\\\\\\\'

 do istate = 1, n_states
  print*,'< HOMO | dE_c,md/dn_alpha+dE_c/dn_beta PBE UEG | HOMO >          =',potential_c_alpha_mo_sr_pbe_ueg(elec_alpha_num,elec_alpha_num,istate) + potential_c_beta_mo_sr_pbe_ueg(elec_alpha_num,elec_alpha_num,istate)
  print*,'< HOMO | dE_c,md/dn_alpha PBE UEG | HOMO >                       =',potential_c_alpha_mo_sr_pbe_ueg(elec_alpha_num,elec_alpha_num,istate)
 enddo
end


!subroutine truc_a_faire_sr_pbe 
!implicit none
!integer :: i,m,istate
!double precision :: r(3),accu_3,mu 
!double precision :: rhoa(N_states),rhob(N_states),grad_rho_a(3,N_states),grad_rho_b(3,N_states)
!double precision :: aos_array(ao_num),grad_aos_array(3,ao_num),mos_array(mo_num)
!double precision :: grad_rho_a_2(N_states),grad_rho_b_2(N_states),grad_rho_a_b(N_states)
!double precision :: rhoc,rhoo,sigmacc,sigmaco,sigmaoo,e_PBE
!double precision :: vsigmacc,vsigmaco,vsigmaoo,vrhoo,vrhoc,vc_rho_a,vc_rho_b
!double precision :: vc_grad_rho_a_2(N_states),vc_grad_rho_b_2(N_states),vc_grad_rho_a_b(N_states)

!double precision :: threshold
!double precision :: homo
!threshold = 1d-15

!grad_rho_a_2 = 0.d0
!grad_rho_b_2 = 0.d0
!grad_rho_a_b = 0.d0

!do istate = 1, n_states

! accu_3 = 0.d0 
! do i = 1, n_points_final_grid
!  r(1) = final_grid_points(1,i)
!  r(2) = final_grid_points(2,i)
!  r(3) = final_grid_points(3,i)
!  mu = mu_of_r_vector(i)

!  rhoa(istate) = one_e_dm_alpha_at_r(i,istate)
!  rhob(istate) = one_e_dm_beta_at_r(i,istate)
!  call give_all_mos_at_r(r,mos_array)
!  call density_and_grad_alpha_beta_and_all_aos_and_grad_aos_at_r(r,rhoa,rhob, grad_rho_a, grad_rho_b, aos_array, grad_aos_array)
!  do m = 1, 3
!   grad_rho_a_2(istate) += grad_rho_a(m,istate)*grad_rho_a(m,istate)
!   grad_rho_b_2(istate) += grad_rho_b(m,istate)*grad_rho_b(m,istate)
!   grad_rho_a_b(istate) += grad_rho_a(m,istate)*grad_rho_b(m,istate)
!  enddo

! ! convertion from (alpha,beta) formalism to (closed, open) formalism

!  call rho_ab_to_rho_oc(rhoa(istate),rhob(istate),rhoo,rhoc)
!  call grad_rho_ab_to_grad_rho_oc(grad_rho_a_2(istate),grad_rho_b_2(istate),grad_rho_a_b(istate),sigmaoo,sigmacc,sigmaco)
!  call ec_pbe_sr(mu,rhoc,rhoo,sigmacc,sigmaco,sigmaoo,ec,vrhoc,vrhoo,vsigmacc,vsigmaco,vsigmaoo)

!  call v_grad_rho_oc_to_v_grad_rho_ab(vsigmaoo,vsigmacc,vsigmaco,vc_grad_rho_a_2(istate),vc_grad_rho_b_2(istate),vc_grad_rho_a_b(istate))

!  homo = mos_in_r_array_transp(i,elec_alpha_num)
!  if(dabs(rhoa(istate)+rhob(istate)).lt.threshold) cycle 
!  accu_3 += sr_vc_a*homo**2.d0*final_weight_at_r_vector(i)

!
! enddo
! print*,'< HOMO | dE_c/dn_alpha PBE | HOMO >                          =',accu_3
!enddo
!nd
