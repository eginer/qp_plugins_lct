
 double precision function f_pbeueg(rs,xi)
 BEGIN_DOC
  !f_pbeueg(rs,xi)= (1-xi^2)g0f(rs) 
 END_DOC
 implicit none        
 double precision rs,xi,g0f
 f_pbeueg=(1d0-xi**2)*g0f(rs)
 return
 end

 double precision function f_pbeueg2(rhoa,rhob)
 BEGIN_DOC
  !f_pbeueg(rs,xi)= (1-xi^2)g0_UEG_mu_inf(rhoa,rhob), same function than f_pbeueg but with a slightly different g0 
 END_DOC
 implicit none        
 double precision rhoa,rhob,xi,g0_UEG_mu_inf
 xi= (rhoa-rhob)/(rhoa+rhob)
 f_pbeueg2=(1d0-xi**2)*g0_UEG_mu_inf(rhoa,rhob)
 return
 end

 double precision function d_f_pbeueg(rs,xi)
 BEGIN_DOC
  !Derivative of (1-xi^2)g0f(rs) with respect to rs
 END_DOC
 implicit none        
 double precision rs,xi,g0d
 d_f_pbeueg=(1d0-xi**2)*g0d(rs)
 return
 end

 double precision function d_xi_f_pbeueg(rs,xi)
 BEGIN_DOC
  !Derivative of (1-xi^2)g0f(rs) with respect to xi
 END_DOC
 implicit none
 double precision rs,xi,g0f
 d_xi_f_pbeueg = -2d0*xi*g0f(rs)
 return
 end


 double precision function d_f_pbeueg_rhoa(rhoa,rhob)
 BEGIN_DOC
  ! derivative of (1-xi^2)g0f(rs) with respect to rho_alpha
 END_DOC
 implicit none
 double precision d_wignerseitz_radius,wignerseitz_radius,d_xi_rhoa,xi,d_f_pbeueg,d_xi_f_pbeueg,rhot,rs,drs,dxi_dna,rhoa,rhob
 rhot=rhoa+rhob
 rs = wignerseitz_radius(rhot)
 drs= d_wignerseitz_radius(rhot)
 xi= (rhoa-rhob)/(rhoa+rhob)
 dxi_dna=d_xi_rhoa(rhoa,rhob)
 d_f_pbeueg_rhoa=d_f_pbeueg(rs,xi)*drs+d_xi_f_pbeueg(rs,xi)*dxi_dna
 return
 end


 double precision function d_f_pbeueg_rhob(rhoa,rhob)
 BEGIN_DOC
  ! derivative of (1-xi^2)g0f(rs) with respect to rho_beta
 END_DOC
 implicit none
 double precision d_wignerseitz_radius,wignerseitz_radius,d_xi_rhob,xi,d_f_pbeueg,d_xi_f_pbeueg,rhot,rs,drs,dxi_dnb,rhoa,rhob
 rhot=rhoa+rhob
 rs = wignerseitz_radius(rhot)
 drs= d_wignerseitz_radius(rhot)
 xi= (rhoa-rhob)/(rhoa+rhob)
 dxi_dnb=d_xi_rhob(rhoa,rhob)
 d_f_pbeueg_rhob=d_f_pbeueg(rs,xi)*drs+d_xi_f_pbeueg(rs,xi)*dxi_dnb
 return
 end



subroutine give_Ec_pbeueg_test(mu,rhoa,rhob,grad_rho_a,grad_rho_b,epsilon_c_pbeueg_bart1,epsilon_c_pbeueg_bart2,beta_test,beta_test2)
 BEGIN_DOC
  ! Epsilon_c,md^srPBE and Beta function (eq 14a and 14b of ref J.Chem.Lett 2019, 2931-2937) test version 
 END_DOC
 implicit none
 double precision, intent(in)  :: rhoa(N_states),rhob(N_states),grad_rho_a(3,N_states),grad_rho_b(3,N_states),mu
 double precision, intent(out) :: epsilon_c_pbeueg_bart1(N_states),epsilon_c_pbeueg_bart2(N_states),beta_test(N_states),beta_test2(N_states)
 integer :: istate,m 
 double precision :: grad_rho_a_2(N_states),grad_rho_b_2(N_states),grad_rho_a_b(N_states),e_PBE(N_states)
 double precision :: wignerseitz_radius,xi,rhot,rs,f_pbeueg,f_pbeueg2,c_beta,pi,rhoc,rhoo,sigmacc,sigmaco,sigmaoo

 pi=dacos(-1.d0) 
 c_beta=3.d0/(2.d0*sqrt(pi)*(1.d0-sqrt(2.d0)))

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
  !convertion from (alpha,beta) formalism to (closed, open) formalism
  rhot=rhoa(istate)+rhob(istate)
  rs = wignerseitz_radius(rhot)
  xi= (rhoa(istate)-rhob(istate))/(rhot)

  call rho_ab_to_rho_oc(rhoa(istate),rhob(istate),rhoo,rhoc)
  call grad_rho_ab_to_grad_rho_oc(grad_rho_a_2(istate),grad_rho_b_2(istate),grad_rho_a_b(istate),sigmaoo,sigmacc,sigmaco)
  call ec_pbe_only(0.d0,rhoc,rhoo,sigmacc,sigmaco,sigmaoo,e_PBE(istate))

  double precision :: denom1,denom2,g0_UEG_mu_inf,g0f

  if(mu == 0.d0) then
   epsilon_c_pbeueg_bart1(istate)=e_PBE(istate)
   epsilon_c_pbeueg_bart2(istate)=e_PBE(istate)
  else

   denom1 = (-2.d0+dsqrt(2d0))*sqrt(2.d0*pi) * 4.d0*rhoa(istate)*rhob(istate)*g0_UEG_mu_inf(rhoa(istate),rhob(istate))

  !denom2 = rhot**2.d0 * f_pbeueg2(rhoa,rhob) / c_beta
   denom2 = rhot**2.d0 * f_pbeueg(rs,xi) / c_beta

   if (dabs(denom1) > 1.d-10) then
    beta_test(istate) = (3.d0*e_PBE(istate))/denom1
    epsilon_c_pbeueg_bart1(istate)=e_PBE(istate)/(1.d0+beta_test(istate)*mu**3)
   else
    beta_test(istate)=0.d0
    epsilon_c_pbeueg_bart1(istate)=0.d0
   endif

   if (dabs(denom2) > 1.d-10) then
    beta_test2(istate) = e_PBE(istate)/denom2
    epsilon_c_pbeueg_bart2(istate)=e_PBE(istate)/(1.d0+beta_test2(istate)*mu**3)
   else
    beta_test2(istate)=0.d0
    epsilon_c_pbeueg_bart2(istate)=0.d0
   endif

  endif
 enddo
 end


subroutine d_beta_pbeueg_rho(rhoa,rhob,grad_rho_a,grad_rho_b,d_beta_rhoa,d_beta_rhob)
 BEGIN_DOC
  !Derivative of the beta function (eq 14b of ref J.Chem.Lett 2019, 2931-2937) with respect to rhoa and rhob 
 END_DOC
 implicit none
 double precision, intent(in)  :: rhoa(N_states),rhob(N_states),grad_rho_a(3,N_states),grad_rho_b(3,N_states)
 double precision, intent(out) :: d_beta_rhoa(N_states),d_beta_rhob(N_states)
 integer :: istate,m 
 double precision :: grad_rho_a_2(N_states),grad_rho_b_2(N_states),grad_rho_a_b(N_states)
 double precision :: wignerseitz_radius,xi,rhot,rs,f_pbeueg,c_beta,pi,rhoc,rhoo,sigmacc,sigmaco,sigmaoo,d_f_pbeueg_rhoa,d_f_pbeueg_rhob,e_PBE
 double precision :: vsigmacc,vsigmaco,vsigmaoo,vrhoo,vrhoc,vc_rho_a,vc_rho_b
 pi=dacos(-1.d0) 
 c_beta=3.d0/(2.d0*sqrt(pi)*(1.d0-sqrt(2.d0)))

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
  rhot=rhoa(istate)+rhob(istate)
  rs = wignerseitz_radius(rhot)
  xi= (rhoa(istate)-rhob(istate))/(rhot)

  call rho_ab_to_rho_oc(rhoa(istate),rhob(istate),rhoo,rhoc)
  call grad_rho_ab_to_grad_rho_oc(grad_rho_a_2(istate),grad_rho_b_2(istate),grad_rho_a_b(istate),sigmaoo,sigmacc,sigmaco)
  call ec_pbe_sr(0.d0,rhoc,rhoo,sigmacc,sigmaco,sigmaoo,e_PBE,vrhoc,vrhoo,vsigmacc,vsigmaco,vsigmaoo)
  call v_rho_oc_to_v_rho_ab(vrhoo,vrhoc,vc_rho_a,vc_rho_b)

  if (dabs((rhot**2.d0*f_pbeueg(rs,xi))**2.d0) > 1.d-12) then
   d_beta_rhoa(istate) = (vc_rho_a*rhot**2.d0*f_pbeueg(rs,xi)-e_PBE*(2.d0*rhot*f_pbeueg(rs,xi)+rhot**2.d0*d_f_pbeueg_rhoa(rhoa(istate),rhob(istate))))/((rhot**2.d0*f_pbeueg(rs,xi))**2.d0)
   d_beta_rhob(istate) = (vc_rho_b*rhot**2.d0*f_pbeueg(rs,xi)-e_PBE*(2.d0*rhot*f_pbeueg(rs,xi)+rhot**2.d0*d_f_pbeueg_rhob(rhoa(istate),rhob(istate))))/((rhot**2.d0*f_pbeueg(rs,xi))**2.d0)
  else
   d_beta_rhoa(istate) = 0.d0 
   d_beta_rhob(istate) = 0.d0 
  endif
 enddo
 d_beta_rhoa=c_beta*d_beta_rhoa
 d_beta_rhob=c_beta*d_beta_rhob
 end


subroutine give_d_Ec_pbeueg_rho(r,mu,d_ec_pbeueg_rhoa,d_ec_pbeueg_rhob)
 BEGIN_DOC
  !Derivative of the epsilon_{c,md}^{sr,PBE} function (eq 14a of ref J.Chem.Lett 2019, 2931-2937) with respect to rhoa and rhob 
 END_DOC
 implicit none
 double precision, intent(in)  :: r(3),mu
 double precision, intent(out) :: d_ec_pbeueg_rhoa(N_states),d_ec_pbeueg_rhob(N_states)
 integer :: istate,m 
 double precision :: rhoa(N_states),rhob(N_states),grad_rho_a(3,N_states),grad_rho_b(3,N_states)
 double precision :: aos_array(ao_num),grad_aos_array(3,ao_num)
 double precision :: grad_rho_a_2(N_states),grad_rho_b_2(N_states),grad_rho_a_b(N_states)
 double precision :: wignerseitz_radius,xi,rhot,rs,f_pbeueg,c_beta,pi,rhoc,rhoo,sigmacc,sigmaco,sigmaoo,d_f_pbeueg_rhoa,d_f_pbeueg_rhob,e_PBE
 double precision :: vsigmacc,vsigmaco,vsigmaoo,vrhoo,vrhoc,vc_rho_a,vc_rho_b
 double precision :: beta(N_states),d_beta_rhoa(N_states),d_beta_rhob(N_states)
 pi=dacos(-1.d0) 
 c_beta=3.d0/(2.d0*sqrt(pi)*(1.d0-sqrt(2.d0)))

 call density_and_grad_alpha_beta_and_all_aos_and_grad_aos_at_r(r,rhoa,rhob, grad_rho_a, grad_rho_b, aos_array, grad_aos_array)

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
  rhot=rhoa(istate)+rhob(istate)
  rs = wignerseitz_radius(rhot)
  xi= (rhoa(istate)-rhob(istate))/(rhot)

  call rho_ab_to_rho_oc(rhoa(istate),rhob(istate),rhoo,rhoc)
  call grad_rho_ab_to_grad_rho_oc(grad_rho_a_2(istate),grad_rho_b_2(istate),grad_rho_a_b(istate),sigmaoo,sigmacc,sigmaco)
  call ec_pbe_sr(0.d0,rhoc,rhoo,sigmacc,sigmaco,sigmaoo,e_PBE,vrhoc,vrhoo,vsigmacc,vsigmaco,vsigmaoo)
  call v_rho_oc_to_v_rho_ab(vrhoo,vrhoc,vc_rho_a,vc_rho_b)
  call d_beta_pbeueg_rho(rhoa,rhob,grad_rho_a,grad_rho_b,d_beta_rhoa,d_beta_rhob)

  beta(istate)=c_beta*e_PBE/(rhot**2.d0*f_pbeueg(rs,xi))

  if ((dabs((1.d0+beta(istate)*mu**3.d0)**2.d0) > 1.d-12) .AND. (dabs((rhot**2.d0*f_pbeueg(rs,xi))) > 1.d-12) ) then
   beta(istate)=c_beta*e_PBE/(rhot**2.d0*f_pbeueg(rs,xi))
   d_ec_pbeueg_rhoa(istate) = (vc_rho_a*(1.d0+beta(istate)*mu**3.d0)-e_PBE*(d_beta_rhoa(istate)*mu**3.d0))/((1.d0+beta(istate)*mu**3.d0)**2.d0)
   d_ec_pbeueg_rhob(istate) = (vc_rho_b*(1.d0+beta(istate)*mu**3.d0)-e_PBE*(d_beta_rhob(istate)*mu**3.d0))/((1.d0+beta(istate)*mu**3.d0)**2.d0) 
  else
   d_ec_pbeueg_rhoa(istate) = 0.d0
   d_ec_pbeueg_rhob(istate) = 0.d0
  endif

 enddo
 end


subroutine d_beta_pbeueg_d_ec_pbe(rhoa,rhob,d_beta_ec)
 BEGIN_DOC
  !Derivative of the epsilon_{c,md}^{sr,PBE} function (eq 14a of ref J.Chem.Lett 2019, 2931-2937) with respect to e_c^pbe 
 END_DOC
 implicit none
 double precision, intent(in)  :: rhoa(N_states),rhob(N_states)
 double precision, intent(out) :: d_beta_ec(N_states)
 integer :: istate,m 
 double precision :: wignerseitz_radius,xi,rhot,rs,f_pbeueg,c_beta,pi

 pi=dacos(-1.d0) 
 c_beta=3.d0/(2.d0*sqrt(pi)*(1.d0-sqrt(2.d0)))

 do istate = 1, N_states
  rhot=rhoa(istate)+rhob(istate)
  rs = wignerseitz_radius(rhot)
  xi= (rhoa(istate)-rhob(istate))/(rhot)

  if (dabs((rhot**2.d0*f_pbeueg(rs,xi))**2.d0) > 1.d-12) then
   d_beta_ec(istate) = c_beta/(rhot**2.d0*f_pbeueg(rs,xi))
  else
   d_beta_ec(istate) = 0.d0 
  endif
 enddo

 end


subroutine give_d_Ec_pbeueg_d_grad_n(r,mu,d_ec_pbeueg_d_ec_pbe,d_ec_pbeueg_d_grad_n_alpha,d_ec_pbeueg_d_grad_n_beta)
 BEGIN_DOC
  !Derivative of the epsilon_{c,md}^{sr,PBE} function (eq 14a of ref J.Chem.Lett 2019, 2931-2937) with respect to e_c^pbe and with respect to grad(n_{alpha}) 
 END_DOC
 implicit none
 double precision, intent(in)  :: r(3),mu
 double precision, intent(out) :: d_ec_pbeueg_d_ec_pbe(N_states),d_ec_pbeueg_d_grad_n_alpha(3,N_states),d_ec_pbeueg_d_grad_n_beta(3,N_states)
 integer :: istate,m 
 double precision :: rhoa(N_states),rhob(N_states),grad_rho_a(3,N_states),grad_rho_b(3,N_states)
 double precision :: aos_array(ao_num),grad_aos_array(3,ao_num)
 double precision :: grad_rho_a_2(N_states),grad_rho_b_2(N_states),grad_rho_a_b(N_states)
 double precision :: wignerseitz_radius,xi,rhot,rs,f_pbeueg,c_beta,pi,rhoc,rhoo,sigmacc,sigmaco,sigmaoo,e_PBE
 double precision :: vsigmacc,vsigmaco,vsigmaoo,vrhoo,vrhoc,vc_rho_a,vc_rho_b
 double precision :: vc_grad_rho_a_2(N_states),vc_grad_rho_b_2(N_states),vc_grad_rho_a_b(N_states)
 double precision :: beta(N_states),d_beta_ec(N_states)
 pi=dacos(-1.d0) 
 c_beta=3.d0/(2.d0*sqrt(pi)*(1.d0-sqrt(2.d0)))

 call density_and_grad_alpha_beta_and_all_aos_and_grad_aos_at_r(r,rhoa,rhob, grad_rho_a, grad_rho_b, aos_array, grad_aos_array)

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
  rhot=rhoa(istate)+rhob(istate)
  rs = wignerseitz_radius(rhot)
  xi= (rhoa(istate)-rhob(istate))/(rhot)

  call rho_ab_to_rho_oc(rhoa(istate),rhob(istate),rhoo,rhoc)
  call grad_rho_ab_to_grad_rho_oc(grad_rho_a_2(istate),grad_rho_b_2(istate),grad_rho_a_b(istate),sigmaoo,sigmacc,sigmaco)
  call ec_pbe_sr(0.d0,rhoc,rhoo,sigmacc,sigmaco,sigmaoo,e_PBE,vrhoc,vrhoo,vsigmacc,vsigmaco,vsigmaoo)

  call v_grad_rho_oc_to_v_grad_rho_ab(vsigmaoo,vsigmacc,vsigmaco,vc_grad_rho_a_2(istate),vc_grad_rho_b_2(istate),vc_grad_rho_a_b(istate))

  beta(istate)=c_beta*e_PBE/(rhot**2.d0*f_pbeueg(rs,xi))

  if ((dabs((1.d0+beta(istate)*mu**3.d0)**2.d0) > 1.d-12) .AND. (dabs((rhot**2.d0*f_pbeueg(rs,xi))) > 1.d-12) ) then
   call d_beta_pbeueg_d_ec_pbe(rhoa,rhob,d_beta_ec) 
   d_ec_pbeueg_d_ec_pbe(istate) = ((1.d0+beta(istate)*mu**3.d0)-e_PBE*(d_beta_ec(istate)*mu**3.d0))/((1.d0+beta(istate)*mu**3.d0)**2.d0)
  else
   d_ec_pbeueg_d_ec_pbe(istate) = 0.d0
  endif

  do m = 1, 3
   d_ec_pbeueg_d_grad_n_alpha(m,istate) = d_ec_pbeueg_d_ec_pbe(istate)*( 2.d0*vc_grad_rho_a_2(istate)*grad_rho_a(m,istate) +  vc_grad_rho_a_b(istate)  * grad_rho_b(m,istate) )
   d_ec_pbeueg_d_grad_n_beta(m,istate)  = d_ec_pbeueg_d_ec_pbe(istate)*( 2.d0*vc_grad_rho_b_2(istate)*grad_rho_b(m,istate) +  vc_grad_rho_a_b(istate)  * grad_rho_a(m,istate) )
  enddo
 enddo
 end

 BEGIN_PROVIDER[double precision, aos_sr_vc_alpha_pbe_ueg_w  , (ao_num,n_points_final_grid,N_states)]
&BEGIN_PROVIDER[double precision, aos_sr_vc_beta_pbe_ueg_w   , (ao_num,n_points_final_grid,N_states)]
&BEGIN_PROVIDER[double precision, aos_dsr_vc_alpha_pbe_ueg_w  , (ao_num,n_points_final_grid,N_states)]
&BEGIN_PROVIDER[double precision, aos_dsr_vc_beta_pbe_ueg_w   ,  (ao_num,n_points_final_grid,N_states)]
 implicit none
 BEGIN_DOC
  ! Intermediate quantities for the calculation of the vc potentials for the PBE UEG functional
 END_DOC
 integer :: istate,i,j,m
 double precision :: r(3)
 double precision :: mu,weight
 double precision :: d_ec_pbeueg_d_ec_pbe(N_states),d_ec_pbeueg_d_grad_n_alpha(3,N_states),d_ec_pbeueg_d_grad_n_beta(3,N_states)
 double precision :: d_ec_pbeueg_rhoa(N_states),d_ec_pbeueg_rhob(N_states)

 aos_dsr_vc_alpha_pbe_ueg_w= 0.d0
 aos_dsr_vc_beta_pbe_ueg_w = 0.d0
 do istate = 1, N_states
  do i = 1, n_points_final_grid
   r(1) = final_grid_points(1,i)
   r(2) = final_grid_points(2,i)
   r(3) = final_grid_points(3,i)
   mu = mu_of_r_vector(i)
   weight = final_weight_at_r_vector(i)

   call give_d_Ec_pbeueg_rho(r,mu,d_ec_pbeueg_rhoa,d_ec_pbeueg_rhob)
   call give_d_Ec_pbeueg_d_grad_n(r,mu,d_ec_pbeueg_d_ec_pbe,d_ec_pbeueg_d_grad_n_alpha,d_ec_pbeueg_d_grad_n_beta)

   do j = 1, ao_num
    aos_sr_vc_alpha_pbe_ueg_w(j,i,istate) = d_ec_pbeueg_rhoa(istate) * weight * aos_in_r_array(j,i)
    aos_sr_vc_beta_pbe_ueg_w (j,i,istate) = d_ec_pbeueg_rhob(istate) * weight * aos_in_r_array(j,i)
   enddo
   do j = 1, ao_num
    do m = 1,3
     aos_dsr_vc_alpha_pbe_ueg_w(j,i,istate) += weight * d_ec_pbeueg_d_grad_n_alpha(m,istate) * aos_grad_in_r_array_transp_xyz(m,j,i)
     aos_dsr_vc_beta_pbe_ueg_w (j,i,istate) += weight * d_ec_pbeueg_d_grad_n_beta(m,istate)  * aos_grad_in_r_array_transp_xyz(m,j,i)
    enddo

   enddo
  enddo
 enddo

 END_PROVIDER


 BEGIN_PROVIDER [double precision, pot_sr_scal_c_alpha_ao_pbe_ueg, (ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, pot_sr_scal_c_beta_ao_pbe_ueg, (ao_num,ao_num,N_states)]
 implicit none
 integer                        :: istate
   BEGIN_DOC
   ! Intermediate quantities for the calculation of the vc potentials related to the scalar part for the PBE UEG functional
   END_DOC
   pot_sr_scal_c_alpha_ao_pbe_ueg = 0.d0
   pot_sr_scal_c_beta_ao_pbe_ueg = 0.d0
   double precision               :: wall_1,wall_2
   call wall_time(wall_1)
   do istate = 1, N_states
     ! correlation alpha
     call dgemm('N','T',ao_num,ao_num,n_points_final_grid,1.d0,                                                                                       &
                 aos_sr_vc_alpha_pbe_ueg_w(1,1,istate),size(aos_sr_vc_alpha_pbe_ueg_w,1),                                                                   &
                 aos_in_r_array,size(aos_in_r_array,1),1.d0,                                                                                          &
                 pot_sr_scal_c_alpha_ao_pbe_ueg(1,1,istate),size(pot_sr_scal_c_alpha_ao_pbe_ueg,1))
     ! correlation beta
     call dgemm('N','T',ao_num,ao_num,n_points_final_grid,1.d0,                                                                                         &
                 aos_sr_vc_beta_pbe_ueg_w(1,1,istate),size(aos_sr_vc_beta_pbe_ueg_w,1),                                                                       &
                 aos_in_r_array,size(aos_in_r_array,1),1.d0,                                                                                            &
                 pot_sr_scal_c_beta_ao_pbe_ueg(1,1,istate),size(pot_sr_scal_c_beta_ao_pbe_ueg,1))
 
   enddo
 call wall_time(wall_2)

END_PROVIDER



 
 BEGIN_PROVIDER [double precision, pot_sr_grad_c_alpha_ao_pbe_ueg,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, pot_sr_grad_c_beta_ao_pbe_ueg,(ao_num,ao_num,N_states)]
   implicit none
   BEGIN_DOC
   ! Intermediate quantity for the calculation of the vc potentials for the PBA UEG functional related to the gradienst of the density and orbitals 
   END_DOC
   integer                        :: istate
   double precision               :: wall_1,wall_2
   call wall_time(wall_1)
   pot_sr_grad_c_alpha_ao_pbe_ueg = 0.d0
   pot_sr_grad_c_beta_ao_pbe_ueg = 0.d0
   do istate = 1, N_states
       ! correlation alpha
       call dgemm('N','N',ao_num,ao_num,n_points_final_grid,1.d0,                                                                                             &
                  aos_dsr_vc_alpha_pbe_ueg_w(1,1,istate),size(aos_dsr_vc_alpha_pbe_ueg_w,1),                                                                      &
                  aos_in_r_array_transp,size(aos_in_r_array_transp,1),1.d0,                                                                                &
                  pot_sr_grad_c_alpha_ao_pbe_ueg(1,1,istate),size(pot_sr_grad_c_alpha_ao_pbe_ueg,1))
       ! correlation beta
       call dgemm('N','N',ao_num,ao_num,n_points_final_grid,1.d0,                                                                                             &
                  aos_dsr_vc_beta_pbe_ueg_w(1,1,istate),size(aos_dsr_vc_beta_pbe_ueg_w,1),                                                                      &
                  aos_in_r_array_transp,size(aos_in_r_array_transp,1),1.d0,                                                                                &
                  pot_sr_grad_c_beta_ao_pbe_ueg(1,1,istate),size(pot_sr_grad_c_beta_ao_pbe_ueg,1))
   enddo
   
 call wall_time(wall_2)

END_PROVIDER

 BEGIN_PROVIDER [double precision, potential_c_alpha_ao_sr_pbe_ueg,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, potential_c_beta_ao_sr_pbe_ueg,(ao_num,ao_num,N_states)]
   implicit none
 BEGIN_DOC
 ! Correlation potential for alpha / beta electrons  with the PBE UEG functional 
 END_DOC 
   integer :: i,j,istate
   do istate = 1, n_states 
    do i = 1, ao_num
     do j = 1, ao_num
      potential_c_alpha_ao_sr_pbe_ueg(j,i,istate) = pot_sr_scal_c_alpha_ao_pbe_ueg(j,i,istate) + pot_sr_grad_c_alpha_ao_pbe_ueg(j,i,istate) + pot_sr_grad_c_alpha_ao_pbe_ueg(i,j,istate)
      potential_c_beta_ao_sr_pbe_ueg(j,i,istate) = pot_sr_scal_c_beta_ao_pbe_ueg(j,i,istate) + pot_sr_grad_c_beta_ao_pbe_ueg(j,i,istate) + pot_sr_grad_c_beta_ao_pbe_ueg(i,j,istate)
     enddo
    enddo
   enddo

END_PROVIDER


 BEGIN_PROVIDER [double precision, potential_c_alpha_mo_sr_pbe_ueg,(mo_num,mo_num,N_states)]
&BEGIN_PROVIDER [double precision, potential_c_beta_mo_sr_pbe_ueg, (mo_num,mo_num,N_states)]
 implicit none
 BEGIN_DOC
! Providers for the alpha/beta correlation potentials of the PBE UEG on the MO basis
 END_DOC
 integer :: istate
 do istate = 1, N_states
    call ao_to_mo(                                                   &
        potential_c_alpha_ao_sr_pbe_ueg(1,1,istate),                                 &
        size(potential_c_alpha_ao_sr_pbe_ueg,1),                                &
        potential_c_alpha_mo_sr_pbe_ueg(1,1,istate),                                 &
        size(potential_c_alpha_mo_sr_pbe_ueg,1)                                 &
        )

    call ao_to_mo(                                                   &
        potential_c_beta_ao_sr_pbe_ueg(1,1,istate),                                  &
        size(potential_c_beta_ao_sr_pbe_ueg,1),                                 &
        potential_c_beta_mo_sr_pbe_ueg(1,1,istate),                                  &
        size(potential_c_beta_mo_sr_pbe_ueg,1)                                  &
        )
 enddo

END_PROVIDER


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Potential used only for tests !!!!!!!!!!!!!!!!

 
subroutine give_d_Ec_pbeueg_rho_grand_mu(r,mu,d_ec_pbeueg_grandmu_rhoa,d_ec_pbeueg_grandmu_rhob)
 BEGIN_DOC
  !Blavblablablablablalbalbla
 END_DOC
 implicit none
 double precision, intent(in)  :: r(3),mu
 double precision, intent(out) :: d_ec_pbeueg_grandmu_rhoa(N_states),d_ec_pbeueg_grandmu_rhob(N_states)
 integer :: istate,m 
 double precision :: rhoa(N_states),rhob(N_states),grad_rho_a(3,N_states),grad_rho_b(3,N_states)
 double precision :: aos_array(ao_num),grad_aos_array(3,ao_num)
 double precision :: wignerseitz_radius,xi,rhot,rs,f_pbeueg,c_beta,pi,d_f_pbeueg_rhoa,d_f_pbeueg_rhob
 pi=dacos(-1.d0) 
 c_beta=3.d0/(2.d0*sqrt(pi)*(1.d0-sqrt(2.d0)))

 call density_and_grad_alpha_beta_and_all_aos_and_grad_aos_at_r(r,rhoa,rhob, grad_rho_a, grad_rho_b, aos_array, grad_aos_array)

 do istate = 1, N_states
  ! convertion from (alpha,beta) formalism to (closed, open) formalism
  rhot=rhoa(istate)+rhob(istate)
  rs = wignerseitz_radius(rhot)
  xi= (rhoa(istate)-rhob(istate))/(rhot)

  d_ec_pbeueg_grandmu_rhoa(istate) =  (2.d0* rhot * f_pbeueg(rs,xi) + rhot**2.d0 * d_f_pbeueg_rhoa(rhoa(istate),rhob(istate)))/(c_beta*mu**3.d0)
  d_ec_pbeueg_grandmu_rhob(istate) =  (2.d0* rhot * f_pbeueg(rs,xi) + rhot**2.d0 * d_f_pbeueg_rhob(rhoa(istate),rhob(istate)))/(c_beta*mu**3.d0)

 enddo
 end



 BEGIN_PROVIDER[double precision, aos_sr_vc_alpha_pbe_ueg_grandmu  , (ao_num,n_points_final_grid,N_states)]
&BEGIN_PROVIDER[double precision, aos_sr_vc_beta_pbe_ueg_grandmu   , (ao_num,n_points_final_grid,N_states)]
 implicit none
 BEGIN_DOC

  !Derivative of the epsilon_{c,md}^{sr,PBE} function (eq 14a of ref J.Chem.Lett 2019, 2931-2937) with respect to e_c^pbe and with respect to grad(n_{alpha}) 
  aos_sr_vxc_alpha_pbe_ueg_w(j,i) = ao_i(r_j) * (v^x_alpha(r_j) + v^c_alpha(r_j)) * W(r_j)
 END_DOC
 integer :: istate,i,j,m
 double precision :: r(3)
 double precision :: mu,weight
 double precision :: d_ec_pbeueg_grandmu_rhoa(N_states),d_ec_pbeueg_grandmu_rhob(N_states)

 aos_sr_vc_alpha_pbe_ueg_grandmu = 0.d0
 aos_sr_vc_beta_pbe_ueg_grandmu  = 0.d0

 do istate = 1, N_states
  do i = 1, n_points_final_grid
   r(1) = final_grid_points(1,i)
   r(2) = final_grid_points(2,i)
   r(3) = final_grid_points(3,i)
   mu = mu_of_r_vector(i)
  !mu = mu_erf_dft 
   weight = final_weight_at_r_vector(i)

   call give_d_Ec_pbeueg_rho_grand_mu(r,mu,d_ec_pbeueg_grandmu_rhoa,d_ec_pbeueg_grandmu_rhob)

   do j = 1, ao_num
    aos_sr_vc_alpha_pbe_ueg_grandmu(j,i,istate) = d_ec_pbeueg_grandmu_rhoa(istate) * weight * aos_in_r_array(j,i)
    aos_sr_vc_beta_pbe_ueg_grandmu (j,i,istate) = d_ec_pbeueg_grandmu_rhob(istate) * weight * aos_in_r_array(j,i)
   enddo

  enddo
 enddo

 END_PROVIDER


 BEGIN_PROVIDER [double precision, pot_sr_scal_c_alpha_ao_pbe_ueg_grandmu, (ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, pot_sr_scal_c_beta_ao_pbe_ueg_grandmu, (ao_num,ao_num,N_states)]
 implicit none
 integer                        :: istate
   BEGIN_DOC
   ! intermediate quantity for the calculation of the vxc potentials for the GGA functionals  related to the scalar part of the potential 
   END_DOC
   pot_sr_scal_c_alpha_ao_pbe_ueg_grandmu = 0.d0
   pot_sr_scal_c_beta_ao_pbe_ueg_grandmu = 0.d0
   double precision               :: wall_1,wall_2
   call wall_time(wall_1)
   do istate = 1, N_states
     ! correlation alpha
     call dgemm('N','T',ao_num,ao_num,n_points_final_grid,1.d0,                                                                                       &
                 aos_sr_vc_alpha_pbe_ueg_grandmu(1,1,istate),size(aos_sr_vc_alpha_pbe_ueg_grandmu,1),                                                                   &
                 aos_in_r_array,size(aos_in_r_array,1),1.d0,                                                                                          &
                 pot_sr_scal_c_alpha_ao_pbe_ueg_grandmu(1,1,istate),size(pot_sr_scal_c_alpha_ao_pbe_ueg_grandmu,1))
     ! correlation beta
     call dgemm('N','T',ao_num,ao_num,n_points_final_grid,1.d0,                                                                                         &
                 aos_sr_vc_beta_pbe_ueg_grandmu(1,1,istate),size(aos_sr_vc_beta_pbe_ueg_grandmu,1),                                                                       &
                 aos_in_r_array,size(aos_in_r_array,1),1.d0,                                                                                            &
                 pot_sr_scal_c_beta_ao_pbe_ueg_grandmu(1,1,istate),size(pot_sr_scal_c_beta_ao_pbe_ueg_grandmu,1))
 
   enddo
 call wall_time(wall_2)

END_PROVIDER

 

 BEGIN_PROVIDER [double precision, potential_c_alpha_ao_sr_pbe_ueg_grandmu,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, potential_c_beta_ao_sr_pbe_ueg_grandmu,(ao_num,ao_num,N_states)]
   implicit none
 BEGIN_DOC
 ! Correlation potential for alpha / beta electrons  with the PBE UEG functional 
 END_DOC 
   integer :: i,j,istate
   do istate = 1, n_states 
    do i = 1, ao_num
     do j = 1, ao_num
      potential_c_alpha_ao_sr_pbe_ueg_grandmu(j,i,istate) = pot_sr_scal_c_alpha_ao_pbe_ueg_grandmu(j,i,istate) 
      potential_c_beta_ao_sr_pbe_ueg_grandmu(j,i,istate) = pot_sr_scal_c_beta_ao_pbe_ueg_grandmu(j,i,istate) 
     enddo
    enddo
   enddo

END_PROVIDER





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 2nd order pot Scalar Part!!!! 
