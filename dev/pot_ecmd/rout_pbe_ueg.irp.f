

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

  if(dabs(rhot**2.d0*f_pbeueg(rs,xi)).gt.1.d-15)then
   beta(istate)=c_beta*e_PBE/(rhot**2.d0*f_pbeueg(rs,xi))
  else
   beta(istate)=1.d+10
  endif

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

  if(dabs(rhot**2.d0*f_pbeueg(rs,xi)).gt.1.d-10)then
   beta(istate)=c_beta*e_PBE/(rhot**2.d0*f_pbeueg(rs,xi))
  else
   beta(istate)=1.d+10
  endif

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

