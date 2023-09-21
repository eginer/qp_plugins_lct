  subroutine exmdsrLDAn(mu,rho_a,rho_b,ex_srmuLDAn,dexdrho_a,dexdrho_b)
  
  implicit none
  BEGIN_DOC
  ! Calculation of exchange energy and chemical potential in LDAn approximation using multideterminantal wave function (short-range part)
  END_DOC
  double precision, intent(in)  :: mu
  double precision, intent(in)  :: rho_a, rho_b
  double precision, intent(out) :: ex_srmuLDAn, dexdrho_a, dexdrho_b
 ! double precision              :: dexdrho, dexdgrad_rho_2
  double precision              :: exLDAn, dexLDAndrho_a, dexLDAndrho_b, dexLDAndrho
  double precision              :: gamma, dgammadrho_a, dgammadrho_b
  double precision              :: delta, ddeltadrho_a, ddeltadrho_b
  double precision              :: denom, ddenomdrho_a, ddenomdrho_b
  double precision              :: pi, a, b, thr
  double precision              :: rho, m  
  double precision              :: n2_UEG, dn2_UEGdrho, dn2_UEGdrho_a, dn2_UEGdrho_b
  double precision              :: n2xc_UEG, dn2xc_UEGdrho, dn2xc_UEGdrho_a, dn2xc_UEGdrho_b
  double precision              :: g0, dg0drho

  if(dabs(rho_a-rho_b)/dabs(rho_a+rho_b) > 1.d-1)then
  print*,dabs(rho_a-rho_b),dabs(rho_a+rho_b),dabs(rho_a-rho_b)/dabs(rho_a+rho_b)
  stop "routine implemented only for closed-shell systems"
  endif 

  pi = dacos(-1.d0)
  rho = rho_a + rho_b
  m = rho_a - rho_b
  thr = 1.d-12

! exchange LDAn standard and on-top pair distribution
 ! call ex_pbe_sr(1.d-12,rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b,exLDAn,dexLDAndrho_a,dexLDAndrho_b,dexLDAndgrad_rho_a_2,dexLDAndgrad_rho_b_2,dexLDAndgrad_rho_a_b)
  call ex_lda(rho_a,rho_b,exLDAn,dexLDAndrho_a,dexLDAndrho_b)
  call g0_dg0(rho, rho_a, rho_b, g0, dg0drho)
  
  if(dabs(exLDAn).lt.thr)then
   exLDAn = 1.d-12
  endif

  if(dabs(rho).lt.thr)then
   rho = 1.d-12
  endif

  if(dabs(g0).lt.thr)then
    g0 = 1.d-12
  endif

! calculation of energy
  a = pi / 2.d0
  b = 2*dsqrt(pi)*(2*dsqrt(2.d0) - 1.d0)/3.d0   
  
  n2_UEG = (rho**2)*g0
  if(dabs(n2_UEG).lt.thr)then
   n2_UEG = 1.d-12
  endif
  
  n2xc_UEG = n2_UEG - rho**2
  if(dabs(n2xc_UEG).lt.thr)then
   n2xc_UEG = 1.d-12
  endif

  gamma = exLDAn / (a*n2xc_UEG)
  if(dabs(gamma).lt.thr)then
   gamma = 1.d-12
  endif

  delta = -(b*n2_UEG*gamma**2) / exLDAn
  if(dabs(delta).lt.thr)then
   delta = 1.d-12
  endif
  
  denom = 1.d0 + delta*mu + gamma*(mu**2)

  ex_srmuLDAn=exLDAn/denom

! calculation of derivatives
  !dex/dn
  dn2_UEGdrho     = 2.d0*rho*g0 + (rho**2)*dg0drho
  dn2_UEGdrho_a   = 2.d0*rho*g0 + (rho**2)*dg0drho
  dn2_UEGdrho_b   = 2.d0*rho*g0 + (rho**2)*dg0drho
  dn2xc_UEGdrho   = dn2_UEGdrho - 2.d0*rho
  dn2xc_UEGdrho_a = dn2xc_UEGdrho
  dn2xc_UEGdrho_b = dn2xc_UEGdrho
  
  dgammadrho_a = (1.d0/(a*n2xc_UEG))*dexLDAndrho_a  - (exLDAn/(a*n2xc_UEG**2))*dn2xc_UEGdrho_a
  ddeltadrho_a = -((b*gamma**2)/exLDAn)*dn2_UEGdrho_a -(b*n2_UEG/exLDAn)*2.d0*gamma*dgammadrho_a + (b*n2_UEG*gamma**2)/(exLDAn**2)*dexLDAndrho_a
  dgammadrho_b = (1.d0/(a*n2xc_UEG))*dexLDAndrho_b  - (exLDAn/(a*n2xc_UEG**2))*dn2xc_UEGdrho_b
  ddeltadrho_b = -((b*gamma**2)/exLDAn)*dn2_UEGdrho_b -(b*n2_UEG/exLDAn)*2.d0*gamma*dgammadrho_b + (b*n2_UEG*gamma**2)/(exLDAn**2)*dexLDAndrho_b

  ddenomdrho_a = ddeltadrho_a * mu + dgammadrho_a * (mu**2)
  ddenomdrho_b = ddeltadrho_b * mu + dgammadrho_b * (mu**2)
  
  dexdrho_a = dexLDAndrho_a/denom - exLDAn*ddenomdrho_a/denom**2
  dexdrho_b = dexLDAndrho_b/denom - exLDAn*ddenomdrho_b/denom**2
  
!  dexdrho = 0.5d0*(dexdrho_a + dexdrho_b)
 
! !dex/d((gradn)^2)
! dgammadgrad_rho_a_2 =dexLDAndgrad_rho_a_2/(a*n2xc_UEG)
! dgammadgrad_rho_b_2 =dexLDAndgrad_rho_b_2/(a*n2xc_UEG)
! dgammadgrad_rho_a_b =dexLDAndgrad_rho_a_b/(a*n2xc_UEG)

! ddeltadgrad_rho_a_2 = ((b*n2_UEG*gamma**2)/(exLDAn**2))*dexLDAndgrad_rho_a_2 - b*(n2_UEG/exLDAn)*2.d0*gamma*dgammadgrad_rho_a_2
! ddeltadgrad_rho_b_2 = ((b*n2_UEG*gamma**2)/(exLDAn**2))*dexLDAndgrad_rho_b_2 - b*(n2_UEG/exLDAn)*2.d0*gamma*dgammadgrad_rho_b_2
! ddeltadgrad_rho_a_b = ((b*n2_UEG*gamma**2)/(exLDAn**2))*dexLDAndgrad_rho_a_b - b*(n2_UEG/exLDAn)*2.d0*gamma*dgammadgrad_rho_a_b

! ddenomdgrad_rho_a_2 = ddeltadgrad_rho_a_2*mu + dgammadgrad_rho_a_2*mu**2  
! ddenomdgrad_rho_b_2 = ddeltadgrad_rho_b_2*mu + dgammadgrad_rho_b_2*mu**2
! ddenomdgrad_rho_a_b =ddeltadgrad_rho_a_b*mu + dgammadgrad_rho_a_b*mu**2
! 
! dexdgrad_rho_a_2 = dexLDAndgrad_rho_a_2/denom - exLDAn*ddenomdgrad_rho_a_2/(denom**2)
! dexdgrad_rho_b_2 = dexLDAndgrad_rho_b_2/denom  - exLDAn*ddenomdgrad_rho_b_2/(denom**2)
! dexdgrad_rho_a_b = dexLDAndgrad_rho_a_b/denom - exLDAn*ddenomdgrad_rho_a_b/(denom**2)

!  dexdgrad_rho_2 = 0.25d0*(dexdgrad_rho_a_2 + dexdgrad_rho_b_2 + dexdgrad_rho_a_b)
  
!  print*, '..................................'
!  print*, 'rhoa                =', rho_a
!  print*, 'rhob                =', rho_b
!  print*, 'gradrho_a_2         =', grad_rho_a_2
!  print*, 'gradrho_b_2         =', grad_rho_b_2
!  print*, 'gradrho_a_b         =', grad_rho_a_b
!  print*, 'exLDAn =', exLDAn, 'ex_srmuLDAn =', ex_srmuLDAn
!  print*, 'dexLDAndrho_a        =', dexLDAndrho_a 
!  print*, 'dexLDAndrho_b        =', dexLDAndrho_b
!  print*, 'dexLDAndgrad_rho_a_2 =', dexLDAndgrad_rho_a_2
!  print*, 'dexLDAndgrad_rho_b_2 =', dexLDAndgrad_rho_b_2
!  print*, 'dexLDAndgrad_rho_a_b= ', dexLDAndgrad_rho_a_b
!  print*, 'dexLDAndgrad_rho_2   =', dexLDAndgrad_rho_2
!  print*, 'dexdrho_a           =', dexdrho_a          
!  print*, 'dexdrho_b           =', dexdrho_b          
!  print*, 'dexdgrad_rho_a_2    =', dexdgrad_rho_a_2   
!  print*, 'dexdgrad_rho_b_2    =', dexdgrad_rho_b_2   
!  print*, 'dexdgrad_rho_a_b    =', dexdgrad_rho_a_b   
!  print*, 'dgammadgrad_rho_a_2 =', dgammadgrad_rho_a_2  
!  print*, 'ddeltadgrad_rho_a_2 =', ddeltadgrad_rho_a_2  
!  print*, 'dgammadgrad_rho_b_2 =', dgammadgrad_rho_b_2  
!  print*, 'ddeltadgrad_rho_b_2 =', ddeltadgrad_rho_b_2  
!  print*, 'dgammadgrad_rho_a_b =', dgammadgrad_rho_a_b  
!  print*, 'ddeltadgrad_rho_a_b =', ddeltadgrad_rho_a_b
!  print*, 'denom               =', denom
!  print*, 'ddenomdgrad_rho_a_2 =', ddenomdgrad_rho_a_2 
!  print*, 'ddenomdgrad_rho_b_2 =', ddenomdgrad_rho_b_2 
!  print*, 'ddenomdgrad_rho_a_b =', ddenomdgrad_rho_a_b 
!  print*, '------> dexdgradrho =', dexdgrad_rho_2   
! print*, '..................................'
  end subroutine exmdsrLDAn
!---------------------------------------------------------------------------------------------------------------------------------------------

  subroutine ecmdsrLDAn(mu,rho_a,rho_b,ec_srmuLDAn,decdrho_a,decdrho_b,d2ecdrho_a2,d2ecdrho_b2)

  implicit none
  BEGIN_DOC
  ! Calculation of correlation energy and chemical potential in LDAn approximation using multideterminantal wave function (short-range part)
  END_DOC
 
  double precision, intent(in)  :: mu
  double precision, intent(in)  :: rho_a, rho_b
  double precision, intent(out) :: ec_srmuLDAn, decdrho_a, decdrho_b
  double precision, intent(out) :: d2ecdrho_a2, d2ecdrho_b2
  double precision              :: decdrho
  double precision              :: ecLDAn, decLDAndrho_a, decLDAndrho_b
  double precision              :: rho_c, rho_o, decLDAndrho_c, decLDAndrho_o
  double precision              :: beta, dbetadrho_a, dbetadrho_b
  double precision              :: denom, ddenomdrho_a, ddenomdrho_b
  double precision              :: pi, c, thr
  double precision              :: rho, m  
  double precision              :: n2_UEG, dn2_UEGdrho_a, dn2_UEGdrho_b
  double precision              :: g0, dg0drho, d2g0drho2

  if(dabs(rho_a-rho_b)/dabs(rho_a+rho_b) > 1.d-1)then
   print*,'rho_a,rho_b        = ',rho_a,rho_b
   print*,'dabs(rho_a-rho_b)  = ',dabs(rho_a-rho_b)
   stop "routine implemented only for closed-shell systems"
  endif 

  pi = dacos(-1.d0)
  rho = rho_a + rho_b

  m = rho_a - rho_b

  thr = 1.d-12
  
! correlation LDAn standard and on-top pair distribution 
  call rho_ab_to_rho_oc(rho_a,rho_b,rho_o,rho_c)

!  call ec_pbe_sr(1.d-12,rho_c,rho_o,grad_rho_c_2,grad_rho_o_2,grad_rho_o_c,ecLDAn,decLDAndrho_c,decLDAndrho_o,decLDAndgrad_rho_c_2,decLDAndgrad_rho_o_2, decLDAndgrad_rho_c_o)

  call ec_lda(rho_a,rho_b,ecLDAn,decLDAndrho_a,decLDAndrho_b)

!  call v_rho_oc_to_v_rho_ab(decLDAndrho_o, decLDAndrho_c, decLDAndrho_a, decLDAndrho_b)

  call g0_dg0(rho, rho_a, rho_b, g0, dg0drho)
!  call g0_dg0_d2g0(rho, rho_a, rho_b, g0, dg0drho, d2g0drho2)

! calculation of energy
  c = 2*dsqrt(pi)*(1.d0 - dsqrt(2.d0))/3.d0
  
  n2_UEG = (rho**2)*g0
 
  if(dabs(n2_UEG).lt.thr)then
   n2_UEG = 1.d-12
  endif  

  beta = ecLDAn/(c*n2_UEG)

  if(dabs(beta).lt.thr)then
   beta = 1.d-12
  endif

  denom = 1.d0 + beta*mu**3

  ec_srmuLDAn=ecLDAn/denom

! calculation of derivatives 
 ! dec/dn
  
 dn2_UEGdrho_a = 2.d0*rho*g0 + (rho**2)*dg0drho
 dn2_UEGdrho_b = 2.d0*rho*g0 + (rho**2)*dg0drho
 
 dbetadrho_a  = decLDAndrho_a/(c*n2_UEG) - (ecLDAn/(c*n2_UEG**2))*dn2_UEGdrho_a
 dbetadrho_b  = decLDAndrho_b/(c*n2_UEG) - (ecLDAn/(c*n2_UEG**2))*dn2_UEGdrho_b

 ddenomdrho_a = dbetadrho_a*mu**3
 ddenomdrho_b = dbetadrho_b*mu**3

 decdrho_a = decLDAndrho_a/denom - ecLDAn*ddenomdrho_a/(denom**2)
 decdrho_b = decLDAndrho_b/denom - ecLDAn*ddenomdrho_b/(denom**2)
 
 decdrho = decdrho_a + decdrho_b
 
  end subroutine ecmdsrLDAn

!---------------------------------------------------------------------------------------------------------------------------------------------

subroutine exc_dexc_md_sr_LDAn(mu,rho_a,rho_b, &
       ex_srmuLDAn,dexdrho_a,dexdrho_b, &
       ec_srmuLDAn,decdrho_a,decdrho_b,d2ecdrho_a2,d2ecdrho_b2)
 
 implicit none
 BEGIN_DOC
 ! Give exchange and correlation energies and chemical potentials
 ! Use qp_plugins_dtraore/sr_md_energies/utils_modified_jt.irp.f's plugins : exmdsrLDAn and ecmdsrLDAn
 END_DOC
 double precision, intent(in)  :: mu
 double precision, intent(in)  :: rho_a,rho_b
 double precision, intent(out) :: ex_srmuLDAn,dexdrho_a,dexdrho_b
 double precision, intent(out) :: ec_srmuLDAn,decdrho_a,decdrho_b,d2ecdrho_a2,d2ecdrho_b2

 call exmdsrLDAn(mu,rho_a,rho_b,ex_srmuLDAn,dexdrho_a,dexdrho_b)

 call ecmdsrLDAn(mu,rho_a,rho_b,ec_srmuLDAn,decdrho_a,decdrho_b,d2ecdrho_a2,d2ecdrho_b2)

 end subroutine excmdsrLDAn

 ---------------------------------------------------------------------------------------------------------------------------------------------
 subroutine fc_LDAUEG(fcmdsrLDAUEG)
 implicit none
 BEGIN_DOC
 END_DOC
 double precision, intent(out) :: fcmdsrLDAUEG(n_points_final_grid)
 
 double precision :: ecLDAn,decLDAndrho_a,decLDAndrho_b,decLDAndrho
 double precision :: ecLDAn_delta_a,decLDAndrho_a_delta_a,decLDAndrho_b_delta_a
 double precision :: ec_srmuLDAn,decdrho_a,decdrho_b,d2ecdrho_a2,d2ecdrho_b2
 double precision, allocatable :: kernelc(:)
 double precision :: thr, pi, c, rho, delta_rho_a
 double precision :: g0, dg0drho, d2g0drho2
 double precision :: n2_UEG, dn2_UEGdrho, d2n2_UEGdrho2 
 double precision :: beta, dbetadrho, d2betadrho2
 double precision :: D, E, F, G, H
 double precision :: A1, A2, dA1, dA2
 double precision :: B1, B2, dB1, dB2

 integer :: i
 double precision :: mu, rho_a, rho_b
 allocate(kernelc(n_points_final_grid))

  pi = dacos(-1.d0)
  c = 2*dsqrt(pi)*(1.d0 - dsqrt(2.d0))/3.d0
  thr = 1.d-12

 fcmdsrLDAUEG = 0.d0

 call kernel_ldac(0,kernelc)
 do i=1,n_points_final_grid
  rho_a = one_e_dm_and_grad_alpha_in_r(4,i,1)
  rho_b = one_e_dm_and_grad_beta_in_r(4,i,1)
  mu = mu_of_r_prov(i,1)
  
  if (dabs(rho_a + rho_b) < 1.d-4)then

    fcmdsrLDAUEG(i) = 0.d0

  else

  if(dabs(rho_a-rho_b)/dabs(rho_a+rho_b) > 1.d-1)then
   print*,'rho_a,rho_b        = ',rho_a,rho_b
   print*,'dabs(rho_a-rho_b)  = ',dabs(rho_a-rho_b)
   stop "routine implemented only for closed-shell systems"
  endif 

  rho = rho_a + rho_b
  delta_rho_a = 1.d-12!*(rho_a + 1.d-6)
  !delta_rho_a = 1.d-3*(rho_a + 1.d-3)

  call ec_lda(rho_a,rho_b,ecLDAn,decLDAndrho_a,decLDAndrho_b)
  call ec_lda(rho_a + delta_rho_a,rho_b,ecLDAn_delta_a,decLDAndrho_a_delta_a,decLDAndrho_b_delta_a)
  call ecmdsrLDAn(mu,rho_a,rho_b,ec_srmuLDAn,decdrho_a,decdrho_b,d2ecdrho_a2,d2ecdrho_b2) 
 ! dv/drho = (dv/drho_a)(drho_a/drho) + (dv/drho_b)(drho_b/drho)
 !         = 0.5*(dv/drho_a + dv/drho_b)
 !         = dv/drho_a (because dv/drho_a = dv/drho_b)
 kernelc(i) = (decLDAndrho_a_delta_a - decLDAndrho_a)/(delta_rho_a)
 ! write(33,'(100(F16.10,X))') dsqrt(final_grid_points(1,i)**2+final_grid_points(2,i)**2+final_grid_points(3,i)**2),rho_a+rho_b, kernelc(i), (decLDAndrho_a_delta_a - decLDAndrho_a)/(delta_rho_a) 

  call g0_dg0_d2g0(rho, rho_a, rho_b, g0, dg0drho, d2g0drho2)

  n2_UEG = (rho**2)*g0

  if(dabs(n2_UEG).lt.thr)then
   n2_UEG = 1.d-12
  endif  

  dn2_UEGdrho = 2.d0*rho*g0 + (rho**2)*dg0drho
  d2n2_UEGdrho2 = 2.d0*g0 + 2.d0*rho*dg0drho + 2.d0*rho*dg0drho + (rho**2)*d2g0drho2
  beta = ecLDAn/(c*n2_UEG)
  dbetadrho  = (decLDAndrho_a + decLDAndrho_b)/(c*n2_UEG) - (ecLDAn/(c*n2_UEG**2))*dn2_UEGdrho

  D = kernelc(i)/(c*n2_UEG)
  E = ((decLDAndrho_a + decLDAndrho_b)*dn2_UEGdrho)/(c*n2_UEG**2)
  F = ((decLDAndrho_a + decLDAndrho_b)*dn2_UEGdrho)/(c*n2_UEG**2)
  G = ecLDAn*d2n2_UEGdrho2/(c*n2_UEG**2)
  H = (2.d0*ecLDAn*(dn2_UEGdrho)**2)/(c*n2_UEG**3)

  d2betadrho2 = D - E -F -G + H

  decLDAndrho = decLDAndrho_a + decLDAndrho_b
  A1 = decLDAndrho*(1.d0 + beta*mu**3)
  A2 = (1.d0 + beta*mu**3)**2

  B1 = ecLDAn*(dbetadrho*mu**3)
  B2 = (1.d0 + beta*mu**3)**2

  dA1 = kernelc(i)*(1.d0 + beta*mu**3) + decLDAndrho*dbetadrho*mu**3
  dA2 = 2.d0*(1.d0 + beta*mu**3)*dbetadrho*mu**3

  dB1 = decLDAndrho*dbetadrho*mu**3 + ecLDAn*d2betadrho2*mu**3
  dB2 = 2.d0*(1.d0 + beta*mu**3)*dbetadrho*mu**3
 
  fcmdsrLDAUEG(i) = (dA1*A2 - dA2*A1)/(A2**2) - (dB1*B2 - dB2*B1)/(B2**2) 
 endif
 enddo

 deallocate(kernelc)
 end subroutine fc_LDAUEG
