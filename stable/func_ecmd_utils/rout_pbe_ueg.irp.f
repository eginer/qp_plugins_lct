  subroutine exmdsrPBE(mu,rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b,ex_srmuPBE,dexdrho_a,dexdrho_b, dexdgrad_rho_a_2,dexdgrad_rho_b_2,dexdgrad_rho_a_b)
  
  implicit none
  BEGIN_DOC
  ! Calculation of exchange energy and chemical potential in PBE approximation using multideterminantal wave function (short-range part)
  END_DOC
  double precision, intent(in)  :: mu
  double precision, intent(in)  :: rho_a, rho_b, grad_rho_a_2, grad_rho_b_2, grad_rho_a_b
  double precision, intent(out) :: ex_srmuPBE, dexdrho_a, dexdrho_b, dexdgrad_rho_a_2, dexdgrad_rho_b_2, dexdgrad_rho_a_b
 ! double precision              :: dexdrho, dexdgrad_rho_2
  double precision              :: exPBE, dexPBEdrho_a, dexPBEdrho_b, dexPBEdrho, dexPBEdgrad_rho_a_2, dexPBEdgrad_rho_b_2, dexPBEdgrad_rho_a_b, dexPBEdgrad_rho_2
  double precision              :: gamma, dgammadrho_a, dgammadrho_b, dgammadgrad_rho_a_2, dgammadgrad_rho_b_2, dgammadgrad_rho_a_b
  double precision              :: delta, ddeltadrho_a, ddeltadrho_b, ddeltadgrad_rho_a_2, ddeltadgrad_rho_b_2, ddeltadgrad_rho_a_b
  double precision              :: denom, ddenomdrho_a, ddenomdrho_b, ddenomdgrad_rho_a_2, ddenomdgrad_rho_b_2, ddenomdgrad_rho_a_b  
  double precision              :: pi, a, b, thr
  double precision              :: rho, m  
  double precision              :: n2_UEG, dn2_UEGdrho, dn2_UEGdrho_a, dn2_UEGdrho_b
  double precision              :: n2xc_UEG, dn2xc_UEGdrho, dn2xc_UEGdrho_a, dn2xc_UEGdrho_b
  double precision              :: g0, dg0drho

  if(dabs(rho_a-rho_b) > 1.d-10)then
  stop "routine implemented only for closed-shell systems"
  endif 

  pi = dacos(-1.d0)
  rho = rho_a + rho_b
  m = rho_a - rho_b
  thr = 1.d-12

! exchange PBE standard and on-top pair distribution
  call ex_pbe_sr(1.d-12,rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b,exPBE,dexPBEdrho_a,dexPBEdrho_b,dexPBEdgrad_rho_a_2,dexPBEdgrad_rho_b_2,dexPBEdgrad_rho_a_b)
  call g0_dg0(rho, rho_a, rho_b, g0, dg0drho)
  
  if(dabs(exPBE).lt.thr)then
   exPBE = 1.d-12
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

  gamma = exPBE / (a*n2xc_UEG)
  if(dabs(gamma).lt.thr)then
   gamma = 1.d-12
  endif

  delta = -(b*n2_UEG*gamma**2) / exPBE
  if(dabs(delta).lt.thr)then
   delta = 1.d-12
  endif
  
  denom = 1.d0 + delta*mu + gamma*(mu**2)

  ex_srmuPBE=exPBE/denom

! calculation of derivatives
  !dex/dn
  dn2_UEGdrho     = 2.d0*rho*g0 + (rho**2)*dg0drho
  dn2_UEGdrho_a   = 2.d0*rho*g0 + (rho**2)*dg0drho
  dn2_UEGdrho_b   = 2.d0*rho*g0 + (rho**2)*dg0drho
  dn2xc_UEGdrho   = dn2_UEGdrho - 2.d0*rho
  dn2xc_UEGdrho_a = dn2xc_UEGdrho
  dn2xc_UEGdrho_b = dn2xc_UEGdrho
  
  dgammadrho_a = (1.d0/(a*n2xc_UEG))*dexPBEdrho_a  - (exPBE/(a*n2xc_UEG**2))*dn2xc_UEGdrho_a
  ddeltadrho_a = -((b*gamma**2)/exPBE)*dn2_UEGdrho_a -(b*n2_UEG/exPBE)*2.d0*gamma*dgammadrho_a + (b*n2_UEG*gamma**2)/(exPBE**2)*dexPBEdrho_a
  dgammadrho_b = (1.d0/(a*n2xc_UEG))*dexPBEdrho_b  - (exPBE/(a*n2xc_UEG**2))*dn2xc_UEGdrho_b
  ddeltadrho_b = -((b*gamma**2)/exPBE)*dn2_UEGdrho_b -(b*n2_UEG/exPBE)*2.d0*gamma*dgammadrho_b + (b*n2_UEG*gamma**2)/(exPBE**2)*dexPBEdrho_b

  ddenomdrho_a = ddeltadrho_a * mu + dgammadrho_a * (mu**2)
  ddenomdrho_b = ddeltadrho_b * mu + dgammadrho_b * (mu**2)
  
  dexdrho_a = dexPBEdrho_a/denom - exPBE*ddenomdrho_a/denom**2
  dexdrho_b = dexPBEdrho_b/denom - exPBE*ddenomdrho_b/denom**2
  
!  dexdrho = 0.5d0*(dexdrho_a + dexdrho_b)
 
  !dex/d((gradn)^2)
  dgammadgrad_rho_a_2 =dexPBEdgrad_rho_a_2/(a*n2xc_UEG)
  dgammadgrad_rho_b_2 =dexPBEdgrad_rho_b_2/(a*n2xc_UEG)
  dgammadgrad_rho_a_b =dexPBEdgrad_rho_a_b/(a*n2xc_UEG)

  ddeltadgrad_rho_a_2 = ((b*n2_UEG*gamma**2)/(exPBE**2))*dexPBEdgrad_rho_a_2 - b*(n2_UEG/exPBE)*2.d0*gamma*dgammadgrad_rho_a_2
  ddeltadgrad_rho_b_2 = ((b*n2_UEG*gamma**2)/(exPBE**2))*dexPBEdgrad_rho_b_2 - b*(n2_UEG/exPBE)*2.d0*gamma*dgammadgrad_rho_b_2
  ddeltadgrad_rho_a_b = ((b*n2_UEG*gamma**2)/(exPBE**2))*dexPBEdgrad_rho_a_b - b*(n2_UEG/exPBE)*2.d0*gamma*dgammadgrad_rho_a_b

  ddenomdgrad_rho_a_2 = ddeltadgrad_rho_a_2*mu + dgammadgrad_rho_a_2*mu**2  
  ddenomdgrad_rho_b_2 = ddeltadgrad_rho_b_2*mu + dgammadgrad_rho_b_2*mu**2
  ddenomdgrad_rho_a_b =ddeltadgrad_rho_a_b*mu + dgammadgrad_rho_a_b*mu**2
  
  dexdgrad_rho_a_2 = dexPBEdgrad_rho_a_2/denom - exPBE*ddenomdgrad_rho_a_2/(denom**2)
  dexdgrad_rho_b_2 = dexPBEdgrad_rho_b_2/denom  - exPBE*ddenomdgrad_rho_b_2/(denom**2)
  dexdgrad_rho_a_b = dexPBEdgrad_rho_a_b/denom - exPBE*ddenomdgrad_rho_a_b/(denom**2)

!  dexdgrad_rho_2 = 0.25d0*(dexdgrad_rho_a_2 + dexdgrad_rho_b_2 + dexdgrad_rho_a_b)
  
!  print*, '..................................'
!  print*, 'rhoa                =', rho_a
!  print*, 'rhob                =', rho_b
!  print*, 'gradrho_a_2         =', grad_rho_a_2
!  print*, 'gradrho_b_2         =', grad_rho_b_2
!  print*, 'gradrho_a_b         =', grad_rho_a_b
!  print*, 'exPBE =', exPBE, 'ex_srmuPBE =', ex_srmuPBE
!  print*, 'dexPBEdrho_a        =', dexPBEdrho_a 
!  print*, 'dexPBEdrho_b        =', dexPBEdrho_b
!  print*, 'dexPBEdgrad_rho_a_2 =', dexPBEdgrad_rho_a_2
!  print*, 'dexPBEdgrad_rho_b_2 =', dexPBEdgrad_rho_b_2
!  print*, 'dexPBEdgrad_rho_a_b= ', dexPBEdgrad_rho_a_b
!  print*, 'dexPBEdgrad_rho_2   =', dexPBEdgrad_rho_2
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
  end subroutine exmdsrPBE
!---------------------------------------------------------------------------------------------------------------------------------------------

  subroutine ecmdsrPBE(mu,rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b,ec_srmuPBE,decdrho_a,decdrho_b, decdgrad_rho_a_2,decdgrad_rho_b_2,decdgrad_rho_a_b)

  implicit none
  BEGIN_DOC
  ! Calculation of correlation energy and chemical potential in PBE approximation using multideterminantal wave function (short-range part)
  END_DOC
 
  double precision, intent(in)  :: mu
  double precision, intent(in)  :: rho_a, rho_b, grad_rho_a_2, grad_rho_b_2, grad_rho_a_b
  double precision, intent(out) :: ec_srmuPBE, decdrho_a, decdrho_b, decdgrad_rho_a_2, decdgrad_rho_b_2, decdgrad_rho_a_b
  ! double precision              :: decdgrad_rho_2, decdrho 
  double precision              :: ecPBE, decPBEdrho_a, decPBEdrho_b, decPBEdgrad_rho_a_2, decPBEdgrad_rho_b_2,decPBEdgrad_rho_a_b
  double precision              :: rho_c, rho_o, grad_rho_c_2, grad_rho_o_2, grad_rho_o_c, decPBEdrho_c, decPBEdrho_o, decPBEdgrad_rho_c_2, decPBEdgrad_rho_o_2, decPBEdgrad_rho_c_o
  double precision              :: beta, dbetadrho_a, dbetadrho_b, dbetadgradrho_a_2, dbetadgradrho_b_2, dbetadgradrho_a_b
  double precision              :: denom, ddenomdrho_a, ddenomdrho_b, ddenomdgrad_rho_a_2, ddenomdgrad_rho_b_2, ddenomdgrad_rho_a_b
  double precision              :: pi, c, thr
  double precision              :: rho, m  
  double precision              :: n2_UEG, dn2_UEGdrho_a, dn2_UEGdrho_b
  double precision              :: g0, dg0drho
 
  if(dabs(rho_a-rho_b)/dabs(rho_a+rho_b) > 1.d-1)then
   print*,'rho_a,rho_b        = ',rho_a,rho_b
   print*,'dabs(rho_a-rho_b)  = ',dabs(rho_a-rho_b)
   stop "routine implemented only for closed-shell systems"
  endif 

  pi = dacos(-1.d0)
  rho = rho_a + rho_b
  m = rho_a - rho_b
  thr = 1.d-12
  
! correlation PBE standard and on-top pair distribution 
  call rho_ab_to_rho_oc(rho_a,rho_b,rho_o,rho_c)
  call grad_rho_ab_to_grad_rho_oc(grad_rho_a_2,grad_rho_b_2,grad_rho_a_b,grad_rho_o_2,grad_rho_c_2,grad_rho_o_c)

  call ec_pbe_sr(1.d-12,rho_c,rho_o,grad_rho_c_2,grad_rho_o_2,grad_rho_o_c,ecPBE,decPBEdrho_c,decPBEdrho_o,decPBEdgrad_rho_c_2,decPBEdgrad_rho_o_2, decPBEdgrad_rho_c_o)

  call v_rho_oc_to_v_rho_ab(decPBEdrho_o, decPBEdrho_c, decPBEdrho_a, decPBEdrho_b)
  call v_grad_rho_oc_to_v_grad_rho_ab(decPBEdgrad_rho_o_2, decPBEdgrad_rho_c_2, decPBEdgrad_rho_c_o, decPBEdgrad_rho_a_2, decPBEdgrad_rho_b_2, decPBEdgrad_rho_a_b)

  call g0_dg0(rho, rho_a, rho_b, g0, dg0drho)

! calculation of energy
  c = 2*dsqrt(pi)*(1.d0 - dsqrt(2.d0))/3.d0
  
  n2_UEG = (rho**2)*g0
  if(dabs(n2_UEG).lt.thr)then
   n2_UEG = 1.d-12
  endif  
  
  beta = ecPBE/(c*n2_UEG)
  if(dabs(beta).lt.thr)then
   beta = 1.d-12
  endif

  denom = 1.d0 + beta*mu**3
  ec_srmuPBE=ecPBE/denom

! calculation of derivatives 
  !dec/dn
  
  !n2_UEG = (na^2 + nb^2 + 2nanb)g0
  !dn2_UEGdrhoa = (2na + 2nb)g0 + (na + nb)^2 * dg0drhoa
  !dn2_UEGdrhoa = 2*rho*g0 + rho^2 *dg0drho
  dn2_UEGdrho_a = 2.d0*rho*g0 + (rho**2)*dg0drho
  dn2_UEGdrho_b = 2.d0*rho*g0 + (rho**2)*dg0drho
  
  dbetadrho_a  = decPBEdrho_a/(c*n2_UEG) - (ecPBE/(c*n2_UEG**2))*dn2_UEGdrho_a
  dbetadrho_b  = decPBEdrho_b/(c*n2_UEG) - (ecPBE/(c*n2_UEG**2))*dn2_UEGdrho_b

  ddenomdrho_a = dbetadrho_a*mu**3
  ddenomdrho_b = dbetadrho_b*mu**3
  decdrho_a = decPBEdrho_a/denom - ecPBE*ddenomdrho_a/(denom**2)
  decdrho_b = decPBEdrho_b/denom - ecPBE*ddenomdrho_b/(denom**2)
  ! decdrho = 0.5d0*(decdrho_a + decdrho_b)
  
  !dec/((dgradn)^2)
  dbetadgradrho_a_2 = decPBEdgrad_rho_a_2/(c*n2_UEG)
  dbetadgradrho_b_2 = decPBEdgrad_rho_b_2/(c*n2_UEG)
  dbetadgradrho_a_b = decPBEdgrad_rho_a_b/(c*n2_UEG)
  
  ddenomdgrad_rho_a_2 = dbetadgradrho_a_2*mu**3
  ddenomdgrad_rho_b_2 = dbetadgradrho_b_2*mu**3
  ddenomdgrad_rho_a_b = dbetadgradrho_a_b*mu**3
  
  decdgrad_rho_a_2 = decPBEdgrad_rho_a_2/denom - ecPBE*ddenomdgrad_rho_a_2/(denom**2)
  decdgrad_rho_b_2 = decPBEdgrad_rho_b_2/denom - ecPBE*ddenomdgrad_rho_b_2/(denom**2)
  decdgrad_rho_a_b = decPBEdgrad_rho_a_b/denom - ecPBE*ddenomdgrad_rho_a_b/(denom**2)

  ! decdgrad_rho_2 = 0.25d0*(decdgrad_rho_a_2 + decdgrad_rho_b_2 + decdgrad_rho_a_b)
  end subroutine ecmdsrPBE

!---------------------------------------------------------------------------------------------------------------------------------------------

subroutine exc_dexc_md_sr_PBE(mu,rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b, &
       ex_srmuPBE,dexdrho_a,dexdrho_b,dexdgrad_rho_a_2,dexdgrad_rho_b_2,dexdgrad_rho_a_b, &
       ec_srmuPBE,decdrho_a,decdrho_b,decdgrad_rho_a_2,decdgrad_rho_b_2,decdgrad_rho_a_b)
 
 implicit none
 BEGIN_DOC
 ! Give exchange and correlation energies and chemical potentials
 ! Use qp_plugins_dtraore/sr_md_energies/utils_modified_jt.irp.f's plugins : exmdsrPBE and ecmdsrPBE
 END_DOC
 double precision, intent(in)  :: mu
 double precision, intent(in)  :: rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b
 double precision, intent(out) :: ex_srmuPBE,dexdrho_a,dexdrho_b,dexdgrad_rho_a_2,dexdgrad_rho_b_2,dexdgrad_rho_a_b
 double precision, intent(out) :: ec_srmuPBE,decdrho_a,decdrho_b,decdgrad_rho_a_2,decdgrad_rho_b_2,decdgrad_rho_a_b

 call exmdsrPBE(mu,rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b,ex_srmuPBE,dexdrho_a,dexdrho_b,dexdgrad_rho_a_2,dexdgrad_rho_b_2,dexdgrad_rho_a_b)

 call ecmdsrPBE(mu,rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b,ec_srmuPBE,decdrho_a,decdrho_b,decdgrad_rho_a_2,decdgrad_rho_b_2,decdgrad_rho_a_b)

 end subroutine excmdsrPBE

