
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine ecmdsrPBEn2(mu_in,rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b,rho2,ec_srmuPBE,decdrho_a,decdrho_b, decdrho, decdgrad_rho_a_2,decdgrad_rho_b_2,decdgrad_rho_a_b, decdgrad_rho_2,decdrho2)

  implicit none
  BEGIN_DOC
  ! Calculation of correlation energy and chemical potential in PBE approximation using multideterminantal wave function (short-range part) with exact on top pair density
  ! The on-top pair density should be normalized to N(N-1)
  END_DOC
 
  double precision, intent(in)  :: mu_in
  double precision, intent(in)  :: rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b, rho2
  double precision, intent(out) :: ec_srmuPBE,decdrho_a,decdrho_b,decdgrad_rho_a_2,decdgrad_rho_b_2,decdgrad_rho_a_b, decdgrad_rho_2, decdrho, decdrho2!, decdrho2_a, decdrho2_b
  double precision              :: mu
  double precision              :: ecPBE,decPBEdrho_a,decPBEdrho_b,decPBEdrho, decPBEdgrad_rho_a_2,decPBEdgrad_rho_b_2,decPBEdgrad_rho_a_b
  double precision              :: rho_c, rho_o,grad_rho_2,grad_rho_o_2,grad_rho_o_c,decPBEdrho_c,decPBEdrho_o,decPBEdgrad_rho_2,decPBEdgrad_rho_o_2, decPBEdgrad_rho_o
  double precision              :: beta, dbetadrho, dbetadgrad_rho_2, denom, ddenomdrho, ddenomdgrad_rho_2, ddenomdrho2
  double precision              :: pi, c, thr
  double precision              :: rho, m  
 
  ec_srmuPBE       = 0.d0 
  decdrho_a        = 0.d0
  decdrho_b        = 0.d0
  decdrho          = 0.d0
  decdgrad_rho_a_2 = 0.d0 
  decdgrad_rho_b_2 = 0.d0
  decdgrad_rho_a_b = 0.d0  
  decdgrad_rho_2   = 0.d0
  decdrho2         = 0.d0
  rho = rho_a + rho_b
  if(rho.lt.1.d-10)then
   return
  else if(rho2/(rho**2) .lt. 1.d-6)then
   return
  endif
  m = rho_a - rho_b

  pi = dacos(-1.d0)
  thr = 1.d-12
  mu = min(mu_in,1.d+10)
  
! correlation PBE standard and on-top pair distribution 
  call rho_ab_to_rho_oc(rho_a,rho_b,rho_o,rho_c)
  call grad_rho_ab_to_grad_rho_oc(grad_rho_a_2,grad_rho_b_2,grad_rho_a_b,grad_rho_o_2,grad_rho_2,grad_rho_o_c)

  call ec_pbe_sr(1.d-12,rho_c,rho_o,grad_rho_2,grad_rho_o_2,grad_rho_o_c,ecPBE,decPBEdrho_c,decPBEdrho_o,decPBEdgrad_rho_2,decPBEdgrad_rho_o_2, decPBEdgrad_rho_o)

  call v_rho_oc_to_v_rho_ab(decPBEdrho_o, decPBEdrho_c, decPBEdrho_a, decPBEdrho_b)
  call v_grad_rho_oc_to_v_grad_rho_ab(decPBEdgrad_rho_o_2, decPBEdgrad_rho_2, decPBEdgrad_rho_o, decPBEdgrad_rho_a_2, decPBEdgrad_rho_b_2, decPBEdgrad_rho_a_b)

! calculation of energy
  c = 2*dsqrt(pi)*(1.d0 - dsqrt(2.d0))/3.d0
   
  if(dabs(rho2).lt.1.d-10)then
   beta = 1.d+10
  else
   beta = ecPBE/(c*rho2)
  endif
  if(dabs(beta).lt.thr)then
   beta = 1.d-12
  endif

  denom = 1.d0 + beta*mu**3
  ec_srmuPBE=ecPBE/denom
  if(isnan(ec_srmuPBE))then
   print*,'stop !!! isnan(ec_srmuPBE)'
   print*,ecPBE,denom
   print*,beta,mu
   print*,rho_c,rho_o 
   print*,grad_rho_2,grad_rho_o_2,grad_rho_o_c
   stop
  endif

! calculation of derivatives 
  !dec/dn
  decPBEdrho = 0.5d0 *(decPBEdrho_a + decPBEdrho_b)

  if(dabs(rho2).lt.1.d-10)then
   dbetadrho = 1.d+10 ! - (ecPBE/(c*rho2**2))*dn2_UEGdrho
  else
   dbetadrho = decPBEdrho/(c*rho2) ! - (ecPBE/(c*rho2**2))*dn2_UEGdrho
  endif
  ddenomdrho = dbetadrho*mu**3

  decdrho = decPBEdrho/denom - ecPBE*ddenomdrho/(denom**2)
  decdrho_a = decdrho
  decdrho_b = decdrho

  !dec/((dgradn)^2)
  decPBEdgrad_rho_2 = 0.25d0 *(decPBEdgrad_rho_a_2 + decPBEdgrad_rho_b_2 + decPBEdgrad_rho_a_b) 
 
  if(dabs(rho2).lt.1.d-10)then
   dbetadgrad_rho_2 = 1.d+10
  else
   dbetadgrad_rho_2 = decPBEdgrad_rho_2/(c*rho2)
  endif
  ddenomdgrad_rho_2 = dbetadgrad_rho_2*mu**3
  
  decdgrad_rho_2 = decPBEdgrad_rho_2/denom - ecPBE*ddenomdgrad_rho_2/(denom**2)
  decdgrad_rho_a_2 = decdgrad_rho_2 ! + decdgrad_n_m + decdgrad_m_2
  decdgrad_rho_b_2 = decdgrad_rho_2 ! - decdgrad_n_m + decdgrad_m_2
  decdgrad_rho_a_b = 2.d0*decdgrad_rho_2 ! - 2.d0*decdgrad_m_2 

  !dec/dn2
  
  if(dabs(rho2).gt.1.d-10)then
   ddenomdrho2 = - (mu**3)* ecPBE/(c*rho2**2)
  else 
   ddenomdrho2 = - (mu**3)* ecPBE * 1.d-20
  endif

  decdrho2 = - ecPBE*ddenomdrho2/(denom**2)
  ! decdrho2_a = decdrho2
  ! decdrho2_b = decdrho2

  end subroutine ecmdsrPBEn2


  subroutine exmdsrPBEn2(mu_in,rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b,rho2,ex_srmuPBE,dexdrho_a,dexdrho_b, dexdrho, dexdgrad_rho_a_2,dexdgrad_rho_b_2,dexdgrad_rho_a_b, dexdgrad_rho_2,dexdrho2)

  implicit none
  BEGIN_DOC
  ! Calculation of exchange energy and chemical potential in PBE approximation using multideterminantal wave function (short-range part) with exact on top pair density
  ! The on-top pair density should be normalized to N(N-1)
  END_DOC
 
  double precision, intent(in)  :: mu_in
  double precision, intent(in)  :: rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b, rho2
  double precision, intent(out) :: ex_srmuPBE, dexdrho_a,dexdrho_b,dexdgrad_rho_a_2,dexdgrad_rho_b_2,dexdgrad_rho_a_b, dexdgrad_rho_2, dexdrho, dexdrho2
  double precision              :: mu
  double precision              :: exPBE,dexPBEdrho_a,dexPBEdrho_b,dexPBEdrho, dexPBEdgrad_rho_a_2,dexPBEdgrad_rho_b_2,dexPBEdgrad_rho_a_b,dexPBEdgrad_rho_2
  double precision              :: denom, ddenomdrho, ddenomdgrad_rho_2, ddenomdrho2
  double precision              :: delta, ddeltadrho, ddeltadgrad_rho_2, ddeltadrho2
  double precision              :: gamma, dgammadrho, dgammadgrad_rho_2, dgammadrho2, dgamma_squaredrho, dgamma_squaredgrad_rho_2
  double precision              :: grad_rho_2
  double precision              :: pi, a, b, thr
  double precision              :: rho, m  
  double precision              :: n2xc, dn2xcdrho
  ex_srmuPBE       = 0.d0 
  rho = rho_a + rho_b
  if(rho.lt.1.d-10)then
   return
  else if(rho2/(rho**2) .lt. 1.d-6)then
   return
  endif
  m = rho_a - rho_b

  pi = dacos(-1.d0)
  thr = 1.d-12
  mu = min(mu_in,1.d+10)
  
! exchange PBE standard and on-top pair distribution 
  call ex_pbe_sr(1.d-12,rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b,exPBE,dexPBEdrho_a,dexPBEdrho_b,dexPBEdgrad_rho_a_2,dexPBEdgrad_rho_b_2,dexPBEdgrad_rho_a_b)

! calculation of energy
  a = pi / 2.d0
  b = 2.d0*dsqrt(pi)*(2.d0*dsqrt(2.d0) - 1.d0)/3.d0   
   
  n2xc = rho2 - rho**2
  if(dabs(n2xc).lt.thr)then
   n2xc = 1.d-12
  endif

  gamma = exPBE / (a*n2xc)
  if(dabs(gamma).lt.thr)then
   gamma = 1.d-12
  endif

  delta = -(b*rho2*gamma**2) / exPBE
  if(dabs(delta).lt.thr)then
   delta = 1.d-12
  endif

  denom = 1.d0 + delta*mu + gamma*(mu**2)
  ex_srmuPBE=exPBE/denom
  if(isnan(ex_srmuPBE))then
   print*,'stop !!! isnan(ex_srmuPBE)'
   print*,exPBE,denom
   print*,gamma,delta,mu
   print*,rho_a,rho_b 
   print*,grad_rho_2,grad_rho_a_2,grad_rho_a_b
   stop
  endif

! calculation of derivatives 
  !dex/dn
  dexPBEdrho = 0.5d0 *(dexPBEdrho_a + dexPBEdrho_b)

  dn2xcdrho = -2.d0*rho 
  if(dabs(rho2).lt.1.d-10)then
   dgammadrho = 1.d+10 ! - (ecPBE/(c*rho2**2))*dn2_UEGdrho
  else
   dgammadrho = dexPBEdrho/(a*n2xc) - (exPBE*dn2xcdrho)/(a*n2xc**2)
  endif

  dgamma_squaredrho= 2.d0*gamma*dgammadrho
  if(dabs(rho2).lt.1.d-10)then
   ddeltadrho = 1.d+10 ! - (ecPBE/(c*rho2**2))*dn2_UEGdrho
  else
   ddeltadrho = -(b*rho2/exPBE)*dgamma_squaredrho + ((b*rho2*gamma**2)/(exPBE**2))*(dexPBEdrho)
  endif

  ddenomdrho = ddeltadrho*mu + dgammadrho*mu**2

  dexdrho = dexPBEdrho/denom - exPBE*ddenomdrho/(denom**2)
  dexdrho_a = dexdrho
  dexdrho_b = dexdrho

  !dec/((dgradn)^2)
  dexPBEdgrad_rho_2 = 0.25d0 *(dexPBEdgrad_rho_a_2 + dexPBEdgrad_rho_b_2 + dexPBEdgrad_rho_a_b) 
 
  if(dabs(rho2).lt.1.d-10)then
   dgammadgrad_rho_2 = 1.d+10
  else
   dgammadgrad_rho_2 = dexPBEdgrad_rho_2/(a*n2xc)
  endif

  dgamma_squaredgrad_rho_2 = 2.d0*gamma*dgammadgrad_rho_2 
  if(dabs(rho2).lt.1.d-10)then
   ddeltadgrad_rho_2 = 1.d+10
  else
   ddeltadgrad_rho_2 = -(b*rho2/exPBE)*dgamma_squaredgrad_rho_2 + ((b*rho2*gamma**2)/(exPBE**2))*(dexPBEdgrad_rho_2)
  endif

  ddenomdgrad_rho_2 = ddeltadgrad_rho_2*mu + dgammadgrad_rho_2*mu**2
  
  dexdgrad_rho_2 = dexPBEdgrad_rho_2/denom - exPBE*ddenomdgrad_rho_2/(denom**2)
  dexdgrad_rho_a_2 = dexdgrad_rho_2 ! + dexdgrad_n_m + dexdgrad_m_2
  dexdgrad_rho_b_2 = dexdgrad_rho_2 ! - dexdgrad_n_m + dexdgrad_m_2
  dexdgrad_rho_a_b = 2.d0*dexdgrad_rho_2 ! - 2.d0*dexdgrad_m_2

  !dec/dn2
  
  dgammadrho2 = exPBE/(a*n2xc**2) 
  ddeltadrho2 = (b*gamma**2)/(exPBE) + (b*dgammadrho2*2.d0*gamma)/(exPBE) 
  if(dabs(rho2).gt.1.d-10)then
   ddenomdrho2 = - (mu)* dgammadrho2 - (mu**2)*ddeltadrho2
  else 
   ddenomdrho2 = - (mu**3)* exPBE * 1.d-20
  endif

  dexdrho2 = - exPBE*ddenomdrho2/(denom**2)
  ! dexdrho2_a = dexdrho2
  ! dexdrho2_b = dexdrho2

!    print*, '..................................'
!    print*, 'rhoa                =', rho_a
!   print*, 'rhob                =', rho_b
!    print*, 'rho2                =', rho2
!    print*, 'mu                  =', mu_in
!    print*, 'gradrho_a_2         =', grad_rho_a_2
!    print*, 'gradrho_b_2         =', grad_rho_b_2
!    print*, 'gradrho_a_b         =', grad_rho_a_b
!    print*, 'exPBE =', exPBE, 'ex_srmuPBE =', ex_srmuPBE
!    print*, 'dexPBEdrho_a        =', dexPBEdrho_a 
!    print*, 'dexPBEdrho_b        =', dexPBEdrho_b
!   ! print*, 'dexPBEdgrad_rho_a_2 =', dexPBEdgrad_rho_a_2
!   ! print*, 'dexPBEdgrad_rho_b_2 =', dexPBEdgrad_rho_b_2
!   ! print*, 'dexPBEdgrad_rho_a_b= ', dexPBEdgrad_rho_a_b
!    print*, 'dexPBEdgrad_rho_2   =', dexPBEdgrad_rho_2
!    print*, 'dexdrho_a           =', dexdrho_a          
!    print*, 'dexdrho_b           =', dexdrho_b          
!    print*, 'dexdgrad_rho_a_2    =', dexdgrad_rho_a_2   
!    print*, 'dexdgrad_rho_b_2    =', dexdgrad_rho_b_2   
!    print*, 'dexdgrad_rho_a_b    =', dexdgrad_rho_a_b   
!    print*, 'denom               =', denom
!    print*, '------> dexdgradrho =', dexdgrad_rho_2   
!    print*, 'dexdrho2            =', dexdrho2
!    print*, '..................................'

  end subroutine exmdsrPBEn2
