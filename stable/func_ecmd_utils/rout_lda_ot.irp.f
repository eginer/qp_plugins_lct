
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine ecmdsrLDAnn2(mu_in,rho_a,rho_b,rho2,ec_srmuLDAn,decdrho_a,decdrho_b, decdrho,decdrho2)

  implicit none
  BEGIN_DOC
  ! Calculation of correlation energy and chemical potential in LDAn approximation using multideterminantal wave function (short-range part) with exact on top pair density
  ! The on-top pair density should be normalized to N(N-1)
  END_DOC
 
  double precision, intent(in)  :: mu_in
  double precision, intent(in)  :: rho_a,rho_b,rho2
  double precision, intent(out) :: ec_srmuLDAn,decdrho_a,decdrho_b, decdrho, decdrho2!, decdrho2_a, decdrho2_b
  double precision              :: mu
  double precision              :: ecLDAn,decLDAndrho_a,decLDAndrho_b,decLDAndrho
  double precision              :: rho_c, rho_o,decLDAndrho_c,decLDAndrho_o
  double precision              :: beta, dbetadrho, denom, ddenomdrho,ddenomdrho2
  double precision              :: pi, c, thr
  double precision              :: rho, m  
 
  ec_srmuLDAn       = 0.d0 
  decdrho_a        = 0.d0
  decdrho_b        = 0.d0
  decdrho          = 0.d0
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
  
! correlation LDAn standard and on-top pair distribution 
  call rho_ab_to_rho_oc(rho_a,rho_b,rho_o,rho_c)
!  call grad_rho_ab_to_grad_rho_oc(grad_rho_a_2,grad_rho_b_2,grad_rho_a_b,grad_rho_o_2,grad_rho_2,grad_rho_o_c)

!  call ec_pbe_sr(1.d-12,rho_c,rho_o,grad_rho_2,grad_rho_o_2,grad_rho_o_c,ecLDAn,decLDAndrho_c,decLDAndrho_o,decLDAndgrad_rho_2,decLDAndgrad_rho_o_2, decLDAndgrad_rho_o)
   call ec_lda(rho_a,rho_b,ecLDAn,decLDAndrho_a, decLDAndrho_b)

  call v_rho_oc_to_v_rho_ab(decLDAndrho_o, decLDAndrho_c, decLDAndrho_a, decLDAndrho_b)
!  call v_grad_rho_oc_to_v_grad_rho_ab(decLDAndgrad_rho_o_2, decLDAndgrad_rho_2, decLDAndgrad_rho_o, decLDAndgrad_rho_a_2, decLDAndgrad_rho_b_2, decLDAndgrad_rho_a_b)

! calculation of energy
  c = 2*dsqrt(pi)*(1.d0 - dsqrt(2.d0))/3.d0
   
  if(dabs(rho2).lt.1.d-10)then
   beta = 1.d+10
  else
   beta = ecLDAn/(c*rho2)
  endif
  if(dabs(beta).lt.thr)then
   beta = 1.d-12
  endif

  denom = 1.d0 + beta*mu**3
  ec_srmuLDAn=ecLDAn/denom
  if(isnan(ec_srmuLDAn))then
   print*,'stop !!! isnan(ec_srmuLDAn)'
   print*,ecLDAn,denom
   print*,beta,mu
   print*,rho_c,rho_o 
!   print*,grad_rho_2,grad_rho_o_2,grad_rho_o_c
   stop
  endif

! calculation of derivatives 
  !dec/dn
  decLDAndrho = 0.5d0 *(decLDAndrho_a + decLDAndrho_b)

  if(dabs(rho2).lt.1.d-10)then
   dbetadrho = 1.d+10 ! - (ecLDAn/(c*rho2**2))*dn2_UEGdrho
  else
   dbetadrho = decLDAndrho/(c*rho2) ! - (ecLDAn/(c*rho2**2))*dn2_UEGdrho
  endif
  ddenomdrho = dbetadrho*mu**3

  decdrho = decLDAndrho/denom - ecLDAn*ddenomdrho/(denom**2)
  decdrho_a = decdrho
  decdrho_b = decdrho

  !dec/((dgradn)^2)
!  decLDAndgrad_rho_2 = 0.25d0 *(decLDAndgrad_rho_a_2 + decLDAndgrad_rho_b_2 + decLDAndgrad_rho_a_b) 
 
! if(dabs(rho2).lt.1.d-10)then
!  dbetadgrad_rho_2 = 1.d+10
! else
!  dbetadgrad_rho_2 = decLDAndgrad_rho_2/(c*rho2)
! endif
! ddenomdgrad_rho_2 = dbetadgrad_rho_2*mu**3
  
! decdgrad_rho_2 = decLDAndgrad_rho_2/denom - ecLDAn*ddenomdgrad_rho_2/(denom**2)
! decdgrad_rho_a_2 = decdgrad_rho_2 ! + decdgrad_n_m + decdgrad_m_2
! decdgrad_rho_b_2 = decdgrad_rho_2 ! - decdgrad_n_m + decdgrad_m_2
! decdgrad_rho_a_b = 2.d0*decdgrad_rho_2 ! - 2.d0*decdgrad_m_2 

  !dec/dn2
  
  if(dabs(rho2).gt.1.d-10)then
   ddenomdrho2 = - (mu**3)* ecLDAn/(c*rho2**2)
  else 
   ddenomdrho2 = - (mu**3)* ecLDAn * 1.d-20
  endif

  decdrho2 = - ecLDAn*ddenomdrho2/(denom**2)
  ! decdrho2_a = decdrho2
  ! decdrho2_b = decdrho2

  end subroutine ecmdsrLDAnn2


  subroutine exmdsrLDAnn2(mu_in,rho_a,rho_b,rho2,ex_srmuLDAn,dexdrho_a,dexdrho_b, dexdrho,dexdrho2)

  implicit none
  BEGIN_DOC
  ! Calculation of exchange energy and chemical potential in LDAn approximation using multideterminantal wave function (short-range part) with exact on top pair density
  ! The on-top pair density should be normalized to N(N-1)
  END_DOC
 
  double precision, intent(in)  :: mu_in
  double precision, intent(in)  :: rho_a,rho_b,rho2
  double precision, intent(out) :: ex_srmuLDAn, dexdrho_a,dexdrho_b, dexdrho, dexdrho2
  double precision              :: mu
  double precision              :: exLDAn,dexLDAndrho_a,dexLDAndrho_b,dexLDAndrho
  double precision              :: denom, ddenomdrho, ddenomdrho2
  double precision              :: delta, ddeltadrho, ddeltadrho2
  double precision              :: gamma, dgammadrho, dgammadrho2, dgamma_squaredrho
  double precision              :: pi, a, b, thr
  double precision              :: rho, m  
  double precision              :: n2xc, dn2xcdrho
  ex_srmuLDAn       = 0.d0 
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
  
! exchange LDAn standard and on-top pair distribution 
!  call ex_pbe_sr(1.d-12,rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b,exLDAn,dexLDAndrho_a,dexLDAndrho_b,dexLDAndgrad_rho_a_2,dexLDAndgrad_rho_b_2,dexLDAndgrad_rho_a_b)
  call ex_lda(rho_a,rho_b,exLDAn,dexLDAndrho_a,dexLDAndrho_b)

! calculation of energy
  a = pi / 2.d0
  b = 2.d0*dsqrt(pi)*(2.d0*dsqrt(2.d0) - 1.d0)/3.d0   
   
  n2xc = rho2 - rho**2
  if(dabs(n2xc).lt.thr)then
   n2xc = 1.d-12
  endif

  gamma = exLDAn / (a*n2xc)
  if(dabs(gamma).lt.thr)then
   gamma = 1.d-12
  endif

  delta = -(b*rho2*gamma**2) / exLDAn
  if(dabs(delta).lt.thr)then
   delta = 1.d-12
  endif

  denom = 1.d0 + delta*mu + gamma*(mu**2)
  ex_srmuLDAn=exLDAn/denom
  if(isnan(ex_srmuLDAn))then
   print*,'stop !!! isnan(ex_srmuLDAn)'
   print*,exLDAn,denom
   print*,gamma,delta,mu
   print*,rho_a,rho_b 
!   print*,grad_rho_2,grad_rho_a_2,grad_rho_a_b
   stop
  endif

! calculation of derivatives 
  !dex/dn
  dexLDAndrho = 0.5d0 *(dexLDAndrho_a + dexLDAndrho_b)

  dn2xcdrho = -2.d0*rho 
  if(dabs(rho2).lt.1.d-10)then
   dgammadrho = 1.d+10 ! - (ecLDAn/(c*rho2**2))*dn2_UEGdrho
  else
   dgammadrho = dexLDAndrho/(a*n2xc) - (exLDAn*dn2xcdrho)/(a*n2xc**2)
  endif

  dgamma_squaredrho= 2.d0*gamma*dgammadrho
  if(dabs(rho2).lt.1.d-10)then
   ddeltadrho = 1.d+10 ! - (ecLDAn/(c*rho2**2))*dn2_UEGdrho
  else
   ddeltadrho = -(b*rho2/exLDAn)*dgamma_squaredrho + ((b*rho2*gamma**2)/(exLDAn**2))*(dexLDAndrho)
  endif

  ddenomdrho = ddeltadrho*mu + dgammadrho*mu**2

  dexdrho = dexLDAndrho/denom - exLDAn*ddenomdrho/(denom**2)
  dexdrho_a = dexdrho
  dexdrho_b = dexdrho

! !dec/((dgradn)^2)
! dexLDAndgrad_rho_2 = 0.25d0 *(dexLDAndgrad_rho_a_2 + dexLDAndgrad_rho_b_2 + dexLDAndgrad_rho_a_b) 
!
! if(dabs(rho2).lt.1.d-10)then
!  dgammadgrad_rho_2 = 1.d+10
! else
!  dgammadgrad_rho_2 = dexLDAndgrad_rho_2/(a*n2xc)
! endif

! dgamma_squaredgrad_rho_2 = 2.d0*gamma*dgammadgrad_rho_2 
! if(dabs(rho2).lt.1.d-10)then
!  ddeltadgrad_rho_2 = 1.d+10
! else
!  ddeltadgrad_rho_2 = -(b*rho2/exLDAn)*dgamma_squaredgrad_rho_2 + ((b*rho2*gamma**2)/(exLDAn**2))*(dexLDAndgrad_rho_2)
! endif

! ddenomdgrad_rho_2 = ddeltadgrad_rho_2*mu + dgammadgrad_rho_2*mu**2
! 
! dexdgrad_rho_2 = dexLDAndgrad_rho_2/denom - exLDAn*ddenomdgrad_rho_2/(denom**2)
! dexdgrad_rho_a_2 = dexdgrad_rho_2 ! + dexdgrad_n_m + dexdgrad_m_2
! dexdgrad_rho_b_2 = dexdgrad_rho_2 ! - dexdgrad_n_m + dexdgrad_m_2
! dexdgrad_rho_a_b = 2.d0*dexdgrad_rho_2 ! - 2.d0*dexdgrad_m_2

  !dec/dn2
  
  dgammadrho2 = exLDAn/(a*n2xc**2) 
  ddeltadrho2 = (b*gamma**2)/(exLDAn) + (b*dgammadrho2*2.d0*gamma)/(exLDAn) 
  if(dabs(rho2).gt.1.d-10)then
   ddenomdrho2 = - (mu)* dgammadrho2 - (mu**2)*ddeltadrho2
  else 
   ddenomdrho2 = - (mu**3)* exLDAn * 1.d-20
  endif

  dexdrho2 = - exLDAn*ddenomdrho2/(denom**2)
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
!    print*, 'exLDAn =', exLDAn, 'ex_srmuLDAn =', ex_srmuLDAn
!    print*, 'dexLDAndrho_a        =', dexLDAndrho_a 
!    print*, 'dexLDAndrho_b        =', dexLDAndrho_b
!   ! print*, 'dexLDAndgrad_rho_a_2 =', dexLDAndgrad_rho_a_2
!   ! print*, 'dexLDAndgrad_rho_b_2 =', dexLDAndgrad_rho_b_2
!   ! print*, 'dexLDAndgrad_rho_a_b= ', dexLDAndgrad_rho_a_b
!    print*, 'dexLDAndgrad_rho_2   =', dexLDAndgrad_rho_2
!    print*, 'dexdrho_a           =', dexdrho_a          
!    print*, 'dexdrho_b           =', dexdrho_b          
!    print*, 'dexdgrad_rho_a_2    =', dexdgrad_rho_a_2   
!    print*, 'dexdgrad_rho_b_2    =', dexdgrad_rho_b_2   
!    print*, 'dexdgrad_rho_a_b    =', dexdgrad_rho_a_b   
!    print*, 'denom               =', denom
!    print*, '------> dexdgradrho =', dexdgrad_rho_2   
!    print*, 'dexdrho2            =', dexdrho2
!    print*, '..................................'

  end subroutine exmdsrLDAnn2
