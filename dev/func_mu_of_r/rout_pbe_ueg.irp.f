!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine g0_dg0(rho, rho_a, rho_b, g0, dg0drho)

  implicit none                                                                                                                              
  BEGIN_DOC
  ! Give the on-top pair distribution function g0 and its derivative according to rho dg0drho
  END_DOC

  double precision, intent (in) :: rho, rho_a, rho_b
  double precision, intent (out) :: g0, dg0drho
  double precision :: pi
  double precision :: g0_UEG_mu_inf, dg0drs
  double precision :: C1, F1, D1, E1, B1, rs

  pi = dacos(-1.d0)
  C1 = 0.0819306d0
  F1 = 0.752411d0
  D1 = -0.0127713d0
  E1 = 0.00185898d0
  B1 = 0.7317d0 - F1
  rs = (3.d0 / (4.d0*pi*rho))**(1.d0/3.d0) 
   
  g0 = g0_UEG_mu_inf(rho_a, rho_b)
  dg0drs = 0.5d0*((-B1 + 2.d0*C1*rs + 3.d0*D1*rs**2 + 4.d0*E1*rs**3)-F1*(1.d0 - B1*rs + C1*rs**2 + D1*rs**3 + E1*rs**4))*exp(-F1*rs)
  dg0drho = -((6.d0*dsqrt(pi)*rho**2)**(-2.d0/3.d0))*dg0drs
 
  end subroutine g0_dg0



  subroutine ecmdsrPBE(mu,rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b,ec_srmuPBE,decdrho_a,decdrho_b, decdrho, decdgrad_rho_a_2,decdgrad_rho_b_2,decdgrad_rho_a_b, decdgrad_rho_2)

  implicit none
  BEGIN_DOC
  ! Calculation of correlation energy and chemical potential in PBE approximation using multideterminantal wave function (short-range part)
  END_DOC
 
  double precision, intent(in)  :: mu
  double precision, intent(in)  :: rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b
  double precision, intent(out) :: ec_srmuPBE,decdrho_a,decdrho_b,decdgrad_rho_a_2,decdgrad_rho_b_2,decdgrad_rho_a_b, decdgrad_rho_2, decdrho
  double precision              :: ecPBE,decPBEdrho_a,decPBEdrho_b,decPBEdgrad_rho_2,decPBEdrho, decPBEdgrad_rho_a_2,decPBEdgrad_rho_b_2,decPBEdgrad_rho_a_b
  double precision              :: rho_c, rho_o,grad_rho_c_2,grad_rho_o_2,grad_rho_o_c,decPBEdrho_c,decPBEdrho_o,decPBEdgrad_rho_c_2,decPBEdgrad_rho_o_2, decPBEdgrad_rho_c_o
  double precision              :: beta, dbetadrho, dbetadgrad_rho_2, denom, ddenomdrho, ddenomdgrad_rho_2
  double precision              :: pi, c, thr
  double precision              :: rho, m  
  double precision              :: n2_UEG, dn2_UEGdrho, n2xc_UEG, dn2xc_UEGdrho, g0, dg0drho
 
  if(abs(rho_a-rho_b) > 1.d-12)then
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
  decPBEdrho = 0.5d0 *(decPBEdrho_a + decPBEdrho_b)

  dn2_UEGdrho = 2.d0*rho*g0 + (rho**2)*dg0drho

  dbetadrho  = decPBEdrho/(c*n2_UEG) - (ecPBE/(c*n2_UEG**2))*dn2_UEGdrho
  ddenomdrho = dbetadrho*mu**3

  decdrho   = decPBEdrho/denom - ecPBE*ddenomdrho/(denom**2)
  decdrho_a = decdrho
  decdrho_b = decdrho

  !dec/((dgradn)^2)
 
  decPBEdgrad_rho_2 = 0.25d0 *(decPBEdgrad_rho_a_2 + decPBEdgrad_rho_b_2 + decPBEdgrad_rho_a_b) 
 
  dbetadgrad_rho_2  = decPBEdgrad_rho_2/(c*n2_UEG)
  ddenomdgrad_rho_2 = dbetadgrad_rho_2*mu**3
  
  decdgrad_rho_2   = decPBEdgrad_rho_2/denom - ecPBE*ddenomdgrad_rho_2/(denom**2)
  decdgrad_rho_a_2 = decdgrad_rho_2 ! + decdgrad_n_m + decdgrad_m_2
  decdgrad_rho_b_2 = decdgrad_rho_2 ! - decdgrad_n_m + decdgrad_m_2
  decdgrad_rho_a_b = 2.d0*decdgrad_rho_2 ! - 2.d0*decdgrad_m_2 
  end subroutine ecmdsrPBE

!---------------------------------------------------------------------------------------------------------------------------------------------

