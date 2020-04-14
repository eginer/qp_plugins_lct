
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine ecmdsrPBEn2(mu,rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b,rho2,ec_srmuPBE,decdrho_a,decdrho_b, decdrho, decdgrad_rho_a_2,decdgrad_rho_b_2,decdgrad_rho_a_b, decdgrad_rho_2,decdrho2)

  implicit none
  BEGIN_DOC
  ! Calculation of correlation energy and chemical potential in PBE approximation using multideterminantal wave function (short-range part) with exact on top pair density
  END_DOC
 
  double precision, intent(in)  :: mu
  double precision, intent(in)  :: rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b, rho2
  double precision, intent(out) :: ec_srmuPBE,decdrho_a,decdrho_b,decdgrad_rho_a_2,decdgrad_rho_b_2,decdgrad_rho_a_b, decdgrad_rho_2, decdrho, decdrho2!, decdrho2_a, decdrho2_b
  double precision              :: ecPBE,decPBEdrho_a,decPBEdrho_b,decPBEdrho, decPBEdgrad_rho_a_2,decPBEdgrad_rho_b_2,decPBEdgrad_rho_a_b
  double precision              :: rho_c, rho_o,grad_rho_2,grad_rho_o_2,grad_rho_o_c,decPBEdrho_c,decPBEdrho_o,decPBEdgrad_rho_2,decPBEdgrad_rho_o_2, decPBEdgrad_rho_o
  double precision              :: beta, dbetadrho, dbetadgrad_rho_2, denom, ddenomdrho, ddenomdgrad_rho_2, ddenomdrho2
  double precision              :: pi, c, thr
  double precision              :: rho, m  
 
  if(abs(rho_a-rho_b) > 1.d-12)then
  stop "routine implemented only for closed-shell systems"
  endif 

  pi = dacos(-1.d0)
  rho = rho_a + rho_b
  m = rho_a - rho_b
  thr = 1.d-12
  
! correlation PBE standard and on-top pair distribution 
  call rho_ab_to_rho_oc(rho_a,rho_b,rho_o,rho_c)
  call grad_rho_ab_to_grad_rho_oc(grad_rho_a_2,grad_rho_b_2,grad_rho_a_b,grad_rho_o_2,grad_rho_2,grad_rho_o_c)

  call ec_pbe_sr(1.d-12,rho_c,rho_o,grad_rho_2,grad_rho_o_2,grad_rho_o_c,ecPBE,decPBEdrho_c,decPBEdrho_o,decPBEdgrad_rho_2,decPBEdgrad_rho_o_2, decPBEdgrad_rho_o)

  call v_rho_oc_to_v_rho_ab(decPBEdrho_o, decPBEdrho_c, decPBEdrho_a, decPBEdrho_b)
  call v_grad_rho_oc_to_v_grad_rho_ab(decPBEdgrad_rho_o_2, decPBEdgrad_rho_2, decPBEdgrad_rho_o, decPBEdgrad_rho_a_2, decPBEdgrad_rho_b_2, decPBEdgrad_rho_a_b)

! calculation of energy
  c = 2*dsqrt(pi)*(1.d0 - dsqrt(2.d0))/3.d0
   
  beta = ecPBE/(c*rho2)
  if(dabs(beta).lt.thr)then
   beta = 1.d-12
  endif

  denom = 1.d0 + beta*mu**3
  ec_srmuPBE=ecPBE/denom

! calculation of derivatives 
  !dec/dn
  decPBEdrho = 0.5d0 *(decPBEdrho_a + decPBEdrho_b)

  dbetadrho = decPBEdrho/(c*rho2) ! - (ecPBE/(c*rho2**2))*dn2_UEGdrho
  ddenomdrho = dbetadrho*mu**3

  decdrho = decPBEdrho/denom - ecPBE*ddenomdrho/(denom**2)
  decdrho_a = decdrho
  decdrho_b = decdrho

  !dec/((dgradn)^2)
  decPBEdgrad_rho_2 = 0.25d0 *(decPBEdgrad_rho_a_2 + decPBEdgrad_rho_b_2 + decPBEdgrad_rho_a_b) 
 
  dbetadgrad_rho_2 = decPBEdgrad_rho_2/(c*rho2)
  ddenomdgrad_rho_2 = dbetadgrad_rho_2*mu**3
  
  decdgrad_rho_2 = decPBEdgrad_rho_2/denom - ecPBE*ddenomdgrad_rho_2/(denom**2)
  decdgrad_rho_a_2 = decdgrad_rho_2 ! + decdgrad_n_m + decdgrad_m_2
  decdgrad_rho_b_2 = decdgrad_rho_2 ! - decdgrad_n_m + decdgrad_m_2
  decdgrad_rho_a_b = 2.d0*decdgrad_rho_2 ! - 2.d0*decdgrad_m_2 

  !dec/dn2
  
  ddenomdrho2 = - (mu**3)* ecPBE/(c*rho2**2)

  decdrho2 = - ecPBE*ddenomdrho2/(denom**2)
  ! decdrho2_a = decdrho2
  ! decdrho2_b = decdrho2

  end subroutine ecmdsrPBE

