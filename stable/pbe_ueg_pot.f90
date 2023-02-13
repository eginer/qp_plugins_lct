  subroutine ecmdpbeueg_autonom(mu,rho_c,rho_o,ecPBE,vec_depbedrho,vec_depbedgradrho, &
        ec_srmuPBE,vec_decmdpbeuegdrho,vec_decmdpbeuegdgradrho)

  implicit none
  ! Computes the ecmdpbeUEG and its derivative with respect to the density and gradients  
  !
  ! input : mu :: value of the range separation parameter 
  !         rho_c, rho_o :: total density (rho_a+rho_b) and spin density (rho_a - rho_b)
  !         ecPBE :: PBE correlation energy (epsilon_PBE) TIMES THE DENSITY 
  !         vec_depbedrho(1) : decpbe/d rho_c, vec_depbedrho(2) : decpbe/d rho_o
  !         vec_depbedgradrho(1) = decpbe/d grad_rho_2, vec_depbedgradrho(2) = decpbe/d grad_rho_b, 
  !         vec_depbedgradrho(3) = decpbe/d (grad_rho_a. grad_rho_b)
  !
  ! output : ec_srmuPBE :: value of the PBE-UEG correlation energy TIMES the density 
  !          vec_decmdpbeuegdrho(1) = d ec_srmuPBE /d rho_c, vec_decmdpbeuegdrho(2) = d ec_srmuPBE /d rho_o
  !          vec_decmdpbeuegdgradrho(1) = d ec_srmuPBE / d grad_rho_a_2
  !          vec_decmdpbeuegdgradrho(2) = d ec_srmuPBE / d grad_rho_b_2
  !          vec_decmdpbeuegdgradrho(3) = d ec_srmuPBE / d (grad_rho_a. grad_rho_b)
 
  double precision, intent(in)  :: mu, rho_c, rho_o
  double precision, intent(in)  :: ecPBE, vec_depbedrho(2), vec_depbedgradrho(3)
  double precision, intent(out) :: ec_srmuPBE, vec_decmdpbeuegdrho(2), vec_decmdpbeuegdgradrho(3)

  double precision  :: decmduegdgrad_rho_c_2, decmduegdgrad_rho_o_2, decmduegdgrad_rho_c_o
  double precision  :: rho_a, rho_b, decPBEdrho_a, decPBEdrho_b, decPBEdgrad_rho_a_2
  double precision  :: decPBEdgrad_rho_b_2, decPBEdgrad_rho_a_b
  double precision  :: beta, dbetadrho_a, dbetadrho_b
  double precision  :: dbetadgradrho_a_2, dbetadgradrho_b_2, dbetadgradrho_a_b
  double precision  :: denom, ddenomdrho_a, ddenomdrho_b, ddenomdgrad_rho_a_2
  double precision  :: ddenomdgrad_rho_b_2, ddenomdgrad_rho_a_b
  double precision  :: pi, c, thr,decPBEdgrad_rho_o_2,decPBEdgrad_rho_c_o  
  double precision  :: n2_UEG, dn2_UEGdrho_a, dn2_UEGdrho_b
  double precision  :: decPBEdrho_c,decPBEdrho_o,decPBEdgrad_rho_c_2
  double precision  :: decdgrad_rho_a_2, decdgrad_rho_b_2, decdgrad_rho_a_b, decdrho_a,decdrho_b
  double precision  :: decmduegdrho_c, decmduegdrho_o
  double precision  :: g0, dg0drho,dn2_UEGdrhoc,dn2_UEGdrhoo
  pi = dacos(-1.d0)
  decPBEdrho_c        = vec_depbedrho(1) 
  decPBEdrho_o        = vec_depbedrho(2)
  decPBEdgrad_rho_c_2 = vec_depbedgradrho(1)
  decPBEdgrad_rho_o_2 = vec_depbedgradrho(2)
  decPBEdgrad_rho_c_o = vec_depbedgradrho(3)
  
  !! closed/open --> alpha-beta
  rho_a= 0.5d0*(rho_c+rho_o)
  rho_b= 0.5d0*(rho_c-rho_o)
  decPBEdrho_a = decPBEdrho_c + decPBEdrho_o
  decPBEdrho_b = decPBEdrho_c - decPBEdrho_o
  decPBEdgrad_rho_a_2 = decPBEdgrad_rho_o_2 + decPBEdgrad_rho_c_2 + decPBEdgrad_rho_c_o
  decPBEdgrad_rho_b_2 = decPBEdgrad_rho_o_2 + decPBEdgrad_rho_c_2 - decPBEdgrad_rho_c_o
  decPBEdgrad_rho_a_b = -2.d0 * decPBEdgrad_rho_o_2 +  2.d0 * decPBEdgrad_rho_c_2
  thr = 1.d-12
  call g0_dg0_coul_UEG(rho_c, g0, dg0drho)

! calculation of energy
  c = 2*dsqrt(pi)*(1.d0 - dsqrt(2.d0))/3.d0
  n2_UEG = (rho_c**2)*g0
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
  
  !n2_UEG = (rho_c**2 - rho_o**2 ) * g0(n) 
  !dn2_UEGdrhoa = (2na + 2nb)g0 + (na + nb)^2 * dg0drhoa
  !dn2_UEGdrhoa = 2*rho*g0 + rho^2 *dg0drho
  dn2_UEGdrhoc = 2.d0 * rho_c * g0 + dg0drho * (rho_c**2 - rho_o**2)
  dn2_UEGdrhoo = -2.d0 * g0 * rho_o
  dn2_UEGdrho_a = dn2_UEGdrhoc + dn2_UEGdrhoo 
  dn2_UEGdrho_b = dn2_UEGdrhoc - dn2_UEGdrhoo 
  
  dbetadrho_a  = decPBEdrho_a/(c*n2_UEG) - ecPBE*dn2_UEGdrho_a/(c*n2_UEG**2)
  dbetadrho_b  = decPBEdrho_b/(c*n2_UEG) - ecPBE*dn2_UEGdrho_b/(c*n2_UEG**2)

  ddenomdrho_a = dbetadrho_a*mu**3
  ddenomdrho_b = dbetadrho_b*mu**3
  decdrho_a = decPBEdrho_a/denom - ecPBE*ddenomdrho_a/(denom**2)
  decdrho_b = decPBEdrho_b/denom - ecPBE*ddenomdrho_b/(denom**2)
  
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
  
  ! convert back to rho_c, rho_o variables 
  decmduegdrho_c = 0.5d0 * (decdrho_a + decdrho_b)
  decmduegdrho_o = 0.5d0 * (decdrho_a - decdrho_b)
  vec_decmdpbeuegdrho(1) = decmduegdrho_c 
  vec_decmdpbeuegdrho(2) = decmduegdrho_o 

  decmduegdgrad_rho_c_2 = 0.25d0 * (decdgrad_rho_a_2 + decdgrad_rho_b_2 + decdgrad_rho_a_b)
  decmduegdgrad_rho_o_2 = 0.25d0 * (decdgrad_rho_a_2 + decdgrad_rho_b_2 - decdgrad_rho_a_b)
  decmduegdgrad_rho_c_o = 0.25d0 * (2d0 * decdgrad_rho_a_2 - 2d0 * decdgrad_rho_b_2       )
  vec_decmdpbeuegdgradrho(1) = decmduegdgrad_rho_c_2
  vec_decmdpbeuegdgradrho(2) = decmduegdgrad_rho_o_2
  vec_decmdpbeuegdgradrho(3) = decmduegdgrad_rho_c_o

  end 

  subroutine g0_dg0_coul_UEG(rho, g0, dg0drho)
  implicit none
  ! Give the on-top pair distribution function g0 and its derivative according to rho dg0drho
  !
  ! for a uniform electron gaz with full coulomb interaction
  !  
  ! the on-top pair density is then n^2(1-xi^2) g0

  double precision, intent (in) :: rho
  double precision, intent (out) :: g0, dg0drho
  double precision :: pi, g0_coul_UEG
  double precision :: g0_UEG_mu_inf, dg0drs
  double precision :: C1, F1, D1, E1, B1, rs, aa, bb

  pi = dacos(-1.d0)
  C1 = 0.0819306d0
  F1 = 0.752411d0
  D1 = -0.0127713d0
  E1 = 0.00185898d0
  B1 = 0.7317d0 - F1
  if(dabs(rho).gt.1.d-20)then
   rs = (3.d0 / (4.d0*pi*rho))**(1.d0/3.d0) 
  else
   rs = (3.d0 / (4.d0*pi*1.d-20))**(1.d0/3.d0) 
  endif
   
  g0 = g0_coul_UEG(rho)
  if(dabs(F1*rs).lt.50.d0)then
   aa = -B1 + 2.d0*C1*rs + 3.d0*D1*rs**2 + 4.d0*E1*rs**3
   bb = 1.d0 - B1*rs + C1*rs**2 + D1*rs**3 + E1*rs**4
   dg0drs = 0.5d0*(aa-F1*bb)*dexp(-F1*rs)
  else
   dg0drs = 0.d0
  endif
  if(dabs(rho).gt.1.d-20)then
   dg0drho = -((6.d0*dsqrt(pi)*rho**2)**(-2.d0/3.d0))*dg0drs
  else
   dg0drho = -((6.d0*dsqrt(pi)*1.d-40)**(-2.d0/3.d0))*dg0drs
  endif
 
  end 

double precision function g0_coul_UEG(rho)
! Pair distribution function g0(n) of the UEG with Coulomb interaction
!
! Taken from Eq. (46)  P. Gori-Giorgi and A. Savin, Phys. Rev. A 73, 032506 (2006).
 implicit none
 double precision, intent(in) :: rho
 double precision :: pi,x
 double precision :: B, C, D, E, d2, rs, ahd
 pi = 4d0 * datan(1d0)
 ahd = -0.36583d0
 d2 = 0.7524d0
 B = -2d0 * ahd - d2
 C = 0.08193d0
 D = -0.01277d0
 E = 0.001859d0
 x = -d2*rs
 if (dabs(rho) > 1.d-20) then
  rs = (3d0 / (4d0*pi*rho))**(1d0/3d0) 
  x = -d2*rs
  if(dabs(x).lt.50.d0)then
   g0_coul_UEG= 0.5d0 * (1d0+ rs* (-B + rs*(C + rs*(D + rs*E))))*dexp(x)
  else
   g0_coul_UEG= 0.d0
  endif
 else
  g0_coul_UEG= 0.d0
 endif
 g0_coul_UEG = max(g0_coul_UEG,1.d-14)

end
