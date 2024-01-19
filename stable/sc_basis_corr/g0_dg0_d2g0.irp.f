!  subroutine g0_dg0_d2g0(rho, rho_a, rho_b, g0, dg0drho, d2g0drho2)
!
!  implicit none
!  BEGIN_DOC
!  ! Give the on-top pair distribution function g0 second derivative according to rho d2g0drho2
!  END_DOC
!
!  double precision, intent (in) :: rho, rho_a, rho_b
!  double precision, intent (out) :: g0, dg0drho, d2g0drho2
!  double precision :: pi
!  double precision :: g0_UEG_mu_inf, dg0drs, d2g0drs2, d2rsdrho2
!  double precision :: C1, F1, D1, E1, B1, rs
!
!  pi = dacos(-1.d0)
!  C1 = 0.0819306d0
!  F1 = 0.752411d0
!  D1 = -0.0127713d0
!  E1 = 0.00185898d0
!  B1 = 0.7317d0 - F1
!  if(dabs(rho).gt.1.d-20)then
!   rs = (3.d0 / (4.d0*pi*rho))**(1.d0/3.d0)
!  else
!   rs = (3.d0 / (4.d0*pi*1.d-20))**(1.d0/3.d0)
!  endif
!
!  g0 = g0_UEG_mu_inf(rho_a, rho_b)
!  if(dabs(F1*rs).lt.50.d0)then
!   dg0drs = 0.5d0*((-B1 + 2.d0*C1*rs + 3.d0*D1*rs**2 + 4.d0*E1*rs**3)-F1*(1.d0 - B1*rs + C1*rs**2 + D1*rs**3 + E1*rs**4))*dexp(-F1*rs)
!   d2g0drs2 = 0.5d0*((2.d0*C1 + 6.d0*D1*rs + 12*E1*rs**2) - 2.d0*F1*(-B1 + 2.d0*C1*rs + 3.d0*D1*rs**2 + 4.d0*E1*rs**3)&
!              &+ (F1**2)*(1.d0 - B1*rs + C1*rs**2 + D1*rs**3 + E1*rs**4))*dexp(-F1*rs)
!  else
!   dg0drs = 0.d0
!   d2g0drs2 = 0.d0
!  endif
!
!  if(dabs(rho).gt.1.d-20)then
!   dg0drho = -((6.d0*dsqrt(pi)*rho**2)**(-2.d0/3.d0))*dg0drs
!   d2rsdrho2 = -8.d0*dsqrt(pi)*rho*(6.d0*dsqrt(pi)*rho**2)**(-5.d0/3.d0)
!   d2g0drho2 = dg0drho*d2rsdrho2 -((6.d0*dsqrt(pi)*rho**2)**(-4.d0/3.d0))*d2g0drs2
!  else
!   dg0drho = -((6.d0*dsqrt(pi)*1.d-40)**(-2.d0/3.d0))*dg0drs
!   d2rsdrho2 = -8.d0*dsqrt(pi)*(1.d-20)*(6.d0*dsqrt(pi)*1.d-40)**(-5.d0/3.d0)
!   d2g0drho2 = dg0drho*d2rsdrho2 - ((6.d0*dsqrt(pi)*1.d-40)**(-4.d0/3.d0))*d2g0drs2
!  endif
!
!  end subroutine 

