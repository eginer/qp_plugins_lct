double precision function cos_x2_int(m)
 implicit none
 BEGIN_DOC
! int(-m*pi;+m*pi) x^2 cos(x) 
 END_DOC
 include 'constants.include.F'
 integer, intent(in) :: m
 double precision :: dm
 dm = dble(m)
! cos_x2_int = 2.d0 * (pi*pi*dm*dm - 2.d0) * dsin(pi*m) + 4.d0 * pi * dm * dcos(pi*m)
 cos_x2_int = 4.d0 * pi * dm * dcos(pi*m)

end

double precision function cos_x4_int(m)
 implicit none
 BEGIN_DOC
! int(-m*pi;+m*pi) x^4 cos(x) 
 END_DOC
 include 'constants.include.F'
 integer, intent(in) :: m
 double precision :: dm
 dm = dble(m)
 cos_x4_int = 8.d0 * pi * dm * (pi*pi*dm*dm - 6.d0) * dcos(pi*m) 
! cos_x2_int = 8.d0 * pi * dm * (pi*pi*dm*dm - 6.d0) * dcos(pi*m) + 2.d0 * (pi * dm * )* dsin(pi*m)

end

double precision function cos_x6_int(m)
 implicit none
 BEGIN_DOC
! int(-m*pi;+m*pi) x^6 cos(x) 
 END_DOC
 include 'constants.include.F'
 integer, intent(in) :: m
 double precision :: dm,dm2,pi2
 dm = dble(m)
 dm2=dm*dm
 pi2=pi*pi
 cos_x6_int = 12.d0 * pi * dm * (pi2*pi2 * dm2*dm2 - 20.d0 *pi2 * dm2 + 120.d0)*dcos(pi*dm)

end

double precision function cos_x8_int(m)
 implicit none
 BEGIN_DOC
! int(-m*pi;+m*pi) x^8 cos(x) 
 END_DOC
 include 'constants.include.F'
 integer, intent(in) :: m
 double precision :: dm,dm2,pi2,pi4,dm4
 dm = dble(m)
 dm2=dm*dm
 pi2=pi*pi
 pi4=pi2*pi2
 dm4=dm2*dm2
 cos_x8_int = 16.d0 * pi * dm * (pi4*pi2 * dm4*dm2 - 42.d0 *pi4 * dm4 + 840.d0 * pi2 * dm2 - 5040.d0)*dcos(pi*dm)

end


double precision function cos_x10_int(m)
 implicit none
 BEGIN_DOC
! int(-m*pi;+m*pi) x^10 cos(x) 
 END_DOC
 include 'constants.include.F'
 integer, intent(in) :: m
 double precision :: dm,dm2,pi2,pi4,dm4
 dm = dble(m)
 dm2=dm*dm
 pi2=pi*pi
 pi4=pi2*pi2
 dm4=dm2*dm2
 cos_x10_int = 20.d0 * pi * dm * (pi4*pi4 * dm4*dm4 - 72.d0 *pi4*pi2* dm4*dm2 + 3024.d0 * pi4 * dm4 - 60480.d0*pi2*dm2 + 362880.d0)*dcos(pi*dm)

end


double precision function cos_x12_int(m)
 implicit none
 BEGIN_DOC
! int(-m*pi;+m*pi) x^12 cos(x) 
 END_DOC
 include 'constants.include.F'
 integer, intent(in) :: m
 double precision :: dm,dm2,pi2,pi4,dm4,pi6,dm6
 dm = dble(m)
 dm2=dm*dm
 pi2=pi*pi
 pi4=pi2*pi2
 dm4=dm2*dm2
 pi6=pi4*pi2
 dm6=dm4*dm2
 cos_x12_int = 24.d0 * pi * dm * (pi6*pi4 * dm6*dm4 - 110.d0 *pi4*pi4* dm4*dm4 + 7920.d0 * pi6 * dm6 - 332640.d0*pi4*dm4 + 6652800.d0*pi2*dm2 -39916800.d0)*dcos(pi*dm)

end


double precision function num_int_cos_pol(m,n,nx)
 implicit none
 integer, intent(in) :: m,nx,n
 BEGIN_DOC
! int(-m*pi;+m*pi) x^n cos(x) 
 END_DOC
 double precision :: domain,dx,dm,x_min,x,dn
 integer :: i
 include 'constants.include.F'
 dm = dble(m)
 dn = dble(n)
 domain = pi*dm
 dx = 2.d0*domain/dble(nx+1)
 x_min = -domain
 x = x_min
 num_int_cos_pol = 0.d0
 do i = 1, nx
  num_int_cos_pol += dcos(x) * x**dn
  x += dx
 enddo
 num_int_cos_pol *= dx
end
