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


double precision function cos_x14_int(m)
 implicit none
 BEGIN_DOC
! int(-m*pi;+m*pi) x^12 cos(x) 
 END_DOC
 include 'constants.include.F'
 integer, intent(in) :: m
 double precision :: dm,dm2,pi2,pi4,dm4,pi6,dm6,pi8,dm8
 dm = dble(m)
 dm2=dm*dm
 pi2=pi*pi
 pi4=pi2*pi2
 dm4=dm2*dm2
 pi6=pi4*pi2
 dm6=dm4*dm2
 pi8=pi4*pi4
 dm8=dm4*dm4
 cos_x14_int = 28.d0 * pi * dm * (pi6*pi6 * dm6*dm6 - 156.d0 *pi6*pi4* dm6*dm4 + 17160.d0 * pi8 * dm8 - 1235520.d0*pi6*dm6 + 51891840.d0*pi4*dm4 - 1037836800.d0 * pi2*dm2 + 6227020800.d0)*dcos(pi*dm)

end


double precision function cos_x16_int(m)
 implicit none
 BEGIN_DOC
! int(-m*pi;+m*pi) x^12 cos(x) 
 END_DOC
 include 'constants.include.F'
 integer, intent(in) :: m
 double precision :: dm,dm2,pi2,pi4,dm4,pi6,dm6,pi8,dm8,pi10,dm10
 dm = dble(m)
 dm2=dm*dm
 pi2=pi*pi
 pi4=pi2*pi2
 dm4=dm2*dm2
 pi6=pi4*pi2
 dm6=dm4*dm2
 pi8=pi4*pi4
 dm8=dm4*dm4
 pi10=pi4*pi6
 dm10=dm4*dm6
 cos_x16_int = 32.d0 * pi * dm * (pi8*pi6 * dm8*dm6 - 210.d0 *pi6*pi6* dm6*dm6 + 32760.d0 * pi10 * dm10 - 3603600.d0*pi8*dm8 + 259459200.d0*pi6*dm6 - 10897286400.d0 * pi4*dm4 + 217945728000.d0 * pi2*dm2 + 1307674368000.d0)*dcos(pi*dm)

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
 print*,'nx =',nx
 dx = 2.d0*domain/dble(nx+1)
 print*,'dx = ',dx
 x_min = -domain
 x = x_min
 num_int_cos_pol = 0.d0
 do i = 1, nx
  num_int_cos_pol += dcos(x) * x**dn
  x += dx
 enddo
 num_int_cos_pol *= dx
end

double precision function int_xn_cos_m(m,n)
 implicit none
 BEGIN_DOC
! int(-m*pi;+m*pi) x^n cos(x) 
 END_DOC
 integer, intent(in) :: m,n
 double precision :: cos_x2_int,cos_x4_int,cos_x6_int,cos_x8_int,cos_x10_int,cos_x12_int,cos_x14_int
 if(n.gt.14)then
  print*,'int_xn_cos_m is defined only for n up to 14'
  print*,'you called it for n = ',n
  stop
 endif
 int_xn_cos_m = 0.d0 
 if (iand(n,1).eq.0)then 
  if(n==0)then
   int_xn_cos_m = 0.d0
  else if(n==2)then
   int_xn_cos_m = cos_x2_int(m)
  else if(n==4)then
   int_xn_cos_m = cos_x4_int(m)
  else if(n==6)then
   int_xn_cos_m = cos_x6_int(m)
  else if(n==8)then
   int_xn_cos_m = cos_x8_int(m)
  else if(n==10)then
   int_xn_cos_m = cos_x10_int(m)
  else if(n==12)then
   int_xn_cos_m = cos_x12_int(m)
  else if(n==14)then
   int_xn_cos_m = cos_x14_int(m)
  endif
 endif

end
