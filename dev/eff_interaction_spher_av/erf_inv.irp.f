 subroutine give_third_order_polynome_solutions(a,b,c,d,x1,x2,x3)
 implicit none
 !Formule found here:
 !https://www.lucaswillems.com/fr/articles/58/equations-troisieme-degre?cache=update
 double precision, intent(in)  :: a 
 double precision, intent(in)  :: b
 double precision, intent(in)  :: c
 double precision, intent(in)  :: d
 double precision, intent(out) :: x1 
 double precision, intent(out) :: x2
 double precision, intent(out) :: x3

 double precision :: t,p,q,delta_1,delta_2
 !x1 solution:
 p = (3.d0 * a * c -b**2.d0)/(3.d0*a**2.d0)
 q = (2.d0 * b**3.d0 - 9.d0 * a * b * c + 27.d0*a**2.d0*d )/(27.d0*a**3.d0)
 
 delta_1 = q**2.d0 +4.d0*p**3/27.d0
 t = ((-q-sqrt(delta_1))/2.d0)**(1.d0/3.d0)+((-q+sqrt(delta_1))/2.d0)**(1.d0/3.d0)  

 x1 = t - b/(3.d0*a)

 !x1,x3 solutions:
 delta_2 = (b+a*x1)**2.d0 - 4.d0*a*(c + (b+a*x1) * x1)

 x2 = (-b -a*x1 -sqrt(delta_2))/(2.d0*a)
 x3 = (-b -a*x1 +sqrt(delta_2))/(2.d0*a)

 end


 subroutine give_mu_r12_looser(eff_inter,r12,mu_1,mu_2,mu_3)
 implicit none
 double precision, intent(in)  :: eff_inter
 double precision, intent(in)  :: r12 
 double precision, intent(out) :: mu_1 
 double precision, intent(out) :: mu_2
 double precision, intent(out) :: mu_3
 include 'utils/constants.include.F' 
 double precision :: a,b,c,d
 a = -2.d0 * r12 **2.d0/(3.d0*sqrt(pi))
 b = 0.d0
 c = 2.d0/sqrt(pi)
 d = -eff_inter

 call give_third_order_polynome_solutions(a,b,c,d,mu_1,mu_2,mu_3)

 end


 subroutine give_mu_r12(eff_inter,r12,mu)
 implicit none
 double precision, intent(in)  :: eff_inter
 double precision, intent(in)  :: r12 
 double precision, intent(out) :: mu 
 double precision :: y,erf_inv,thr
 y = eff_inter*r12 
 thr = y*1.d-4
!print*,'eff_inter =',eff_inter
!print*,'eff_inter*r12 =',y
 mu= erf_inv(y,thr)/r12 
 
 end

 double precision function erf_inv(y,thr)
  implicit none
  BEGIN_DOC
 ! inverse function of erf
 ! erf( erf_inv(y) ) = y
 ! thr = threshold on y
  END_DOC
  double precision, intent(in) :: y,thr
  double precision :: x_plus,x_plus_before,x_moins
  x_plus = 50.d0
  x_moins = 0.d0
  x_plus_before = x_plus
  if(y.ge.1.d0)then
   print*,'y must lt than 1'
   print*,y
   stop
  endif
  do while(dabs(derf(0.5d0 *( x_plus+x_moins)) - y).gt.thr)
   if(y.le.derf(x_plus))then
    x_plus_before = x_plus
    x_plus = ( x_plus + x_moins ) * 0.5d0
   else
    x_moins = x_plus
    x_plus = x_plus_before
   endif
!  print*,'x_moins,x_plus',x_moins,x_plus 
  enddo
  erf_inv = 0.5d0 * (x_moins + x_plus)
 end
