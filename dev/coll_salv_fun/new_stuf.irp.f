double precision function u0(t,a,b)
 implicit none
 include 'constants.include.F'
 double precision, intent(in) :: t,a,b
 double precision :: cst
 cst = 0.d0

u0 = 1.d0/(16.d0 * b**(1.5d0)) * dexp(-t * (a + b * t)) *  (2.d0 * dsqrt(b) * (-1 + 4.d0*b * dexp(t * (a +b*t)) & 
 * (t + 2.d0*cst)) + dexp((a + 2.d0*b*t)**2/(4.d0*b)) * dsqrt(pi) * (2.d0*b*t*derf(a/(2.d0 * dsqrt(b))) &
 - (a + 2.d0*b*t) * derf((a + 2.d0*b*t)/(2.d0* dsqrt(b)))))

end

