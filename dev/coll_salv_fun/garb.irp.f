
double precision function j_factor_slat(x,g0,g1,g2)
 implicit none
 double precision, intent(in) :: x,g0,g1,g2
 BEGIN_DOC
! simple jastrow factor 
!
! g0 is the depth of the hole at x = 0
!
! g1 is the first derivative in x WHICH MUST BE POSITIVE
!
! g2 is the second derivative 
 END_DOC
 double precision :: slater_fit_gam,gam,alpha,beta
 gam = 1.d0/(1.d0 - g0)
 alpha = gam * g1
 beta = (alpha**2 - 2.d0 * g2 * gam) / (2.d0 * gam)
 j_factor_slat = 1.d0 - 1.d0/gam * slater_fit_gam(x,alpha)* dexp(-beta*x*x )
end

