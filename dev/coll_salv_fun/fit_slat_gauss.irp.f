 BEGIN_PROVIDER [integer, n_max_fit_slat]
 implicit none
 BEGIN_DOC
! number of gaussian to fit exp(-x)
!
! I took the 8 gaussians of the ANO-RCC basis set of the H atom
 END_DOC
 n_max_fit_slat = 8
 END_PROVIDER

 BEGIN_PROVIDER [double precision, coef_fit_slat_gauss, (n_max_fit_slat)]
&BEGIN_PROVIDER [double precision, expo_fit_slat_gauss, (n_max_fit_slat)]
 implicit none
  include 'constants.include.F'
 BEGIN_DOC
 ! fit the exp(-x) as 
 !
 ! \sum_{i = 1, 8} coef_fit_slat_gauss(i) * exp(-expo_fit_slat_gauss(i) * x**2)
 !
 ! The coefficient are taken from the ano-rcc expansion of the hydrogen 1s orbital
 END_DOC

 coef_fit_slat_gauss(1) =  2.6032524673598918d-002
 coef_fit_slat_gauss(2) =  4.8538671362083070d-002
 coef_fit_slat_gauss(3) =  8.1756487143886281d-002
 coef_fit_slat_gauss(4) = 0.1220544746390489d0    
 coef_fit_slat_gauss(5) = 0.1443695761749570d0    
 coef_fit_slat_gauss(6) = 0.1056239657977568d0    
 coef_fit_slat_gauss(7) =  2.3962067692955683d-002
 coef_fit_slat_gauss(8) =  1.0571415716647108d-004

 coef_fit_slat_gauss *= sqpi ! for normalization

 expo_fit_slat_gauss(1)  = 188.61445d0    
 expo_fit_slat_gauss(2)  =  28.276596d0  
 expo_fit_slat_gauss(3)  =   6.4248300d0  
 expo_fit_slat_gauss(4)  =   1.8150410d0 
 expo_fit_slat_gauss(5)  =    .59106300d0 
 expo_fit_slat_gauss(6)  =    .21214900d0 
 expo_fit_slat_gauss(7)  =    .07989100d0 
 expo_fit_slat_gauss(8)  =    .02796200d0 

END_PROVIDER 

double precision function slater_fit_gam(x,gam)
 implicit none
 double precision, intent(in) :: x,gam
 BEGIN_DOC
! fit of the function exp(-gam * x) with gaussian functions 
 END_DOC
 integer :: i
 slater_fit_gam = 0.d0
 do i = 1, n_max_fit_slat
  slater_fit_gam += coef_fit_slat_gauss(i) * dexp(-expo_fit_slat_gauss(i) * gam * gam * x * x)
 enddo
end

subroutine expo_fit_slater_gam(gam,expos)
 implicit none
 BEGIN_DOC
! returns the array of the exponents of the gaussians to fit exp(-gam*x)
 END_DOC
 double precision, intent(in)  :: gam
 double precision, intent(out) :: expos(n_max_fit_slat)
 integer :: i
 do i = 1, n_max_fit_slat
  expos(i) = expo_fit_slat_gauss(i) * gam * gam
 enddo
end

