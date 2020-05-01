
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

 BEGIN_PROVIDER [integer, n_max_fit_ten_no_slat]
 implicit none
 n_max_fit_ten_no_slat = 6
 END_PROVIDER

 BEGIN_PROVIDER [double precision, coef_fit_ten_no_slat_gauss, (n_max_fit_ten_no_slat)]
&BEGIN_PROVIDER [double precision, expo_fit_ten_no_slat_gauss, (n_max_fit_ten_no_slat)]
 implicit none
 BEGIN_DOC
 ! fit the exp(-x) as 
 !
 ! \sum_{i = 1, 8} coef_fit_ten_no_slat_gauss(i) * exp(-expo_fit_ten_no_slat_gauss(i) * x**2)
 !
 ! The coefficient are taken from the ano-rcc expansion of the hydrogen 1s orbital
 END_DOC
 coef_fit_ten_no_slat_gauss(1) = 0.078215d0
 coef_fit_ten_no_slat_gauss(2) = 0.132037d0
 coef_fit_ten_no_slat_gauss(3) = 0.068633d0
 coef_fit_ten_no_slat_gauss(4) = 0.029047d0
 coef_fit_ten_no_slat_gauss(5) = 0.012063d0
 coef_fit_ten_no_slat_gauss(6) = 0.004346d0

 expo_fit_ten_no_slat_gauss(1) = 0.621698d0
 expo_fit_ten_no_slat_gauss(2) = 3.371717d0
 expo_fit_ten_no_slat_gauss(3) = 14.27116d0
 expo_fit_ten_no_slat_gauss(4) = 82.76522d0
 expo_fit_ten_no_slat_gauss(5) = 605.5295d0
 expo_fit_ten_no_slat_gauss(6) = 6596.808d0

END_PROVIDER 


double precision function slater_fit_ten_no(x)
 implicit none
 double precision, intent(in) :: x
 integer :: i
 slater_fit_ten_no = 0.d0
 do i = 1, n_max_fit_ten_no_slat
  slater_fit_ten_no += coef_fit_ten_no_slat_gauss(i) * dexp(-expo_fit_ten_no_slat_gauss(i) * x * x)
 enddo
end

double precision function slater_ten_no(x,gam)
 implicit none
 double precision, intent(in) :: x,gam
 integer :: i
 slater_ten_no = dexp(-0.5d0/gam * dexp(- gam * x))
end
