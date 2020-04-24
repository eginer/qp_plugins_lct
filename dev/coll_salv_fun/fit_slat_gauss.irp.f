 BEGIN_PROVIDER [double precision, coef_fit_slat_gauss, (8)]
&BEGIN_PROVIDER [double precision, expo_fit_slat_gauss, (8)]
 implicit none
 BEGIN_DOC
 ! fit the exp(-x) as 
 !
 ! \sum_{i = 1, 8} coef_fit_slat_gauss(i) * exp(-expo_fit_slat_gauss(i) * x**2)
 !
 ! The coefficient are taken from the ano-rcc expansion of the hydrogen 1s orbital
 END_DOC
 coef_fit_slat_gauss(1)  = 0.00096385
 coef_fit_slat_gauss(2)  = 0.00749196
 coef_fit_slat_gauss(3)  = 0.03759541
 coef_fit_slat_gauss(4)  = 0.14339498
 coef_fit_slat_gauss(5)  = 0.34863630
 coef_fit_slat_gauss(6)  = 0.43829736
 coef_fit_slat_gauss(7)  = 0.16510661
 coef_fit_slat_gauss(8)  = 0.02102287               

 expo_fit_slat_gauss(1)  = 188.61445    
 expo_fit_slat_gauss(2)  =  28.276596   
 expo_fit_slat_gauss(3)  =   6.4248300  
 expo_fit_slat_gauss(4)  =   1.8150410  
 expo_fit_slat_gauss(5)  =    .59106300 
 expo_fit_slat_gauss(6)  =    .21214900 
 expo_fit_slat_gauss(7)  =    .07989100 
 expo_fit_slat_gauss(8)  =    .02796200 

END_PROVIDER 

 BEGIN_PROVIDER [double precision, norm_fit_slat]
&BEGIN_PROVIDER [double precision, inv_sq_norm_fit_slat]
 implicit none
 integer :: n,l, dim1
 double precision :: A_center(3), B_center(3), alpha,beta
 double precision :: overlap_x,overlap_y,overlap_z,overlap, c 
 integer :: power_A(3), power_B(3)
 dim1=100
 A_center = 0.d0
 power_A  = 0
 B_center = 0.d0
 power_B  = 0
 norm_fit_slat = 0.d0
 do n = 1, 8
  alpha = expo_fit_slat_gauss(n)
  do l = 1, 8
   beta = expo_fit_slat_gauss(l)
   call overlap_gaussian_xyz(A_center,B_center,alpha,beta,power_A,power_B,overlap_x,overlap_y,overlap_z,overlap,dim1)
   c = coef_fit_slat_gauss(n) * coef_fit_slat_gauss(l)
   norm_fit_slat += c * overlap
  enddo
 enddo
 inv_sq_norm_fit_slat = 1.d0/dsqrt(norm_fit_slat)

END_PROVIDER 

double precision function slater_fit(x)
 implicit none
 double precision, intent(in) :: x
 integer :: i
 slater_fit = 0.d0
 do i = 1, 8
  slater_fit += coef_fit_slat_gauss(i) * dexp(-expo_fit_slat_gauss(i) * x * x)
 enddo
 slater_fit *= inv_sq_norm_fit_slat
end
