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
 coef_fit_slat_gauss(1) = 7.1767030247679102d-004
 coef_fit_slat_gauss(2) = 5.5540200969662898d-003
 coef_fit_slat_gauss(3) = 2.8426009700541901d-002
 coef_fit_slat_gauss(4) = 0.10951630526849800     
 coef_fit_slat_gauss(5) = 0.30049708190292101     
 coef_fit_slat_gauss(6) = 0.47410183437598802     
 coef_fit_slat_gauss(7) = 0.22373854296317000     
 coef_fit_slat_gauss(8) = 2.1691864646938498d-003

 expo_fit_slat_gauss(1)  = 188.61445d0    
 expo_fit_slat_gauss(2)  =  28.276596d0  
 expo_fit_slat_gauss(3)  =   6.4248300d0  
 expo_fit_slat_gauss(4)  =   1.8150410d0 
 expo_fit_slat_gauss(5)  =    .59106300d0 
 expo_fit_slat_gauss(6)  =    .21214900d0 
 expo_fit_slat_gauss(7)  =    .07989100d0 
 expo_fit_slat_gauss(8)  =    .02796200d0 

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
