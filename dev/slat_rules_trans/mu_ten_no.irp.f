BEGIN_PROVIDER [double precision, mu_ten_no]
 implicit none
 double precision :: slater_fit_ten_no,c0
 c0 = slater_fit_ten_no(0.d0)
 mu_ten_no = -1.d0/(2.d0 * dsqrt(dacos(-1.d0))*c0 )

END_PROVIDER 
