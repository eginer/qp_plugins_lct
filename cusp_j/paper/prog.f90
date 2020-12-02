
double precision function slater_fit_ten_no(x)
 implicit none
 double precision, intent(in) :: x
 integer :: i
 double precision :: coef_fit_ten_no_slat_gauss(6)
 double precision :: expo_fit_ten_no_slat_gauss(6)
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

 slater_fit_ten_no = 0.d0
 do i = 1, 6
  slater_fit_ten_no = slater_fit_ten_no - coef_fit_ten_no_slat_gauss(i) * dexp(-expo_fit_ten_no_slat_gauss(i) * x * x)
 enddo                                                                                                                                       
end




double precision function shank3_f(array,n,nmax)
 implicit none
 integer, intent(in) :: n,nmax
 double precision, intent(in) :: array(0:nmax) ! array of the partial sums
 integer :: ntmp
 double precision :: shank1(0:nmax),shank2(0:nmax),shank3(0:nmax)
 ntmp = n
 call shank(array,ntmp,nmax,shank1)
 ntmp = ntmp - 2
 call shank(shank1,ntmp,nmax,shank2)
 ntmp = ntmp - 2
 call shank(shank2,ntmp,nmax,shank3)
 ntmp = ntmp - 2
 shank3_f = shank3(ntmp)
end


double precision function shank2_f(array,n,nmax)
 implicit none
 integer, intent(in) :: n,nmax
 double precision, intent(in) :: array(0:nmax) ! array of the partial sums
 integer :: ntmp
 double precision :: shank1(0:nmax),shank2(0:nmax)
 ntmp = n
 call shank(array,ntmp,nmax,shank1)
 ntmp = ntmp - 2
 call shank(shank1,ntmp,nmax,shank2)
 ntmp = ntmp - 2
 shank2_f = shank2(ntmp)
end



subroutine shank(array,n,nmax,shank1)
 implicit none
 integer, intent(in) :: n,nmax
 double precision, intent(in)  :: array(0:nmax)
 double precision, intent(out) :: shank1(0:nmax)
 integer :: i,j
 double precision :: shank_function
 do i = 1, n-1
  shank1(i-1) = shank_function(array,i,nmax)
 enddo
end

double precision function shank_function(array,i,n)
 implicit none
 integer, intent(in) :: i,n
 double precision, intent(in) :: array(0:n)
 double precision :: b_n, b_n1
 b_n = array(i) - array(i-1)
 b_n1 = array(i+1) - array(i)
 if(dabs(b_n1-b_n).lt.1.d-12)then
  shank_function = array(i+1) 
 else
  shank_function = array(i+1) - b_n1*b_n1/(b_n1-b_n)
 endif
end


double precision function f0(x,mu)
 implicit none
 double precision, intent(in) :: x,mu
 double precision :: pi
 pi = dacos(-1.d0)
 f0 = 0.5d0 * x * (1.d0 - derf(mu*x)) - 0.5d0/(dsqrt(pi)*mu)
end

double precision function f(x,mu)
 implicit none
 double precision, intent(in) :: x,mu
 double precision :: pi
 pi = dacos(-1.d0)
 f = 0.5d0 * x * (1.d0 - derf(mu*x)) - 0.5d0/(dsqrt(pi)*mu) * dexp(-(x*mu)**2.d0)
end

double precision function g_exp(x,mu)
 implicit none
 double precision, intent(in) :: x,mu
 double precision :: g
 g_exp = dexp(g(x,mu))
end

double precision function j(x,a,b)
 implicit none
 double precision, intent(in) :: x,a,b
 j = dexp(-a/b) * dexp(a*x/(1.d0 + b*x))
end

double precision function g(x,mu)
 implicit none
 double precision, intent(in) :: x,mu
 double precision :: pi
 pi = dacos(-1.d0)
 g = dsqrt(pi) * mu * x * (1.d0 - derf(mu*x)) - dexp(-(mu*x)**2.d0)
end

double precision function s(x,alpha,beta)
 implicit none
 double precision, intent(in) :: x,alpha,beta
 s = -dexp(-alpha*x -beta*x*x)
end

double precision function s_exp(x,alpha,beta)
 implicit none
 double precision, intent(in) :: x,alpha,beta
 double precision :: s
 s_exp = dexp(s(x,alpha,beta))
end

subroutine give_parameters(mu,alpha,beta)
 implicit none
 double precision, intent(in) :: mu
 double precision, intent(out):: alpha,beta
 double precision :: pi
 pi = dacos(-1.d0)
 alpha = dsqrt(pi) * mu
 beta = -mu*mu + 0.5d0 * alpha*alpha
end

double precision function dl_exp(x,n)
 implicit none
 double precision, intent(in) :: x
 integer, intent(in) :: n
 integer :: i
 double precision :: fact
 dl_exp = 1.d0
 fact = 1.d0
 do i = 1, n
  fact = fact * dble(i)
  dl_exp = dl_exp + 1.d0/fact * x**dble(i)
 enddo
end

subroutine dl_exp_rout(x,n,nmax,array)
 implicit none
 double precision, intent(in) :: x
 integer, intent(in) :: n,nmax
 double precision, intent(out) :: array(0:nmax)
 integer :: i
 double precision :: fact
 fact = 1.d0
 array(0) = 1.d0
 do i = 1, n
  fact = fact * dble(i)
  array(i) = array(i-1) + 1.d0/fact * x**dble(i)
 enddo
end


program pouet
 implicit none
 integer :: i,nx
 double precision :: x, dx,xmax,f,mu0,mu1,mu2,mu3,mu4,mu5,mu6,shank3_f
 double precision :: mu7, mu8, mu9,mu10,mu11, mu12, mu13, mu14
 double precision :: j,a,b
 double precision :: pi,f0,g,s,alpha,beta
 double precision :: alpha0,beta0
 double precision :: alpha1,beta1
 double precision :: alpha2,beta2
 double precision :: alpha3,beta3
 double precision :: alpha4,beta4
 double precision :: alpha5,beta5
 double precision :: alpha6,beta6
 double precision :: alpha7,beta7
 double precision :: alpha8,beta8
 double precision :: s_exp,g_exp,dl_exp,shank_exp,shank2_f
 double precision :: array(0:100), shank1(0:100),shank2(0:100),array_tmp(0:100)
 double precision :: contrib(0:100),shank3(0:100)
 double precision :: sum,mu,slater_fit_ten_no
 integer :: n,nmax,k
 pi = dacos(-1.d0)
 xmax = 2.d0
 nx = 10000
 dx = xmax/dble(nx)
 print*,'0.5**0.25', 0.5d0**0.25d0
 mu0 = 0.05d0
 mu1 = 0.10d0
 mu2 = 0.20d0
 mu3 = 0.3d0
 mu4 = 0.4d0
 mu5 = 0.5d0
 mu6 = 0.6d0
 mu7 = 0.7071d0
 mu8 = 0.8d0
 mu9 = 0.9d0
 mu10= 1.0d0
 mu11= 1.5d0
 mu12= 3.0d0
 mu13= 5.0d0
! mu15= 10.0d0

! mu0 = 1.00d0
! mu1 = 1.50d0
! mu2 = 2.00d0
! mu3 = 2.5d0
! mu4 = 3.0d0
! mu5 = 3.5d0
! mu6 = 4.0d0
! mu7 = 4.5000d0
! mu8 = 5.0d0
! mu9 = 5.5d0
! mu10= 6.0d0

 write(33,'(A5,12X,100(F16.10,X))')"#  ", mu0, mu1, mu2, mu3, mu4, mu5, mu6,mu7, mu8, mu9, mu10
 x = dx
 do i = 1, nx
  write(33,'(100(F16.10,X))')x,dexp(f(x,mu0)),dexp(f(x,mu1)),dexp(f(x,mu2))& 
                              ,dexp(f(x,mu3)),dexp(f(x,mu4)),dexp(f(x,mu5)),dexp(f(x,mu6))&
                              ,dexp(f(x,mu7)),dexp(f(x,mu8)),dexp(f(x,mu9)),dexp(f(x,mu10))& 
                              ,dexp(f(x,mu11)),dexp(f(x,mu12)),dexp(f(x,mu13)),dexp(slater_fit_ten_no(x))
                              
  write(34,'(100(F16.10,X))')x,(f(x,mu0)),(f(x,mu1)),(f(x,mu2))& 
                              ,(f(x,mu3)),(f(x,mu4)),(f(x,mu5)),(f(x,mu6))&
                              ,(f(x,mu7)),(f(x,mu8)),(f(x,mu9)),(f(x,mu10))&
                              ,(f(x,mu11)),(f(x,mu12)),(f(x,mu13)),slater_fit_ten_no(x)
  x = x + dx
 enddo

 stop

 xmax = 5.d0
 dx = xmax/dble(nx)

 mu0 = 1.0d0
 call give_parameters(mu0,alpha0,beta0)
 mu1 = 1.5d0
 call give_parameters(mu1,alpha1,beta1)
 mu2 = 2.0d0
 call give_parameters(mu2,alpha2,beta2)
 mu3 = 2.5d0
 call give_parameters(mu3,alpha3,beta3)
 mu4 = 3.0d0
 call give_parameters(mu4,alpha4,beta4)
 mu5 = 3.5d0
 call give_parameters(mu5,alpha5,beta5)
 mu6 = 4.0d0
 call give_parameters(mu6,alpha6,beta6)
 mu7 = 4.5d0
 call give_parameters(mu7,alpha7,beta7)
 mu8 = 5.0d0
 call give_parameters(mu7,alpha8,beta8)
 x = dx
! do i = 1, nx
!  write(34,'(100(F16.10,X))')x,g(x,mu0),s(x,alpha0,beta0)& 
!                              ,g(x,mu1),s(x,alpha1,beta1)& 
!                              ,g(x,mu2),s(x,alpha2,beta2)& 
!                              ,g(x,mu3),s(x,alpha3,beta3)& 
!                              ,g(x,mu4),s(x,alpha4,beta4)& 
!                              ,g(x,mu5),s(x,alpha5,beta5)& 
!                              ,g(x,mu6),s(x,alpha6,beta6) 
!  x = x + dx
! enddo

 mu0 = 0.1d0
 call give_parameters(mu0,alpha0,beta0)
 mu1 = 0.2d0
 call give_parameters(mu1,alpha1,beta1)
 mu2 = 0.3d0
 call give_parameters(mu2,alpha2,beta2)
 mu3 = 0.4d0
 call give_parameters(mu3,alpha3,beta3)
 mu4 = 0.5d0
 call give_parameters(mu4,alpha4,beta4)
 mu5 = 0.6d0
 call give_parameters(mu5,alpha5,beta5)
 mu6 = 0.7d0
 call give_parameters(mu6,alpha6,beta6)
 mu7 = 0.8d0
 call give_parameters(mu7,alpha7,beta7)
 mu8 = 0.9d0
 call give_parameters(mu7,alpha8,beta8)
 x = dx

 mu0 = 1.0d0
 call give_parameters(mu0,alpha0,beta0)
 mu1 = 1.5d0
 call give_parameters(mu1,alpha1,beta1)
 mu2 = 2.0d0
 call give_parameters(mu2,alpha2,beta2)
 mu3 = 2.5d0
 call give_parameters(mu3,alpha3,beta3)
 mu4 = 3.0d0
 call give_parameters(mu4,alpha4,beta4)
 mu5 = 3.5d0
 call give_parameters(mu5,alpha5,beta5)
 mu6 = 4.0d0
 call give_parameters(mu6,alpha6,beta6)
 mu7 = 4.5d0
 call give_parameters(mu7,alpha7,beta7)
 mu8 = 5.0d0
 call give_parameters(mu7,alpha8,beta8)
 x = dx
 n = 6
 nmax = 100
 x = dx
 mu = 0.01d0
 do i = 1, nx
  call dl_exp_rout(g(x,mu),n,nmax,array)
  shank_exp = shank2_f(array,n,nmax)
  write(35,'(100(F16.10,X))')x,g_exp(x,mu),dl_exp(g(x,mu),n),shank_exp
  !,dl_exp(s(x,alpha0,beta0),n)
!                              ,dl_exp(g(x,mu1),n),dl_exp(s(x,alpha1,beta1),n)& 
!                              ,dl_exp(g(x,mu2),n),dl_exp(s(x,alpha2,beta2),n)& 
!                              ,dl_exp(g(x,mu3),n),dl_exp(s(x,alpha3,beta3),n)& 
!                              ,dl_exp(g(x,mu4),n),dl_exp(s(x,alpha4,beta4),n)& 
!                              ,dl_exp(g(x,mu5),n),dl_exp(s(x,alpha5,beta5),n)& 
!                              ,dl_exp(g(x,mu6),n),dl_exp(s(x,alpha6,beta6),n) 
  x = x + dx
 enddo
 stop 
 n = 10
 x = dx
 do i = 1, nx
  call dl_exp_rout(-x,n,nmax,array)
  shank_exp = shank3_f(array,n,nmax)
  write(36,'(100(F16.10,X))')x,dexp(-x),dl_exp(-x,n), shank_exp
  x = x + dx
 enddo
 
 sum = 0.d0
 n = 12
 do i = 0, n
  contrib(i) = 4.d0 * (-1.d0)**dble(i)/(2.d0*dble(i)+1.d0)
  sum = sum + contrib(i)
 enddo
 array = 0.d0
 do i = 0, n
  do k = 0, i
   array(i) = array(i) + contrib(k)
  enddo
 enddo
 print*,'array(n) = ',array(n)
 print*,'' 

 nmax = 100
 print*,'sum = ',sum
 print*,'First Shank'
 call shank(array,n,nmax,shank1)
 do i = 0, n-2
  print*,i,'shank1 = ',shank1(i)
 enddo
 print*,'Second Shank'
 n = n - 2
 call shank(shank1,n,nmax,shank2)
 do i = 0, n-2
  print*,i,'shank2 = ',shank2(i)
 enddo
 print*,'Third Shank'
 n = n - 2
 call shank(shank2,n,nmax,shank3)
 do i = 0, n-2
  print*,i,'shank3 = ',shank3(i)
 enddo
 n = 12
 sum = shank3_f(array,n,nmax)
 print*,'sum = ',sum
 print*,'pi  = ',pi
end

