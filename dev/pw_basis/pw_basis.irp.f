program pw_basis
  implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC
!  provide pw_eigval
!  call routine_2
!   provide n_delta_k_uniq
 double precision :: integral_exact,integral_num
 double precision :: cos_x2_int,cos_x4_int,cos_x6_int,cos_x8_int,cos_x10_int,cos_x12_int,num_int_cos_pol
 integer :: n,i,nx,m
 include 'constants.include.F'
 print*,''
 print*,''
 m = 1
 n = 12
 print*,'m,n',m,n
 do i = 5,7
  nx = 10**i
  integral_num   = num_int_cos_pol(m,n,nx)
  integral_exact = cos_x12_int(m)
  print*,''
  print*,nx
  print*,integral_exact,integral_num,dabs(integral_exact-integral_num)/dabs(integral_exact)
 enddo

 print*,''
 print*,''
 m = 2
 n = 12
 print*,'m,n',m,n
 do i = 5,7
  nx = 10**i
  integral_num   = num_int_cos_pol(m,n,nx)
  integral_exact = cos_x12_int(m)
  print*,''
  print*,nx
  print*,integral_exact,integral_num,dabs(integral_exact-integral_num)/dabs(integral_exact)
 enddo

 print*,''
 print*,''
 m = 3
 n = 12
 print*,'m,n',m,n
 do i = 5,7
  nx = 10**i
  integral_num   = num_int_cos_pol(m,n,nx)
  integral_exact = cos_x12_int(m)
  print*,''
  print*,nx
  print*,integral_exact,integral_num,dabs(integral_exact-integral_num)/dabs(integral_exact)
 enddo

 print*,''
 print*,''
 m = 4
 n = 12
 print*,'m,n',m,n
 do i = 5,7
  nx = 10**i
  integral_num   = num_int_cos_pol(m,n,nx)
  integral_exact = cos_x12_int(m)
  print*,''
  print*,nx
  print*,integral_exact,integral_num,dabs(integral_exact-integral_num)/dabs(integral_exact)
 enddo
! complex*8 :: z
! z = dcmplx(1.d0,1.d0)
! print*,erf(1.d0)
! print*,''
! print*,'erf(z) = ',erf(z)
! n = 10
! alpha = 2.d0
! do i = 1, 7
!  dx = 10.d0**(-dble(i))
!  call int_gauss_cos(n,alpha,nx,dx,integral)
!  print*,nx,dx,integral
! enddo
end
subroutine routine
  implicit none
  double precision :: r(3),dx,ao_value,aos,val_1,val_2,val_3,mos_array(mo_num)
  complex*8 :: pw_x,pw_y,pw_z,pw
  integer :: n,i,i_pw,i_vec
  n = 100
  dx = pw_box_x/dble(n)
  r = 0.d0 
  r(1) = -0.5d0 * pw_box_x
  print*,'delta_k_pw(1)',delta_k_pw(1)
  do i = 1, n+1
   call give_pw_i_at_r(i_pw,r,pw_x,pw_y,pw_z,pw)
   aos = ao_value(1,r)
   call give_all_mos_at_r(r,mos_array)
   i_vec = 1
   call give_complex_function_at_r(pw_eigvec(1,i_vec),r,pw)
   val_1 = dsqrt(IMAG(pw)**2.d0 +  REAL(pw)**2.d0 )
   i_vec = 300
   call give_complex_function_at_r(pw_eigvec(1,i_vec),r,pw)
   val_2 = dsqrt(IMAG(pw)**2.d0 +  REAL(pw)**2.d0 )
   i_vec = 15625
   call give_complex_function_at_r(pw_eigvec(1,i_vec),r,pw)
   val_3 = dsqrt(IMAG(pw)**2.d0 +  REAL(pw)**2.d0 )
   write(35,'(F10.5,X,10(F20.16,X))')r(1),val_1,val_2,val_3,mos_array(77)
   r(1) += dx
  enddo
end

subroutine routine_2
 implicit none
 integer :: i,nx,n
 double precision :: Lx,alpha,omega
 complex*8 :: tf

 n = 0

 alpha = 1.d0
 omega = delta_k_pw(1)
 Lx = pw_box_x
 print*,'omega = ',omega
 print*,''
 do i = 1,6
  nx = 10**i
  call numerical_four_series_pol_gauss(alpha,n,omega,Lx,nx,tf)
  print*,nx,'tf',tf
 enddo
 stop
 print*,''
 print*,''
 omega = dble(pw_nx-1) * delta_k_pw(1)
 print*,'omega = ',omega
 print*,''
 do i = 1, 6
  nx = 10*i
  call numerical_four_series_pol_gauss(alpha,n,omega,Lx,nx,tf)
  print*,nx,'tf',tf
 enddo

end

double precision function int_cos(Lx,beta,dx)
 implicit none
 double precision, intent(in) :: Lx,beta,dx
 double precision :: x
 integer :: i , nx
 nx = int(Lx/dx)+1
 print*,'nx = ',nx
 int_cos = 0.d0
 x = 0.d0
 do i = 1, nx
  int_cos += dcos(beta*x)
  x += dx
 enddo
 int_cos *= dx

end
