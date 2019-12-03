
subroutine give_pw_i_at_r(i,r,pw_x,pw_y,pw_z,pw)
 implicit none
 integer, intent(in) :: i
 double precision, intent(in) :: r(3)
 complex*8, intent(out) :: pw_x, pw_y, pw_z
 complex*8, intent(out) :: pw
 double precision :: a,b
 include 'constants.include.F'
 a = dcos(d_list_pw(1,i) * delta_k_pw(1) * r(1))
 b = dsin(d_list_pw(1,i) * delta_k_pw(1) * r(1))
 pw_x = dcmplx(norm_box_pw(1)*a,norm_box_pw(1)*b)

 a = dcos(d_list_pw(2,i) * delta_k_pw(2) * r(2))
 b = dsin(d_list_pw(2,i) * delta_k_pw(2) * r(2))
 pw_y = dcmplx(norm_box_pw(2)*a,norm_box_pw(2)*b)
 
 a = dcos(d_list_pw(3,i) * delta_k_pw(3) * r(3))
 b = dsin(d_list_pw(3,i) * delta_k_pw(3) * r(3))
 pw_z = dcmplx(norm_box_pw(3)*a,norm_box_pw(3)*b)

 pw = pw_x * pw_y * pw_z
 
end

subroutine give_complex_function_at_r(array_coef,r,pw)
 implicit none
 double precision, intent(in) :: array_coef(n_pw_basis),r(3)
 complex*8, intent(out) :: pw
 complex*8 :: pw_x, pw_y, pw_z,pw_tmp
 pw = dcmplx(0.d0,0.d0)
 integer :: i
 do i = 1, n_pw_basis
  call give_pw_i_at_r(i,r,pw_x,pw_y,pw_z,pw_tmp)
  pw += pw_tmp * array_coef(i)
 enddo
end
